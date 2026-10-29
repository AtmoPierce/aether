#include <cuda_runtime.h>
#include <cuda_fp16.h>
#include <cuda_fp4.h>
#include <cuda_fp6.h>
#include <cuda_fp8.h>

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <csignal>
#include <cmath>
#include <cctype>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <utility>
#include <vector>

namespace {

struct Range {
    double lo;
    double hi;
};

struct Args {
    std::vector<std::string> ops;
    std::vector<std::string> modes;
    std::vector<std::string> precisions;
    std::string csv_path = "target/cuda_perf_raw.csv";
    int runs = 10;
    int launches = 100;
    std::size_t elements = 1u << 20;
    int block_size = 256;
    int power_sample_ms = 0;
    int power_poll_ms = 20;
};

struct Record {
    std::string op;
    std::string mode;
    std::string precision;
    int run;
    std::size_t elements;
    std::size_t launches;
    double time_elapsed_s;
    double ns_per_launch;
    double ns_per_element;
    double elements_per_second;
    double sample_output;
    double average_power_watts;
    std::string power_source;
};

struct PowerSamplerProcess {
    pid_t pid = -1;
    std::filesystem::path output_path;
    std::string power_source;
};

struct PowerMeasurement {
    double average_power_watts = std::numeric_limits<double>::quiet_NaN();
    std::string power_source;
};

[[noreturn]] void fail(const std::string& message) {
    std::cerr << "cuda_perf_main: " << message << '\n';
    std::exit(1);
}

void cuda_check(cudaError_t status, const char* context) {
    if (status != cudaSuccess) {
        std::ostringstream oss;
        oss << context << ": " << cudaGetErrorString(status);
        fail(oss.str());
    }
}

std::vector<std::string> split_csv_list(const std::string& value) {
    std::vector<std::string> out;
    std::stringstream ss(value);
    std::string item;
    while (std::getline(ss, item, ',')) {
        auto first = item.find_first_not_of(" \t\n\r");
        if (first == std::string::npos) {
            continue;
        }
        auto last = item.find_last_not_of(" \t\n\r");
        out.push_back(item.substr(first, last - first + 1));
    }
    return out;
}

bool contains_value(const std::vector<std::string>& values, const std::string& needle) {
    return std::find(values.begin(), values.end(), needle) != values.end();
}

std::string trim_copy(std::string value) {
    const auto first = value.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) {
        return "";
    }
    const auto last = value.find_last_not_of(" \t\n\r");
    return value.substr(first, last - first + 1);
}

bool nvidia_smi_available() {
    static const bool available = []() {
        return std::system("command -v nvidia-smi >/dev/null 2>&1") == 0;
    }();
    return available;
}

std::string gpu_selector_for_nvidia_smi() {
    const char* visible_devices = std::getenv("CUDA_VISIBLE_DEVICES");
    if (visible_devices == nullptr) {
        return "0";
    }

    std::string selector = visible_devices;
    const auto comma = selector.find(',');
    if (comma != std::string::npos) {
        selector.erase(comma);
    }

    selector = trim_copy(selector);
    return selector.empty() ? "0" : selector;
}

bool power_sampling_enabled(const Args& args) {
    return args.power_sample_ms > 0 && nvidia_smi_available();
}

std::size_t checked_mul(std::size_t lhs, std::size_t rhs, const char* label) {
    if (rhs != 0 && lhs > std::numeric_limits<std::size_t>::max() / rhs) {
        fail(std::string(label) + " overflow");
    }
    return lhs * rhs;
}

bool parse_power_draw_line(const std::string& line, double* watts) {
    const std::string trimmed = trim_copy(line);
    if (trimmed.empty()) {
        return false;
    }

    char* end = nullptr;
    errno = 0;
    const double parsed = std::strtod(trimmed.c_str(), &end);
    if (end == trimmed.c_str() || errno != 0) {
        return false;
    }

    while (*end != '\0' && std::isspace(static_cast<unsigned char>(*end)) != 0) {
        ++end;
    }

    if (*end != '\0' || !std::isfinite(parsed) || parsed < 0.0) {
        return false;
    }

    *watts = parsed;
    return true;
}

PowerSamplerProcess start_power_sampler(const Args& args) {
    if (!power_sampling_enabled(args)) {
        return {};
    }

    char output_template[] = "/tmp/fpx_cuda_power_XXXXXX";
    const int output_fd = mkstemp(output_template);
    if (output_fd < 0) {
        return {};
    }

    const pid_t pid = fork();
    if (pid < 0) {
        close(output_fd);
        std::filesystem::remove(output_template);
        return {};
    }

    if (pid == 0) {
        if (dup2(output_fd, STDOUT_FILENO) < 0) {
            _exit(127);
        }

        const int null_fd = open("/dev/null", O_WRONLY);
        if (null_fd >= 0) {
            dup2(null_fd, STDERR_FILENO);
            close(null_fd);
        }

        close(output_fd);

        const std::string gpu_selector = gpu_selector_for_nvidia_smi();
        const std::string loop_arg = "--loop-ms=" + std::to_string(args.power_poll_ms);
        execlp(
            "nvidia-smi",
            "nvidia-smi",
            "-i",
            gpu_selector.c_str(),
            "--query-gpu=power.draw",
            "--format=csv,noheader,nounits",
            loop_arg.c_str(),
            static_cast<char*>(nullptr));
        _exit(127);
    }

    close(output_fd);
    return {pid, output_template, "nvidia-smi:power.draw"};
}

std::uintmax_t power_sampler_offset(const PowerSamplerProcess& sampler) {
    if (sampler.pid <= 0) {
        return 0;
    }

    std::error_code error;
    return std::filesystem::file_size(sampler.output_path, error);
}

PowerMeasurement read_power_samples_between(
    const PowerSamplerProcess& sampler,
    std::uintmax_t start_offset,
    std::uintmax_t end_offset) {
    if (sampler.pid <= 0 || end_offset <= start_offset) {
        return {};
    }

    std::ifstream input(sampler.output_path);
    if (!input) {
        return {};
    }

    std::string line;
    std::uintmax_t line_offset = 0;
    double power_sum = 0.0;
    std::size_t samples = 0;
    while (std::getline(input, line)) {
        const std::uintmax_t next_offset = line_offset + line.size() + 1;
        if (line_offset >= start_offset && line_offset < end_offset) {
            double watts = 0.0;
            if (parse_power_draw_line(line, &watts)) {
                power_sum += watts;
                ++samples;
            }
        }
        line_offset = next_offset;
        if (line_offset >= end_offset) {
            break;
        }
    }

    if (samples == 0) {
        return {};
    }

    return {
        power_sum / static_cast<double>(samples),
        sampler.power_source,
    };
}

void stop_power_sampler(const PowerSamplerProcess& sampler) {
    if (sampler.pid <= 0) {
        return;
    }

    if (kill(sampler.pid, SIGTERM) != 0 && errno != ESRCH) {
        waitpid(sampler.pid, nullptr, 0);
        std::filesystem::remove(sampler.output_path);
        return;
    }

    waitpid(sampler.pid, nullptr, 0);
    std::filesystem::remove(sampler.output_path);
}

int timed_graph_repeats(const Args& args, float graph_elapsed_ms) {
    if (!power_sampling_enabled(args) || graph_elapsed_ms <= 0.0f) {
        return 1;
    }

    const double repeats =
        std::ceil(static_cast<double>(args.power_sample_ms) / static_cast<double>(graph_elapsed_ms));
    return std::max(1, static_cast<int>(repeats));
}

void print_usage(const char* argv0) {
    std::cout
        << "usage: " << argv0 << " [options]\n\n"
        << "Options:\n"
        << "  --only <op1,op2,...>          Limit benchmark to specific ops\n"
        << "  --modes <m1,m2,...>           Limit benchmark to specific modes\n"
        << "  --precisions <p1,p2,...>      Limit benchmark to specific precisions\n"
        << "  --runs <n>                    Number of measured runs per op/precision (default: 10)\n"
        << "  --launches <n>                Kernel launches per measured run (default: 100)\n"
        << "  --elements <n>                Elements processed per kernel launch (default: 1048576)\n"
        << "  --block-size <n>              CUDA block size (default: 256)\n"
        << "  --power-sample-ms <n>         Target minimum timed interval per row for nvidia-smi power sampling (default: 0, disabled)\n"
        << "  --power-poll-ms <n>           nvidia-smi power polling interval in milliseconds (default: 20)\n"
        << "  --csv <path>                  Output raw performance CSV (default: target/cuda_perf_raw.csv)\n"
        << "  --help                        Show this help message\n\n"
        << "Reduced-domain modes (all ops): fpx_kernel_horner,nat_kernel\n"
        << "Full-kernel modes (log only): fpx,nat\n"
        << "Supported precisions: f16,f32,f64,fp4_e2m1,fp6_e3m2,fp8_e4m3,fp8_e5m2\n";
}

Args parse_args(int argc, char** argv) {
    Args args;
    for (int index = 1; index < argc; ++index) {
        const std::string flag = argv[index];
        auto require_value = [&](const char* name) -> std::string {
            if (index + 1 >= argc) {
                fail(std::string("missing value for ") + name);
            }
            return argv[++index];
        };

        if (flag == "--only") {
            args.ops = split_csv_list(require_value("--only"));
        } else if (flag == "--modes") {
            args.modes = split_csv_list(require_value("--modes"));
        } else if (flag == "--precisions") {
            args.precisions = split_csv_list(require_value("--precisions"));
        } else if (flag == "--runs") {
            args.runs = std::stoi(require_value("--runs"));
        } else if (flag == "--launches") {
            args.launches = std::stoi(require_value("--launches"));
        } else if (flag == "--elements") {
            args.elements = static_cast<std::size_t>(std::stoull(require_value("--elements")));
        } else if (flag == "--block-size") {
            args.block_size = std::stoi(require_value("--block-size"));
        } else if (flag == "--power-sample-ms") {
            args.power_sample_ms = std::stoi(require_value("--power-sample-ms"));
        } else if (flag == "--power-poll-ms") {
            args.power_poll_ms = std::stoi(require_value("--power-poll-ms"));
        } else if (flag == "--csv") {
            args.csv_path = require_value("--csv");
        } else if (flag == "--help" || flag == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        } else {
            fail(std::string("unknown argument: ") + flag);
        }
    }

    if (args.runs <= 0) {
        fail("--runs must be positive");
    }
    if (args.launches <= 0) {
        fail("--launches must be positive");
    }
    if (args.elements == 0) {
        fail("--elements must be positive");
    }
    if (args.block_size <= 0) {
        fail("--block-size must be positive");
    }
    if (args.power_sample_ms < 0) {
        fail("--power-sample-ms must be non-negative");
    }
    if (args.power_poll_ms <= 0) {
        fail("--power-poll-ms must be positive");
    }

    return args;
}

Range full_domain_for_op(const std::string& op) {
    if (op == "log") {
        return {0.25, 16.0};
    }
    fail(std::string("no full-kernel benchmark domain configured for op ") + op);
}

template <typename T>
struct ScalarCodec;

template <>
struct ScalarCodec<float> {
    __host__ __device__ static float from_double(double value) { return static_cast<float>(value); }
    __host__ __device__ static float from_float(float value) { return value; }
    __host__ __device__ static double to_double(float value) { return static_cast<double>(value); }
};

template <>
struct ScalarCodec<double> {
    __host__ __device__ static double from_double(double value) { return value; }
    __host__ __device__ static double from_float(float value) { return static_cast<double>(value); }
    __host__ __device__ static double to_double(double value) { return value; }
};

template <>
struct ScalarCodec<__half> {
    __host__ __device__ static __half from_double(double value) { return __half(value); }
    __host__ __device__ static __half from_float(float value) { return __half(value); }
    __host__ __device__ static double to_double(__half value) { return static_cast<double>(static_cast<float>(value)); }
};

template <>
struct ScalarCodec<__nv_fp4_e2m1> {
    __host__ __device__ static __nv_fp4_e2m1 from_double(double value) { return __nv_fp4_e2m1(value); }
    __host__ __device__ static __nv_fp4_e2m1 from_float(float value) { return __nv_fp4_e2m1(value); }
    __host__ __device__ static double to_double(__nv_fp4_e2m1 value) { return static_cast<double>(static_cast<float>(value)); }
};

template <>
struct ScalarCodec<__nv_fp6_e3m2> {
    __host__ __device__ static __nv_fp6_e3m2 from_double(double value) { return __nv_fp6_e3m2(value); }
    __host__ __device__ static __nv_fp6_e3m2 from_float(float value) { return __nv_fp6_e3m2(value); }
    __host__ __device__ static double to_double(__nv_fp6_e3m2 value) { return static_cast<double>(static_cast<float>(value)); }
};

template <>
struct ScalarCodec<__nv_fp8_e4m3> {
    __host__ __device__ static __nv_fp8_e4m3 from_double(double value) { return __nv_fp8_e4m3(value); }
    __host__ __device__ static __nv_fp8_e4m3 from_float(float value) { return __nv_fp8_e4m3(value); }
    __host__ __device__ static double to_double(__nv_fp8_e4m3 value) { return static_cast<double>(static_cast<float>(value)); }
};

template <>
struct ScalarCodec<__nv_fp8_e5m2> {
    __host__ __device__ static __nv_fp8_e5m2 from_double(double value) { return __nv_fp8_e5m2(value); }
    __host__ __device__ static __nv_fp8_e5m2 from_float(float value) { return __nv_fp8_e5m2(value); }
    __host__ __device__ static double to_double(__nv_fp8_e5m2 value) { return static_cast<double>(static_cast<float>(value)); }
};

template <typename T>
using LaunchFn = void (*)(const T*, T*, std::size_t, dim3, dim3, cudaStream_t);

template <typename T>
struct BenchSpec {
    const char* op;
    const char* mode;
    Range range;
    LaunchFn<T> launch;
};

#define FPX_CUDA_OPS(X) \
    X(acos)             \
    X(asin)             \
    X(atan)             \
    X(cos)              \
    X(cosh)             \
    X(exp)              \
    X(expm1)            \
    X(log)              \
    X(log1p)            \
    X(sin)              \
    X(sinh)             \
    X(tan)              \
    X(tanh)

#define DECLARE_REDUCED_KERNELS(op)                                                                \
    extern "C" __global__ void op##_reduced_kernel_double(const double*, double*, std::size_t); \
    extern "C" __global__ void op##_reduced_kernel_float(const float*, float*, std::size_t);    \
    extern "C" __global__ void op##_reduced_kernel_f16(const __half*, __half*, std::size_t);    \
    extern "C" __global__ void op##_reduced_kernel_fp4_e2m1(                                   \
        const __nv_fp4_e2m1*, __nv_fp4_e2m1*, std::size_t);                                      \
    extern "C" __global__ void op##_reduced_kernel_fp6_e3m2(                                   \
        const __nv_fp6_e3m2*, __nv_fp6_e3m2*, std::size_t);                                      \
    extern "C" __global__ void op##_reduced_kernel_fp8_e4m3(                                   \
        const __nv_fp8_e4m3*, __nv_fp8_e4m3*, std::size_t);                                      \
    extern "C" __global__ void op##_reduced_kernel_fp8_e5m2(                                   \
        const __nv_fp8_e5m2*, __nv_fp8_e5m2*, std::size_t);                                      \
    extern "C" double op##_reduced_domain_lo();                                                 \
    extern "C" double op##_reduced_domain_hi();

FPX_CUDA_OPS(DECLARE_REDUCED_KERNELS)

#define DECLARE_FULL_KERNELS(op)                                                               \
    extern "C" __global__ void op##_kernel_double(const double*, double*, std::size_t);      \
    extern "C" __global__ void op##_kernel_float(const float*, float*, std::size_t);         \
    extern "C" __global__ void op##_kernel_f16(const __half*, __half*, std::size_t);         \
    extern "C" __global__ void op##_kernel_fp4_e2m1(                                         \
        const __nv_fp4_e2m1*, __nv_fp4_e2m1*, std::size_t);                                   \
    extern "C" __global__ void op##_kernel_fp6_e3m2(                                         \
        const __nv_fp6_e3m2*, __nv_fp6_e3m2*, std::size_t);                                   \
    extern "C" __global__ void op##_kernel_fp8_e4m3(                                         \
        const __nv_fp8_e4m3*, __nv_fp8_e4m3*, std::size_t);                                   \
    extern "C" __global__ void op##_kernel_fp8_e5m2(                                         \
        const __nv_fp8_e5m2*, __nv_fp8_e5m2*, std::size_t);

DECLARE_FULL_KERNELS(log)

Range reduced_domain_for_op(const std::string& op) {
#define RETURN_REDUCED_DOMAIN(candidate)                                      \
    if (op == #candidate) {                                                   \
        return {candidate##_reduced_domain_lo(), candidate##_reduced_domain_hi()}; \
    }
    FPX_CUDA_OPS(RETURN_REDUCED_DOMAIN)
#undef RETURN_REDUCED_DOMAIN
    fail(std::string("no reduced kernel domain configured for op ") + op);
}

#define DEFINE_NATIVE_DEVICE_FUNCS(op, float_fn, double_fn)                                            \
    __device__ __forceinline__ float native_##op(float value) { return float_fn(value); }            \
    __device__ __forceinline__ double native_##op(double value) { return double_fn(value); }         \
    __device__ __forceinline__ __half native_##op(__half value) {                                     \
        return ScalarCodec<__half>::from_float(float_fn(static_cast<float>(value)));                  \
    }                                                                                                  \
    __device__ __forceinline__ __nv_fp4_e2m1 native_##op(__nv_fp4_e2m1 value) {                       \
        return ScalarCodec<__nv_fp4_e2m1>::from_float(float_fn(static_cast<float>(value)));           \
    }                                                                                                  \
    __device__ __forceinline__ __nv_fp6_e3m2 native_##op(__nv_fp6_e3m2 value) {                       \
        return ScalarCodec<__nv_fp6_e3m2>::from_float(float_fn(static_cast<float>(value)));           \
    }                                                                                                  \
    __device__ __forceinline__ __nv_fp8_e4m3 native_##op(__nv_fp8_e4m3 value) {                       \
        return ScalarCodec<__nv_fp8_e4m3>::from_float(float_fn(static_cast<float>(value)));           \
    }                                                                                                  \
    __device__ __forceinline__ __nv_fp8_e5m2 native_##op(__nv_fp8_e5m2 value) {                       \
        return ScalarCodec<__nv_fp8_e5m2>::from_float(float_fn(static_cast<float>(value)));           \
    }

DEFINE_NATIVE_DEVICE_FUNCS(acos, acosf, acos)
DEFINE_NATIVE_DEVICE_FUNCS(asin, asinf, asin)
DEFINE_NATIVE_DEVICE_FUNCS(atan, atanf, atan)
DEFINE_NATIVE_DEVICE_FUNCS(cos, cosf, cos)
DEFINE_NATIVE_DEVICE_FUNCS(cosh, coshf, cosh)
DEFINE_NATIVE_DEVICE_FUNCS(exp, expf, exp)
DEFINE_NATIVE_DEVICE_FUNCS(expm1, expm1f, expm1)
DEFINE_NATIVE_DEVICE_FUNCS(log, logf, log)
DEFINE_NATIVE_DEVICE_FUNCS(log1p, log1pf, log1p)
DEFINE_NATIVE_DEVICE_FUNCS(sin, sinf, sin)
DEFINE_NATIVE_DEVICE_FUNCS(sinh, sinhf, sinh)
DEFINE_NATIVE_DEVICE_FUNCS(tan, tanf, tan)
DEFINE_NATIVE_DEVICE_FUNCS(tanh, tanhf, tanh)

#define DEFINE_REDUCED_LAUNCHERS(op)                                                                   \
    static void launch_##op##_f64(const double* xs,                                                   \
                                  double* ys,                                                         \
                                  std::size_t n,                                                      \
                                  dim3 grid,                                                          \
                                  dim3 blk,                                                           \
                                  cudaStream_t stream) {                                              \
        op##_reduced_kernel_double<<<grid, blk, 0, stream>>>(xs, ys, n);                             \
    }                                                                                                  \
    static void launch_##op##_f32(const float* xs,                                                    \
                                  float* ys,                                                          \
                                  std::size_t n,                                                      \
                                  dim3 grid,                                                          \
                                  dim3 blk,                                                           \
                                  cudaStream_t stream) {                                              \
        op##_reduced_kernel_float<<<grid, blk, 0, stream>>>(xs, ys, n);                              \
    }                                                                                                  \
    static void launch_##op##_f16(const __half* xs,                                                   \
                                  __half* ys,                                                         \
                                  std::size_t n,                                                      \
                                  dim3 grid,                                                          \
                                  dim3 blk,                                                           \
                                  cudaStream_t stream) {                                              \
        op##_reduced_kernel_f16<<<grid, blk, 0, stream>>>(xs, ys, n);                                \
    }                                                                                                  \
    static void launch_##op##_fp4(const __nv_fp4_e2m1* xs,                                            \
                                  __nv_fp4_e2m1* ys,                                                  \
                                  std::size_t n,                                                      \
                                  dim3 grid,                                                          \
                                  dim3 blk,                                                           \
                                  cudaStream_t stream) {                                              \
        op##_reduced_kernel_fp4_e2m1<<<grid, blk, 0, stream>>>(xs, ys, n);                           \
    }                                                                                                  \
    static void launch_##op##_fp6(const __nv_fp6_e3m2* xs,                                            \
                                  __nv_fp6_e3m2* ys,                                                  \
                                  std::size_t n,                                                      \
                                  dim3 grid,                                                          \
                                  dim3 blk,                                                           \
                                  cudaStream_t stream) {                                              \
        op##_reduced_kernel_fp6_e3m2<<<grid, blk, 0, stream>>>(xs, ys, n);                           \
    }                                                                                                  \
    static void launch_##op##_fp8e4(const __nv_fp8_e4m3* xs,                                          \
                                    __nv_fp8_e4m3* ys,                                                \
                                    std::size_t n,                                                    \
                                    dim3 grid,                                                        \
                                    dim3 blk,                                                         \
                                    cudaStream_t stream) {                                            \
        op##_reduced_kernel_fp8_e4m3<<<grid, blk, 0, stream>>>(xs, ys, n);                           \
    }                                                                                                  \
    static void launch_##op##_fp8e5(const __nv_fp8_e5m2* xs,                                          \
                                    __nv_fp8_e5m2* ys,                                                \
                                    std::size_t n,                                                    \
                                    dim3 grid,                                                        \
                                    dim3 blk,                                                         \
                                    cudaStream_t stream) {                                            \
        op##_reduced_kernel_fp8_e5m2<<<grid, blk, 0, stream>>>(xs, ys, n);                           \
    }

FPX_CUDA_OPS(DEFINE_REDUCED_LAUNCHERS)

#define DEFINE_FULL_LAUNCHERS(op)                                                                      \
    static void launch_full_##op##_f64(const double* xs,                                             \
                                       double* ys,                                                   \
                                       std::size_t n,                                                \
                                       dim3 grid,                                                    \
                                       dim3 blk,                                                     \
                                       cudaStream_t stream) {                                        \
        op##_kernel_double<<<grid, blk, 0, stream>>>(xs, ys, n);                                    \
    }                                                                                                 \
    static void launch_full_##op##_f32(const float* xs,                                              \
                                       float* ys,                                                    \
                                       std::size_t n,                                                \
                                       dim3 grid,                                                    \
                                       dim3 blk,                                                     \
                                       cudaStream_t stream) {                                        \
        op##_kernel_float<<<grid, blk, 0, stream>>>(xs, ys, n);                                     \
    }                                                                                                 \
    static void launch_full_##op##_f16(const __half* xs,                                             \
                                       __half* ys,                                                   \
                                       std::size_t n,                                                \
                                       dim3 grid,                                                    \
                                       dim3 blk,                                                     \
                                       cudaStream_t stream) {                                        \
        op##_kernel_f16<<<grid, blk, 0, stream>>>(xs, ys, n);                                       \
    }                                                                                                 \
    static void launch_full_##op##_fp4(const __nv_fp4_e2m1* xs,                                     \
                                       __nv_fp4_e2m1* ys,                                           \
                                       std::size_t n,                                                \
                                       dim3 grid,                                                    \
                                       dim3 blk,                                                     \
                                       cudaStream_t stream) {                                        \
        op##_kernel_fp4_e2m1<<<grid, blk, 0, stream>>>(xs, ys, n);                                  \
    }                                                                                                 \
    static void launch_full_##op##_fp6(const __nv_fp6_e3m2* xs,                                     \
                                       __nv_fp6_e3m2* ys,                                           \
                                       std::size_t n,                                                \
                                       dim3 grid,                                                    \
                                       dim3 blk,                                                     \
                                       cudaStream_t stream) {                                        \
        op##_kernel_fp6_e3m2<<<grid, blk, 0, stream>>>(xs, ys, n);                                  \
    }                                                                                                 \
    static void launch_full_##op##_fp8e4(const __nv_fp8_e4m3* xs,                                   \
                                         __nv_fp8_e4m3* ys,                                         \
                                         std::size_t n,                                              \
                                         dim3 grid,                                                  \
                                         dim3 blk,                                                   \
                                         cudaStream_t stream) {                                      \
        op##_kernel_fp8_e4m3<<<grid, blk, 0, stream>>>(xs, ys, n);                                  \
    }                                                                                                 \
    static void launch_full_##op##_fp8e5(const __nv_fp8_e5m2* xs,                                   \
                                         __nv_fp8_e5m2* ys,                                         \
                                         std::size_t n,                                              \
                                         dim3 grid,                                                  \
                                         dim3 blk,                                                   \
                                         cudaStream_t stream) {                                      \
        op##_kernel_fp8_e5m2<<<grid, blk, 0, stream>>>(xs, ys, n);                                  \
    }

DEFINE_FULL_LAUNCHERS(log)

#define DEFINE_NATIVE_LAUNCHERS(op)                                                                    \
    __global__ void native_##op##_kernel_double(const double* xs, double* ys, std::size_t n) {      \
        const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;   \
        if (index < n) {                                                                               \
            ys[index] = native_##op(xs[index]);                                                        \
        }                                                                                              \
    }                                                                                                  \
    __global__ void native_##op##_kernel_float(const float* xs, float* ys, std::size_t n) {          \
        const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;   \
        if (index < n) {                                                                               \
            ys[index] = native_##op(xs[index]);                                                        \
        }                                                                                              \
    }                                                                                                  \
    __global__ void native_##op##_kernel_f16(const __half* xs, __half* ys, std::size_t n) {          \
        const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;   \
        if (index < n) {                                                                               \
            ys[index] = native_##op(xs[index]);                                                        \
        }                                                                                              \
    }                                                                                                  \
    __global__ void native_##op##_kernel_fp4_e2m1(const __nv_fp4_e2m1* xs,                            \
                                                   __nv_fp4_e2m1* ys,                                  \
                                                   std::size_t n) {                                    \
        const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;   \
        if (index < n) {                                                                               \
            ys[index] = native_##op(xs[index]);                                                        \
        }                                                                                              \
    }                                                                                                  \
    __global__ void native_##op##_kernel_fp6_e3m2(const __nv_fp6_e3m2* xs,                            \
                                                   __nv_fp6_e3m2* ys,                                  \
                                                   std::size_t n) {                                    \
        const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;   \
        if (index < n) {                                                                               \
            ys[index] = native_##op(xs[index]);                                                        \
        }                                                                                              \
    }                                                                                                  \
    __global__ void native_##op##_kernel_fp8_e4m3(const __nv_fp8_e4m3* xs,                            \
                                                   __nv_fp8_e4m3* ys,                                  \
                                                   std::size_t n) {                                    \
        const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;   \
        if (index < n) {                                                                               \
            ys[index] = native_##op(xs[index]);                                                        \
        }                                                                                              \
    }                                                                                                  \
    __global__ void native_##op##_kernel_fp8_e5m2(const __nv_fp8_e5m2* xs,                            \
                                                   __nv_fp8_e5m2* ys,                                  \
                                                   std::size_t n) {                                    \
        const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;   \
        if (index < n) {                                                                               \
            ys[index] = native_##op(xs[index]);                                                        \
        }                                                                                              \
    }                                                                                                  \
    static void launch_native_##op##_f64(const double* xs,                                            \
                                         double* ys,                                                  \
                                         std::size_t n,                                               \
                                         dim3 grid,                                                   \
                                         dim3 blk,                                                    \
                                         cudaStream_t stream) {                                       \
        native_##op##_kernel_double<<<grid, blk, 0, stream>>>(xs, ys, n);                            \
    }                                                                                                  \
    static void launch_native_##op##_f32(const float* xs,                                             \
                                         float* ys,                                                   \
                                         std::size_t n,                                               \
                                         dim3 grid,                                                   \
                                         dim3 blk,                                                    \
                                         cudaStream_t stream) {                                       \
        native_##op##_kernel_float<<<grid, blk, 0, stream>>>(xs, ys, n);                             \
    }                                                                                                  \
    static void launch_native_##op##_f16(const __half* xs,                                            \
                                         __half* ys,                                                  \
                                         std::size_t n,                                               \
                                         dim3 grid,                                                   \
                                         dim3 blk,                                                    \
                                         cudaStream_t stream) {                                       \
        native_##op##_kernel_f16<<<grid, blk, 0, stream>>>(xs, ys, n);                               \
    }                                                                                                  \
    static void launch_native_##op##_fp4(const __nv_fp4_e2m1* xs,                                     \
                                         __nv_fp4_e2m1* ys,                                           \
                                         std::size_t n,                                               \
                                         dim3 grid,                                                   \
                                         dim3 blk,                                                    \
                                         cudaStream_t stream) {                                       \
        native_##op##_kernel_fp4_e2m1<<<grid, blk, 0, stream>>>(xs, ys, n);                          \
    }                                                                                                  \
    static void launch_native_##op##_fp6(const __nv_fp6_e3m2* xs,                                     \
                                         __nv_fp6_e3m2* ys,                                           \
                                         std::size_t n,                                               \
                                         dim3 grid,                                                   \
                                         dim3 blk,                                                    \
                                         cudaStream_t stream) {                                       \
        native_##op##_kernel_fp6_e3m2<<<grid, blk, 0, stream>>>(xs, ys, n);                          \
    }                                                                                                  \
    static void launch_native_##op##_fp8e4(const __nv_fp8_e4m3* xs,                                   \
                                           __nv_fp8_e4m3* ys,                                         \
                                           std::size_t n,                                             \
                                           dim3 grid,                                                 \
                                           dim3 blk,                                                  \
                                           cudaStream_t stream) {                                     \
        native_##op##_kernel_fp8_e4m3<<<grid, blk, 0, stream>>>(xs, ys, n);                          \
    }                                                                                                  \
    static void launch_native_##op##_fp8e5(const __nv_fp8_e5m2* xs,                                   \
                                           __nv_fp8_e5m2* ys,                                         \
                                           std::size_t n,                                             \
                                           dim3 grid,                                                 \
                                           dim3 blk,                                                  \
                                           cudaStream_t stream) {                                     \
        native_##op##_kernel_fp8_e5m2<<<grid, blk, 0, stream>>>(xs, ys, n);                          \
    }

FPX_CUDA_OPS(DEFINE_NATIVE_LAUNCHERS)

template <typename T>
std::vector<BenchSpec<T>> specs_for_precision();

#define SPEC_F64(op)                                                                                  \
    BenchSpec<double>{#op, "fpx_kernel_horner", reduced_domain_for_op(#op), launch_##op##_f64},   \
        BenchSpec<double>{#op, "nat_kernel", reduced_domain_for_op(#op), launch_native_##op##_f64},
template <>
std::vector<BenchSpec<double>> specs_for_precision<double>() {
    return {
        FPX_CUDA_OPS(SPEC_F64)
        BenchSpec<double>{"log", "fpx", full_domain_for_op("log"), launch_full_log_f64},
        BenchSpec<double>{"log", "nat", full_domain_for_op("log"), launch_native_log_f64},
    };
}

#define SPEC_F32(op)                                                                                 \
    BenchSpec<float>{#op, "fpx_kernel_horner", reduced_domain_for_op(#op), launch_##op##_f32},   \
        BenchSpec<float>{#op, "nat_kernel", reduced_domain_for_op(#op), launch_native_##op##_f32},
template <>
std::vector<BenchSpec<float>> specs_for_precision<float>() {
    return {
        FPX_CUDA_OPS(SPEC_F32)
        BenchSpec<float>{"log", "fpx", full_domain_for_op("log"), launch_full_log_f32},
        BenchSpec<float>{"log", "nat", full_domain_for_op("log"), launch_native_log_f32},
    };
}

#define SPEC_F16(op)                                                                                  \
    BenchSpec<__half>{#op, "fpx_kernel_horner", reduced_domain_for_op(#op), launch_##op##_f16},  \
        BenchSpec<__half>{#op, "nat_kernel", reduced_domain_for_op(#op), launch_native_##op##_f16},
template <>
std::vector<BenchSpec<__half>> specs_for_precision<__half>() {
    return {
        FPX_CUDA_OPS(SPEC_F16)
        BenchSpec<__half>{"log", "fpx", full_domain_for_op("log"), launch_full_log_f16},
        BenchSpec<__half>{"log", "nat", full_domain_for_op("log"), launch_native_log_f16},
    };
}

#define SPEC_FP4(op)                                                                                         \
    BenchSpec<__nv_fp4_e2m1>{#op, "fpx_kernel_horner", reduced_domain_for_op(#op), launch_##op##_fp4},   \
        BenchSpec<__nv_fp4_e2m1>{#op, "nat_kernel", reduced_domain_for_op(#op), launch_native_##op##_fp4},
template <>
std::vector<BenchSpec<__nv_fp4_e2m1>> specs_for_precision<__nv_fp4_e2m1>() {
    return {
        FPX_CUDA_OPS(SPEC_FP4)
        BenchSpec<__nv_fp4_e2m1>{"log", "fpx", full_domain_for_op("log"), launch_full_log_fp4},
        BenchSpec<__nv_fp4_e2m1>{"log", "nat", full_domain_for_op("log"), launch_native_log_fp4},
    };
}

#define SPEC_FP6(op)                                                                                         \
    BenchSpec<__nv_fp6_e3m2>{#op, "fpx_kernel_horner", reduced_domain_for_op(#op), launch_##op##_fp6},   \
        BenchSpec<__nv_fp6_e3m2>{#op, "nat_kernel", reduced_domain_for_op(#op), launch_native_##op##_fp6},
template <>
std::vector<BenchSpec<__nv_fp6_e3m2>> specs_for_precision<__nv_fp6_e3m2>() {
    return {
        FPX_CUDA_OPS(SPEC_FP6)
        BenchSpec<__nv_fp6_e3m2>{"log", "fpx", full_domain_for_op("log"), launch_full_log_fp6},
        BenchSpec<__nv_fp6_e3m2>{"log", "nat", full_domain_for_op("log"), launch_native_log_fp6},
    };
}

#define SPEC_FP8E4(op)                                                                                           \
    BenchSpec<__nv_fp8_e4m3>{#op, "fpx_kernel_horner", reduced_domain_for_op(#op), launch_##op##_fp8e4},      \
        BenchSpec<__nv_fp8_e4m3>{#op, "nat_kernel", reduced_domain_for_op(#op), launch_native_##op##_fp8e4},
template <>
std::vector<BenchSpec<__nv_fp8_e4m3>> specs_for_precision<__nv_fp8_e4m3>() {
    return {
        FPX_CUDA_OPS(SPEC_FP8E4)
        BenchSpec<__nv_fp8_e4m3>{"log", "fpx", full_domain_for_op("log"), launch_full_log_fp8e4},
        BenchSpec<__nv_fp8_e4m3>{"log", "nat", full_domain_for_op("log"), launch_native_log_fp8e4},
    };
}

#define SPEC_FP8E5(op)                                                                                           \
    BenchSpec<__nv_fp8_e5m2>{#op, "fpx_kernel_horner", reduced_domain_for_op(#op), launch_##op##_fp8e5},      \
        BenchSpec<__nv_fp8_e5m2>{#op, "nat_kernel", reduced_domain_for_op(#op), launch_native_##op##_fp8e5},
template <>
std::vector<BenchSpec<__nv_fp8_e5m2>> specs_for_precision<__nv_fp8_e5m2>() {
    return {
        FPX_CUDA_OPS(SPEC_FP8E5)
        BenchSpec<__nv_fp8_e5m2>{"log", "fpx", full_domain_for_op("log"), launch_full_log_fp8e5},
        BenchSpec<__nv_fp8_e5m2>{"log", "nat", full_domain_for_op("log"), launch_native_log_fp8e5},
    };
}

std::vector<std::string> supported_ops() {
#define OP_NAME(op) #op,
    return {FPX_CUDA_OPS(OP_NAME)};
#undef OP_NAME
}

std::vector<std::string> supported_precisions() {
    return {"f16", "f32", "f64", "fp4_e2m1", "fp6_e3m2", "fp8_e4m3", "fp8_e5m2"};
}

std::vector<std::string> supported_modes() {
    return {"fpx_kernel_horner", "nat_kernel", "fpx", "nat"};
}

void validate_filters(const Args& args) {
    const auto ops = supported_ops();
    const auto modes = supported_modes();
    const auto precisions = supported_precisions();

    for (const auto& op : args.ops) {
        if (!contains_value(ops, op)) {
            fail(std::string("unsupported op filter: ") + op);
        }
    }

    for (const auto& mode : args.modes) {
        if (!contains_value(modes, mode)) {
            fail(std::string("unsupported mode filter: ") + mode);
        }
    }

    for (const auto& precision : args.precisions) {
        if (!contains_value(precisions, precision)) {
            fail(std::string("unsupported precision filter: ") + precision);
        }
    }

    const bool reduced_mode_selected = args.modes.empty() ||
        contains_value(args.modes, "fpx_kernel_horner") || contains_value(args.modes, "nat_kernel");
    const bool log_selected = args.ops.empty() || contains_value(args.ops, "log");
    if (!reduced_mode_selected && !log_selected) {
        fail("full-kernel modes fpx and nat are available only for log");
    }
}

template <typename T>
Record benchmark_one(const BenchSpec<T>& spec,
                     const std::string& precision,
                     const Args& args,
                     int run,
                     const PowerSamplerProcess& power_sampler) {
    std::vector<T> host_input(args.elements);
    const double span = spec.range.hi - spec.range.lo;
    for (std::size_t index = 0; index < args.elements; ++index) {
        const double alpha = args.elements == 1
            ? 0.0
            : static_cast<double>(index) / static_cast<double>(args.elements - 1);
        host_input[index] = ScalarCodec<T>::from_double(spec.range.lo + span * alpha);
    }

    T* device_input = nullptr;
    T* device_output = nullptr;
    cudaStream_t stream = nullptr;
    cudaGraph_t graph = nullptr;
    cudaGraphExec_t graph_exec = nullptr;
    cuda_check(cudaMalloc(&device_input, sizeof(T) * args.elements), "cudaMalloc(device_input)");
    cuda_check(cudaMalloc(&device_output, sizeof(T) * args.elements), "cudaMalloc(device_output)");
    cuda_check(cudaMemcpy(device_input, host_input.data(), sizeof(T) * args.elements, cudaMemcpyHostToDevice),
               "cudaMemcpy(host->device)");
    cuda_check(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking), "cudaStreamCreateWithFlags");

    const dim3 block(static_cast<unsigned int>(args.block_size));
    const dim3 grid(static_cast<unsigned int>((args.elements + static_cast<std::size_t>(args.block_size) - 1) /
                                              static_cast<std::size_t>(args.block_size)));

    spec.launch(device_input, device_output, args.elements, grid, block, stream);
    cuda_check(cudaGetLastError(), "warmup launch");
    cuda_check(cudaStreamSynchronize(stream), "warmup synchronize");

    cuda_check(cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal), "cudaStreamBeginCapture");
    for (int launch = 0; launch < args.launches; ++launch) {
        spec.launch(device_input, device_output, args.elements, grid, block, stream);
    }
    cuda_check(cudaStreamEndCapture(stream, &graph), "cudaStreamEndCapture");
    cuda_check(cudaGraphInstantiate(&graph_exec, graph, 0), "cudaGraphInstantiate");

    cuda_check(cudaGraphLaunch(graph_exec, stream), "cudaGraphLaunch(warmup)");
    cuda_check(cudaStreamSynchronize(stream), "cudaGraphLaunch(warmup) synchronize");

    cudaEvent_t start = nullptr;
    cudaEvent_t stop = nullptr;
    cuda_check(cudaEventCreate(&start), "cudaEventCreate(start)");
    cuda_check(cudaEventCreate(&stop), "cudaEventCreate(stop)");

    cuda_check(cudaEventRecord(start, stream), "cudaEventRecord(start pilot)");
    cuda_check(cudaGraphLaunch(graph_exec, stream), "cudaGraphLaunch(pilot)");
    cuda_check(cudaEventRecord(stop, stream), "cudaEventRecord(stop pilot)");
    cuda_check(cudaEventSynchronize(stop), "cudaEventSynchronize(stop pilot)");

    float graph_elapsed_ms = 0.0f;
    cuda_check(cudaEventElapsedTime(&graph_elapsed_ms, start, stop), "cudaEventElapsedTime(pilot)");

    const int graph_repeats = timed_graph_repeats(args, graph_elapsed_ms);
    const std::size_t total_launches = checked_mul(
        static_cast<std::size_t>(args.launches),
        static_cast<std::size_t>(graph_repeats),
        "launch count");

    const std::uintmax_t power_sample_start = power_sampler_offset(power_sampler);

    cuda_check(cudaEventRecord(start, stream), "cudaEventRecord(start)");
    for (int repeat = 0; repeat < graph_repeats; ++repeat) {
        cuda_check(cudaGraphLaunch(graph_exec, stream), "cudaGraphLaunch(timed)");
    }
    cuda_check(cudaEventRecord(stop, stream), "cudaEventRecord(stop)");
    cuda_check(cudaEventSynchronize(stop), "cudaEventSynchronize(stop)");

    float elapsed_ms = 0.0f;
    cuda_check(cudaEventElapsedTime(&elapsed_ms, start, stop), "cudaEventElapsedTime");
    const std::uintmax_t power_sample_end = power_sampler_offset(power_sampler);
    const PowerMeasurement power =
        read_power_samples_between(power_sampler, power_sample_start, power_sample_end);

    T sample{};
    cuda_check(cudaMemcpy(&sample,
                          device_output + (args.elements / 2),
                          sizeof(T),
                          cudaMemcpyDeviceToHost),
               "cudaMemcpy(device->host sample)");

    cuda_check(cudaEventDestroy(start), "cudaEventDestroy(start)");
    cuda_check(cudaEventDestroy(stop), "cudaEventDestroy(stop)");
    cuda_check(cudaGraphExecDestroy(graph_exec), "cudaGraphExecDestroy");
    cuda_check(cudaGraphDestroy(graph), "cudaGraphDestroy");
    cuda_check(cudaStreamDestroy(stream), "cudaStreamDestroy");
    cuda_check(cudaFree(device_input), "cudaFree(device_input)");
    cuda_check(cudaFree(device_output), "cudaFree(device_output)");

    const double elapsed_s = static_cast<double>(elapsed_ms) / 1000.0;
    const double total_elements = static_cast<double>(args.elements) * static_cast<double>(total_launches);
    const double ns_per_launch = elapsed_s > 0.0 ? (elapsed_s * 1.0e9) / static_cast<double>(total_launches) : 0.0;
    const double ns_per_element = total_elements > 0.0 ? (elapsed_s * 1.0e9) / total_elements : 0.0;
    const double elements_per_second = elapsed_s > 0.0 ? total_elements / elapsed_s : 0.0;

    return {
        spec.op,
        spec.mode,
        precision,
        run,
        args.elements,
        total_launches,
        elapsed_s,
        ns_per_launch,
        ns_per_element,
        elements_per_second,
        ScalarCodec<T>::to_double(sample),
        power.average_power_watts,
        power.power_source,
    };
}

template <typename T>
void run_precision_suite(const Args& args,
                         const std::string& precision,
                         std::vector<Record>* records,
                         const PowerSamplerProcess& power_sampler) {
    if (!args.precisions.empty() && !contains_value(args.precisions, precision)) {
        return;
    }

    const auto specs = specs_for_precision<T>();
    std::vector<const BenchSpec<T>*> selected_specs;
    for (const auto& spec : specs) {
        if (!args.ops.empty() && !contains_value(args.ops, spec.op)) {
            continue;
        }
        if (!args.modes.empty() && !contains_value(args.modes, spec.mode)) {
            continue;
        }
        selected_specs.push_back(&spec);
    }

    for (const auto* selected_spec : selected_specs) {
        const auto& spec = *selected_spec;
        for (int run = 0; run < args.runs; ++run) {
            const auto record = benchmark_one(spec, precision, args, run, power_sampler);
            std::cout << "run=" << std::setw(2) << run
                      << " op=" << std::setw(6) << spec.op
                      << " mode=" << std::setw(4) << spec.mode
                      << " precision=" << std::setw(10) << precision
                      << " launches=" << record.launches
                      << " elements=" << args.elements
                      << " ns/launch=" << std::fixed << std::setprecision(3) << record.ns_per_launch
                      << " ns/element=" << record.ns_per_element
                      << " sample=" << record.sample_output;
            if (std::isfinite(record.average_power_watts)) {
                std::cout << " power(W)=" << record.average_power_watts;
            }
            std::cout
                      << '\n';
            records->push_back(record);
        }
    }
}

void write_csv(const std::string& path, const std::vector<Record>& records) {
    const std::filesystem::path output(path);
    if (output.has_parent_path()) {
        std::filesystem::create_directories(output.parent_path());
    }

    std::ofstream out(path);
    if (!out) {
        fail(std::string("could not open CSV for writing: ") + path);
    }

    out << "op,mode,precision,run,elements,launches,time_elapsed_s,ns_per_launch,ns_per_element,elements_per_second,sample_output,average_power_watts,power_source\n";
    out << std::setprecision(17);
    for (const auto& record : records) {
        out << record.op << ','
            << record.mode << ','
            << record.precision << ','
            << record.run << ','
            << record.elements << ','
            << record.launches << ','
            << record.time_elapsed_s << ','
            << record.ns_per_launch << ','
            << record.ns_per_element << ','
            << record.elements_per_second << ','
            << record.sample_output << ',';
        if (std::isfinite(record.average_power_watts)) {
            out << record.average_power_watts;
        }
        out << ',' << record.power_source << '\n';
    }
}

double mean(const std::vector<double>& values) {
    if (values.empty()) {
        return 0.0;
    }
    double sum = 0.0;
    for (double value : values) {
        sum += value;
    }
    return sum / static_cast<double>(values.size());
}

double stdev(const std::vector<double>& values, double avg) {
    if (values.size() < 2) {
        return 0.0;
    }
    double sum_sq = 0.0;
    for (double value : values) {
        const double delta = value - avg;
        sum_sq += delta * delta;
    }
    return std::sqrt(sum_sq / static_cast<double>(values.size() - 1));
}

void print_summary(const std::vector<Record>& records) {
    using Key = std::tuple<std::string, std::string, std::string>;
    std::map<Key, std::vector<const Record*>> grouped;
    for (const auto& record : records) {
        grouped[{record.op, record.mode, record.precision}].push_back(&record);
    }

    std::cout << "\nSummary (mean +/- stdev ns/element):\n";
    for (const auto& entry : grouped) {
        std::vector<double> ns_per_element;
        std::vector<double> throughput;
        for (const auto* record : entry.second) {
            ns_per_element.push_back(record->ns_per_element);
            throughput.push_back(record->elements_per_second);
        }

        const double avg_ns = mean(ns_per_element);
        const double sd_ns = stdev(ns_per_element, avg_ns);
        const double avg_eps = mean(throughput);
        std::cout << "  op=" << std::setw(6) << std::get<0>(entry.first)
              << " mode=" << std::setw(17) << std::get<1>(entry.first)
              << " precision=" << std::setw(10) << std::get<2>(entry.first)
                  << " ns/element=" << std::fixed << std::setprecision(3) << avg_ns
                  << " +/- " << sd_ns
                  << " elems/s=" << std::setprecision(2) << avg_eps
                  << '\n';
    }
}

}  // namespace

int main(int argc, char** argv) {
    const Args args = parse_args(argc, argv);
    validate_filters(args);
    const PowerSamplerProcess power_sampler = start_power_sampler(args);

    if (args.power_sample_ms > 0) {
        if (power_sampling_enabled(args)) {
            std::cout << "Per-row power sampling: nvidia-smi enabled (target="
                      << args.power_sample_ms
                      << "ms, poll="
                      << args.power_poll_ms
                      << "ms, gpu="
                      << gpu_selector_for_nvidia_smi()
                      << ")\n";
        } else {
            std::cout << "Per-row power sampling: nvidia-smi unavailable; power columns will remain empty\n";
        }
    }

    cuda_check(cudaSetDevice(0), "cudaSetDevice(0)");

    std::vector<Record> records;
    run_precision_suite<__half>(args, "f16", &records, power_sampler);
    run_precision_suite<float>(args, "f32", &records, power_sampler);
    run_precision_suite<double>(args, "f64", &records, power_sampler);
    run_precision_suite<__nv_fp4_e2m1>(args, "fp4_e2m1", &records, power_sampler);
    run_precision_suite<__nv_fp6_e3m2>(args, "fp6_e3m2", &records, power_sampler);
    run_precision_suite<__nv_fp8_e4m3>(args, "fp8_e4m3", &records, power_sampler);
    run_precision_suite<__nv_fp8_e5m2>(args, "fp8_e5m2", &records, power_sampler);
    stop_power_sampler(power_sampler);

    if (records.empty()) {
        fail("no benchmark combinations matched the requested filters");
    }

    write_csv(args.csv_path, records);
    std::cout << "\nWrote " << args.csv_path << '\n';
    print_summary(records);

    cuda_check(cudaDeviceReset(), "cudaDeviceReset()");
    return 0;
}
