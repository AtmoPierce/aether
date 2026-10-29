#include <cuda_runtime.h>
#include <cuda_fp16.h>
#include <cuda_fp4.h>
#include <cuda_fp6.h>
#include <cuda_fp8.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
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
    std::string csv_path = "target/cuda_precision_summary.csv";
    std::size_t samples = 1u << 20;
    int block_size = 256;
};

struct PrecisionRecord {
    std::string op;
    std::string mode;
    std::string precision;
    std::string test;
    std::size_t samples;
    double max_abs;
    std::string max_ulp;
    double max_abs_x;
    double max_abs_ref;
    double max_abs_out;
    std::string max_ulp_x;
    std::string max_ulp_ref;
    std::string max_ulp_out;
};

[[noreturn]] void fail(const std::string& message) {
    std::cerr << "cuda_precision_main: " << message << '\n';
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
        const auto first = item.find_first_not_of(" \t\n\r");
        if (first == std::string::npos) {
            continue;
        }
        const auto last = item.find_last_not_of(" \t\n\r");
        out.push_back(item.substr(first, last - first + 1));
    }
    return out;
}

bool contains_value(const std::vector<std::string>& values, const std::string& needle) {
    return std::find(values.begin(), values.end(), needle) != values.end();
}

std::string format_double(double value) {
    std::ostringstream oss;
    oss << std::setprecision(17) << value;
    return oss.str();
}

void print_usage(const char* argv0) {
    std::cout
        << "usage: " << argv0 << " [options]\n\n"
        << "Options:\n"
        << "  --only <op1,op2,...>          Limit precision sweep to specific ops\n"
        << "  --modes <m1,m2,...>           Limit precision sweep to specific modes\n"
        << "  --precisions <p1,p2,...>      Limit precision sweep to specific precisions\n"
        << "  --samples <n>                 Samples per op/mode/precision over the configured domain\n"
        << "                                (default: 1048576)\n"
        << "  --block-size <n>              CUDA block size (default: 256)\n"
        << "  --csv <path>                  Output summary CSV (default: target/cuda_precision_summary.csv)\n"
        << "  --help                        Show this help message\n\n"
        << "Reduced-domain modes (all ops): fpx_kernel_horner,nat_kernel\n"
        << "Full-kernel modes (log only): fpx,nat\n"
        << "Supported precisions: f16,f32,f64,fp4_e2m1,fp6_e3m2,fp8_e4m3,fp8_e5m2\n";
}

Args parse_args(int argc, char** argv) {
    Args args;
    for (int index = 1; index < argc; ++index) {
        const std::string flag = argv[index];
        const auto require_value = [&](const char* name) -> std::string {
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
        } else if (flag == "--samples") {
            args.samples = static_cast<std::size_t>(std::stoull(require_value("--samples")));
        } else if (flag == "--block-size") {
            args.block_size = std::stoi(require_value("--block-size"));
        } else if (flag == "--csv") {
            args.csv_path = require_value("--csv");
        } else if (flag == "--help" || flag == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        } else {
            fail(std::string("unknown argument: ") + flag);
        }
    }

    if (args.samples == 0) {
        fail("--samples must be positive");
    }
    if (args.block_size <= 0) {
        fail("--block-size must be positive");
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
using LaunchFn = void (*)(const T*, T*, std::size_t, dim3, dim3);

template <typename T>
using RefFn = T (*)(T);

template <typename T>
struct BenchSpec {
    const char* op;
    const char* mode;
    Range range;
    LaunchFn<T> launch;
    RefFn<T> reference;
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

#define DEFINE_NATIVE_DEVICE_FUNCS(op, float_fn, double_fn)                                      \
    __device__ __forceinline__ float native_##op(float value) { return float_fn(value); }       \
    __device__ __forceinline__ double native_##op(double value) { return double_fn(value); }    \
    __device__ __forceinline__ __half native_##op(__half value) {                               \
        return ScalarCodec<__half>::from_float(float_fn(static_cast<float>(value)));            \
    }                                                                                            \
    __device__ __forceinline__ __nv_fp4_e2m1 native_##op(__nv_fp4_e2m1 value) {                 \
        return ScalarCodec<__nv_fp4_e2m1>::from_float(float_fn(static_cast<float>(value)));     \
    }                                                                                            \
    __device__ __forceinline__ __nv_fp6_e3m2 native_##op(__nv_fp6_e3m2 value) {                 \
        return ScalarCodec<__nv_fp6_e3m2>::from_float(float_fn(static_cast<float>(value)));     \
    }                                                                                            \
    __device__ __forceinline__ __nv_fp8_e4m3 native_##op(__nv_fp8_e4m3 value) {                 \
        return ScalarCodec<__nv_fp8_e4m3>::from_float(float_fn(static_cast<float>(value)));     \
    }                                                                                            \
    __device__ __forceinline__ __nv_fp8_e5m2 native_##op(__nv_fp8_e5m2 value) {                 \
        return ScalarCodec<__nv_fp8_e5m2>::from_float(float_fn(static_cast<float>(value)));     \
    }

DEFINE_NATIVE_DEVICE_FUNCS(acos, ::acosf, ::acos)
DEFINE_NATIVE_DEVICE_FUNCS(asin, ::asinf, ::asin)
DEFINE_NATIVE_DEVICE_FUNCS(atan, ::atanf, ::atan)
DEFINE_NATIVE_DEVICE_FUNCS(cos, ::cosf, ::cos)
DEFINE_NATIVE_DEVICE_FUNCS(cosh, ::coshf, ::cosh)
DEFINE_NATIVE_DEVICE_FUNCS(exp, ::expf, ::exp)
DEFINE_NATIVE_DEVICE_FUNCS(expm1, ::expm1f, ::expm1)
DEFINE_NATIVE_DEVICE_FUNCS(log, ::logf, ::log)
DEFINE_NATIVE_DEVICE_FUNCS(log1p, ::log1pf, ::log1p)
DEFINE_NATIVE_DEVICE_FUNCS(sin, ::sinf, ::sin)
DEFINE_NATIVE_DEVICE_FUNCS(sinh, ::sinhf, ::sinh)
DEFINE_NATIVE_DEVICE_FUNCS(tan, ::tanf, ::tan)
DEFINE_NATIVE_DEVICE_FUNCS(tanh, ::tanhf, ::tanh)

#define DEFINE_REDUCED_LAUNCHERS(op)                                                                   \
    static void launch_##op##_f64(const double* xs, double* ys, std::size_t n, dim3 grid, dim3 blk) { \
        op##_reduced_kernel_double<<<grid, blk>>>(xs, ys, n);                                        \
    }                                                                                                  \
    static void launch_##op##_f32(const float* xs, float* ys, std::size_t n, dim3 grid, dim3 blk) {   \
        op##_reduced_kernel_float<<<grid, blk>>>(xs, ys, n);                                         \
    }                                                                                                  \
    static void launch_##op##_f16(const __half* xs, __half* ys, std::size_t n, dim3 grid, dim3 blk) { \
        op##_reduced_kernel_f16<<<grid, blk>>>(xs, ys, n);                                           \
    }                                                                                                  \
    static void launch_##op##_fp4(const __nv_fp4_e2m1* xs,                                            \
                                  __nv_fp4_e2m1* ys,                                                  \
                                  std::size_t n,                                                      \
                                  dim3 grid,                                                          \
                                  dim3 blk) {                                                         \
        op##_reduced_kernel_fp4_e2m1<<<grid, blk>>>(xs, ys, n);                                      \
    }                                                                                                  \
    static void launch_##op##_fp6(const __nv_fp6_e3m2* xs,                                            \
                                  __nv_fp6_e3m2* ys,                                                  \
                                  std::size_t n,                                                      \
                                  dim3 grid,                                                          \
                                  dim3 blk) {                                                         \
        op##_reduced_kernel_fp6_e3m2<<<grid, blk>>>(xs, ys, n);                                      \
    }                                                                                                  \
    static void launch_##op##_fp8e4(const __nv_fp8_e4m3* xs,                                          \
                                    __nv_fp8_e4m3* ys,                                                \
                                    std::size_t n,                                                    \
                                    dim3 grid,                                                        \
                                    dim3 blk) {                                                       \
        op##_reduced_kernel_fp8_e4m3<<<grid, blk>>>(xs, ys, n);                                      \
    }                                                                                                  \
    static void launch_##op##_fp8e5(const __nv_fp8_e5m2* xs,                                          \
                                    __nv_fp8_e5m2* ys,                                                \
                                    std::size_t n,                                                    \
                                    dim3 grid,                                                        \
                                    dim3 blk) {                                                       \
        op##_reduced_kernel_fp8_e5m2<<<grid, blk>>>(xs, ys, n);                                      \
    }

FPX_CUDA_OPS(DEFINE_REDUCED_LAUNCHERS)

#define DEFINE_FULL_LAUNCHERS(op)                                                                   \
    static void launch_full_##op##_f64(                                                            \
        const double* xs, double* ys, std::size_t n, dim3 grid, dim3 blk) {                         \
        op##_kernel_double<<<grid, blk>>>(xs, ys, n);                                              \
    }                                                                                               \
    static void launch_full_##op##_f32(                                                            \
        const float* xs, float* ys, std::size_t n, dim3 grid, dim3 blk) {                           \
        op##_kernel_float<<<grid, blk>>>(xs, ys, n);                                               \
    }                                                                                               \
    static void launch_full_##op##_f16(                                                            \
        const __half* xs, __half* ys, std::size_t n, dim3 grid, dim3 blk) {                         \
        op##_kernel_f16<<<grid, blk>>>(xs, ys, n);                                                 \
    }                                                                                               \
    static void launch_full_##op##_fp4(const __nv_fp4_e2m1* xs,                                   \
                                       __nv_fp4_e2m1* ys,                                         \
                                       std::size_t n,                                              \
                                       dim3 grid,                                                  \
                                       dim3 blk) {                                                 \
        op##_kernel_fp4_e2m1<<<grid, blk>>>(xs, ys, n);                                           \
    }                                                                                               \
    static void launch_full_##op##_fp6(const __nv_fp6_e3m2* xs,                                   \
                                       __nv_fp6_e3m2* ys,                                         \
                                       std::size_t n,                                              \
                                       dim3 grid,                                                  \
                                       dim3 blk) {                                                 \
        op##_kernel_fp6_e3m2<<<grid, blk>>>(xs, ys, n);                                           \
    }                                                                                               \
    static void launch_full_##op##_fp8e4(const __nv_fp8_e4m3* xs,                                 \
                                         __nv_fp8_e4m3* ys,                                       \
                                         std::size_t n,                                            \
                                         dim3 grid,                                                \
                                         dim3 blk) {                                               \
        op##_kernel_fp8_e4m3<<<grid, blk>>>(xs, ys, n);                                           \
    }                                                                                               \
    static void launch_full_##op##_fp8e5(const __nv_fp8_e5m2* xs,                                 \
                                         __nv_fp8_e5m2* ys,                                       \
                                         std::size_t n,                                            \
                                         dim3 grid,                                                \
                                         dim3 blk) {                                               \
        op##_kernel_fp8_e5m2<<<grid, blk>>>(xs, ys, n);                                           \
    }

DEFINE_FULL_LAUNCHERS(log)

#define DEFINE_NATIVE_LAUNCHERS(op)                                                                     \
    __global__ void native_##op##_kernel_double(const double* xs, double* ys, std::size_t n) {       \
        const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;    \
        if (index < n) {                                                                                \
            ys[index] = native_##op(xs[index]);                                                         \
        }                                                                                               \
    }                                                                                                   \
    __global__ void native_##op##_kernel_float(const float* xs, float* ys, std::size_t n) {           \
        const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;    \
        if (index < n) {                                                                                \
            ys[index] = native_##op(xs[index]);                                                         \
        }                                                                                               \
    }                                                                                                   \
    __global__ void native_##op##_kernel_f16(const __half* xs, __half* ys, std::size_t n) {           \
        const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;    \
        if (index < n) {                                                                                \
            ys[index] = native_##op(xs[index]);                                                         \
        }                                                                                               \
    }                                                                                                   \
    __global__ void native_##op##_kernel_fp4_e2m1(const __nv_fp4_e2m1* xs,                             \
                                                   __nv_fp4_e2m1* ys,                                   \
                                                   std::size_t n) {                                     \
        const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;    \
        if (index < n) {                                                                                \
            ys[index] = native_##op(xs[index]);                                                         \
        }                                                                                               \
    }                                                                                                   \
    __global__ void native_##op##_kernel_fp6_e3m2(const __nv_fp6_e3m2* xs,                             \
                                                   __nv_fp6_e3m2* ys,                                   \
                                                   std::size_t n) {                                     \
        const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;    \
        if (index < n) {                                                                                \
            ys[index] = native_##op(xs[index]);                                                         \
        }                                                                                               \
    }                                                                                                   \
    __global__ void native_##op##_kernel_fp8_e4m3(const __nv_fp8_e4m3* xs,                             \
                                                   __nv_fp8_e4m3* ys,                                   \
                                                   std::size_t n) {                                     \
        const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;    \
        if (index < n) {                                                                                \
            ys[index] = native_##op(xs[index]);                                                         \
        }                                                                                               \
    }                                                                                                   \
    __global__ void native_##op##_kernel_fp8_e5m2(const __nv_fp8_e5m2* xs,                             \
                                                   __nv_fp8_e5m2* ys,                                   \
                                                   std::size_t n) {                                     \
        const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;    \
        if (index < n) {                                                                                \
            ys[index] = native_##op(xs[index]);                                                         \
        }                                                                                               \
    }                                                                                                   \
    static void launch_native_##op##_f64(const double* xs, double* ys, std::size_t n, dim3 grid, dim3 blk) { \
        native_##op##_kernel_double<<<grid, blk>>>(xs, ys, n);                                         \
    }                                                                                                   \
    static void launch_native_##op##_f32(const float* xs, float* ys, std::size_t n, dim3 grid, dim3 blk) { \
        native_##op##_kernel_float<<<grid, blk>>>(xs, ys, n);                                          \
    }                                                                                                   \
    static void launch_native_##op##_f16(const __half* xs, __half* ys, std::size_t n, dim3 grid, dim3 blk) { \
        native_##op##_kernel_f16<<<grid, blk>>>(xs, ys, n);                                            \
    }                                                                                                   \
    static void launch_native_##op##_fp4(const __nv_fp4_e2m1* xs,                                      \
                                         __nv_fp4_e2m1* ys,                                            \
                                         std::size_t n,                                                \
                                         dim3 grid,                                                    \
                                         dim3 blk) {                                                   \
        native_##op##_kernel_fp4_e2m1<<<grid, blk>>>(xs, ys, n);                                       \
    }                                                                                                   \
    static void launch_native_##op##_fp6(const __nv_fp6_e3m2* xs,                                      \
                                         __nv_fp6_e3m2* ys,                                            \
                                         std::size_t n,                                                \
                                         dim3 grid,                                                    \
                                         dim3 blk) {                                                   \
        native_##op##_kernel_fp6_e3m2<<<grid, blk>>>(xs, ys, n);                                       \
    }                                                                                                   \
    static void launch_native_##op##_fp8e4(const __nv_fp8_e4m3* xs,                                    \
                                           __nv_fp8_e4m3* ys,                                          \
                                           std::size_t n,                                              \
                                           dim3 grid,                                                  \
                                           dim3 blk) {                                                 \
        native_##op##_kernel_fp8_e4m3<<<grid, blk>>>(xs, ys, n);                                       \
    }                                                                                                   \
    static void launch_native_##op##_fp8e5(const __nv_fp8_e5m2* xs,                                    \
                                           __nv_fp8_e5m2* ys,                                          \
                                           std::size_t n,                                              \
                                           dim3 grid,                                                  \
                                           dim3 blk) {                                                 \
        native_##op##_kernel_fp8_e5m2<<<grid, blk>>>(xs, ys, n);                                       \
    }

FPX_CUDA_OPS(DEFINE_NATIVE_LAUNCHERS)

#define DEFINE_REF_FUNCS(op, float_fn, double_fn)                                                        \
    static float ref_##op##_f32(float value) { return float_fn(value); }                                \
    static double ref_##op##_f64(double value) { return double_fn(value); }                              \
    static __half ref_##op##_f16(__half value) {                                                          \
        return ScalarCodec<__half>::from_float(float_fn(static_cast<float>(value)));                    \
    }                                                                                                     \
    static __nv_fp4_e2m1 ref_##op##_fp4(__nv_fp4_e2m1 value) {                                            \
        return ScalarCodec<__nv_fp4_e2m1>::from_float(float_fn(static_cast<float>(value)));              \
    }                                                                                                     \
    static __nv_fp6_e3m2 ref_##op##_fp6(__nv_fp6_e3m2 value) {                                            \
        return ScalarCodec<__nv_fp6_e3m2>::from_float(float_fn(static_cast<float>(value)));              \
    }                                                                                                     \
    static __nv_fp8_e4m3 ref_##op##_fp8e4(__nv_fp8_e4m3 value) {                                          \
        return ScalarCodec<__nv_fp8_e4m3>::from_float(float_fn(static_cast<float>(value)));              \
    }                                                                                                     \
    static __nv_fp8_e5m2 ref_##op##_fp8e5(__nv_fp8_e5m2 value) {                                          \
        return ScalarCodec<__nv_fp8_e5m2>::from_float(float_fn(static_cast<float>(value)));              \
    }

DEFINE_REF_FUNCS(acos, ::acosf, ::acos)
DEFINE_REF_FUNCS(asin, ::asinf, ::asin)
DEFINE_REF_FUNCS(atan, ::atanf, ::atan)
DEFINE_REF_FUNCS(cos, ::cosf, ::cos)
DEFINE_REF_FUNCS(cosh, ::coshf, ::cosh)
DEFINE_REF_FUNCS(exp, ::expf, ::exp)
DEFINE_REF_FUNCS(expm1, ::expm1f, ::expm1)
DEFINE_REF_FUNCS(log, ::logf, ::log)
DEFINE_REF_FUNCS(log1p, ::log1pf, ::log1p)
DEFINE_REF_FUNCS(sin, ::sinf, ::sin)
DEFINE_REF_FUNCS(sinh, ::sinhf, ::sinh)
DEFINE_REF_FUNCS(tan, ::tanf, ::tan)
DEFINE_REF_FUNCS(tanh, ::tanhf, ::tanh)

template <typename T>
std::vector<BenchSpec<T>> specs_for_precision();

#define SPEC_F64(op)                                                                                                 \
    BenchSpec<double>{#op, "fpx_kernel_horner", reduced_domain_for_op(#op), launch_##op##_f64, ref_##op##_f64},   \
        BenchSpec<double>{                                                                                           \
            #op, "nat_kernel", reduced_domain_for_op(#op), launch_native_##op##_f64, ref_##op##_f64},
template <>
std::vector<BenchSpec<double>> specs_for_precision<double>() {
    return {
        FPX_CUDA_OPS(SPEC_F64)
        BenchSpec<double>{
            "log", "fpx", full_domain_for_op("log"), launch_full_log_f64, ref_log_f64},
        BenchSpec<double>{
            "log", "nat", full_domain_for_op("log"), launch_native_log_f64, ref_log_f64},
    };
}

#define SPEC_F32(op)                                                                                               \
    BenchSpec<float>{#op, "fpx_kernel_horner", reduced_domain_for_op(#op), launch_##op##_f32, ref_##op##_f32},  \
        BenchSpec<float>{                                                                                         \
            #op, "nat_kernel", reduced_domain_for_op(#op), launch_native_##op##_f32, ref_##op##_f32},
template <>
std::vector<BenchSpec<float>> specs_for_precision<float>() {
    return {
        FPX_CUDA_OPS(SPEC_F32)
        BenchSpec<float>{
            "log", "fpx", full_domain_for_op("log"), launch_full_log_f32, ref_log_f32},
        BenchSpec<float>{
            "log", "nat", full_domain_for_op("log"), launch_native_log_f32, ref_log_f32},
    };
}

#define SPEC_F16(op)                                                                                                 \
    BenchSpec<__half>{#op, "fpx_kernel_horner", reduced_domain_for_op(#op), launch_##op##_f16, ref_##op##_f16},  \
        BenchSpec<__half>{                                                                                           \
            #op, "nat_kernel", reduced_domain_for_op(#op), launch_native_##op##_f16, ref_##op##_f16},
template <>
std::vector<BenchSpec<__half>> specs_for_precision<__half>() {
    return {
        FPX_CUDA_OPS(SPEC_F16)
        BenchSpec<__half>{
            "log", "fpx", full_domain_for_op("log"), launch_full_log_f16, ref_log_f16},
        BenchSpec<__half>{
            "log", "nat", full_domain_for_op("log"), launch_native_log_f16, ref_log_f16},
    };
}

#define SPEC_FP4(op)                                                                                                       \
    BenchSpec<__nv_fp4_e2m1>{                                                                                             \
        #op, "fpx_kernel_horner", reduced_domain_for_op(#op), launch_##op##_fp4, ref_##op##_fp4},                     \
        BenchSpec<__nv_fp4_e2m1>{                                                                                         \
            #op, "nat_kernel", reduced_domain_for_op(#op), launch_native_##op##_fp4, ref_##op##_fp4},
template <>
std::vector<BenchSpec<__nv_fp4_e2m1>> specs_for_precision<__nv_fp4_e2m1>() {
    return {
        FPX_CUDA_OPS(SPEC_FP4)
        BenchSpec<__nv_fp4_e2m1>{
            "log", "fpx", full_domain_for_op("log"), launch_full_log_fp4, ref_log_fp4},
        BenchSpec<__nv_fp4_e2m1>{
            "log", "nat", full_domain_for_op("log"), launch_native_log_fp4, ref_log_fp4},
    };
}

#define SPEC_FP6(op)                                                                                                       \
    BenchSpec<__nv_fp6_e3m2>{                                                                                             \
        #op, "fpx_kernel_horner", reduced_domain_for_op(#op), launch_##op##_fp6, ref_##op##_fp6},                     \
        BenchSpec<__nv_fp6_e3m2>{                                                                                         \
            #op, "nat_kernel", reduced_domain_for_op(#op), launch_native_##op##_fp6, ref_##op##_fp6},
template <>
std::vector<BenchSpec<__nv_fp6_e3m2>> specs_for_precision<__nv_fp6_e3m2>() {
    return {
        FPX_CUDA_OPS(SPEC_FP6)
        BenchSpec<__nv_fp6_e3m2>{
            "log", "fpx", full_domain_for_op("log"), launch_full_log_fp6, ref_log_fp6},
        BenchSpec<__nv_fp6_e3m2>{
            "log", "nat", full_domain_for_op("log"), launch_native_log_fp6, ref_log_fp6},
    };
}

#define SPEC_FP8E4(op)                                                                                                       \
    BenchSpec<__nv_fp8_e4m3>{                                                                                               \
        #op, "fpx_kernel_horner", reduced_domain_for_op(#op), launch_##op##_fp8e4, ref_##op##_fp8e4},                   \
        BenchSpec<__nv_fp8_e4m3>{                                                                                           \
            #op, "nat_kernel", reduced_domain_for_op(#op), launch_native_##op##_fp8e4, ref_##op##_fp8e4},
template <>
std::vector<BenchSpec<__nv_fp8_e4m3>> specs_for_precision<__nv_fp8_e4m3>() {
    return {
        FPX_CUDA_OPS(SPEC_FP8E4)
        BenchSpec<__nv_fp8_e4m3>{
            "log", "fpx", full_domain_for_op("log"), launch_full_log_fp8e4, ref_log_fp8e4},
        BenchSpec<__nv_fp8_e4m3>{
            "log", "nat", full_domain_for_op("log"), launch_native_log_fp8e4, ref_log_fp8e4},
    };
}

#define SPEC_FP8E5(op)                                                                                                       \
    BenchSpec<__nv_fp8_e5m2>{                                                                                               \
        #op, "fpx_kernel_horner", reduced_domain_for_op(#op), launch_##op##_fp8e5, ref_##op##_fp8e5},                   \
        BenchSpec<__nv_fp8_e5m2>{                                                                                           \
            #op, "nat_kernel", reduced_domain_for_op(#op), launch_native_##op##_fp8e5, ref_##op##_fp8e5},
template <>
std::vector<BenchSpec<__nv_fp8_e5m2>> specs_for_precision<__nv_fp8_e5m2>() {
    return {
        FPX_CUDA_OPS(SPEC_FP8E5)
        BenchSpec<__nv_fp8_e5m2>{
            "log", "fpx", full_domain_for_op("log"), launch_full_log_fp8e5, ref_log_fp8e5},
        BenchSpec<__nv_fp8_e5m2>{
            "log", "nat", full_domain_for_op("log"), launch_native_log_fp8e5, ref_log_fp8e5},
    };
}

std::vector<std::string> supported_ops() {
#define OP_NAME(op) #op,
    return {FPX_CUDA_OPS(OP_NAME)};
#undef OP_NAME
}

std::vector<std::string> supported_modes() {
    return {"fpx_kernel_horner", "nat_kernel", "fpx", "nat"};
}

std::vector<std::string> supported_precisions() {
    return {"f16", "f32", "f64", "fp4_e2m1", "fp6_e3m2", "fp8_e4m3", "fp8_e5m2"};
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

template <typename U>
U abs_diff_unsigned(U lhs, U rhs) {
    return lhs > rhs ? static_cast<U>(lhs - rhs) : static_cast<U>(rhs - lhs);
}

template <typename T>
bool is_finite_value(T value) {
    return std::isfinite(ScalarCodec<T>::to_double(value));
}

inline uint64_t ordered_u64_from_f64(double value) {
    union {
        double f64;
        uint64_t u64;
    } bits{};
    bits.f64 = value;
    if ((bits.u64 >> 63) != 0) {
        return ~bits.u64;
    }
    return bits.u64 | (uint64_t{1} << 63);
}

inline uint32_t ordered_u32_from_f32(float value) {
    union {
        float f32;
        uint32_t u32;
    } bits{};
    bits.f32 = value;
    if ((bits.u32 >> 31) != 0) {
        return ~bits.u32;
    }
    return bits.u32 | (uint32_t{1} << 31);
}

inline uint16_t half_bits(__half value) {
    const __half_raw raw = value;
    return raw.x;
}

inline __half half_from_float(float value) {
    return ScalarCodec<__half>::from_float(value);
}

inline uint16_t ordered_u16_from_f16(__half value) {
    const uint16_t bits = half_bits(value);
    if ((bits >> 15) != 0) {
        return static_cast<uint16_t>(~bits);
    }
    return static_cast<uint16_t>(bits | (uint16_t{1} << 15));
}

bool valid_ulp_distance(double lhs, double rhs, uint64_t* out) {
    if (!std::isfinite(lhs) || !std::isfinite(rhs)) {
        return false;
    }

    if (lhs == 0.0) {
        lhs = 0.0;
    }
    if (rhs == 0.0) {
        rhs = 0.0;
    }

    if (std::signbit(lhs) != std::signbit(rhs)) {
        const uint64_t zero = ordered_u64_from_f64(0.0);
        const uint64_t dl = abs_diff_unsigned(ordered_u64_from_f64(std::abs(lhs)), zero);
        const uint64_t dr = abs_diff_unsigned(ordered_u64_from_f64(std::abs(rhs)), zero);
        *out = dl + dr;
    } else {
        *out = abs_diff_unsigned(ordered_u64_from_f64(lhs), ordered_u64_from_f64(rhs));
    }
    return true;
}

bool valid_ulp_distance(float lhs, float rhs, uint64_t* out) {
    if (!std::isfinite(lhs) || !std::isfinite(rhs)) {
        return false;
    }

    if (lhs == 0.0f) {
        lhs = 0.0f;
    }
    if (rhs == 0.0f) {
        rhs = 0.0f;
    }

    if (std::signbit(lhs) != std::signbit(rhs)) {
        const uint32_t zero = ordered_u32_from_f32(0.0f);
        const uint32_t dl = abs_diff_unsigned(ordered_u32_from_f32(std::fabs(lhs)), zero);
        const uint32_t dr = abs_diff_unsigned(ordered_u32_from_f32(std::fabs(rhs)), zero);
        *out = static_cast<uint64_t>(dl) + static_cast<uint64_t>(dr);
    } else {
        *out = abs_diff_unsigned(ordered_u32_from_f32(lhs), ordered_u32_from_f32(rhs));
    }
    return true;
}

bool valid_ulp_distance(__half lhs, __half rhs, uint64_t* out) {
    if (!is_finite_value(lhs) || !is_finite_value(rhs)) {
        return false;
    }

    if (ScalarCodec<__half>::to_double(lhs) == 0.0) {
        lhs = half_from_float(0.0f);
    }
    if (ScalarCodec<__half>::to_double(rhs) == 0.0) {
        rhs = half_from_float(0.0f);
    }

    const bool lhs_negative = (half_bits(lhs) >> 15) != 0;
    const bool rhs_negative = (half_bits(rhs) >> 15) != 0;

    if (lhs_negative != rhs_negative) {
        const uint16_t zero = ordered_u16_from_f16(half_from_float(0.0f));
        const uint16_t dl = abs_diff_unsigned(ordered_u16_from_f16(half_from_float(std::fabs(static_cast<float>(lhs)))), zero);
        const uint16_t dr = abs_diff_unsigned(ordered_u16_from_f16(half_from_float(std::fabs(static_cast<float>(rhs)))), zero);
        *out = static_cast<uint64_t>(dl) + static_cast<uint64_t>(dr);
    } else {
        *out = abs_diff_unsigned(ordered_u16_from_f16(lhs), ordered_u16_from_f16(rhs));
    }
    return true;
}

template <typename T>
bool valid_ulp_distance(T, T, uint64_t*) {
    return false;
}

template <typename T>
PrecisionRecord benchmark_one(const BenchSpec<T>& spec, const std::string& precision, const Args& args) {
    std::vector<T> host_input(args.samples);
    for (std::size_t index = 0; index < args.samples; ++index) {
        const double alpha = args.samples == 1
            ? 0.0
            : static_cast<double>(index) / static_cast<double>(args.samples - 1);
        const double value = spec.range.lo + (spec.range.hi - spec.range.lo) * alpha;
        host_input[index] = ScalarCodec<T>::from_double(value);
    }

    T* device_input = nullptr;
    T* device_output = nullptr;
    cuda_check(cudaMalloc(&device_input, sizeof(T) * args.samples), "cudaMalloc(device_input)");
    cuda_check(cudaMalloc(&device_output, sizeof(T) * args.samples), "cudaMalloc(device_output)");
    cuda_check(
        cudaMemcpy(device_input, host_input.data(), sizeof(T) * args.samples, cudaMemcpyHostToDevice),
        "cudaMemcpy(host->device)"
    );

    const dim3 block(static_cast<unsigned int>(args.block_size));
    const dim3 grid(static_cast<unsigned int>((args.samples + static_cast<std::size_t>(args.block_size) - 1) /
                                              static_cast<std::size_t>(args.block_size)));

    spec.launch(device_input, device_output, args.samples, grid, block);
    cuda_check(cudaGetLastError(), "precision launch");
    cuda_check(cudaDeviceSynchronize(), "precision synchronize");

    std::vector<T> host_output(args.samples);
    cuda_check(
        cudaMemcpy(host_output.data(), device_output, sizeof(T) * args.samples, cudaMemcpyDeviceToHost),
        "cudaMemcpy(device->host)"
    );

    cuda_check(cudaFree(device_input), "cudaFree(device_input)");
    cuda_check(cudaFree(device_output), "cudaFree(device_output)");

    bool saw_finite = false;
    double max_abs = -1.0;
    double max_abs_x = 0.0;
    double max_abs_ref = 0.0;
    double max_abs_out = 0.0;

    bool has_max_ulp = false;
    uint64_t max_ulp = 0;
    double max_ulp_x = 0.0;
    double max_ulp_ref = 0.0;
    double max_ulp_out = 0.0;

    for (std::size_t index = 0; index < args.samples; ++index) {
        const T x = host_input[index];
        const T reference_value = spec.reference(x);
        const T output_value = host_output[index];

        const double x_value = ScalarCodec<T>::to_double(x);
        const double reference_double = ScalarCodec<T>::to_double(reference_value);
        const double output_double = ScalarCodec<T>::to_double(output_value);

        if (std::isfinite(x_value) && std::isfinite(reference_double) && std::isfinite(output_double)) {
            saw_finite = true;
            const double abs_error = std::fabs(reference_double - output_double);
            if (abs_error > max_abs || (abs_error == max_abs && x_value < max_abs_x)) {
                max_abs = abs_error;
                max_abs_x = x_value;
                max_abs_ref = reference_double;
                max_abs_out = output_double;
            }
        }

        uint64_t ulp_distance = 0;
        if (valid_ulp_distance(reference_value, output_value, &ulp_distance)) {
            if (!has_max_ulp || ulp_distance > max_ulp || (ulp_distance == max_ulp && x_value < max_ulp_x)) {
                has_max_ulp = true;
                max_ulp = ulp_distance;
                max_ulp_x = x_value;
                max_ulp_ref = reference_double;
                max_ulp_out = output_double;
            }
        }
    }

    if (!saw_finite) {
        fail(std::string("no finite comparisons for op=") + spec.op + " mode=" + spec.mode + " precision=" + precision);
    }

    return {
        spec.op,
        spec.mode,
        precision,
        "uniform_grid",
        args.samples,
        max_abs,
        has_max_ulp ? std::to_string(max_ulp) : "NA",
        max_abs_x,
        max_abs_ref,
        max_abs_out,
        has_max_ulp ? format_double(max_ulp_x) : "NA",
        has_max_ulp ? format_double(max_ulp_ref) : "NA",
        has_max_ulp ? format_double(max_ulp_out) : "NA",
    };
}

template <typename T>
void run_precision_suite(const Args& args,
                         const std::string& precision,
                         std::vector<PrecisionRecord>* records) {
    if (!args.precisions.empty() && !contains_value(args.precisions, precision)) {
        return;
    }

    const auto specs = specs_for_precision<T>();
    for (const auto& spec : specs) {
        if (!args.ops.empty() && !contains_value(args.ops, spec.op)) {
            continue;
        }
        if (!args.modes.empty() && !contains_value(args.modes, spec.mode)) {
            continue;
        }

        const auto record = benchmark_one(spec, precision, args);
        std::cout << "op=" << std::setw(6) << spec.op
                  << " mode=" << std::setw(17) << spec.mode
                  << " precision=" << std::setw(10) << precision
                  << " samples=" << args.samples
                  << " max_abs=" << std::fixed << std::setprecision(9) << record.max_abs
                  << " max_ulp=" << record.max_ulp
                  << '\n';
        records->push_back(record);
    }
}

void write_csv(const std::string& path, const std::vector<PrecisionRecord>& records) {
    const std::filesystem::path output(path);
    if (output.has_parent_path()) {
        std::filesystem::create_directories(output.parent_path());
    }

    std::ofstream out(path);
    if (!out) {
        fail(std::string("could not open CSV for writing: ") + path);
    }

    out << "op,mode,precision,test,samples,max_abs,max_ulp,max_abs_x,max_abs_ref,max_abs_out,max_ulp_x,max_ulp_ref,max_ulp_out\n";
    out << std::setprecision(17);
    for (const auto& record : records) {
        out << record.op << ','
            << record.mode << ','
            << record.precision << ','
            << record.test << ','
            << record.samples << ','
            << record.max_abs << ','
            << record.max_ulp << ','
            << record.max_abs_x << ','
            << record.max_abs_ref << ','
            << record.max_abs_out << ','
            << record.max_ulp_x << ','
            << record.max_ulp_ref << ','
            << record.max_ulp_out << '\n';
    }
}

}  // namespace

int main(int argc, char** argv) {
    const Args args = parse_args(argc, argv);
    validate_filters(args);

    cuda_check(cudaSetDevice(0), "cudaSetDevice(0)");

    std::vector<PrecisionRecord> records;
    run_precision_suite<__half>(args, "f16", &records);
    run_precision_suite<float>(args, "f32", &records);
    run_precision_suite<double>(args, "f64", &records);
    run_precision_suite<__nv_fp4_e2m1>(args, "fp4_e2m1", &records);
    run_precision_suite<__nv_fp6_e3m2>(args, "fp6_e3m2", &records);
    run_precision_suite<__nv_fp8_e4m3>(args, "fp8_e4m3", &records);
    run_precision_suite<__nv_fp8_e5m2>(args, "fp8_e5m2", &records);

    if (records.empty()) {
        fail("no precision combinations matched the requested filters");
    }

    write_csv(args.csv_path, records);
    std::cout << "\nWrote " << args.csv_path << '\n';

    cuda_check(cudaDeviceReset(), "cudaDeviceReset()");
    return 0;
}
