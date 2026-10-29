#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../../../../.." && pwd)"
X86_TARGET_DIR="${REPO_ROOT}/target/asm-simd-x86_64"
AARCH64_TARGET_DIR="${REPO_ROOT}/target/asm-simd-aarch64"
REPORT_DIR="${REPO_ROOT}/outputs/asm"
REPORT_PATH="${REPORT_DIR}/simd_instruction_forms.csv"
FULL_PATH_DIR="${REPORT_DIR}/fpx_paths"
DOC_PATH_DIR="${REPORT_DIR}/excerpts"

usage() {
    echo "Usage: $0 [all|x86_64|aarch64|analyze|paths]" >&2
    exit 2
}

extract_function() {
    local assembly="$1"
    local function_name="$2"
    local output_path="$3"
    awk -v function_name="${function_name}" '
        $0 == function_name ":" { active = 1 }
        active { print }
        active && /^\.Lfunc_end/ { exit }
    ' "${assembly}" > "${output_path}"
}

extract_x86_excerpt() {
    local assembly="$1"
    local function_name="$2"
    local path="$3"
    local output_path="$4"
    awk -v function_name="${function_name}" -v path="${path}" '
        BEGIN {
            print "# Compiler-produced arithmetic excerpt from " function_name
            print function_name ":"
        }
        $0 == function_name ":" { active = 1; next }
        active && /^\.Lfunc_end/ { exit }
        active && path == "scalar" && $1 ~ /^(vmulss|vaddss)$/ {
            if (arithmetic++ < 6) print
        }
        active && path == "fma" && $1 ~ /^vfmadd.*ss$/ {
            if (arithmetic++ < 4) print
        }
        active && path == "simd" {
            if (!loaded && $1 == "vmovups" && $0 ~ /\([^)]*\), %ymm[0-9]+$/) {
                print
                loaded = 1
            } else if (loaded && fmas < 4 && $1 ~ /^vfmadd.*ps$/) {
                print
                fmas++
            } else if (fmas == 4 && !stored && $1 == "vmovups" && $0 ~ /%ymm[0-9]+, .*\(/) {
                print "\t# remaining stages repeat the packed FMA"
                print
                stored = 1
                exit
            }
        }
        END {
            if (path == "scalar" && arithmetic < 6) exit 11
            if (path == "fma" && arithmetic < 4) exit 12
            if (path == "simd" && (!loaded || fmas < 4 || !stored)) exit 13
        }
    ' "${assembly}" > "${output_path}"
}

extract_aarch64_excerpt() {
    local assembly="$1"
    local function_name="$2"
    local path="$3"
    local output_path="$4"
    awk -v function_name="${function_name}" -v path="${path}" '
        BEGIN {
            print "// Compiler-produced arithmetic excerpt from " function_name
            print function_name ":"
        }
        $0 == function_name ":" { active = 1; next }
        active && /^\.Lfunc_end/ { exit }
        active && path == "scalar" && $1 ~ /^(fmul|fadd)$/ && $2 ~ /^s[0-9]+,/ {
            if (arithmetic++ < 6) print
        }
        active && path == "fma" && $1 == "fmadd" && $2 ~ /^s[0-9]+,/ {
            if (arithmetic++ < 4) print
        }
        active && path == "simd" {
            if (!loaded && $1 == "ldr" && $2 ~ /^q[0-9]+,/) {
                print
                loaded = 1
            } else if (loaded && fmas < 4 && $1 == "fmla" && $2 ~ /^v[0-9]+\.4s,/) {
                print
                fmas++
            } else if (fmas == 4 && !stored && $1 == "str" && $2 ~ /^q[0-9]+,/) {
                print "\t// remaining stages repeat the packed FMA"
                print
                stored = 1
                exit
            }
        }
        END {
            if (path == "scalar" && arithmetic < 6) exit 21
            if (path == "fma" && arithmetic < 4) exit 22
            if (path == "simd" && (!loaded || fmas < 4 || !stored)) exit 23
        }
    ' "${assembly}" > "${output_path}"
}

generate_path_dumps() {
    local x86_assembly
    local aarch64_assembly
    local path function_name
    x86_assembly="$(find_dump_assembly "${X86_TARGET_DIR}")"
    aarch64_assembly="$(find_dump_assembly "${AARCH64_TARGET_DIR}")"
    mkdir -p "${FULL_PATH_DIR}" "${DOC_PATH_DIR}"

    for path in scalar fma simd; do
        case "${path}" in
            scalar) function_name="dump_sin_raw_scalar_f32" ;;
            fma) function_name="dump_sin_raw_fma_f32" ;;
            simd) function_name="dump_sin_raw_simd_horner_f32" ;;
        esac
        extract_function "${x86_assembly}" "${function_name}" \
            "${FULL_PATH_DIR}/x86_64_${path}_f32_full.s"
        extract_function "${aarch64_assembly}" "${function_name}" \
            "${FULL_PATH_DIR}/aarch64_${path}_f32_full.s"
        extract_x86_excerpt "${x86_assembly}" "${function_name}" "${path}" \
            "${DOC_PATH_DIR}/x86_64_${path}_f32.s"
        extract_aarch64_excerpt "${aarch64_assembly}" "${function_name}" "${path}" \
            "${DOC_PATH_DIR}/aarch64_${path}_f32.s"
    done

    echo "Full FPX path functions: ${FULL_PATH_DIR}"
    echo "Assembly excerpts:       ${DOC_PATH_DIR}"
}

build_x86_64() {
    echo "Building x86-64 assembly (Zen 2, AVX2 + FMA)..."
    (
        cd "${REPO_ROOT}"
        env CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_RUNNER=env \
            cargo rustc --release -p aether-fpx-research --no-default-features --features simd \
            --bin dump \
            --target x86_64-unknown-linux-gnu \
            --target-dir "${X86_TARGET_DIR}" -- \
            --emit=asm \
            -C target-cpu=znver2 \
            -C target-feature=+avx2,+fma
    )
}

build_aarch64() {
    echo "Building AArch64 assembly (Cortex-A76, NEON)..."
    (
        cd "${REPO_ROOT}"
        cargo rustc --release -p aether-fpx-research --no-default-features --features simd \
            --bin dump \
            --target aarch64-unknown-linux-musl \
            --target-dir "${AARCH64_TARGET_DIR}" -- \
            --emit=asm \
            -C target-cpu=cortex-a76 \
            -C target-feature=+neon \
            -C linker=rust-lld
    )
}

find_dump_assembly() {
    local target_dir="$1"
    local assembly
    assembly="$(
        find "${target_dir}" -type f -name 'dump.dump.*.s' -printf '%T@ %p\n' \
            | sort -nr \
            | sed -n '1p' \
            | cut -d ' ' -f 2-
    )"
    if [[ -z "${assembly}" ]]; then
        echo "No dump assembly found under ${target_dir}; build that target first." >&2
        exit 1
    fi
    printf '%s\n' "${assembly}"
}

analyze_x86_64() {
    local assembly="$1"
    awk -v arch="x86_64" '
        /^dump_[a-z0-9_]+_(horner|estrin)_f(32|64):$/ {
            fn = $0
            sub(/:$/, "", fn)
            active = 1
            scalar_fma = scalar_mul = scalar_addsub = 0
            packed_fma = packed_mul = packed_addsub = 0
            next
        }
        active {
            mnemonic = $1
            if (mnemonic ~ /^vfm(add|sub).*s[sd]$/) scalar_fma++
            else if (mnemonic ~ /^vfm(add|sub).*p[sd]$/) packed_fma++
            else if (mnemonic ~ /^vmuls[sd]$/) scalar_mul++
            else if (mnemonic ~ /^vmulp[sd]$/) packed_mul++
            else if (mnemonic ~ /^v(add|sub)s[sd]$/) scalar_addsub++
            else if (mnemonic ~ /^v(add|sub)p[sd]$/) packed_addsub++
        }
        active && /^\.Lfunc_end/ {
            printf "%s,%s,%d,%d,%d,%d,%d,%d\n",
                arch, fn, scalar_fma, scalar_mul, scalar_addsub,
                packed_fma, packed_mul, packed_addsub
            active = 0
        }
    ' "${assembly}"
}

analyze_aarch64() {
    local assembly="$1"
    awk -v arch="aarch64" '
        /^dump_[a-z0-9_]+_(horner|estrin)_f(32|64):$/ {
            fn = $0
            sub(/:$/, "", fn)
            active = 1
            scalar_fma = scalar_mul = scalar_addsub = 0
            packed_fma = packed_mul = packed_addsub = 0
            next
        }
        active {
            mnemonic = $1
            operand = $2
            if (mnemonic ~ /^fm(add|sub)$/ && operand ~ /^[sd][0-9]+,/) scalar_fma++
            else if (mnemonic ~ /^(fm(add|sub)|fml[as])$/ && operand ~ /^v[0-9]+\./) packed_fma++
            else if (mnemonic == "fmul" && operand ~ /^[sd][0-9]+,/) scalar_mul++
            else if (mnemonic == "fmul" && operand ~ /^v[0-9]+\./) packed_mul++
            else if (mnemonic ~ /^f(add|sub)$/ && operand ~ /^[sd][0-9]+,/) scalar_addsub++
            else if (mnemonic ~ /^f(add|sub)$/ && operand ~ /^v[0-9]+\./) packed_addsub++
        }
        active && /^\.Lfunc_end/ {
            printf "%s,%s,%d,%d,%d,%d,%d,%d\n",
                arch, fn, scalar_fma, scalar_mul, scalar_addsub,
                packed_fma, packed_mul, packed_addsub
            active = 0
        }
    ' "${assembly}"
}

analyze() {
    local x86_assembly
    local aarch64_assembly
    x86_assembly="$(find_dump_assembly "${X86_TARGET_DIR}")"
    aarch64_assembly="$(find_dump_assembly "${AARCH64_TARGET_DIR}")"
    mkdir -p "${REPORT_DIR}"
    {
        echo "architecture,function,scalar_fma,scalar_mul,scalar_addsub,packed_fma,packed_mul,packed_addsub"
        analyze_x86_64 "${x86_assembly}"
        analyze_aarch64 "${aarch64_assembly}"
    } > "${REPORT_PATH}"

    echo
    echo "Assembly inputs:"
    echo "  x86_64:  ${x86_assembly}"
    echo "  aarch64: ${aarch64_assembly}"
    echo
    awk -F, '
        NR > 1 {
            method = ($2 ~ /_horner_/) ? "horner" : "estrin"
            kind = ($2 ~ /_raw_simd_/) ? "explicit_simd" : "scalar_api"
            key = $1 SUBSEP kind SUBSEP method
            anchors[key]++
            scalar_fma[key] += $3
            scalar_mul[key] += $4
            scalar_addsub[key] += $5
            packed[key] += $6 + $7 + $8
            scalar_total += $3 + $4 + $5
            packed_total += $6 + $7 + $8
        }
        END {
            simd_ok = 1
            print "architecture  api            method  anchors  scalar_fma  scalar_mul  scalar_addsub  packed_total"
            keys[1] = "x86_64" SUBSEP "scalar_api" SUBSEP "horner"
            keys[2] = "x86_64" SUBSEP "scalar_api" SUBSEP "estrin"
            keys[3] = "x86_64" SUBSEP "explicit_simd" SUBSEP "horner"
            keys[4] = "x86_64" SUBSEP "explicit_simd" SUBSEP "estrin"
            keys[5] = "aarch64" SUBSEP "scalar_api" SUBSEP "horner"
            keys[6] = "aarch64" SUBSEP "scalar_api" SUBSEP "estrin"
            keys[7] = "aarch64" SUBSEP "explicit_simd" SUBSEP "horner"
            keys[8] = "aarch64" SUBSEP "explicit_simd" SUBSEP "estrin"
            for (i = 1; i <= 8; i++) {
                split(keys[i], parts, SUBSEP)
                key = keys[i]
                printf "%-12s  %-13s  %-6s  %7d  %10d  %10d  %13d  %12d\n",
                    parts[1], parts[2], parts[3], anchors[key], scalar_fma[key],
                    scalar_mul[key], scalar_addsub[key], packed[key]
                if (parts[2] == "explicit_simd" && (anchors[key] != 2 || packed[key] == 0)) {
                    simd_ok = 0
                }
            }
            printf "\nPolynomial arithmetic totals: scalar=%d, packed=%d\n",
                scalar_total, packed_total
            if (!simd_ok) {
                print "Result: ERROR: at least one architecture/method lacks packed arithmetic."
                exit 1
            }
            print "Result: packed arithmetic is present for both methods on both architectures."
        }
    ' "${REPORT_PATH}"
    echo "Report: ${REPORT_PATH}"
    generate_path_dumps
}

case "${1:-all}" in
    all)
        build_x86_64
        build_aarch64
        analyze
        ;;
    x86_64)
        build_x86_64
        ;;
    aarch64)
        build_aarch64
        ;;
    analyze)
        analyze
        ;;
    paths)
        generate_path_dumps
        ;;
    *)
        usage
        ;;
esac
