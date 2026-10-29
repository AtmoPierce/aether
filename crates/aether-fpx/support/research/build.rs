use std::env;

fn main() {
    println!("cargo:rustc-check-cfg=cfg(fpx_x87)");
    println!("cargo:rerun-if-env-changed=CARGO_FEATURE_X87");
    println!("cargo:rerun-if-changed=c/x87/sin_x87.c");
    println!("cargo:rerun-if-changed=c/x87/cos_x87.c");
    println!("cargo:rerun-if-changed=c/x87/tan_x87.c");
    println!("cargo:rerun-if-changed=c/x87/sqrt_x87.c");

    if env::var_os("CARGO_FEATURE_X87").is_none() {
        return;
    }

    let target_arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap_or_default();
    let target_env = env::var("CARGO_CFG_TARGET_ENV").unwrap_or_default();
    let supported_arch = matches!(target_arch.as_str(), "x86" | "x86_64");
    let supported_environment = target_env != "msvc";

    if !supported_arch || !supported_environment {
        println!(
            "cargo:warning=x87 benchmarks disabled for target architecture {target_arch} ({target_env})"
        );
        return;
    }

    let mut build = cc::Build::new();

    build
        .file("c/x87/tan_x87.c")
        .file("c/x87/sin_x87.c")
        .file("c/x87/cos_x87.c")
        .file("c/x87/sqrt_x87.c")
        .flag("-std=c11")
        .flag_if_supported("-mfpmath=387")
        .flag("-fno-builtin")
        .flag("-ffp-contract=off")
        // determinism
        .flag("-fno-fast-math")
        .flag("-frounding-math")
        .flag("-fexcess-precision=standard")
        // debug / research
        .flag("-O2")
        .flag("-fverbose-asm")
        .flag("-save-temps=obj");

    build.compile("fpx_x87");
    println!("cargo:rustc-cfg=fpx_x87");
}
