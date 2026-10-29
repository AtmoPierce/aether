use std::{env, path::PathBuf};

fn main() {
    let manifest_dir =
        PathBuf::from(env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR not set"));
    let memory_dir = manifest_dir.join("memory");

    println!(
        "cargo:rerun-if-changed={}",
        memory_dir.join("pico2_arm.X").display()
    );
    println!(
        "cargo:rerun-if-changed={}",
        memory_dir.join("pico2_hazard3.X").display()
    );

    let target = env::var("TARGET").unwrap_or_default();
    let arm_enabled = env::var_os("CARGO_FEATURE_PICO2_ARM").is_some();
    let hazard3_enabled = env::var_os("CARGO_FEATURE_PICO2_HAZARD3").is_some();

    if arm_enabled && target.starts_with("thumbv8m") {
        let arm_script = memory_dir.join("pico2_arm.X");
        println!("cargo:rustc-link-arg-bins=-T{}", arm_script.display());
        println!("cargo:rustc-link-arg-benches=-T{}", arm_script.display());
    }

    if hazard3_enabled && target.starts_with("riscv32") {
        let hazard3_script = memory_dir.join("pico2_hazard3.X");
        println!("cargo:rustc-link-arg-bins=-T{}", hazard3_script.display());
        println!(
            "cargo:rustc-link-arg-benches=-T{}",
            hazard3_script.display()
        );
    }
}
