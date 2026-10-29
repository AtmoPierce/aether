use aether_fpx_gen::generators::*;

fn main() {
    println!("cargo:rerun-if-changed=../../src/core/");
    println!("cargo:rerun-if-changed=../../src/core/trig");
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=fpx-remez.toml");
    println!("cargo:rerun-if-changed=./inputs");
    println!("cargo:rerun-if-changed=./templates");

    let force = std::env::var("FPX_GEN").is_ok();

    if force {
        println!("cargo:warning=FPX_GEN set: forcing codegen");
    }

    generate_all_from_config();
}
