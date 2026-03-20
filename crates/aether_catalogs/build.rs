use std::{
    fs,
    path::Path,
    process::{Command, Stdio},
};

fn main() {
    usno_gnc_catalog();
}

fn usno_gnc_catalog() {
    download_if_missing(
        "catalog/gnc_v1_1_mar_2_2023.csv",
        "https://crf.usno.navy.mil/data_products/ICRF/GNC/2023/gnc_v1_1_mar_2_2023.csv",
    );
}

fn download_if_missing<P: AsRef<Path>>(local_path: P, url: &str) {
    let path = local_path.as_ref();

    if path.exists() {
        return;
    }

    let dir = path.parent().unwrap();
    fs::create_dir_all(dir).unwrap();

    println!("cargo:warning=Starting download: {}", url);

    let mut cmd = Command::new("curl");
    cmd.args([
        "-L",
        "-C",
        "-",
        "--connect-timeout",
        "5",
        "--max-time",
        "600",
        "--retry",
        "3",
        "--retry-delay",
        "2",
        "--retry-connrefused",
        "--progress-bar",
        "-o",
        path.to_str().unwrap(),
        url,
    ]);

    // Try to show live output if terminal allows
    cmd.stdout(Stdio::inherit()).stderr(Stdio::inherit());

    let status = cmd.status().expect("Failed to run curl");

    if !status.success() {
        panic!("Download failed for: {}", url);
    }

    println!("cargo:warning=Finished download: {}", url);
}