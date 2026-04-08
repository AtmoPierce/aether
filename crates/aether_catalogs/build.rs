use std::{
    env,
    fs,
    path::{Path, PathBuf},
    process::{Command, Stdio},
};

const WMM_MAX_DEGREE: usize = 12;
const WMM_RECORD_COUNT: usize = 90;
const WMM_TRIANGULAR_COUNT: usize = (WMM_MAX_DEGREE + 1) * (WMM_MAX_DEGREE + 2) / 2;
const GENERATED_WMM_RELATIVE_PATH: &str = "src/wmm/wmm_generated.rs";

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed={GENERATED_WMM_RELATIVE_PATH}");
    println!("cargo:rerun-if-env-changed=CARGO_FEATURE_DOWNLOAD_GNC_CATALOG");
    println!("cargo:rerun-if-env-changed=CARGO_FEATURE_DOWNLOAD_WORLD_MAGNETIC_MODEL");
    println!("cargo:rerun-if-env-changed=CARGO_FEATURE_DOWNLOAD_CATALOGS");

    if feature_enabled("DOWNLOAD_GNC_CATALOG") || feature_enabled("DOWNLOAD_CATALOGS") {
        usno_gnc_catalog();
    }

    if feature_enabled("DOWNLOAD_WORLD_MAGNETIC_MODEL") || feature_enabled("DOWNLOAD_CATALOGS") {
        world_magnetic_model();
    }

    if feature_enabled("WORLD_MAGNETIC_MODEL")
        || feature_enabled("DOWNLOAD_WORLD_MAGNETIC_MODEL")
        || feature_enabled("DOWNLOAD_CATALOGS")
    {
        generate_embedded_wmm();
    }
}

fn usno_gnc_catalog() {
    download_if_missing(
        "catalog/gnc_v1_1_mar_2_2023.csv",
        "https://crf.usno.navy.mil/data_products/ICRF/GNC/2023/gnc_v1_1_mar_2_2023.csv",
    );
}

fn world_magnetic_model(){
    let zip_path = Path::new("catalog/wmm_coefficients.zip");
    let extracted_dir = Path::new("catalog/wmm");

    download_if_missing(
        zip_path,
        "https://www.ncei.noaa.gov/sites/default/files/2024-12/WMM2025COF.zip",
    );
    extract_zip_if_needed(zip_path, extracted_dir);
}

fn feature_enabled(feature: &str) -> bool {
    env::var_os(format!("CARGO_FEATURE_{feature}")).is_some()
}

fn generate_embedded_wmm() {
    let output_path = PathBuf::from(GENERATED_WMM_RELATIVE_PATH);
    if output_path.is_file() {
        return;
    }

    let source_path = find_wmm_source_file().expect("unable to locate a WMM coefficient file to embed");
    println!("cargo:rerun-if-changed={}", source_path.display());

    let contents = fs::read_to_string(&source_path)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", source_path.display()));

    let generated = render_wmm_include(&contents, &source_path);
    if let Some(parent) = output_path.parent() {
        fs::create_dir_all(parent)
            .unwrap_or_else(|error| panic!("failed to create {}: {error}", parent.display()));
    }
    fs::write(&output_path, generated)
        .unwrap_or_else(|error| panic!("failed to write {}: {error}", output_path.display()));
    println!("cargo:warning=Generated {} from {}", output_path.display(), source_path.display());
}

fn find_wmm_source_file() -> Option<PathBuf> {
    let candidates = [
        PathBuf::from("catalog/wmm/WMM2025COF/WMM.COF"),
        PathBuf::from("catalog/wmm/WMM2025COF/WMM2025.COF"),
        PathBuf::from("catalog/wmm/WMM.COF"),
        PathBuf::from("catalog/wmm/WMM2025.COF"),
    ];

    candidates.into_iter().find(|path| path.is_file())
}

fn render_wmm_include(contents: &str, source_path: &Path) -> String {
    let mut non_empty_lines = contents.lines().filter(|line| !line.trim().is_empty());
    let header_line = non_empty_lines
        .next()
        .unwrap_or_else(|| panic!("{} is missing the WMM header", source_path.display()));
    let header_tokens: Vec<_> = header_line.split_whitespace().collect();
    assert!(header_tokens.len() >= 3, "invalid WMM header in {}", source_path.display());

    let epoch = header_tokens[0]
        .parse::<f64>()
        .unwrap_or_else(|error| panic!("invalid WMM epoch in {}: {error}", source_path.display()));
    let model_name = header_tokens[1];
    let release_date = header_tokens[2];

    let mut records = Vec::with_capacity(WMM_RECORD_COUNT);
    let mut g = [0.0_f64; WMM_TRIANGULAR_COUNT];
    let mut h = [0.0_f64; WMM_TRIANGULAR_COUNT];
    let mut g_dot = [0.0_f64; WMM_TRIANGULAR_COUNT];
    let mut h_dot = [0.0_f64; WMM_TRIANGULAR_COUNT];
    let mut max_degree = 0usize;

    for (line_number, raw_line) in contents.lines().enumerate().skip(1) {
        let line = raw_line.trim();
        if line.is_empty() || !line.chars().next().unwrap_or_default().is_ascii_digit() {
            continue;
        }

        let tokens: Vec<_> = line.split_whitespace().collect();
        if tokens.len() != 6 {
            continue;
        }

        let degree = tokens[0]
            .parse::<usize>()
            .unwrap_or_else(|error| panic!("invalid degree on line {} in {}: {error}", line_number + 1, source_path.display()));
        let order = tokens[1]
            .parse::<usize>()
            .unwrap_or_else(|error| panic!("invalid order on line {} in {}: {error}", line_number + 1, source_path.display()));
        assert!(order <= degree, "invalid n/m pair on line {} in {}", line_number + 1, source_path.display());

        let g_nm = tokens[2]
            .parse::<f64>()
            .unwrap_or_else(|error| panic!("invalid g_nm on line {} in {}: {error}", line_number + 1, source_path.display()));
        let h_nm = tokens[3]
            .parse::<f64>()
            .unwrap_or_else(|error| panic!("invalid h_nm on line {} in {}: {error}", line_number + 1, source_path.display()));
        let g_dot_nm = tokens[4]
            .parse::<f64>()
            .unwrap_or_else(|error| panic!("invalid g_dot_nm on line {} in {}: {error}", line_number + 1, source_path.display()));
        let h_dot_nm = tokens[5]
            .parse::<f64>()
            .unwrap_or_else(|error| panic!("invalid h_dot_nm on line {} in {}: {error}", line_number + 1, source_path.display()));

        let index = triangular_index(degree, order);
        g[index] = g_nm;
        h[index] = h_nm;
        g_dot[index] = g_dot_nm;
        h_dot[index] = h_dot_nm;
        max_degree = max_degree.max(degree);
        records.push((degree, order, g_nm, h_nm, g_dot_nm, h_dot_nm));
    }

    assert_eq!(max_degree, WMM_MAX_DEGREE, "unexpected WMM max degree in {}", source_path.display());
    assert_eq!(records.len(), WMM_RECORD_COUNT, "unexpected WMM record count in {}", source_path.display());

    let mut output = String::new();
    output.push_str("pub const GENERATED_WMM_HEADER: WmmHeader = WmmHeader {\n");
    output.push_str(&format!("    epoch: {epoch:?},\n"));
    output.push_str(&format!("    model_name: {model_name:?},\n"));
    output.push_str(&format!("    release_date: {release_date:?},\n"));
    output.push_str("};\n\n");

    output.push_str(&format!(
        "pub const GENERATED_WMM_RECORDS: [WmmCoefficientRecord; {WMM_RECORD_COUNT}] = [\n"
    ));
    for (degree, order, g_nm, h_nm, g_dot_nm, h_dot_nm) in records {
        output.push_str(&format!(
            "    WmmCoefficientRecord {{ degree: {degree}, order: {order}, g_nm_nt: {g_nm:?}, h_nm_nt: {h_nm:?}, g_dot_nm_nt_per_year: {g_dot_nm:?}, h_dot_nm_nt_per_year: {h_dot_nm:?} }},\n"
        ));
    }
    output.push_str("];\n\n");

    write_array(&mut output, "GENERATED_WMM_G", &g);
    write_array(&mut output, "GENERATED_WMM_H", &h);
    write_array(&mut output, "GENERATED_WMM_G_DOT", &g_dot);
    write_array(&mut output, "GENERATED_WMM_H_DOT", &h_dot);
    output
}

fn triangular_index(degree: usize, order: usize) -> usize {
    degree * (degree + 1) / 2 + order
}

fn write_array(output: &mut String, name: &str, values: &[f64; WMM_TRIANGULAR_COUNT]) {
    output.push_str(&format!(
        "pub const {name}: [f64; {WMM_TRIANGULAR_COUNT}] = [\n"
    ));
    for value in values {
        output.push_str(&format!("    {value:?},\n"));
    }
    output.push_str("];\n\n");
}

fn extract_zip_if_needed<P: AsRef<Path>, Q: AsRef<Path>>(archive_path: P, output_dir: Q) {
    let archive_path = archive_path.as_ref();
    let output_dir = output_dir.as_ref();

    if output_dir.exists() && output_dir.read_dir().map(|mut entries| entries.next().is_some()).unwrap_or(false) {
        return;
    }

    fs::create_dir_all(output_dir).unwrap();

    println!(
        "cargo:warning=Extracting archive: {} -> {}",
        archive_path.display(),
        output_dir.display()
    );

    let status = Command::new("unzip")
        .args(["-o", archive_path.to_str().unwrap(), "-d", output_dir.to_str().unwrap()])
        .stdout(Stdio::inherit())
        .stderr(Stdio::inherit())
        .status()
        .expect("Failed to run unzip");

    if !status.success() {
        panic!("Extraction failed for: {}", archive_path.display());
    }

    println!("cargo:warning=Finished extracting: {}", archive_path.display());
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

    // show live output if terminal allows
    cmd.stdout(Stdio::inherit()).stderr(Stdio::inherit());

    let status = cmd.status().expect("Failed to run curl");

    if !status.success() {
        panic!("Download failed for: {}", url);
    }

    println!("cargo:warning=Finished download: {}", url);
}