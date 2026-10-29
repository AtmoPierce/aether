use std::{
    fs,
    path::PathBuf,
    process::Command,
    time::{SystemTime, UNIX_EPOCH},
};

use super::config::PbsConfig;

#[derive(Debug)]
pub struct SubmittedJob {
    pub id: String,
}

pub fn submit(config: &PbsConfig, script: &str) -> Result<SubmittedJob, String> {
    fs::create_dir_all(&config.script_directory).map_err(|error| {
        format!(
            "failed to create PBS script directory {}: {error}",
            config.script_directory.display()
        )
    })?;

    let script_path = create_script_path(config);

    fs::write(&script_path, script).map_err(|error| {
        format!(
            "failed to write PBS script {}: {error}",
            script_path.display()
        )
    })?;

    let output = Command::new("qsub")
        .arg("-V")
        .arg(&script_path)
        .output()
        .map_err(|error| format!("failed to execute qsub: {error}"))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);

        return Err(format!("qsub failed:\n{stderr}"));
    }

    let stdout = String::from_utf8_lossy(&output.stdout);

    let job_id = stdout.trim();

    if job_id.is_empty() {
        return Err("qsub succeeded but returned no job ID".to_string());
    }

    // qsub has already consumed the script by this point.
    let _ = fs::remove_file(&script_path);

    Ok(SubmittedJob {
        id: job_id.to_string(),
    })
}

fn create_script_path(config: &PbsConfig) -> PathBuf {
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();

    config
        .script_directory
        .join(format!("cargo-{timestamp}.pbs"))
}
