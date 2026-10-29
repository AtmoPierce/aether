use std::{env, path::PathBuf};

#[derive(Debug, Clone)]
pub struct PbsConfig {
    pub queue: Option<String>,
    pub account: Option<String>,

    pub ncpus: usize,
    pub memory: String,
    pub walltime: String,

    pub target_resource: String,

    pub script_directory: PathBuf,
}

impl PbsConfig {
    pub fn from_env() -> Self {
        let queue = env::var("PBS_RUNNER_QUEUE").ok();
        let account = env::var("PBS_RUNNER_ACCOUNT").ok();

        let ncpus = env::var("PBS_RUNNER_NCPUS")
            .ok()
            .and_then(|value| value.parse().ok())
            .unwrap_or(1);

        let memory = env::var("PBS_RUNNER_MEMORY").unwrap_or_else(|_| "2gb".to_string());

        let walltime = env::var("PBS_RUNNER_WALLTIME").unwrap_or_else(|_| "10:00:00".to_string());

        let target_resource =
            env::var("PBS_RUNNER_TARGET_RESOURCE").unwrap_or_else(|_| "cpu_arch".to_string());

        let script_directory = env::var("PBS_RUNNER_SCRIPT_DIRECTORY")
            .map(PathBuf::from)
            .unwrap_or_else(|_| env::temp_dir().join("aether-runner"));

        Self {
            queue,
            account,
            ncpus,
            memory,
            walltime,
            target_resource,
            script_directory,
        }
    }
}
