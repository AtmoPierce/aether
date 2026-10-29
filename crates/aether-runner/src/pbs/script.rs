use std::path::Path;

use super::{config::PbsConfig, target::PbsTarget};

pub struct PbsJob<'a> {
    pub config: &'a PbsConfig,
    pub target: &'a PbsTarget,

    pub executable: &'a Path,
    pub arguments: &'a [String],
}

impl PbsJob<'_> {
    pub fn render(&self) -> Result<String, String> {
        let executable = self.executable.canonicalize().map_err(|error| {
            format!(
                "failed to canonicalize executable {}: {error}",
                self.executable.display()
            )
        })?;

        let job_name = executable
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("cargo-run");

        let working_directory = std::env::current_dir()
            .map_err(|error| format!("failed to determine working directory: {error}"))?;

        let mut script = String::new();

        script.push_str("#!/bin/bash\n");

        script.push_str(&format!("#PBS -N {}\n", sanitize_job_name(job_name)));

        if let Some(queue) = &self.config.queue {
            script.push_str(&format!("#PBS -q {}\n", queue));
        }

        if let Some(account) = &self.config.account {
            script.push_str(&format!("#PBS -A {}\n", account));
        }

        script.push_str(&format!(
            "#PBS -l select=1:ncpus={}:mem={}:{}={}\n",
            self.config.ncpus,
            self.config.memory,
            self.config.target_resource,
            self.target.resource_value(),
        ));

        script.push_str(&format!("#PBS -l walltime={}\n", self.config.walltime));

        script.push_str("#PBS -j oe\n");

        script.push('\n');

        script.push_str("set -e\n\n");

        script.push_str("cd ");
        script.push_str(&shell_escape(working_directory.to_string_lossy().as_ref()));
        script.push('\n');

        script.push_str("aether_runner_original_executable=");
        script.push_str(&shell_escape(executable.to_string_lossy().as_ref()));
        script.push('\n');
        script.push_str("aether_runner_tmpdir=\"${PBS_TMPDIR:-${TMPDIR:-/tmp}}\"\n");
        script.push_str("mkdir -p \"$aether_runner_tmpdir\"\n");
        script.push_str(
            "aether_runner_executable_name=$(basename \"$aether_runner_original_executable\")\n",
        );
        script.push_str("aether_runner_staged_executable=\"$aether_runner_tmpdir/aether-runner-${PBS_JOBID:-$$}-$aether_runner_executable_name\"\n");
        script.push_str(
            "cp \"$aether_runner_original_executable\" \"$aether_runner_staged_executable\"\n",
        );
        script.push_str("chmod +x \"$aether_runner_staged_executable\"\n");
        script.push_str("trap 'rm -f \"$aether_runner_staged_executable\"' EXIT\n");

        script.push_str("\"$aether_runner_staged_executable\"");

        for argument in self.arguments {
            script.push(' ');
            script.push_str(&shell_escape(argument));
        }

        script.push('\n');

        Ok(script)
    }
}

fn sanitize_job_name(name: &str) -> String {
    name.chars()
        .map(|character| {
            if character.is_ascii_alphanumeric() || character == '_' || character == '-' {
                character
            } else {
                '_'
            }
        })
        .collect()
}

fn shell_escape(value: &str) -> String {
    if value
        .chars()
        .all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '-' | '.' | '/' | ':'))
    {
        return value.to_string();
    }

    format!("'{}'", value.replace('\'', "'\"'\"'"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pbs::{config::PbsConfig, target::PbsTarget};
    use std::{
        fs,
        time::{SystemTime, UNIX_EPOCH},
    };

    #[test]
    fn pbs_script_stages_executable_before_running_it() {
        let executable = temp_file("aether-runner-script-test");
        fs::write(&executable, "#!/bin/sh\n").unwrap();
        let target = PbsTarget::from_rust_target("aarch64-unknown-linux-gnu").unwrap();
        let config = PbsConfig {
            queue: None,
            account: None,
            ncpus: 1,
            memory: "1gb".to_string(),
            walltime: "10:00:00".to_string(),
            target_resource: "cpu_arch".to_string(),
            script_directory: std::env::temp_dir(),
        };
        let arguments = vec!["--flag".to_string(), "value with spaces".to_string()];
        let job = PbsJob {
            config: &config,
            target: &target,
            executable: &executable,
            arguments: &arguments,
        };

        let script = job.render().unwrap();

        assert!(script.contains("aether_runner_original_executable="));
        assert!(script.contains(
            "cp \"$aether_runner_original_executable\" \"$aether_runner_staged_executable\""
        ));
        assert!(script.contains("chmod +x \"$aether_runner_staged_executable\""));
        assert!(script.contains("\"$aether_runner_staged_executable\" --flag 'value with spaces'"));

        let _ = fs::remove_file(executable);
    }

    fn temp_file(name: &str) -> std::path::PathBuf {
        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();

        std::env::temp_dir().join(format!("{name}-{timestamp}"))
    }
}
