mod config;
mod script;
mod submit;
mod target;

use std::path::Path;

use config::PbsConfig;
use script::PbsJob;
use target::PbsTarget;

pub fn run(args: &[String]) -> Result<(), String> {
    if matches!(args.first().map(String::as_str), Some("-h" | "--help")) {
        println!("usage: aether-runner pbs --target <rust-target> <executable> [arguments...]");
        return Ok(());
    }

    let [target_option, rust_target, executable, arguments @ ..] = args else {
        return Err(
            "usage: aether-runner pbs --target <rust-target> <executable> [arguments...]"
                .to_string(),
        );
    };

    if target_option != "--target" {
        return Err(format!(
            "expected `--target`, found `{target_option}`\nusage: aether-runner pbs --target <rust-target> <executable> [arguments...]"
        ));
    }

    let target = PbsTarget::from_rust_target(rust_target)?;
    let config = PbsConfig::from_env();
    let job = PbsJob {
        config: &config,
        target: &target,
        executable: Path::new(executable),
        arguments,
    };
    let script = job.render()?;
    let submitted_job = submit::submit(&config, &script)?;

    println!("submitted PBS job {}", submitted_job.id);
    Ok(())
}
