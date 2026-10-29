mod forgejo;
mod local;
mod pbs;

use std::{env, process::ExitCode};

fn main() -> ExitCode {
    let mut args = env::args().skip(1);

    let Some(runner) = args.next() else {
        eprintln!("usage: aether-runner <local|pbs|forgejo> [arguments...]");
        return ExitCode::FAILURE;
    };

    let args: Vec<String> = args.collect();

    let result: Result<u8, String> = match runner.as_str() {
        "local" => local::run(&args),
        "pbs" => pbs::run(&args).map(|()| 0),
        "forgejo" => forgejo::run(&args).map(|()| 0),

        unknown => {
            eprintln!("unknown runner: {unknown}");
            return ExitCode::FAILURE;
        }
    };

    match result {
        Ok(code) => ExitCode::from(code),

        Err(error) => {
            eprintln!("runner error: {error}");
            ExitCode::FAILURE
        }
    }
}
