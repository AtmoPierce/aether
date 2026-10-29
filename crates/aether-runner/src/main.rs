mod forgejo;
mod pbs;

use std::{env, process::ExitCode};

fn main() -> ExitCode {
    let mut args = env::args().skip(1);

    let Some(runner) = args.next() else {
        eprintln!("usage: aether-runner <runner> [arguments...]");
        return ExitCode::FAILURE;
    };

    let args: Vec<String> = args.collect();

    let result = match runner.as_str() {
        "pbs" => pbs::run(&args),
        "forgejo" => forgejo::run(&args),

        unknown => {
            eprintln!("unknown runner: {unknown}");
            return ExitCode::FAILURE;
        }
    };

    match result {
        Ok(()) => ExitCode::SUCCESS,

        Err(error) => {
            eprintln!("runner error: {error}");
            ExitCode::FAILURE
        }
    }
}
