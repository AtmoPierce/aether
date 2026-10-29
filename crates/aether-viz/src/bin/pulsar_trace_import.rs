use std::error::Error;
use std::io::{Error as IoError, ErrorKind};
use std::path::PathBuf;

use aether_viz::rtos::{import_probe_log, load_rtos_run};

fn usage_error() -> Box<dyn Error> {
    IoError::new(
        ErrorKind::InvalidInput,
        "usage: pulsar_trace_import <probe-output.log> <run-directory>",
    )
    .into()
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut arguments = std::env::args_os().skip(1);
    let log = PathBuf::from(arguments.next().ok_or_else(usage_error)?);
    let output = PathBuf::from(arguments.next().ok_or_else(usage_error)?);
    if arguments.next().is_some() {
        return Err(usage_error());
    }

    import_probe_log(&log, &output)?;
    let run = load_rtos_run(&output)?;
    println!(
        "imported {} records for run {} ({})",
        run.records.len(),
        run.run_id,
        run.metadata.get("data_origin").unwrap_or("unknown origin")
    );
    for warning in &run.warnings {
        eprintln!("warning: {warning}");
    }
    Ok(())
}
