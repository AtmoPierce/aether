pub mod utils;
use utils::*;

#[cfg(feature = "repl")]
fn main() {
    Repl::run();
}

#[cfg(not(feature = "repl"))]
fn main() {
    // Optionally, do nothing or print a message
    println!("Build with `--features cli` to enable the interactive shell.");
}
