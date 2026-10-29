#![no_std]
#![no_main]
#![feature(f16, f128)]

use panic_halt as _;

#[unsafe(no_mangle)]
pub extern "C" fn main() -> ! {
    aether_fpx_boards::probe::run_trig_probe()
}
