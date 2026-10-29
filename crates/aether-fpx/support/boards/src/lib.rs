#![cfg_attr(feature = "trig-probe", feature(f16, f128))]
#![no_std]

#[cfg(all(feature = "pico2-arm", feature = "pico2-hazard3"))]
compile_error!("Enable only one of 'pico2-arm' or 'pico2-hazard3' at a time.");

#[cfg(feature = "trig-probe")]
pub mod probe;

#[cfg(feature = "trig-bench-rtt")]
pub mod rtt_bench;

#[cfg(feature = "pico2-arm")]
pub mod pico2_arm {
    pub const BOARD_NAME: &str = "Raspberry Pi Pico 2";
    pub const MCU_NAME: &str = "RP2350";
    pub const TARGET_TRIPLE: &str = "thumbv8m.main-none-eabihf";

    #[inline]
    pub const fn cortex_m33() -> bool {
        true
    }
}

#[cfg(feature = "pico2-hazard3")]
pub mod pico2_hazard3 {
    pub const BOARD_NAME: &str = "Raspberry Pi Pico 2";
    pub const MCU_NAME: &str = "RP2350";
    pub const TARGET_TRIPLE: &str = "riscv32imac-unknown-none-elf";

    #[inline]
    pub const fn hazard3() -> bool {
        true
    }
}
