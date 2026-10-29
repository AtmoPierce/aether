mod remez_tables {
    include!(concat!(env!("OUT_DIR"), "/fpx_sin.rs"));
    include!(concat!(env!("OUT_DIR"), "/fpx_cos.rs"));
    include!(concat!(env!("OUT_DIR"), "/fpx_tan.rs"));
    include!(concat!(env!("OUT_DIR"), "/fpx_exp.rs"));
}

pub use remez_tables::*;
