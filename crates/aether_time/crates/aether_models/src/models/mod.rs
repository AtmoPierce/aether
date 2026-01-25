// Bridge module to match expected `crate::models::...` paths

pub mod terrestial {
    pub use crate::terrestial::*;
}

pub mod celestial {
    pub mod constants {
        pub use crate::celestial::constants::*;
    }
}
