// Bridge module to match expected `crate::models::...` paths

pub mod terrestrial {
    pub use crate::terrestrial::*;
}

pub mod celestial {
    pub mod constants {
        pub use crate::celestial::constants::*;
    }
}
