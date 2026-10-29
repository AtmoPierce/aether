pub mod euler;
pub mod rk4;
pub mod trapezoidal;
pub use euler::EulerIntegrate;
pub use rk4::Rk4Integrate;
pub use trapezoidal::TrapezoidalIntegrate;
