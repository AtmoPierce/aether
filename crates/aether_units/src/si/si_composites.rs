#[derive(Copy, Clone, Debug, PartialEq)]
pub enum MeterPerSecond {}              // velocity
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum MeterPerSecondSquared {}       // acceleration
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum RadianPerSecond {}             // angular rate
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum RadianPerSecondSquared {}      // angular accel

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum NewtonMeter {}                 // torque
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum NewtonSecond {}                // impulse

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum CubicMeterPerSecond {}         // volume flow
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum KilogramPerCubicMeter {}       // density
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum PascalSecond {}                // dynamic viscosity
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum SquareMeterPerSecond {}        // kinematic viscosity
