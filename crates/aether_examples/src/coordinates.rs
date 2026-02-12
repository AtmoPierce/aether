use aether_core::attitude::DirectionCosineMatrix;
use aether_core::coordinate::Cartesian;
use aether_core::reference_frame::{Body, NED, ITRF};
use aether_core::reference_frame::transforms::{body_to_ned, itrf_to_ned};

fn main() {
    // Example 1: gravity measured in the BODY frame
    //
    // Imagine an IMU on a vehicle. In body coordinates, it "sees" gravity
    // along +X (nose-forward) due to pitch attitude.
    let gravity_body: Cartesian<f64, Body<f64>> = Cartesian::new(9.8, 0.0, 0.0);

    // Current attitude: roll, pitch, yaw (rad)
    let roll  = 0.0_f64.to_radians();
    let pitch = -90.0_f64.to_radians(); // nose straight down
    let yaw   = 0.0_f64.to_radians();

    // Build the Body -> NED direction cosine matrix from Euler angles.
    // body_to_ned() returns a DirectionCosineMatrix<Body, NED>.
    let dcm_body_to_ned: DirectionCosineMatrix<f64, Body<f64>, NED<f64>> =
        body_to_ned(roll, pitch, yaw);

    // Rotate the measurement into NED coordinates.
    let gravity_ned: Cartesian<f64, NED<f64>> = dcm_body_to_ned * gravity_body;

    println!("Gravity in body frame: {}", gravity_body);
    println!("Gravity in NED frame:  {}", gravity_ned);
    // Result:
    // Gravity in body frame: Cartesian [x,y,z]: [9.8, 0, 0]
    // Gravity in NED frame:  Cartesian [x,y,z]: [0.0000000000000006000769315822031, 0, -9.8]

    // At this point, gravity_ned is tagged as NED<f64>, so downstream code
    // cannot accidentally treat it as if it's already in Earth-fixed or inertial
    // coordinates. The frame is part of the type.

    // ---------------------------------------------------------
    // Example 2: connect local navigation to Earth-fixed
    // ---------------------------------------------------------

    // Suppose we know where we are on Earth:
    let latitude  = 40.0_f64.to_radians();
    let longitude = -90.0_f64.to_radians();

    // Build the ECEF(ITRF) -> NED transform at that geodetic location.
    // itrf_to_ned() returns DirectionCosineMatrix<ITRF, NED>.
    let dcm_itrf_to_ned: DirectionCosineMatrix<f64, ITRF<f64>, NED<f64>> =
        itrf_to_ned(latitude, longitude);

    // If we had some global force/velocity/etc. in ITRF/ECEF, we could move it into
    // local navigation frame with dcm_itrf_to_ned. Conversely, its inverse moves local NED vectors back to Earth-fixed coordinates.
    // inverse/transpose is done by the user, and not automatically stored.

    // This gives you an auditable pipeline like:
    // Body --> NED --> ITRF --> ICRF
    // with no chance of "oops I mixed frames" bugs.

    
}