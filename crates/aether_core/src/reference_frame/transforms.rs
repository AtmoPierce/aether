use crate::attitude::DirectionCosineMatrix;
use crate::coordinate::Cartesian;
use crate::math::{Matrix, Vector};
use crate::matrix;
use crate::reference_frame::{Body, ReferenceFrame, ICRF, ITRF, NED};
use crate::real::Real;

pub fn angular_rate_dcm(roll: f64, pitch: f64, yaw: f64) -> Matrix<f64, 3, 3> {
    let rotation_matrix = matrix![
        0.0, -yaw, pitch;
        yaw, 0.0, -roll;
        -pitch, roll, 0.0;
    ];
    return rotation_matrix;
}

pub fn aircraft_wind_to_body(angle_of_attack: f64, side_slip: f64) -> Matrix<f64, 3, 3> {
    let t_bs = matrix![
        angle_of_attack.cos(), 0.0, -angle_of_attack.sin();
        0.0, 1.0, 0.0;
        angle_of_attack.sin(), 0.0, angle_of_attack.cos();
    ];

    let t_ws = matrix![
        side_slip.cos(), side_slip.sin(), 0.0;
        -side_slip.sin(), side_slip.cos(), 0.0;
        0.0, 0.0, 1.0;
    ];

    let t_wb = t_ws * t_bs.transpose();
    return t_wb;
}

pub fn aeroballistic_wind_to_body(
    angle_of_attack: f64,
    aerodynamic_roll_angle: f64,
) -> Matrix<f64, 3, 3> {
    let t_ab = matrix![
        angle_of_attack.cos(),  angle_of_attack.sin()*aerodynamic_roll_angle.sin(), angle_of_attack.sin()*aerodynamic_roll_angle.cos();
        0.0,                    aerodynamic_roll_angle.cos(),                       -aerodynamic_roll_angle.sin();
        -angle_of_attack.sin(), angle_of_attack.cos()*aerodynamic_roll_angle.sin(), angle_of_attack.cos()*aerodynamic_roll_angle.cos();
    ];
    return t_ab;
}

pub fn flight_path_to_geographic(heading_angle: f64, flight_path_angle: f64) -> Matrix<f64, 3, 3> {
    let rotation_matrix = matrix![
        flight_path_angle.cos()*heading_angle.cos(),    flight_path_angle.cos()*heading_angle.sin(),    -flight_path_angle.sin();
        -heading_angle.sin(),                           heading_angle.cos(),                            0.0;
        flight_path_angle.sin()*heading_angle.cos(),    flight_path_angle.sin()*heading_angle.sin(),    flight_path_angle.cos();
    ];
    return rotation_matrix;
}

pub fn body_to_ned(
    roll: f64,
    pitch: f64,
    yaw: f64,
) -> DirectionCosineMatrix<f64, Body<f64>, NED<f64>> {
    let phi = roll;
    let tht = pitch;
    let psi = yaw;
    let rotation_matrix = matrix![
        psi.cos()*tht.cos(),                                  psi.sin()*tht.cos(),                                     -tht.sin();
        psi.cos()*tht.sin()*phi.sin()-psi.sin()*phi.cos(),    psi.sin()*tht.sin()*phi.sin()+psi.cos()*phi.cos(),  tht.cos()*phi.sin();
        psi.cos()*tht.sin()*phi.cos()+psi.sin()*phi.sin(),    psi.sin()*tht.sin()*phi.cos()-psi.cos()*phi.sin(),  tht.cos()*phi.cos();
    ];
    let dcm: DirectionCosineMatrix<f64, Body<f64>, NED<f64>> =
        DirectionCosineMatrix::from_matrix(rotation_matrix);
    return dcm;
}

pub fn itrf_to_ned(
    latitude: f64,
    longitude: f64,
) -> DirectionCosineMatrix<f64, ITRF<f64>, NED<f64>> {
    let rotation_matrix = matrix![
        -latitude.sin()*longitude.cos(),    -latitude.sin()*longitude.sin(),    latitude.cos();
        -longitude.sin(),                   longitude.cos(),                    0.0;
        -latitude.cos()*longitude.cos(),    -latitude.cos()*longitude.sin(),    -latitude.sin();
    ];
    let dcm = DirectionCosineMatrix::from_matrix(rotation_matrix);
    return dcm;
}

pub fn icrf_to_itrf(
    time: f64,
    rotational_velocity: Cartesian<f64, ITRF<f64>>,
) -> DirectionCosineMatrix<f64, ICRF<f64>, ITRF<f64>> {
    return DirectionCosineMatrix::new(
        (rotational_velocity.z() * time).cos(),
        (rotational_velocity.z() * time).sin(),
        0.0,
        (-rotational_velocity.z() * time).sin(),
        (rotational_velocity.z() * time).cos(),
        0.0,
        0.0,
        0.0,
        1.0,
    );
}
pub fn itrf_to_icrf(
    time: f64,
    rotational_velocity: Cartesian<f64, ITRF<f64>>,
) -> DirectionCosineMatrix<f64, ITRF<f64>, ICRF<f64>> {
    return DirectionCosineMatrix::new(
        (rotational_velocity.z() * time).cos(),
        (-rotational_velocity.z() * time).sin(),
        0.0,
        (rotational_velocity.z() * time).sin(),
        (rotational_velocity.z() * time).cos(),
        0.0,
        0.0,
        0.0,
        1.0,
    );
}

#[cfg(test)]
mod tests {
    use crate::coordinate::Spherical;

    use super::*;
    use approx::relative_eq;
    #[test]
    fn test_aircraft_wind_to_body() {
        let wind = Vector::new([1.0, 1.0, 1.0]);
        let wind_magnitude: f64 = wind.magnitude();

        let angle_of_attack = wind[1].atan2(wind[0]);
        let aircraft_side_slip = (wind[1] / wind_magnitude).asin();

        let result = aircraft_wind_to_body(angle_of_attack, aircraft_side_slip);
        let aoa_result = result[0][2].atan2(result[0][0]);
        let sideslip_result = result[0][2].asin();

        let wut = relative_eq!(angle_of_attack, aoa_result, epsilon = 1e-3);
        assert_eq!(true, wut);
        let wut = relative_eq!(aircraft_side_slip, sideslip_result, epsilon = 1e-3);
        assert_eq!(true, wut);
    }

    #[test]
    fn test_aeroballistic_wind_to_body() {
        let wind = Vector::new([1.0, 1.0, 1.0]);
        let wind_magnitude: f64 = wind.magnitude();

        let angle_of_attack = wind[2].atan2(wind_magnitude);
        let roll_angle = wind[0].atan2(wind[2]);

        let result = aeroballistic_wind_to_body(angle_of_attack, roll_angle);
        let r1 = -result[0][2].asin();
        let r2 = result[0][1].atan2(result[0][0]);
        let r3 = result[1][2].atan2(result[2][2]);

        println!("Angles: {},{}", angle_of_attack, roll_angle);
        println!("R: {},{},{}", r1, r2, r3);

        // let aoa_result =
        // let roll_angle_result = result.m22.atan2(result.m21);
        // println!("AoA: {},{}", angle_of_attack, aoa_result);
        // println!("Roll: {},{}", roll_angle, roll_angle_result);

        // let wut = relative_eq!(angle_of_attack, aoa_result, epsilon = 1e-3);
        // assert_eq!(true, wut);
        // let wut = relative_eq!(roll_angle, roll_angle_result, epsilon = 1e-3);
        // assert_eq!(true, wut);
    }

    #[test]
    fn test_flight_path_to_geographic() {
        let geographic_velocity = Vector::new([100.0, 100.0, 0.0]);

        let heading = 45.0_f64.to_radians();
        let flight_path = 45.0_f64.to_radians();

        let result = flight_path_to_geographic(heading, flight_path);
        let velocity_coordinates = result * geographic_velocity;

        let heading_result = geographic_velocity[1].atan2(geographic_velocity[0]);
        println!("Heading: {}, {}", heading, heading_result);

        let flight_path_result = velocity_coordinates[0]
            .atan2((velocity_coordinates[0].powf(2.0) + velocity_coordinates[1].powf(2.0)).sqrt());
        println!("FP: {}, {}", flight_path, flight_path_result);

        let wut = relative_eq!(heading, heading_result, epsilon = 1e-3);
        assert_eq!(true, wut);
        let wut = relative_eq!(flight_path, flight_path_result, epsilon = 1e-3);
        assert_eq!(true, wut);
    }

    #[test]
    fn test_body_to_ned() {
        let gravity_in_body: Cartesian<f64, Body<f64>> = Cartesian::new(9.8_f64, 0.0_f64, 0.0_f64);

        let roll: f64 = 0.0_f64.to_radians();
        let pitch: f64 = -90.0_f64.to_radians();
        let yaw: f64 = 0.0_f64.to_radians();
        let result: DirectionCosineMatrix<f64, Body<f64>, NED<f64>> = body_to_ned(roll, pitch, yaw);

        let roll_result = result.m23().atan2(result.m33());
        let pitch_result = -result.m13().asin();
        let yaw_result = result.m12().atan2(result.m11());

        let grav_positive_in_ned = result * gravity_in_body;
        println!(
            "Gravity: body - {} | ned - {}",
            gravity_in_body, grav_positive_in_ned
        );
        assert_eq!(-gravity_in_body[0], grav_positive_in_ned[2]);
        assert_eq!(roll_result, roll);
        assert_eq!(pitch_result, pitch);
        assert_eq!(yaw_result, yaw);
    }

    #[test]
    fn test_itrf_to_ned() {
        // let d2r = 3.14159265359/180.0;
        // let latitude = 45.0 * d2r;
        // let longitude = 30.0 * d2r;

        // let x_e_s = Spherical::new(1.25, 0.0, longitude);
        // let x_e: Cartesian<f64, ITRF<f64>> = Cartesian::from(&x_e_s);
        // println!("x_e: {}", x_e.x());

        // let y_e = DirectionCosineMatrix::rotate_x(-90.0*d2r) * x_e;
        // println!("y_e: {}", y_e);

        // let n_origin_s = Spherical::new(1.25, latitude, longitude+15.0*d2r);
        // let n_origin = Cartesian::from(&n_origin_s);
        // println!("N Origin{n_origin}");

        // let d = -1.0 * n_origin.magnitude();
        // println!("D: {}", d+n_origin);

        // let z = Vector::new([0.0, 0.0, 1.25]);

        // let e = d.cross(&z);
        // let n = e.cross(&d);
        // println!("NED: {}, {}, {}", n_origin+n,n_origin+e,n_origin+d);
    }
}
