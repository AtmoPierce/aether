use crate::attitude::DirectionCosineMatrix;
use crate::coordinate::Cartesian;
use crate::reference_frame::{Body, ITRF, ICRF, NED};
use crate::math::Matrix;
use crate::matrix;

pub fn angular_rate_dcm(roll: f64, pitch: f64, yaw: f64)-> Matrix<f64, 3, 3>{
    let rotation_matrix = matrix![
        0.0, -yaw, pitch;
        yaw, 0.0, -roll;
        -pitch, roll, 0.0;
    ];
    return rotation_matrix;
}

pub fn aircraft_wind_to_body(angle_of_attack: f64, side_slip: f64)->Matrix<f64, 3, 3>{
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

pub fn aeroballistic_wind_to_body(angle_of_attack: f64, aerodynamic_roll_angle: f64)->Matrix<f64, 3, 3>{
    let t_ab = matrix![
        angle_of_attack.cos(),  angle_of_attack.sin()*aerodynamic_roll_angle.sin(), angle_of_attack.sin()*aerodynamic_roll_angle.cos();
        0.0,                    aerodynamic_roll_angle.cos(),                       -aerodynamic_roll_angle.sin();
        -angle_of_attack.sin(), angle_of_attack.cos()*aerodynamic_roll_angle.sin(), angle_of_attack.cos()*aerodynamic_roll_angle.cos();
    ];
    return t_ab;
}

pub fn flight_path_to_geographic(heading_angle: f64, flight_path_angle: f64)->Matrix<f64, 3, 3>{
    let rotation_matrix = matrix![
        flight_path_angle.cos()*heading_angle.cos(),    flight_path_angle.cos()*heading_angle.sin(),    -flight_path_angle.sin();
        -heading_angle.sin(),                           heading_angle.cos(),                            0.0;
        flight_path_angle.sin()*heading_angle.cos(),    flight_path_angle.sin()*heading_angle.sin(),    flight_path_angle.cos();
    ];
    return rotation_matrix;
}

pub fn body_to_ned(roll: f64, pitch: f64, yaw: f64)->DirectionCosineMatrix<f64, Body<f64>, NED<f64>>{
    let phi = roll;
    let tht = pitch;
    let psi = yaw;
    let rotation_matrix = matrix![
        psi.cos()*tht.cos(),                                  psi.sin()*tht.cos(),                                     -tht.sin();
        psi.cos()*tht.sin()*phi.sin()-psi.sin()*phi.cos(),    psi.sin()*tht.sin()*phi.sin()+psi.cos()*phi.cos(),  tht.cos()*phi.sin();
        psi.cos()*tht.sin()*phi.cos()+psi.sin()*phi.sin(),    psi.sin()*tht.sin()*phi.cos()-psi.cos()*phi.sin(),  tht.cos()*phi.cos();
    ];
    let dcm: DirectionCosineMatrix<f64, Body<f64>, NED<f64>> = DirectionCosineMatrix::from_matrix(rotation_matrix);
    return dcm;
}

pub fn ecef_to_ned(latitude: f64, longitude: f64)->DirectionCosineMatrix<f64, ITRF<f64>, NED<f64>>{
    let rotation_matrix = matrix![
        -latitude.sin()*longitude.cos(),    -latitude.sin()*longitude.sin(),    latitude.cos();
        -longitude.sin(),                   longitude.cos(),                    0.0;
        -latitude.cos()*longitude.cos(),    -latitude.cos()*longitude.sin(),    -latitude.sin();
    ];
    let dcm = DirectionCosineMatrix::from_matrix(rotation_matrix);
    return dcm;
}

pub fn geocentric_to_ecef(latitude: f64, longitude: f64, altitude: f64)->Cartesian<f64, ITRF<f64>>{
    return wgs84::transforms::transforms::geocentric_to_ecef(latitude, longitude, altitude);
}

pub fn ecef_to_geocentric_ferrari(x: f64, y: f64, z: f64) -> Cartesian<f64, ITRF<f64>>{
    return wgs84::transforms::transforms::ecef_to_geocentric_ferrari(x, y, z);
}

pub fn ecef_to_geocentric(x: f64, y: f64, z: f64) -> Cartesian<f64, ITRF<f64>> {
    return wgs84::transforms::transforms::ecef_to_geocentric(x, y, z);
}

pub fn eci_to_ecef(time: f64, rotational_velocity: Cartesian<f64, ITRF<f64>>)->DirectionCosineMatrix<f64, ICRF<f64>, ITRF<f64>>{
    return DirectionCosineMatrix::new(
        (rotational_velocity.z()*time).cos(), (rotational_velocity.z()*time).sin(), 0.0,
        (-rotational_velocity.z()*time).sin(), (rotational_velocity.z()*time).cos(), 0.0,
        0.0 ,0.0 ,1.0
    );
}
pub fn ecef_to_eci(time: f64, rotational_velocity: Cartesian<f64, ITRF<f64>>)->DirectionCosineMatrix<f64, ITRF<f64>, ICRF<f64>>{
    return DirectionCosineMatrix::new(
        (rotational_velocity.z()*time).cos(), (-rotational_velocity.z()*time).sin(), 0.0,
        (rotational_velocity.z()*time).sin(), (rotational_velocity.z()*time).cos(), 0.0,
        0.0 ,0.0 ,1.0
    );
}

#[cfg(test)]
mod tests{
    use super::*;
    use approx::relative_eq;
    use wgs84::transforms::transforms::ecef_to_geocentric_ferrari;
    #[test]
    fn test_aircraft_wind_to_body(){
        let wind = Vector::new([1.0, 1.0, 1.0]);
        let wind_magnitude: f64 = wind.magnitude();

        let angle_of_attack = wind.z.atan2(wind.x);
        let aircraft_side_slip = (wind.y/wind_magnitude).asin();

        let result = aircraft_wind_to_body(angle_of_attack, aircraft_side_slip);
        let aoa_result = result.m12.atan2(result.m11);
        let sideslip_result = result.m13.asin();

        let wut = relative_eq!(angle_of_attack, aoa_result, epsilon = 1e-3);
        assert_eq!(true, wut);
        let wut = relative_eq!(aircraft_side_slip, sideslip_result, epsilon = 1e-3);
        assert_eq!(true, wut);
    }


    #[test]
    fn test_aeroballistic_wind_to_body(){
        let wind = Vector::new([1.0, 1.0, 1.0]);
        let wind_magnitude: f64 = wind.magnitude();

        let angle_of_attack = wind.z.atan2(wind_magnitude);
        let roll_angle = wind.y.atan2(wind.z);

        let result = aeroballistic_wind_to_body(angle_of_attack, roll_angle);
        let r1 = -result.m13.asin();
        let r2 = result.m12.atan2(result.m11);
        let r3 = result.m23.atan2(result.m33);



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
    fn test_flight_path_to_geographic(){
        let geographic_velocity = Vector::new([100.0, 100.0, 0.0]);

        let heading = 45.0_f64.to_radians();
        let flight_path = 45.0_f64.to_radians();

        let result = flight_path_to_geographic(heading, flight_path);
        let velocity_coordinates = result * geographic_velocity;

        let heading_result = geographic_velocity.y.atan2(geographic_velocity.x);
        println!("Heading: {}, {}", heading, heading_result);

        let flight_path_result = velocity_coordinates.z.atan2((velocity_coordinates.x.powf(2.0) + velocity_coordinates.y.powf(2.0)).sqrt());
        println!("FP: {}, {}", flight_path, flight_path_result);

        let wut = relative_eq!(heading, heading_result, epsilon = 1e-3);
        assert_eq!(true, wut);
        let wut = relative_eq!(flight_path, flight_path_result, epsilon = 1e-3);
        assert_eq!(true, wut);
    }

    #[test]
    fn test_body_to_ned(){
        let gravity_in_body = Matrixx1::new([-9.8, 0.0, 0.0]);

        let roll: f64 = 0.0_f64.to_radians();
        let pitch: f64 = 90.0_f64.to_radians();
        let yaw: f64 = 0.0_f64.to_radians();
        let result = body_to_ned(roll, pitch, yaw);

        let roll_result = result.m23.atan2(result.m33);
        let pitch_result = -result.m13.asin();
        let yaw_result = result.m12.atan2(result.m11);

        let grav_positive_in_ned =  result * gravity_in_body;
        assert_eq!(gravity_in_body.x, grav_positive_in_ned.z);
        assert_eq!(roll_result, roll);
        assert_eq!(pitch_result, pitch);
        assert_eq!(yaw_result, yaw);
    }

    #[test]
    fn test_ecef_to_ned(){
        let d2r = 3.14159265359/180.0;
        let latitude = 45.0 * d2r;
        let longitude = 30.0 * d2r;

        let x_e = spherical_to_rectangular(1.25, 0.0, longitude);
        println!("x_e: {}", x_e);

        let y_e = r3(-90.0*d2r) * x_e;
        println!("y_e: {}", y_e);

        let n_origin = spherical_to_rectangular(1.25, latitude, longitude+15.0*d2r);
        println!("N Origin{n_origin}");

        let d = -1.0 * n_origin.normalize();
        println!("D: {}", d+n_origin);            

        let z = Matrixx1::new(0.0, 0.0, 1.25);

        let e = d.cross(&z);
        let n = e.cross(&d);
        println!("NED: {}, {}, {}", n_origin+n,n_origin+e,n_origin+d);

    }
}