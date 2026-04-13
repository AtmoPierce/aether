use crate::models::terrestrial::wgs84::constants::{a, b, e2};
use crate::{
    attitude::DirectionCosineMatrix,
    coordinate::Cartesian,
    math::Vector,
    real::Real,
    reference_frame::{ICRF, ITRF},
};

pub fn icrf_to_itrf<T: Real>(
    time: T,
    rotational_velocity: Cartesian<T, ITRF<T>>,
) -> DirectionCosineMatrix<T, ICRF<T>, ITRF<T>> {
    return DirectionCosineMatrix::new(
        (rotational_velocity.z() * time).cos(),
        (rotational_velocity.z() * time).sin(),
        T::ZERO,
        (-rotational_velocity.z() * time).sin(),
        (rotational_velocity.z() * time).cos(),
        T::ZERO,
        T::ZERO,
        T::ZERO,
        T::ONE,
    );
}
pub fn itrf_to_icrf<T: Real>(
    time: T,
    rotational_velocity: Cartesian<T, ITRF<T>>,
) -> DirectionCosineMatrix<T, ITRF<T>, ICRF<T>> {
    return DirectionCosineMatrix::new(
        (rotational_velocity.z() * time).cos(),
        (-rotational_velocity.z() * time).sin(),
        T::ZERO,
        (rotational_velocity.z() * time).sin(),
        (rotational_velocity.z() * time).cos(),
        T::ZERO,
        T::ZERO,
        T::ZERO,
        T::ONE,
    );
}

pub fn geocentric_to_ecef<T: Real>(
    latitude: T,
    longitude: T,
    altitude: T,
) -> Cartesian<T, ITRF<T>> {
    let semi_major_axis = T::from_f64(a);
    let eccentricity_squared = T::from_f64(e2);
    let n = semi_major_axis / (T::ONE - eccentricity_squared * latitude.sin() * latitude.sin()).sqrt();
    let x = (n + altitude) * latitude.cos() * longitude.cos();
    let y = (n + altitude) * latitude.cos() * longitude.sin();
    let z = (n * (T::ONE - eccentricity_squared) + altitude) * latitude.sin();
    Cartesian::new(x, y, z)
}

pub fn ecef_to_geocentric_ferrari<T: Real>(x: T, y: T, z: T) -> Vector<T, 3> {
    let semi_major_axis = T::from_f64(a);
    let semi_minor_axis = T::from_f64(b);
    let a2 = semi_major_axis * semi_major_axis;
    let b2 = semi_minor_axis * semi_minor_axis;
    let f = (semi_major_axis - semi_minor_axis) / semi_major_axis;
    let two = T::ONE + T::ONE;
    let e_2 = two * f - f * f;
    let ep = (a2 / b2 - T::ONE).sqrt();
    let z2 = z * z;
    let r = (x * x + y * y).sqrt();
    let r2 = r * r;
    let e2_2 = a2 - b2;
    let ff = T::from_f64(54.0) * b2 * z2;
    let g = r2 + (T::ONE - e_2) * z2 - e_2 * e2_2;
    let c = e_2.powf(two) * ff * r2 / g.powf(T::from_f64(3.0));
    let s = (T::ONE + c + (c.powf(two) + two * c)).powf(T::ONE / T::from_f64(3.0));
    let pp = ff / (T::from_f64(3.0) * (s + T::ONE / s + T::ONE).powf(two) * g.powf(two));
    let q = (T::ONE + two * e_2.powf(two) * pp).sqrt();

    // NaNs are gross. -> This needs fixing. !TODO
    let mut r0_op = T::ONE / two * a2 * (T::ONE + T::ONE / q)
        - (pp * (T::ONE - e_2) * z2) / (q * (T::ONE + q))
        - T::ONE / two * pp * r2;
    if !((r0_op > T::ZERO) && (r0_op == r0_op)) {
        r0_op = T::ZERO;
    }

    let r0 = -pp * e_2 * r / (T::ONE + q) + r0_op.sqrt();
    let u = ((r - e_2 * r0).powf(two) + z2).sqrt();
    let v = ((r - e_2 * r0).powf(two) + (T::ONE - e_2) * z2).sqrt();
    let z0 = b2 * z / (semi_major_axis * v);

    let altitude = u * (T::ONE - (b2 / (semi_major_axis * v)));
    let latitude = (z + ep * ep * z0).atan2(r);
    let longitude = y.atan2(x);
    Vector::new([latitude, longitude, altitude])
}

pub fn ecef_to_geocentric<T: Real>(x: T, y: T, z: T) -> Vector<T, 3> {
    // Compute intermediate quantities
    let eps = T::EPSILON * T::from_f64(1.0e2);
    let eccentricity_squared = T::from_f64(e2);
    let semi_major_axis = T::from_f64(a);
    let rho2 = x * x + y * y;
    let mut dz = eccentricity_squared * z;
    let mut n = T::ZERO;

    // Iterative refine coordinate estimate
    let mut iter = 0;
    while iter < 101 {
        let zdz = z + dz;
        let nh = (rho2 + zdz * zdz).sqrt();
        let sinphi = zdz / nh;
        n = semi_major_axis / (T::ONE - eccentricity_squared * sinphi * sinphi).sqrt();
        let dz_new = n * eccentricity_squared * sinphi;

        // Check convergence requirement
        if (dz - dz_new).abs() < eps {
            break;
        }

        dz = dz_new;
        iter += 1;
    }

    if iter == 100 {
        let zdz = z + dz;
        let lon = y.atan2(x);
        let lat = zdz.atan2(rho2.sqrt());
        let alt = (rho2 + zdz * zdz).sqrt() - n;

        Vector::new([lat, lon, alt])
    } else {
        let zdz = z + dz;
        let lon = y.atan2(x);
        let lat = zdz.atan2(rho2.sqrt());
        let alt = (rho2 + zdz * zdz).sqrt() - n;

        Vector::new([lat, lon, alt])
    }
}
