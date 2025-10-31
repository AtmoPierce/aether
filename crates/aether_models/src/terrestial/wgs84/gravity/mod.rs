use crate::models::terrestrial::wgs84::{
    constants::{a, b, e2, g_e, k, w, E, GM},
    transforms::ecef_to_geocentric_ferrari,
};
use crate::{
    coordinate::{coordinate::Coordinate, Cartesian},
    math::{Matrix, Vector},
    reference_frame::ITRF,
};
use core::marker::PhantomData;

pub fn gravity_normal(latitude: f64) -> f64 {
    return g_e * (1.0 + k * latitude.sin().sin()) / (1.0 - e2 * latitude.sin().sin());
}

fn ellipsoidal_beta(x: f64, y: f64, z: f64, u: f64) -> f64 {
    let beta = ((z * (u.powf(2.0) + E.powf(2.0)).sqrt())
        / (u * (x.powf(2.0) + y.powf(2.0)).sqrt()))
    .atan();
    return beta;
}

fn ellipsoidal_height(x: f64, y: f64, z: f64) -> f64 {
    let t_share = x.powf(2.0) + y.powf(2.0) + z.powf(2.0) - E.powf(2.0);
    let t1 = 0.5 * t_share;
    let t2 = 1.0 + (1.0 + (4.0 * E.powf(2.0) * z.powf(2.0)) / (t_share.powf(2.0))).sqrt();
    let u = (t1 * t2).sqrt();
    return u;
}

fn ellipsoidal_w(u: f64, beta: f64) -> f64 {
    let ellip_w =
        ((u.powf(2.0) + E.powf(2.0) * (beta.sin().sin())) / (u.powf(2.0) + E.powf(2.0))).sqrt();
    return ellip_w;
}

pub fn gravity_ellipsoidal(x: f64, y: f64, z: f64) -> Vector<f64, 3> {
    let u = ellipsoidal_height(x, y, z);
    let beta = ellipsoidal_beta(x, y, z, u);
    let ellip_w = ellipsoidal_w(u, beta);

    let q = 0.5 * ((1.0 + 3.0 * u.powf(2.0) / E.powf(2.0)) * (E / u).atan() - 3.0 * u / E);
    let q_0 = 0.5 * ((1.0 + 3.0 * b.powf(2.0) / E.powf(2.0)) * (E / b).atan() - 3.0 * b / E);
    let qprime = 3.0 * (1.0 + u.powf(2.0) / E.powf(2.0)) * (1.0 - u / E * (E / u).atan()) - 1.0;

    // let t1 = -1.0 / ellip_w;
    // let t2 = (GM/(u.powf(2.0)+E.powf(2.0)));
    // let t3 = (w.powf(2.0)*a.powf(2.0)*E)/(u.powf(2.0)+E.powf(2.0));
    // let t4 = (qprime/q_0)*(0.5 * beta.sin().sin()-1.0/6.0);

    let g_u = -1.0 / ellip_w
        * ((GM / (u.powf(2.0) + E.powf(2.0)))
            + (w.powf(2.0) * a.powf(2.0) * E) / (u.powf(2.0) + E.powf(2.0))
                * (qprime / q_0)
                * (0.5 * beta.sin().sin() - 1.0 / 6.0))
        + 1.0 / ellip_w * w.powf(2.0) * u * beta.cos().cos();
    let g_b = 1.0 / ellip_w * (w.powf(2.0) * a.powf(2.0) / (u.powf(2.0) + E.powf(2.0))) * q / q_0
        * beta.sin()
        * beta.cos()
        - 1.0 / ellip_w
            * w.powf(2.0)
            * (u.powf(2.0) + E.powf(2.0)).sqrt()
            * beta.sin()
            * beta.cos();
    let g_l = 0.0;
    return Vector::new([g_u, g_b, g_l]);
}
pub fn gravity_rectangular(x: f64, y: f64, z: f64) -> Cartesian<f64, ITRF<f64>> {
    let u = ellipsoidal_height(x, y, z);
    let beta = ellipsoidal_beta(x, y, z, u);
    let ellip_w = ellipsoidal_w(u, beta);
    let lla = ecef_to_geocentric_ferrari(x, y, z);
    let g_ellip = gravity_ellipsoidal(x, y, z);

    let long = lla[1];

    let low = u / (ellip_w * (u.powf(2.0) + E.powf(2.0)).sqrt()) * beta.cos();

    let m11 = low * lla[1].cos();
    let m12 = -1.0 / ellip_w * beta.sin() * long.cos();
    let m13 = -(long.sin());
    let m21 = low * long.sin();
    let m22 = -1.0 / ellip_w * beta.sin() * long.sin();
    let m23 = long.cos();
    let m31 = 1.0 / ellip_w * beta.sin();
    let m32 = low;
    let m33 = 0.0;
    let r1 = Matrix::new([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]]);
    let g_rect = r1 * g_ellip;
    return Cartesian {
        data: g_rect,
        _reference_frame: PhantomData::<ITRF<f64>>,
    };
}

#[cfg(test)]
#[path = "tests/mod.rs"]
mod tests;
