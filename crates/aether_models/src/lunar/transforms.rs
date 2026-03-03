use crate::attitude::DirectionCosineMatrix;
use crate::reference_frame::ICRF;

use super::{
    constants::MOON_SIDEREAL_ROTATION_RATE_RAD_S,
    frames::{LIRF, LTRF},
};

pub fn icrf_to_lirf() -> DirectionCosineMatrix<f64, ICRF<f64>, LIRF> {
    DirectionCosineMatrix::new(
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    )
}

pub fn lirf_to_icrf() -> DirectionCosineMatrix<f64, LIRF, ICRF<f64>> {
    DirectionCosineMatrix::new(
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    )
}

pub fn lirf_to_ltrf(time: f64) -> DirectionCosineMatrix<f64, LIRF, LTRF> {
    let theta = MOON_SIDEREAL_ROTATION_RATE_RAD_S * time;
    DirectionCosineMatrix::new(
        theta.cos(),
        theta.sin(),
        0.0,
        -theta.sin(),
        theta.cos(),
        0.0,
        0.0,
        0.0,
        1.0,
    )
}

pub fn ltrf_to_lirf(time: f64) -> DirectionCosineMatrix<f64, LTRF, LIRF> {
    let theta = MOON_SIDEREAL_ROTATION_RATE_RAD_S * time;
    DirectionCosineMatrix::new(
        theta.cos(),
        -theta.sin(),
        0.0,
        theta.sin(),
        theta.cos(),
        0.0,
        0.0,
        0.0,
        1.0,
    )
}

pub fn icrf_to_ltrf(time: f64) -> DirectionCosineMatrix<f64, ICRF<f64>, LTRF> {
    lirf_to_ltrf(time) * icrf_to_lirf()
}

pub fn ltrf_to_icrf(time: f64) -> DirectionCosineMatrix<f64, LTRF, ICRF<f64>> {
    lirf_to_icrf() * ltrf_to_lirf(time)
}

pub fn icrf_to_ltrf_at_epoch(
    time_since_epoch: f64,
    icrf_to_lirf_at_epoch: DirectionCosineMatrix<f64, ICRF<f64>, LIRF>,
) -> DirectionCosineMatrix<f64, ICRF<f64>, LTRF> {
    lirf_to_ltrf(time_since_epoch) * icrf_to_lirf_at_epoch
}

pub fn ltrf_to_icrf_at_epoch(
    time_since_epoch: f64,
    icrf_to_lirf_at_epoch: DirectionCosineMatrix<f64, ICRF<f64>, LIRF>,
) -> DirectionCosineMatrix<f64, LTRF, ICRF<f64>> {
    icrf_to_lirf_at_epoch.transpose() * ltrf_to_lirf(time_since_epoch)
}

#[cfg(test)]
mod tests {
    use super::{
        icrf_to_lirf,
        icrf_to_ltrf,
        icrf_to_ltrf_at_epoch,
        lirf_to_icrf,
        lirf_to_ltrf,
        ltrf_to_icrf,
        ltrf_to_icrf_at_epoch,
        ltrf_to_lirf,
    };
    use crate::{
        attitude::DirectionCosineMatrix,
        coordinate::Cartesian,
        lunar::frames::LIRF,
        reference_frame::ICRF,
    };

    #[test]
    fn transforms_are_inverse() {
        let t = 12_345.0;
        let r_l2t = lirf_to_ltrf(t);
        let r_t2l = ltrf_to_lirf(t);

        let v_l = Cartesian::new(1.0, 2.0, 3.0);
        let v_t = r_l2t * v_l;
        let v_l_back = r_t2l * v_t;

        assert!((v_l_back.x() - 1.0).abs() < 1.0e-12);
        assert!((v_l_back.y() - 2.0).abs() < 1.0e-12);
        assert!((v_l_back.z() - 3.0).abs() < 1.0e-12);
    }

    #[test]
    fn zero_time_is_identity() {
        let r_l2t = lirf_to_ltrf(0.0);
        let r_t2l = ltrf_to_lirf(0.0);

        let v = Cartesian::new(4.0, -5.0, 6.0);
        let a = r_l2t * v;
        assert!((a.x() - 4.0).abs() < 1.0e-12);
        assert!((a.y() + 5.0).abs() < 1.0e-12);
        assert!((a.z() - 6.0).abs() < 1.0e-12);

        let v2 = Cartesian::new(-4.0, 5.0, -6.0);
        let c = r_t2l * v2;
        assert!((c.x() + 4.0).abs() < 1.0e-12);
        assert!((c.y() - 5.0).abs() < 1.0e-12);
        assert!((c.z() + 6.0).abs() < 1.0e-12);
    }

    #[test]
    fn icrf_lirf_are_aligned() {
        let r_i2l = icrf_to_lirf();
        let r_l2i = lirf_to_icrf();

        let v_i = Cartesian::<f64, ICRF<f64>>::new(7.0, -8.0, 9.0);
        let v_l = r_i2l * v_i;
        assert!((v_l.x() - 7.0).abs() < 1.0e-12);
        assert!((v_l.y() + 8.0).abs() < 1.0e-12);
        assert!((v_l.z() - 9.0).abs() < 1.0e-12);

        let v_i_back = r_l2i * v_l;
        assert!((v_i_back.x() - 7.0).abs() < 1.0e-12);
        assert!((v_i_back.y() + 8.0).abs() < 1.0e-12);
        assert!((v_i_back.z() - 9.0).abs() < 1.0e-12);
    }

    #[test]
    fn icrf_ltrf_round_trip() {
        let t = 9876.5;
        let r_i2t = icrf_to_ltrf(t);
        let r_t2i = ltrf_to_icrf(t);

        let v_i = Cartesian::<f64, ICRF<f64>>::new(3.0, 2.0, 1.0);
        let v_t = r_i2t * v_i;
        let v_i_back = r_t2i * v_t;

        assert!((v_i_back.x() - 3.0).abs() < 1.0e-12);
        assert!((v_i_back.y() - 2.0).abs() < 1.0e-12);
        assert!((v_i_back.z() - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn central_rotation_can_be_non_identity() {
        let theta: f64 = 0.25;
        let central = DirectionCosineMatrix::<f64, ICRF<f64>, LIRF>::new(
            theta.cos(), theta.sin(), 0.0,
            -theta.sin(), theta.cos(), 0.0,
            0.0, 0.0, 1.0,
        );

        let i2l = central;
        let l2i = central.transpose();

        let v_i = Cartesian::<f64, ICRF<f64>>::new(1.0, 0.0, 0.0);
        let v_l = i2l * v_i;
        assert!((v_l.x() - theta.cos()).abs() < 1.0e-12);
        assert!((v_l.y() + theta.sin()).abs() < 1.0e-12);

        let v_i_back = l2i * v_l;
        assert!((v_i_back.x() - 1.0).abs() < 1.0e-12);
        assert!(v_i_back.y().abs() < 1.0e-12);
        assert!(v_i_back.z().abs() < 1.0e-12);

        let t = 1000.0;
        let i2t = icrf_to_ltrf_at_epoch(t, central);
        let t2i = ltrf_to_icrf_at_epoch(t, central);
        let v_t = i2t * v_i;
        let v_i_round_trip = t2i * v_t;
        assert!((v_i_round_trip.x() - 1.0).abs() < 1.0e-12);
        assert!(v_i_round_trip.y().abs() < 1.0e-12);
        assert!(v_i_round_trip.z().abs() < 1.0e-12);
    }
}
