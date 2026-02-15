use crate::attitude::{DirectionCosineMatrix, Euler};
use crate::math::{Matrix, Vector};
use crate::reference_frame::ReferenceFrame;
use crate::real::Real;
use crate::coordinate::Cartesian;

use core::marker::PhantomData;
use core::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Debug, Clone, Copy, PartialEq, Default)]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
pub struct Quaternion<T: Real, From: ReferenceFrame, To: ReferenceFrame> {
    pub data: Vector<T, 4>, // [w, i, j, k]
    _from: PhantomData<From>,
    _to: PhantomData<To>,
}

impl<T: Real, From: ReferenceFrame, To: ReferenceFrame> Quaternion<T, From, To> {
    #[inline]
    pub fn new(w: T, i: T, j: T, k: T) -> Self {
        Self {
            data: Vector { data: [w, i, j, k] },
            _from: PhantomData,
            _to: PhantomData,
        }
    }

    #[inline]
    pub fn identity() -> Self {
        Self::new(T::ONE, T::ZERO, T::ZERO, T::ZERO)
    }

    #[inline] pub fn w(&self) -> T { self.data[0] }
    #[inline] pub fn i(&self) -> T { self.data[1] }
    #[inline] pub fn j(&self) -> T { self.data[2] }
    #[inline] pub fn k(&self) -> T { self.data[3] }

    /// Build a quaternion from scalar + Euler components (x,y,z).
    #[inline]
    pub fn from_parts(w: T, v: Cartesian<T, From>) -> Self {
        Self::new(w, v.x(), v.y(), v.z())
    }

    #[inline]
    pub fn norm2(&self) -> T {
        self.data.dot(&self.data)
    }

    #[inline]
    pub fn norm(&self) -> T {
        self.norm2().sqrt()
    }

    #[inline]
    pub fn normalized(&self) -> Self {
        let n = self.norm();
        Self::new(
            self.w() / n,
            self.i() / n,
            self.j() / n,
            self.k() / n,
        )
    }

    #[inline]
    pub fn dot(&self, rhs: &Self) -> T {
        self.data.dot(&rhs.data)
    }

    /// Conjugate (frames unchanged)
    #[inline]
    pub fn conjugate(&self) -> Self {
        Self::new(
            self.w(),
            -self.i(),
            -self.j(),
            -self.k(),
        )
    }

    /// Inverse swaps frames
    #[inline]
    pub fn inverse(&self) -> Quaternion<T, To, From> {
        let n2 = self.norm2();
        Quaternion::<T, To, From>::new(
            self.w() / n2,
            -self.i() / n2,
            -self.j() / n2,
            -self.k() / n2,
        )
    }

    /// Passive frame transform: From -> To
    #[inline]
    pub fn rotate_vector(&self, v_from: Vector<T, 3>) -> Vector<T, 3> {
        // Passive rotation: q* ⊗ (0, v) ⊗ q
        let q = self.normalized();
        let qc = q.conjugate();

        let px = v_from[0];
        let py = v_from[1];
        let pz = v_from[2];

        // qc ⊗ p
        let aw = -(qc.i() * px + qc.j() * py + qc.k() * pz);
        let ax =  qc.w() * px + qc.j() * pz - qc.k() * py;
        let ay =  qc.w() * py + qc.k() * px - qc.i() * pz;
        let az =  qc.w() * pz + qc.i() * py - qc.j() * px;

        // (qc ⊗ p) ⊗ q
        Vector::new([
            aw * q.i() + ax * q.w() + ay * q.k() - az * q.j(),
            aw * q.j() - ax * q.k() + ay * q.w() + az * q.i(),
            aw * q.k() + ax * q.j() - ay * q.i() + az * q.w(),
        ])
    }

    /// Passive frame transform: Cartesian<From> -> Cartesian<To>
    #[inline]
    pub fn rotate_cartesian(
        &self,
        v_from: Cartesian<T, From>,
    ) -> Cartesian<T, To> {
        Cartesian {
            data: self.rotate_vector(v_from.data),
            _reference_frame: PhantomData,
        }
    }

    /// Convert to DCM (From -> To)
    #[inline]
    pub fn to_dcm(&self) -> DirectionCosineMatrix<T, From, To> {
        let q = self.normalized();
        let (w, x, y, z) = (q.w(), q.i(), q.j(), q.k());

        let one = T::ONE;
        let two = one + one;

        let ww = w * w;
        let xx = x * x;
        let yy = y * y;
        let zz = z * z;

        let wx = w * x;
        let wy = w * y;
        let wz = w * z;
        let xy = x * y;
        let xz = x * z;
        let yz = y * z;

        let m = Matrix::<T, 3, 3>::new([
            [ ww + xx - yy - zz, two * (xy - wz),     two * (xz + wy) ],
            [ two * (xy + wz),   ww - xx + yy - zz,  two * (yz - wx) ],
            [ two * (xz - wy),   two * (yz + wx),    ww - xx - yy + zz ],
        ]);

        DirectionCosineMatrix::from_matrix(m)
    }

    #[inline]
    fn mul_raw(&self, rhs: &Self) -> Self {
        let [w1, x1, y1, z1] = self.data.data;
        let [w2, x2, y2, z2] = rhs.data.data;

        Self::new(
            w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
            w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
        )
    }

    /// Integrate quaternion using angular velocity/acceleration (Taylor expansion).
    #[inline]
    pub fn integrate(self, ang_vel: Cartesian<T, From>, ang_acc: Cartesian<T, From>, dt: T) -> Self {
        let q = self.normalized();

        let one = T::ONE;
        let two = one + one;
        let three = two + one;
        let four = two + two;
        let nine = three * three;

        let half = one / two;
        let quarter = one / four;
        let ninth = one / nine;

        let ang_acc_quat = Quaternion::from_parts(one, ang_acc);
        let ang_vel_quat = Quaternion::from_parts(one, ang_vel);

        let q_dot = q.mul_raw(&ang_vel_quat) * half;
        let q_ddot = q.mul_raw(&ang_vel_quat).mul_raw(&ang_vel_quat) * quarter
            + q.mul_raw(&ang_acc_quat) * half;
        let q_dddot = q.mul_raw(&ang_vel_quat).mul_raw(&ang_vel_quat).mul_raw(&ang_vel_quat) * (one / (three + three))
            + q.mul_raw(&ang_acc_quat).mul_raw(&ang_vel_quat) * quarter
            + q.mul_raw(&ang_vel_quat).mul_raw(&ang_acc_quat) * half;

        let dt2 = dt * dt;
        let dt3 = dt2 * dt;

        let mut new_q = q + q_dot * dt;
        new_q = new_q + q_ddot * (quarter * dt2);
        new_q = new_q + q_dddot * (ninth * dt3);

        new_q.normalized()
    }
}

// ===== Quaternion composition (PASSIVE, matches DCM chaining) =====
//
// Convention: Quaternion<T, From, To> is a passive frame transform From -> To.
//
// Composition rule (rightmost acts first):
//   (Mid -> To) * (From -> Mid) = (From -> To)
//
// This matches DCM rule: R_ac = R_bc * R_ab
//

/// (Mid -> To) * (From -> Mid) = (From -> To)
impl<T: Real, From: ReferenceFrame, Mid: ReferenceFrame, To: ReferenceFrame>
    Mul<&Quaternion<T, From, Mid>> for &Quaternion<T, Mid, To>
{
    type Output = Quaternion<T, From, To>;

    #[inline]
    fn mul(self, rhs: &Quaternion<T, From, Mid>) -> Self::Output {
        // Hamilton product: self x rhs
        let [w1, x1, y1, z1] = self.data.data; // Mid -> To
        let [w2, x2, y2, z2] = rhs.data.data;  // From -> Mid

        Quaternion::new(
            w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
            w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
        )
    }
}

impl<T: Real, From: ReferenceFrame, Mid: ReferenceFrame, To: ReferenceFrame>
    Mul<Quaternion<T, From, Mid>> for Quaternion<T, Mid, To>
{
    type Output = Quaternion<T, From, To>;

    #[inline]
    fn mul(self, rhs: Quaternion<T, From, Mid>) -> Self::Output {
        (&self) * (&rhs)
    }
}

impl<T: Real, From: ReferenceFrame, To: ReferenceFrame>
    Mul<&Cartesian<T, From>> for &Quaternion<T, From, To>
{
    type Output = Cartesian<T, To>;

    fn mul(self, rhs: &Cartesian<T, From>) -> Self::Output {
        Cartesian {
            data: self.rotate_vector(rhs.data),
            _reference_frame: PhantomData,
        }
    }
}

impl<T: Real, From: ReferenceFrame, To: ReferenceFrame>
    Mul<Cartesian<T, From>> for Quaternion<T, From, To>
{
    type Output = Cartesian<T, To>;

    fn mul(self, rhs: Cartesian<T, From>) -> Self::Output {
        (&self) * (&rhs)
    }
}

impl<T: Real, From: ReferenceFrame, To: ReferenceFrame>
    Mul<Cartesian<T, From>> for &Quaternion<T, From, To>
{
    type Output = Cartesian<T, To>;

    fn mul(self, rhs: Cartesian<T, From>) -> Self::Output {
        self * (&rhs)
    }
}

impl<T: Real, From: ReferenceFrame, To: ReferenceFrame>
    Mul<&Cartesian<T, From>> for Quaternion<T, From, To>
{
    type Output = Cartesian<T, To>;

    fn mul(self, rhs: &Cartesian<T, From>) -> Self::Output {
        (&self) * rhs
    }
}

/// Scalar multiply
impl<T: Real, From: ReferenceFrame, To: ReferenceFrame> Mul<T>
    for Quaternion<T, From, To>
{
    type Output = Self;
    fn mul(self, rhs: T) -> Self {
        Self {
            data: self.data * rhs,
            _from: PhantomData,
            _to: PhantomData,
        }
    }
}

/// Addition
impl<T: Real, From: ReferenceFrame, To: ReferenceFrame> Add
    for Quaternion<T, From, To>
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            data: self.data + rhs.data,
            _from: PhantomData,
            _to: PhantomData,
        }
    }
}

/// Subtraction
impl<T: Real, From: ReferenceFrame, To: ReferenceFrame> Sub
    for Quaternion<T, From, To>
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self {
            data: self.data - rhs.data,
            _from: PhantomData,
            _to: PhantomData,
        }
    }
}

/// Negation
impl<T: Real, From: ReferenceFrame, To: ReferenceFrame> Neg
    for Quaternion<T, From, To>
{
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            data: -self.data,
            _from: PhantomData,
            _to: PhantomData,
        }
    }
}

/// Scalar division
impl<T: Real, From: ReferenceFrame, To: ReferenceFrame> Div<T>
    for Quaternion<T, From, To>
{
    type Output = Self;
    fn div(self, rhs: T) -> Self {
        Self {
            data: self.data / rhs,
            _from: PhantomData,
            _to: PhantomData,
        }
    }
}

impl<T: Real, A: ReferenceFrame, B: ReferenceFrame>
    TryFrom<&DirectionCosineMatrix<T, A, B>> for Quaternion<T, A, B>
{
    type Error = ();

    fn try_from(dcm: &DirectionCosineMatrix<T, A, B>) -> Result<Self, Self::Error> {
        // Norm Constraint
        // https://motoq.github.io/doc/tnotes/dcmq.pdf

        let m = &dcm.as_matrix().data;

        let one = T::ONE;
        let two = one + one;
        let four = two + two;
        let k = one / four;

        let c1 = one + m[0][0] + m[1][1] + m[2][2];
        let c2 = one + m[0][0] - m[1][1] - m[2][2];
        let c3 = one - m[0][0] + m[1][1] - m[2][2];
        let c4 = one - m[0][0] - m[1][1] + m[2][2];

        if c1 > k {
            let qs = (c1 / four).sqrt();
            let qs4 = qs * four;
            Ok(Quaternion::new(
                qs,
                (m[1][2] - m[2][1]) / qs4,
                (m[2][0] - m[0][2]) / qs4,
                (m[0][1] - m[1][0]) / qs4,
            ))
        } else if c2 > k {
            let qi = (c2 / four).sqrt();
            let qi4 = qi * four;
            Ok(Quaternion::new(
                (m[1][2] - m[2][1]) / qi4,
                qi,
                (m[0][1] + m[1][0]) / qi4,
                (m[2][0] + m[0][2]) / qi4,
            ))
        } else if c3 > k {
            let qj = (c3 / four).sqrt();
            let qj4 = qj * four;
            Ok(Quaternion::new(
                (m[2][0] - m[0][2]) / qj4,
                (m[0][1] + m[1][0]) / qj4,
                qj,
                (m[1][2] + m[2][1]) / qj4,
            ))
        } else if c4 > k {
            let qk = (c4 / four).sqrt();
            let qk4 = qk * four;
            Ok(Quaternion::new(
                (m[0][1] - m[1][0]) / qk4,
                (m[2][0] + m[0][2]) / qk4,
                (m[1][2] + m[2][1]) / qk4,
                qk,
            ))
        } else {
            Err(())
        }
    }
}

impl<T: Real, From: ReferenceFrame, To: ReferenceFrame>
    core::convert::From<&Euler<T, From, To>> for Quaternion<T, From, To>
{
    fn from(e: &Euler<T, From, To>) -> Self {
        let [roll, pitch, yaw] = e.data.data;

        let half = T::ONE / (T::ONE + T::ONE);
        let (cr, sr) = ((roll * half).cos(), (roll * half).sin());
        let (cp, sp) = ((pitch * half).cos(), (pitch * half).sin());
        let (cy, sy) = ((yaw * half).cos(), (yaw * half).sin());

        Quaternion::new(
            cr * cp * cy + sr * sp * sy,
            sr * cp * cy - cr * sp * sy,
            cr * sp * cy + sr * cp * sy,
            cr * cp * sy - sr * sp * cy,
        )
    }
}

#[cfg(feature = "std")]
impl<T, From, To> core::fmt::Display for Quaternion<T, From, To>
where
    T: Real + core::fmt::Display,
    From: ReferenceFrame,
    To: ReferenceFrame,
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "q=[{}, {}, {}, {}]",
            self.w(),
            self.i(),
            self.j(),
            self.k(),
        )
    }
}




#[cfg(test)]
#[path = "tests/quaternion_tests.rs"]
mod quaternion_tests;
