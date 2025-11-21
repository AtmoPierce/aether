use aether_core::attitude::DirectionCosineMatrix;
use aether_core::coordinate::Cartesian;
use aether_core::math::{Matrix, Vector};
use aether_core::reference_frame::ReferenceFrame;
use aether_core::real::Real;

pub fn look_at<T: Real, F: ReferenceFrame>(
    eye: &Cartesian<T, F>,
    center: &Cartesian<T, F>,
    up: &Cartesian<T, F>,
) -> Matrix<T, 4, 4> {
    let f = (center - eye).normalize().data;
    let r = up.data.cross(&f).normalize();
    let u = f.cross(&r);

    Matrix::new([
        [r[0], r[1], r[2], -r.dot(&eye.data)],
        [u[0], u[1], u[2], -u.dot(&eye.data)],
        [-f[0], -f[1], -f[2], f.dot(&eye.data)],
        [T::ZERO, T::ZERO, T::ZERO, T::ONE],
    ])
}

pub fn perspective<T: Real>(aspect: T, fov_y_rad: T, near: T, far: T) -> Matrix<T, 4, 4> {
    // assert!(aspect > T::ZERO, "Aspect ratio must be > 0");
    // assert!(fov_y_rad > T::ZERO, "Field of view must be > 0");
    // assert!(near != far, "Near and far planes must not be equal");

    let f = T::ONE / (fov_y_rad / T::from_f32(2.0)).tan();
    let nf = T::ONE / (near - far);

    Matrix {
        data: [
            [f / aspect, T::ZERO, T::ZERO, T::ZERO],
            [T::ZERO, f, T::ZERO, T::ZERO],
            [T::ZERO, T::ZERO, (far + near) * nf, -T::ONE],
            [
                T::ZERO,
                T::ZERO,
                (T::from_f32(2.0) * far * near) * nf,
                T::ZERO,
            ],
        ],
    }
}

pub fn look_at_perspective<T: Real, F: ReferenceFrame>(
    eye: &Cartesian<T, F>,
    center: &Cartesian<T, F>,
    up: &Cartesian<T, F>,
    aspect: T,
    fov_y_rad: T,
    near: T,
    far: T,
) -> (Matrix<T, 4, 4>, Matrix<T, 4, 4>) {
    let view = look_at(eye, center, up);
    let projection = perspective(aspect, fov_y_rad, near, far);
    (view, projection)
}

pub fn translate<T: Real>(mat: &Matrix<T, 4, 4>, trans: &Vector<T, 3>) -> Matrix<T, 4, 4> {
    let mut result = *mat;

    // Row-major translation: update the last column (row 0–2, col 3)
    result.data[0][3] = result.data[0][0] * trans[0]
        + result.data[0][1] * trans[1]
        + result.data[0][2] * trans[2]
        + result.data[0][3];

    result.data[1][3] = result.data[1][0] * trans[0]
        + result.data[1][1] * trans[1]
        + result.data[1][2] * trans[2]
        + result.data[1][3];

    result.data[2][3] = result.data[2][0] * trans[0]
        + result.data[2][1] * trans[1]
        + result.data[2][2] * trans[2]
        + result.data[2][3];

    result
}
