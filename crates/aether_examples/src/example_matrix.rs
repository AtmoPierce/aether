use aether_core::math::*;

fn main() {
    // Define a 3×3 rotation matrix (example: 45° rotation about Z-axis)
    let theta = std::f64::consts::FRAC_PI_4;
    let rot_z = Matrix::<f64, 3, 3>::new([
        [ theta.cos(), -theta.sin(), 0.0 ],
        [ theta.sin(),  theta.cos(), 0.0 ],
        [ 0.0,          0.0,         1.0 ],
    ]);

    // Define a 3×1 vector (position, velocity, or generic state)
    let v = Vector::<f64, 3>::new([1.0, 0.0, 0.0]);

    // Apply the rotation
    let v_rot = rot_z * v;

    println!("Rotated vector: {:?}", v_rot);
}