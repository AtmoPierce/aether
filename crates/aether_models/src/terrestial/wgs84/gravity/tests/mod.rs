use crate::models::terrestrial::wgs84::gravity::gravity_rectangular;
#[test]
fn gravity_test() {
    let x = 6378137.0;
    let y = 0.0;
    let z = 0.0;
    let g_r = gravity_rectangular(x, y, z);

    println!("Gravity: {}", g_r);
}
