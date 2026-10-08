pub mod xor_shift;
pub use xor_shift::*;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reexport_available() {
        let mut r = XorShift64Star::new(1);
        let _ = r.next_u64();
    }
}
