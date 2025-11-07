use num_traits::Float;
pub struct NewtonRaphson<T: Float> {
    tolerance: T,
    max_iterations: usize,
}

impl <T: Float> NewtonRaphson<T> {
    pub fn new(tolerance: T, max_iterations: usize) -> Self {
        Self {
            tolerance,
            max_iterations,
        }
    }

    pub fn solve<F, DF>(&self, mut guess: T, f: F, df: DF) -> Result<T, &'static str>
    where
        F: Fn(T) -> T,
        DF: Fn(T) -> T{
        for _ in 0..self.max_iterations {
            let f_val = f(guess);
            let df_val = df(guess);

            if df_val.is_zero() {
                return Err("Derivative is zero; no solution found.");
            }

            let next_guess = guess - f_val / df_val;

            if (next_guess - guess).abs() < self.tolerance {
                return Ok(next_guess);
            }

            guess = next_guess;
        }
        Err("Maximum iterations reached; no solution found.")
    }
}