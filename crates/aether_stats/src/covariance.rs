use aether_core::{
    math::{Matrix, Vector},
    real::Real,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CovarianceError {
    EmptySamples,
    InsufficientSamples,
    SampleCountMismatch,
    ZeroVarianceAxis,
}

pub fn mean<T: Real, const N: usize>(samples: &[Vector<T, N>]) -> Result<Vector<T, N>, CovarianceError> {
    if samples.is_empty() {
        return Err(CovarianceError::EmptySamples);
    }

    let mut mean = Vector::new([T::ZERO; N]);
    for sample in samples {
        mean = mean + *sample;
    }

    Ok(mean / T::from_usize(samples.len()))
}

pub fn sample_covariance<T: Real, const N: usize>(
    samples: &[Vector<T, N>],
) -> Result<Matrix<T, N, N>, CovarianceError> {
    if samples.is_empty() {
        return Err(CovarianceError::EmptySamples);
    }
    if samples.len() < 2 {
        return Err(CovarianceError::InsufficientSamples);
    }

    let mean = mean(samples)?;
    let scale = T::ONE / T::from_usize(samples.len() - 1);
    let mut covariance = Matrix::<T, N, N>::zeros();

    for sample in samples {
        let centered = *sample - mean;
        for row in 0..N {
            for col in 0..N {
                covariance.data[row][col] = covariance.data[row][col] + centered[row] * centered[col];
            }
        }
    }

    Ok(covariance * scale)
}

pub fn sample_cross_covariance<T: Real, const N: usize, const M: usize>(
    xs: &[Vector<T, N>],
    ys: &[Vector<T, M>],
) -> Result<Matrix<T, N, M>, CovarianceError> {
    if xs.len() != ys.len() {
        return Err(CovarianceError::SampleCountMismatch);
    }
    if xs.is_empty() {
        return Err(CovarianceError::EmptySamples);
    }
    if xs.len() < 2 {
        return Err(CovarianceError::InsufficientSamples);
    }

    let mean_x = mean(xs)?;
    let mean_y = mean(ys)?;
    let scale = T::ONE / T::from_usize(xs.len() - 1);
    let mut covariance = Matrix::<T, N, M>::zeros();

    for (x, y) in xs.iter().zip(ys.iter()) {
        let centered_x = *x - mean_x;
        let centered_y = *y - mean_y;
        for row in 0..N {
            for col in 0..M {
                covariance.data[row][col] = covariance.data[row][col] + centered_x[row] * centered_y[col];
            }
        }
    }

    Ok(covariance * scale)
}

pub fn standard_deviations_from_covariance<T: Real, const N: usize>(
    covariance: &Matrix<T, N, N>,
) -> Result<Vector<T, N>, CovarianceError> {
    let mut standard_deviations = [T::ZERO; N];
    for (index, standard_deviation) in standard_deviations.iter_mut().enumerate() {
        let variance = covariance.data[index][index];
        if variance <= T::EPSILON {
            return Err(CovarianceError::ZeroVarianceAxis);
        }
        *standard_deviation = variance.sqrt();
    }

    Ok(Vector::new(standard_deviations))
}

pub fn correlation_from_covariance<T: Real, const N: usize>(
    covariance: &Matrix<T, N, N>,
) -> Result<Matrix<T, N, N>, CovarianceError> {
    let standard_deviations = standard_deviations_from_covariance(covariance)?;
    let mut correlation = Matrix::<T, N, N>::zeros();

    for row in 0..N {
        for col in 0..N {
            correlation.data[row][col] = covariance.data[row][col]
                / (standard_deviations[row] * standard_deviations[col]);
        }
    }

    Ok(correlation)
}

pub fn covariance_from_correlation<T: Real, const N: usize>(
    correlation: &Matrix<T, N, N>,
    standard_deviations: &Vector<T, N>,
) -> Matrix<T, N, N> {
    let mut covariance = Matrix::<T, N, N>::zeros();

    for row in 0..N {
        for col in 0..N {
            covariance.data[row][col] = correlation.data[row][col]
                * standard_deviations[row]
                * standard_deviations[col];
        }
    }

    covariance
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn computes_sample_covariance() {
        let samples = [
            Vector::new([1.0_f64, 2.0_f64]),
            Vector::new([2.0_f64, 1.0_f64]),
            Vector::new([3.0_f64, 0.0_f64]),
        ];

        let covariance = sample_covariance(&samples).expect("covariance should exist");

        assert!((covariance.data[0][0] - 1.0).abs() < 1.0e-12);
        assert!((covariance.data[1][1] - 1.0).abs() < 1.0e-12);
        assert!((covariance.data[0][1] + 1.0).abs() < 1.0e-12);
        assert!((covariance.data[1][0] + 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn computes_sample_cross_covariance() {
        let xs = [
            Vector::new([1.0_f64]),
            Vector::new([2.0_f64]),
            Vector::new([3.0_f64]),
        ];
        let ys = [
            Vector::new([2.0_f64]),
            Vector::new([4.0_f64]),
            Vector::new([6.0_f64]),
        ];

        let covariance = sample_cross_covariance(&xs, &ys).expect("cross covariance should exist");

        assert!((covariance.data[0][0] - 2.0).abs() < 1.0e-12);
    }

    #[test]
    fn correlation_and_covariance_round_trip() {
        let covariance = Matrix::new([
            [4.0_f64, 3.0_f64],
            [3.0_f64, 9.0_f64],
        ]);

        let standard_deviations = standard_deviations_from_covariance(&covariance)
            .expect("standard deviations should exist");
        let correlation = correlation_from_covariance(&covariance)
            .expect("correlation should exist");
        let rebuilt_covariance = covariance_from_correlation(&correlation, &standard_deviations);

        for row in 0..2 {
            for col in 0..2 {
                assert!((rebuilt_covariance.data[row][col] - covariance.data[row][col]).abs() < 1.0e-12);
            }
        }
    }
}
