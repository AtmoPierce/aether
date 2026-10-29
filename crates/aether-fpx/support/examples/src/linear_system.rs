use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use aether_core::math::{Matrix, Vector};
use aether_core::real::Real;
use aether_fpx::ocp::fp4::Fp4E2M1;
use aether_fpx::ocp::fp6::Fp6E3M2;
use aether_fpx::ocp::fp8::{Fp8E4M3, Fp8E5M2};

const N: usize = 3;
const PRECISION_COUNT: usize = 8;

// A, b, and the exact solution are representable in every tested format.
const A_REFERENCE: [[f64; N]; N] = [[1.0, 0.5, 0.0], [0.5, 1.0, 0.5], [0.0, 0.5, 1.5]];
const B_REFERENCE: [f64; N] = [0.5, 0.0, 1.0];
const X_REFERENCE: [f64; N] = [1.0, -1.0, 1.0];

const OUTPUT_DIR: &str = "outputs/linear_system";
const CSV_PATH: &str = "outputs/linear_system/linear_system.csv";

const LABELS: [&str; PRECISION_COUNT] = [
    "FP4 E2M1", "FP6 E3M2", "FP8 E4M3", "FP8 E5M2", "f16", "f32", "f64", "f128",
];

const CSV_NAMES: [&str; PRECISION_COUNT] = [
    "fpx_fp4_e2m1",
    "fpx_fp6_e3m2",
    "fpx_fp8_e4m3",
    "fpx_fp8_e5m2",
    "fpx_f16",
    "fpx_f32",
    "fpx_f64",
    "fpx_f128",
];

#[derive(Clone, Copy, Debug)]
struct ResultRow {
    inverse_available: bool,
    inverse_solution: [f64; N],
    inverse_solution_error: f64,
    inverse_residual: f64,
    inverse_identity_defect: f64,
    lu_available: bool,
    lu_solution: [f64; N],
    lu_solution_error: f64,
    lu_residual: f64,
}

fn vector_norm(values: [f64; N]) -> f64 {
    f64::sqrt(values.into_iter().map(|value| value * value).sum())
}

fn solution_error(solution: [f64; N]) -> f64 {
    vector_norm(std::array::from_fn(|index| {
        solution[index] - X_REFERENCE[index]
    }))
}

fn residual(solution: [f64; N]) -> f64 {
    vector_norm(std::array::from_fn(|row| {
        let ax = (0..N)
            .map(|column| A_REFERENCE[row][column] * solution[column])
            .sum::<f64>();
        ax - B_REFERENCE[row]
    }))
}

fn identity_defect(inverse: [[f64; N]; N]) -> f64 {
    f64::sqrt(
        (0..N)
            .flat_map(|row| {
                (0..N).map(move |column| {
                    let value = (0..N)
                        .map(|inner| A_REFERENCE[row][inner] * inverse[inner][column])
                        .sum::<f64>();
                    let error = value - f64::from(row == column);
                    error * error
                })
            })
            .sum(),
    )
}

fn evaluate<T: Real + Copy>() -> ResultRow {
    let a = Matrix::<T, N, N>::new(std::array::from_fn(|row| {
        std::array::from_fn(|column| T::from_f64(A_REFERENCE[row][column]))
    }));
    let b = Vector::<T, N>::new(std::array::from_fn(|index| T::from_f64(B_REFERENCE[index])));

    let inverse = a.inverse_gauss_jordan();
    let (
        inverse_available,
        inverse_solution,
        inverse_solution_error,
        inverse_residual,
        inverse_identity_defect,
    ) = if let Some(inverse) = inverse {
        let solution = &inverse * &b;
        let widened_solution = std::array::from_fn(|index| solution[index].to_f64());
        let widened_inverse = std::array::from_fn(|row| {
            std::array::from_fn(|column| inverse[(row, column)].to_f64())
        });
        (
            true,
            widened_solution,
            solution_error(widened_solution),
            residual(widened_solution),
            identity_defect(widened_inverse),
        )
    } else {
        (false, [f64::NAN; N], f64::NAN, f64::NAN, f64::NAN)
    };

    let lu_solution = a.solve(&b);
    let (lu_available, lu_solution, lu_solution_error, lu_residual) =
        if let Some(solution) = lu_solution {
            let widened_solution = std::array::from_fn(|index| solution[index].to_f64());
            (
                true,
                widened_solution,
                solution_error(widened_solution),
                residual(widened_solution),
            )
        } else {
            (false, [f64::NAN; N], f64::NAN, f64::NAN)
        };

    ResultRow {
        inverse_available,
        inverse_solution,
        inverse_solution_error,
        inverse_residual,
        inverse_identity_defect,
        lu_available,
        lu_solution,
        lu_solution_error,
        lu_residual,
    }
}

fn evaluate_all() -> [ResultRow; PRECISION_COUNT] {
    [
        evaluate::<Fp4E2M1>(),
        evaluate::<Fp6E3M2>(),
        evaluate::<Fp8E4M3>(),
        evaluate::<Fp8E5M2>(),
        evaluate::<f16>(),
        evaluate::<f32>(),
        evaluate::<f64>(),
        evaluate::<f128>(),
    ]
}

fn write_csv(results: &[ResultRow; PRECISION_COUNT], path: &Path) -> std::io::Result<()> {
    let mut output = BufWriter::new(File::create(path)?);
    writeln!(
        output,
        "format,inverse_available,inverse_x0,inverse_x1,inverse_x2,inverse_solution_error,inverse_residual,inverse_identity_defect,lu_available,lu_x0,lu_x1,lu_x2,lu_solution_error,lu_residual"
    )?;
    for (index, result) in results.iter().enumerate() {
        writeln!(
            output,
            "{},{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            CSV_NAMES[index],
            result.inverse_available,
            result.inverse_solution[0],
            result.inverse_solution[1],
            result.inverse_solution[2],
            result.inverse_solution_error,
            result.inverse_residual,
            result.inverse_identity_defect,
            result.lu_available,
            result.lu_solution[0],
            result.lu_solution[1],
            result.lu_solution[2],
            result.lu_solution_error,
            result.lu_residual,
        )?;
    }
    output.flush()
}

pub(crate) fn run() -> Result<(), Box<dyn std::error::Error>> {
    fs::create_dir_all(OUTPUT_DIR)?;
    let results = evaluate_all();
    write_csv(&results, Path::new(CSV_PATH))?;

    println!("wrote {CSV_PATH}");
    println!("A = {:?}", A_REFERENCE);
    println!("b = {:?}; exact x = {:?}", B_REFERENCE, X_REFERENCE);
    for (label, result) in LABELS.iter().zip(results) {
        let lu_status = if result.lu_available {
            format!("residual {:.6e}", result.lu_residual)
        } else {
            "unavailable (singularity guard)".to_owned()
        };
        println!(
            "{label:>9}: inverse x = [{:>8.5}, {:>8.5}, {:>8.5}], error = {:.6e}, residual = {:.6e}, AA^-1-I = {:.6e}; LU {lu_status}",
            result.inverse_solution[0],
            result.inverse_solution[1],
            result.inverse_solution[2],
            result.inverse_solution_error,
            result.inverse_residual,
            result.inverse_identity_defect,
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_precision_solves_the_linear_system_with_finite_metrics() {
        for (label, result) in LABELS.iter().zip(evaluate_all()) {
            assert!(result.inverse_available, "{label} inversion failed");
            assert!(
                result.inverse_solution.into_iter().all(f64::is_finite)
                    && result.inverse_solution_error.is_finite()
                    && result.inverse_residual.is_finite()
                    && result.inverse_identity_defect.is_finite(),
                "{label} produced a non-finite result",
            );
            if result.lu_available {
                assert!(
                    result.lu_solution.into_iter().all(f64::is_finite)
                        && result.lu_solution_error.is_finite()
                        && result.lu_residual.is_finite(),
                    "{label} produced a non-finite LU result",
                );
            }
        }
    }
}
