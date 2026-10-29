use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use aether_core::math::Vector;
use aether_core::real::Real;
use aether_ml::classification::perceptron::Perceptron;
use aether_fpx::ocp::fp4::Fp4E2M1;
use aether_fpx::ocp::fp6::Fp6E3M2;
use aether_fpx::ocp::fp8::{Fp8E4M3, Fp8E5M2};

const FEATURE_COUNT: usize = 2;
const MAX_EPOCHS: usize = 64;
const LEARNING_RATE: f64 = 0.5;
const TRAIN_GRID_STEPS: usize = 16;
const TRAIN_GRID_SPACING: f64 = 0.25;
const TRAIN_MARGIN: f64 = 0.08;
const VALIDATION_GRID_STEPS: usize = 80;
const VALIDATION_GRID_SPACING: f64 = 0.05;

const OUTPUT_DIR: &str = "outputs/perceptron";
const CSV_PATH: &str = "outputs/perceptron/perceptron.csv";

const LABELS: [&str; 8] = [
    "FPX fp4 e2m1",
    "FPX fp6 e3m2",
    "FPX fp8 e4m3",
    "FPX fp8 e5m2",
    "FPX f16",
    "FPX f32",
    "FPX f64",
    "FPX f128",
];

const CSV_NAMES: [&str; 8] = [
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
struct TrainingSample {
    features: [f64; FEATURE_COUNT],
    label: i8,
}

#[derive(Clone, Copy, Debug)]
struct PerceptronResult {
    epochs: usize,
    stable: bool,
    converged: bool,
    training_accuracy: f64,
    validation_accuracy: f64,
    weights: [f64; FEATURE_COUNT],
    bias: f64,
    cast_weights: [f64; FEATURE_COUNT],
    cast_bias: f64,
    cast_validation_accuracy: f64,
    cast_reference_agreement: f64,
    weight_only_reference_agreement: f64,
}

fn analytic_score(x: f64, y: f64) -> f64 {
    0.73 * x - 0.41 * y + 0.17
}

fn analytic_label(x: f64, y: f64) -> i8 {
    if analytic_score(x, y) >= 0.0 {
        1
    } else {
        -1
    }
}

fn training_samples() -> Vec<TrainingSample> {
    let mut samples = Vec::new();
    for x_index in 0..=TRAIN_GRID_STEPS {
        let x = -2.0 + x_index as f64 * TRAIN_GRID_SPACING;
        for y_index in 0..=TRAIN_GRID_STEPS {
            let y = -2.0 + y_index as f64 * TRAIN_GRID_SPACING;
            let score = analytic_score(x, y);
            if score.abs() >= TRAIN_MARGIN {
                samples.push(TrainingSample {
                    features: [x, y],
                    label: if score >= 0.0 { 1 } else { -1 },
                });
            }
        }
    }
    samples
}

fn accuracy<T: Real + Copy>(
    model: &Perceptron<T, FEATURE_COUNT>,
    samples: &[TrainingSample],
) -> f64 {
    let correct = samples
        .iter()
        .filter(|sample| {
            let input = Vector::new([
                T::from_f64(sample.features[0]),
                T::from_f64(sample.features[1]),
            ]);
            model.predict(&input) == sample.label
        })
        .count();
    correct as f64 / samples.len() as f64
}

fn validation_accuracy<T: Real + Copy>(model: &Perceptron<T, FEATURE_COUNT>) -> f64 {
    let mut correct = 0usize;
    let mut total = 0usize;
    for x_index in 0..=VALIDATION_GRID_STEPS {
        let x = -2.0 + x_index as f64 * VALIDATION_GRID_SPACING;
        for y_index in 0..=VALIDATION_GRID_STEPS {
            let y = -2.0 + y_index as f64 * VALIDATION_GRID_SPACING;
            let input = Vector::new([T::from_f64(x), T::from_f64(y)]);
            correct += usize::from(model.predict(&input) == analytic_label(x, y));
            total += 1;
        }
    }
    correct as f64 / total as f64
}

fn validation_agreement<T: Real + Copy>(
    model: &Perceptron<T, FEATURE_COUNT>,
    reference: &Perceptron<f64, FEATURE_COUNT>,
) -> f64 {
    let mut matches = 0usize;
    let mut total = 0usize;
    for x_index in 0..=VALIDATION_GRID_STEPS {
        let x = -2.0 + x_index as f64 * VALIDATION_GRID_SPACING;
        for y_index in 0..=VALIDATION_GRID_STEPS {
            let y = -2.0 + y_index as f64 * VALIDATION_GRID_SPACING;
            let input = Vector::new([T::from_f64(x), T::from_f64(y)]);
            let reference_input = Vector::new([x, y]);
            matches += usize::from(model.predict(&input) == reference.predict(&reference_input));
            total += 1;
        }
    }
    matches as f64 / total as f64
}

fn train_model<T: Real + Copy>(
    samples: &[TrainingSample],
) -> (Perceptron<T, FEATURE_COUNT>, usize, bool) {
    let features = samples
        .iter()
        .map(|sample| {
            Vector::new([
                T::from_f64(sample.features[0]),
                T::from_f64(sample.features[1]),
            ])
        })
        .collect::<Vec<_>>();
    let labels = samples
        .iter()
        .map(|sample| sample.label)
        .collect::<Vec<_>>();

    let mut model = Perceptron::<T, FEATURE_COUNT>::new();
    model.lr = T::from_f64(LEARNING_RATE);
    let mut stable = false;
    let mut epochs = MAX_EPOCHS;

    for epoch in 1..=MAX_EPOCHS {
        let previous_weights = model.weights;
        let previous_bias = model.bias;
        model.train_epoch(&features, &labels);
        if model.weights == previous_weights && model.bias == previous_bias {
            stable = true;
            epochs = epoch;
            break;
        }
    }

    (model, epochs, stable)
}

fn cast_reference_model<T: Real + Copy>(
    reference: &Perceptron<f64, FEATURE_COUNT>,
) -> Perceptron<T, FEATURE_COUNT> {
    let mut model = Perceptron::<T, FEATURE_COUNT>::new();
    model.weights = Vector::new([
        T::from_f64(reference.weights[0]),
        T::from_f64(reference.weights[1]),
    ]);
    model.bias = T::from_f64(reference.bias);
    model.lr = T::from_f64(reference.lr);
    model
}

fn summarize<T: Real + Copy>(
    model: &Perceptron<T, FEATURE_COUNT>,
    samples: &[TrainingSample],
    epochs: usize,
    stable: bool,
    reference: &Perceptron<f64, FEATURE_COUNT>,
) -> PerceptronResult {
    let training_accuracy = accuracy(&model, samples);
    let cast_model = cast_reference_model::<T>(reference);
    let mut weight_only_model = Perceptron::<f64, FEATURE_COUNT>::new();
    weight_only_model.weights = Vector::new([
        cast_model.weights[0].to_f64(),
        cast_model.weights[1].to_f64(),
    ]);
    weight_only_model.bias = cast_model.bias.to_f64();
    weight_only_model.lr = reference.lr;

    PerceptronResult {
        epochs,
        stable,
        converged: stable && training_accuracy == 1.0,
        training_accuracy,
        validation_accuracy: validation_accuracy(&model),
        weights: [model.weights[0].to_f64(), model.weights[1].to_f64()],
        bias: model.bias.to_f64(),
        cast_weights: [
            cast_model.weights[0].to_f64(),
            cast_model.weights[1].to_f64(),
        ],
        cast_bias: cast_model.bias.to_f64(),
        cast_validation_accuracy: validation_accuracy(&cast_model),
        cast_reference_agreement: validation_agreement(&cast_model, reference),
        weight_only_reference_agreement: validation_agreement(&weight_only_model, reference),
    }
}

fn evaluate<T: Real + Copy>(
    samples: &[TrainingSample],
    reference: &Perceptron<f64, FEATURE_COUNT>,
) -> PerceptronResult {
    let (model, epochs, stable) = train_model::<T>(samples);
    summarize(&model, samples, epochs, stable, reference)
}

fn evaluate_all(samples: &[TrainingSample]) -> [PerceptronResult; 8] {
    let (reference, reference_epochs, reference_stable) = train_model::<f64>(samples);
    let reference_result = summarize(
        &reference,
        samples,
        reference_epochs,
        reference_stable,
        &reference,
    );
    [
        evaluate::<Fp4E2M1>(samples, &reference),
        evaluate::<Fp6E3M2>(samples, &reference),
        evaluate::<Fp8E4M3>(samples, &reference),
        evaluate::<Fp8E5M2>(samples, &reference),
        evaluate::<f16>(samples, &reference),
        evaluate::<f32>(samples, &reference),
        reference_result,
        evaluate::<f128>(samples, &reference),
    ]
}

fn write_csv(
    samples: &[TrainingSample],
    results: &[PerceptronResult; 8],
    path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "precision,training_samples,epochs,stable,converged,native_training_accuracy,native_validation_accuracy,native_weight_x,native_weight_y,native_bias,cast_f64_weight_x,cast_f64_weight_y,cast_f64_bias,cast_f64_validation_accuracy,cast_f64_reference_agreement,weight_only_reference_agreement"
    )?;
    for (index, result) in results.iter().enumerate() {
        writeln!(
            writer,
            "{},{},{},{},{},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            CSV_NAMES[index],
            samples.len(),
            result.epochs,
            result.stable,
            result.converged,
            result.training_accuracy,
            result.validation_accuracy,
            result.weights[0],
            result.weights[1],
            result.bias,
            result.cast_weights[0],
            result.cast_weights[1],
            result.cast_bias,
            result.cast_validation_accuracy,
            result.cast_reference_agreement,
            result.weight_only_reference_agreement,
        )?;
    }
    writer.flush()?;
    Ok(())
}

pub(crate) fn run() -> Result<(), Box<dyn std::error::Error>> {
    fs::create_dir_all(OUTPUT_DIR)?;
    let samples = training_samples();
    let results = evaluate_all(&samples);
    write_csv(&samples, &results, Path::new(CSV_PATH))?;

    println!("wrote {CSV_PATH}");
    println!("training samples = {}", samples.len());
    println!(
        "validation grid  = {} x {}",
        VALIDATION_GRID_STEPS + 1,
        VALIDATION_GRID_STEPS + 1
    );
    println!();
    for (index, result) in results.iter().enumerate() {
        let status = if result.converged {
            "converged"
        } else if result.stable {
            "stalled"
        } else {
            "epoch limit"
        };
        println!(
            "{:>13}: native val = {:7.3}%, cast val = {:7.3}%, cast/f64 = {:7.3}%, weights-only/f64 = {:7.3}%, epoch = {:>2}, status = {:>9}",
            LABELS[index],
            100.0 * result.validation_accuracy,
            100.0 * result.cast_validation_accuracy,
            100.0 * result.cast_reference_agreement,
            100.0 * result.weight_only_reference_agreement,
            result.epochs,
            status,
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_fpx_perceptrons_produce_finite_results() {
        let samples = training_samples();
        let results = evaluate_all(&samples);
        assert!(!samples.is_empty());
        for (label, result) in LABELS.iter().zip(results) {
            assert!(
                result.training_accuracy.is_finite()
                    && (0.0..=1.0).contains(&result.training_accuracy),
                "{label} produced an invalid training accuracy",
            );
            assert!(
                result.validation_accuracy.is_finite()
                    && (0.0..=1.0).contains(&result.validation_accuracy),
                "{label} produced an invalid validation accuracy",
            );
            assert!(
                result.cast_validation_accuracy.is_finite()
                    && (0.0..=1.0).contains(&result.cast_validation_accuracy),
                "{label} produced an invalid cast-f64 validation accuracy",
            );
            assert!(
                result.cast_reference_agreement.is_finite()
                    && (0.0..=1.0).contains(&result.cast_reference_agreement),
                "{label} produced an invalid cast-f64 agreement",
            );
            assert!(
                result.weight_only_reference_agreement.is_finite()
                    && (0.0..=1.0).contains(&result.weight_only_reference_agreement),
                "{label} produced an invalid weight-only agreement",
            );
            assert!(
                result.weights.into_iter().all(f64::is_finite)
                    && result.bias.is_finite()
                    && result.cast_weights.into_iter().all(f64::is_finite)
                    && result.cast_bias.is_finite(),
                "{label} produced non-finite parameters",
            );
        }
    }
}
