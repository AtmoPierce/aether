//! A small, trainable convolutional image classifier.
//!
//! The network uses one two-dimensional convolution, a ReLU activation,
//! global-average pooling, and a dense softmax classifier. Images and kernels
//! use Aether's fixed-size [`Matrix`] type, while pooled features and class
//! scores use [`Vector`]. Callers supply the convolution output dimensions as
//! const parameters; constructors and forward passes assert that they match the
//! dimensions implied by the input, kernel, stride, and padding.

use aether_core::math::{Matrix, Vector};
use aether_core::real::Real;

/// A channels-first image backed by Aether matrices.
pub type Image<F, const CHANNELS: usize, const HEIGHT: usize, const WIDTH: usize> =
    [Matrix<F, HEIGHT, WIDTH>; CHANNELS];

/// A trainable two-dimensional convolution layer.
///
/// `kernels[out_channel][in_channel]` is a `KERNEL x KERNEL` matrix.
#[derive(Clone, Copy, Debug)]
pub struct Conv2d<
    F: Real + Copy,
    const IN_CHANNELS: usize,
    const OUT_CHANNELS: usize,
    const KERNEL: usize,
> {
    pub kernels: [[Matrix<F, KERNEL, KERNEL>; IN_CHANNELS]; OUT_CHANNELS],
    pub bias: Vector<F, OUT_CHANNELS>,
    pub stride: usize,
    pub padding: usize,
}

impl<F: Real + Copy, const IN_CHANNELS: usize, const OUT_CHANNELS: usize, const KERNEL: usize>
    Conv2d<F, IN_CHANNELS, OUT_CHANNELS, KERNEL>
{
    /// Creates a convolution layer with deterministic He-style weights.
    pub fn new(stride: usize, padding: usize) -> Self {
        Self::with_seed(stride, padding, 0xA37E_4D91_2C6B_5F08)
    }

    /// Creates a convolution layer using a reproducible initialization seed.
    pub fn with_seed(stride: usize, padding: usize, seed: u64) -> Self {
        assert!(
            IN_CHANNELS > 0,
            "Conv2d requires at least one input channel"
        );
        assert!(
            OUT_CHANNELS > 0,
            "Conv2d requires at least one output channel"
        );
        assert!(KERNEL > 0, "Conv2d kernel size must be nonzero");
        assert!(stride > 0, "Conv2d stride must be nonzero");

        let zero_kernel = Matrix::new([[F::ZERO; KERNEL]; KERNEL]);
        let mut kernels = [[zero_kernel; IN_CHANNELS]; OUT_CHANNELS];
        let fan_in = IN_CHANNELS
            .checked_mul(KERNEL)
            .and_then(|count| count.checked_mul(KERNEL))
            .expect("Conv2d fan-in overflows usize");
        let scale = (F::from_f64(2.0) / F::from_usize(fan_in)).sqrt();
        let mut state = normalized_seed(seed);

        for output_kernels in &mut kernels {
            for kernel in output_kernels {
                for row in 0..KERNEL {
                    for column in 0..KERNEL {
                        kernel[(row, column)] = F::from_f64(next_symmetric(&mut state)) * scale;
                    }
                }
            }
        }

        Self {
            kernels,
            bias: Vector::new([F::ZERO; OUT_CHANNELS]),
            stride,
            padding,
        }
    }

    /// Returns the output height and width for an input matrix shape.
    pub fn output_shape<const HEIGHT: usize, const WIDTH: usize>(&self) -> (usize, usize) {
        (
            output_extent(HEIGHT, KERNEL, self.stride, self.padding, "height"),
            output_extent(WIDTH, KERNEL, self.stride, self.padding, "width"),
        )
    }

    /// Applies the convolution and returns fixed-size feature-map matrices.
    ///
    /// `OUTPUT_HEIGHT` and `OUTPUT_WIDTH` must equal
    /// `(input + 2 * padding - KERNEL) / stride + 1`.
    pub fn forward<
        const HEIGHT: usize,
        const WIDTH: usize,
        const OUTPUT_HEIGHT: usize,
        const OUTPUT_WIDTH: usize,
    >(
        &self,
        input: &Image<F, IN_CHANNELS, HEIGHT, WIDTH>,
    ) -> [Matrix<F, OUTPUT_HEIGHT, OUTPUT_WIDTH>; OUT_CHANNELS] {
        let actual_shape = self.output_shape::<HEIGHT, WIDTH>();
        assert_eq!(
            actual_shape,
            (OUTPUT_HEIGHT, OUTPUT_WIDTH),
            "Conv2d output constants must match the dimensions implied by the input, kernel, stride, and padding"
        );
        let zero_output = Matrix::new([[F::ZERO; OUTPUT_WIDTH]; OUTPUT_HEIGHT]);
        let mut output = [zero_output; OUT_CHANNELS];

        for out_channel in 0..OUT_CHANNELS {
            for out_y in 0..OUTPUT_HEIGHT {
                for out_x in 0..OUTPUT_WIDTH {
                    let mut sum = self.bias[out_channel];
                    for in_channel in 0..IN_CHANNELS {
                        for kernel_y in 0..KERNEL {
                            let padded_y = out_y * self.stride + kernel_y;
                            if padded_y < self.padding {
                                continue;
                            }
                            let input_y = padded_y - self.padding;
                            if input_y >= HEIGHT {
                                continue;
                            }

                            for kernel_x in 0..KERNEL {
                                let padded_x = out_x * self.stride + kernel_x;
                                if padded_x < self.padding {
                                    continue;
                                }
                                let input_x = padded_x - self.padding;
                                if input_x >= WIDTH {
                                    continue;
                                }

                                sum = sum
                                    + input[in_channel][(input_y, input_x)]
                                        * self.kernels[out_channel][in_channel]
                                            [(kernel_y, kernel_x)];
                            }
                        }
                    }
                    output[out_channel][(out_y, out_x)] = sum;
                }
            }
        }

        output
    }
}

/// A compact convolutional neural network for image classification.
///
/// The architecture is `Conv2d -> ReLU -> global average pool -> dense ->
/// softmax`. Training uses batch gradient descent and cross-entropy loss.
#[derive(Clone, Copy, Debug)]
pub struct ConvolutionalNeuralNetwork<
    F: Real + Copy,
    const CHANNELS: usize,
    const HEIGHT: usize,
    const WIDTH: usize,
    const FILTERS: usize,
    const KERNEL: usize,
    const OUTPUT_HEIGHT: usize,
    const OUTPUT_WIDTH: usize,
    const CLASSES: usize,
> {
    pub convolution: Conv2d<F, CHANNELS, FILTERS, KERNEL>,
    /// Dense weights in `[class, convolution_filter]` order.
    pub classifier_weights: Matrix<F, CLASSES, FILTERS>,
    pub classifier_bias: Vector<F, CLASSES>,
    pub learning_rate: F,
}

/// Short name for [`ConvolutionalNeuralNetwork`].
pub type Cnn<
    F,
    const CHANNELS: usize,
    const HEIGHT: usize,
    const WIDTH: usize,
    const FILTERS: usize,
    const KERNEL: usize,
    const OUTPUT_HEIGHT: usize,
    const OUTPUT_WIDTH: usize,
    const CLASSES: usize,
> = ConvolutionalNeuralNetwork<
    F,
    CHANNELS,
    HEIGHT,
    WIDTH,
    FILTERS,
    KERNEL,
    OUTPUT_HEIGHT,
    OUTPUT_WIDTH,
    CLASSES,
>;

impl<
        F: Real + Copy,
        const CHANNELS: usize,
        const HEIGHT: usize,
        const WIDTH: usize,
        const FILTERS: usize,
        const KERNEL: usize,
        const OUTPUT_HEIGHT: usize,
        const OUTPUT_WIDTH: usize,
        const CLASSES: usize,
    >
    ConvolutionalNeuralNetwork<
        F,
        CHANNELS,
        HEIGHT,
        WIDTH,
        FILTERS,
        KERNEL,
        OUTPUT_HEIGHT,
        OUTPUT_WIDTH,
        CLASSES,
    >
{
    /// Creates a CNN with stride one, no padding, and deterministic weights.
    pub fn new(learning_rate: F) -> Self {
        Self::with_seed(learning_rate, 0xC011_7A2D_D15C_A11E)
    }

    /// Creates a CNN using a reproducible initialization seed.
    pub fn with_seed(learning_rate: F, seed: u64) -> Self {
        Self::with_options(learning_rate, 1, 0, seed)
    }

    /// Creates a CNN with explicit convolution stride, padding, and seed.
    pub fn with_options(learning_rate: F, stride: usize, padding: usize, seed: u64) -> Self {
        assert!(CHANNELS > 0, "CNN requires at least one input channel");
        assert!(HEIGHT > 0, "CNN input height must be nonzero");
        assert!(WIDTH > 0, "CNN input width must be nonzero");
        assert!(FILTERS > 0, "CNN requires at least one convolution filter");
        assert!(CLASSES > 1, "CNN requires at least two classes");

        let convolution = Conv2d::with_seed(stride, padding, seed);
        let actual_shape = convolution.output_shape::<HEIGHT, WIDTH>();
        assert_eq!(
            actual_shape,
            (OUTPUT_HEIGHT, OUTPUT_WIDTH),
            "CNN output constants must match the dimensions implied by the input, kernel, stride, and padding"
        );

        let mut classifier_weights = Matrix::new([[F::ZERO; FILTERS]; CLASSES]);
        let scale = (F::ONE / F::from_usize(FILTERS)).sqrt();
        let mut state = normalized_seed(seed ^ 0x9E37_79B9_7F4A_7C15);
        for class in 0..CLASSES {
            for filter in 0..FILTERS {
                classifier_weights[(class, filter)] =
                    F::from_f64(next_symmetric(&mut state)) * scale;
            }
        }

        Self {
            convolution,
            classifier_weights,
            classifier_bias: Vector::new([F::ZERO; CLASSES]),
            learning_rate,
        }
    }

    /// Returns the pooled feature vector produced by the convolutional stage.
    pub fn extract_features(
        &self,
        input: &Image<F, CHANNELS, HEIGHT, WIDTH>,
    ) -> Vector<F, FILTERS> {
        self.forward_cache(input).features
    }

    /// Returns the unnormalized class scores for one image.
    pub fn logits(&self, input: &Image<F, CHANNELS, HEIGHT, WIDTH>) -> Vector<F, CLASSES> {
        self.forward_cache(input).logits
    }

    /// Returns softmax class probabilities for one image.
    pub fn predict_probabilities(
        &self,
        input: &Image<F, CHANNELS, HEIGHT, WIDTH>,
    ) -> Vector<F, CLASSES> {
        self.forward_cache(input).probabilities
    }

    /// Returns the class with the greatest softmax probability.
    pub fn predict(&self, input: &Image<F, CHANNELS, HEIGHT, WIDTH>) -> usize {
        let probabilities = self.predict_probabilities(input);
        let mut best_class = 0;
        for class in 1..CLASSES {
            if probabilities[class] > probabilities[best_class] {
                best_class = class;
            }
        }
        best_class
    }

    /// Returns cross-entropy loss for one labeled image.
    pub fn loss(&self, input: &Image<F, CHANNELS, HEIGHT, WIDTH>, label: usize) -> F {
        assert!(label < CLASSES, "CNN label {label} is out of range");
        let probabilities = self.predict_probabilities(input);
        -probabilities[label].max(F::EPSILON).ln()
    }

    /// Trains one batch-gradient epoch and returns its mean pre-update loss.
    pub fn train_epoch(
        &mut self,
        inputs: &[Image<F, CHANNELS, HEIGHT, WIDTH>],
        labels: &[usize],
    ) -> F {
        assert_eq!(
            inputs.len(),
            labels.len(),
            "CNN inputs and labels must have equal lengths"
        );
        if inputs.is_empty() {
            return F::ZERO;
        }

        let spatial_count = OUTPUT_HEIGHT * OUTPUT_WIDTH;
        let spatial_scale = F::ONE / F::from_usize(spatial_count);
        let zero_kernel = Matrix::new([[F::ZERO; KERNEL]; KERNEL]);
        let mut convolution_kernel_gradient = [[zero_kernel; CHANNELS]; FILTERS];
        let mut convolution_bias_gradient = Vector::new([F::ZERO; FILTERS]);
        let mut classifier_weight_gradient = Matrix::new([[F::ZERO; FILTERS]; CLASSES]);
        let mut classifier_bias_gradient = Vector::new([F::ZERO; CLASSES]);
        let mut loss = F::ZERO;

        for (input, &label) in inputs.iter().zip(labels.iter()) {
            assert!(label < CLASSES, "CNN label {label} is out of range");
            let cache = self.forward_cache(input);
            loss = loss - cache.probabilities[label].max(F::EPSILON).ln();

            let mut logits_gradient = cache.probabilities;
            logits_gradient[label] = logits_gradient[label] - F::ONE;
            let mut feature_gradient = Vector::new([F::ZERO; FILTERS]);

            for class in 0..CLASSES {
                classifier_bias_gradient[class] =
                    classifier_bias_gradient[class] + logits_gradient[class];
                for filter in 0..FILTERS {
                    classifier_weight_gradient[(class, filter)] = classifier_weight_gradient
                        [(class, filter)]
                        + logits_gradient[class] * cache.features[filter];
                    feature_gradient[filter] = feature_gradient[filter]
                        + self.classifier_weights[(class, filter)] * logits_gradient[class];
                }
            }

            for out_channel in 0..FILTERS {
                let pooled_gradient = feature_gradient[out_channel] * spatial_scale;
                for out_y in 0..OUTPUT_HEIGHT {
                    for out_x in 0..OUTPUT_WIDTH {
                        if cache.convolution_output[out_channel][(out_y, out_x)] <= F::ZERO {
                            continue;
                        }
                        convolution_bias_gradient[out_channel] =
                            convolution_bias_gradient[out_channel] + pooled_gradient;

                        for in_channel in 0..CHANNELS {
                            for kernel_y in 0..KERNEL {
                                let padded_y = out_y * self.convolution.stride + kernel_y;
                                if padded_y < self.convolution.padding {
                                    continue;
                                }
                                let input_y = padded_y - self.convolution.padding;
                                if input_y >= HEIGHT {
                                    continue;
                                }

                                for kernel_x in 0..KERNEL {
                                    let padded_x = out_x * self.convolution.stride + kernel_x;
                                    if padded_x < self.convolution.padding {
                                        continue;
                                    }
                                    let input_x = padded_x - self.convolution.padding;
                                    if input_x >= WIDTH {
                                        continue;
                                    }

                                    let gradient = &mut convolution_kernel_gradient[out_channel]
                                        [in_channel][(kernel_y, kernel_x)];
                                    *gradient = *gradient
                                        + pooled_gradient * input[in_channel][(input_y, input_x)];
                                }
                            }
                        }
                    }
                }
            }
        }

        let batch_size = F::from_usize(inputs.len());
        let step = self.learning_rate / batch_size;
        for out_channel in 0..FILTERS {
            for in_channel in 0..CHANNELS {
                for kernel_y in 0..KERNEL {
                    for kernel_x in 0..KERNEL {
                        self.convolution.kernels[out_channel][in_channel][(kernel_y, kernel_x)] =
                            self.convolution.kernels[out_channel][in_channel][(kernel_y, kernel_x)]
                                - step
                                    * convolution_kernel_gradient[out_channel][in_channel]
                                        [(kernel_y, kernel_x)];
                    }
                }
            }
            self.convolution.bias[out_channel] =
                self.convolution.bias[out_channel] - step * convolution_bias_gradient[out_channel];
        }
        for class in 0..CLASSES {
            for filter in 0..FILTERS {
                self.classifier_weights[(class, filter)] = self.classifier_weights[(class, filter)]
                    - step * classifier_weight_gradient[(class, filter)];
            }
            self.classifier_bias[class] =
                self.classifier_bias[class] - step * classifier_bias_gradient[class];
        }

        loss / batch_size
    }

    /// Runs several training epochs and returns the final epoch's loss.
    pub fn fit(
        &mut self,
        inputs: &[Image<F, CHANNELS, HEIGHT, WIDTH>],
        labels: &[usize],
        epochs: usize,
    ) -> F {
        let mut loss = F::ZERO;
        for _ in 0..epochs {
            loss = self.train_epoch(inputs, labels);
        }
        loss
    }

    fn forward_cache(
        &self,
        input: &Image<F, CHANNELS, HEIGHT, WIDTH>,
    ) -> ForwardCache<F, FILTERS, OUTPUT_HEIGHT, OUTPUT_WIDTH, CLASSES> {
        let convolution_output: [Matrix<F, OUTPUT_HEIGHT, OUTPUT_WIDTH>; FILTERS] =
            self.convolution.forward(input);
        let spatial_count = OUTPUT_HEIGHT * OUTPUT_WIDTH;
        let spatial_scale = F::ONE / F::from_usize(spatial_count);
        let mut features = Vector::new([F::ZERO; FILTERS]);

        for out_channel in 0..FILTERS {
            let mut sum = F::ZERO;
            for out_y in 0..OUTPUT_HEIGHT {
                for out_x in 0..OUTPUT_WIDTH {
                    sum = sum + convolution_output[out_channel][(out_y, out_x)].max(F::ZERO);
                }
            }
            features[out_channel] = sum * spatial_scale;
        }

        let mut logits = &self.classifier_weights * &features;
        for class in 0..CLASSES {
            logits[class] = logits[class] + self.classifier_bias[class];
        }
        let probabilities = softmax(&logits);

        ForwardCache {
            convolution_output,
            features,
            logits,
            probabilities,
        }
    }
}

#[derive(Debug)]
struct ForwardCache<
    F,
    const FILTERS: usize,
    const OUTPUT_HEIGHT: usize,
    const OUTPUT_WIDTH: usize,
    const CLASSES: usize,
> {
    convolution_output: [Matrix<F, OUTPUT_HEIGHT, OUTPUT_WIDTH>; FILTERS],
    features: Vector<F, FILTERS>,
    logits: Vector<F, CLASSES>,
    probabilities: Vector<F, CLASSES>,
}

fn softmax<F: Real + Copy, const N: usize>(logits: &Vector<F, N>) -> Vector<F, N> {
    debug_assert!(N > 0);
    let mut maximum = logits[0];
    for index in 1..N {
        if logits[index] > maximum {
            maximum = logits[index];
        }
    }

    let mut denominator = F::ZERO;
    let mut probabilities = Vector::new([F::ZERO; N]);
    for index in 0..N {
        probabilities[index] = (logits[index] - maximum).exp();
        denominator = denominator + probabilities[index];
    }
    for index in 0..N {
        probabilities[index] = probabilities[index] / denominator;
    }
    probabilities
}

fn output_extent(
    input: usize,
    kernel: usize,
    stride: usize,
    padding: usize,
    dimension_name: &str,
) -> usize {
    let double_padding = padding
        .checked_mul(2)
        .expect("Conv2d padding overflows usize");
    let padded_input = input
        .checked_add(double_padding)
        .expect("Conv2d padded input overflows usize");
    assert!(
        padded_input >= kernel,
        "Conv2d kernel is larger than padded input {dimension_name}"
    );
    (padded_input - kernel) / stride + 1
}

#[inline]
fn normalized_seed(seed: u64) -> u64 {
    if seed == 0 {
        0xD1B5_4A32_D192_ED03
    } else {
        seed
    }
}

#[inline]
fn next_symmetric(state: &mut u64) -> f64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    let unit = (*state >> 11) as f64 / ((1_u64 << 53) as f64);
    unit * 2.0 - 1.0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1.0e-10,
            "expected {expected}, got {actual}"
        );
    }

    #[test]
    fn conv2d_computes_known_valid_convolution() {
        let mut convolution = Conv2d::<f64, 1, 1, 2>::new(1, 0);
        convolution.kernels[0][0] = Matrix::new([[1.0, 1.0], [1.0, 1.0]]);
        convolution.bias[0] = 0.0;
        let input = [Matrix::new([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ])];

        let output: [Matrix<f64, 2, 2>; 1] = convolution.forward(&input);

        for row in 0..2 {
            for column in 0..2 {
                let expected = [[12.0, 16.0], [24.0, 28.0]][row][column];
                assert_close(output[0][(row, column)], expected);
            }
        }
    }

    #[test]
    #[should_panic(expected = "CNN output constants must match")]
    fn cnn_rejects_incorrect_output_dimensions() {
        let _ = Cnn::<f64, 1, 4, 4, 2, 2, 2, 2, 2>::new(0.05);
    }

    #[test]
    fn probabilities_are_normalized() {
        let network = Cnn::<f64, 1, 4, 4, 3, 2, 3, 3, 3>::new(0.05);
        let input = [Matrix::new([[0.5; 4]; 4])];
        let probabilities = network.predict_probabilities(&input);
        let sum: f64 = probabilities.data.iter().sum();

        assert_close(sum, 1.0);
        assert!(probabilities
            .data
            .iter()
            .all(|&probability| probability > 0.0));
    }

    #[test]
    fn cnn_learns_horizontal_and_vertical_lines() {
        let images = [
            [Matrix::new([
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
            ])],
            [Matrix::new([
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
            ])],
            [Matrix::new([
                [0.0, 0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0, 1.0],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
            ])],
            [Matrix::new([
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0, 1.0],
                [0.0, 0.0, 0.0, 0.0],
            ])],
        ];
        let labels = [0, 0, 1, 1];
        let mut network = Cnn::<f64, 1, 4, 4, 2, 2, 3, 3, 2>::new(0.2);

        network.convolution.kernels[0][0] = Matrix::new([[1.0, -1.0], [1.0, -1.0]]);
        network.convolution.kernels[1][0] = Matrix::new([[1.0, 1.0], [-1.0, -1.0]]);
        network.convolution.bias = Vector::new([0.0; 2]);
        network.classifier_weights = Matrix::new([[0.0; 2]; 2]);
        network.classifier_bias = Vector::new([0.0; 2]);

        let initial_loss: f64 = images
            .iter()
            .zip(labels.iter())
            .map(|(image, &label)| network.loss(image, label))
            .sum::<f64>()
            / images.len() as f64;
        network.fit(&images, &labels, 300);
        let final_loss: f64 = images
            .iter()
            .zip(labels.iter())
            .map(|(image, &label)| network.loss(image, label))
            .sum::<f64>()
            / images.len() as f64;

        assert!(final_loss < initial_loss * 0.25);
        let predictions: [usize; 4] = core::array::from_fn(|index| network.predict(&images[index]));
        assert_eq!(predictions, labels);
    }
}
