use aether_core::math::Matrix;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use nalgebra::{SMatrix, SVector};

fn aether_matmul<const N: usize>(a: Matrix<f64, N, N>, b: Matrix<f64, N, N>) -> Matrix<f64, N, N> {
    a * b
}

fn nalgebra_matmul<const N: usize>(a: SMatrix<f64, N, N>, b: SMatrix<f64, N, N>) -> SMatrix<f64, N, N> {
    a * b
}

fn aether_matvec<const N: usize>(a: Matrix<f64, N, N>, v: [f64; N]) -> [f64; N] {
    let v = aether_core::math::Vector::new(v);
    let out = a * v;
    out.data
}

fn nalgebra_matvec<const N: usize>(a: SMatrix<f64, N, N>, v: [f64; N]) -> [f64; N] {
    let v = SVector::<f64, N>::from_row_slice(&v);
    let out = a * v;
    out.as_slice().try_into().unwrap()
}

fn bench_matmul<const N: usize>(c: &mut Criterion) {
    let mut group = c.benchmark_group(format!("matmul_{N}x{N}"));
    group.throughput(Throughput::Elements((N * N * N) as u64));

    let aether_a = Matrix::new([[1.001_f64; N]; N]);
    let aether_b = Matrix::new([[0.999_f64; N]; N]);
    let na_a = SMatrix::<f64, N, N>::from_element(1.001_f64);
    let na_b = SMatrix::<f64, N, N>::from_element(0.999_f64);

    group.bench_function(BenchmarkId::new("aether", N), |b| {
        b.iter(|| aether_matmul::<N>(black_box(aether_a), black_box(aether_b)))
    });

    group.bench_function(BenchmarkId::new("nalgebra", N), |b| {
        b.iter(|| nalgebra_matmul::<N>(black_box(na_a), black_box(na_b)))
    });

    group.finish();
}

fn bench_matvec<const N: usize>(c: &mut Criterion) {
    let mut group = c.benchmark_group(format!("matvec_{N}x{N}"));
    group.throughput(Throughput::Elements((N * N) as u64));

    let aether_a = Matrix::new([[1.001_f64; N]; N]);
    let na_a = SMatrix::<f64, N, N>::from_element(1.001_f64);
    let vec = [0.999_f64; N];

    group.bench_function(BenchmarkId::new("aether", N), |b| {
        b.iter(|| aether_matvec::<N>(black_box(aether_a), black_box(vec)))
    });

    group.bench_function(BenchmarkId::new("nalgebra", N), |b| {
        b.iter(|| nalgebra_matvec::<N>(black_box(na_a), black_box(vec)))
    });

    group.finish();
}

fn throughput_benches(c: &mut Criterion) {
    bench_matmul::<3>(c);
    bench_matmul::<6>(c);
    bench_matvec::<3>(c);
    bench_matvec::<6>(c);
}

criterion_group!(benches, throughput_benches);
criterion_main!(benches);
