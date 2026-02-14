use aether_core::attitude::{DirectionCosineMatrix, Euler, Quaternion};
use aether_core::coordinate::Cartesian;
use aether_core::math::{Matrix, Vector};
use aether_core::reference_frame::Unknown;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use nalgebra::{DMatrix, DVector, SMatrix, SVector};

#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::*;

fn aether_matmul<const N: usize>(a: &Matrix<f64, N, N>, b: &Matrix<f64, N, N>) -> Matrix<f64, N, N> {
    a * b
}

fn aether_matmul_f32<const N: usize>(a: &Matrix<f32, N, N>, b: &Matrix<f32, N, N>) -> Matrix<f32, N, N> {
    a * b
}

fn nalgebra_matmul<const N: usize>(a: &SMatrix<f64, N, N>, b: &SMatrix<f64, N, N>) -> SMatrix<f64, N, N> {
    a * b
}

fn nalgebra_matmul_f32<const N: usize>(a: &SMatrix<f32, N, N>, b: &SMatrix<f32, N, N>) -> SMatrix<f32, N, N> {
    a * b
}

fn aether_matvec<const N: usize>(a: &Matrix<f64, N, N>, v: &Vector<f64, N>) -> [f64; N] {
    let out = a * v;
    out.data
}

fn aether_matvec_f32<const N: usize>(a: &Matrix<f32, N, N>, v: &Vector<f32, N>) -> [f32; N] {
    let out = a * v;
    out.data
}

fn aether_matvec6_simd(a: &Matrix<f64, 6, 6>, v: &Vector<f64, 6>) -> [f64; 6] {
    a.mul_vec6_simd(v).data
}

fn aether_matvec4_simd(a: &Matrix<f64, 4, 4>, v: &Vector<f64, 4>) -> [f64; 4] {
    a.mul_vec4_simd(v).data
}

fn aether_matvec6_simd_f32(a: &Matrix<f32, 6, 6>, v: &Vector<f32, 6>) -> [f32; 6] {
    a.mul_vec6_simd(v).data
}

fn aether_matvec4_simd_f32(a: &Matrix<f32, 4, 4>, v: &Vector<f32, 4>) -> [f32; 4] {
    a.mul_vec4_simd(v).data
}

fn nalgebra_matvec<const N: usize>(a: &SMatrix<f64, N, N>, v: &SVector<f64, N>) -> [f64; N] {
    let out = a * v;
    out.as_slice().try_into().unwrap()
}

fn nalgebra_matvec_f32<const N: usize>(a: &SMatrix<f32, N, N>, v: &SVector<f32, N>) -> [f32; N] {
    let out = a * v;
    out.as_slice().try_into().unwrap()
}

fn aether_dot4_simd(a: &Vector<f64, 4>, b: &Vector<f64, 4>) -> f64 {
    a.dot4_simd(b)
}

fn aether_dot4_scalar(a: &Vector<f64, 4>, b: &Vector<f64, 4>) -> f64 {
    a.dot(b)
}

fn aether_dot4_simd_f32(a: &Vector<f32, 4>, b: &Vector<f32, 4>) -> f32 {
    a.dot4_simd(b)
}

fn aether_dot4_scalar_f32(a: &Vector<f32, 4>, b: &Vector<f32, 4>) -> f32 {
    a.dot(b)
}

fn nalgebra_dot4(a: &SVector<f64, 4>, b: &SVector<f64, 4>) -> f64 {
    a.dot(b)
}

fn nalgebra_dot4_f32(a: &SVector<f32, 4>, b: &SVector<f32, 4>) -> f32 {
    a.dot(b)
}

fn native_vec3_dot(a: &Vector<f64, 3>, b: &Vector<f64, 3>) -> f64 {
    a.dot(b)
}

fn native_vec3_cross(a: &Vector<f64, 3>, b: &Vector<f64, 3>) -> Vector<f64, 3> {
    a.cross(b)
}

fn native_mat3_vec3(m: &Matrix<f64, 3, 3>, v: &Vector<f64, 3>) -> Vector<f64, 3> {
    m * v
}

fn native_mat3_mul(a: &Matrix<f64, 3, 3>, b: &Matrix<f64, 3, 3>) -> Matrix<f64, 3, 3> {
    a * b
}

fn native_mat3_mul_owned(a: Matrix<f64, 3, 3>, b: Matrix<f64, 3, 3>) -> Matrix<f64, 3, 3> {
    a * b
}

fn native_quat_rotate_vector_equivalent(
    q: &Quaternion<f64, Unknown, Unknown>,
    v_from: &Vector<f64, 3>,
) -> Vector<f64, 3> {
    let qn = q.normalized();
    let qc = qn.conjugate();

    let px = v_from[0];
    let py = v_from[1];
    let pz = v_from[2];

    let aw = -(qc.i() * px + qc.j() * py + qc.k() * pz);
    let ax = qc.w() * px + qc.j() * pz - qc.k() * py;
    let ay = qc.w() * py + qc.k() * px - qc.i() * pz;
    let az = qc.w() * pz + qc.i() * py - qc.j() * px;

    Vector::new([
        aw * qn.i() + ax * qn.w() + ay * qn.k() - az * qn.j(),
        aw * qn.j() - ax * qn.k() + ay * qn.w() + az * qn.i(),
        aw * qn.k() + ax * qn.j() - ay * qn.i() + az * qn.w(),
    ])
}

fn matvec_naive_rowmajor(a: &[f64], x: &[f64], y: &mut [f64], n: usize) {
    for i in 0..n {
        let row = &a[i * n..(i + 1) * n];
        let mut sum = 0.0;
        for j in 0..n {
            sum += row[j] * x[j];
        }
        y[i] = sum;
    }
}

fn matvec_naive_rowmajor_f32(a: &[f32], x: &[f32], y: &mut [f32], n: usize) {
    for i in 0..n {
        let row = &a[i * n..(i + 1) * n];
        let mut sum = 0.0_f32;
        for j in 0..n {
            sum += row[j] * x[j];
        }
        y[i] = sum;
    }
}

fn matmul_naive_rowmajor(a: &[f64], b: &[f64], c: &mut [f64], n: usize) {
    for i in 0..n {
        for j in 0..n {
            let mut sum = 0.0_f64;
            for k in 0..n {
                sum += a[i * n + k] * b[k * n + j];
            }
            c[i * n + j] = sum;
        }
    }
}

fn matmul_ikj_rowmajor(a: &[f64], b: &[f64], c: &mut [f64], n: usize) {
    c.fill(0.0_f64);
    for i in 0..n {
        let a_row = &a[i * n..(i + 1) * n];
        let c_row = &mut c[i * n..(i + 1) * n];
        for k in 0..n {
            let aik = a_row[k];
            let b_row = &b[k * n..(k + 1) * n];
            for j in 0..n {
                c_row[j] += aik * b_row[j];
            }
        }
    }
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn matmul_simd_neon_rowmajor(a: &[f64], b: &[f64], c: &mut [f64], n: usize) {
    c.fill(0.0_f64);
    for i in 0..n {
        let a_row = &a[i * n..(i + 1) * n];
        let c_row = &mut c[i * n..(i + 1) * n];
        for k in 0..n {
            let aik = a_row[k];
            let aik_vec = vdupq_n_f64(aik);
            let b_row = &b[k * n..(k + 1) * n];

            let mut j = 0;
            while j + 2 <= n {
                let c_vec = vld1q_f64(c_row.as_ptr().add(j));
                let b_vec = vld1q_f64(b_row.as_ptr().add(j));
                let out = vmlaq_f64(c_vec, b_vec, aik_vec);
                vst1q_f64(c_row.as_mut_ptr().add(j), out);
                j += 2;
            }

            while j < n {
                c_row[j] += aik * b_row[j];
                j += 1;
            }
        }
    }
}

fn matmul_naive_rowmajor_f32(a: &[f32], b: &[f32], c: &mut [f32], n: usize) {
    for i in 0..n {
        for j in 0..n {
            let mut sum = 0.0_f32;
            for k in 0..n {
                sum += a[i * n + k] * b[k * n + j];
            }
            c[i * n + j] = sum;
        }
    }
}

fn matmul_ikj_rowmajor_f32(a: &[f32], b: &[f32], c: &mut [f32], n: usize) {
    c.fill(0.0_f32);
    for i in 0..n {
        let a_row = &a[i * n..(i + 1) * n];
        let c_row = &mut c[i * n..(i + 1) * n];
        for k in 0..n {
            let aik = a_row[k];
            let b_row = &b[k * n..(k + 1) * n];
            for j in 0..n {
                c_row[j] += aik * b_row[j];
            }
        }
    }
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn matmul_simd_neon_rowmajor_f32(a: &[f32], b: &[f32], c: &mut [f32], n: usize) {
    c.fill(0.0_f32);
    for i in 0..n {
        let a_row = &a[i * n..(i + 1) * n];
        let c_row = &mut c[i * n..(i + 1) * n];
        for k in 0..n {
            let aik = a_row[k];
            let aik_vec = vdupq_n_f32(aik);
            let b_row = &b[k * n..(k + 1) * n];

            let mut j = 0;
            while j + 4 <= n {
                let c_vec = vld1q_f32(c_row.as_ptr().add(j));
                let b_vec = vld1q_f32(b_row.as_ptr().add(j));
                let out = vmlaq_f32(c_vec, b_vec, aik_vec);
                vst1q_f32(c_row.as_mut_ptr().add(j), out);
                j += 4;
            }

            while j < n {
                c_row[j] += aik * b_row[j];
                j += 1;
            }
        }
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx")]
unsafe fn matvec_simd_avx_rowmajor(a: &[f64], x: &[f64], y: &mut [f64], n: usize) {
    for i in 0..n {
        let row = &a[i * n..(i + 1) * n];
        let mut acc = _mm256_setzero_pd();

        let mut j = 0;
        while j + 4 <= n {
            let r = unsafe { _mm256_loadu_pd(row.as_ptr().add(j)) };
            let v = unsafe { _mm256_loadu_pd(x.as_ptr().add(j)) };
            let p = _mm256_mul_pd(r, v);
            acc = _mm256_add_pd(acc, p);
            j += 4;
        }

        let mut lanes = [0.0_f64; 4];
        unsafe { _mm256_storeu_pd(lanes.as_mut_ptr(), acc) };
        let mut sum = lanes[0] + lanes[1] + lanes[2] + lanes[3];

        while j < n {
            sum += row[j] * x[j];
            j += 1;
        }

        y[i] = sum;
    }
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn matvec_simd_neon_rowmajor(a: &[f64], x: &[f64], y: &mut [f64], n: usize) {
    for i in 0..n {
        let row = &a[i * n..(i + 1) * n];
        let mut acc = vdupq_n_f64(0.0_f64);

        let mut j = 0;
        while j + 2 <= n {
            let r = vld1q_f64(row.as_ptr().add(j));
            let v = vld1q_f64(x.as_ptr().add(j));
            acc = vaddq_f64(acc, vmulq_f64(r, v));
            j += 2;
        }

        let mut sum = vgetq_lane_f64(acc, 0) + vgetq_lane_f64(acc, 1);

        while j < n {
            sum += row[j] * x[j];
            j += 1;
        }

        y[i] = sum;
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx,fma")]
unsafe fn matvec_fma_avx_rowmajor(a: &[f64], x: &[f64], y: &mut [f64], n: usize) {
    for i in 0..n {
        let row = &a[i * n..(i + 1) * n];
        let mut acc = _mm256_setzero_pd();

        let mut j = 0;
        while j + 4 <= n {
            let r = unsafe { _mm256_loadu_pd(row.as_ptr().add(j)) };
            let v = unsafe { _mm256_loadu_pd(x.as_ptr().add(j)) };
            acc = _mm256_fmadd_pd(r, v, acc);
            j += 4;
        }

        let mut lanes = [0.0_f64; 4];
        unsafe { _mm256_storeu_pd(lanes.as_mut_ptr(), acc) };
        let mut sum = lanes[0] + lanes[1] + lanes[2] + lanes[3];

        while j < n {
            sum += row[j] * x[j];
            j += 1;
        }

        y[i] = sum;
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx")]
unsafe fn matvec_simd_avx_rowmajor_f32(a: &[f32], x: &[f32], y: &mut [f32], n: usize) {
    for i in 0..n {
        let row = &a[i * n..(i + 1) * n];
        let mut acc = _mm256_setzero_ps();

        let mut j = 0;
        while j + 8 <= n {
            let r = unsafe { _mm256_loadu_ps(row.as_ptr().add(j)) };
            let v = unsafe { _mm256_loadu_ps(x.as_ptr().add(j)) };
            let p = _mm256_mul_ps(r, v);
            acc = _mm256_add_ps(acc, p);
            j += 8;
        }

        let mut lanes = [0.0_f32; 8];
        unsafe { _mm256_storeu_ps(lanes.as_mut_ptr(), acc) };
        let mut sum = lanes[0] + lanes[1] + lanes[2] + lanes[3] + lanes[4] + lanes[5] + lanes[6] + lanes[7];

        while j < n {
            sum += row[j] * x[j];
            j += 1;
        }

        y[i] = sum;
    }
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn matvec_simd_neon_rowmajor_f32(a: &[f32], x: &[f32], y: &mut [f32], n: usize) {
    for i in 0..n {
        let row = &a[i * n..(i + 1) * n];
        let mut acc = vdupq_n_f32(0.0_f32);

        let mut j = 0;
        while j + 4 <= n {
            let r = vld1q_f32(row.as_ptr().add(j));
            let v = vld1q_f32(x.as_ptr().add(j));
            acc = vaddq_f32(acc, vmulq_f32(r, v));
            j += 4;
        }

        let mut lanes = [0.0_f32; 4];
        vst1q_f32(lanes.as_mut_ptr(), acc);
        let mut sum = lanes[0] + lanes[1] + lanes[2] + lanes[3];

        while j < n {
            sum += row[j] * x[j];
            j += 1;
        }

        y[i] = sum;
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx,fma")]
unsafe fn matvec_fma_avx_rowmajor_f32(a: &[f32], x: &[f32], y: &mut [f32], n: usize) {
    for i in 0..n {
        let row = &a[i * n..(i + 1) * n];
        let mut acc = _mm256_setzero_ps();

        let mut j = 0;
        while j + 8 <= n {
            let r = unsafe { _mm256_loadu_ps(row.as_ptr().add(j)) };
            let v = unsafe { _mm256_loadu_ps(x.as_ptr().add(j)) };
            acc = _mm256_fmadd_ps(r, v, acc);
            j += 8;
        }

        let mut lanes = [0.0_f32; 8];
        unsafe { _mm256_storeu_ps(lanes.as_mut_ptr(), acc) };
        let mut sum = lanes[0] + lanes[1] + lanes[2] + lanes[3] + lanes[4] + lanes[5] + lanes[6] + lanes[7];

        while j < n {
            sum += row[j] * x[j];
            j += 1;
        }

        y[i] = sum;
    }
}

fn avx_available() -> bool {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        std::is_x86_feature_detected!("avx")
    }

    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        false
    }
}

fn neon_available() -> bool {
    #[cfg(target_arch = "aarch64")]
    {
        std::arch::is_aarch64_feature_detected!("neon")
    }

    #[cfg(not(target_arch = "aarch64"))]
    {
        false
    }
}

fn simd_available() -> bool {
    avx_available() || neon_available()
}

fn avx_fma_available() -> bool {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        std::is_x86_feature_detected!("avx") && std::is_x86_feature_detected!("fma")
    }

    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        false
    }
}

fn bench_matvec_large(c: &mut Criterion, n: usize) {
    let mut a = vec![0.0_f64; n * n];
    for i in 0..n {
        for j in 0..n {
            a[i * n + j] = 1.001 + ((i + j) % 7) as f64 * 1e-3;
        }
    }

    let x = vec![0.999_f64; n];

    let na_a = DMatrix::<f64>::from_row_slice(n, n, &a);
    let na_x = DVector::<f64>::from_row_slice(&x);
    let mut na_y = DVector::<f64>::zeros(n);

    c.bench_function(&format!("matvec_{}_nalgebra", n), |b| {
        b.iter(|| {
            na_y.gemv(1.0, black_box(&na_a), black_box(&na_x), 0.0);
            black_box(na_y[0])
        })
    });
}

fn bench_matvec_large_f32(c: &mut Criterion, n: usize) {
    let mut a = vec![0.0_f32; n * n];
    for i in 0..n {
        for j in 0..n {
            a[i * n + j] = 1.001_f32 + ((i + j) % 7) as f32 * 1e-3_f32;
        }
    }

    let x = vec![0.999_f32; n];

    let na_a = DMatrix::<f32>::from_row_slice(n, n, &a);
    let na_x = DVector::<f32>::from_row_slice(&x);
    let mut na_y = DVector::<f32>::zeros(n);

    c.bench_function(&format!("matvec_{}_nalgebra_f32", n), |b| {
        b.iter(|| {
            na_y.gemv(1.0_f32, black_box(&na_a), black_box(&na_x), 0.0_f32);
            black_box(na_y[0])
        })
    });
}

fn bench_matmul_large(c: &mut Criterion, n: usize) {
    let mut a = vec![0.0_f64; n * n];
    let mut b = vec![0.0_f64; n * n];
    for i in 0..n {
        for j in 0..n {
            a[i * n + j] = 1.001 + ((i + j) % 7) as f64 * 1e-3;
            b[i * n + j] = 0.999 + ((i + 2 * j) % 5) as f64 * 1e-3;
        }
    }

    #[cfg(target_arch = "aarch64")]
    let mut c_simd = vec![0.0_f64; n * n];

    let na_a = DMatrix::<f64>::from_row_slice(n, n, &a);
    let na_b = DMatrix::<f64>::from_row_slice(n, n, &b);
    let mut na_c = DMatrix::<f64>::zeros(n, n);

    #[cfg(target_arch = "aarch64")]
    c.bench_function(&format!("matmul_{}_simd", n), |bch| {
        bch.iter(|| {
            unsafe {
                matmul_simd_neon_rowmajor(black_box(&a), black_box(&b), black_box(&mut c_simd), n);
            }
            black_box(c_simd[0])
        })
    });

    c.bench_function(&format!("matmul_{}_nalgebra", n), |bch| {
        bch.iter(|| {
            na_c.gemm(1.0_f64, black_box(&na_a), black_box(&na_b), 0.0_f64);
            black_box(na_c[(0, 0)])
        })
    });
}

fn bench_matmul_large_f32(c: &mut Criterion, n: usize) {
    let mut a = vec![0.0_f32; n * n];
    let mut b = vec![0.0_f32; n * n];
    for i in 0..n {
        for j in 0..n {
            a[i * n + j] = 1.001_f32 + ((i + j) % 7) as f32 * 1e-3_f32;
            b[i * n + j] = 0.999_f32 + ((i + 2 * j) % 5) as f32 * 1e-3_f32;
        }
    }

    #[cfg(target_arch = "aarch64")]
    let mut c_simd = vec![0.0_f32; n * n];

    let na_a = DMatrix::<f32>::from_row_slice(n, n, &a);
    let na_b = DMatrix::<f32>::from_row_slice(n, n, &b);
    let mut na_c = DMatrix::<f32>::zeros(n, n);

    #[cfg(target_arch = "aarch64")]
    c.bench_function(&format!("matmul_{}_simd_f32", n), |bch| {
        bch.iter(|| {
            unsafe {
                matmul_simd_neon_rowmajor_f32(black_box(&a), black_box(&b), black_box(&mut c_simd), n);
            }
            black_box(c_simd[0])
        })
    });

    c.bench_function(&format!("matmul_{}_nalgebra_f32", n), |bch| {
        bch.iter(|| {
            na_c.gemm(1.0_f32, black_box(&na_a), black_box(&na_b), 0.0_f32);
            black_box(na_c[(0, 0)])
        })
    });
}

fn bench_matmul<const N: usize>(c: &mut Criterion) {
    let aether_a = Matrix::new([[1.001_f64; N]; N]);
    let aether_b = Matrix::new([[0.999_f64; N]; N]);
    let na_a = SMatrix::<f64, N, N>::from_element(1.001_f64);
    let na_b = SMatrix::<f64, N, N>::from_element(0.999_f64);

    c.bench_function(&format!("matmul_{}_aether", N), |b| {
        b.iter(|| aether_matmul::<N>(black_box(&aether_a), black_box(&aether_b)))
    });

    c.bench_function(&format!("matmul_{}_nalgebra", N), |b| {
        b.iter(|| nalgebra_matmul::<N>(black_box(&na_a), black_box(&na_b)))
    });
}

fn bench_matmul_f32<const N: usize>(c: &mut Criterion) {
    let aether_a = Matrix::new([[1.001_f32; N]; N]);
    let aether_b = Matrix::new([[0.999_f32; N]; N]);
    let na_a = SMatrix::<f32, N, N>::from_element(1.001_f32);
    let na_b = SMatrix::<f32, N, N>::from_element(0.999_f32);

    c.bench_function(&format!("matmul_{}_aether_f32", N), |b| {
        b.iter(|| aether_matmul_f32::<N>(black_box(&aether_a), black_box(&aether_b)))
    });

    c.bench_function(&format!("matmul_{}_nalgebra_f32", N), |b| {
        b.iter(|| nalgebra_matmul_f32::<N>(black_box(&na_a), black_box(&na_b)))
    });
}

fn bench_matvec<const N: usize>(c: &mut Criterion) {
    let aether_a = Matrix::new([[1.001_f64; N]; N]);
    let na_a = SMatrix::<f64, N, N>::from_element(1.001_f64);
    let aether_v = Vector::new([0.999_f64; N]);
    let na_v = SVector::<f64, N>::from_element(0.999_f64);

    c.bench_function(&format!("matvec_{}_aether", N), |b| {
        b.iter(|| aether_matvec::<N>(black_box(&aether_a), black_box(&aether_v)))
    });

    c.bench_function(&format!("matvec_{}_nalgebra", N), |b| {
        b.iter(|| nalgebra_matvec::<N>(black_box(&na_a), black_box(&na_v)))
    });
}

fn bench_matvec_f32<const N: usize>(c: &mut Criterion) {
    let aether_a = Matrix::new([[1.001_f32; N]; N]);
    let na_a = SMatrix::<f32, N, N>::from_element(1.001_f32);
    let aether_v = Vector::new([0.999_f32; N]);
    let na_v = SVector::<f32, N>::from_element(0.999_f32);

    c.bench_function(&format!("matvec_{}_aether_f32", N), |b| {
        b.iter(|| aether_matvec_f32::<N>(black_box(&aether_a), black_box(&aether_v)))
    });

    c.bench_function(&format!("matvec_{}_nalgebra_f32", N), |b| {
        b.iter(|| nalgebra_matvec_f32::<N>(black_box(&na_a), black_box(&na_v)))
    });
}

fn bench_matvec_6(c: &mut Criterion) {
    let aether_a = Matrix::new([[1.001_f64; 6]; 6]);
    let na_a = SMatrix::<f64, 6, 6>::from_element(1.001_f64);
    let aether_v = Vector::new([0.999_f64; 6]);
    let na_v = SVector::<f64, 6>::from_element(0.999_f64);

    c.bench_function("matvec_6_aether", |b| {
        b.iter(|| aether_matvec6_simd(black_box(&aether_a), black_box(&aether_v)))
    });

    c.bench_function("matvec_6_nalgebra", |b| {
        b.iter(|| nalgebra_matvec::<6>(black_box(&na_a), black_box(&na_v)))
    });
}

fn bench_matvec_4(c: &mut Criterion) {
    let aether_a = Matrix::new([[1.001_f64; 4]; 4]);
    let na_a = SMatrix::<f64, 4, 4>::from_element(1.001_f64);
    let aether_v = Vector::new([0.999_f64; 4]);
    let na_v = SVector::<f64, 4>::from_element(0.999_f64);

    c.bench_function("matvec_4_aether", |b| {
        b.iter(|| aether_matvec4_simd(black_box(&aether_a), black_box(&aether_v)))
    });

    c.bench_function("matvec_4_nalgebra", |b| {
        b.iter(|| nalgebra_matvec::<4>(black_box(&na_a), black_box(&na_v)))
    });
}

fn bench_matvec_4_f32(c: &mut Criterion) {
    let aether_a = Matrix::new([[1.001_f32; 4]; 4]);
    let na_a = SMatrix::<f32, 4, 4>::from_element(1.001_f32);
    let aether_v = Vector::new([0.999_f32; 4]);
    let na_v = SVector::<f32, 4>::from_element(0.999_f32);

    c.bench_function("matvec_4_aether_f32", |b| {
        b.iter(|| aether_matvec4_simd_f32(black_box(&aether_a), black_box(&aether_v)))
    });

    c.bench_function("matvec_4_nalgebra_f32", |b| {
        b.iter(|| nalgebra_matvec_f32::<4>(black_box(&na_a), black_box(&na_v)))
    });
}

fn bench_dot_4(c: &mut Criterion) {
    let aether_a = Vector::new([1.001_f64; 4]);
    let aether_b = Vector::new([0.999_f64; 4]);
    let na_a = SVector::<f64, 4>::from_element(1.001_f64);
    let na_b = SVector::<f64, 4>::from_element(0.999_f64);

    c.bench_function("dot_4_aether", |b| {
        b.iter(|| aether_dot4_simd(black_box(&aether_a), black_box(&aether_b)))
    });

    c.bench_function("dot_4_scalar", |b| {
        b.iter(|| aether_dot4_scalar(black_box(&aether_a), black_box(&aether_b)))
    });

    c.bench_function("dot_4_nalgebra", |b| {
        b.iter(|| nalgebra_dot4(black_box(&na_a), black_box(&na_b)))
    });
}

fn bench_dot_4_f32(c: &mut Criterion) {
    let aether_a = Vector::new([1.001_f32; 4]);
    let aether_b = Vector::new([0.999_f32; 4]);
    let na_a = SVector::<f32, 4>::from_element(1.001_f32);
    let na_b = SVector::<f32, 4>::from_element(0.999_f32);

    c.bench_function("dot_4_aether_f32", |b| {
        b.iter(|| aether_dot4_simd_f32(black_box(&aether_a), black_box(&aether_b)))
    });

    c.bench_function("dot_4_scalar_f32", |b| {
        b.iter(|| aether_dot4_scalar_f32(black_box(&aether_a), black_box(&aether_b)))
    });

    c.bench_function("dot_4_nalgebra_f32", |b| {
        b.iter(|| nalgebra_dot4_f32(black_box(&na_a), black_box(&na_b)))
    });
}

fn bench_attitude_coordinate_abstractions(c: &mut Criterion) {
    type RF = Unknown;

    let cart_a = Cartesian::<f64, RF>::new(1.1, -2.2, 3.3);
    let cart_b = Cartesian::<f64, RF>::new(-0.9, 4.2, 0.7);
    let vec_a = cart_a.data;
    let vec_b = cart_b.data;

    let eul_a = Euler::<f64, RF, RF>::new(0.11, -0.22, 0.33);
    let eul_b = Euler::<f64, RF, RF>::new(-0.27, 0.19, -0.41);

    let dcm_a: DirectionCosineMatrix<f64, RF, RF> = DirectionCosineMatrix::from(eul_a);
    let dcm_b: DirectionCosineMatrix<f64, RF, RF> = DirectionCosineMatrix::from(eul_b);

    let mat_a = *dcm_a.as_matrix();
    let mat_b = *dcm_b.as_matrix();

    let quat_a: Quaternion<f64, RF, RF> = Quaternion::from(&eul_a);
    let quat_b: Quaternion<f64, RF, RF> = Quaternion::from(&eul_b);
    let quat_a_dcm: DirectionCosineMatrix<f64, RF, RF> = quat_a.to_dcm();
    let quat_b_dcm: DirectionCosineMatrix<f64, RF, RF> = quat_b.to_dcm();
    let quat_a_mat = *quat_a_dcm.as_matrix();
    let quat_b_mat = *quat_b_dcm.as_matrix();

    c.bench_function("cartesian_dot_abstraction", |bch| {
        bch.iter(|| black_box(black_box(&cart_a).dot(black_box(&cart_b))))
    });

    c.bench_function("cartesian_dot_native_vector", |bch| {
        bch.iter(|| black_box(native_vec3_dot(black_box(&vec_a), black_box(&vec_b))))
    });

    c.bench_function("cartesian_cross_abstraction", |bch| {
        bch.iter(|| {
            let out = black_box(&cart_a).cross(black_box(&cart_b));
            black_box(out.x())
        })
    });

    c.bench_function("cartesian_cross_native_vector", |bch| {
        bch.iter(|| {
            let out = native_vec3_cross(black_box(&vec_a), black_box(&vec_b));
            black_box(out.data[0])
        })
    });

    c.bench_function("dcm_rotate_cartesian_abstraction", |bch| {
        bch.iter(|| {
            let out = black_box(&dcm_a) * black_box(cart_a);
            black_box(out.x())
        })
    });

    c.bench_function("dcm_rotate_native_matrix_vector", |bch| {
        bch.iter(|| {
            let out = native_mat3_vec3(black_box(&mat_a), black_box(&vec_a));
            black_box(out.data[0])
        })
    });

    c.bench_function("dcm_compose_abstraction", |bch| {
        bch.iter(|| {
            let out = black_box(dcm_a) * black_box(dcm_b);
            black_box(out.m11())
        })
    });

    c.bench_function("dcm_compose_native_matrix", |bch| {
        bch.iter(|| {
            let out = native_mat3_mul_owned(black_box(mat_a), black_box(mat_b));
            black_box(out[(0, 0)])
        })
    });

    c.bench_function("quaternion_rotate_cartesian_abstraction", |bch| {
        bch.iter(|| {
            let out = black_box(&quat_a) * black_box(cart_a);
            black_box(out.x())
        })
    });

    c.bench_function("quaternion_rotate_native_matrix_vector", |bch| {
        bch.iter(|| {
            let out = native_quat_rotate_vector_equivalent(black_box(&quat_a), black_box(&vec_a));
            black_box(out.data[0])
        })
    });

    c.bench_function("quaternion_compose_abstraction", |bch| {
        bch.iter(|| {
            let out = black_box(&quat_a) * black_box(&quat_b);
            black_box(out.w())
        })
    });

    c.bench_function("quaternion_compose_native_matrix", |bch| {
        bch.iter(|| {
            let out = native_mat3_mul(black_box(&quat_a_mat), black_box(&quat_b_mat));
            black_box(out[(0, 0)])
        })
    });
}

fn bench_matvec_6_f32(c: &mut Criterion) {
    let aether_a = Matrix::new([[1.001_f32; 6]; 6]);
    let na_a = SMatrix::<f32, 6, 6>::from_element(1.001_f32);
    let aether_v = Vector::new([0.999_f32; 6]);
    let na_v = SVector::<f32, 6>::from_element(0.999_f32);

    c.bench_function("matvec_6_aether_f32", |b| {
        b.iter(|| aether_matvec6_simd_f32(black_box(&aether_a), black_box(&aether_v)))
    });

    c.bench_function("matvec_6_nalgebra_f32", |b| {
        b.iter(|| nalgebra_matvec_f32::<6>(black_box(&na_a), black_box(&na_v)))
    });
}

fn throughput_benches(c: &mut Criterion) {
    bench_matmul::<3>(c);
    bench_matmul::<4>(c);
    bench_matmul::<6>(c);
    bench_matmul_f32::<3>(c);
    bench_matmul_f32::<4>(c);
    bench_matmul_f32::<6>(c);

    bench_matvec::<3>(c);
    bench_matvec_4(c);
    bench_matvec_6(c);
    bench_matvec_f32::<3>(c);
    bench_matvec_4_f32(c);
    bench_matvec_6_f32(c);
    bench_dot_4(c);
    bench_dot_4_f32(c);
    bench_attitude_coordinate_abstractions(c);

    for n in [32, 64, 128, 256, 512, 1024] {
        bench_matvec_large(c, n);
        bench_matvec_large_f32(c, n);
    }

    for n in [32, 64, 128, 256] {
        bench_matmul_large(c, n);
        bench_matmul_large_f32(c, n);
    }
}

criterion_group!(benches, throughput_benches);
criterion_main!(benches);
