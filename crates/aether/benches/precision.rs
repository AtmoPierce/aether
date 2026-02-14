use aether_core::math::{Matrix, Vector};
use criterion::{criterion_group, criterion_main, Criterion};
use nalgebra::{DMatrix, DVector, SMatrix, SVector};

#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

fn max_abs_rel_slices(a: &[f64], b: &[f64]) -> (f64, f64) {
    let mut max_abs = 0.0;
    let mut max_rel = 0.0;

    for (&lhs, &rhs) in a.iter().zip(b.iter()) {
        let abs = (lhs - rhs).abs();
        let rel = abs / rhs.abs().max(1e-12);
        if abs > max_abs {
            max_abs = abs;
        }
        if rel > max_rel {
            max_rel = rel;
        }
    }

    (max_abs, max_rel)
}

fn max_abs_rel_slices_f32(a: &[f32], b: &[f32]) -> (f32, f32) {
    let mut max_abs = 0.0_f32;
    let mut max_rel = 0.0_f32;

    for (&lhs, &rhs) in a.iter().zip(b.iter()) {
        let abs = (lhs - rhs).abs();
        let rel = abs / rhs.abs().max(1e-12_f32);
        if abs > max_abs {
            max_abs = abs;
        }
        if rel > max_rel {
            max_rel = rel;
        }
    }

    (max_abs, max_rel)
}

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

fn report_precision(_c: &mut Criterion) {
    let aether3_a = Matrix::new([[1.001_f64; 3]; 3]);
    let aether3_b = Matrix::new([[0.999_f64; 3]; 3]);
    let na3_a = SMatrix::<f64, 3, 3>::from_element(1.001_f64);
    let na3_b = SMatrix::<f64, 3, 3>::from_element(0.999_f64);
    let out3_a = aether_matmul::<3>(&aether3_a, &aether3_b);
    let out3_n = nalgebra_matmul::<3>(&na3_a, &na3_b);
    let mut max_abs = 0.0;
    let mut max_rel = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            let abs = (out3_a.data[i][j] - out3_n[(i, j)]).abs();
            let rel = abs / out3_n[(i, j)].abs().max(1e-12);
            if abs > max_abs {
                max_abs = abs;
            }
            if rel > max_rel {
                max_rel = rel;
            }
        }
    }
    println!("precision matmul_3_aether vs nalgebra: max_abs={:.3e} max_rel={:.3e}", max_abs, max_rel);

    let aether4_a = Matrix::new([[1.001_f64; 4]; 4]);
    let aether4_b = Matrix::new([[0.999_f64; 4]; 4]);
    let na4_a = SMatrix::<f64, 4, 4>::from_element(1.001_f64);
    let na4_b = SMatrix::<f64, 4, 4>::from_element(0.999_f64);
    let out4_a = aether_matmul::<4>(&aether4_a, &aether4_b);
    let out4_n = nalgebra_matmul::<4>(&na4_a, &na4_b);
    let mut max_abs4 = 0.0;
    let mut max_rel4 = 0.0;
    for i in 0..4 {
        for j in 0..4 {
            let abs = (out4_a.data[i][j] - out4_n[(i, j)]).abs();
            let rel = abs / out4_n[(i, j)].abs().max(1e-12);
            if abs > max_abs4 {
                max_abs4 = abs;
            }
            if rel > max_rel4 {
                max_rel4 = rel;
            }
        }
    }
    println!("precision matmul_4_aether vs nalgebra: max_abs={:.3e} max_rel={:.3e}", max_abs4, max_rel4);

    let aether6_a = Matrix::new([[1.001_f64; 6]; 6]);
    let aether6_b = Matrix::new([[0.999_f64; 6]; 6]);
    let na6_a = SMatrix::<f64, 6, 6>::from_element(1.001_f64);
    let na6_b = SMatrix::<f64, 6, 6>::from_element(0.999_f64);
    let out6_a = aether_matmul::<6>(&aether6_a, &aether6_b);
    let out6_n = nalgebra_matmul::<6>(&na6_a, &na6_b);
    let mut max_abs6 = 0.0;
    let mut max_rel6 = 0.0;
    for i in 0..6 {
        for j in 0..6 {
            let abs = (out6_a.data[i][j] - out6_n[(i, j)]).abs();
            let rel = abs / out6_n[(i, j)].abs().max(1e-12);
            if abs > max_abs6 {
                max_abs6 = abs;
            }
            if rel > max_rel6 {
                max_rel6 = rel;
            }
        }
    }
    println!("precision matmul_6_aether vs nalgebra: max_abs={:.3e} max_rel={:.3e}", max_abs6, max_rel6);

    let aether_v3 = Vector::new([0.999_f64; 3]);
    let na_v3 = SVector::<f64, 3>::from_element(0.999_f64);
    let outv3_a = aether_matvec::<3>(&aether3_a, &aether_v3);
    let outv3_n = nalgebra_matvec::<3>(&na3_a, &na_v3);
    let (v3_abs, v3_rel) = max_abs_rel_slices(&outv3_a, &outv3_n);
    println!("precision matvec_3_aether vs nalgebra: max_abs={:.3e} max_rel={:.3e}", v3_abs, v3_rel);

    let aether_v4 = Vector::new([0.999_f64; 4]);
    let na_v4 = SVector::<f64, 4>::from_element(0.999_f64);
    let outv4_a = aether_matvec4_simd(&aether4_a, &aether_v4);
    let outv4_n = nalgebra_matvec::<4>(&na4_a, &na_v4);
    let (v4_abs, v4_rel) = max_abs_rel_slices(&outv4_a, &outv4_n);
    println!("precision matvec_4_aether vs nalgebra: max_abs={:.3e} max_rel={:.3e}", v4_abs, v4_rel);

    let dot_a4 = Vector::new([1.001_f64; 4]);
    let dot_b4 = Vector::new([0.999_f64; 4]);
    let na_dot_a4 = SVector::<f64, 4>::from_element(1.001_f64);
    let na_dot_b4 = SVector::<f64, 4>::from_element(0.999_f64);
    let aether_dot4 = aether_dot4_simd(&dot_a4, &dot_b4);
    let nalgebra_dot4_ref = nalgebra_dot4(&na_dot_a4, &na_dot_b4);
    let dot4_abs = (aether_dot4 - nalgebra_dot4_ref).abs();
    let dot4_rel = dot4_abs / nalgebra_dot4_ref.abs().max(1e-12);
    println!("precision dot_4_aether vs nalgebra: max_abs={:.3e} max_rel={:.3e}", dot4_abs, dot4_rel);

    let aether_dot4_scalar_v = aether_dot4_scalar(&dot_a4, &dot_b4);
    let dot4_scalar_abs = (aether_dot4_scalar_v - nalgebra_dot4_ref).abs();
    let dot4_scalar_rel = dot4_scalar_abs / nalgebra_dot4_ref.abs().max(1e-12);
    println!("precision dot_4_scalar vs nalgebra: max_abs={:.3e} max_rel={:.3e}", dot4_scalar_abs, dot4_scalar_rel);

    let aether_v6 = Vector::new([0.999_f64; 6]);
    let na_v6 = SVector::<f64, 6>::from_element(0.999_f64);
    let outv6_a = aether_matvec6_simd(&aether6_a, &aether_v6);
    let outv6_n = nalgebra_matvec::<6>(&na6_a, &na_v6);
    let (v6_abs, v6_rel) = max_abs_rel_slices(&outv6_a, &outv6_n);
    println!("precision matvec_6_aether vs nalgebra: max_abs={:.3e} max_rel={:.3e}", v6_abs, v6_rel);

    let aether3_a_f32 = Matrix::new([[1.001_f32; 3]; 3]);
    let aether3_b_f32 = Matrix::new([[0.999_f32; 3]; 3]);
    let na3_a_f32 = SMatrix::<f32, 3, 3>::from_element(1.001_f32);
    let na3_b_f32 = SMatrix::<f32, 3, 3>::from_element(0.999_f32);
    let out3_a_f32 = aether_matmul_f32::<3>(&aether3_a_f32, &aether3_b_f32);
    let out3_n_f32 = nalgebra_matmul_f32::<3>(&na3_a_f32, &na3_b_f32);
    let mut max_abs_f32 = 0.0_f32;
    let mut max_rel_f32 = 0.0_f32;
    for i in 0..3 {
        for j in 0..3 {
            let abs = (out3_a_f32.data[i][j] - out3_n_f32[(i, j)]).abs();
            let rel = abs / out3_n_f32[(i, j)].abs().max(1e-12_f32);
            if abs > max_abs_f32 {
                max_abs_f32 = abs;
            }
            if rel > max_rel_f32 {
                max_rel_f32 = rel;
            }
        }
    }
    println!("precision matmul_3_aether_f32 vs nalgebra_f32: max_abs={:.3e} max_rel={:.3e}", max_abs_f32, max_rel_f32);

    let aether4_a_f32 = Matrix::new([[1.001_f32; 4]; 4]);
    let aether4_b_f32 = Matrix::new([[0.999_f32; 4]; 4]);
    let na4_a_f32 = SMatrix::<f32, 4, 4>::from_element(1.001_f32);
    let na4_b_f32 = SMatrix::<f32, 4, 4>::from_element(0.999_f32);
    let out4_a_f32 = aether_matmul_f32::<4>(&aether4_a_f32, &aether4_b_f32);
    let out4_n_f32 = nalgebra_matmul_f32::<4>(&na4_a_f32, &na4_b_f32);
    let mut max_abs4_f32 = 0.0_f32;
    let mut max_rel4_f32 = 0.0_f32;
    for i in 0..4 {
        for j in 0..4 {
            let abs = (out4_a_f32.data[i][j] - out4_n_f32[(i, j)]).abs();
            let rel = abs / out4_n_f32[(i, j)].abs().max(1e-12_f32);
            if abs > max_abs4_f32 {
                max_abs4_f32 = abs;
            }
            if rel > max_rel4_f32 {
                max_rel4_f32 = rel;
            }
        }
    }
    println!("precision matmul_4_aether_f32 vs nalgebra_f32: max_abs={:.3e} max_rel={:.3e}", max_abs4_f32, max_rel4_f32);

    let aether6_a_f32 = Matrix::new([[1.001_f32; 6]; 6]);
    let aether6_b_f32 = Matrix::new([[0.999_f32; 6]; 6]);
    let na6_a_f32 = SMatrix::<f32, 6, 6>::from_element(1.001_f32);
    let na6_b_f32 = SMatrix::<f32, 6, 6>::from_element(0.999_f32);
    let out6_a_f32 = aether_matmul_f32::<6>(&aether6_a_f32, &aether6_b_f32);
    let out6_n_f32 = nalgebra_matmul_f32::<6>(&na6_a_f32, &na6_b_f32);
    let mut max_abs6_f32 = 0.0_f32;
    let mut max_rel6_f32 = 0.0_f32;
    for i in 0..6 {
        for j in 0..6 {
            let abs = (out6_a_f32.data[i][j] - out6_n_f32[(i, j)]).abs();
            let rel = abs / out6_n_f32[(i, j)].abs().max(1e-12_f32);
            if abs > max_abs6_f32 {
                max_abs6_f32 = abs;
            }
            if rel > max_rel6_f32 {
                max_rel6_f32 = rel;
            }
        }
    }
    println!("precision matmul_6_aether_f32 vs nalgebra_f32: max_abs={:.3e} max_rel={:.3e}", max_abs6_f32, max_rel6_f32);

    let aether_v3_f32 = Vector::new([0.999_f32; 3]);
    let na_v3_f32 = SVector::<f32, 3>::from_element(0.999_f32);
    let outv3_a_f32 = aether_matvec_f32::<3>(&aether3_a_f32, &aether_v3_f32);
    let outv3_n_f32 = nalgebra_matvec_f32::<3>(&na3_a_f32, &na_v3_f32);
    let (v3_abs_f32, v3_rel_f32) = max_abs_rel_slices_f32(&outv3_a_f32, &outv3_n_f32);
    println!("precision matvec_3_aether_f32 vs nalgebra_f32: max_abs={:.3e} max_rel={:.3e}", v3_abs_f32, v3_rel_f32);

    let aether_v4_f32 = Vector::new([0.999_f32; 4]);
    let na_v4_f32 = SVector::<f32, 4>::from_element(0.999_f32);
    let outv4_a_f32 = aether_matvec4_simd_f32(&aether4_a_f32, &aether_v4_f32);
    let outv4_n_f32 = nalgebra_matvec_f32::<4>(&na4_a_f32, &na_v4_f32);
    let (v4_abs_f32, v4_rel_f32) = max_abs_rel_slices_f32(&outv4_a_f32, &outv4_n_f32);
    println!("precision matvec_4_aether_f32 vs nalgebra_f32: max_abs={:.3e} max_rel={:.3e}", v4_abs_f32, v4_rel_f32);

    let dot_a4_f32 = Vector::new([1.001_f32; 4]);
    let dot_b4_f32 = Vector::new([0.999_f32; 4]);
    let na_dot_a4_f32 = SVector::<f32, 4>::from_element(1.001_f32);
    let na_dot_b4_f32 = SVector::<f32, 4>::from_element(0.999_f32);
    let aether_dot4_f32 = aether_dot4_simd_f32(&dot_a4_f32, &dot_b4_f32);
    let nalgebra_dot4_ref_f32 = nalgebra_dot4_f32(&na_dot_a4_f32, &na_dot_b4_f32);
    let dot4_abs_f32 = (aether_dot4_f32 - nalgebra_dot4_ref_f32).abs();
    let dot4_rel_f32 = dot4_abs_f32 / nalgebra_dot4_ref_f32.abs().max(1e-12_f32);
    println!("precision dot_4_aether_f32 vs nalgebra_f32: max_abs={:.3e} max_rel={:.3e}", dot4_abs_f32, dot4_rel_f32);

    let aether_dot4_scalar_f32_v = aether_dot4_scalar_f32(&dot_a4_f32, &dot_b4_f32);
    let dot4_scalar_abs_f32 = (aether_dot4_scalar_f32_v - nalgebra_dot4_ref_f32).abs();
    let dot4_scalar_rel_f32 = dot4_scalar_abs_f32 / nalgebra_dot4_ref_f32.abs().max(1e-12_f32);
    println!("precision dot_4_scalar_f32 vs nalgebra_f32: max_abs={:.3e} max_rel={:.3e}", dot4_scalar_abs_f32, dot4_scalar_rel_f32);

    let aether_v6_f32 = Vector::new([0.999_f32; 6]);
    let na_v6_f32 = SVector::<f32, 6>::from_element(0.999_f32);
    let outv6_a_f32 = aether_matvec6_simd_f32(&aether6_a_f32, &aether_v6_f32);
    let outv6_n_f32 = nalgebra_matvec_f32::<6>(&na6_a_f32, &na_v6_f32);
    let (v6_abs_f32, v6_rel_f32) = max_abs_rel_slices_f32(&outv6_a_f32, &outv6_n_f32);
    println!("precision matvec_6_aether_f32 vs nalgebra_f32: max_abs={:.3e} max_rel={:.3e}", v6_abs_f32, v6_rel_f32);

    for n in [64_usize, 128, 256, 512, 1024] {
        let mut a = vec![0.0_f64; n * n];
        for i in 0..n {
            for j in 0..n {
                a[i * n + j] = 1.001 + ((i + j) % 7) as f64 * 1e-3;
            }
        }
        let x = vec![0.999_f64; n];

        let mut y_naive = vec![0.0_f64; n];
        let mut y_native = vec![0.0_f64; n];
        let mut y_fma = vec![0.0_f64; n];
        let mut y_simd = vec![0.0_f64; n];

        let na_a = DMatrix::<f64>::from_row_slice(n, n, &a);
        let na_x = DVector::<f64>::from_row_slice(&x);
        let na_ref = &na_a * &na_x;

        matvec_naive_rowmajor(&a, &x, &mut y_native, n);
        let (native_abs, native_rel) = max_abs_rel_slices(&y_native, na_ref.as_slice());
        println!(
            "precision matvec_{}_native vs nalgebra: max_abs={:.3e} max_rel={:.3e}",
            n, native_abs, native_rel
        );

        matvec_naive_rowmajor(&a, &x, &mut y_naive, n);
        let (naive_abs, naive_rel) = max_abs_rel_slices(&y_naive, na_ref.as_slice());
        println!(
            "precision matvec_{}_naive vs nalgebra: max_abs={:.3e} max_rel={:.3e}",
            n, naive_abs, naive_rel
        );

        if avx_fma_available() {
            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
            unsafe {
                matvec_fma_avx_rowmajor(&a, &x, &mut y_fma, n);
            }
            let (fma_abs, fma_rel) = max_abs_rel_slices(&y_fma, na_ref.as_slice());
            println!(
                "precision matvec_{}_fma vs nalgebra: max_abs={:.3e} max_rel={:.3e}",
                n, fma_abs, fma_rel
            );
        }

        if avx_available() {
            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
            unsafe {
                matvec_simd_avx_rowmajor(&a, &x, &mut y_simd, n);
            }
        } else {
            matvec_naive_rowmajor(&a, &x, &mut y_simd, n);
        }
        let (simd_abs, simd_rel) = max_abs_rel_slices(&y_simd, na_ref.as_slice());
        println!(
            "precision matvec_{}_simd vs nalgebra: max_abs={:.3e} max_rel={:.3e}",
            n, simd_abs, simd_rel
        );

        let mut a_f32 = vec![0.0_f32; n * n];
        for i in 0..n {
            for j in 0..n {
                a_f32[i * n + j] = 1.001_f32 + ((i + j) % 7) as f32 * 1e-3_f32;
            }
        }
        let x_f32 = vec![0.999_f32; n];

        let mut y_naive_f32 = vec![0.0_f32; n];
        let mut y_native_f32 = vec![0.0_f32; n];
        let mut y_fma_f32 = vec![0.0_f32; n];
        let mut y_simd_f32 = vec![0.0_f32; n];

        let na_a_f32 = DMatrix::<f32>::from_row_slice(n, n, &a_f32);
        let na_x_f32 = DVector::<f32>::from_row_slice(&x_f32);
        let na_ref_f32 = &na_a_f32 * &na_x_f32;

        matvec_naive_rowmajor_f32(&a_f32, &x_f32, &mut y_native_f32, n);
        let (native_abs_f32, native_rel_f32) = max_abs_rel_slices_f32(&y_native_f32, na_ref_f32.as_slice());
        println!(
            "precision matvec_{}_native_f32 vs nalgebra_f32: max_abs={:.3e} max_rel={:.3e}",
            n, native_abs_f32, native_rel_f32
        );

        matvec_naive_rowmajor_f32(&a_f32, &x_f32, &mut y_naive_f32, n);
        let (naive_abs_f32, naive_rel_f32) = max_abs_rel_slices_f32(&y_naive_f32, na_ref_f32.as_slice());
        println!(
            "precision matvec_{}_naive_f32 vs nalgebra_f32: max_abs={:.3e} max_rel={:.3e}",
            n, naive_abs_f32, naive_rel_f32
        );

        if avx_fma_available() {
            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
            unsafe {
                matvec_fma_avx_rowmajor_f32(&a_f32, &x_f32, &mut y_fma_f32, n);
            }
            let (fma_abs_f32, fma_rel_f32) = max_abs_rel_slices_f32(&y_fma_f32, na_ref_f32.as_slice());
            println!(
                "precision matvec_{}_fma_f32 vs nalgebra_f32: max_abs={:.3e} max_rel={:.3e}",
                n, fma_abs_f32, fma_rel_f32
            );
        }

        if avx_available() {
            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
            unsafe {
                matvec_simd_avx_rowmajor_f32(&a_f32, &x_f32, &mut y_simd_f32, n);
            }
        } else {
            matvec_naive_rowmajor_f32(&a_f32, &x_f32, &mut y_simd_f32, n);
        }
        let (simd_abs_f32, simd_rel_f32) = max_abs_rel_slices_f32(&y_simd_f32, na_ref_f32.as_slice());
        println!(
            "precision matvec_{}_simd_f32 vs nalgebra_f32: max_abs={:.3e} max_rel={:.3e}",
            n, simd_abs_f32, simd_rel_f32
        );
    }
}

criterion_group!(benches, report_precision);
criterion_main!(benches);
