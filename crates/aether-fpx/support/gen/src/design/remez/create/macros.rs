/* ------------------------------------------------------------------------------------------------
 * Macros
 * ------------------------------------------------------------------------------------------------ */

#[macro_export]
macro_rules! remez_create {
    (
        name        = $name:ident,
        N           = $N:expr,
        K           = $K:expr,
        func        = $f:path,
        deriv       = $df:path,
        interval    = [$a:expr, $b:expr],
        max_error   = $max_err:expr,
        samples     = $samples:expr,
        title       = $title:expr,
        file_prefix = $prefix:expr
        $(, boundary_tol = $boundary_tol:expr, enforce_endpoints = $enforce_endpoints:expr)?
    ) => {
        pub fn $name() {
            type Rz = $crate::RemezApprox<f128, $N, $K>;

            let a: f128 = ($a) as f128;
            let b: f128 = ($b) as f128;
            let a64 = a as f64;
            let b64 = b as f64;

            let mut rz = Rz::new(a, b);
            rz.solve($f, $df);

            // ULP acceptance in every emitted IEEE format. Keep the f128
            // interval and coefficient evaluation native to f128.
            let coeffs_f128: [f128; $N] = rz.coefficients;
            let (_ulp_poly128, ulp_func128) =
                $crate::design::remez::create::metrics::max_ulps_for_coeffs_f128(
                    &rz, &coeffs_f128, &$f, a, b, $samples
                );

            assert!(
                ulp_func128 <= $crate::design::remez::create::metrics::MAX_ULP_ALLOWED_F128,
                "max ulp(f128) = {} over [{:?},{:?}] (allowed {})",
                ulp_func128, a, b,
                $crate::design::remez::create::metrics::MAX_ULP_ALLOWED_F128
            );

            let coeffs_f64: [f64; $N] = core::array::from_fn(|i| rz.coefficients[i] as f64);
            let (_ulp_poly64, ulp_func64) =
                $crate::design::remez::create::accept::max_ulps_for_coeffs_f64(
                    &rz, &coeffs_f64, &$f, a, b, $samples
                );

            assert!(
                ulp_func64 <= $crate::design::remez::create::accept::MAX_ULP_ALLOWED_F64,
                "max ulp(f64) = {} over [{:.6},{:.6}] (allowed {})",
                ulp_func64, a64, b64, $crate::design::remez::create::accept::MAX_ULP_ALLOWED_F64
            );

            let coeffs_f32: [f32; $N] = core::array::from_fn(|i| rz.coefficients[i] as f32);
            let (_ulp_poly32, ulp_func32) =
                $crate::design::remez::create::accept::max_ulps_for_coeffs_f32(
                    &rz, &coeffs_f32, &$f, a, b, $samples
                );

            assert!(
                ulp_func32 <= $crate::design::remez::create::accept::MAX_ULP_ALLOWED_F32,
                "max ulp(f32) = {} over [{:.6},{:.6}] (allowed {})",
                ulp_func32, a64, b64, $crate::design::remez::create::accept::MAX_ULP_ALLOWED_F32
            );

            // abs max error for reporting / exports
            let max_e128 = $crate::max_error_on_grid(&rz, $f, a, b, $samples);

            $(
            if $enforce_endpoints {
                let fa128 = $f(a);
                let fb128 = $f(b);
                let pa128 = rz.eval_poly(a);
                let pb128 = rz.eval_poly(b);

                let ea128 = (pa128 - fa128).abs();
                let eb128 = (pb128 - fb128).abs();

                let tol128: f128 = $boundary_tol as f128;

                assert!(
                    ea128 <= tol128 && eb128 <= tol128,
                    "boundary error too large: |P(a)-f(a)|={:.3e}, |P(b)-f(b)|={:.3e}, tol={:.3e}",
                    ea128 as f64, eb128 as f64, tol128 as f64
                );
            }
            )?

            $crate::design::remez::create::export::write_coeffs_to_json::<$N, $K, _>(
                &rz, stringify!($f), $f, (a, b), max_e128, $prefix
            );
            $crate::design::remez::create::export::write_coeffs_to_c::<$N, $K, _>(
                &rz, stringify!($f), $f, (a, b), max_e128, $prefix
            );
            $crate::design::remez::create::export::write_coeffs_to_rs::<$N, $K, _>(
                &rz, stringify!($f), $f, (a, b), max_e128, $prefix
            );

            // Call plotter
            {
                $crate::design::remez::create::plot::remez_generate_plots::<$N, $K, _>(
                    &rz, $f, a, b, $samples, $title, $prefix
                );
            }
        }
    };
}

#[macro_export]
macro_rules! remez_create_opt {
    (
        name        = $name:ident,
        N           = $N:expr,
        K           = $K:expr,
        func        = $f:path,
        deriv       = $df:path,
        interval    = [$a:expr, $b:expr],
        max_error   = $max_err:expr,
        samples     = $samples:expr,
        title       = $title:expr,
        file_prefix = $prefix:expr
        $(, boundary_tol = $boundary_tol:expr, enforce_endpoints = $enforce_endpoints:expr)?
    ) => {
        $crate::remez_create! {
            name        = $name,
            N           = $N,
            K           = $K,
            func        = $f,
            deriv       = $df,
            interval    = [$a, $b],
            max_error   = $max_err,
            samples     = $samples,
            title       = $title,
            file_prefix = $prefix
            $(, boundary_tol = $boundary_tol,
               enforce_endpoints = $enforce_endpoints)?
        }
    };
}

#[macro_export]
macro_rules! remez_piecewise {
    (
        name  = $name:ident,
        func  = $f:path,
        deriv = $df:path,
        pieces = [
            $(
                {
                    seg_name     = $seg_name:ident,
                    N            = $N:expr,
                    K            = $K:expr,
                    interval     = [$a:expr, $b:expr],
                    max_error    = $max_err:expr,
                    samples      = $samples:expr,
                    file_prefix  = $prefix:expr
                    $(, boundary_tol = $boundary_tol:expr,
                       enforce_endpoints = $enforce_endpoints:expr)?
                }
            ),+ $(,)?
        ]
    ) => {
        $(
            $crate::remez_create_opt! {
                name        = $seg_name,
                N           = $N,
                K           = $K,
                func        = $f,
                deriv       = $df,
                interval    = [$a, $b],
                max_error   = $max_err,
                samples     = $samples,
                title       = concat!(stringify!($f), " on [", stringify!($a), ", ", stringify!($b), "]"),
                file_prefix = $prefix
                $(, boundary_tol = $boundary_tol,
                   enforce_endpoints = $enforce_endpoints)?
            }
        )+

        pub fn $name() {
            $(
                $seg_name();
            )+
        }
    };
}

#[macro_export]
macro_rules! remez_create_quadrant_kernels {
    (
        name        = $name:ident,
        sin_N       = $NSIN:expr,
        sin_K       = $KSIN:expr,
        cos_N       = $NCOS:expr,
        cos_K       = $KCOS:expr,
        interval    = [$a:expr, $b:expr],
        samples     = $samples:expr,
        file_prefix = $prefix:expr
    ) => {
        pub fn $name(out_dir: &std::path::Path) {
            {
                std::env::set_var("FPX_REMEZ_OUT_DIR", out_dir);

                $crate::design::remez::create::quadrant::write_quadrant_kernels_to_rs::<
                    $NSIN,
                    $KSIN,
                    $NCOS,
                    $KCOS,
                >(($a, $b), $samples, $prefix);
            }
        }
    };
}

// plotting
use aether_viz::*;
#[macro_export]
macro_rules! plots { ($($tt:tt)*) => { $($tt)* } }
