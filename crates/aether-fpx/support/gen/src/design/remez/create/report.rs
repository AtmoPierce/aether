pub fn print_ulp_report_for_plan(
    plan: &crate::design::remez::create::coefficients::FunctionPlan,
    f_x: fn(f128) -> f128,
    samples: usize,
    label: &str,
) {
    use crate::design::remez::create::coefficients::{eval_plan_f32, eval_plan_f64};
    use aether_fpx_bench::ulp::{ulp_dist_f32, ulp_dist_f64};

    #[derive(Copy, Clone, Debug)]
    struct Row64 {
        x: f64,
        pq: f64,
        fx: f64,
        ulp: u64,
        abs: f64,
        seg: usize,
    }

    #[derive(Copy, Clone, Debug)]
    struct Row32 {
        x: f32,
        pq: f32,
        fx: f32,
        ulp: u32,
        abs: f32,
        seg: usize,
    }

    let mut worst64 = Row64 {
        x: 0.0,
        pq: 0.0,
        fx: 0.0,
        ulp: 0,
        abs: 0.0,
        seg: 0,
    };
    let mut worst32 = Row32 {
        x: 0.0,
        pq: 0.0,
        fx: 0.0,
        ulp: 0,
        abs: 0.0,
        seg: 0,
    };

    const TOP: usize = 12;
    let mut top64 = [worst64; TOP];
    let mut top32 = [worst32; TOP];

    let mut push_top64 = |r: Row64| {
        for i in 0..TOP {
            if r.ulp > top64[i].ulp {
                for j in (i + 1..TOP).rev() {
                    top64[j] = top64[j - 1];
                }
                top64[i] = r;
                break;
            }
        }
    };

    let mut push_top32 = |r: Row32| {
        for i in 0..TOP {
            if r.ulp > top32[i].ulp {
                for j in (i + 1..TOP).rev() {
                    top32[j] = top32[j - 1];
                }
                top32[i] = r;
                break;
            }
        }
    };

    // Sample each segment uniformly in x-domain (exact runtime behavior)
    for (si, seg) in plan.segments.iter().enumerate() {
        let a64 = if si == 0 {
            0.0
        } else {
            plan.segments[si - 1].b64
        };
        let b64 = seg.b64;

        for i in 0..=samples {
            let t = (i as f64) / (samples as f64);
            let x = a64 + (b64 - a64) * t;

            // f64 path
            let pq = eval_plan_f64(plan, x);
            let fx = f_x(x as f128) as f64;
            let ulp = ulp_dist_f64(pq, fx);
            let abs = (pq - fx).abs();

            let r64 = Row64 {
                x,
                pq,
                fx,
                ulp,
                abs,
                seg: si,
            };
            if r64.ulp > worst64.ulp {
                worst64 = r64;
            }
            push_top64(r64);

            // f32 path
            let x32 = x as f32;
            let pq32 = eval_plan_f32(plan, x32);
            let fx32 = f_x(x32 as f128) as f32;
            let ulp32 = ulp_dist_f32(pq32, fx32);
            let abs32 = (pq32 - fx32).abs();

            let r32 = Row32 {
                x: x32,
                pq: pq32,
                fx: fx32,
                ulp: ulp32,
                abs: abs32,
                seg: si,
            };
            if r32.ulp > worst32.ulp {
                worst32 = r32;
            }
            push_top32(r32);
        }
    }

    println!("\n[aether-fpx-gen][ulp-report] {label}");
    println!(
        "  segments={}  samples/seg={}",
        plan.segments.len(),
        samples
    );

    println!(
        "\n  WORST f64: ulp={}  seg={}  x={:.17}",
        worst64.ulp, worst64.seg, worst64.x
    );
    println!("    pq={:.17e}", worst64.pq);
    println!("    fx={:.17e}", worst64.fx);
    println!("    abs_err={:.3e}", worst64.abs);

    println!(
        "\n  WORST f32: ulp={}  seg={}  x={:.9}",
        worst32.ulp, worst32.seg, worst32.x
    );
    println!("    pq={:.9e}", worst32.pq);
    println!("    fx={:.9e}", worst32.fx);
    println!("    abs_err={:.3e}", worst32.abs);

    println!("\n  Top f64 (by ULP):");
    println!(
        "    {:>3}  {:>4}  {:>22}  {:>6}  {:>12}",
        "#", "seg", "x", "ulp", "abs_err"
    );
    for (i, r) in top64.iter().enumerate() {
        println!(
            "    {:>3}  {:>4}  {:>22.17}  {:>6}  {:>12.3e}",
            i, r.seg, r.x, r.ulp, r.abs
        );
    }

    println!("\n  Top f32 (by ULP):");
    println!(
        "    {:>3}  {:>4}  {:>14}  {:>6}  {:>12}",
        "#", "seg", "x", "ulp", "abs_err"
    );
    for (i, r) in top32.iter().enumerate() {
        println!(
            "    {:>3}  {:>4}  {:>14.9}  {:>6}  {:>12.3e}",
            i, r.seg, r.x, r.ulp, r.abs
        );
    }
}
