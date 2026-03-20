pub mod gaia_color {
    /// BP−RP -> effective temperature (K), rough MS prior.
    pub fn bp_rp_to_teff(bp_rp: f32) -> f32 {
        const KNOTS: &[(f32, f32)] = &[
            (-0.2, 15000.0),
            (0.0, 10000.0),
            (0.3, 8000.0),
            (0.65, 5800.0), // ~Sun
            (1.0, 5000.0),  // K
            (1.5, 4000.0),
            (2.0, 3500.0),
            (2.5, 3200.0),
            (3.0, 3000.0), // late-M
        ];
        let x = bp_rp.clamp(KNOTS[0].0, KNOTS[KNOTS.len() - 1].0);
        for w in KNOTS.windows(2) {
            let (x0, y0) = (w[0].0, w[0].1);
            let (x1, y1) = (w[1].0, w[1].1);
            if x <= x1 {
                let t = (x - x0) / (x1 - x0);
                return y0 + t * (y1 - y0);
            }
        }
        KNOTS[KNOTS.len() - 1].1
    }

    /// Blackbody Teff -> RGB [0..1] (approx along Planckian locus).
    pub fn blackbody_to_rgb(teff: f32) -> [f32; 3] {
        let t = teff.clamp(1000.0, 40000.0) / 100.0;

        let r = if t <= 66.0 {
            1.0
        } else {
            (329.698727446 * (t - 60.0).powf(-0.1332047592) / 255.0).clamp(0.0, 1.0)
        };

        let g = if t <= 66.0 {
            (99.4708025861 * t.ln() - 161.1195681661)
                .max(0.0)
                .min(255.0)
                / 255.0
        } else {
            (288.1221695283 * (t - 60.0).powf(-0.0755148492) / 255.0).clamp(0.0, 1.0)
        };

        let b = if t >= 66.0 {
            1.0
        } else if t <= 19.0 {
            0.0
        } else {
            (138.5177312231 * (t - 10.0).ln() - 305.0447927307)
                .max(0.0)
                .min(255.0)
                / 255.0
        };

        [r as f32, g as f32, b as f32]
    }

    /// Optional sRGB->linear (keep if your FB is linear and you *don’t* enable FRAMEBUFFER_SRGB).
    pub fn srgb_to_linear(c: [f32; 3]) -> [f32; 3] {
        fn f(u: f32) -> f32 {
            if u <= 0.04045 {
                u / 12.92
            } else {
                ((u + 0.055) / 1.055).powf(2.4)
            }
        }
        [f(c[0]), f(c[1]), f(c[2])]
    }

    /// Final: BP−RP -> RGB. Set `linear_out` based on your framebuffer.
    pub fn bp_rp_to_rgb(bp_rp: f32, linear_out: bool) -> [f32; 3] {
        let teff = bp_rp_to_teff(bp_rp);
        let mut c = blackbody_to_rgb(teff);

        // mild saturation boost so the Milky Way doesn’t go all pastel
        let sat = 1.15_f32;
        let avg = (c[0] + c[1] + c[2]) / 3.0;
        c = [
            (avg + (c[0] - avg) * sat).clamp(0.0, 1.0),
            (avg + (c[1] - avg) * sat).clamp(0.0, 1.0),
            (avg + (c[2] - avg) * sat).clamp(0.0, 1.0),
        ];

        if linear_out {
            srgb_to_linear(c)
        } else {
            c
        }
    }
}
