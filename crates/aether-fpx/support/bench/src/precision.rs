use core::f64::consts::PI;
use core::fmt::Write;
use core::{f128, f16};
use rayon::prelude::*;
use std::cmp::Ordering;
use std::fs;
use std::path::Path;

use crate::ulp::{valid_ulp_dist_f128, valid_ulp_dist_f16, valid_ulp_dist_f32, valid_ulp_dist_f64};

#[derive(Clone, Debug)]
pub struct PrecisionSummary {
    pub name: String,
    pub precision: String,
    pub test: String,
    pub max_abs: String,
    pub max_ulp: String,
    pub max_abs_x: String,
    pub max_abs_ref: String,
    pub max_abs_fpx: String,
    pub max_ulp_x: String,
    pub max_ulp_ref: String,
    pub max_ulp_fpx: String,
}

#[derive(Clone, Debug)]
pub struct PrecisionReport {
    pub summary: PrecisionSummary,
    pub report: String,
}

impl PrecisionReport {
    pub fn print(&self) {
        print!("{}", self.report);
    }
}

fn format_optional_x_f32(x: Option<f32>) -> String {
    x.map(|value| format!("{value:.9e}"))
        .unwrap_or_else(|| "NA".to_string())
}

fn format_optional_x_f16(x: Option<f16>) -> String {
    x.map(|value| format!("{:.5e}", value as f32))
        .unwrap_or_else(|| "NA".to_string())
}

fn format_optional_x_f64(x: Option<f64>) -> String {
    x.map(|value| format!("{value:.17e}"))
        .unwrap_or_else(|| "NA".to_string())
}

fn format_optional_x_f128(x: Option<f128>) -> String {
    x.map(|value| format!("{value:?}"))
        .unwrap_or_else(|| "NA".to_string())
}

fn format_value_f32(x: f32) -> String {
    format!("{x:.9e}")
}

fn format_value_f16(x: f16) -> String {
    format!("{:.5e}", x as f32)
}

fn format_value_f64(x: f64) -> String {
    format!("{x:.17e}")
}

fn format_value_f128(x: f128) -> String {
    format!("{x:?}")
}

fn format_optional_value_f32(x: Option<f32>) -> String {
    x.map(format_value_f32).unwrap_or_else(|| "NA".to_string())
}

fn format_optional_value_f16(x: Option<f16>) -> String {
    x.map(format_value_f16).unwrap_or_else(|| "NA".to_string())
}

fn format_optional_value_f64(x: Option<f64>) -> String {
    x.map(format_value_f64).unwrap_or_else(|| "NA".to_string())
}

fn format_optional_value_f128(x: Option<f128>) -> String {
    x.map(format_value_f128).unwrap_or_else(|| "NA".to_string())
}

fn parse_max_ulp(value: &str) -> Option<u128> {
    if value == "NA" {
        None
    } else {
        value.parse::<u128>().ok()
    }
}

pub fn write_precision_csv<P: AsRef<Path>>(
    path: P,
    rows: &[PrecisionSummary],
) -> Result<(), String> {
    let mut out = String::new();
    out.push_str(
        "name,precision,test,max_abs,max_ulp,max_abs_x,max_abs_ref,max_abs_fpx,max_ulp_x,max_ulp_ref,max_ulp_fpx\n",
    );

    let mut sorted_rows: Vec<&PrecisionSummary> = rows.iter().collect();
    sorted_rows.sort_by(|left, right| {
        match (parse_max_ulp(&left.max_ulp), parse_max_ulp(&right.max_ulp)) {
            (Some(left_ulp), Some(right_ulp)) => right_ulp
                .cmp(&left_ulp)
                .then_with(|| left.precision.cmp(&right.precision))
                .then_with(|| left.name.cmp(&right.name))
                .then_with(|| left.test.cmp(&right.test)),
            (Some(_), None) => Ordering::Less,
            (None, Some(_)) => Ordering::Greater,
            (None, None) => left
                .precision
                .cmp(&right.precision)
                .then_with(|| left.name.cmp(&right.name))
                .then_with(|| left.test.cmp(&right.test)),
        }
    });

    for r in sorted_rows {
        out.push_str(&format!(
            "{},{},{},{},{},{},{},{},{},{},{}\n",
            r.name,
            r.precision,
            r.test,
            r.max_abs,
            r.max_ulp,
            r.max_abs_x,
            r.max_abs_ref,
            r.max_abs_fpx,
            r.max_ulp_x,
            r.max_ulp_ref,
            r.max_ulp_fpx,
        ));
    }

    fs::write(path, out).map_err(|e| e.to_string())
}

/* ============================================================
 * Ordered-float utilities
 *
 * These map IEEE floats into monotonically increasing unsigned
 * integers. They are suitable for exhaustive numeric sweeps over
 * intervals that may be negative or cross zero.
 * ============================================================ */

#[inline]
fn ordered_f32_bits(x: f32) -> u32 {
    let bits = x.to_bits();
    if (bits >> 31) != 0 {
        !bits
    } else {
        bits | (1u32 << 31)
    }
}

#[inline]
fn f32_from_ordered_bits(ord: u32) -> f32 {
    let bits = if (ord >> 31) != 0 {
        ord & !(1u32 << 31)
    } else {
        !ord
    };
    f32::from_bits(bits)
}

#[inline]
fn ordered_f64_bits(x: f64) -> u64 {
    let bits = x.to_bits();
    if (bits >> 63) != 0 {
        !bits
    } else {
        bits | (1u64 << 63)
    }
}

#[inline]
fn f64_from_ordered_bits(ord: u64) -> f64 {
    let bits = if (ord >> 63) != 0 {
        ord & !(1u64 << 63)
    } else {
        !ord
    };
    f64::from_bits(bits)
}

#[inline]
fn ordered_f16_bits(x: f16) -> u16 {
    let bits = x.to_bits();
    if (bits >> 15) != 0 {
        !bits
    } else {
        bits | (1u16 << 15)
    }
}

#[inline]
fn f16_from_ordered_bits(ord: u16) -> f16 {
    let bits = if (ord >> 15) != 0 {
        ord & !(1u16 << 15)
    } else {
        !ord
    };
    f16::from_bits(bits)
}

#[inline]
fn ordered_f128_bits(x: f128) -> u128 {
    let bits = x.to_bits();
    if (bits >> 127) != 0 {
        !bits
    } else {
        bits | (1u128 << 127)
    }
}

#[inline]
fn f128_from_ordered_bits(ord: u128) -> f128 {
    let bits = if (ord >> 127) != 0 {
        ord & !(1u128 << 127)
    } else {
        !ord
    };
    f128::from_bits(bits)
}

#[inline]
fn ordered_range_f32(lo: f32, hi: f32) -> Option<(u32, u32)> {
    if !lo.is_finite() || !hi.is_finite() {
        return None;
    }
    let a = ordered_f32_bits(lo);
    let b = ordered_f32_bits(hi);
    Some((a.min(b), a.max(b)))
}

#[inline]
fn ordered_range_f16(lo: f16, hi: f16) -> Option<(u16, u16)> {
    if !lo.is_finite() || !hi.is_finite() {
        return None;
    }
    let a = ordered_f16_bits(lo);
    let b = ordered_f16_bits(hi);
    Some((a.min(b), a.max(b)))
}

#[inline]
fn finite_ulp_dist_f32(a: f32, b: f32) -> Option<u32> {
    if a.is_finite() && b.is_finite() {
        valid_ulp_dist_f32(a, b)
    } else {
        None
    }
}

#[inline]
fn finite_ulp_dist_f16(a: f16, b: f16) -> Option<u16> {
    if a.is_finite() && b.is_finite() {
        valid_ulp_dist_f16(a, b)
    } else {
        None
    }
}

#[inline]
fn finite_ulp_dist_f64(a: f64, b: f64) -> Option<u64> {
    if a.is_finite() && b.is_finite() {
        valid_ulp_dist_f64(a, b)
    } else {
        None
    }
}

#[inline]
fn finite_ulp_dist_f128(a: f128, b: f128) -> Option<u128> {
    if a.is_finite() && b.is_finite() {
        valid_ulp_dist_f128(a, b)
    } else {
        None
    }
}

#[inline]
fn f32_neighbors(x: f32) -> Vec<f32> {
    if !x.is_finite() {
        return Vec::new();
    }

    let ord = ordered_f32_bits(x);
    let mut out = Vec::with_capacity(3);

    if ord > 0 {
        let v = f32_from_ordered_bits(ord - 1);
        if v.is_finite() {
            out.push(v);
        }
    }

    out.push(x);

    if ord < u32::MAX {
        let v = f32_from_ordered_bits(ord + 1);
        if v.is_finite() {
            out.push(v);
        }
    }

    out
}

#[inline]
fn f64_neighbors(x: f64) -> Vec<f64> {
    if !x.is_finite() {
        return Vec::new();
    }

    let ord = ordered_f64_bits(x);
    let mut out = Vec::with_capacity(3);

    if ord > 0 {
        let v = f64_from_ordered_bits(ord - 1);
        if v.is_finite() {
            out.push(v);
        }
    }

    out.push(x);

    if ord < u64::MAX {
        let v = f64_from_ordered_bits(ord + 1);
        if v.is_finite() {
            out.push(v);
        }
    }

    out
}

#[inline]
fn f16_neighbors(x: f16) -> Vec<f16> {
    if !x.is_finite() {
        return Vec::new();
    }

    let ord = ordered_f16_bits(x);
    let mut out = Vec::with_capacity(3);

    if ord > 0 {
        let v = f16_from_ordered_bits(ord - 1);
        if v.is_finite() {
            out.push(v);
        }
    }

    out.push(x);

    if ord < u16::MAX {
        let v = f16_from_ordered_bits(ord + 1);
        if v.is_finite() {
            out.push(v);
        }
    }

    out
}

#[inline]
fn f128_neighbors(x: f128) -> Vec<f128> {
    if !x.is_finite() {
        return Vec::new();
    }

    let ord = ordered_f128_bits(x);
    let mut out = Vec::with_capacity(3);

    if ord > 0 {
        let v = f128_from_ordered_bits(ord - 1);
        if v.is_finite() {
            out.push(v);
        }
    }

    out.push(x);

    if ord < u128::MAX {
        let v = f128_from_ordered_bits(ord + 1);
        if v.is_finite() {
            out.push(v);
        }
    }

    out
}

/* ============================================================
 * Sweep states
 * ============================================================ */

#[derive(Clone, Copy)]
struct SweepStateF32 {
    max_abs: f32,
    max_abs_x: f32,
    max_abs_ref: f32,
    max_abs_fpx: f32,
    max_ulp: Option<u32>,
    max_ulp_x: Option<f32>,
    max_ulp_ref: Option<f32>,
    max_ulp_fpx: Option<f32>,
}

impl SweepStateF32 {
    #[inline]
    fn new() -> Self {
        Self {
            max_abs: 0.0,
            max_abs_x: 0.0,
            max_abs_ref: 0.0,
            max_abs_fpx: 0.0,
            max_ulp: None,
            max_ulp_x: None,
            max_ulp_ref: None,
            max_ulp_fpx: None,
        }
    }

    #[inline]
    fn observe(&mut self, x: f32, y0: f32, y1: f32) {
        if !x.is_finite() || !y0.is_finite() || !y1.is_finite() {
            return;
        }

        let abs = (y0 - y1).abs();
        if abs > self.max_abs
            || (abs == self.max_abs && ordered_f32_bits(x) < ordered_f32_bits(self.max_abs_x))
        {
            self.max_abs = abs;
            self.max_abs_x = x;
            self.max_abs_ref = y0;
            self.max_abs_fpx = y1;
        }

        if let Some(ulp) = finite_ulp_dist_f32(y0, y1) {
            let should_update = match (self.max_ulp, self.max_ulp_x) {
                (None, _) => true,
                (Some(current), Some(current_x)) => {
                    ulp > current
                        || (ulp == current && ordered_f32_bits(x) < ordered_f32_bits(current_x))
                }
                (Some(current), None) => ulp > current,
            };
            if should_update {
                self.max_ulp = Some(ulp);
                self.max_ulp_x = Some(x);
                self.max_ulp_ref = Some(y0);
                self.max_ulp_fpx = Some(y1);
            }
        }
    }

    #[inline]
    fn merge(self, other: Self) -> Self {
        let mut merged = self;
        if other.max_abs > merged.max_abs
            || (other.max_abs == merged.max_abs
                && ordered_f32_bits(other.max_abs_x) < ordered_f32_bits(merged.max_abs_x))
        {
            merged.max_abs = other.max_abs;
            merged.max_abs_x = other.max_abs_x;
            merged.max_abs_ref = other.max_abs_ref;
            merged.max_abs_fpx = other.max_abs_fpx;
        }

        let other_is_better = match (
            merged.max_ulp,
            merged.max_ulp_x,
            other.max_ulp,
            other.max_ulp_x,
        ) {
            (_, _, None, _) => false,
            (None, _, Some(_), _) => true,
            (Some(current), Some(current_x), Some(candidate), Some(candidate_x)) => {
                candidate > current
                    || (candidate == current
                        && ordered_f32_bits(candidate_x) < ordered_f32_bits(current_x))
            }
            (Some(current), None, Some(candidate), _) => candidate > current,
            _ => false,
        };
        if other_is_better {
            merged.max_ulp = other.max_ulp;
            merged.max_ulp_x = other.max_ulp_x;
            merged.max_ulp_ref = other.max_ulp_ref;
            merged.max_ulp_fpx = other.max_ulp_fpx;
        }
        merged
    }
}

#[derive(Clone, Copy)]
struct SweepStateF16 {
    max_abs: f16,
    max_abs_x: f16,
    max_abs_ref: f16,
    max_abs_fpx: f16,
    max_ulp: Option<u16>,
    max_ulp_x: Option<f16>,
    max_ulp_ref: Option<f16>,
    max_ulp_fpx: Option<f16>,
}

impl SweepStateF16 {
    #[inline]
    fn new() -> Self {
        let zero = f16::from_bits(0);
        Self {
            max_abs: zero,
            max_abs_x: zero,
            max_abs_ref: zero,
            max_abs_fpx: zero,
            max_ulp: None,
            max_ulp_x: None,
            max_ulp_ref: None,
            max_ulp_fpx: None,
        }
    }

    #[inline]
    fn observe(&mut self, x: f16, y0: f16, y1: f16) {
        if !x.is_finite() || !y0.is_finite() || !y1.is_finite() {
            return;
        }

        let abs = (y0 - y1).abs();
        if abs > self.max_abs
            || (abs == self.max_abs && ordered_f16_bits(x) < ordered_f16_bits(self.max_abs_x))
        {
            self.max_abs = abs;
            self.max_abs_x = x;
            self.max_abs_ref = y0;
            self.max_abs_fpx = y1;
        }

        if let Some(ulp) = finite_ulp_dist_f16(y0, y1) {
            let should_update = match (self.max_ulp, self.max_ulp_x) {
                (None, _) => true,
                (Some(current), Some(current_x)) => {
                    ulp > current
                        || (ulp == current && ordered_f16_bits(x) < ordered_f16_bits(current_x))
                }
                (Some(current), None) => ulp > current,
            };
            if should_update {
                self.max_ulp = Some(ulp);
                self.max_ulp_x = Some(x);
                self.max_ulp_ref = Some(y0);
                self.max_ulp_fpx = Some(y1);
            }
        }
    }

    #[inline]
    fn merge(self, other: Self) -> Self {
        let mut merged = self;
        if other.max_abs > merged.max_abs
            || (other.max_abs == merged.max_abs
                && ordered_f16_bits(other.max_abs_x) < ordered_f16_bits(merged.max_abs_x))
        {
            merged.max_abs = other.max_abs;
            merged.max_abs_x = other.max_abs_x;
            merged.max_abs_ref = other.max_abs_ref;
            merged.max_abs_fpx = other.max_abs_fpx;
        }

        let other_is_better = match (
            merged.max_ulp,
            merged.max_ulp_x,
            other.max_ulp,
            other.max_ulp_x,
        ) {
            (_, _, None, _) => false,
            (None, _, Some(_), _) => true,
            (Some(current), Some(current_x), Some(candidate), Some(candidate_x)) => {
                candidate > current
                    || (candidate == current
                        && ordered_f16_bits(candidate_x) < ordered_f16_bits(current_x))
            }
            (Some(current), None, Some(candidate), _) => candidate > current,
            _ => false,
        };
        if other_is_better {
            merged.max_ulp = other.max_ulp;
            merged.max_ulp_x = other.max_ulp_x;
            merged.max_ulp_ref = other.max_ulp_ref;
            merged.max_ulp_fpx = other.max_ulp_fpx;
        }
        merged
    }
}

#[derive(Clone, Copy)]
struct SweepStateF64 {
    max_abs: f64,
    max_abs_x: f64,
    max_abs_ref: f64,
    max_abs_fpx: f64,
    max_ulp: Option<u64>,
    max_ulp_x: Option<f64>,
    max_ulp_ref: Option<f64>,
    max_ulp_fpx: Option<f64>,
}

impl SweepStateF64 {
    #[inline]
    fn new() -> Self {
        Self {
            max_abs: 0.0,
            max_abs_x: 0.0,
            max_abs_ref: 0.0,
            max_abs_fpx: 0.0,
            max_ulp: None,
            max_ulp_x: None,
            max_ulp_ref: None,
            max_ulp_fpx: None,
        }
    }

    #[inline]
    fn observe(&mut self, x: f64, y0: f64, y1: f64) {
        if !x.is_finite() || !y0.is_finite() || !y1.is_finite() {
            return;
        }

        let abs = (y0 - y1).abs();
        if abs > self.max_abs
            || (abs == self.max_abs && ordered_f64_bits(x) < ordered_f64_bits(self.max_abs_x))
        {
            self.max_abs = abs;
            self.max_abs_x = x;
            self.max_abs_ref = y0;
            self.max_abs_fpx = y1;
        }

        if let Some(ulp) = finite_ulp_dist_f64(y0, y1) {
            let should_update = match (self.max_ulp, self.max_ulp_x) {
                (None, _) => true,
                (Some(current), Some(current_x)) => {
                    ulp > current
                        || (ulp == current && ordered_f64_bits(x) < ordered_f64_bits(current_x))
                }
                (Some(current), None) => ulp > current,
            };
            if should_update {
                self.max_ulp = Some(ulp);
                self.max_ulp_x = Some(x);
                self.max_ulp_ref = Some(y0);
                self.max_ulp_fpx = Some(y1);
            }
        }
    }

    #[inline]
    fn merge(self, other: Self) -> Self {
        let mut merged = self;
        if other.max_abs > merged.max_abs
            || (other.max_abs == merged.max_abs
                && ordered_f64_bits(other.max_abs_x) < ordered_f64_bits(merged.max_abs_x))
        {
            merged.max_abs = other.max_abs;
            merged.max_abs_x = other.max_abs_x;
            merged.max_abs_ref = other.max_abs_ref;
            merged.max_abs_fpx = other.max_abs_fpx;
        }

        let other_is_better = match (
            merged.max_ulp,
            merged.max_ulp_x,
            other.max_ulp,
            other.max_ulp_x,
        ) {
            (_, _, None, _) => false,
            (None, _, Some(_), _) => true,
            (Some(current), Some(current_x), Some(candidate), Some(candidate_x)) => {
                candidate > current
                    || (candidate == current
                        && ordered_f64_bits(candidate_x) < ordered_f64_bits(current_x))
            }
            (Some(current), None, Some(candidate), _) => candidate > current,
            _ => false,
        };
        if other_is_better {
            merged.max_ulp = other.max_ulp;
            merged.max_ulp_x = other.max_ulp_x;
            merged.max_ulp_ref = other.max_ulp_ref;
            merged.max_ulp_fpx = other.max_ulp_fpx;
        }
        merged
    }
}

#[derive(Clone, Copy)]
struct SweepStateF128 {
    max_abs: f128,
    max_abs_x: f128,
    max_abs_ref: f128,
    max_abs_fpx: f128,
    max_ulp: Option<u128>,
    max_ulp_x: Option<f128>,
    max_ulp_ref: Option<f128>,
    max_ulp_fpx: Option<f128>,
}

impl SweepStateF128 {
    #[inline]
    fn new() -> Self {
        Self {
            max_abs: 0.0_f128,
            max_abs_x: 0.0_f128,
            max_abs_ref: 0.0_f128,
            max_abs_fpx: 0.0_f128,
            max_ulp: None,
            max_ulp_x: None,
            max_ulp_ref: None,
            max_ulp_fpx: None,
        }
    }

    #[inline]
    fn observe(&mut self, x: f128, y0: f128, y1: f128) {
        if !x.is_finite() || !y0.is_finite() || !y1.is_finite() {
            return;
        }

        let abs = (y0 - y1).abs();
        if abs > self.max_abs
            || (abs == self.max_abs && ordered_f128_bits(x) < ordered_f128_bits(self.max_abs_x))
        {
            self.max_abs = abs;
            self.max_abs_x = x;
            self.max_abs_ref = y0;
            self.max_abs_fpx = y1;
        }

        if let Some(ulp) = finite_ulp_dist_f128(y0, y1) {
            let should_update = match (self.max_ulp, self.max_ulp_x) {
                (None, _) => true,
                (Some(current), Some(current_x)) => {
                    ulp > current
                        || (ulp == current && ordered_f128_bits(x) < ordered_f128_bits(current_x))
                }
                (Some(current), None) => ulp > current,
            };
            if should_update {
                self.max_ulp = Some(ulp);
                self.max_ulp_x = Some(x);
                self.max_ulp_ref = Some(y0);
                self.max_ulp_fpx = Some(y1);
            }
        }
    }

    #[inline]
    fn merge(self, other: Self) -> Self {
        let mut merged = self;
        if other.max_abs > merged.max_abs
            || (other.max_abs == merged.max_abs
                && ordered_f128_bits(other.max_abs_x) < ordered_f128_bits(merged.max_abs_x))
        {
            merged.max_abs = other.max_abs;
            merged.max_abs_x = other.max_abs_x;
            merged.max_abs_ref = other.max_abs_ref;
            merged.max_abs_fpx = other.max_abs_fpx;
        }

        let other_is_better = match (
            merged.max_ulp,
            merged.max_ulp_x,
            other.max_ulp,
            other.max_ulp_x,
        ) {
            (_, _, None, _) => false,
            (None, _, Some(_), _) => true,
            (Some(current), Some(current_x), Some(candidate), Some(candidate_x)) => {
                candidate > current
                    || (candidate == current
                        && ordered_f128_bits(candidate_x) < ordered_f128_bits(current_x))
            }
            (Some(current), None, Some(candidate), _) => candidate > current,
            _ => false,
        };
        if other_is_better {
            merged.max_ulp = other.max_ulp;
            merged.max_ulp_x = other.max_ulp_x;
            merged.max_ulp_ref = other.max_ulp_ref;
            merged.max_ulp_fpx = other.max_ulp_fpx;
        }
        merged
    }
}

/* ============================================================
 * Neighbor dump helpers
 * ============================================================ */

#[inline]
fn dump_f64_triplet(
    out: &mut String,
    label: &str,
    f_ref: fn(f64) -> f64,
    f_fpx: fn(f64) -> f64,
    x: f64,
) {
    let _ = writeln!(out, "  {label}:");
    for xi in f64_neighbors(x) {
        let y0 = f_ref(xi);
        let y1 = f_fpx(xi);
        let abs = if y0.is_finite() && y1.is_finite() {
            (y0 - y1).abs().to_string()
        } else {
            "NA".to_string()
        };
        let ulp = finite_ulp_dist_f64(y0, y1)
            .map(|value| value.to_string())
            .unwrap_or_else(|| "NA".to_string());
        let _ = writeln!(out, "    x={:+.17e}  abs={}  ulp={}", xi, abs, ulp);
    }
}

#[inline]
fn dump_f32_triplet(
    out: &mut String,
    label: &str,
    f_ref: fn(f32) -> f32,
    f_fpx: fn(f32) -> f32,
    x: f32,
) {
    let _ = writeln!(out, "  {label}:");
    for xi in f32_neighbors(x) {
        let y0 = f_ref(xi);
        let y1 = f_fpx(xi);
        let abs = if y0.is_finite() && y1.is_finite() {
            (y0 - y1).abs().to_string()
        } else {
            "NA".to_string()
        };
        let ulp = finite_ulp_dist_f32(y0, y1)
            .map(|value| value.to_string())
            .unwrap_or_else(|| "NA".to_string());
        let _ = writeln!(out, "    x={:+.9e}  abs={}  ulp={}", xi, abs, ulp);
    }
}

#[inline]
fn dump_f16_triplet(
    out: &mut String,
    label: &str,
    f_ref: fn(f16) -> f16,
    f_fpx: fn(f16) -> f16,
    x: f16,
) {
    let _ = writeln!(out, "  {label}:");
    for xi in f16_neighbors(x) {
        let y0 = f_ref(xi);
        let y1 = f_fpx(xi);
        let abs = if y0.is_finite() && y1.is_finite() {
            format!("{}", (y0 - y1).abs() as f32)
        } else {
            "NA".to_string()
        };
        let ulp = finite_ulp_dist_f16(y0, y1)
            .map(|value| value.to_string())
            .unwrap_or_else(|| "NA".to_string());
        let _ = writeln!(out, "    x={:+.5e}  abs={}  ulp={}", xi as f32, abs, ulp);
    }
}

#[inline]
fn dump_f128_triplet(
    out: &mut String,
    label: &str,
    f_ref: fn(f128) -> f128,
    f_fpx: fn(f128) -> f128,
    x: f128,
) {
    let _ = writeln!(out, "  {label}:");
    for xi in f128_neighbors(x) {
        let y0 = f_ref(xi);
        let y1 = f_fpx(xi);
        let abs = if y0.is_finite() && y1.is_finite() {
            format!("{:?}", (y0 - y1).abs())
        } else {
            "NA".to_string()
        };
        let ulp = finite_ulp_dist_f128(y0, y1)
            .map(|value| value.to_string())
            .unwrap_or_else(|| "NA".to_string());
        let _ = writeln!(out, "    x={:?}  abs={}  ulp={}", xi, abs, ulp);
    }
}

/* ============================================================
 * Exhaustive f32 bit-pattern test
 * ============================================================ */

pub fn test_f32_exhaustive(
    name: &str,
    f_ref: fn(f32) -> f32,
    f_fpx: fn(f32) -> f32,
    lo: f32,
    hi: f32,
) -> PrecisionReport {
    let state = if let Some((lo_ord, hi_ord)) = ordered_range_f32(lo, hi) {
        (lo_ord..=hi_ord)
            .into_par_iter()
            .fold(SweepStateF32::new, |mut state, ord| {
                let x = f32_from_ordered_bits(ord);
                if !x.is_finite() {
                    return state;
                }

                let y0 = f_ref(x);
                let y1 = f_fpx(x);
                state.observe(x, y0, y1);
                state
            })
            .reduce(SweepStateF32::new, SweepStateF32::merge)
    } else {
        SweepStateF32::new()
    };

    let max_abs = state.max_abs;
    let max_ulp = state.max_ulp;
    let max_abs_x = state.max_abs_x;
    let max_abs_ref = state.max_abs_ref;
    let max_abs_fpx = state.max_abs_fpx;
    let max_ulp_x = state.max_ulp_x;
    let max_ulp_ref = state.max_ulp_ref;
    let max_ulp_fpx = state.max_ulp_fpx;

    let mut report = String::new();
    let _ = writeln!(report, "{name} f32 exhaustive");
    let _ = writeln!(report, "  max abs : {:e} @ x={}", max_abs, max_abs_x);
    dump_f32_triplet(&mut report, "abs neighbors", f_ref, f_fpx, max_abs_x);

    match (max_ulp, max_ulp_x) {
        (Some(max_ulp), Some(max_ulp_x)) => {
            let _ = writeln!(report, "  max ulp : {} @ x={}", max_ulp, max_ulp_x);
            dump_f32_triplet(&mut report, "ulp neighbors", f_ref, f_fpx, max_ulp_x);
        }
        _ => {
            let _ = writeln!(report, "  max ulp : NA");
        }
    }

    PrecisionReport {
        report,
        summary: PrecisionSummary {
            name: name.to_string(),
            precision: "f32".to_string(),
            test: "exhaustive".to_string(),
            max_abs: format!("{:.9e}", max_abs),
            max_ulp: max_ulp
                .map(|value| value.to_string())
                .unwrap_or_else(|| "NA".to_string()),
            max_abs_x: format_value_f32(max_abs_x),
            max_abs_ref: format_value_f32(max_abs_ref),
            max_abs_fpx: format_value_f32(max_abs_fpx),
            max_ulp_x: format_optional_x_f32(max_ulp_x),
            max_ulp_ref: format_optional_value_f32(max_ulp_ref),
            max_ulp_fpx: format_optional_value_f32(max_ulp_fpx),
        },
    }
}

/* ============================================================
 * Stratified f32 bit-pattern test
 * ============================================================ */

pub fn test_f32_stratified_bits(
    name: &str,
    f_ref: fn(f32) -> f32,
    f_fpx: fn(f32) -> f32,
    lo: f32,
    hi: f32,
    exponents: &[u8],
    mantissas: &[u32],
) -> PrecisionReport {
    let state = exponents
        .par_iter()
        .copied()
        .map(|exp| {
            let mut state = SweepStateF32::new();
            for &mant in mantissas {
                for sign in [0_u32, 1_u32 << 31] {
                    let bits = sign | ((exp as u32) << 23) | (mant & ((1_u32 << 23) - 1));
                    let x = f32::from_bits(bits);
                    if !x.is_finite() || x < lo || x > hi {
                        continue;
                    }
                    state.observe(x, f_ref(x), f_fpx(x));
                }
            }
            state
        })
        .reduce(SweepStateF32::new, SweepStateF32::merge);

    let mut report = String::new();
    let _ = writeln!(report, "{name} f32 full-domain stratified bits");
    let _ = writeln!(
        report,
        "  max abs : {:e} @ x={}",
        state.max_abs, state.max_abs_x
    );
    match (state.max_ulp, state.max_ulp_x) {
        (Some(max_ulp), Some(max_ulp_x)) => {
            let _ = writeln!(report, "  max ulp : {} @ x={}", max_ulp, max_ulp_x);
            dump_f32_triplet(&mut report, "ulp neighbors", f_ref, f_fpx, max_ulp_x);
        }
        _ => {
            let _ = writeln!(report, "  max ulp : NA");
        }
    }

    PrecisionReport {
        report,
        summary: PrecisionSummary {
            name: name.to_string(),
            precision: "f32".to_string(),
            test: "full_domain_stratified_bits".to_string(),
            max_abs: format!("{:.9e}", state.max_abs),
            max_ulp: state
                .max_ulp
                .map(|v| v.to_string())
                .unwrap_or_else(|| "NA".to_string()),
            max_abs_x: format_value_f32(state.max_abs_x),
            max_abs_ref: format_value_f32(state.max_abs_ref),
            max_abs_fpx: format_value_f32(state.max_abs_fpx),
            max_ulp_x: format_optional_x_f32(state.max_ulp_x),
            max_ulp_ref: format_optional_value_f32(state.max_ulp_ref),
            max_ulp_fpx: format_optional_value_f32(state.max_ulp_fpx),
        },
    }
}

/* ============================================================
 * Exhaustive f16 bit-pattern test
 * ============================================================ */

pub fn test_f16_exhaustive(
    name: &str,
    f_ref: fn(f16) -> f16,
    f_fpx: fn(f16) -> f16,
    lo: f16,
    hi: f16,
) -> PrecisionReport {
    let state = if let Some((lo_ord, hi_ord)) = ordered_range_f16(lo, hi) {
        (lo_ord..=hi_ord)
            .into_par_iter()
            .fold(SweepStateF16::new, |mut state, ord| {
                let x = f16_from_ordered_bits(ord);
                if !x.is_finite() {
                    return state;
                }

                let y0 = f_ref(x);
                let y1 = f_fpx(x);
                state.observe(x, y0, y1);
                state
            })
            .reduce(SweepStateF16::new, SweepStateF16::merge)
    } else {
        SweepStateF16::new()
    };

    let max_abs = state.max_abs;
    let max_ulp = state.max_ulp;
    let max_abs_x = state.max_abs_x;
    let max_abs_ref = state.max_abs_ref;
    let max_abs_fpx = state.max_abs_fpx;
    let max_ulp_x = state.max_ulp_x;
    let max_ulp_ref = state.max_ulp_ref;
    let max_ulp_fpx = state.max_ulp_fpx;

    let mut report = String::new();
    let _ = writeln!(report, "{name} f16 exhaustive");
    let _ = writeln!(
        report,
        "  max abs : {:e} @ x={}",
        max_abs as f32, max_abs_x as f32
    );
    dump_f16_triplet(&mut report, "abs neighbors", f_ref, f_fpx, max_abs_x);

    match (max_ulp, max_ulp_x) {
        (Some(max_ulp), Some(max_ulp_x)) => {
            let _ = writeln!(report, "  max ulp : {} @ x={}", max_ulp, max_ulp_x as f32);
            dump_f16_triplet(&mut report, "ulp neighbors", f_ref, f_fpx, max_ulp_x);
        }
        _ => {
            let _ = writeln!(report, "  max ulp : NA");
        }
    }

    PrecisionReport {
        report,
        summary: PrecisionSummary {
            name: name.to_string(),
            precision: "f16".to_string(),
            test: "exhaustive".to_string(),
            max_abs: format!("{:.5e}", max_abs as f32),
            max_ulp: max_ulp
                .map(|value| value.to_string())
                .unwrap_or_else(|| "NA".to_string()),
            max_abs_x: format_value_f16(max_abs_x),
            max_abs_ref: format_value_f16(max_abs_ref),
            max_abs_fpx: format_value_f16(max_abs_fpx),
            max_ulp_x: format_optional_x_f16(max_ulp_x),
            max_ulp_ref: format_optional_value_f16(max_ulp_ref),
            max_ulp_fpx: format_optional_value_f16(max_ulp_fpx),
        },
    }
}

/* ============================================================
 * Stratified f64 bit-pattern test
 * ============================================================ */

pub fn test_f64_stratified_bits(
    name: &str,
    f_ref: fn(f64) -> f64,
    f_fpx: fn(f64) -> f64,
    lo: f64,
    hi: f64,
    exponents: &[u16],
    mantissas: &[u64],
) -> PrecisionReport {
    let state = exponents
        .par_iter()
        .copied()
        .map(|exp| {
            let mut state = SweepStateF64::new();
            for &mant in mantissas {
                for sign in [0_u64, 1_u64 << 63] {
                    let bits = sign | ((exp as u64) << 52) | (mant & ((1u64 << 52) - 1));
                    let x = f64::from_bits(bits);

                    if !x.is_finite() || x < lo || x > hi {
                        continue;
                    }

                    let y0 = f_ref(x);
                    let y1 = f_fpx(x);
                    state.observe(x, y0, y1);
                }
            }
            state
        })
        .reduce(SweepStateF64::new, SweepStateF64::merge);

    let max_abs = state.max_abs;
    let max_ulp = state.max_ulp;
    let max_abs_x = state.max_abs_x;
    let max_abs_ref = state.max_abs_ref;
    let max_abs_fpx = state.max_abs_fpx;
    let max_ulp_x = state.max_ulp_x;
    let max_ulp_ref = state.max_ulp_ref;
    let max_ulp_fpx = state.max_ulp_fpx;

    let mut report = String::new();
    let _ = writeln!(report, "{name} f64 stratified bits");
    let _ = writeln!(report, "  max abs : {:e} @ x={}", max_abs, max_abs_x);
    dump_f64_triplet(&mut report, "abs neighbors", f_ref, f_fpx, max_abs_x);

    match (max_ulp, max_ulp_x) {
        (Some(max_ulp), Some(max_ulp_x)) => {
            let _ = writeln!(report, "  max ulp : {} @ x={}", max_ulp, max_ulp_x);
            dump_f64_triplet(&mut report, "ulp neighbors", f_ref, f_fpx, max_ulp_x);
        }
        _ => {
            let _ = writeln!(report, "  max ulp : NA");
        }
    }

    PrecisionReport {
        report,
        summary: PrecisionSummary {
            name: name.to_string(),
            precision: "f64".to_string(),
            test: "full_domain_stratified_bits".to_string(),
            max_abs: format!("{:.17e}", max_abs),
            max_ulp: max_ulp
                .map(|value| value.to_string())
                .unwrap_or_else(|| "NA".to_string()),
            max_abs_x: format_value_f64(max_abs_x),
            max_abs_ref: format_value_f64(max_abs_ref),
            max_abs_fpx: format_value_f64(max_abs_fpx),
            max_ulp_x: format_optional_x_f64(max_ulp_x),
            max_ulp_ref: format_optional_value_f64(max_ulp_ref),
            max_ulp_fpx: format_optional_value_f64(max_ulp_fpx),
        },
    }
}

/* ============================================================
 * Stratified f128 bit-pattern test
 * ============================================================ */

pub fn test_f128_stratified_bits(
    name: &str,
    f_ref: fn(f128) -> f128,
    f_fpx: fn(f128) -> f128,
    lo: f128,
    hi: f128,
    exponents: &[u16],
    mantissas: &[u128],
) -> PrecisionReport {
    let mant_mask: u128 = (1u128 << 112) - 1;

    let state = exponents
        .par_iter()
        .copied()
        .map(|exp| {
            let mut state = SweepStateF128::new();
            for &mant in mantissas {
                for sign in [0_u128, 1_u128 << 127] {
                    let bits = sign | ((exp as u128) << 112) | (mant & mant_mask);
                    let x = f128::from_bits(bits);

                    if !x.is_finite() || x < lo || x > hi {
                        continue;
                    }

                    let y0 = f_ref(x);
                    let y1 = f_fpx(x);
                    state.observe(x, y0, y1);
                }
            }
            state
        })
        .reduce(SweepStateF128::new, SweepStateF128::merge);

    let max_abs = state.max_abs;
    let max_ulp = state.max_ulp;
    let max_abs_x = state.max_abs_x;
    let max_abs_ref = state.max_abs_ref;
    let max_abs_fpx = state.max_abs_fpx;
    let max_ulp_x = state.max_ulp_x;
    let max_ulp_ref = state.max_ulp_ref;
    let max_ulp_fpx = state.max_ulp_fpx;

    let mut report = String::new();
    let _ = writeln!(report, "{name} f128 stratified bits");
    let _ = writeln!(report, "  max abs : {:?} @ x={:?}", max_abs, max_abs_x);
    dump_f128_triplet(&mut report, "abs neighbors", f_ref, f_fpx, max_abs_x);

    match (max_ulp, max_ulp_x) {
        (Some(max_ulp), Some(max_ulp_x)) => {
            let _ = writeln!(report, "  max ulp : {} @ x={:?}", max_ulp, max_ulp_x);
            dump_f128_triplet(&mut report, "ulp neighbors", f_ref, f_fpx, max_ulp_x);
        }
        _ => {
            let _ = writeln!(report, "  max ulp : NA");
        }
    }

    PrecisionReport {
        report,
        summary: PrecisionSummary {
            name: name.to_string(),
            precision: "f128".to_string(),
            test: "full_domain_stratified_bits".to_string(),
            max_abs: format!("{:?}", max_abs),
            max_ulp: max_ulp
                .map(|value| value.to_string())
                .unwrap_or_else(|| "NA".to_string()),
            max_abs_x: format_value_f128(max_abs_x),
            max_abs_ref: format_value_f128(max_abs_ref),
            max_abs_fpx: format_value_f128(max_abs_fpx),
            max_ulp_x: format_optional_x_f128(max_ulp_x),
            max_ulp_ref: format_optional_value_f128(max_ulp_ref),
            max_ulp_fpx: format_optional_value_f128(max_ulp_fpx),
        },
    }
}

/* ============================================================
 * PI-fraction refinement, using ordered-bit neighborhoods.
 * This avoids raw to_bits().wrapping_sub() around +0.0.
 * ============================================================ */

#[inline]
fn add_ordered_ulp_f64(center: f64, k: i64) -> Option<f64> {
    if !center.is_finite() {
        return None;
    }

    let ord = ordered_f64_bits(center) as i128;
    let next = ord + k as i128;

    if !(0..=(u64::MAX as i128)).contains(&next) {
        return None;
    }

    let x = f64_from_ordered_bits(next as u64);
    x.is_finite().then_some(x)
}

pub fn test_f64_pi_refinement(
    name: &str,
    f_ref: fn(f64) -> f64,
    f_fpx: fn(f64) -> f64,
    ulp_radius: i64,
) -> PrecisionReport {
    let centers = [0.0, PI / 6.0, PI / 4.0, PI / 3.0, PI / 2.0];

    let state = centers
        .par_iter()
        .copied()
        .map(|c| {
            let mut state = SweepStateF64::new();
            for k in -ulp_radius..=ulp_radius {
                let Some(x) = add_ordered_ulp_f64(c, k) else {
                    continue;
                };

                let y0 = f_ref(x);
                let y1 = f_fpx(x);
                state.observe(x, y0, y1);
            }
            state
        })
        .reduce(SweepStateF64::new, SweepStateF64::merge);

    let max_abs = state.max_abs;
    let max_ulp = state.max_ulp;
    let max_abs_x = state.max_abs_x;
    let max_abs_ref = state.max_abs_ref;
    let max_abs_fpx = state.max_abs_fpx;
    let max_ulp_x = state.max_ulp_x;
    let max_ulp_ref = state.max_ulp_ref;
    let max_ulp_fpx = state.max_ulp_fpx;

    let mut report = String::new();
    let _ = writeln!(report, "{name} f64 PI-refinement");
    let _ = writeln!(report, "  max abs : {:e} @ x={}", max_abs, max_abs_x);
    dump_f64_triplet(&mut report, "abs neighbors", f_ref, f_fpx, max_abs_x);

    match (max_ulp, max_ulp_x) {
        (Some(max_ulp), Some(max_ulp_x)) => {
            let _ = writeln!(report, "  max ulp : {} @ x={}", max_ulp, max_ulp_x);
            dump_f64_triplet(&mut report, "ulp neighbors", f_ref, f_fpx, max_ulp_x);
        }
        _ => {
            let _ = writeln!(report, "  max ulp : NA");
        }
    }

    PrecisionReport {
        report,
        summary: PrecisionSummary {
            name: name.to_string(),
            precision: "f64".to_string(),
            test: format!("pi_refinement_{}", ulp_radius),
            max_abs: format!("{:.17e}", max_abs),
            max_ulp: max_ulp
                .map(|value| value.to_string())
                .unwrap_or_else(|| "NA".to_string()),
            max_abs_x: format_value_f64(max_abs_x),
            max_abs_ref: format_value_f64(max_abs_ref),
            max_abs_fpx: format_value_f64(max_abs_fpx),
            max_ulp_x: format_optional_x_f64(max_ulp_x),
            max_ulp_ref: format_optional_value_f64(max_ulp_ref),
            max_ulp_fpx: format_optional_value_f64(max_ulp_fpx),
        },
    }
}

/* ============================================================
 * Endpoint + boundary testing
 * ============================================================ */

pub fn test_f64_boundaries(
    name: &str,
    f_ref: fn(f64) -> f64,
    f_fpx: fn(f64) -> f64,
    boundaries: &[f64],
) -> PrecisionReport {
    let state = boundaries
        .par_iter()
        .copied()
        .map(|b| {
            let mut state = SweepStateF64::new();
            for x in f64_neighbors(b) {
                let y0 = f_ref(x);
                let y1 = f_fpx(x);
                state.observe(x, y0, y1);
            }
            state
        })
        .reduce(SweepStateF64::new, SweepStateF64::merge);

    let max_abs = state.max_abs;
    let max_ulp = state.max_ulp;
    let max_abs_x = state.max_abs_x;
    let max_abs_ref = state.max_abs_ref;
    let max_abs_fpx = state.max_abs_fpx;
    let max_ulp_x = state.max_ulp_x;
    let max_ulp_ref = state.max_ulp_ref;
    let max_ulp_fpx = state.max_ulp_fpx;

    let mut report = String::new();
    let _ = writeln!(report, "{name} f64 boundaries");
    let _ = writeln!(report, "  max abs : {:e} @ x={}", max_abs, max_abs_x);
    dump_f64_triplet(&mut report, "abs neighbors", f_ref, f_fpx, max_abs_x);

    match (max_ulp, max_ulp_x) {
        (Some(max_ulp), Some(max_ulp_x)) => {
            let _ = writeln!(report, "  max ulp : {} @ x={}", max_ulp, max_ulp_x);
            dump_f64_triplet(&mut report, "ulp neighbors", f_ref, f_fpx, max_ulp_x);
        }
        _ => {
            let _ = writeln!(report, "  max ulp : NA");
        }
    }

    PrecisionReport {
        report,
        summary: PrecisionSummary {
            name: name.to_string(),
            precision: "f64".to_string(),
            test: "boundaries".to_string(),
            max_abs: format!("{:.17e}", max_abs),
            max_ulp: max_ulp
                .map(|value| value.to_string())
                .unwrap_or_else(|| "NA".to_string()),
            max_abs_x: format_value_f64(max_abs_x),
            max_abs_ref: format_value_f64(max_abs_ref),
            max_abs_fpx: format_value_f64(max_abs_fpx),
            max_ulp_x: format_optional_x_f64(max_ulp_x),
            max_ulp_ref: format_optional_value_f64(max_ulp_ref),
            max_ulp_fpx: format_optional_value_f64(max_ulp_fpx),
        },
    }
}
