/* ---------------- formatting ---------------- */

pub fn fmt_i64(v: Option<f64>) -> String {
    match v {
        Some(x) => format!("{}", x as i64),
        None => "NA".into(),
    }
}

pub fn fmt_f64(v: Option<f64>) -> String {
    match v {
        Some(x) => format!("{x:.6}"),
        None => "NA".into(),
    }
}

pub fn opt_f64(v: Option<f64>) -> String {
    match v {
        Some(x) => format!("{x}"),
        None => String::new(),
    }
}
