use matrslab::{Matrix, Vector};
use regex::Regex;

fn parse_vector_dynamic(s: &str) -> Option<Vec<f64>> {
    let inside = s.trim().trim_start_matches('[').trim_end_matches(']');
    let nums: Vec<f64> = inside
        .split(|c| c == ' ' || c == ',')
        .filter(|x| !x.trim().is_empty())
        .map(|x| x.trim().parse::<f64>())
        .collect::<Result<_, _>>()
        .ok()?;
    if !nums.is_empty() {
        Some(nums)
    } else {
        None
    }
}

fn format_vector(v: &[f64]) -> String {
    let s = v
        .iter()
        .map(|x| format!("{}", x))
        .collect::<Vec<_>>()
        .join(" ");
    format!("[{}]", s)
}

fn format_matrix(rows: usize, cols: usize, data: &[f64]) -> String {
    let mut out = String::from("[");
    for r in 0..rows {
        if r > 0 {
            out.push_str("; ");
        }
        let row = (0..cols)
            .map(|c| format!("{}", data[r * cols + c]))
            .collect::<Vec<_>>()
            .join(" ");
        out.push_str(&row);
    }
    out.push(']');
    out
}

/// Parse a 2x2 or 3x3 matrix from [1 2; 3 4] or [1,2;3,4] syntax
fn parse_matrix(s: &str) -> Option<(usize, usize, Vec<f64>)> {
    let inside = s.trim().trim_start_matches('[').trim_end_matches(']');
    let rows: Vec<&str> = inside.split(';').collect();
    let matrix: Vec<Vec<f64>> = rows
        .iter()
        .map(|row| {
            row.split(|c| c == ' ' || c == ',')
                .filter(|x| !x.trim().is_empty())
                .map(|x| x.trim().parse::<f64>())
                .collect::<Result<Vec<f64>, _>>()
        })
        .collect::<Result<Vec<_>, _>>()
        .ok()?;
    let nrows = matrix.len();
    let ncols = if nrows > 0 { matrix[0].len() } else { 0 };
    let flat: Vec<f64> = matrix.into_iter().flat_map(|row| row.into_iter()).collect();
    Some((nrows, ncols, flat))
}

fn parse_scalar(s: &str) -> Option<f64> {
    s.trim().parse().ok()
}

fn vector_elementwise_op<F>(a: &[f64], b: &[f64], f: F) -> Result<Vec<f64>, String>
where
    F: Fn(f64, f64) -> f64,
{
    if a.len() != b.len() {
        Err(format!(
            "Dimension mismatch: left has {} elements, right has {}.",
            a.len(),
            b.len()
        ))
    } else {
        Ok(a.iter().zip(b.iter()).map(|(&x, &y)| f(x, y)).collect())
    }
}

pub fn parse_and_eval(expr: &str) -> Result<String, String> {
    let expr = expr.trim();

    // -------- Matrix constructor support ---------
    let eye_re = Regex::new(r"^(eye|identity)\((\d+)\)$").unwrap();
    if let Some(caps) = eye_re.captures(expr) {
        let n: usize = caps[2].parse().unwrap();
        return match n {
            2 => Ok(format_matrix(
                2,
                2,
                &Matrix::<f64, 2, 2>::identity().data.concat(),
            )),
            3 => Ok(format_matrix(
                3,
                3,
                &Matrix::<f64, 3, 3>::identity().data.concat(),
            )),
            4 => Ok(format_matrix(
                4,
                4,
                &Matrix::<f64, 4, 4>::identity().data.concat(),
            )),
            _ => Err("Only n=2,3,4 supported".to_string()),
        };
    }

    let diag_re = Regex::new(r"^diag\(\s*(\[.*\])\s*\)$").unwrap();
    if let Some(caps) = diag_re.captures(expr) {
        let v = parse_vector_dynamic(&caps[1]).ok_or("Failed to parse diag vector")?;
        return match v.len() {
            2 => Ok(format_matrix(
                2,
                2,
                &Matrix::<f64, 2, 2>::diag(&[v[0], v[1]]).data.concat(),
            )),
            3 => Ok(format_matrix(
                3,
                3,
                &Matrix::<f64, 3, 3>::diag(&[v[0], v[1], v[2]]).data.concat(),
            )),
            4 => Ok(format_matrix(
                4,
                4,
                &Matrix::<f64, 4, 4>::diag(&[v[0], v[1], v[2], v[3]])
                    .data
                    .concat(),
            )),
            _ => Err("Only diag with 2, 3, or 4 elements supported.".to_string()),
        };
    }

    let zeros_re = Regex::new(r"^zeros\(\s*(\d+)\s*,\s*(\d+)\s*\)$").unwrap();
    if let Some(caps) = zeros_re.captures(expr) {
        let m: usize = caps[1].parse().unwrap();
        let n: usize = caps[2].parse().unwrap();
        return match (m, n) {
            (2, 2) => Ok(format_matrix(
                2,
                2,
                &Matrix::<f64, 2, 2>::zeros().data.concat(),
            )),
            (3, 3) => Ok(format_matrix(
                3,
                3,
                &Matrix::<f64, 3, 3>::zeros().data.concat(),
            )),
            (4, 4) => Ok(format_matrix(
                4,
                4,
                &Matrix::<f64, 4, 4>::zeros().data.concat(),
            )),
            _ => Err("Only zeros of 2x2, 3x3, or 4x4 supported.".to_string()),
        };
    }

    let ones_re = Regex::new(r"^ones\(\s*(\d+)\s*,\s*(\d+)\s*\)$").unwrap();
    if let Some(caps) = ones_re.captures(expr) {
        let m: usize = caps[1].parse().unwrap();
        let n: usize = caps[2].parse().unwrap();
        return match (m, n) {
            (2, 2) => Ok(format_matrix(
                2,
                2,
                &Matrix::<f64, 2, 2>::ones().data.concat(),
            )),
            (3, 3) => Ok(format_matrix(
                3,
                3,
                &Matrix::<f64, 3, 3>::ones().data.concat(),
            )),
            (4, 4) => Ok(format_matrix(
                4,
                4,
                &Matrix::<f64, 4, 4>::ones().data.concat(),
            )),
            _ => Err("Only ones of 2x2, 3x3, or 4x4 supported.".to_string()),
        };
    }
    // -------- End constructor support ---------

    // Operator handling (unchanged from your version)
    let ops = ['+', '-', '*', '/'];
    let mut op_idx = None;
    let mut op_char = ' ';
    let mut depth = 0;
    for (i, c) in expr.chars().enumerate() {
        match c {
            '[' => depth += 1,
            ']' => depth -= 1,
            _ if ops.contains(&c) && depth == 0 => {
                op_idx = Some(i);
                op_char = c;
                break;
            }
            _ => (),
        }
    }
    let idx = op_idx.ok_or("No operator found. Use +, -, *, or /.")?;
    let left = &expr[..idx];
    let right = &expr[idx + 1..];

    // Dynamic N-vectors
    let left_vec = parse_vector_dynamic(left);
    let right_vec = parse_vector_dynamic(right);
    let left_scalar = parse_scalar(left);
    let right_scalar = parse_scalar(right);
    let left_mat = parse_matrix(left);
    let right_mat = parse_matrix(right);

    match (
        left_vec.as_ref(),
        right_vec.as_ref(),
        left_scalar,
        right_scalar,
        &left_mat,
        &right_mat,
        op_char,
    ) {
        // [n] + [n] (N-dim vector addition)
        (Some(a), Some(b), _, _, _, _, '+') => {
            vector_elementwise_op(a, b, |x, y| x + y).map(|v| format_vector(&v))
        }
        (Some(a), Some(b), _, _, _, _, '-') => {
            vector_elementwise_op(a, b, |x, y| x - y).map(|v| format_vector(&v))
        }
        (Some(a), _, _, Some(s), _, _, '*') => {
            Ok(format_vector(&a.iter().map(|x| x * s).collect::<Vec<_>>()))
        }
        (Some(a), _, _, Some(s), _, _, '/') => {
            if s == 0.0 {
                Err("Division by zero!".to_string())
            } else {
                Ok(format_vector(&a.iter().map(|x| x / s).collect::<Vec<_>>()))
            }
        }
        (_, Some(b), Some(s), _, _, _, '*') => {
            Ok(format_vector(&b.iter().map(|x| x * s).collect::<Vec<_>>()))
        }
        // Matrix * scalar (e.g., [1 2; 3 4] * 2)
        (_, _, _, Some(s), Some((r, c, m)), _, '*')
            if (*r == 2 && *c == 2) || (*r == 3 && *c == 3) =>
        {
            match (*r, *c) {
                (2, 2) => {
                    let a = Matrix::<f64, 2, 2> {
                        data: [[m[0], m[1]], [m[2], m[3]]],
                    };
                    let ab = a * s;
                    let flat = [ab.data[0][0], ab.data[0][1], ab.data[1][0], ab.data[1][1]];
                    Ok(format_matrix(2, 2, &flat))
                }
                (3, 3) => {
                    let a = Matrix::<f64, 3, 3> {
                        data: [[m[0], m[1], m[2]], [m[3], m[4], m[5]], [m[6], m[7], m[8]]],
                    };
                    let ab = a * s;
                    let flat = [
                        ab.data[0][0],
                        ab.data[0][1],
                        ab.data[0][2],
                        ab.data[1][0],
                        ab.data[1][1],
                        ab.data[1][2],
                        ab.data[2][0],
                        ab.data[2][1],
                        ab.data[2][2],
                    ];
                    Ok(format_matrix(3, 3, &flat))
                }
                _ => Err("Matrix size not supported".to_string()),
            }
        }
        // Scalar * Matrix (e.g., 2 * [1 2; 3 4])
        (_, _, Some(s), _, _, Some((r, c, m)), '*')
            if (*r == 2 && *c == 2) || (*r == 3 && *c == 3) =>
        {
            match (*r, *c) {
                (2, 2) => {
                    let a = Matrix::<f64, 2, 2> {
                        data: [[m[0], m[1]], [m[2], m[3]]],
                    };
                    let ab = a * s;
                    let flat = [ab.data[0][0], ab.data[0][1], ab.data[1][0], ab.data[1][1]];
                    Ok(format_matrix(2, 2, &flat))
                }
                (3, 3) => {
                    let a = Matrix::<f64, 3, 3> {
                        data: [[m[0], m[1], m[2]], [m[3], m[4], m[5]], [m[6], m[7], m[8]]],
                    };
                    let ab = a * s;
                    let flat = [
                        ab.data[0][0],
                        ab.data[0][1],
                        ab.data[0][2],
                        ab.data[1][0],
                        ab.data[1][1],
                        ab.data[1][2],
                        ab.data[2][0],
                        ab.data[2][1],
                        ab.data[2][2],
                    ];
                    Ok(format_matrix(3, 3, &flat))
                }
                _ => Err("Matrix size not supported".to_string()),
            }
        }
        // Matrix + Matrix ([1 2;3 4] + [5 6;7 8])
        (_, _, _, _, Some((r1, c1, m1)), Some((r2, c2, m2)), '+')
            if r1 == r2 && c1 == c2 && ((*r1 == 2 && *c1 == 2) || (*r1 == 3 && *c1 == 3)) =>
        {
            match (*r1, *c1) {
                (2, 2) => {
                    let a = Matrix::<f64, 2, 2> {
                        data: [[m1[0], m1[1]], [m1[2], m1[3]]],
                    };
                    let b = Matrix::<f64, 2, 2> {
                        data: [[m2[0], m2[1]], [m2[2], m2[3]]],
                    };
                    let ab = a + b;
                    let flat = [ab.data[0][0], ab.data[0][1], ab.data[1][0], ab.data[1][1]];
                    Ok(format_matrix(2, 2, &flat))
                }
                (3, 3) => {
                    let a = Matrix::<f64, 3, 3> {
                        data: [
                            [m1[0], m1[1], m1[2]],
                            [m1[3], m1[4], m1[5]],
                            [m1[6], m1[7], m1[8]],
                        ],
                    };
                    let b = Matrix::<f64, 3, 3> {
                        data: [
                            [m2[0], m2[1], m2[2]],
                            [m2[3], m2[4], m2[5]],
                            [m2[6], m2[7], m2[8]],
                        ],
                    };
                    let ab = a + b;
                    let flat = [
                        ab.data[0][0],
                        ab.data[0][1],
                        ab.data[0][2],
                        ab.data[1][0],
                        ab.data[1][1],
                        ab.data[1][2],
                        ab.data[2][0],
                        ab.data[2][1],
                        ab.data[2][2],
                    ];
                    Ok(format_matrix(3, 3, &flat))
                }
                _ => Err("Matrix size not supported".to_string()),
            }
        }
        // Matrix * Matrix ([a] * [b])
        (_, _, _, _, Some((r1, c1, m1)), Some((r2, c2, m2)), '*')
            if *c1 == *r2
                && ((*r1 == 2 && *c1 == 2 && *r2 == 2 && *c2 == 2)
                    || (*r1 == 3 && *c1 == 3 && *r2 == 3 && *c2 == 3)) =>
        {
            match (*r1, *c1, *r2, *c2) {
                (2, 2, 2, 2) => {
                    let a = Matrix::<f64, 2, 2> {
                        data: [[m1[0], m1[1]], [m1[2], m1[3]]],
                    };
                    let b = Matrix::<f64, 2, 2> {
                        data: [[m2[0], m2[1]], [m2[2], m2[3]]],
                    };
                    let ab = a * b;
                    let flat = [ab.data[0][0], ab.data[0][1], ab.data[1][0], ab.data[1][1]];
                    Ok(format_matrix(2, 2, &flat))
                }
                (3, 3, 3, 3) => {
                    let a = Matrix::<f64, 3, 3> {
                        data: [
                            [m1[0], m1[1], m1[2]],
                            [m1[3], m1[4], m1[5]],
                            [m1[6], m1[7], m1[8]],
                        ],
                    };
                    let b = Matrix::<f64, 3, 3> {
                        data: [
                            [m2[0], m2[1], m2[2]],
                            [m2[3], m2[4], m2[5]],
                            [m2[6], m2[7], m2[8]],
                        ],
                    };
                    let ab = a * b;
                    let flat = [
                        ab.data[0][0],
                        ab.data[0][1],
                        ab.data[0][2],
                        ab.data[1][0],
                        ab.data[1][1],
                        ab.data[1][2],
                        ab.data[2][0],
                        ab.data[2][1],
                        ab.data[2][2],
                    ];
                    Ok(format_matrix(3, 3, &flat))
                }
                _ => Err("Matrix size not supported".to_string()),
            }
        }
        _ => Err("Couldn't parse or match any supported operation or dimension.".to_string()),
    }
}
