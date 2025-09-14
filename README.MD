
**Aether** is a Rust-native, no_std-friendly linear algebra library and interactive “lab” for vectors and matrices.  
It provides statically-typed, stack-allocated types, operator overloading, and an optional interactive REPL (Read-Eval-Print Loop) for hands-on computations and experimentation.

![Aether](docs/images/aether.png)

---

## Features

- **Statically-sized, stack-allocated vectors and matrices:**  
  - Generic types: `Vector<T, N>`, `Matrix<T, R, C>`
- **No heap allocation:** Efficient and suitable for embedded or `no_std` environments.
- **Operator overloading:** Supports `+`, `-`, `*`, `/` (elementwise and matrix-matrix)
- **Space-separated, bracketed syntax:**  
  - Vectors: `[1 2 3]`
  - Matrices: `[1 2; 3 4]`
- **Interactive REPL:** Try out calculations instantly from the command line.
- **Clear, readable output:**  
  `[2 4 6]`, `[6 8; 10 12]`
- **Flexible usage:** Use as a library or as an interactive calculator.

> **Note:**  
> The interactive REPL supports a subset of the available operations and types for demonstration and convenience.  
> **All operations are available and fully supported when using Aether as a Rust library.**

---

## Library Example

```rust
use aether::{Vector, Matrix};

fn main() {
    let v1 = Vector::<f64, 3> { data: [1.0, 2.0, 3.0] };
    let v2 = Vector::<f64, 3> { data: [4.0, 5.0, 6.0] };
    let sum = v1 + v2;
    println!("{:?}", sum); // Output: Vector { data: [5.0, 7.0, 9.0] }

    let m1 = Matrix::<f64, 2, 2> { data: [[1.0, 2.0], [3.0, 4.0]] };
    let m2 = Matrix::<f64, 2, 2> { data: [[5.0, 6.0], [7.0, 8.0]] };
    let prod = m1 * m2;
    println!("{:?}", prod); // Output: Matrix { data: [[19.0, 22.0], [43.0, 50.0]] }
}
```
---
## Example REPL Session
```
>> [1 2 3] + [4 5 6]
[5 7 9]
>> [1 2; 3 4] * [5 6; 7 8]
[19 22; 43 50]
>> [1 2 3] * 2
[2 4 6]
```

**Planned and Proposed Features:**
- **Parallel and SIMD-enabled operations**  
  [ ] Leverage Rust’s parallel iterators and/or explicit SIMD for large statically-sized types (when supported by hardware)
- **Support for larger matrices and vectors**  
  [ ] Ergonomic helpers for generating and handling NxM matrices (with compile-time or run-time dimensions)
- **Advanced algebra:**  
  [ ] Determinants, matrix inversion, eigenvalues/eigenvectors for supported types
  [ ] QR, LU, and SVD decomposition
- **Flexible output and formatting:**  
  [ ] Optional pretty-print for CLI/repl and library display
  [ ] Support for exporting to CSV/JSON/text
- **Extended CLI syntax:**  
  [ ] Chained expressions, parentheses for grouping, transposition (`'`), slicing/indexing
  [ ] User-defined variables in the REPL session
  [ ] Saving of variables as a memento
- **Interoperability:**  
  [ ] Bindings or bridges for FFI to C, CUDA kernels
- **Parallel CLI:**  
  [ ] parallelize bulk calculations or REPL commands (e.g., batch evaluation)
- **Extensive tests and documentation:**  
  [ ] Property-based tests, benchmarks, and `no_std`/embedded examples