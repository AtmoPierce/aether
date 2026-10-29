# fpx - Floating-Point eXtensions

`fpx` is a `no_std` / core-only Rust library focused on constructing **minimax polynomial approximations** of standard mathematical functions using **Chebyshev** and **Remez** techniques.
The crate is designed around transparent numerical methods: explicit interval selection, boundary-value handling, and static, allocation-free coefficient storage.

Although no functions are implemented yet, the project provides the mathematical and structural foundation for building accurate, deterministic, and portable math kernels entirely atop **core Rust**.

The library is made public for complete auditability of the process, and finetuning for platform specificity.

Serializable/Deserializable coefficient tables and boundaries shall also be implemented to enable extension to other languages.

## Features (in progress)

- `no_std` / core compatibility  
- IEEE-754 bit-level utilities  
- Deterministic polynomial evaluation (Horner, Clenshaw)  
- Minimax coefficient generation (Chebyshev / Remez)  
- Static `const` coefficient storage (no allocation)  
- Clear numerical boundaries and interval reduction strategy  

## Planned Functionality

### Trigonometric
- [ ] sin(x)
- [ ] cos(x)
- [ ] tan(x)
- [ ] asin(x)
- [ ] acos(x)
- [ ] atan(x)
- [ ] atan2(y, x)
- [ ] ...

### Exponential & Logarithmic
- [ ] exp(x)
- [ ] exp2(x)
- [ ] exp10(x)
- [ ] ln(x)
- [ ] log2(x)
- [ ] log10(x)
- [ ] ...

### Hyperbolic
- [ ] sinh(x)
- [ ] cosh(x)
- [ ] tanh(x)
- [ ] ...

### Reciprocal & Roots
- [ ] 1/x
- [ ] sqrt(x)
- [ ] cbrt(x)
- [ ] ...

### Error & Special Functions
- [ ] erf(x)
- [ ] erfc(x)
- [ ] gamma(x)
- [ ] lgamma(x)
- [ ] ...

## Numerical Foundations

### Primary References
- William J. Cody & William Waite, *Software Manual for the Elementary Functions*, Prentence-Hall, 1980.
- J. F. Hart et al., *Computer Approximations*, Wiley, 1968.
- E. Y. Remez, *Sur un procédé convergent d'approximations successives pour déterminer les polynômes d'approximation*, 1934.

### Additional References
- Hart, Cheney, Lawson, Maehly, Selfridge, Wheeler, *Table of Chebyshev Approximations*, 1968.
- Cody, "Chebyshev Approximations for the Error Function" (NIST/NASA)
- Demanet & Ying, *On Chebyshev Interpolation of Analytic Functions*, 2010.

## Status
The current version defines the mathematical scaffold (Chebyshev basis, Remez infrastructure, coefficient layout), with function implementations planned.

## Cortex-M33 inline assembly module

`aether-fpx` includes a Cortex-M33 assembly probe under `crates/aether-fpx/src/bin`.

Feature flags:

- `m33-inline-asm`: enables raw instruction emission helpers (`.inst.w`)
- `m33-inline-asm-mnemonics`: enables mnemonic wrappers (currently a `CX1` starter wrapper)

Example build:

```bash
cargo check -p aether-fpx --features m33-inline-asm
```

On Cortex-M targets, use `m33::emit_inst32::<OPCODE>()` for raw FPX opcodes and grow typed wrappers from there.

## Optional x87 research baseline

The `aether-fpx-research` crate includes an `x87` feature for its hardware-specific
binary64 benchmark wrappers. It is enabled by default on supported x86 and
x86-64 targets, but the build script automatically skips it on other
architectures, including Raspberry Pi targets.

To disable x87 explicitly, build the research crate without its default
features:

```bash
cargo build -p aether-fpx-research --no-default-features
```

Optional features can be added independently; for example, the visualization
tools without x87 are built with:

```bash
cargo build -p aether-fpx-research --no-default-features --features viz
```
