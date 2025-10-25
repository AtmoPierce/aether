# Aether Workspace

> **Aether** — an open, strongly-typed scientific computing framework for aerospace and physics simulation in Rust.

This workspace contains the family of crates that make up **Aether**, a modular library suite for high-fidelity dynamics, control, and visualization.  
Aether was designed to enable rigorous, reference-frame-aware computation for spacecraft, aircraft, and robotic systems — without sacrificing performance or clarity. Please feel free to contribute your own models, mathematics, and thoughts on this work.

---

## Scientific Intent

Aether was originally developed as part of my internal research ecosystem for **Guidance, Navigation, and Control (GNC)**, **physics-based simulation**, and **model-predictive optimization**.  
The intent behind publishing these crates is to **advance open, verifiable research in scientific computing** and **promote reproducibility in aerospace and control theory**.

Our goals for the open library family are:

1. **Transparency and reproducibility** — scientific methods should be inspectable, not black boxes.
2. **Type-safe physics** — coordinate frames, units, and state representations should be encoded in the type system.
3. **Performance through purity** — numerics should compile down to optimal code paths, suitable for HPC and embedded environments alike.
4. **Interdisciplinary reusability** — the same primitives should serve orbital dynamics, machine learning, and structural analysis.
5. **Open acceleration of science** — researchers, students, and engineers should be able to build upon a common foundation of reliable math and simulation code.

Aether's open core is therefore made public to contribute to a verifiable, type-safe, and hardware-accelerated ecosystem for computational physics.

---

## Workspace Structure

This workspace collects the following crates:

| Crate | Description |
|-------|--------------|
| **aether** | Umbrella crate aggregating core modules for convenience. |
| **aether_core** | Strongly-typed math foundation — matrices, vectors, quaternions, and reference frame abstractions. |
| **aether_models** | Physical and dynamical models (rigid bodies, atmosphere, gravitation, etc.). |
| **aether_fluids** | Fluid dynamics and continuum-mechanics primitives. |
| **aether_shapes** | Geometric primitives for collision, inertia, and volumetric modeling. |
| **aether_graphics** | Rendering and visualization utilities for simulation and analysis. |
| **aether_ml** | Lightweight scientific machine learning and regression utilities integrated with Aether math types. |
| **aether_opt** | Optimization and control algorithms (gradient descent, Riccati solvers, MPC scaffolding). |
| **aether_rand** | Deterministic RNGs and sampling utilities for simulations. |
| **aether_stats** | Statistical analysis, regression, and signal-processing tools. |
| **aether_benchmark** | Performance tests and HPC kernels for benchmarking Aether numerics. |

---

## Philosophy

Aether is built on three core principles:

1. **Strong Types for Strong Science**  
   Physical correctness is encoded at compile-time using the Rust type system.  
   Reference frames, units, and coordinate systems are not strings — they are types.

2. **Numerical Transparency**  
   Every algorithm is written with clear linear-algebraic intent, leveraging static matrices, explicit arithmetic, and traceable computation steps.

3. **Unified Continuum**  
   Whether simulating a spacecraft’s attitude, a fluid field, or a statistical regression, all Aether crates share the same mathematical substrate.

---

## Build and Usage

To build the entire workspace:

```bash
git clone https://github.com/atmopierce/aether.git
cd aether
cargo build --release
```

To include a singular package via git
```toml
[dependencies]
aether_core = { git = "https://github.com/nuntius-aerospace/aether.git", package = "aether_core" }
```

## Citation
Michael Angeles. Aether: A Strongly-Typed Scientific Computing Framework for Simulation in Rust. 2025.
https://github.com/atmopierce/aether

## Contributing
While Aether is primarily developed in support of my research, community contributions are welcome — especially those improving mathematical clarity, documentation, or academic integration (e.g., external bindings, notebooks, examples).

Please open an issue or pull request if you wish to contribute.