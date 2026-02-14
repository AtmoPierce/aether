#!/usr/bin/env python3
import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CRITERION_DIR = ROOT / "target" / "criterion"
PATTERN = re.compile(r"^(matmul|matvec|dot)_(\d+)_([a-z0-9_]+)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare Criterion benchmark means against nalgebra and optionally generate speedup plots."
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Disable PNG plot generation.",
    )
    parser.add_argument(
        "--plot-dir",
        type=Path,
        default=CRITERION_DIR / "plots",
        help="Directory where plot PNG files are written.",
    )
    return parser.parse_args()


def load_means() -> dict[tuple[str, int, str], float]:
    means: dict[tuple[str, int, str], float] = {}

    if not CRITERION_DIR.exists():
        return means

    for estimates_file in CRITERION_DIR.glob("*/new/estimates.json"):
        bench_name = estimates_file.parent.parent.name
        match = PATTERN.match(bench_name)
        if not match:
            continue

        op, n_str, impl = match.groups()
        with estimates_file.open("r", encoding="utf-8") as handle:
            estimates = json.load(handle)

        mean_ns = float(estimates["mean"]["point_estimate"])
        means[(op, int(n_str), impl)] = mean_ns

    return means


def format_ns(value: float) -> str:
    if value >= 1_000_000_000.0:
        return f"{value / 1_000_000_000.0:.3f} s"
    if value >= 1_000_000.0:
        return f"{value / 1_000_000.0:.3f} ms"
    if value >= 1_000.0:
        return f"{value / 1_000.0:.3f} µs"
    return f"{value:.3f} ns"


def collect_rows(means: dict[tuple[str, int, str], float]) -> list[dict[str, object]]:
    keys = sorted({(op, n) for (op, n, _impl) in means})
    rows: list[dict[str, object]] = []

    for op, n in keys:
        baselines = {
            "f64": means.get((op, n, "nalgebra")),
            "f32": means.get((op, n, "nalgebra_f32")),
        }

        for dtype, nalgebra in baselines.items():
            if nalgebra is None:
                continue

            if dtype == "f64":
                impls = sorted(
                    impl
                    for (op_k, n_k, impl) in means
                    if op_k == op and n_k == n and impl != "nalgebra" and not impl.endswith("_f32")
                )
            else:
                impls = sorted(
                    impl
                    for (op_k, n_k, impl) in means
                    if op_k == op and n_k == n and impl != "nalgebra_f32" and impl.endswith("_f32")
                )

            for impl in impls:
                impl_time = means[(op, n, impl)]
                ratio = impl_time / nalgebra
                rows.append(
                    {
                        "op": op,
                        "n": n,
                        "dtype": dtype,
                        "impl": impl,
                        "impl_time": impl_time,
                        "nalgebra_time": nalgebra,
                        "ratio": ratio,
                        "delta_pct": (ratio - 1.0) * 100.0,
                        "speedup": nalgebra / impl_time,
                    }
                )

            aether_candidates: list[tuple[str, float]] = []
            for impl in impls:
                lower = impl.lower()
                if "nalgebra" in lower:
                    continue
                if "naive" in lower:
                    continue
                if "fma" in lower:
                    continue
                aether_candidates.append((impl, means[(op, n, impl)]))

            if aether_candidates:
                best_impl, best_time = min(aether_candidates, key=lambda item: item[1])
                ratio = best_time / nalgebra
                rows.append(
                    {
                        "op": op,
                        "n": n,
                        "dtype": dtype,
                        "impl": "aether_baseline",
                        "impl_time": best_time,
                        "nalgebra_time": nalgebra,
                        "ratio": ratio,
                        "delta_pct": (ratio - 1.0) * 100.0,
                        "speedup": nalgebra / best_time,
                        "source_impl": best_impl,
                    }
                )

    return rows


def generate_plots(rows: list[dict[str, object]], plot_dir: Path) -> tuple[int, str | None]:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return 0, "matplotlib is not installed (pip install matplotlib)"

    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["op"]), str(row["dtype"]))].append(row)

    def impl_suffix(dtype: str) -> str:
        return "_f32" if dtype == "f32" else ""

    def pick_present(candidates: list[str], available: set[str]) -> str | None:
        for name in candidates:
            if name in available:
                return name
        return None

    plot_dir.mkdir(parents=True, exist_ok=True)
    created = 0

    for (op, dtype), items in grouped.items():
        impl_to_n_to_time: dict[str, dict[int, float]] = defaultdict(dict)
        nalgebra_to_n: dict[int, float] = {}
        nalgebra_name = "nalgebra" if dtype == "f64" else "nalgebra_f32"

        for item in items:
            n = int(item["n"])
            impl = str(item["impl"])
            impl_to_n_to_time[impl][n] = float(item["impl_time"])
            nalgebra_to_n[n] = float(item["nalgebra_time"])

        if not nalgebra_to_n:
            continue

        ns = sorted(nalgebra_to_n.keys())
        all_impls = set(impl_to_n_to_time.keys())

        suffix = impl_suffix(dtype)
        native_impl = pick_present([f"native{suffix}", f"scalar{suffix}", f"aether{suffix}"], all_impls)
        simd_candidates = [f"simd{suffix}"]
        if op == "dot":
            simd_candidates = [f"simd{suffix}", f"aether{suffix}"]
        simd_impl = pick_present(simd_candidates, all_impls)

        selected_impls: list[tuple[str, str]] = []
        if native_impl is not None:
            selected_impls.append((native_impl, "aether_native"))
        if simd_impl is not None and simd_impl != native_impl:
            selected_impls.append((simd_impl, "aether_simd"))

        if not selected_impls:
            continue

        fig, ax = plt.subplots(figsize=(max(10, len(ns) * 1.4), 5))

        group_width = 0.82
        bar_labels = [nalgebra_name] + [alias for _, alias in selected_impls]
        bar_width = group_width / max(1, len(bar_labels))
        x_positions = list(range(len(ns)))

        for idx, label in enumerate(bar_labels):
            offset = (idx - (len(bar_labels) - 1) / 2.0) * bar_width
            xs = [x + offset for x in x_positions]

            if label == nalgebra_name:
                ys = [nalgebra_to_n[n] for n in ns]
                color = "#808080"
            else:
                source_impl = next(src for src, alias in selected_impls if alias == label)
                ys = [impl_to_n_to_time[source_impl].get(n, float("nan")) for n in ns]
                color = None

            ax.bar(xs, ys, width=bar_width * 0.95, label=label, color=color)

        ax.set_title(f"{op} runtime histogram ({dtype}) — lower is better")
        ax.set_xlabel("Problem size (N)")
        ax.set_ylabel("Mean time (ns)")
        ax.set_xticks(x_positions)
        ax.set_xticklabels([str(n) for n in ns])
        ax.grid(True, axis="y", alpha=0.25)
        ax.legend(loc="best")

        out_file = plot_dir / f"{op}_{dtype}_time_hist.png"
        fig.tight_layout()
        fig.savefig(out_file, dpi=150)
        plt.close(fig)
        created += 1

        # Add an explicit summary view: best Aether implementation vs Nalgebra.
        best_aether_to_n: dict[int, float] = {}
        if "aether_baseline" in impl_to_n_to_time:
            best_aether_to_n = dict(impl_to_n_to_time["aether_baseline"])
        else:
            for n in ns:
                candidates = [
                    impl_to_n_to_time[impl][n]
                    for impl, _alias in selected_impls
                    if n in impl_to_n_to_time[impl]
                ]
                if candidates:
                    best_aether_to_n[n] = min(candidates)

        if best_aether_to_n:
            ns_best = [n for n in ns if n in best_aether_to_n]
            fig2, ax2 = plt.subplots(figsize=(max(10, len(ns_best) * 1.25), 5))

            x2 = list(range(len(ns_best)))
            width2 = 0.36
            nalgebra_vals = [nalgebra_to_n[n] for n in ns_best]
            aether_vals = [best_aether_to_n[n] for n in ns_best]

            ax2.bar([x - width2 / 2 for x in x2], nalgebra_vals, width=width2, label=nalgebra_name, color="#808080")
            ax2.bar([x + width2 / 2 for x in x2], aether_vals, width=width2, label="aether_best")

            ax2.set_title(f"{op} runtime ({dtype}) — aether best vs nalgebra (lower is better)")
            ax2.set_xlabel("Problem size (N)")
            ax2.set_ylabel("Mean time (ns)")
            ax2.set_xticks(x2)
            ax2.set_xticklabels([str(n) for n in ns_best])
            ax2.grid(True, axis="y", alpha=0.25)
            ax2.legend(loc="best")

            out_file2 = plot_dir / f"{op}_{dtype}_aether_vs_nalgebra_time.png"
            fig2.tight_layout()
            fig2.savefig(out_file2, dpi=150)
            plt.close(fig2)
            created += 1

    return created, None


def main() -> int:
    args = parse_args()
    means = load_means()
    if not means:
        print("No matching Criterion results found in target/criterion.")
        print("Run: cargo bench -p aether --bench throughput")
        return 1

    rows = collect_rows(means)

    print("Implementations vs Nalgebra (mean time)")
    print("-" * 96)
    print(f"{'Benchmark':<20} {'Type':<6} {'Impl':<12} {'Impl Time':>14} {'Nalgebra':>14} {'Ratio':>10} {'Delta':>12}")
    print("-" * 96)

    for row in rows:
        label = f"{row['op']}_{row['n']}"
        print(
            f"{label:<20} {row['dtype']:<6} {row['impl']:<12} {format_ns(float(row['impl_time'])):>14} {format_ns(float(row['nalgebra_time'])):>14} {float(row['ratio']):>10.3f}x {float(row['delta_pct']):>+11.2f}%"
        )

    if not rows:
        print("No paired implementation/Nalgebra benchmarks were found.")
        print("Expected names like: matvec_1000_basic + matvec_1000_nalgebra (or *_f32 variants)")
        return 1

    if not args.no_plots:
        created, warning = generate_plots(rows, args.plot_dir)
        if warning is not None:
            print(f"\nPlot generation skipped: {warning}")
        else:
            print(f"\nWrote {created} histogram plot(s) to: {args.plot_dir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
