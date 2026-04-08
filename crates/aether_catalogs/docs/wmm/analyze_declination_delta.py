#!/usr/bin/env python3
"""Report maximum adjacent magnetic-north deltas from a WMM CSV grid."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


@dataclass(frozen=True)
class GridPoint:
    latitude_deg: float
    longitude_deg: float
    declination_deg: float
    horizontal_intensity_nt: float


@dataclass(frozen=True)
class DeltaResult:
    delta_deg: float
    axis: str
    start: GridPoint
    end: GridPoint


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze adjacent declination changes from a WMM CSV grid.")
    parser.add_argument("csv_path", type=Path, help="Path to the WMM CSV file")
    parser.add_argument(
        "--h-min",
        type=float,
        nargs="*",
        default=[0.0, 1000.0, 2000.0, 5000.0, 10000.0],
        help="Minimum horizontal field thresholds in nT for filtered reports",
    )
    parser.add_argument("--plot", action="store_true", help="Render and save a 2D declination-delta heatmap")
    parser.add_argument(
        "--plot-h-min",
        type=float,
        default=10000.0,
        help="Horizontal-field threshold in nT used for the heatmap plot",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional output image path for the heatmap plot",
    )
    return parser.parse_args()


def angular_difference_deg(a_deg: float, b_deg: float) -> float:
    return abs(((a_deg - b_deg + 180.0) % 360.0) - 180.0)


def read_grid(csv_path: Path):
    with csv_path.open("r", newline="") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise ValueError(f"{csv_path} contains no rows")

    latitudes = sorted({float(row["latitude_deg"]) for row in rows})
    longitudes = sorted({float(row["longitude_deg"]) for row in rows})

    grid: dict[tuple[float, float], GridPoint] = {}
    for row in rows:
        point = GridPoint(
            latitude_deg=float(row["latitude_deg"]),
            longitude_deg=float(row["longitude_deg"]),
            declination_deg=float(row["declination_deg"]),
            horizontal_intensity_nt=float(row["h_nt"]),
        )
        grid[(point.latitude_deg, point.longitude_deg)] = point

    return latitudes, longitudes, grid, rows[0]


def scan_axis(
    latitudes: list[float],
    longitudes: list[float],
    grid: dict[tuple[float, float], GridPoint],
    h_min_nt: float,
) -> tuple[DeltaResult | None, DeltaResult | None]:
    best_lat: DeltaResult | None = None
    for i in range(len(latitudes) - 1):
        for longitude_deg in longitudes:
            start = grid[(latitudes[i], longitude_deg)]
            end = grid[(latitudes[i + 1], longitude_deg)]
            if min(start.horizontal_intensity_nt, end.horizontal_intensity_nt) < h_min_nt:
                continue
            delta_deg = angular_difference_deg(end.declination_deg, start.declination_deg)
            candidate = DeltaResult(delta_deg, "lat", start, end)
            if best_lat is None or candidate.delta_deg > best_lat.delta_deg:
                best_lat = candidate

    best_lon: DeltaResult | None = None
    for latitude_deg in latitudes:
        for j in range(len(longitudes) - 1):
            start = grid[(latitude_deg, longitudes[j])]
            end = grid[(latitude_deg, longitudes[j + 1])]
            if min(start.horizontal_intensity_nt, end.horizontal_intensity_nt) < h_min_nt:
                continue
            delta_deg = angular_difference_deg(end.declination_deg, start.declination_deg)
            candidate = DeltaResult(delta_deg, "lon", start, end)
            if best_lon is None or candidate.delta_deg > best_lon.delta_deg:
                best_lon = candidate

    return best_lat, best_lon


def build_delta_grids(
    latitudes: list[float],
    longitudes: list[float],
    grid: dict[tuple[float, float], GridPoint],
    h_min_nt: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lat_delta = np.full((len(latitudes), len(longitudes)), np.nan, dtype=float)
    lon_delta = np.full((len(latitudes), len(longitudes)), np.nan, dtype=float)

    for i in range(len(latitudes) - 1):
        for j, longitude_deg in enumerate(longitudes):
            start = grid[(latitudes[i], longitude_deg)]
            end = grid[(latitudes[i + 1], longitude_deg)]
            if min(start.horizontal_intensity_nt, end.horizontal_intensity_nt) < h_min_nt:
                continue
            lat_delta[i, j] = angular_difference_deg(end.declination_deg, start.declination_deg)

    for i, latitude_deg in enumerate(latitudes):
        for j in range(len(longitudes) - 1):
            start = grid[(latitude_deg, longitudes[j])]
            end = grid[(latitude_deg, longitudes[j + 1])]
            if min(start.horizontal_intensity_nt, end.horizontal_intensity_nt) < h_min_nt:
                continue
            lon_delta[i, j] = angular_difference_deg(end.declination_deg, start.declination_deg)

    lat_filled = np.where(np.isfinite(lat_delta), lat_delta, -np.inf)
    lon_filled = np.where(np.isfinite(lon_delta), lon_delta, -np.inf)
    max_delta = np.maximum(lat_filled, lon_filled)
    max_delta[max_delta == -np.inf] = np.nan
    return lat_delta, lon_delta, max_delta


def default_plot_path(csv_path: Path, h_min_nt: float) -> Path:
    suffix = f"_declination_delta_hmin_{int(round(h_min_nt))}.png"
    return csv_path.with_name(f"{csv_path.stem}{suffix}")


def plot_delta_grids(
    csv_path: Path,
    output_path: Path,
    latitudes: list[float],
    longitudes: list[float],
    lat_delta: np.ndarray,
    lon_delta: np.ndarray,
    max_delta: np.ndarray,
    h_min_nt: float,
    metadata: dict[str, str],
) -> None:
    longitude_grid_deg, latitude_grid_deg = np.meshgrid(longitudes, latitudes)
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True, sharex=True, sharey=True)

    valid = np.concatenate([
        lat_delta[np.isfinite(lat_delta)],
        lon_delta[np.isfinite(lon_delta)],
        max_delta[np.isfinite(max_delta)],
    ])
    vmax = float(np.nanmax(valid)) if valid.size else 1.0

    panels = [
        (lat_delta, "Latitude-adjacent |ΔD|"),
        (lon_delta, "Longitude-adjacent |ΔD|"),
        (max_delta, "Max adjacent |ΔD|"),
    ]

    mesh = None
    for ax, (data, title) in zip(axes, panels):
        mesh = ax.pcolormesh(
            longitude_grid_deg,
            latitude_grid_deg,
            data,
            shading="auto",
            cmap="magma",
            vmin=0.0,
            vmax=vmax,
        )
        ax.set_title(title)
        ax.set_xlabel("Longitude (deg)")
        ax.grid(True, linestyle=":", linewidth=0.6, alpha=0.6)

    axes[0].set_ylabel("Latitude (deg)")
    colorbar = fig.colorbar(mesh, ax=axes, shrink=0.92, pad=0.02)
    colorbar.set_label("Adjacent declination delta (deg)")

    year = float(metadata["decimal_year"])
    altitude_km = float(metadata["altitude_km"])
    fig.suptitle(
        f"WMM2025 adjacent magnetic-north delta | year={year:.2f} | altitude={altitude_km:.0f} km | h_min={h_min_nt:.0f} nT",
        weight="bold",
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {output_path}")


def format_result(label: str, result: DeltaResult | None) -> str:
    if result is None:
        return f"{label}: no qualifying samples"

    return (
        f"{label}: {result.delta_deg:.6f} deg"
        f" | axis={result.axis}"
        f" | start=({result.start.latitude_deg:.1f}, {result.start.longitude_deg:.1f})"
        f" | end=({result.end.latitude_deg:.1f}, {result.end.longitude_deg:.1f})"
        f" | H=({result.start.horizontal_intensity_nt:.3f}, {result.end.horizontal_intensity_nt:.3f}) nT"
    )


def main() -> None:
    args = parse_args()
    latitudes, longitudes, grid, metadata = read_grid(args.csv_path)

    lat_step_deg = latitudes[1] - latitudes[0] if len(latitudes) > 1 else float("nan")
    lon_step_deg = longitudes[1] - longitudes[0] if len(longitudes) > 1 else float("nan")

    print(f"csv: {args.csv_path}")
    print(
        f"grid: {len(latitudes)} lat x {len(longitudes)} lon"
        f" | lat_step={lat_step_deg:.3f} deg"
        f" | lon_step={lon_step_deg:.3f} deg"
        f" | year={float(metadata['decimal_year']):.3f}"
        f" | altitude_km={float(metadata['altitude_km']):.3f}"
    )
    print()

    for h_min_nt in args.h_min:
        best_lat, best_lon = scan_axis(latitudes, longitudes, grid, h_min_nt)
        overall = max(
            [result for result in (best_lat, best_lon) if result is not None],
            key=lambda result: result.delta_deg,
            default=None,
        )

        print(f"h_min >= {h_min_nt:.1f} nT")
        print(format_result("  max_lat_adjacent", best_lat))
        print(format_result("  max_lon_adjacent", best_lon))
        print(format_result("  max_any_adjacent", overall))
        print()

    if args.plot:
        lat_delta, lon_delta, max_delta = build_delta_grids(latitudes, longitudes, grid, args.plot_h_min)
        output_path = args.output or default_plot_path(args.csv_path, args.plot_h_min)
        plot_delta_grids(
            args.csv_path,
            output_path,
            latitudes,
            longitudes,
            lat_delta,
            lon_delta,
            max_delta,
            args.plot_h_min,
            metadata,
        )


if __name__ == "__main__":
    main()
