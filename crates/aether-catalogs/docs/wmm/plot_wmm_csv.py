#!/usr/bin/env python3
"""Render a 3D globe from a WMM CSV grid."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm, colors
from matplotlib.lines import Line2D


REQUIRED_COLUMNS = {
    "decimal_year",
    "altitude_km",
    "latitude_deg",
    "longitude_deg",
    "x_nt",
    "y_nt",
    "z_nt",
    "h_nt",
    "f_nt",
}

MIN_3D_ARROW_STRIDE = 10


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot a 3D WMM globe from a generated CSV grid.")
    parser.add_argument("csv_path", type=Path, help="Path to the WMM CSV file")
    parser.add_argument("--output", type=Path, default=None, help="Optional output image path")
    parser.add_argument("--show", action="store_true", help="Display the figure interactively after rendering")
    parser.add_argument(
        "--view",
        choices=["2d", "3d"],
        default="2d",
        help="Render either a 2D longitude/latitude map or a 3D globe",
    )
    parser.add_argument("--field", choices=["f_nt", "h_nt", "z_nt"], default="f_nt", help="Scalar field for globe coloring")
    parser.add_argument("--cmap", default="jet", help="Matplotlib colormap name")
    parser.add_argument("--quiver-step", type=int, default=4, help="Grid subsampling for horizontal field arrows")
    parser.add_argument("--quiver-scale", type=float, default=0.14, help="Arrow length scale on the unit sphere")
    parser.add_argument("--elev", type=float, default=25.0, help="3D camera elevation in degrees")
    parser.add_argument("--azim", type=float, default=-120.0, help="3D camera azimuth in degrees")
    parser.add_argument(
        "--saa-threshold",
        type=float,
        default=None,
        help="Optional total-field threshold in nT for highlighting low-field regions",
    )
    return parser.parse_args()


def default_output_path(csv_path: Path, field_name: str) -> Path:
    return csv_path.with_name(f"{csv_path.stem}_{field_name}_3d.png")


def has_display() -> bool:
    return bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))


def horizontal_components_to_lon_lat(lat_rad, north_nt, east_nt):
    u = east_nt
    v = north_nt
    scale = np.cos(lat_rad)
    scale = np.where(np.abs(scale) < 1.0e-6, 1.0e-6, scale)
    return u / scale, v


def grid_spacing_degrees(values):
    if len(values) < 2:
        return None
    deltas = np.diff(values)
    positive_deltas = deltas[deltas > 0.0]
    if positive_deltas.size == 0:
        return None
    return float(np.min(positive_deltas))


def quiver_stride(latitudes_deg, longitudes_deg, requested_step, minimum_spacing_deg=None):
    lat_step = max(1, requested_step)
    lon_step = max(1, requested_step)

    if minimum_spacing_deg is not None:
        lat_spacing = grid_spacing_degrees(latitudes_deg)
        lon_spacing = grid_spacing_degrees(longitudes_deg)

        if lat_spacing is not None and lat_spacing > 0.0:
            lat_step = max(lat_step, int(np.ceil(minimum_spacing_deg / lat_spacing)))
        if lon_spacing is not None and lon_spacing > 0.0:
            lon_step = max(lon_step, int(np.ceil(minimum_spacing_deg / lon_spacing)))

    return lat_step, lon_step


def quiver_stride_3d(latitudes_deg, longitudes_deg, requested_step):
    lat_step, lon_step = quiver_stride(latitudes_deg, longitudes_deg, requested_step)
    return max(lat_step, MIN_3D_ARROW_STRIDE), max(lon_step, MIN_3D_ARROW_STRIDE)


def render_2d(
    latitudes_deg,
    longitudes_deg,
    latitude_grid_deg,
    longitude_grid_deg,
    lat_rad,
    fields,
    args,
):
    scalar = fields[args.field]
    cmap = plt.get_cmap(args.cmap)
    norm = colors.Normalize(vmin=np.nanmin(scalar), vmax=np.nanmax(scalar))

    fig, ax = plt.subplots(figsize=(14, 7), constrained_layout=True)
    mesh = ax.pcolormesh(longitudes_deg, latitudes_deg, scalar, shading="auto", cmap=cmap, norm=norm)

    lat_step, lon_step = quiver_stride(latitudes_deg, longitudes_deg, args.quiver_step)
    u, v = horizontal_components_to_lon_lat(
        lat_rad[::lat_step, ::lon_step],
        fields["x_nt"][::lat_step, ::lon_step],
        fields["y_nt"][::lat_step, ::lon_step],
    )
    magnitude = np.sqrt(u * u + v * v)
    magnitude[magnitude == 0.0] = 1.0
    u = u / magnitude
    v = v / magnitude
    ax.quiver(
        longitude_grid_deg[::lat_step, ::lon_step],
        latitude_grid_deg[::lat_step, ::lon_step],
        u,
        v,
        color="white",
        scale=40.0 / max(args.quiver_scale, 1.0e-6),
        width=0.0022,
        headwidth=3.0,
        headlength=4.0,
        headaxislength=3.5,
        pivot="mid",
    )

    if args.saa_threshold is not None:
        ax.contour(
            longitude_grid_deg,
            latitude_grid_deg,
            fields["f_nt"],
            levels=[args.saa_threshold],
            colors="black",
            linewidths=2.0,
        )
        ax.legend(
            handles=[Line2D([0], [0], color="black", linewidth=2.0, label=f"F = {args.saa_threshold:.0f} nT")],
            loc="lower left",
        )

    colorbar = fig.colorbar(mesh, ax=ax, shrink=0.92, pad=0.02)
    colorbar.set_label(f"{args.field} (nT)")

    decimal_year = fields["decimal_year"][0, 0]
    altitude_km = fields["altitude_km"][0, 0]
    ax.set_title(
        f"WMM2025 {args.field} at {altitude_km:.0f} km, year {decimal_year:.2f}",
        weight="bold",
    )
    ax.set_xlabel("Longitude (deg)")
    ax.set_ylabel("Latitude (deg)")
    ax.set_xlim(longitudes_deg.min(), longitudes_deg.max())
    ax.set_ylim(latitudes_deg.min(), latitudes_deg.max())
    ax.set_aspect("auto")
    ax.grid(True, linestyle=":", linewidth=0.7, alpha=0.65)
    return fig


def render_3d(latitudes_deg, longitudes_deg, lat_rad, lon_rad, x, y, z, fields, args):
    scalar = fields[args.field]
    cmap = plt.get_cmap(args.cmap)
    norm = colors.Normalize(vmin=np.nanmin(scalar), vmax=np.nanmax(scalar))
    facecolors = cmap(norm(scalar))

    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(
        x,
        y,
        z,
        facecolors=facecolors,
        rstride=1,
        cstride=1,
        linewidth=0,
        antialiased=False,
        shade=False,
        alpha=0.98,
    )

    lat_step, lon_step = quiver_stride_3d(latitudes_deg, longitudes_deg, args.quiver_step)
    qx = x[::lat_step, ::lon_step]
    qy = y[::lat_step, ::lon_step]
    qz = z[::lat_step, ::lon_step]
    vx, vy, vz = ned_horizontal_to_ecef(
        lat_rad[::lat_step, ::lon_step],
        lon_rad[::lat_step, ::lon_step],
        fields["x_nt"][::lat_step, ::lon_step],
        fields["y_nt"][::lat_step, ::lon_step],
    )
    magnitude = np.sqrt(vx * vx + vy * vy + vz * vz)
    magnitude[magnitude == 0.0] = 1.0
    vx = args.quiver_scale * vx / magnitude
    vy = args.quiver_scale * vy / magnitude
    vz = args.quiver_scale * vz / magnitude
    ax.quiver(qx, qy, qz, vx, vy, vz, color="white", linewidth=0.8, normalize=False)

    if args.saa_threshold is not None:
        saa_mask = fields["f_nt"] <= args.saa_threshold
        if np.any(saa_mask):
            ax.scatter(
                x[saa_mask] * 1.01,
                y[saa_mask] * 1.01,
                z[saa_mask] * 1.01,
                color="black",
                s=10,
                alpha=0.45,
                label=f"F <= {args.saa_threshold:.0f} nT",
            )
            ax.legend(loc="upper left")

    scalar_mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar_mappable.set_array([])
    colorbar = fig.colorbar(scalar_mappable, ax=ax, shrink=0.68, pad=0.08)
    colorbar.set_label(f"{args.field} (nT)")

    decimal_year = fields["decimal_year"][0, 0]
    altitude_km = fields["altitude_km"][0, 0]
    ax.set_title(f"WMM2025 {args.field} at {altitude_km:.0f} km, year {decimal_year:.2f}", pad=18, weight="bold")
    ax.set_box_aspect((1.0, 1.0, 1.0))
    ax.view_init(elev=args.elev, azim=args.azim)
    ax.set_axis_off()
    return fig


def read_grid(csv_path: Path):
    with csv_path.open("r", newline="") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise ValueError(f"{csv_path} contains no rows")

    missing = REQUIRED_COLUMNS.difference(rows[0].keys())
    if missing:
        raise ValueError(f"missing required columns: {sorted(missing)}")

    latitudes = sorted({float(row["latitude_deg"]) for row in rows})
    longitudes = sorted({float(row["longitude_deg"]) for row in rows})
    grid_shape = (len(latitudes), len(longitudes))

    if len(rows) != grid_shape[0] * grid_shape[1]:
        raise ValueError("CSV is not a complete latitude/longitude grid")

    lat_index = {value: idx for idx, value in enumerate(latitudes)}
    lon_index = {value: idx for idx, value in enumerate(longitudes)}
    fields = {
        key: np.full(grid_shape, np.nan, dtype=float)
        for key in ["x_nt", "y_nt", "z_nt", "h_nt", "f_nt", "decimal_year", "altitude_km"]
    }

    for row in rows:
        i = lat_index[float(row["latitude_deg"])]
        j = lon_index[float(row["longitude_deg"])]
        for key in fields:
            fields[key][i, j] = float(row[key])

    return np.array(latitudes), np.array(longitudes), fields


def ned_horizontal_to_ecef(lat_rad, lon_rad, north_nt, east_nt):
    north_unit_x = -np.sin(lat_rad) * np.cos(lon_rad)
    north_unit_y = -np.sin(lat_rad) * np.sin(lon_rad)
    north_unit_z = np.cos(lat_rad)

    east_unit_x = -np.sin(lon_rad)
    east_unit_y = np.cos(lon_rad)
    east_unit_z = np.zeros_like(lat_rad)

    vx = north_nt * north_unit_x + east_nt * east_unit_x
    vy = north_nt * north_unit_y + east_nt * east_unit_y
    vz = north_nt * north_unit_z + east_nt * east_unit_z
    return vx, vy, vz


def main() -> None:
    args = parse_args()
    latitudes_deg, longitudes_deg, fields = read_grid(args.csv_path)
    output_path = args.output or default_output_path(args.csv_path, args.field)
    if args.output is None:
        suffix = f"_{args.field}_{args.view}.png"
        output_path = args.csv_path.with_name(f"{args.csv_path.stem}{suffix}")

    latitude_grid_deg, longitude_grid_deg = np.meshgrid(latitudes_deg, longitudes_deg, indexing="ij")
    lat_rad = np.deg2rad(latitude_grid_deg)
    lon_rad = np.deg2rad(longitude_grid_deg)

    x = np.cos(lat_rad) * np.cos(lon_rad)
    y = np.cos(lat_rad) * np.sin(lon_rad)
    z = np.sin(lat_rad)

    if args.view == "2d":
        fig = render_2d(
            latitudes_deg,
            longitudes_deg,
            latitude_grid_deg,
            longitude_grid_deg,
            lat_rad,
            fields,
            args,
        )
    else:
        fig = render_3d(latitudes_deg, longitudes_deg, lat_rad, lon_rad, x, y, z, fields, args)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    print(f"wrote {output_path}")

    if args.show and has_display():
        plt.show()
    elif args.show:
        print("display not available; saved image instead")
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
