# plot_q_guidance.py
# Usage:
#   python plot_q_guidance.py path/to/q_guidance_250.csv
# Or plot many:
#   python plot_q_guidance.py path/to/q_guidance_*.csv

import sys
import glob
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

def plot_one(csv_path: Path):
    df = pd.read_csv(csv_path)

    # Basic sanity
    required = [
        "time","t_go","pos_x","pos_y","pos_z",
        "vel_x","vel_y","vel_z","q11","q22","q33","q12","q13","q23"
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{csv_path}: missing columns {missing}")

    # Title suffix (e.g., t_f from filename)
    title_suffix = f" ({csv_path.stem})"

    # 1) Position vs time
    plt.figure()
    plt.plot(df["time"], df["pos_x"], label="pos_x")
    plt.plot(df["time"], df["pos_y"], label="pos_y")
    plt.plot(df["time"], df["pos_z"], label="pos_z")
    plt.xlabel("time [s]")
    plt.ylabel("position [m]")
    plt.title("Position vs Time" + title_suffix)
    plt.grid(True, alpha=0.3)
    plt.legend()
    out = csv_path.with_suffix(".pos.png")
    plt.savefig(out, dpi=160, bbox_inches="tight")
    print(f"Saved {out}")

    # 2) Velocity vs time
    plt.figure()
    plt.plot(df["time"], df["vel_x"], label="vel_x")
    plt.plot(df["time"], df["vel_y"], label="vel_y")
    plt.plot(df["time"], df["vel_z"], label="vel_z")
    plt.xlabel("time [s]")
    plt.ylabel("velocity [m/s]")
    plt.title("Velocity vs Time" + title_suffix)
    plt.grid(True, alpha=0.3)
    plt.legend()
    out = csv_path.with_suffix(".vel.png")
    plt.savefig(out, dpi=160, bbox_inches="tight")
    print(f"Saved {out}")

    # 3) Q diagonals vs time
    plt.figure()
    plt.plot(df["time"], df["q11"], label="q11 (Q_xx)")
    plt.plot(df["time"], df["q22"], label="q22 (Q_yy)")
    plt.plot(df["time"], df["q33"], label="q33 (Q_zz)")
    plt.xlabel("time [s]")
    plt.ylabel("Q diagonal [1/s]")
    plt.title("Q Diagonals vs Time" + title_suffix)
    plt.grid(True, alpha=0.3)
    plt.legend()
    out = csv_path.with_suffix(".qdiag.png")
    plt.savefig(out, dpi=160, bbox_inches="tight")
    print(f"Saved {out}")

    # 4) Q off-diagonals vs time
    plt.figure()
    plt.plot(df["time"], df["q12"], label="q12 (Q_xy)")
    plt.plot(df["time"], df["q13"], label="q13 (Q_xz)")
    plt.plot(df["time"], df["q23"], label="q23 (Q_yz)")
    plt.xlabel("time [s]")
    plt.ylabel("Q off-diagonal [1/s]")
    plt.title("Q Off-Diagonals vs Time" + title_suffix)
    plt.grid(True, alpha=0.3)
    plt.legend()
    out = csv_path.with_suffix(".qoff.png")
    plt.savefig(out, dpi=160, bbox_inches="tight")
    print(f"Saved {out}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Please provide a CSV path or glob (e.g., q_guidance_*.csv)")
        sys.exit(2)

    patterns = sys.argv[1:]
    files = []
    for pat in patterns:
        files.extend(glob.glob(pat))
    if not files:
        print("No files matched.")
        sys.exit(1)

    for f in files:
        plot_one(Path(f))
