import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.patches import Rectangle, Circle
import argparse

def animate_cartpole(csv_path: str,
                     pole_length: float = 1.0,
                     fps: int = 50,
                     save_gif: str | None = None,
                     save_mp4: str | None = None):
    """
    Animate a cart–pole from CSV.
    Track & axes are STATIC; only the cart/pole move horizontally.
    CSV columns: t,x,xdot,theta,thetadot,u,cart_KE,pole_KE,potential_E,total_E
    theta = 0 -> upright (y up).
    """
    df = pd.read_csv(csv_path)
    t     = df["t"].to_numpy()
    x     = df["x"].to_numpy()
    theta = df["theta"].to_numpy()
    frames = len(t)

    # Geometry (meters)
    L = float(pole_length)
    pole_draw_len = 2.0 * L            # draw full pole if L was COM distance
    cart_w, cart_h = 0.6, 0.3
    wheel_r = 0.08
    track_y = 0.0
    cart_top_y = track_y + cart_h/2.0

    # Pole tip (upright convention)
    px = x + pole_draw_len * np.sin(theta)
    py = cart_top_y + pole_draw_len * np.cos(theta)

    # Static view bounds (do NOT update during animation)
    pad = 2.0
    x_lo, x_hi = x.min() - pad, x.max() + pad
    y_lo, y_hi = -0.6, cart_top_y + pole_draw_len + 0.6

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(x_lo, x_hi)     # fixed limits -> track is static
    ax.set_ylim(y_lo, y_hi)
    ax.set_xlabel("Horizontal position  x  [m]")
    ax.set_ylabel("Vertical position  y  [m]")
    ax.set_title("Cart–Pole MPC (swing-up + translate)")

    # Static track line (extend beyond view so it looks infinite)
    track_margin = 50.0
    track_line, = ax.plot([x_lo - track_margin, x_hi + track_margin],
                          [track_y, track_y], "k-", lw=1)

    # Cart body (rectangle)
    cart = Rectangle((x[0] - cart_w/2.0, track_y - cart_h/2.0), cart_w, cart_h,
                     facecolor="tab:blue", edgecolor="k", zorder=3)
    ax.add_patch(cart)

    # Wheels (purely visual, move with cart)
    wheel_left  = Circle((x[0] - cart_w/3.0, track_y - cart_h/2.0 - wheel_r), wheel_r,
                         facecolor="k", edgecolor="k", zorder=4)
    wheel_right = Circle((x[0] + cart_w/3.0, track_y - cart_h/2.0 - wheel_r), wheel_r,
                         facecolor="k", edgecolor="k", zorder=4)
    ax.add_patch(wheel_left)
    ax.add_patch(wheel_right)

    # Pole
    pole_line, = ax.plot([x[0], px[0]], [cart_top_y, py[0]],
                         color="tab:red", lw=2.5, solid_capstyle="round", zorder=2)

    # Time label (inside axes; no clipping)
    time_text = ax.text(
        0.01, 0.98, "", transform=ax.transAxes,
        ha="left", va="top",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7),
        clip_on=False, zorder=5
    )

    def init():
        cart.set_x(x[0] - cart_w/2.0)
        wheel_left.center  = (x[0] - cart_w/3.0, track_y - cart_h/2.0 - wheel_r)
        wheel_right.center = (x[0] + cart_w/3.0, track_y - cart_h/2.0 - wheel_r)
        pole_line.set_data([x[0], px[0]], [cart_top_y, py[0]])
        time_text.set_text(f"t = {t[0]:.2f}s")
        return cart, wheel_left, wheel_right, pole_line, time_text, track_line

    ani = None  # set after animate is defined

    def animate(i):
        # Move cart body and wheels horizontally; track/axes do NOT move
        cart.set_x(x[i] - cart_w/2.0)
        wheel_left.center  = (x[i] - cart_w/3.0, track_y - cart_h/2.0 - wheel_r)
        wheel_right.center = (x[i] + cart_w/3.0, track_y - cart_h/2.0 - wheel_r)

        # Update pole geometry
        pole_line.set_data([x[i], px[i]], [cart_top_y, py[i]])

        # Update time text
        time_text.set_text(f"t = {t[i]:.2f}s")

        # Hard-stop timer at final frame so saves finish
        if i == frames - 1 and ani is not None:
            ani.event_source.stop()

        return cart, wheel_left, wheel_right, pole_line, time_text, track_line

    interval_ms = int(1000 / fps)
    ani = animation.FuncAnimation(
        fig, animate, init_func=init,
        frames=frames, interval=interval_ms,
        blit=True, repeat=False  # do not loop
    )

    # Optional save (non-looping)
    if save_gif:
        ani.save(save_gif, writer=animation.PillowWriter(fps=fps, metadata={"loop": 1}))
        print(f"Saved GIF: {save_gif}")

    if save_mp4:
        ani.save(
            save_mp4,
            writer=animation.FFMpegWriter(
                fps=fps, bitrate=3000, codec="libx264",
                extra_args=["-preset", "veryfast", "-movflags", "+faststart"]
            )
        )
        print(f"Saved MP4: {save_mp4}")

    plt.tight_layout()
    plt.show()
    plt.close(fig)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", default="cartpole_mpc.csv")
    parser.add_argument("--L", type=float, default=1.0, help="Pole COM length; drawn pole ≈ 2*L")
    parser.add_argument("--fps", type=int, default=50)
    parser.add_argument("--gif", default=None)
    parser.add_argument("--mp4", default=None)
    args = parser.parse_args()

    animate_cartpole(
        csv_path=args.csv,
        pole_length=args.L,
        fps=args.fps,
        save_gif=args.gif,
        save_mp4=args.mp4,
    )
