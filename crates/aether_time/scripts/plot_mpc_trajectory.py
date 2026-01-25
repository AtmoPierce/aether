import csv
import matplotlib.pyplot as plt

path = "mpc_trajectory.csv"  # or set MPC_CSV before running `cargo test`

t, pos0, vel0, acc0, u0 = [], [], [], [], []
with open(path, newline="") as f:
    rdr = csv.DictReader(f)
    for row in rdr:
        t.append(float(row["t"]))
        pos0.append(float(row["pos_0"]))
        vel0.append(float(row["vel_0"]))
        acc0.append(float(row["acc_0"]))
        # last row has empty u's; guard it
        u0.append(float(row["u_0"])) if row["u_0"] else u0.append(float("nan"))

plt.figure(); plt.plot(t, pos0); plt.title("pos_0"); plt.xlabel("t [s]"); plt.ylabel("m")
plt.figure(); plt.plot(t, vel0); plt.title("vel_0"); plt.xlabel("t [s]"); plt.ylabel("m/s")
plt.figure(); plt.plot(t, acc0); plt.title("acc_0"); plt.xlabel("t [s]"); plt.ylabel("m/s²")
plt.figure(); plt.plot(t, u0);   plt.title("u_0");   plt.xlabel("t [s]"); plt.ylabel("cmd")

plt.show()
