import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import root_scalar
from pathlib import Path

# ----------------------------
# Base hourglass half-profile
# This is the HALF shape you defined
# ----------------------------
def Nsc_hourglass_half(x):
    return -0.1*(x/1.2)**2 + 5.1*(x/1.2) + 5.0

# ----------------------------
# Hourglass -> cylinder transition
# t = 0 : hourglass half-profile
# t = 1 : cylinder y = 40
# ----------------------------
def make_Nsc_transition(t):
    def Nsc(x):
        return (1.0 - t) * Nsc_hourglass_half(x) + t * 40.0
    return Nsc

# ----------------------------
# Solve D_list for one HALF profile
# D^3 = alpha^2 * integral_x^{x+D} Nsc(s) ds
# ----------------------------
def solve_D_list_for_profile(Nsc_func, x_start0, x_max, alpha):
    def make_equation(x_start):
        def f(D):
            val, _ = quad(Nsc_func, x_start, x_start + D)
            return D**3 - alpha**2 * val
        return f

    x_start = float(x_start0)
    D_list = []
    x_list = []

    while x_start < x_max:
        f = make_equation(x_start)

        try:
            sol = root_scalar(f, bracket=[0.1, 100.0], method="brentq")
        except ValueError:
            break

        if not sol.converged:
            break

        D_i = sol.root

        # stop if next bead goes beyond half-domain
        if x_start + D_i > x_max:
            break

        x_list.append(x_start)
        D_list.append(D_i)
        x_start += D_i

    if len(D_list) == 0:
        x_end = x_start0
    else:
        x_end = x_list[-1] + D_list[-1]

    return np.array(x_list), np.array(D_list), x_end

# ----------------------------
# Mirror half list to full list
# If you want a central duplicated bead:
# [a,b,c] -> [a,b,c,c,b,a]
# ----------------------------
def mirror_D_list(D_list_half):
    if len(D_list_half) == 0:
        return np.array([])
    return np.concatenate([D_list_half, D_list_half[::-1]])

# ----------------------------
# Insert middle ISC beads
# [d1,d2,d3] -> [d1,(d1+d2)/2,d2,(d2+d3)/2,d3]
# ----------------------------
def build_ISC_list(D_list):
    if len(D_list) == 0:
        return np.array([])

    Tot_ISC = 2 * len(D_list) - 1
    ISC_list = []

    for i in range(1, Tot_ISC + 1):
        if i % 2 != 0:
            D_ISC = D_list[i // 2]
        else:
            D_ISC = 0.5 * (D_list[(i - 1) // 2] + D_list[(i + 1) // 2])
        ISC_list.append(D_ISC)

    return np.array(ISC_list)

# ----------------------------
# Build x positions for plotting circles
# starting from center and expanding to both sides
# ----------------------------
def build_symmetric_x_list(D_list_half, x_center=0.0):
    if len(D_list_half) == 0:
        return np.array([]), np.array([])

    # right half
    x_right = []
    pos = x_center
    for D in D_list_half:
        x_right.append(pos)
        pos += D
    x_right = np.array(x_right)

    # left half: mirrored placement
    x_left = []
    pos = x_center
    for D in D_list_half:
        pos -= D
        x_left.append(pos)
    x_left = np.array(x_left[::-1])  # left to center order

    D_full = np.concatenate([D_list_half[::-1], D_list_half])
    x_full = np.concatenate([x_left, x_right])

    return x_full, D_full

# ----------------------------
# Parameters
# ----------------------------
alpha = 3.578
x_half_max = 70.0
nchange = 10
transition_list = np.linspace(0.0, 1.0, nchange)

base_path = Path("/Users/haisukang/remote3/Workspace/SCMF_BB/shape/v3_evenNisc_longbb/hourglass_to_cylinder_full")
base_path.mkdir(parents=True, exist_ok=True)

# ----------------------------
# Main loop
# ----------------------------
for nth, t in enumerate(transition_list[0:1]):
    Nsc_func = make_Nsc_transition(t)

    plt.figure(figsize=(6, 6))
    ax = plt.subplot(1, 1, 1)

    x_start0 = -1

    # solve only HALF
    x_list_half, D_list_half, x_end = solve_D_list_for_profile(
        Nsc_func=Nsc_func,
        x_start0=x_start0,
        x_max=x_half_max,
        alpha=alpha
    )

    # optional adjustment
    for _ in range(8):
        if len(D_list_half) == 0:
            break

        if x_end > x_half_max:
            x_start0_adj = x_start0 - 1.0

            x_list_half_new, D_list_half_new, x_end_new = solve_D_list_for_profile(
                Nsc_func=Nsc_func,
                x_start0=x_start0_adj,
                x_max=x_half_max,
                alpha=alpha
            )

            if abs(x_half_max - x_end_new) < abs(x_half_max - x_end):
                x_start0 = x_start0_adj
                x_list_half, D_list_half, x_end = x_list_half_new, D_list_half_new, x_end_new
            else:
                break
        else:
            break

    # mirror to full list
    print(D_list_half)
    D_list_half[-1] = D_list_half[-1]-1.5
    D_list_half[0] = D_list_half[0]-4
    D_list_half[1] = D_list_half[1]+1.45
    D_list = mirror_D_list(D_list_half)

    # build ISC from full list
    ISC_list = build_ISC_list(D_list)

    print(f"shape {nth:02d}, t={t:.3f}, Nbead_half={len(D_list_half)}, Nbead_full={len(D_list)}, Nisc={len(ISC_list)}")

    # save directory
    mother_path = base_path / f"shape_{nth:02d}_t{t:.3f}"
    mother_path.mkdir(parents=True, exist_ok=True)

    # save D_cfgs
    file_path = mother_path / f"D_cfg_t{t:.3f}.txt"
    with open(file_path, 'w') as file:
        for i in range(len(ISC_list)):
            file.write("isc%d %.2lf %.4lf\n" % (i+1, ISC_list[i], 2.2534*ISC_list[i] - 6.7449))

    # profile for plotting: mirrored half-profile
    x_vals_half = np.linspace(0.0, x_half_max, 400)
    y_vals_half = Nsc_func(x_vals_half)

    x_vals_full = np.concatenate([-x_vals_half[::-1], x_vals_half])
    y_vals_full = np.concatenate([y_vals_half[::-1], y_vals_half])

    ax.plot(x_vals_full, y_vals_full, '--', color='black', label=r'$N_{sc}(x)$')
    ax.plot(x_vals_full, alpha*np.sqrt(y_vals_full), '-', color='purple', label=r'$\alpha N_{sc}^{1/2}$')

    # draw circles
    x_plot, D_plot = build_symmetric_x_list(D_list_half, x_center=0.0)
    for x, D in zip(x_plot, D_plot):
        circle = plt.Circle((x + D/2, D/2), D/2, color='lightblue', ec='purple', fill=False)
        ax.add_patch(circle)

    ax.set_xlim(-x_half_max - 5, x_half_max + 5)
    ax.set_ylim(-5, 80)
    ax.set_xlabel(r"$N_{bb}$")
    ax.set_ylabel(r"$N_{sc}$ or $D$")
    ax.set_aspect('equal')
    ax.legend()
    ax.set_title(f"hourglass → cylinder, t={t:.2f}")

    plt.savefig(mother_path / f"ISC_t{t:.3f}.pdf", bbox_inches="tight")
    plt.show()
    plt.close()