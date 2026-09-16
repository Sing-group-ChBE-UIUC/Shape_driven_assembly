import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import root_scalar
from pathlib import Path

# ----------------------------
# Base football profile
# ----------------------------
def Nsc_football(x):
    return -0.0245*(x/1.2)**2 + 2.4735*(x/1.2) + 7.551

# ----------------------------
# Football -> cylinder transition
# t = 0 : football
# t = 1 : cylinder y = 40
# ----------------------------
def make_Nsc_transition(t):
    def Nsc(x):
        return (1.0 - t) * Nsc_football(x) + t * 40.0
    return Nsc

# ----------------------------
# Solve D_list for one profile
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

        # stop if the next bead would go beyond x_max
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
# Parameters
# ----------------------------
alpha = 3.578
x_max = 130.0
nchange = 10
transition_list = np.linspace(0.0, 1.0, nchange)

base_path = Path("/Users/haisukang/remote3/Workspace/SCMF_BB/shape/v3_evenNisc_longbb/football_to_cylinder_full")
base_path.mkdir(parents=True, exist_ok=True)

# ----------------------------
# Main loop
# ----------------------------
for nth, t in enumerate(transition_list[4:]):
    Nsc_func = make_Nsc_transition(t)

    plt.figure(figsize=(6, 6))
    ax = plt.subplot(1, 1, 1)

    x_start0 = -7
    # first solve
    x_list, D_list, x_end = solve_D_list_for_profile(
        Nsc_func=Nsc_func,
        x_start0=x_start0,
        x_max=x_max,
        alpha=alpha
    )

    # optional adjustment to improve how close final end gets to x_max
    for _ in range(8):
        if len(D_list) == 0:
            break

        if x_end > x_max:
            # dev = x_max - x_end
            # x_start0_adj = x_start0 + 0.5 * dev
            x_start0_adj = x_start0 -1

            x_list_new, D_list_new, x_end_new = solve_D_list_for_profile(
                Nsc_func=Nsc_func,
                x_start0=x_start0_adj,
                x_max=x_max,
                alpha=alpha
            )

            if abs(x_max - x_end_new) < abs(x_max - x_end):
                x_start0 = x_start0_adj
                x_list, D_list, x_end = x_list_new, D_list_new, x_end_new
            else:
                break
        else:
            break

    ISC_list = build_ISC_list(D_list)
 
    # print(f"shape {nth+1:02d}, t={t:.3f}, Nbead={len(D_list)}, Nisc={len(ISC_list)}, x_end={x_end:.4f}")

    # save directory
    mother_path = base_path / f"shape_{nth:02d}_t{t:.3f}"
    mother_path.mkdir(parents=True, exist_ok=True)

    # save D_cfgs
    file_path = mother_path / f"D_cfg_t{t:.3f}.txt"
    with open(file_path, 'w') as file:
        for i in range(len(ISC_list)):
            file.write("isc%d %.2lf %.4lf\n" % (i+1, ISC_list[i], 2.2534*ISC_list[i] - 6.7449))

    # plot contour
    x_vals = np.linspace(-7, x_max, 400)
    y_vals = Nsc_func(x_vals)

    ax.plot(x_vals, y_vals, '--', color='black', label=r'$N_{sc}(x)$')
    ax.plot(x_vals, alpha*np.sqrt(y_vals), '-', color='purple', label=r'$\alpha N_{sc}^{1/2}$')

    # draw circles
    for x, D in zip(x_list, D_list):
        circle = plt.Circle((x + D/2, D/2), D/2, color='lightblue', ec='purple', fill=False)
        ax.add_patch(circle)

    ax.set_xlim(-7, x_max + 2)
    ax.set_ylim(-5, 80)
    ax.set_xlabel(r"$N_{bb}$")
    ax.set_ylabel(r"$N_{sc}$ or $D$")
    ax.set_aspect('equal')
    ax.legend()
    ax.set_title(f"football → cylinder, t={t:.2f}")

    plt.savefig(mother_path / f"ISC_t{t:.3f}.pdf", bbox_inches="tight")
    plt.show()
    plt.close()