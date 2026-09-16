#!/usr/bin/env python3
# plot_sq_Bonly_triblock_rpa.py
import sys
import numpy as np
import matplotlib.pyplot as plt
# -------------------------- simulation curve (B only) --------------------------
# Expect 2 columns: q  SBB(q)
# q_sim, SBB_sim = np.loadtxt(sys.argv[1], unpack=True)
q_sim, SBB_sim, S_stderr1, n_modes1 = np.loadtxt(sys.argv[1], skiprows=1, unpack=True)
data = np.loadtxt(sys.argv[2], dtype={'names': ('index', 'D_i', 'kuhn'), 'formats': ('U10', 'f8', 'f8')})
D_i = data['D_i']
mon_v_i = (D_i/2)**3
# fb = (np.sum(mon_v_i[1:3])+np.sum(mon_v_i[10:]))/np.sum(mon_v_i)
fb = (np.sum(mon_v_i[5:8]))/np.sum(mon_v_i)
fa = (1-fb)/2
plt.figure(figsize=(5, 5))
plt.plot(q_sim[1:], (SBB_sim[1:] - 0.44)*0.18, "-o", ms=3, lw=1.2,
         label=f"Simulation)") # bowtie
# plt.plot(q_sim[1:], (SBB_sim[1:] - 0.37)*0.42, "-o", ms=3, lw=1.2,
#          label=f"Simulation)") # football
# plt.plot(q_sim[1:], (SBB_sim[1:] - 0.35)*0.45, "-o", ms=3, lw=1.2,
#          label=f"Simulation)") # dia


# -----------------------------
# ABC triblock RPA (Gaussian chains) + incompressibility + scattering intensity
# Based on Cochran et al. (Macromolecules 2003): build ideal-gas S0_ij(q),
# then RPA: S = (S0^{-1} + chi_matrix)^{-1}, then enforce incompressibility,
# then I(q) = c^T S c.
# -----------------------------

def S0_triblock(q, N, f, b, nu0=1.0):
    """
    Ideal-gas correlation matrix S0_ij(q) for a linear ABC triblock.

    q   : array
    N   : total DP (or "RPA N")
    f   : [fA,fB,fC] (sum to 1)
    b   : [bA,bB,bC] statistical segment lengths
    nu0 : reference monomer volume (often set 1 in reduced units)

    Returns: (len(q), 3, 3) array
    """
    f = np.asarray(f, float)
    b = np.asarray(b, float)
    assert len(f) == 3 and len(b) == 3
    assert abs(f.sum() - 1.0) < 1e-10

    q = np.asarray(q, float)
    nq = q.size
    S0 = np.zeros((nq, 3, 3), float)

    # alpha_i(q) = q^2 * N * b_i^2 / 6
    alpha = (q[:, None] ** 2) * (N * (b[None, :] ** 2) / 6.0)  # (nq,3)

    # Diagonal (i=i):  2*(exp(-f_i*alpha_i) + f_i*alpha_i - 1)/alpha_i^2
    eps = 1e-300
    for i in range(3):
        ai = np.maximum(alpha[:, i], eps)
        fi = f[i]
        gii = 2.0 * (np.exp(-fi * ai) + fi * ai - 1.0) / (ai * ai)
        S0[:, i, i] = (N / nu0) * gii

    # Off-diagonal (i<j):
    # gij = exp(-gap) * (1-exp(-f_i alpha_i))*(1-exp(-f_j alpha_j)) / (alpha_i alpha_j)
    # gap = q^2 * N/6 * sum_{k=i+1}^{j-1} f_k b_k^2
    def gap_between(i, j):
        if j <= i + 1:
            return 0.0
        mid = np.arange(i + 1, j)
        return np.sum(f[mid] * b[mid] ** 2)

    for i in range(3):
        for j in range(i + 1, 3):
            ai = np.maximum(alpha[:, i], eps)
            aj = np.maximum(alpha[:, j], eps)
            fi, fj = f[i], f[j]

            gap_const = gap_between(i, j)  # scalar
            gap = (q ** 2) * (N / 6.0) * gap_const  # (nq,)

            gij = np.exp(-gap) * (1.0 - np.exp(-fi * ai)) * (1.0 - np.exp(-fj * aj)) / (ai * aj)
            S0[:, i, j] = (N / nu0) * gij
            S0[:, j, i] = (N / nu0) * gij

    return S0


def rpa_S(q, N, f, b, chi, nu0=1.0, enforce_incompressibility=True):
    """
    Build RPA structure-factor matrix S(q) for ABC.

    chi: 3x3 symmetric matrix (zeros on diagonal), in whatever convention you want.
         If you want the "2*chi" convention, just pass 2*chi here.
    """
    chi = np.asarray(chi, float)
    assert chi.shape == (3, 3)

    S0 = S0_triblock(q, N, f, b, nu0=nu0)
    nq = len(q)
    S = np.zeros_like(S0)

    for k in range(nq):
        S0k = S0[k]
        Mk = np.linalg.inv(S0k) + chi
        Sk = np.linalg.inv(Mk)

        if enforce_incompressibility:
            # Enforce constraint u^T δφ = 0 with u = (1,1,1)
            # Conditional covariance:
            # S_c = S - S u (u^T S u)^{-1} u^T S
            u = np.ones(3)
            Su = Sk @ u
            denom = float(u @ Su)
            Sk = Sk - np.outer(Su, Su) / denom

        S[k] = Sk

    return S


def intensity_from_contrast(S, c):
    """
    I(q) = c^T S(q) c
    S: (nq,3,3), c: length-3 vector
    """
    c = np.asarray(c, float).reshape(3, 1)
    I = np.einsum("kij,ia,jb->k", S, c, c).reshape(-1)
    return I


if __name__ == "__main__":
    # ---- Example usage ----
    q = np.linspace(1e-3, 2.0, 4000)

    # triblock composition
    f = [fa, fb, fa]          # fA,fB,fC
    N = 12
    b=22.6# total DP (RPA N)
    b_list = [b, b, b]          # statistical segment lengths (can differ!)

    # Interaction matrix chi_ij (symmetric, diagonal=0)
    # Put your temperature-dependent chis here.
    chi_i = float(sys.argv[3])
    chiAB, chiAC, chiBC = chi_i, 0.0, chi_i
    chi = np.array([[0.0,   chiAB, chiAC],
                    [chiAB, 0.0,   chiBC],
                    [chiAC, chiBC, 0.0 ]])

    # If you prefer the AB-diblock convention S^{-1} = S0^{-1} - 2*chi*...,
    # just replace with chi_eff = 2*chi (sign depends on your definition).
    chi_eff = chi

    S = rpa_S(q, N, f, b_list, chi_eff, nu0=1.0, enforce_incompressibility=True)

    # Contrast vectors c = (cA,cB,cC): choose what your scattering "sees"
    # Example 1: emphasize B vs (A,C)
    c1 = np.array([0, 1, 0])

    # Example 2: neutron-like contrast A - C (B invisible)
    c2 = np.array([1.0, 0.0, -1.0])

    I1 = intensity_from_contrast(S, c1)
    I2 = intensity_from_contrast(S, c2)
    plt.plot(q, I1, label="c=(0,1,0)  (B-only)")
    # plt.plot(q, I2, label="c=(1,0,-1) (A-C)")


# # diblock RPA no matrixN = 12
# chi = 0.885  # IMPORTANT: interpret this as chiN (not bare chi)
# N = 12
# b = 24
# # b = np.mean(D_i[:6])
# f = 0.5
# # chilist = np.linspace(0, 0.9, 10)
# q = np.linspace(1e-2, 10.0, 20000)   # 100 is overkill; start smaller
# Rg2 = N * b * b / 6.0
# x = np.maximum(q*q*Rg2, 1e-12)
# rho = 1.42

# # Avoid divide-by-zero issues at very small x
# # x = np.maximum(x, 1e-12)

# # Gaussian diblock intrachain correlators (dimensionless, Leibler)
# wAA = 2.0 * (np.exp(-f * x) - 1.0 + f * x) / (x * x)
# wBB = 2.0 * (np.exp(-(1.0 - f) * x) - 1.0 + (1.0 - f) * x) / (x * x)
# wAB = (1.0 - np.exp(-f * x)) * (1.0 - np.exp(-(1.0 - f) * x)) / (x * x)

# # Leibler kernel F(x,f)
# F = (wAA + wBB + 2.0 * wAB) / (wAA * wBB - wAB * wAB)
# # for i in range(len(chilist)):
# # If you define chi as chiN, just do:
# # Sq = rho*1.0 / ((F / N) - 2.0 * chilist[i])/4
# Sq = rho*1.0 / ((F) - 2.0 * chi*N)

# # plt.plot(q, Sq, label=rf'$\chi N = {chilist[i]*N:.1f}$')
# plt.plot(q, Sq, label=rf'RPA $\chi N = {chi*N:.1f}$')

# -------------------------- formatting --------------------------
plt.xlim(0, 1.2)
# plt.ylim(0, 4)
plt.text(0.5, 0.70, rf"$f_a={fa:.2f}$", fontsize=14, transform=plt.gca().transAxes)
plt.text(0.5, 0.65, rf"$N={N}$", fontsize=14, transform=plt.gca().transAxes)
plt.text(0.5, 0.60, rf"$b={b:.2f}$", fontsize=14, transform=plt.gca().transAxes)
plt.text(0.5, 0.55, rf"$\chi={chi_i}$", fontsize=14, transform=plt.gca().transAxes)

plt.xlabel(r"$q[b^{-1}]$")
plt.ylabel(r"$S(q)$")
plt.tight_layout()
plt.legend(frameon=False, loc="upper right", ncol=1)
plt.savefig("Sq_triblock_RPA_re2.pdf")
plt.show()