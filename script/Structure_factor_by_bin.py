#!/usr/bin/env python3
# plot_sq_shellavg.py
import sys, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fname = sys.argv[1] 

# choose ONE of these: set dq directly, or nbins
dq    = float(sys.argv[2]) if len(sys.argv) > 3 else None   # e.g. dq = 0.05
nbins = 200 if dq is None else None
sq_list = int(sys.argv[3])
# --- load ---
q, S = np.loadtxt(fname, usecols=(0,1), unpack=True)
m = np.isfinite(q) & np.isfinite(S)
q, S = q[m], S[m]

# --- make |q|-shell bins ---
qmin, qmax = q.min(), q.max()
if dq is None:
    edges = np.linspace(qmin, qmax, nbins+1)
else:
    edges = np.arange(qmin, qmax + dq, dq)

bin_id = np.digitize(q, edges) - 1
centers, means, stderrs, counts = [], [], [], []

for b in range(len(edges)-1):
    sel = (bin_id == b)
    if not np.any(sel): 
        continue
    centers.append(0.5*(edges[b]+edges[b+1]))
    Sb = S[sel]
    means.append(Sb.mean())
    # stderr = std / sqrt(N) (clip to avoid log(0))
    stderrs.append(Sb.std(ddof=1)/np.sqrt(Sb.size) if Sb.size>1 else 0.0)
    counts.append(Sb.size)

centers  = np.array(centers)
means    = np.array(means)
stderrs  = np.array(stderrs)
counts   = np.array(counts)

# --- plot ---
plt.figure(figsize=(7,3.5))
# plt.plot(q, S, lw=0.5, alpha=0.3, label='raw modes')  # optional
plt.errorbar(centers[1:], means[1:], fmt='-o', ms=3, lw=1.2) #label='|q|-binned avg')
# plt.ylim(0.1, 20)
plt.xscale('log'); plt.yscale('log')
# plt.axhline(1.0, ls='--', lw=0.8, color='gray')
plt.xlabel(r'$q[b^{-1}]$'); plt.ylabel(r'$S(q)$')
# plt.tight_layout(); plt.legend(frameon=False)
plt.savefig("Sq_%d.pdf"%(sq_list))
#plt.show()
plt.close()

# If you want to save the binned data:
np.savetxt("S_q%d_bin.txt"%(sq_list), np.c_[centers, means, stderrs, counts],
           header="q  S_avg  S_stderr  n_modes")
