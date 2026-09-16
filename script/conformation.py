import numpy as np
import matplotlib.pyplot as plt

# ----- Data -----
shapes = ["Cylinder", "Tapered", "Diamond", "Bowtie", "Football", "Hourglass"]

ree = np.array([169.559636, 176.351045, 180.669608, 170.435317, 183.194756, 174.420206])
rg  = np.array([78.5394542, 81.5139244, 86.747016, 77.8279196, 89.1120146, 80.2875703])

reestd = np.array([1.39391127, 1.59475098, 1.41606885, 1.26387602, 0.93857814, 1.18161152])
rgstd  = np.array([0.59101864, 0.49159709, 0.47267525, 0.47939683, 0.31809061, 0.49941822])

x = np.arange(len(shapes))

# ----- Plot -----
fig, ax1 = plt.subplots(figsize=(10,8))

# Left axis (Ree)
ax1.errorbar(x, ree, yerr=reestd,
             fmt='o', markersize=15,
             capsize=4, label=r'$\langle R_{ee}^{2} \rangle ^{1/2}$',color = 'Dodgerblue')
ax1.set_ylabel(r'$\langle R_{ee}^{2} \rangle ^{1/2}$', fontsize=25)
ax1.set_xticks(x)
ax1.set_xticklabels(shapes, rotation=30)

# Right axis (Rg)
ax2 = ax1.twinx()
ax2.errorbar(x, rg, yerr=rgstd,
             fmt='o', markersize=15,
             capsize=4, label=r'$\langle R_{g}^{2} \rangle ^{1/2}$',color='Purple')
ax2.set_ylabel(r'$\langle R_{g}^{2} \rangle ^{1/2}$', fontsize=25)
ax1.tick_params(axis='both', which='major', length=6, direction='in', width=1.6, top=True, labelsize=20)
ax1.tick_params(axis='both', which='minor', length=4, direction='in', width=1.6, top=True, labelsize=20)
ax1.xaxis.set_tick_params(pad=10)
ax1.yaxis.set_tick_params(pad=10)
ax2.tick_params(axis='both', which='major', length=6, direction='in', width=1.6, top=True, labelsize=20)
ax2.tick_params(axis='both', which='minor', length=4, direction='in', width=1.6, top=True, labelsize=20)
ax2.xaxis.set_tick_params(pad=10)
ax2.yaxis.set_tick_params(pad=10)
ax1.set_ylim(160, 200)
ax2.set_ylim(60, 100)
# Clean layout
# ax1.set_xlabel("Shape")
# ax1.grid(alpha=0.3)

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1+h2 , l1+l2, ncol=3, frameon=False, fontsize=24, loc="upper right")
fig.tight_layout()
plt.savefig("conformation.pdf")
plt.show()