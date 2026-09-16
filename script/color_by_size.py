import numpy as np
import sys

def isc_to_rgb(isc_values):
    isc_values = np.array(isc_values)

    # Normalize isc values to 0–1
    isc_min = 28.43
    isc_max = 12.35
    norm = (isc_values - isc_min) / (isc_max - isc_min)

    # Color endpoints
    R = 255
    B = 102
    G_start = 102
    G_end = 255

    # Scale green channel
    G = G_start + norm * (G_end - G_start)

    # Build RGB list
    rgb_colors = [(R, g, B) for g in G]

    return rgb_colors

data = np.loadtxt(sys.argv[1], dtype={'names': ('index', 'D_i', 'kuhn'), 'formats': ('U10', 'f8', 'f8')})
D_i = data['D_i']

print(isc_to_rgb(D_i))