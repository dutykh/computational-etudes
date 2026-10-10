#!/usr/bin/env python3
"""
polar_grid_geometry.py - Grid geometry and the cost of clustering

Computational Étude 12.1: Compares naive (r ∈ [0,1]) and doubled
(r ∈ [-1,1]) Chebyshev–Fourier grids on the unit disk. Visualises
grid clustering, area distribution, and CFL time-step scaling.

Chapter 12: Spectral Methods in Polar Coordinates
Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
Last modified: February 2026
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'ch07'))
from cheb_matrix import cheb_matrix

# Output directory
fig_dir = os.path.join(os.path.dirname(__file__),
                       '..', '..', '..', 'textbook', 'figures', 'ch12', 'python')
os.makedirs(fig_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Figure 12.1: Naive vs doubled polar grids on the disk
# ---------------------------------------------------------------------------

Nr = 25  # must be odd for doubled grid
M = 20   # angular points

# Chebyshev points on [-1, 1]
_, x_cheb = cheb_matrix(Nr)

# --- Naive grid: map x ∈ [-1,1] to r ∈ [0,1] via r = (x+1)/2 ---
r_naive = (x_cheb + 1) / 2  # includes r=0 and r=1
theta = np.linspace(0, 2 * np.pi, M, endpoint=False)
RR_n, TT_n = np.meshgrid(r_naive, theta)
XX_n = RR_n * np.cos(TT_n)
YY_n = RR_n * np.sin(TT_n)

# --- Doubled grid: r ∈ [-1,1], use only r > 0, odd Nr so no r=0 ---
r_doubled = x_cheb  # full Chebyshev grid on [-1,1]
# Keep only positive-r points (interior, excluding r=1 boundary)
N2 = (Nr - 1) // 2
r_pos = r_doubled[1:N2+1]  # r > 0 interior points
RR_d, TT_d = np.meshgrid(r_pos, theta)
XX_d = RR_d * np.cos(TT_d)
YY_d = RR_d * np.sin(TT_d)

# Median radius (circle containing half the points)
r_naive_interior = r_naive[1:-1]  # exclude endpoints
r_median_naive = np.median(r_naive_interior)
r_median_doubled = np.median(r_pos)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.8))

# Draw disk boundary
circle = plt.Circle((0, 0), 1.0, fill=False, color='k', linewidth=1.2)
ax1.add_patch(circle)
ax1.plot(XX_n.flatten(), YY_n.flatten(), 'k.', markersize=2.5)
# Median circle
circ_med = plt.Circle((0, 0), r_median_naive, fill=False, color='C3',
                       linewidth=1.0, linestyle='--')
ax1.add_patch(circ_med)
area_frac_naive = r_median_naive**2 * 100
ax1.set_title(f'Naive grid ($r \\in [0, 1]$)\n'
              f'Median circle: {area_frac_naive:.0f}% of area',
              fontsize=10)
ax1.set_xlim(-1.15, 1.15)
ax1.set_ylim(-1.15, 1.15)
ax1.set_aspect('equal')
ax1.axis('off')

circle2 = plt.Circle((0, 0), 1.0, fill=False, color='k', linewidth=1.2)
ax2.add_patch(circle2)
ax2.plot(XX_d.flatten(), YY_d.flatten(), 'k.', markersize=2.5)
circ_med2 = plt.Circle((0, 0), r_median_doubled, fill=False, color='C3',
                        linewidth=1.0, linestyle='--')
ax2.add_patch(circ_med2)
area_frac_doubled = r_median_doubled**2 * 100
ax2.set_title(f'Doubled grid ($r \\in [-1, 1]$)\n'
              f'Median circle: {area_frac_doubled:.0f}% of area',
              fontsize=10)
ax2.set_xlim(-1.15, 1.15)
ax2.set_ylim(-1.15, 1.15)
ax2.set_aspect('equal')
ax2.axis('off')

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'polar_grids.pdf'), bbox_inches='tight')
plt.close()

# ---------------------------------------------------------------------------
# Figure 12.2: Area distribution per ring and CFL scaling
# ---------------------------------------------------------------------------

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2))

# --- Left: area per radial ring ---
# Naive grid
r_sorted_naive = np.sort(r_naive)
areas_naive = []
for i in range(len(r_sorted_naive) - 1):
    a = np.pi * (r_sorted_naive[i+1]**2 - r_sorted_naive[i]**2)
    areas_naive.append(a)
midpoints_naive = 0.5 * (r_sorted_naive[:-1] + r_sorted_naive[1:])

# Doubled grid: use sorted positive-r points, plus r=0 and r=1 as boundaries
r_sorted_doubled = np.sort(np.concatenate(([0], r_pos, [1])))
areas_doubled = []
for i in range(len(r_sorted_doubled) - 1):
    a = np.pi * (r_sorted_doubled[i+1]**2 - r_sorted_doubled[i]**2)
    areas_doubled.append(a)
midpoints_doubled = 0.5 * (r_sorted_doubled[:-1] + r_sorted_doubled[1:])

width = 0.012
ax1.bar(midpoints_naive - width, areas_naive, width=2*width, alpha=0.7,
        label='Naive', color='C0')
ax1.bar(midpoints_doubled + width, areas_doubled, width=2*width, alpha=0.7,
        label='Doubled', color='C3')
ax1.set_xlabel('$r$')
ax1.set_ylabel('Area of radial annulus')
ax1.legend(fontsize=9)
ax1.set_title('Area per ring', fontsize=10)

# --- Right: CFL time step vs Nr ---
Nr_values = np.arange(7, 52, 2)  # odd values only
dt_naive_list = []
dt_doubled_list = []

for N in Nr_values:
    _, x = cheb_matrix(N)
    # Naive: r = (x+1)/2 ∈ [0,1]
    r_n = (x + 1) / 2
    h_min_naive = np.min(np.abs(np.diff(r_n)))
    dt_naive_list.append(h_min_naive**2)  # CFL ~ h_min^2

    # Doubled: use positive interior points
    n2 = (N - 1) // 2
    r_p = x[1:n2+1]
    h_min_doubled = np.min(np.abs(np.diff(np.sort(r_p))))
    dt_doubled_list.append(h_min_doubled**2)

ax2.loglog(Nr_values, dt_naive_list, 'o-', markersize=3.5,
           label='Naive', color='C0')
ax2.loglog(Nr_values, dt_doubled_list, 's-', markersize=3.5,
           label='Doubled', color='C3')

# Reference slopes
Nr_ref = Nr_values.astype(float)
ax2.loglog(Nr_ref, 5e2 * Nr_ref**(-4), 'k--', linewidth=0.8, alpha=0.5,
           label='$O(N^{-4})$')
ax2.loglog(Nr_ref, 2e1 * Nr_ref**(-2), 'k:', linewidth=0.8, alpha=0.5,
           label='$O(N^{-2})$')

ax2.set_xlabel('$N_r$')
ax2.set_ylabel('$\\Delta t_{\\max} \\propto h_{\\min}^2$')
ax2.legend(fontsize=8)
ax2.set_title('CFL time-step scaling', fontsize=10)

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'polar_area_cfl.pdf'), bbox_inches='tight')
plt.close()

print("Figures saved to:", fig_dir)
print("  polar_grids.pdf")
print("  polar_area_cfl.pdf")
