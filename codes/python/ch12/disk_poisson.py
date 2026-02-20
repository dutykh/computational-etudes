#!/usr/bin/env python3
"""
disk_poisson.py - Poisson equation on the unit disk

Computational Étude 12.4: Solves -Δu = f on the disk with u = 0 at r = 1,
using an off-centre Gaussian source. Also compares radially symmetric vs
asymmetric forcing to verify symmetry preservation.

Chapter 12: Spectral Methods in Polar Coordinates
Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
Last modified: February 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from laplacian_polar import laplacian_polar

# Output directory
fig_dir = os.path.join(os.path.dirname(__file__),
                       '..', '..', '..', 'textbook', 'figures', 'ch12', 'python')
os.makedirs(fig_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Build the Laplacian
# ---------------------------------------------------------------------------

Nr = 31   # radial points (odd)
M = 40    # angular points (even, more than eigenmodes for smoother plots)
N2 = (Nr - 1) // 2

L, r_pos, theta = laplacian_polar(Nr, M)

# Grid in (r, θ) space
RR, TT = np.meshgrid(r_pos, theta, indexing='ij')
XX = RR * np.cos(TT)
YY = RR * np.sin(TT)

# ---------------------------------------------------------------------------
# Poisson solve with off-centre Gaussian forcing
# ---------------------------------------------------------------------------

# Source: localized load at (x0, y0) = (0.4, 0.2)
x0, y0 = 0.4, 0.2
beta = 40.0
f_vals = np.exp(-beta * ((XX - x0)**2 + (YY - y0)**2))

# Flatten and solve: L u = f (note: L is the Laplacian, so -L u = f → u = -L \ f)
# Since L u = Δu, and we want -Δu = f, we need L u = -f → u = L \ (-f)
f_vec = -f_vals.flatten()
u_vec = np.linalg.solve(L, f_vec)
U = u_vec.reshape(N2, M)

# Extend for plotting: add boundary (r=1, u=0) and wrap in θ
U_ext = np.zeros((N2 + 1, M + 1))
U_ext[1:, :M] = U
U_ext[:, M] = U_ext[:, 0]

r_ext = np.concatenate(([1.0], r_pos))
theta_ext = np.concatenate((theta, [2 * np.pi]))
RR_ext, TT_ext = np.meshgrid(r_ext, theta_ext, indexing='ij')
XX_ext = RR_ext * np.cos(TT_ext)
YY_ext = RR_ext * np.sin(TT_ext)

# ---------------------------------------------------------------------------
# Figure 12.6: 3D surface plot of the Poisson solution
# ---------------------------------------------------------------------------

fig = plt.figure(figsize=(7, 5.5))
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(XX_ext, YY_ext, U_ext, cmap=cm.viridis, linewidth=0.2,
                edgecolor='k', alpha=0.9, rstride=1, cstride=1)

# Draw disk boundary circle at z=0
circle_th = np.linspace(0, 2 * np.pi, 200)
ax.plot(np.cos(circle_th), np.sin(circle_th),
        np.zeros(200), 'k-', linewidth=1.0)

ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$u$')
ax.view_init(elev=25, azim=-55)
ax.set_title('Poisson equation on the disk: off-centre Gaussian load',
             fontsize=10, pad=12)

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'poisson_disk_surface.pdf'),
            bbox_inches='tight')
plt.close()

# ---------------------------------------------------------------------------
# Figure 12.7: Filled contour plot on the disk
# ---------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(5.5, 5))

# Filled contour
cf = ax.contourf(XX_ext, YY_ext, U_ext, levels=20, cmap='viridis')
ax.contour(XX_ext, YY_ext, U_ext, levels=20, colors='k',
           linewidths=0.3, alpha=0.5)

# Disk boundary
ax.plot(np.cos(circle_th), np.sin(circle_th), 'k-', linewidth=1.2)

# Mark source location
ax.plot(x0, y0, 'r*', markersize=10, markeredgecolor='k', markeredgewidth=0.5)
ax.annotate('source', (x0, y0), textcoords="offset points",
            xytext=(8, 8), fontsize=9, color='r')

ax.set_xlim(-1.15, 1.15)
ax.set_ylim(-1.15, 1.15)
ax.set_aspect('equal')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_title('Poisson solution: filled contour', fontsize=10)
plt.colorbar(cf, ax=ax, shrink=0.8, label='$u(x,y)$')

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'poisson_disk_contour.pdf'),
            bbox_inches='tight')
plt.close()

# ---------------------------------------------------------------------------
# Figure 12.8: Radial symmetry test
# ---------------------------------------------------------------------------

# Case 1: radially symmetric Gaussian centred at origin
f_sym = np.exp(-beta * (XX**2 + YY**2))
f_vec_sym = -f_sym.flatten()
u_vec_sym = np.linalg.solve(L, f_vec_sym)
U_sym = u_vec_sym.reshape(N2, M)

# Case 2: off-centre Gaussian (already computed as U)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2))

# Left: radially symmetric case — profiles should collapse
n_angles = 8
angle_indices = np.linspace(0, M-1, n_angles, dtype=int, endpoint=False)
for k in angle_indices:
    ax1.plot(r_pos, U_sym[:, k], 'C0-', linewidth=0.8, alpha=0.7)
ax1.set_xlabel('$r$')
ax1.set_ylabel('$u(r, \\theta_k)$')
ax1.set_title('Centred Gaussian: radial profiles\n(all angles coincide)',
              fontsize=10)
ax1.invert_xaxis()  # r decreasing from left (r=1) to right (r≈0)

# Right: off-centre case — profiles differ
for k in angle_indices:
    label = f'$\\theta = {theta[k]:.2f}$' if k in angle_indices[:3] else None
    ax2.plot(r_pos, U[:, k], linewidth=0.8, alpha=0.7, label=label)
ax2.set_xlabel('$r$')
ax2.set_ylabel('$u(r, \\theta_k)$')
ax2.set_title('Off-centre Gaussian: radial profiles\n(angles differ)',
              fontsize=10)
ax2.legend(fontsize=8, loc='upper right')
ax2.invert_xaxis()

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'radial_symmetry_test.pdf'),
            bbox_inches='tight')
plt.close()

print(f"Figures saved to: {fig_dir}")
print("  poisson_disk_surface.pdf")
print("  poisson_disk_contour.pdf")
print("  radial_symmetry_test.pdf")
print(f"\nMax deflection: {np.max(U):.6f} at off-centre load")
print(f"Max deflection: {np.max(U_sym):.6f} at centred load")
