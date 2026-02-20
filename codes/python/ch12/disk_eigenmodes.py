#!/usr/bin/env python3
"""
disk_eigenmodes.py - Eigenmodes of the Laplacian on the unit disk

Computational Étude 12.2: Computes eigenmodes of -Δu = λu on the disk
with u = 0 on the boundary. Compares eigenvalues with Bessel function
zeros to demonstrate spectral accuracy.

Chapter 12: Spectral Methods in Polar Coordinates
Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
Last modified: February 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.special import jn_zeros
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from laplacian_polar import laplacian_polar

# Output directory
fig_dir = os.path.join(os.path.dirname(__file__),
                       '..', '..', '..', 'textbook', 'figures', 'ch12', 'python')
os.makedirs(fig_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Build the Laplacian and compute eigenmodes
# ---------------------------------------------------------------------------

Nr = 25  # radial points (odd)
M = 20   # angular points (even)
N2 = (Nr - 1) // 2

L, r_pos, theta = laplacian_polar(Nr, M)

# Solve the eigenvalue problem -L v = λ v
eigenvalues, eigenvectors = np.linalg.eig(-L)
eigenvalues = np.real(eigenvalues)

# Sort by eigenvalue
idx = np.argsort(eigenvalues)
eigenvalues = eigenvalues[idx]
eigenvectors = eigenvectors[:, idx]

# ---------------------------------------------------------------------------
# Figure 12.3: Six representative eigenmodes as surface plots
# ---------------------------------------------------------------------------

# We show modes 1, 2, 6, 9 (matching Trefethen's selection) plus 3 and 5
# for a richer display. These correspond to (m,n) = (0,1), (1,1), (2,1),
# (0,2), (3,1), (1,2).
mode_indices = [0, 1, 3, 5, 6, 8]  # 0-based indices

fig = plt.figure(figsize=(12, 8))

for panel, mi in enumerate(mode_indices):
    ax = fig.add_subplot(2, 3, panel + 1, projection='3d')

    # Reshape eigenvector: the solution vector is ordered as
    # [u(r_1,θ_1), u(r_1,θ_2), ..., u(r_1,θ_M), u(r_2,θ_1), ...]
    v = np.real(eigenvectors[:, mi])
    U = v.reshape(N2, M)

    # Normalise
    U = U / np.max(np.abs(U))

    # Extend to include boundary (r=1, u=0) and center region
    # Add a row of zeros at r=1 (top) and extend periodically in θ
    U_ext = np.zeros((N2 + 1, M + 1))
    U_ext[1:, :M] = U
    U_ext[:, M] = U_ext[:, 0]  # periodic wrap

    r_ext = np.concatenate(([1.0], r_pos))
    theta_ext = np.concatenate((theta, [2 * np.pi]))

    RR, TT = np.meshgrid(r_ext, theta_ext, indexing='ij')
    XX = RR * np.cos(TT)
    YY = RR * np.sin(TT)

    ax.plot_surface(XX, YY, U_ext, cmap=cm.RdBu_r, linewidth=0.3,
                    edgecolor='k', alpha=0.85, rstride=1, cstride=1,
                    antialiased=True)
    ax.set_zlim(-1.1, 1.1)
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 1.1)
    ax.view_init(elev=25, azim=-60)
    ax.set_axis_off()

    lam_scaled = np.sqrt(eigenvalues[mi] / eigenvalues[0])
    ax.set_title(f'Mode {mi+1}\n$\\lambda = {lam_scaled:.10f}$', fontsize=9)

plt.suptitle('Eigenmodes of the Laplacian on the unit disk', fontsize=12, y=0.98)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'disk_eigenmodes.pdf'), bbox_inches='tight')
plt.close()

# ---------------------------------------------------------------------------
# Figure 12.4: Eigenvalue convergence vs Bessel zeros
# ---------------------------------------------------------------------------

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2))

# --- Left: computed vs exact eigenvalues ---
n_show = 25
lam_computed = eigenvalues[:n_show]

# Collect exact Bessel eigenvalues and sort
bessel_lam = []
for m in range(10):
    z = jn_zeros(m, 5)
    mult = 1 if m == 0 else 2
    for zval in z:
        bessel_lam.extend([zval**2] * mult)
bessel_lam = np.sort(bessel_lam)[:n_show]

ax1.plot(range(1, n_show+1), lam_computed, 'o', markersize=6,
         markerfacecolor='none', markeredgecolor='C0', label='Computed')
ax1.plot(range(1, n_show+1), bessel_lam, 'x', markersize=5,
         color='C3', label='Bessel zeros $j_{m,n}^2$')
ax1.set_xlabel('Mode index')
ax1.set_ylabel('Eigenvalue $\\lambda$')
ax1.legend(fontsize=9)
ax1.set_title('Eigenvalue comparison', fontsize=10)

# --- Right: convergence of first eigenvalue with Nr ---
Nr_values = [9, 13, 17, 21, 25, 29, 33]
M_fixed = 20
lam_exact_1 = jn_zeros(0, 1)[0]**2  # j_{0,1}^2

errors = []
for N in Nr_values:
    L_test, _, _ = laplacian_polar(N, M_fixed)
    eigs = np.sort(np.real(np.linalg.eigvals(-L_test)))
    errors.append(abs(eigs[0] - lam_exact_1))

ax2.semilogy(Nr_values, errors, 'o-', markersize=5, color='C0')
ax2.set_xlabel('$N_r$ (radial points)')
ax2.set_ylabel('$|\\lambda_1 - j_{0,1}^2|$')
ax2.set_title('Exponential convergence', fontsize=10)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'eigenvalue_convergence.pdf'),
            bbox_inches='tight')
plt.close()

# Print eigenvalue comparison table
print(f"{'Mode':>5s}  {'Computed':>16s}  {'Exact (Bessel)':>16s}  {'Error':>12s}")
print("-" * 58)
for i in range(min(15, n_show)):
    err = abs(lam_computed[i] - bessel_lam[i])
    print(f"{i+1:5d}  {lam_computed[i]:16.10f}  {bessel_lam[i]:16.10f}  {err:12.2e}")

print(f"\nFigures saved to: {fig_dir}")
print("  disk_eigenmodes.pdf")
print("  eigenvalue_convergence.pdf")
