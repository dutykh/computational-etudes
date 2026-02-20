#!/usr/bin/env python3
"""
disk_nodal_rotation.py - Eigenvalue degeneracy and nodal rotations

Computational Étude 12.3: For a doubly-degenerate eigenvalue (m ≥ 1),
the two eigenvectors span a 2D eigenspace. Forming u_φ = cos(φ)·u₁ +
sin(φ)·u₂ produces continuously rotating nodal line patterns.

Chapter 12: Spectral Methods in Polar Coordinates
Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
Last modified: February 2026
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from laplacian_polar import laplacian_polar

# Output directory
fig_dir = os.path.join(os.path.dirname(__file__),
                       '..', '..', '..', 'textbook', 'figures', 'ch12', 'python')
os.makedirs(fig_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Build Laplacian and compute the first degenerate pair
# ---------------------------------------------------------------------------

Nr = 25
M = 20
N2 = (Nr - 1) // 2

L, r_pos, theta = laplacian_polar(Nr, M)

eigenvalues, eigenvectors = np.linalg.eig(-L)
eigenvalues = np.real(eigenvalues)
idx = np.argsort(eigenvalues)
eigenvalues = eigenvalues[idx]
eigenvectors = eigenvectors[:, idx]

# The first degenerate pair is modes 2 and 3 (eigenvalue j_{1,1}^2 ≈ 14.682)
u1 = np.real(eigenvectors[:, 1])
u2 = np.real(eigenvectors[:, 2])

# Build Cartesian grid for contour plotting (fine interpolation grid)
n_fine = 200
r_fine = np.linspace(0, 1, n_fine)
theta_fine = np.linspace(0, 2 * np.pi, n_fine)
RR_fine, TT_fine = np.meshgrid(r_fine, theta_fine)
XX_fine = RR_fine * np.cos(TT_fine)
YY_fine = RR_fine * np.sin(TT_fine)

# Build the coarse grid for nodal line plotting
# Extend to include boundary at r=1 and wrap in theta
U1 = u1.reshape(N2, M)
U2 = u2.reshape(N2, M)

# Add boundary (r=1: u=0) row at top
U1_ext = np.zeros((N2 + 1, M + 1))
U2_ext = np.zeros((N2 + 1, M + 1))
U1_ext[1:, :M] = U1
U2_ext[1:, :M] = U2
U1_ext[:, M] = U1_ext[:, 0]
U2_ext[:, M] = U2_ext[:, 0]

r_ext = np.concatenate(([1.0], r_pos))
theta_ext = np.concatenate((theta, [2 * np.pi]))
RR, TT = np.meshgrid(r_ext, theta_ext, indexing='ij')
XX = RR * np.cos(TT)
YY = RR * np.sin(TT)

# ---------------------------------------------------------------------------
# Figure 12.5: Nodal line rotation for 6 values of φ
# ---------------------------------------------------------------------------

phi_values = np.linspace(0, 5 * np.pi / 6, 6)  # 0, π/6, π/3, π/2, 2π/3, 5π/6

fig, axes = plt.subplots(2, 3, figsize=(10, 7))

for idx_p, (phi, ax) in enumerate(zip(phi_values, axes.flatten())):
    # Linear combination
    U_phi = np.cos(phi) * U1_ext + np.sin(phi) * U2_ext

    # Draw disk boundary
    circle_theta = np.linspace(0, 2 * np.pi, 200)
    ax.plot(np.cos(circle_theta), np.sin(circle_theta), 'k-', linewidth=1.0)

    # Plot nodal lines (zero contour)
    ax.contour(XX, YY, U_phi, levels=[0], colors='C0', linewidths=1.5)

    # Light filled contour for context
    ax.contourf(XX, YY, U_phi, levels=20, cmap='RdBu_r', alpha=0.3)

    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 1.1)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(f'$\\varphi = {idx_p}\\pi/6$' if idx_p > 0
                 else '$\\varphi = 0$', fontsize=10)

plt.suptitle('Nodal line rotation in a degenerate eigenspace '
             f'($\\lambda \\approx {eigenvalues[1]:.4f}$)',
             fontsize=11, y=0.99)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'nodal_rotation.pdf'), bbox_inches='tight')
plt.close()

print(f"Figure saved to: {fig_dir}/nodal_rotation.pdf")
print(f"Degenerate eigenvalue: λ = {eigenvalues[1]:.10f} (modes 2 and 3)")
