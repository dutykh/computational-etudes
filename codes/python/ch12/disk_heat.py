#!/usr/bin/env python3
"""
disk_heat.py - Heat equation on the unit disk

Computational Étude 12.5: Solves the heat equation u_t = α Δu + S(r,θ)
on the unit disk with Crank–Nicolson time stepping. Models heat dissipation
on a semiconductor wafer with a localised laser source.

Chapter 12: Spectral Methods in Polar Coordinates
Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
Last modified: February 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.linalg import lu_factor, lu_solve
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from laplacian_polar import laplacian_polar

# Output directory
fig_dir = os.path.join(os.path.dirname(__file__),
                       '..', '..', '..', 'textbook', 'figures', 'ch12', 'python')
os.makedirs(fig_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

Nr = 25       # radial points (odd)
M = 20        # angular points (even)
alpha = 0.1   # thermal diffusivity
dt = 0.01     # time step
t_final = 20.0
n_steps = int(t_final / dt)

N2 = (Nr - 1) // 2
ndof = N2 * M

# Build the Laplacian
L, r_pos, theta = laplacian_polar(Nr, M)

# Grid
RR, TT = np.meshgrid(r_pos, theta, indexing='ij')
XX = RR * np.cos(TT)
YY = RR * np.sin(TT)

# Source: localised laser spot at (0.3, 0.0)
x_src, y_src = 0.3, 0.0
S_vals = 2.0 * np.exp(-30.0 * ((XX - x_src)**2 + (YY - y_src)**2))
S_vec = S_vals.flatten()

# ---------------------------------------------------------------------------
# Crank–Nicolson time stepping
# ---------------------------------------------------------------------------

# (I - α dt/2 L) u^{n+1} = (I + α dt/2 L) u^n + dt S
I_mat = np.eye(ndof)
A = I_mat - alpha * dt / 2 * L     # LHS matrix
B = I_mat + alpha * dt / 2 * L     # RHS matrix part

# LU factorise once (reuse at every step)
lu, piv = lu_factor(A)

# Initial condition: u = 0 (cold start)
u = np.zeros(ndof)

# Store snapshots at selected times
snapshot_times = [0.0, 0.5, 2.0, 5.0, 10.0, 20.0]
snapshots = {}
max_temp = []
times = []

snapshots[0.0] = u.copy()
max_temp.append(0.0)
times.append(0.0)

t = 0.0
for step in range(1, n_steps + 1):
    rhs = B @ u + dt * S_vec
    u = lu_solve((lu, piv), rhs)
    t = step * dt

    max_temp.append(np.max(u))
    times.append(t)

    # Check if this is a snapshot time (within tolerance)
    for ts in snapshot_times:
        if abs(t - ts) < dt / 2 and ts not in snapshots:
            snapshots[ts] = u.copy()

# Ensure final snapshot is stored
if t_final not in snapshots:
    snapshots[t_final] = u.copy()

# ---------------------------------------------------------------------------
# Figure 12.9: Heat equation snapshots at 6 time levels
# ---------------------------------------------------------------------------

# Extend grid for plotting
r_ext = np.concatenate(([1.0], r_pos))
theta_ext = np.concatenate((theta, [2 * np.pi]))
RR_ext, TT_ext = np.meshgrid(r_ext, theta_ext, indexing='ij')
XX_ext = RR_ext * np.cos(TT_ext)
YY_ext = RR_ext * np.sin(TT_ext)

# Global color scale
u_max = max(np.max(np.abs(s)) for s in snapshots.values())
if u_max < 1e-10:
    u_max = 1.0

fig, axes = plt.subplots(2, 3, figsize=(11, 7.5))
circle_th = np.linspace(0, 2 * np.pi, 200)

for idx, (ts, ax) in enumerate(zip(snapshot_times, axes.flatten())):
    U = snapshots[ts].reshape(N2, M)

    # Extend for plotting
    U_ext = np.zeros((N2 + 1, M + 1))
    U_ext[1:, :M] = U
    U_ext[:, M] = U_ext[:, 0]

    cf = ax.contourf(XX_ext, YY_ext, U_ext, levels=20,
                     cmap='hot', vmin=0, vmax=u_max)
    ax.plot(np.cos(circle_th), np.sin(circle_th), 'k-', linewidth=1.0)
    ax.set_xlim(-1.15, 1.15)
    ax.set_ylim(-1.15, 1.15)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(f'$t = {ts:.1f}$', fontsize=10)

# Single colorbar
fig.subplots_adjust(right=0.88)
cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
fig.colorbar(cf, cax=cbar_ax, label='Temperature $u$')

plt.suptitle('Heat dissipation on a disk with localised source', fontsize=11,
             y=0.98)
plt.savefig(os.path.join(fig_dir, 'heat_disk_snapshots.pdf'),
            bbox_inches='tight')
plt.close()

# ---------------------------------------------------------------------------
# Figure 12.10: Temperature evolution and steady-state comparison
# ---------------------------------------------------------------------------

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2))

# Left: maximum temperature vs time
ax1.plot(times, max_temp, 'C0-', linewidth=1.0)
ax1.set_xlabel('Time $t$')
ax1.set_ylabel('$\\max_{{x,y}} \\, u(x,y,t)$')
ax1.set_title('Maximum temperature evolution', fontsize=10)
ax1.grid(True, alpha=0.3)

# Right: compare steady state with direct Poisson solve
# Steady state satisfies α Δu + S = 0, i.e., Δu = -S/α, i.e., u = L \ (-S_vec/α)
u_poisson = np.linalg.solve(L, -S_vec / alpha)

# Plot radial profiles at θ = 0
angle_idx = 0
U_final = snapshots[t_final].reshape(N2, M)
U_poisson = u_poisson.reshape(N2, M)

ax2.plot(r_pos, U_final[:, angle_idx], 'C0-', linewidth=1.5,
         label='CN steady state')
ax2.plot(r_pos, U_poisson[:, angle_idx], 'C3--', linewidth=1.5,
         label='Direct Poisson solve')
ax2.set_xlabel('$r$')
ax2.set_ylabel('$u(r, 0)$')
ax2.set_title('Steady state vs Poisson solution ($\\theta = 0$)', fontsize=10)
ax2.legend(fontsize=9)
ax2.invert_xaxis()

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'heat_disk_energy.pdf'),
            bbox_inches='tight')
plt.close()

# Print summary
err_ss = np.max(np.abs(snapshots[t_final] - u_poisson))
print(f"Steady-state vs Poisson max error: {err_ss:.2e}")
print(f"Max temperature at t={t_final}: {max_temp[-1]:.6f}")
print(f"\nFigures saved to: {fig_dir}")
print("  heat_disk_snapshots.pdf")
print("  heat_disk_energy.pdf")
