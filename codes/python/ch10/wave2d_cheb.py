#!/usr/bin/env python3
"""
wave2d_cheb.py
Chapter 10: Spectral PDE Solvers with Chebyshev and Fourier Grids

2D wave equation solver on Chebyshev tensor grid with leapfrog time stepping.

PDE: u_tt = c^2 (u_xx + u_yy), -1 < x, y < 1, t > 0
BCs: u = 0 on boundary (Dirichlet)
ICs: u(x, y, 0) = exp(-30*((x-0.2)^2 + (y+0.3)^2)), u_t(x, y, 0) = 0

Generates Figure 10.4: Snapshots of wave evolution on 2D domain.

Author: Dr. Denys Dutykh
Date: February 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pathlib import Path

from chebfft import chebfft

# -----------------------------------------------------------------------------
# Publication-quality matplotlib configuration
# -----------------------------------------------------------------------------
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman', 'CMU Serif', 'DejaVu Serif']
rcParams['mathtext.fontset'] = 'cm'
rcParams['axes.labelsize'] = 11
rcParams['axes.titlesize'] = 12
rcParams['xtick.labelsize'] = 10
rcParams['ytick.labelsize'] = 10
rcParams['legend.fontsize'] = 10
rcParams['figure.dpi'] = 150
rcParams['savefig.dpi'] = 300
rcParams['text.usetex'] = False
rcParams['axes.linewidth'] = 0.8

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
NAVY = '#142D6E'

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch10' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'wave2d_snapshots.pdf'


def laplacian_chebfft(U, N):
    """Compute the 2D Laplacian using chebfft along rows and columns."""
    Lap = np.zeros_like(U)
    for i in range(N + 1):
        Lap[i, :] += chebfft(chebfft(U[i, :]))   # d^2/dx^2
    for j in range(N + 1):
        Lap[:, j] += chebfft(chebfft(U[:, j]))   # d^2/dy^2
    return Lap


def wave2d_cheb(N=32, c=1.0, tmax=2.0, alpha=3.0):
    """
    Solve 2D wave equation on Chebyshev tensor grid.

    Parameters:
        N : int
            Number of Chebyshev intervals in each direction
        c : float
            Wave speed
        tmax : float
            Final time
        alpha : float
            CFL parameter (dt = alpha / N^2)

    Returns:
        xx, yy : arrays of shape (N+1, N+1)
            Meshgrid of Chebyshev points
        snapshots : list of (time, U) tuples
            Solution snapshots at selected times
    """
    # Chebyshev grid
    x = np.cos(np.pi * np.arange(N + 1) / N)
    xx, yy = np.meshgrid(x, x)

    # Initial data: offset Gaussian
    U = np.exp(-30 * ((xx - 0.2)**2 + (yy + 0.3)**2))
    U_prev = U.copy()

    # Time stepping
    dt = alpha / N**2
    nsteps = int(np.ceil(tmax / dt))

    # Times for snapshots
    t_snap = [0.0, 0.3, 0.7, 1.2]
    snapshots = [(0.0, U.copy())]
    current_snap_idx = 1

    # Taylor start
    Lap_U = laplacian_chebfft(U, N)
    U_prev = U - 0.5 * (c * dt)**2 * Lap_U

    t = 0.0
    for n in range(nsteps):
        # 2D Laplacian via chebfft along rows and columns
        Lap_U = laplacian_chebfft(U, N)

        # Leapfrog step
        U_new = 2 * U - U_prev + (c * dt)**2 * Lap_U

        # Boundary conditions
        U_new[0, :] = 0
        U_new[-1, :] = 0
        U_new[:, 0] = 0
        U_new[:, -1] = 0

        U_prev = U
        U = U_new
        t += dt

        # Save snapshot if we're at a target time
        if current_snap_idx < len(t_snap) and t >= t_snap[current_snap_idx]:
            snapshots.append((t, U.copy()))
            current_snap_idx += 1

    return xx, yy, snapshots


def main():
    # Solve wave equation
    N = 40
    xx, yy, snapshots = wave2d_cheb(N=N, tmax=1.5)

    # Create figure with 4 panels
    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    axes = axes.flatten()

    # Plot each snapshot
    vmin, vmax = -0.5, 1.0
    for i, (t, U) in enumerate(snapshots[:4]):
        ax = axes[i]

        # Contourf plot
        levels = np.linspace(vmin, vmax, 30)
        cf = ax.contourf(xx, yy, U, levels=levels, cmap='RdBu_r', extend='both')

        # Add contour lines
        ax.contour(xx, yy, U, levels=10, colors='k', linewidths=0.3, alpha=0.5)

        # Mark the initial center
        if i == 0:
            ax.plot(0.2, -0.3, 'k+', markersize=10, markeredgewidth=2)

        ax.set_xlabel(r'$x$', fontsize=11)
        ax.set_ylabel(r'$y$', fontsize=11)
        ax.set_title(f'$t = {t:.2f}$', fontsize=12)
        ax.set_aspect('equal')

        # Add colorbar to each panel
        cbar = plt.colorbar(cf, ax=ax, shrink=0.8)
        if i % 2 == 1:
            cbar.set_label(r'$u(x, y, t)$', fontsize=10)

    plt.suptitle('2D Wave Equation: Vibrating Membrane', fontsize=14, y=1.02)
    plt.tight_layout()

    # Save figure
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_FILE, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(OUTPUT_FILE.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {OUTPUT_FILE.resolve()}')

    plt.close(fig)


if __name__ == '__main__':
    main()
