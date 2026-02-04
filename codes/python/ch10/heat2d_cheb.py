#!/usr/bin/env python3
"""
heat2d_cheb.py
Chapter 10: Spectral PDE Solvers with Chebyshev and Fourier Grids

2D heat equation solver on Chebyshev tensor grid with backward Euler.

PDE: u_t = u_xx + u_yy, -1 < x, y < 1, t > 0
BCs: u = 0 on boundary (Dirichlet)
ICs: u(x, y, 0) = exp(-25*((x+0.3)^2 + (y-0.1)^2))

Generates Figure 10.6: Heat diffusion evolution on 2D domain.

Author: Dr. Denys Dutykh
Date: February 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pathlib import Path

from chebfft import cheb_matrix

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
OUTPUT_FILE = OUTPUT_DIR / 'heat2d_snapshots.pdf'


def heat2d_cheb(N=24, tmax=0.5, dt=0.005):
    """
    Solve 2D heat equation on Chebyshev tensor grid using backward Euler.

    The discrete Laplacian uses Kronecker sum:
        L = D2_int ⊗ I + I ⊗ D2_int

    The system is:
        (I - dt * L) vec(U^{n+1}) = vec(U^n)

    Parameters:
        N : int
            Number of Chebyshev intervals in each direction
        tmax : float
            Final time
        dt : float
            Time step

    Returns:
        xx, yy : arrays of shape (N+1, N+1)
            Meshgrid of Chebyshev points
        snapshots : list of (time, U) tuples
            Solution snapshots
    """
    # Build differentiation matrices
    D, x = cheb_matrix(N)
    D2 = D @ D
    xx, yy = np.meshgrid(x, x)

    # Interior block
    D2_int = D2[1:N, 1:N]
    I_int = np.eye(N - 1)

    # Kronecker sum Laplacian for interior
    L = np.kron(D2_int, I_int) + np.kron(I_int, D2_int)
    A = np.eye((N - 1)**2) - dt * L

    # Initial condition
    U = np.exp(-25 * ((xx + 0.3)**2 + (yy - 0.1)**2))
    U[0, :] = U[-1, :] = U[:, 0] = U[:, -1] = 0

    # Times for snapshots
    t_snap = [0.0, 0.05, 0.15, 0.5]
    snapshots = [(0.0, U.copy())]
    current_snap_idx = 1

    # Time stepping
    nsteps = int(np.ceil(tmax / dt))
    t = 0.0

    for n in range(nsteps):
        # Extract interior and solve
        u_int = U[1:N, 1:N].flatten()
        u_int_new = np.linalg.solve(A, u_int)

        # Reconstruct
        U_new = np.zeros_like(U)
        U_new[1:N, 1:N] = u_int_new.reshape(N - 1, N - 1)
        U = U_new

        t += dt

        # Save snapshot
        if current_snap_idx < len(t_snap) and t >= t_snap[current_snap_idx]:
            snapshots.append((t, U.copy()))
            current_snap_idx += 1

    return xx, yy, snapshots


def main():
    # Solve heat equation
    N = 32
    xx, yy, snapshots = heat2d_cheb(N=N, tmax=0.6, dt=0.002)

    # Create figure with 4 panels
    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    axes = axes.flatten()

    # Plot each snapshot with individual color scaling
    for i, (t, U) in enumerate(snapshots[:4]):
        ax = axes[i]

        # Individual color scaling for each panel
        vmin = 0
        vmax = np.max(U)
        levels = np.linspace(vmin, vmax, 30)

        # Contourf plot
        cf = ax.contourf(xx, yy, U, levels=levels, cmap='hot')

        # Add contour lines
        ax.contour(xx, yy, U, levels=8, colors='k', linewidths=0.3, alpha=0.5)

        # Mark the initial center on first panel
        if i == 0:
            ax.plot(-0.3, 0.1, 'w+', markersize=10, markeredgewidth=2)

        ax.set_xlabel(r'$x$', fontsize=11)
        ax.set_ylabel(r'$y$', fontsize=11)
        ax.set_title(f'$t = {t:.3f}$', fontsize=12)
        ax.set_aspect('equal')

        # Add colorbar
        cbar = plt.colorbar(cf, ax=ax, shrink=0.8)
        if i % 2 == 1:
            cbar.set_label(r'$u(x, y, t)$', fontsize=10)

    plt.suptitle('2D Heat Equation: Diffusion on a Square', fontsize=14, y=1.02)
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
