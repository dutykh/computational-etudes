#!/usr/bin/env python3
"""
heat1d_cheb.py
Chapter 10: Spectral PDE Solvers with Chebyshev and Fourier Grids

1D heat equation solver on Chebyshev grid with Crank--Nicolson time stepping.

PDE: u_t = kappa * u_xx, -1 < x < 1, t > 0
BCs: u(-1, t) = u(1, t) = 0 (Dirichlet)
ICs: u(x, 0) = exp(-20*(x+0.5)^2) + 0.5*exp(-30*(x-0.4)^2) (two bumps)

Generates Figure 10.5: Heat diffusion evolution.

Author: Dr. Denys Dutykh
Date: February 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D
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
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'
PURPLE = '#9B59B6'

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch10' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'heat1d_evolution.pdf'


def heat1d_cheb(N=64, kappa=1.0, tmax=1.0, dt=0.01, n_snapshots=100):
    """
    Solve 1D heat equation using Crank--Nicolson time stepping.

    The semi-discrete system is:
        (I - kappa*dt/2 * D2_int) u^{n+1} = (I + kappa*dt/2 * D2_int) u^n

    Parameters:
        N : int
            Number of Chebyshev intervals
        kappa : float
            Diffusion coefficient
        tmax : float
            Final time
        dt : float
            Time step
        n_snapshots : int
            Number of snapshots to save

    Returns:
        x : array of shape (N+1,)
            Chebyshev grid
        t_save : array
            Times of snapshots
        U_save : array
            Solution snapshots
    """
    # Build differentiation matrices
    D, x = cheb_matrix(N)
    D2 = D @ D

    # Interior indices (Dirichlet BCs: u[0] = u[N] = 0)
    ii = slice(1, N)
    D2_int = D2[np.ix_(range(1, N), range(1, N))]
    I = np.eye(N - 1)

    # Crank--Nicolson matrices
    A = I - 0.5 * kappa * dt * D2_int
    B = I + 0.5 * kappa * dt * D2_int

    # Initial condition: two bumps
    u_full = np.exp(-20 * (x + 0.5)**2) + 0.5 * np.exp(-30 * (x - 0.4)**2)
    u = u_full[ii]

    # Time stepping
    nsteps = int(np.ceil(tmax / dt))
    save_every = max(1, nsteps // n_snapshots)

    U_save = [u_full.copy()]
    t_save = [0.0]

    t = 0.0
    for n in range(nsteps):
        # Solve the CN system
        u = np.linalg.solve(A, B @ u)
        t += dt

        if (n + 1) % save_every == 0:
            # Reconstruct full solution
            u_full = np.zeros(N + 1)
            u_full[ii] = u
            U_save.append(u_full.copy())
            t_save.append(t)

    return x, np.array(t_save), np.array(U_save)


def main():
    # Solve heat equation
    N = 128
    tmax = 0.5
    x, t_save, U_save = heat1d_cheb(N=N, tmax=tmax, dt=0.005, n_snapshots=80)

    # Create figure with main plot and insets
    fig = plt.figure(figsize=(12, 6))

    # Main plot: 3D surface
    ax_main = fig.add_subplot(1, 2, 1, projection='3d')

    X, T = np.meshgrid(x, t_save)
    surf = ax_main.plot_surface(X, T, U_save, cmap='hot',
                                linewidth=0, antialiased=True,
                                alpha=0.9)

    ax_main.set_xlabel(r'$x$', fontsize=11)
    ax_main.set_ylabel(r'$t$', fontsize=11)
    ax_main.set_zlabel(r'$u(x, t)$', fontsize=11)
    ax_main.set_title('Heat Diffusion: Space-Time View', fontsize=12)
    ax_main.view_init(elev=25, azim=-60)

    # Right side: snapshots at selected times
    ax_snap = fig.add_subplot(1, 2, 2)

    # Select specific times for snapshots
    t_targets = [0.0, 0.05, 0.1, 0.2, 0.5]
    colors = [NAVY, SKY, TEAL, CORAL, PURPLE]

    for t_target, color in zip(t_targets, colors):
        idx = np.argmin(np.abs(t_save - t_target))
        ax_snap.plot(x, U_save[idx], '-', color=color, linewidth=1.5,
                     label=f'$t = {t_save[idx]:.2f}$')

    ax_snap.set_xlabel(r'$x$', fontsize=11)
    ax_snap.set_ylabel(r'$u(x, t)$', fontsize=11)
    ax_snap.set_title('Heat Diffusion: Time Snapshots', fontsize=12)
    ax_snap.legend(loc='upper right', fontsize=10)

    # Clean styling
    ax_snap.spines['top'].set_visible(False)
    ax_snap.spines['right'].set_visible(False)
    ax_snap.set_xlim(-1, 1)
    ax_snap.set_ylim(0, 1.1)

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
