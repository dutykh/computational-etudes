#!/usr/bin/env python3
"""
wave1d_cheb.py
Chapter 10: Spectral PDE Solvers with Chebyshev and Fourier Grids

1D wave equation solver on Chebyshev grid with leapfrog time stepping.

PDE: u_tt = c^2 u_xx, -1 < x < 1, t > 0
BCs: u(-1, t) = u(1, t) = 0 (Dirichlet)
ICs: u(x, 0) = sin(pi/2 * (1 + x)), u_t(x, 0) = 0

Generates Figure 10.3: Waterfall plot of wave evolution.

Author: Dr. Denys Dutykh
Date: February 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D
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
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch10' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'wave1d_waterfall.pdf'


def wave1d_cheb(N=80, c=1.0, tmax=6.0, alpha=4.0, n_snapshots=100):
    """
    Solve 1D wave equation on Chebyshev grid using leapfrog time stepping.

    Parameters:
        N : int
            Number of Chebyshev intervals
        c : float
            Wave speed
        tmax : float
            Final time
        alpha : float
            CFL parameter (dt = alpha / N^2)
        n_snapshots : int
            Number of time snapshots to store

    Returns:
        x : array of shape (N+1,)
            Chebyshev grid points
        t_save : array of shape (n_snapshots,)
            Time values for snapshots
        U_save : array of shape (n_snapshots, N+1)
            Solution snapshots
    """
    # Chebyshev grid
    x = np.cos(np.pi * np.arange(N + 1) / N)

    # Initial data: half-sine
    u = np.sin(0.5 * np.pi * (1 + x))

    # Time stepping
    dt = alpha / N**2
    nsteps = int(np.ceil(tmax / dt))

    # Store snapshots
    save_every = max(1, nsteps // n_snapshots)
    U_save = [u.copy()]
    t_save = [0.0]

    # Taylor start for leapfrog (zero initial velocity)
    u_xx = chebfft(chebfft(u))
    u_xx[0] = u_xx[N] = 0  # Enforce BCs
    u_prev = u - 0.5 * (c * dt)**2 * u_xx

    t = 0.0
    for n in range(nsteps):
        # Second derivative via FFT
        u_xx = chebfft(chebfft(u))
        u_xx[0] = u_xx[N] = 0  # Enforce BCs

        # Leapfrog step
        u_new = 2 * u - u_prev + (c * dt)**2 * u_xx

        # Boundary conditions
        u_new[0] = 0
        u_new[N] = 0

        u_prev = u
        u = u_new
        t += dt

        # Save snapshot
        if (n + 1) % save_every == 0:
            U_save.append(u.copy())
            t_save.append(t)

    return x, np.array(t_save), np.array(U_save)


def main():
    # Solve wave equation
    N = 80
    tmax = 6.0
    x, t_save, U_save = wave1d_cheb(N=N, tmax=tmax, n_snapshots=80)

    # Create waterfall plot
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Create meshgrid
    X, T = np.meshgrid(x, t_save)

    # Plot surface
    surf = ax.plot_surface(X, T, U_save, cmap='coolwarm',
                           linewidth=0, antialiased=True,
                           alpha=0.9, rcount=80, ccount=80)

    # Add wireframe for better visualization
    ax.plot_wireframe(X, T, U_save, color='k', linewidth=0.2,
                      rstride=4, cstride=4, alpha=0.3)

    # Labels
    ax.set_xlabel(r'$x$', fontsize=12)
    ax.set_ylabel(r'$t$', fontsize=12)
    ax.set_zlabel(r'$u(x, t)$', fontsize=12)
    ax.set_title('1D Wave Equation: Vibrating String', fontsize=12, pad=10)

    # Adjust view angle
    ax.view_init(elev=25, azim=-60)

    # Add colorbar
    cbar = fig.colorbar(surf, ax=ax, shrink=0.6, aspect=15, pad=0.1)
    cbar.set_label(r'$u(x, t)$', fontsize=11)

    plt.tight_layout()

    # Save figure
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_FILE, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(OUTPUT_FILE.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {OUTPUT_FILE.resolve()}')

    # Also create a 2D waterfall-style plot
    fig2, ax2 = plt.subplots(figsize=(10, 6))

    # Select subset of times for waterfall
    n_lines = 40
    indices = np.linspace(0, len(t_save) - 1, n_lines, dtype=int)

    for i, idx in enumerate(indices):
        offset = t_save[idx] * 0.15
        ax2.plot(x, U_save[idx] + offset, color=plt.cm.viridis(i / n_lines),
                 linewidth=0.8)
        ax2.fill_between(x, offset, U_save[idx] + offset,
                         color=plt.cm.viridis(i / n_lines), alpha=0.1)

    ax2.set_xlabel(r'$x$', fontsize=11)
    ax2.set_ylabel(r'$u(x, t)$ (offset by $t$)', fontsize=11)
    ax2.set_title('1D Wave Equation: Waterfall View', fontsize=12)

    # Clean styling
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    plt.tight_layout()

    OUTPUT_FILE2 = OUTPUT_DIR / 'wave1d_waterfall_2d.pdf'
    fig2.savefig(OUTPUT_FILE2, bbox_inches='tight', pad_inches=0.05)

    print(f'Additional figure saved to: {OUTPUT_FILE2.resolve()}')

    plt.close('all')


if __name__ == '__main__':
    main()
