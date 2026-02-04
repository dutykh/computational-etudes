#!/usr/bin/env python3
"""
transport_variable.py
Chapter 10: Spectral PDE Solvers with Chebyshev and Fourier Grids

Variable coefficient transport equation with periodic BCs using Fourier.

PDE: u_t + c(x) * u_x = 0, 0 < x < 2*pi, t > 0 (periodic)
Speed: c(x) = 0.3 + sin^2(x - 1)
ICs: u(x, 0) = exp(-20*(x - pi)^2)

Uses FFT for spatial differentiation and leapfrog for time stepping.

Generates Figure 10.8: Variable speed transport waterfall plot.

Author: Dr. Denys Dutykh
Date: February 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path

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
OUTPUT_FILE = OUTPUT_DIR / 'transport_variable.pdf'


def fourier_diff(v):
    """
    Compute the derivative of a periodic function using FFT.

    Parameters:
        v : array of length N
            Function values at equispaced points on [0, 2*pi)

    Returns:
        w : array of length N
            Approximate derivative values
    """
    N = len(v)
    v_hat = np.fft.fft(v)

    # Wavenumbers
    if N % 2 == 0:
        k = np.concatenate([np.arange(N // 2), [0], np.arange(1 - N // 2, 0)])
    else:
        k = np.concatenate([np.arange((N - 1) // 2 + 1),
                           np.arange(-(N - 1) // 2, 0)])

    w = np.fft.ifft(1j * k * v_hat).real
    return w


def transport_variable(N=256, tmax=12.0, n_snapshots=100):
    """
    Solve variable coefficient transport equation using Fourier spectral.

    Uses RK4 time stepping for stability.

    Parameters:
        N : int
            Number of grid points
        tmax : float
            Final time
        n_snapshots : int
            Number of snapshots to save

    Returns:
        x : array of shape (N,)
            Grid points on [0, 2*pi)
        t_save : array
            Times of snapshots
        U_save : array
            Solution snapshots
    """
    # Grid
    h = 2 * np.pi / N
    x = h * np.arange(N)

    # Variable speed: c(x) = 0.3 + sin^2(x - 1)
    c = 0.3 + np.sin(x - 1)**2

    # Initial condition: Gaussian pulse
    u = np.exp(-20 * (x - np.pi)**2)

    # CFL condition for variable coefficient
    c_max = np.max(c)
    dt = 0.4 * h / c_max
    nsteps = int(np.ceil(tmax / dt))

    # Store snapshots
    save_every = max(1, nsteps // n_snapshots)
    U_save = [u.copy()]
    t_save = [0.0]

    # RHS function for RK4
    def rhs(u):
        return -c * fourier_diff(u)

    t = 0.0
    for n in range(nsteps):
        # RK4 time stepping
        k1 = rhs(u)
        k2 = rhs(u + 0.5 * dt * k1)
        k3 = rhs(u + 0.5 * dt * k2)
        k4 = rhs(u + dt * k3)
        u = u + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)

        t += dt

        # Save snapshot
        if (n + 1) % save_every == 0:
            U_save.append(u.copy())
            t_save.append(t)

    return x, np.array(t_save), np.array(U_save)


def main():
    # Solve transport equation
    N = 256
    tmax = 8.0
    x, t_save, U_save = transport_variable(N=N, tmax=tmax, n_snapshots=120)

    # Create figure with two panels
    fig, (ax_main, ax_water) = plt.subplots(1, 2, figsize=(12, 5))

    # Left panel: space-time contour plot
    X, T = np.meshgrid(x, t_save)
    pcm = ax_main.pcolormesh(X, T, U_save, cmap='coolwarm', shading='gouraud')
    ax_main.set_xlabel(r'$x$', fontsize=11)
    ax_main.set_ylabel(r'$t$', fontsize=11)
    ax_main.set_title('Variable Coefficient Transport', fontsize=12)
    fig.colorbar(pcm, ax=ax_main, shrink=0.8, label=r'$u(x, t)$')

    # Right panel: waterfall view
    n_lines = 40
    indices = np.linspace(0, len(t_save) - 1, n_lines, dtype=int)

    for i, idx in enumerate(indices):
        offset = t_save[idx] * 0.05
        color = plt.cm.viridis(i / n_lines)
        ax_water.plot(x, U_save[idx] + offset, color=color, linewidth=0.8)

    ax_water.set_xlabel(r'$x$', fontsize=11)
    ax_water.set_ylabel(r'$u(x, t)$ (offset by $t$)', fontsize=11)
    ax_water.set_title('Waterfall View: Pulse in Variable Speed', fontsize=12)

    # Add inset showing the speed profile
    ax_inset = ax_water.inset_axes([0.6, 0.7, 0.35, 0.25])
    c = 0.3 + np.sin(x - 1)**2
    ax_inset.plot(x, c, color=CORAL, linewidth=1.5)
    ax_inset.set_xlabel(r'$x$', fontsize=8)
    ax_inset.set_ylabel(r'$c(x)$', fontsize=8)
    ax_inset.set_title('Speed Profile', fontsize=9)
    ax_inset.tick_params(labelsize=7)

    # Clean styling
    ax_water.spines['top'].set_visible(False)
    ax_water.spines['right'].set_visible(False)

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
