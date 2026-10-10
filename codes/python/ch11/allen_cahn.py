#!/usr/bin/env python3
"""
allen_cahn.py
Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs

Allen-Cahn equation with periodic BCs using Fourier spectral IMEX scheme.

PDE: u_t = eps^2 * u_xx + u(1 - u^2), 0 < x < 2*pi, t > 0 (periodic)

The parameter eps controls the interface width between phases u = +1 and u = -1.
The diffusion term eps^2 * u_xx is treated implicitly and the nonlinear reaction
u(1 - u^2) is treated explicitly (IMEX):
    (I - eps^2 * dt * D_xx) u^{n+1} = u^n + dt * (u^n - (u^n)^3)

In Fourier space this becomes diagonal:
    u_hat_k^{n+1} = (u_hat_k^n + dt * FFT(u^n - (u^n)^3)) / (1 + eps^2 * k^2 * dt)

Generates Figure 11.x: Allen-Cahn phase separation via IMEX spectral method.

Author: Dr. Denys Dutykh
Date: February 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
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
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch11' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'allen_cahn.pdf'


def allen_cahn_imex(N=512, eps=0.05, tmax=100.0, dt=0.01, n_snapshots=400):
    """
    Solve the Allen-Cahn equation using Fourier spectral IMEX scheme.

    Diffusion (eps^2 * u_xx) is treated implicitly in Fourier space.
    Reaction (u - u^3) is treated explicitly.

    Parameters:
        N : int
            Number of grid points
        eps : float
            Interface width parameter
        tmax : float
            Final time
        dt : float
            Time step size
        n_snapshots : int
            Number of snapshots to save

    Returns:
        x : array of shape (N,)
            Grid points on [0, 2*pi)
        t_save : array
            Times of snapshots
        U_save : array
            Solution snapshots (n_snapshots x N)
    """
    # Grid
    h = 2 * np.pi / N
    x = h * np.arange(N)

    # Wavenumbers for FFT
    if N % 2 == 0:
        k = np.concatenate([np.arange(N // 2), [0], np.arange(1 - N // 2, 0)])
    else:
        k = np.concatenate([np.arange((N - 1) // 2 + 1),
                           np.arange(-(N - 1) // 2, 0)])

    # Initial condition: smooth random perturbation (fixed seed)
    np.random.seed(42)
    u = 0.1 * np.random.randn(N)
    # Low-pass filter to smooth the initial data
    u_hat = np.fft.fft(u)
    filter_mask = np.abs(k) <= 10
    u_hat *= filter_mask
    u = np.fft.ifft(u_hat).real

    # IMEX multiplier: 1 / (1 + eps^2 * k^2 * dt)
    # This is the implicit part for the diffusion operator
    imex_denom = 1.0 + eps**2 * k**2 * dt

    # Time stepping
    nsteps = int(np.ceil(tmax / dt))

    # Store snapshots
    save_every = max(1, nsteps // n_snapshots)
    U_save = [u.copy()]
    t_save = [0.0]

    t = 0.0
    for n in range(nsteps):
        # Explicit reaction term: u - u^3 = u(1 - u^2)
        reaction = u - u**3

        # Transform to Fourier space
        u_hat = np.fft.fft(u)
        reaction_hat = np.fft.fft(reaction)

        # IMEX step in Fourier space:
        # u_hat^{n+1} = (u_hat^n + dt * reaction_hat) / (1 + k^2 * dt)
        u_hat = (u_hat + dt * reaction_hat) / imex_denom

        # Transform back to physical space
        u = np.fft.ifft(u_hat).real

        t += dt

        # Save snapshot
        if (n + 1) % save_every == 0:
            U_save.append(u.copy())
            t_save.append(t)

    return x, np.array(t_save), np.array(U_save)


def main():
    # Solve Allen-Cahn equation
    N = 512
    eps = 0.05
    tmax = 100.0
    dt = 0.01
    x, t_save, U_save = allen_cahn_imex(N=N, eps=eps, tmax=tmax, dt=dt,
                                         n_snapshots=400)

    # Create figure with two panels
    fig, (ax_main, ax_snap) = plt.subplots(1, 2, figsize=(12, 5))

    # Left panel: space-time pcolormesh plot
    X, T = np.meshgrid(x, t_save)
    pcm = ax_main.pcolormesh(X, T, U_save, cmap='RdBu_r', shading='gouraud',
                              vmin=-1.0, vmax=1.0)
    ax_main.set_xlabel(r'$x$', fontsize=11)
    ax_main.set_ylabel(r'$t$', fontsize=11)
    ax_main.set_title('Allen-Cahn: Phase Separation', fontsize=12)
    fig.colorbar(pcm, ax=ax_main, shrink=0.8, label=r'$u(x, t)$')

    # Right panel: snapshots at selected times
    snapshot_times = [0.0, 2.0, 5.0, 10.0, 30.0, 100.0]
    colors = [SKY, TEAL, '#2ECC71', '#F39C12', CORAL, NAVY]

    for t_target, color in zip(snapshot_times, colors):
        # Find closest snapshot
        idx = np.argmin(np.abs(np.array(t_save) - t_target))
        t_actual = t_save[idx]
        ax_snap.plot(x, U_save[idx], color=color, linewidth=1.5,
                     label=f'$t = {t_actual:.1f}$')

    ax_snap.set_xlabel(r'$x$', fontsize=11)
    ax_snap.set_ylabel(r'$u(x, t)$', fontsize=11)
    ax_snap.set_title('Snapshots at Selected Times', fontsize=12)
    ax_snap.legend(loc='upper right', framealpha=0.9, fontsize=9)
    ax_snap.set_xlim([0, 2 * np.pi])
    ax_snap.set_ylim([-1.3, 1.3])
    ax_snap.axhline(y=1.0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_snap.axhline(y=-1.0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

    # Clean styling
    ax_snap.spines['top'].set_visible(False)
    ax_snap.spines['right'].set_visible(False)

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
