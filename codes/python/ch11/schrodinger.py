#!/usr/bin/env python3
"""
schrodinger.py
Chapter 10: Spectral PDE Solvers with Chebyshev and Fourier Grids

Time-dependent Schrodinger equation with harmonic potential using Strang splitting.

PDE: i u_t = -u_xx + x^2 u,  0 < x < 2*pi,  t > 0 (periodic)
ICs: u(x, 0) = exp(-5*(x - pi)^2) * exp(3ix)  (Gaussian wavepacket)

Method: Strang splitting (second-order, unitary)
  1. Half potential step (physical space): u -> u * exp(-i x^2 dt/2)
  2. Full kinetic step (Fourier space):   u_hat -> u_hat * exp(-i k^2 dt)
  3. Half potential step (physical space): u -> u * exp(-i x^2 dt/2)

Generates Figure 10.X: Schrodinger equation probability density evolution.

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
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch10' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'schrodinger.pdf'


def schrodinger_strang(N=256, tmax=3.0, n_snapshots=200):
    """
    Solve the time-dependent Schrodinger equation with harmonic potential
    using Strang splitting and FFT.

    The Strang splitting preserves the L^2 norm (unitarity) and is
    second-order accurate in time.

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
            Snapshots of |u(x,t)|^2 (probability density)
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

    # Potential: V(x) = x^2
    V = x**2

    # Initial condition: Gaussian wavepacket with momentum
    u = np.exp(-5.0 * (x - np.pi)**2) * np.exp(3j * x)

    # Time stepping
    dt = 0.5 * h  # Time step (stable for splitting)
    nsteps = int(np.ceil(tmax / dt))
    dt = tmax / nsteps  # Adjust to hit tmax exactly

    # Precompute exponentials for splitting operators
    half_potential = np.exp(-0.5j * V * dt)    # Half potential step
    full_kinetic = np.exp(-1j * k**2 * dt)     # Full kinetic step

    # Store snapshots
    save_every = max(1, nsteps // n_snapshots)
    U_save = [np.abs(u)**2]
    t_save = [0.0]

    # Record initial L^2 norm for unitarity check
    norm0 = np.sqrt(h * np.sum(np.abs(u)**2))

    t = 0.0
    for n in range(nsteps):
        # Strang splitting: potential/2 -> kinetic -> potential/2
        u = half_potential * u                    # Half potential step
        u_hat = np.fft.fft(u)
        u_hat = full_kinetic * u_hat              # Full kinetic step
        u = np.fft.ifft(u_hat)
        u = half_potential * u                    # Half potential step

        t += dt

        # Save snapshot
        if (n + 1) % save_every == 0:
            U_save.append(np.abs(u)**2)
            t_save.append(t)

    # Check unitarity preservation
    norm_final = np.sqrt(h * np.sum(np.abs(u)**2))
    print(f"Initial L2 norm: {norm0:.10f}")
    print(f"Final   L2 norm: {norm_final:.10f}")
    print(f"Relative change: {abs(norm_final - norm0) / norm0:.2e}")

    return x, np.array(t_save), np.array(U_save)


def main():
    # Solve Schrodinger equation
    N = 256
    tmax = 3.0
    x, t_save, U_save = schrodinger_strang(N=N, tmax=tmax, n_snapshots=200)

    # Create figure with two panels
    fig, (ax_main, ax_water) = plt.subplots(1, 2, figsize=(12, 5))

    # Left panel: space-time plot of |u|^2
    X, T = np.meshgrid(x, t_save)
    pcm = ax_main.pcolormesh(X, T, U_save, cmap='inferno', shading='gouraud')
    ax_main.set_xlabel(r'$x$', fontsize=11)
    ax_main.set_ylabel(r'$t$', fontsize=11)
    ax_main.set_title(r'Probability Density $|\psi(x,t)|^2$', fontsize=12)
    fig.colorbar(pcm, ax=ax_main, shrink=0.8, label=r'$|\psi|^2$')

    # Right panel: waterfall view
    n_lines = 40
    indices = np.linspace(0, len(t_save) - 1, n_lines, dtype=int)

    for i, idx in enumerate(indices):
        offset = t_save[idx] * 0.3
        color = plt.cm.viridis(i / n_lines)
        ax_water.plot(x, U_save[idx] + offset, color=color, linewidth=0.8)

    ax_water.set_xlabel(r'$x$', fontsize=11)
    ax_water.set_ylabel(r'$|\psi(x,t)|^2$ (offset by $t$)', fontsize=11)
    ax_water.set_title(r'Waterfall View: $|\psi|^2$ Evolution', fontsize=12)

    # Add inset showing the potential V(x) = x^2
    ax_inset = ax_water.inset_axes([0.6, 0.7, 0.35, 0.25])
    ax_inset.plot(x, x**2, color=CORAL, linewidth=1.5)
    ax_inset.set_xlabel(r'$x$', fontsize=8)
    ax_inset.set_ylabel(r'$V(x)$', fontsize=8)
    ax_inset.set_title('Harmonic Potential', fontsize=9)
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
