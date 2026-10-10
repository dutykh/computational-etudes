#!/usr/bin/env python3
"""
nls_recurrence.py
Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs

Split-step Fourier solver for the focusing nonlinear Schrodinger equation:
    i u_t + u_xx + 2 |u|^2 u = 0 on [0, 2*pi), periodic.

Demonstrates modulation instability and approximate recurrence using
Strang splitting (linear-nonlinear-linear half-steps).

Initial condition: u(x, 0) = 1 + 0.1 * cos(x)

Generates figure: nls_recurrence.pdf

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
OUTPUT_FILE = OUTPUT_DIR / 'nls_recurrence.pdf'


def nls_split_step(N=256, t_final=5.0, dt=1e-3):
    """
    Solve i u_t + u_xx + 2|u|^2 u = 0 on [0, 2*pi) with periodic BCs
    using the Strang split-step Fourier method.

    Linear step:  u_hat_k <- exp(-i k^2 dt/2) * u_hat_k
    Nonlinear step: u <- exp(2i |u|^2 dt) * u

    Returns:
        x, t_save, U_save (|u(x,t)|)
    """
    x = 2.0 * np.pi * np.arange(N) / N
    k = np.fft.fftfreq(N, d=1.0 / N).astype(float)

    # Initial condition: modulated plane wave
    u = (1.0 + 0.1 * np.cos(x)).astype(complex)

    # Precompute linear propagator (half-step)
    linear_half = np.exp(-1j * k**2 * dt / 2.0)

    n_steps = int(np.ceil(t_final / dt))
    dt = t_final / n_steps

    # Recompute propagator with adjusted dt
    linear_half = np.exp(-1j * k**2 * dt / 2.0)

    n_save = 400
    save_every = max(1, n_steps // n_save)
    U_save = [np.abs(u).copy()]
    L2_save = [np.sqrt(np.mean(np.abs(u)**2))]
    t_save = [0.0]

    t = 0.0
    for n in range(n_steps):
        # Strang splitting: L/2 -> N -> L/2
        u_hat = np.fft.fft(u)
        u_hat *= linear_half
        u = np.fft.ifft(u_hat)

        # Full nonlinear step
        u *= np.exp(2j * np.abs(u)**2 * dt)

        # Second linear half-step
        u_hat = np.fft.fft(u)
        u_hat *= linear_half
        u = np.fft.ifft(u_hat)

        t += dt

        if (n + 1) % save_every == 0:
            U_save.append(np.abs(u).copy())
            L2_save.append(np.sqrt(np.mean(np.abs(u)**2)))
            t_save.append(t)

    return x, np.array(t_save), np.array(U_save), np.array(L2_save)


def main():
    N = 256
    t_final = 5.0
    dt = 5e-4

    x, t_save, U_save, L2_save = nls_split_step(N=N, t_final=t_final, dt=dt)

    # -------------------------------------------------------------------------
    # Figure: space-time diagram + L2 norm conservation
    # -------------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5),
                                    gridspec_kw={'width_ratios': [2, 1]})

    # Left: space-time diagram of |u|
    X, T = np.meshgrid(x, t_save)
    pcm = ax1.pcolormesh(X, T, U_save, cmap='inferno', shading='gouraud')
    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$t$')
    ax1.set_title(r'NLS: $|u(x, t)|$ (Modulation Instability)')
    fig.colorbar(pcm, ax=ax1, shrink=0.8, label=r'$|u|$')

    # Right: L2 norm conservation
    L2_0 = L2_save[0]
    ax2.plot(t_save, np.abs(L2_save - L2_0) / L2_0, '-', color=NAVY,
             linewidth=1.0)
    ax2.set_xlabel(r'$t$')
    ax2.set_ylabel(r'$|\|u\|_2 - \|u_0\|_2| / \|u_0\|_2$')
    ax2.set_title(r'$L^2$ Norm Conservation')
    ax2.set_yscale('log')
    ax2.grid(True, alpha=0.3, which='both')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    fig.tight_layout()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_FILE, bbox_inches='tight', pad_inches=0.05)
    print(f'Figure saved: {OUTPUT_FILE.resolve()}')
    plt.close(fig)


if __name__ == '__main__':
    main()
