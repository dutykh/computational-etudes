#!/usr/bin/env python3
"""
advection_fourier.py
Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs

Fourier pseudospectral solver for the linear advection equation:
    u_t + c * u_x = 0 on [0, 2*pi), periodic BCs.

Compares Fourier pseudospectral + RK4 with second-order upwind FD + RK4
to demonstrate the dispersion-free nature of spectral advection.

Initial condition: u(x, 0) = exp(-10*(x - pi)^2)

Generates figure: advection_comparison.pdf

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
OUTPUT_FILE = OUTPUT_DIR / 'advection_comparison.pdf'


def rhs_fourier(u, k, c):
    """Right-hand side for Fourier pseudospectral advection: -c * u_x."""
    u_hat = np.fft.fft(u)
    ux = np.fft.ifft(1j * k * u_hat).real
    return -c * ux


def rhs_fd_upwind2(u, dx, c):
    """Right-hand side for second-order upwind FD advection: -c * u_x."""
    N = len(u)
    ux = np.zeros(N)
    if c > 0:
        # Second-order upwind (backward biased)
        for j in range(N):
            ux[j] = (3.0 * u[j] - 4.0 * u[(j - 1) % N] + u[(j - 2) % N]) / (2.0 * dx)
    else:
        for j in range(N):
            ux[j] = (-3.0 * u[j] + 4.0 * u[(j + 1) % N] - u[(j + 2) % N]) / (2.0 * dx)
    return -c * ux


def rk4_step(u, dt, rhs_func):
    """One step of classical RK4."""
    k1 = rhs_func(u)
    k2 = rhs_func(u + 0.5 * dt * k1)
    k3 = rhs_func(u + 0.5 * dt * k2)
    k4 = rhs_func(u + dt * k3)
    return u + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def main():
    # Parameters
    N = 128
    c = 1.0
    L = 2.0 * np.pi
    t_final = 2.0 * np.pi / c  # One full period

    # Grid
    x = L * np.arange(N) / N
    dx = L / N

    # Wavenumbers
    k = np.fft.fftfreq(N, d=1.0 / N).astype(float)
    if N % 2 == 0:
        k[N // 2] = 0  # Zero Nyquist mode for odd derivatives

    # Initial condition (smooth localised bump, periodised)
    u0 = np.exp(-10.0 * (x - np.pi)**2)

    # Time step (CFL-based)
    k_max = N / 2
    dt = 0.5 / (c * k_max)  # CFL ~ 0.5
    n_steps = int(np.ceil(t_final / dt))
    dt = t_final / n_steps

    # Solve with Fourier pseudospectral
    u_fourier = u0.copy()
    for _ in range(n_steps):
        u_fourier = rk4_step(u_fourier, dt, lambda u: rhs_fourier(u, k, c))

    # Solve with second-order upwind FD (same time step and RK4)
    u_fd = u0.copy()
    for _ in range(n_steps):
        u_fd = rk4_step(u_fd, dt, lambda u: rhs_fd_upwind2(u, dx, c))

    # Exact solution (one full period = initial condition)
    u_exact = u0.copy()

    # -------------------------------------------------------------------------
    # Figure
    # -------------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    # Left: overlay of solutions
    ax1.plot(x, u_exact, '-', color='gray', linewidth=2.5, alpha=0.5,
             label='Exact (initial)')
    ax1.plot(x, u_fourier, '-', color=NAVY, linewidth=1.5,
             label='Fourier PS + RK4')
    ax1.plot(x, u_fd, '--', color=CORAL, linewidth=1.5,
             label='Upwind FD2 + RK4')
    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$u(x, t)$')
    ax1.set_title(f'After one period ($t = 2\\pi / c$, $N = {N}$)')
    ax1.legend(loc='upper right')
    ax1.set_xlim([0, L])
    ax1.grid(True, alpha=0.3)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Right: error profiles
    err_fourier = np.abs(u_fourier - u_exact)
    err_fd = np.abs(u_fd - u_exact)

    ax2.semilogy(x, err_fd, '-', color=CORAL, linewidth=1.5,
                 label=f'Upwind FD2 (max err: {np.max(err_fd):.2e})')
    ax2.semilogy(x, np.maximum(err_fourier, 1e-16), '-', color=NAVY, linewidth=1.5,
                 label=f'Fourier PS (max err: {np.max(err_fourier):.2e})')
    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r'$|u^{\mathrm{num}} - u^{\mathrm{exact}}|$')
    ax2.set_title('Pointwise Error After One Period')
    ax2.legend(loc='upper right')
    ax2.set_xlim([0, L])
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
