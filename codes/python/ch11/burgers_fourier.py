#!/usr/bin/env python3
"""
burgers_fourier.py
Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs

Fourier pseudospectral solver for the viscous Burgers equation:
    u_t + u * u_x = nu * u_xx on [0, 2*pi), periodic BCs.

Initial condition: u(x, 0) = sin(x)
Features: RK4 time stepping, 2/3 dealiasing rule.

Generates two figures:
  - burgers_evolution.pdf: solution evolution showing front steepening
  - burgers_dealiasing.pdf: comparison of energy spectra with/without dealiasing

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


def burgers_fourier(N=256, nu=1e-2, t_final=1.0, CFL=0.4, dealias=True):
    """
    Fourier pseudospectral solver for viscous Burgers equation
        u_t + u u_x = nu u_xx
    on [0, 2*pi) with periodic boundary conditions.

    Parameters:
        N : int - number of grid points (even)
        nu : float - viscosity
        t_final : float - final simulation time
        CFL : float - CFL-like safety factor
        dealias : bool - apply 2/3 dealiasing rule

    Returns:
        x, t_save, U_save : grid, snapshot times, solution snapshots
    """
    # Grid and wavenumbers
    x = 2.0 * np.pi * np.arange(N) / N
    k = np.fft.fftfreq(N, d=1.0 / N).astype(float)

    # Initial condition
    u = np.sin(x)

    # Time step from linearised stability
    k_max = N / 2.0
    dt_adv = CFL / k_max
    dt_diff = 0.4 / (nu * k_max**2) if nu > 0 else np.inf
    dt = min(dt_adv, dt_diff)
    n_steps = int(np.ceil(t_final / dt))
    dt = t_final / n_steps

    # Diffusive symbol
    diff_symbol = -nu * k**2

    # 2/3 dealiasing mask
    dealias_mask = np.ones(N, dtype=bool)
    if dealias:
        cutoff = N // 3
        if cutoff > 0:
            dealias_mask[cutoff + 1:N - cutoff] = False

    def rhs(u_phys):
        u_hat = np.fft.fft(u_phys)
        ux = np.fft.ifft(1j * k * u_hat).real

        # Nonlinear term in physical space, then dealiase
        f_nl = -u_phys * ux
        f_nl_hat = np.fft.fft(f_nl)
        f_nl_hat[~dealias_mask] = 0.0

        # Diffusive term in Fourier space
        f_diff_hat = diff_symbol * u_hat

        return np.fft.ifft(f_nl_hat + f_diff_hat).real

    # Time integration with RK4
    n_save = min(200, n_steps)
    save_every = max(1, n_steps // n_save)
    U_save = [u.copy()]
    t_save = [0.0]

    t = 0.0
    for n in range(n_steps):
        k1 = rhs(u)
        k2 = rhs(u + 0.5 * dt * k1)
        k3 = rhs(u + 0.5 * dt * k2)
        k4 = rhs(u + dt * k3)
        u = u + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        t += dt

        if (n + 1) % save_every == 0:
            U_save.append(u.copy())
            t_save.append(t)

    return x, np.array(t_save), np.array(U_save)


def main():
    N = 256
    nu = 0.01
    t_final = 1.5

    # -------------------------------------------------------------------------
    # Figure 1: Solution evolution
    # -------------------------------------------------------------------------
    x, t_save, U_save = burgers_fourier(N=N, nu=nu, t_final=t_final, dealias=True)

    fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Left: space-time plot
    X, T = np.meshgrid(x, t_save)
    pcm = ax1.pcolormesh(X, T, U_save, cmap='RdBu_r', shading='gouraud')
    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$t$')
    ax1.set_title(f'Burgers Equation ($\\nu = {nu}$)')
    fig1.colorbar(pcm, ax=ax1, shrink=0.8, label=r'$u(x, t)$')

    # Right: snapshots
    snapshot_times = [0.0, 0.3, 0.6, 1.0, 1.5]
    colors = [SKY, TEAL, '#2ECC71', '#F39C12', NAVY]

    for t_target, color in zip(snapshot_times, colors):
        idx = np.argmin(np.abs(t_save - t_target))
        ax2.plot(x, U_save[idx], color=color, linewidth=1.5,
                 label=f'$t = {t_save[idx]:.1f}$')

    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r'$u(x, t)$')
    ax2.set_title('Snapshots: Front Steepening')
    ax2.legend(loc='upper right', framealpha=0.9, fontsize=9)
    ax2.set_xlim([0, 2 * np.pi])
    ax2.grid(True, alpha=0.3)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    fig1.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig1.savefig(OUTPUT_DIR / 'burgers_evolution.pdf',
                 bbox_inches='tight', pad_inches=0.05)
    print(f'Figure saved: {(OUTPUT_DIR / "burgers_evolution.pdf").resolve()}')
    plt.close(fig1)

    # -------------------------------------------------------------------------
    # Figure 2: Dealiasing comparison
    # -------------------------------------------------------------------------
    t_compare = 1.0

    _, t1, U1 = burgers_fourier(N=N, nu=nu, t_final=t_compare, dealias=True)
    _, t2, U2 = burgers_fourier(N=N, nu=nu, t_final=t_compare, dealias=False)

    # Energy spectra at final time
    u_hat_deal = np.fft.fft(U1[-1])
    u_hat_nodeal = np.fft.fft(U2[-1])
    spectrum_deal = np.abs(u_hat_deal[:N // 2])**2 / N**2
    spectrum_nodeal = np.abs(u_hat_nodeal[:N // 2])**2 / N**2
    k_pos = np.arange(N // 2)

    fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    # Left: energy spectra
    ax1.semilogy(k_pos[1:], spectrum_deal[1:], '-', color=NAVY, linewidth=1.5,
                 label='With 2/3 dealiasing')
    ax1.semilogy(k_pos[1:], spectrum_nodeal[1:], '--', color=CORAL, linewidth=1.5,
                 label='Without dealiasing')
    ax1.set_xlabel(r'Wavenumber $k$')
    ax1.set_ylabel(r'$|\hat{u}_k|^2 / N^2$')
    ax1.set_title(f'Energy Spectrum at $t = {t_compare}$')
    ax1.legend()
    ax1.grid(True, alpha=0.3, which='both')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Right: physical space comparison
    ax2.plot(x, U1[-1], '-', color=NAVY, linewidth=1.5,
             label='With dealiasing')
    ax2.plot(x, U2[-1], '--', color=CORAL, linewidth=1.5,
             label='Without dealiasing')
    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r'$u(x, t)$')
    ax2.set_title(f'Physical Space at $t = {t_compare}$')
    ax2.legend()
    ax2.set_xlim([0, 2 * np.pi])
    ax2.grid(True, alpha=0.3)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    fig2.tight_layout()
    fig2.savefig(OUTPUT_DIR / 'burgers_dealiasing.pdf',
                 bbox_inches='tight', pad_inches=0.05)
    print(f'Figure saved: {(OUTPUT_DIR / "burgers_dealiasing.pdf").resolve()}')
    plt.close(fig2)


if __name__ == '__main__':
    main()
