#!/usr/bin/env python3
"""
kdv_soliton.py
Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs

Fourier pseudospectral solver for the KdV equation with integrating factor RK4:
    u_t + 6 u u_x + u_xxx = 0 on [0, L), periodic.

Simulates the elastic collision of two solitons of different amplitudes.
Soliton solution: u(x,t) = 2*a^2 * sech^2(a*(x - x0 - 4*a^2*t))

Generates figure: kdv_snapshots.pdf (snapshots before, during, after collision)

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
OUTPUT_FILE = OUTPUT_DIR / 'kdv_snapshots.pdf'


def soliton(x, a, x0):
    """KdV soliton: 2*a^2 * sech^2(a*(x - x0))."""
    return 2.0 * a**2 / np.cosh(a * (x - x0))**2


def kdv_ifrk4(N=512, L=80.0, a1=1.0, a2=0.5, x1=20.0, x2=40.0,
              t_final=25.0, dt=0.005, dealias=True):
    """
    Integrating factor RK4 solver for KdV:
        u_t + 6 u u_x + u_xxx = 0
    on [0, L) with periodic BCs.

    Parameters:
        N : int - grid points
        L : float - domain length
        a1, a2 : float - soliton parameters (amplitude = 2*a^2, speed = 4*a^2)
        x1, x2 : float - initial soliton positions
        t_final : float - final time
        dt : float - time step
        dealias : bool - apply 2/3 dealiasing

    Returns:
        x, t_save, U_save
    """
    # Grid and wavenumbers
    x = L * np.arange(N) / N
    k = (2.0 * np.pi / L) * np.fft.fftfreq(N, d=1.0 / N).astype(float)

    # Linear symbol
    lam = 1j * k**3

    # Dealiasing mask
    cutoff = N // 3
    mask = np.ones(N, dtype=bool)
    if dealias and cutoff > 0:
        mask[cutoff + 1:N - cutoff] = False

    # Initial condition: two solitons
    u0 = soliton(x, a1, x1) + soliton(x, a2, x2)
    v_hat = np.fft.fft(u0).copy()

    def G(t_now, v_hat_now):
        """RHS in integrating factor variables."""
        E = np.exp(lam * t_now)
        u_hat = E * v_hat_now
        u = np.fft.ifft(u_hat).real
        ux = np.fft.ifft(1j * k * u_hat).real
        f_nl = -6.0 * u * ux
        f_nl_hat = np.fft.fft(f_nl)
        f_nl_hat[~mask] = 0.0
        return np.exp(-lam * t_now) * f_nl_hat

    n_steps = int(np.ceil(t_final / dt))
    dt = t_final / n_steps

    n_save = 200
    save_every = max(1, n_steps // n_save)
    snapshots = [u0.copy()]
    snap_times = [0.0]

    t = 0.0
    for n in range(n_steps):
        k1 = G(t, v_hat)
        k2 = G(t + 0.5 * dt, v_hat + 0.5 * dt * k1)
        k3 = G(t + 0.5 * dt, v_hat + 0.5 * dt * k2)
        k4 = G(t + dt, v_hat + dt * k3)
        v_hat = v_hat + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        t += dt

        if (n + 1) % save_every == 0:
            u_hat = np.exp(lam * t) * v_hat
            u = np.fft.ifft(u_hat).real
            snapshots.append(u.copy())
            snap_times.append(t)

    return x, np.array(snap_times), np.array(snapshots)


def main():
    N = 512
    L = 80.0
    a1, a2 = 1.0, 0.5
    x1, x2 = 20.0, 40.0  # Taller soliton behind
    t_final = 25.0

    x, t_save, U_save = kdv_ifrk4(
        N=N, L=L, a1=a1, a2=a2, x1=x1, x2=x2, t_final=t_final, dt=0.002
    )

    # -------------------------------------------------------------------------
    # Figure: snapshots before, during, and after collision
    # -------------------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # Select snapshot times
    t_targets = [0.0, 8.0, 12.0, t_final]
    titles = ['Initial condition ($t = 0$)',
              'Before collision',
              'During collision',
              f'After collision ($t = {t_final:.0f}$)']
    colors = [NAVY, TEAL, CORAL, NAVY]

    for ax, t_target, title, color in zip(axes.flat, t_targets, titles, colors):
        idx = np.argmin(np.abs(t_save - t_target))
        ax.plot(x, U_save[idx], color=color, linewidth=1.5)
        ax.set_title(title)
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$u(x, t)$')
        ax.set_xlim([0, L])
        ax.set_ylim([-0.2, 2.5])
        ax.grid(True, alpha=0.3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle(
        f'KdV Two-Soliton Collision ($a_1 = {a1}$, $a_2 = {a2}$, $N = {N}$)',
        fontsize=13, y=1.01
    )

    fig.tight_layout()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_FILE, bbox_inches='tight', pad_inches=0.05)
    print(f'Figure saved: {OUTPUT_FILE.resolve()}')
    plt.close(fig)


if __name__ == '__main__':
    main()
