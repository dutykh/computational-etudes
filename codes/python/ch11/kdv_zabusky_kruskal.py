#!/usr/bin/env python3
"""
kdv_zabusky_kruskal.py
Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs

Fourier pseudospectral solver for the KdV equation with integrating factor RK4:
    u_t + 6 u u_x + u_xxx = 0 on [0, 2), periodic.

Replicates the Zabusky-Kruskal (1965) experiment: a cosine initial condition
breaks into a train of solitons, demonstrating the FPU recurrence phenomenon.

Initial condition: u(x, 0) = cos(pi * x)

Generates figure: kdv_zabusky_kruskal.pdf (space-time diagram)

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
OUTPUT_FILE = OUTPUT_DIR / 'kdv_zabusky_kruskal.pdf'


def kdv_zabusky_kruskal(N=256, L=2.0, t_final=3.6 / np.pi, dt=1e-4):
    """
    Solve u_t + u u_x + delta^2 u_xxx = 0 on [0, L), periodic,
    with u(x, 0) = cos(pi * x) and delta^2 = 0.022^2 (Zabusky-Kruskal form).

    Rewritten as: u_t + 6*u*u_x + u_xxx = 0 on [0, L) after rescaling.
    We use the original ZK normalisation for historical fidelity.

    Uses integrating factor RK4 with 2/3 dealiasing.
    """
    delta2 = 0.022**2

    # Grid and wavenumbers
    x = L * np.arange(N) / N
    k = (2.0 * np.pi / L) * np.fft.fftfreq(N, d=1.0 / N).astype(float)

    # Linear symbol: lambda_k = -i * delta^2 * k^3
    lam = -1j * delta2 * k**3

    # 2/3 dealiasing mask
    cutoff = N // 3
    mask = np.ones(N, dtype=bool)
    if cutoff > 0:
        mask[cutoff + 1:N - cutoff] = False

    # Initial condition
    u0 = np.cos(np.pi * x)
    v_hat = np.fft.fft(u0).copy()

    def G(t_now, v_hat_now):
        """RHS in integrating factor variables."""
        E = np.exp(lam * t_now)
        u_hat = E * v_hat_now
        u = np.fft.ifft(u_hat).real
        ux = np.fft.ifft(1j * k * u_hat).real
        f_nl = -u * ux  # Nonlinear term: -u * u_x
        f_nl_hat = np.fft.fft(f_nl)
        f_nl_hat[~mask] = 0.0
        return np.exp(-lam * t_now) * f_nl_hat

    # Time stepping
    n_steps = int(np.ceil(t_final / dt))
    dt = t_final / n_steps

    n_save = 300
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
    N = 256
    t_final = 3.6 / np.pi  # Standard ZK experiment time

    x, t_save, U_save = kdv_zabusky_kruskal(N=N, t_final=t_final, dt=5e-5)

    # -------------------------------------------------------------------------
    # Figure: Space-time diagram + snapshots
    # -------------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Left: space-time plot
    X, T = np.meshgrid(x, t_save)
    pcm = ax1.pcolormesh(X, T, U_save, cmap='RdBu_r', shading='gouraud')
    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$t$')
    ax1.set_title('Zabusky--Kruskal Experiment')
    fig.colorbar(pcm, ax=ax1, shrink=0.8, label=r'$u(x, t)$')

    # Right: snapshots at selected times
    t_targets = [0.0, 0.3, 0.6, 0.9, t_final]
    colors = [SKY, TEAL, '#2ECC71', '#F39C12', NAVY]

    for t_target, color in zip(t_targets, colors):
        idx = np.argmin(np.abs(t_save - t_target))
        ax2.plot(x, U_save[idx], color=color, linewidth=1.3,
                 label=f'$t = {t_save[idx]:.2f}$')

    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r'$u(x, t)$')
    ax2.set_title('Soliton Train Formation')
    ax2.legend(loc='upper right', framealpha=0.9, fontsize=9)
    ax2.set_xlim([0, 2.0])
    ax2.grid(True, alpha=0.3)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    fig.tight_layout()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_FILE, bbox_inches='tight', pad_inches=0.05)
    print(f'Figure saved: {OUTPUT_FILE.resolve()}')
    plt.close(fig)


if __name__ == '__main__':
    main()
