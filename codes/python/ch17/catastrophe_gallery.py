#!/usr/bin/env python3
"""
catastrophe_gallery.py
Chapter 17: Stability of Spectral Methods

Computational Etude 17.1: Gallery of Deliberate Instabilities.

Three-panel figure showing what happens when CFL / stability limits are
deliberately violated in spectral discretisations:

  (a) Fourier advection u_t + u_x = 0 with forward Euler, dt 20 % above CFL.
  (b) Chebyshev wave equation u_tt = u_xx (Dirichlet) with leapfrog,
      dt 50 % above CFL.  Boundary-localised instability develops.
  (c) Chebyshev heat equation u_t = u_xx (Dirichlet) with forward Euler,
      dt doubled above the stability limit.  Immediate blow-up.

Generates Figure 17.1: catastrophe_gallery.pdf

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pathlib import Path

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman', 'CMU Serif', 'DejaVu Serif']
rcParams['mathtext.fontset'] = 'cm'
rcParams['axes.labelsize'] = 11
rcParams['axes.titlesize'] = 12
rcParams['xtick.labelsize'] = 10
rcParams['ytick.labelsize'] = 10
rcParams['legend.fontsize'] = 9
rcParams['figure.dpi'] = 150
rcParams['savefig.dpi'] = 300
rcParams['text.usetex'] = False
rcParams['axes.linewidth'] = 0.8
rcParams['xtick.major.width'] = 0.8
rcParams['ytick.major.width'] = 0.8

NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'
AMBER = '#E67E22'

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch17' / 'python'


def cheb_diff_matrix(N):
    """Chebyshev differentiation matrix of size (N+1) x (N+1)."""
    if N == 0:
        return np.array([[0.0]]), np.array([1.0])
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.ones(N + 1)
    c[0] = 2.0
    c[N] = 2.0
    c *= (-1.0) ** np.arange(N + 1)
    X = np.outer(x, np.ones(N + 1))
    dX = X - X.T + np.eye(N + 1)
    D = np.outer(c, 1.0 / c) / dX
    D -= np.diag(D.sum(axis=1))
    return D, x


def fourier_advection_instability():
    """Forward Euler on Fourier advection, dt above CFL."""
    N = 64
    L = 2.0 * np.pi
    x = np.linspace(0, L, N, endpoint=False)
    k = np.fft.fftfreq(N, d=L / (2.0 * np.pi * N))
    # Wavenumbers
    k = np.fft.fftfreq(N, d=1.0 / N)
    k = 2.0 * np.pi * k / L

    # CFL for forward Euler on advection: dt <= 2 / (max|ik|) doesn't apply
    # Forward Euler on purely imaginary eigenvalues is unconditionally
    # unstable, but with RK4 it is stable up to dt*k_max ~ 2.83.
    # For forward Euler, any dt gives instability on imaginary axis.
    # We use a very small dt that would be "nearly stable" and one 20% larger.
    # Actually, forward Euler on du/dt = iku is |1 + i*k*dt| > 1 for all k!=0.
    # So we just show the growth at a moderate dt.

    u0 = np.exp(-10.0 * (x - np.pi) ** 2)
    dt_stable = 0.5 * L / (N * 1.0)  # reference "CFL" ~ dx/c
    dt = 1.2 * dt_stable  # 20% above

    snapshots_stable = []
    snapshots_unstable = []
    snap_times = [0, 50, 100, 150]

    # "Stable" run with RK4 for reference
    u_hat = np.fft.fft(u0)
    ik = 1j * k
    for n in range(max(snap_times) + 1):
        if n in snap_times:
            snapshots_stable.append(np.real(np.fft.ifft(u_hat)).copy())
        # RK4 step
        k1 = -ik * u_hat
        k2 = -ik * (u_hat + 0.5 * dt * k1)
        k3 = -ik * (u_hat + 0.5 * dt * k2)
        k4 = -ik * (u_hat + dt * k3)
        u_hat = u_hat + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)

    # Unstable run with forward Euler
    u_hat = np.fft.fft(u0)
    for n in range(max(snap_times) + 1):
        if n in snap_times:
            snapshots_unstable.append(np.real(np.fft.ifft(u_hat)).copy())
        u_hat = u_hat + dt * (-ik * u_hat)

    return x, snapshots_stable, snapshots_unstable, snap_times, dt


def chebyshev_wave_instability():
    """Leapfrog on Chebyshev wave equation, dt above CFL."""
    N = 40
    D, x = cheb_diff_matrix(N)
    D2 = D @ D
    # Impose Dirichlet BCs: remove first and last rows/columns
    D2_int = D2[1:N, 1:N]
    x_int = x[1:N]

    # Critical dt for leapfrog on wave equation: dt < 2/rho(D2_int)
    # but rho(D2_int) is very large ~ O(N^4)
    eigs = np.linalg.eigvals(D2_int)
    rho = np.max(np.abs(eigs))
    dt_crit = 2.0 / np.sqrt(rho)
    dt = 1.5 * dt_crit  # 50% above

    # Initial condition: smooth bump
    u = np.exp(-20.0 * x_int ** 2)
    u_prev = u.copy()  # zero velocity start

    n_steps = 60
    snapshots = [np.zeros(N + 1)]
    snap_full = np.zeros(N + 1)
    snap_full[1:N] = u
    snapshots[0] = snap_full.copy()

    snap_times = [0, 20, 40, 60]

    for n in range(1, n_steps + 1):
        if n == 1:
            # First step: use forward Euler for u^1
            u_new = u + 0.5 * dt ** 2 * (D2_int @ u)
        else:
            u_new = 2.0 * u - u_prev + dt ** 2 * (D2_int @ u)
        u_prev = u.copy()
        u = u_new.copy()
        if n in snap_times:
            snap_full = np.zeros(N + 1)
            snap_full[1:N] = u
            snapshots.append(snap_full.copy())

    return x, snapshots, snap_times, dt


def chebyshev_heat_instability():
    """Forward Euler on Chebyshev heat equation, dt above stability limit."""
    N = 20
    D, x = cheb_diff_matrix(N)
    D2 = D @ D
    D2_int = D2[1:N, 1:N]
    x_int = x[1:N]

    eigs = np.linalg.eigvals(D2_int)
    rho = np.max(np.abs(np.real(eigs)))
    dt_crit = 2.0 / rho
    # Use 3x the critical dt -- the amplification factor for the worst
    # eigenvalue is |1 + dt*lambda_min| = |1 - 3*2| = 5 per step,
    # so the outlier mode grows by 5^n and dominates within a few steps.
    dt = 3.0 * dt_crit

    # Initial condition: Gaussian bump (has small but nonzero projection
    # onto the boundary-localised outlier eigenvectors)
    u = np.exp(-5.0 * x_int**2)
    snapshots = []
    snap_times = [0, 3, 5, 6]

    for n in range(max(snap_times) + 1):
        if n in snap_times:
            snap_full = np.zeros(N + 1)
            snap_full[1:N] = u
            snapshots.append(snap_full.copy())
        u = u + dt * (D2_int @ u)

    return x, snapshots, snap_times, dt


def main():
    x_f, snaps_s, snaps_u, times_f, dt_f = fourier_advection_instability()
    x_w, snaps_w, times_w, dt_w = chebyshev_wave_instability()
    x_h, snaps_h, times_h, dt_h = chebyshev_heat_instability()

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 4.2))

    # --- Panel (a): Fourier advection ---
    colors_s = [NAVY, SKY, SKY, SKY]
    colors_u = [NAVY, CORAL, CORAL, CORAL]
    alphas = [1.0, 0.5, 0.7, 1.0]
    for i in range(len(times_f)):
        lbl_s = f'RK4, $n={times_f[i]}$' if i > 0 else 'initial'
        lbl_u = f'FE, $n={times_f[i]}$' if i > 0 else None
        if i == 0:
            ax1.plot(x_f, snaps_s[i], color=NAVY, lw=1.5, label='initial')
        else:
            ax1.plot(x_f, snaps_s[i], color=SKY, lw=1.0, alpha=alphas[i],
                     linestyle='--')
            # Clip unstable solution for display
            u_clip = np.clip(snaps_u[i], -3, 3)
            ax1.plot(x_f, u_clip, color=CORAL, lw=1.0, alpha=alphas[i],
                     label=lbl_u)

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$u(x, t)$')
    ax1.set_title(r'(a) Fourier advection, forward Euler')
    ax1.set_ylim(-3, 3)
    ax1.legend(loc='upper right', fontsize=8)
    ax1.grid(True, alpha=0.3, linewidth=0.5)

    # --- Panel (b): Chebyshev wave equation ---
    colors_w = [NAVY, SKY, AMBER, CORAL]
    for i, snap in enumerate(snaps_w):
        t_idx = [0] + times_w[1:]
        lbl = f'$n = {([0] + times_w[1:])[i]}$'
        snap_clip = np.clip(snap, -5, 5)
        ax2.plot(x_w, snap_clip, color=colors_w[i], lw=1.2, label=lbl)

    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r'$u(x, t)$')
    ax2.set_title(r'(b) Chebyshev wave, leapfrog')
    ax2.set_ylim(-5, 5)
    ax2.legend(loc='upper left', fontsize=8)
    ax2.grid(True, alpha=0.3, linewidth=0.5)

    # --- Panel (c): Chebyshev heat equation ---
    colors_h = [NAVY, SKY, AMBER, CORAL]
    # Find a sensible y-range that shows the blow-up clearly
    all_vals = np.concatenate([s for s in snaps_h])
    ymax_h = min(np.max(np.abs(all_vals)) * 1.1, 50)
    for i, snap in enumerate(snaps_h):
        lbl = f'$n = {times_h[i]}$'
        snap_clip = np.clip(snap, -ymax_h, ymax_h)
        ax3.plot(x_h, snap_clip, color=colors_h[i], lw=1.2, label=lbl)

    ax3.set_xlabel(r'$x$')
    ax3.set_ylabel(r'$u(x, t)$')
    ax3.set_title(r'(c) Chebyshev heat, forward Euler')
    ax3.set_ylim(-ymax_h, ymax_h)
    ax3.legend(loc='upper right', fontsize=8)
    ax3.grid(True, alpha=0.3, linewidth=0.5)

    plt.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'catastrophe_gallery.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'catastrophe_gallery.png', bbox_inches='tight')
    print(f"Figure saved to {OUTPUT_DIR / 'catastrophe_gallery.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
