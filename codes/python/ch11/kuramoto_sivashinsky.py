#!/usr/bin/env python3
"""
kuramoto_sivashinsky.py
Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs

Fourier pseudospectral solver for the Kuramoto-Sivashinsky equation:
    u_t + u u_x + u_xx + u_xxxx = 0 on [0, L), periodic.

Uses integrating factor RK4 to handle the stiff linear terms exactly.
The linear symbol is L(k) = k^2 - k^4 (energy production at low k,
damping at high k).

Generates two figures:
  - ks_spacetime.pdf: space-time "fingerprint" plot
  - ks_spectrum.pdf: time-averaged energy spectrum

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


def kuramoto_sivashinsky(N=256, L=32.0 * np.pi, t_final=200.0, dt=0.05):
    """
    Solve u_t + u u_x + u_xx + u_xxxx = 0 on [0, L), periodic,
    using integrating factor RK4 with 2/3 dealiasing.

    Linear part: L(k) = k^2 - k^4 (treated exactly via IF)
    Nonlinear part: N(u) = -u u_x (treated explicitly)

    In conservation form: N(u) = -(1/2)(u^2)_x so N_hat = -ik/2 * FFT(u^2).
    """
    # Grid and wavenumbers
    x = L * np.arange(N) / N
    k = (2.0 * np.pi / L) * np.fft.fftfreq(N, d=1.0 / N).astype(float)

    # Linear symbol: lambda_k = k^2 - k^4
    lam = k**2 - k**4

    # 2/3 dealiasing mask
    cutoff = N // 3
    mask = np.ones(N, dtype=bool)
    if cutoff > 0:
        mask[cutoff + 1:N - cutoff] = False

    # Initial condition: small random perturbation
    np.random.seed(42)
    u0 = 0.01 * np.random.randn(N)
    # Low-pass filter
    u_hat = np.fft.fft(u0)
    u_hat[np.abs(k) > 5.0 * (2.0 * np.pi / L)] = 0.0
    u0 = np.fft.ifft(u_hat).real

    # Per-step integrating factor (Trefethen's approach, avoids overflow)
    # E = exp(L*dt/2), E2 = exp(L*dt) applied per step, not cumulatively
    u_hat = np.fft.fft(u0).copy()

    n_steps = int(np.ceil(t_final / dt))
    dt = t_final / n_steps

    E = np.exp(lam * dt / 2.0)   # half-step propagator
    E2 = E**2                     # full-step propagator

    def Nhat(u_hat_now):
        """Nonlinear term in Fourier space: -ik/2 * FFT(u^2)."""
        u = np.fft.ifft(u_hat_now).real
        f_hat = -0.5j * k * np.fft.fft(u**2)
        f_hat[~mask] = 0.0
        return f_hat

    n_save = 500
    save_every = max(1, n_steps // n_save)
    snapshots = [u0.copy()]
    snap_times = [0.0]

    t = 0.0
    for n in range(n_steps):
        # IF-RK4 per step (Trefethen SMIM style)
        Na = dt * Nhat(u_hat)
        Nb = dt * Nhat(E * u_hat + Na / 2.0)
        Nc = dt * Nhat(E * u_hat + Nb / 2.0)
        Nd = dt * Nhat(E2 * u_hat + E * Nc)
        u_hat = E2 * u_hat + (E2 * Na + 2.0 * E * (Nb + Nc) + Nd) / 6.0
        t += dt

        if (n + 1) % save_every == 0:
            u = np.fft.ifft(u_hat).real
            snapshots.append(u.copy())
            snap_times.append(t)

    return x, np.array(snap_times), np.array(snapshots)


def main():
    N = 256
    L = 32.0 * np.pi
    t_final = 200.0

    x, t_save, U_save = kuramoto_sivashinsky(N=N, L=L, t_final=t_final, dt=0.05)

    # -------------------------------------------------------------------------
    # Figure 1: Space-time "fingerprint"
    # -------------------------------------------------------------------------
    fig1, ax1 = plt.subplots(figsize=(8, 6))

    X, T = np.meshgrid(x, t_save)
    vmax = np.percentile(np.abs(U_save[len(U_save) // 4:]), 98)
    pcm = ax1.pcolormesh(X, T, U_save, cmap='RdBu_r', shading='gouraud',
                          vmin=-vmax, vmax=vmax)
    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$t$')
    ax1.set_title(f'Kuramoto--Sivashinsky: Spatiotemporal Chaos ($L = {L / np.pi:.0f}\\pi$)')
    fig1.colorbar(pcm, ax=ax1, shrink=0.8, label=r'$u(x, t)$')

    fig1.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig1.savefig(OUTPUT_DIR / 'ks_spacetime.pdf',
                 bbox_inches='tight', pad_inches=0.05)
    print(f'Figure saved: {(OUTPUT_DIR / "ks_spacetime.pdf").resolve()}')
    plt.close(fig1)

    # -------------------------------------------------------------------------
    # Figure 2: Time-averaged energy spectrum
    # -------------------------------------------------------------------------
    # Average over the second half (statistically steady state)
    n_start = len(U_save) // 2
    spectra = np.zeros(N // 2)
    for i in range(n_start, len(U_save)):
        u_hat = np.fft.fft(U_save[i])
        spectra += np.abs(u_hat[:N // 2])**2 / N**2
    spectra /= (len(U_save) - n_start)

    k_pos = (2.0 * np.pi / L) * np.arange(N // 2)

    fig2, ax2 = plt.subplots(figsize=(6, 4.5))
    ax2.semilogy(k_pos[1:], spectra[1:], '-', color=NAVY, linewidth=1.5)
    ax2.set_xlabel(r'Wavenumber $k$')
    ax2.set_ylabel(r'$\langle |\hat{u}_k|^2 \rangle$')
    ax2.set_title('Time-Averaged Energy Spectrum (KS)')
    ax2.grid(True, alpha=0.3, which='both')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    fig2.tight_layout()
    fig2.savefig(OUTPUT_DIR / 'ks_spectrum.pdf',
                 bbox_inches='tight', pad_inches=0.05)
    print(f'Figure saved: {(OUTPUT_DIR / "ks_spectrum.pdf").resolve()}')
    plt.close(fig2)


if __name__ == '__main__':
    main()
