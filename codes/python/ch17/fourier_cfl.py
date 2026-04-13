#!/usr/bin/env python3
"""
fourier_cfl.py
Chapter 17: Stability of Spectral Methods

Computational Etude 17.3: Fourier CFL Scaling.

Two-panel figure:
  (a) Eigenvalues of the Fourier differentiation matrix for N=32,
      plotted on the imaginary axis from -16i to 16i.
  (b) Log-log plot of critical dt vs N for RK4 applied to Fourier
      advection u_t + u_x = 0.  Binary search for dt_crit.
      Reference line 2.83 / (N/2) shown as dashed.

Generates Figure 17.3: fourier_cfl_scaling.pdf

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

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch17' / 'python'


def rk4_amplification(z):
    """RK4 amplification factor g(z) = 1 + z + z^2/2 + z^3/6 + z^4/24."""
    return 1.0 + z + z ** 2 / 2.0 + z ** 3 / 6.0 + z ** 4 / 24.0


def rk4_max_amplification(dt, N):
    """Maximum |g(dt * lambda)| over all Fourier eigenvalues lambda = -ik.
    For advection u_t + u_x = 0 on [0, 2*pi), eigenvalues are -ik for
    k = -N/2+1, ..., N/2."""
    k = np.fft.fftfreq(N, d=1.0 / N)
    z = -1j * k * dt  # dt * lambda_k
    return np.max(np.abs(rk4_amplification(z)))


def find_dt_crit(N, tol=1e-8):
    """Binary search for critical dt where max|g| = 1."""
    dt_lo = 1e-10
    dt_hi = 1.0
    # Find unstable upper bound
    while rk4_max_amplification(dt_hi, N) <= 1.0:
        dt_hi *= 2
        if dt_hi > 100:
            return dt_hi

    for _ in range(80):
        dt_mid = 0.5 * (dt_lo + dt_hi)
        if rk4_max_amplification(dt_mid, N) <= 1.0:
            dt_lo = dt_mid
        else:
            dt_hi = dt_mid
        if dt_hi - dt_lo < tol * dt_lo:
            break
    return 0.5 * (dt_lo + dt_hi)


def main():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # --- Panel (a): Fourier eigenvalues for N=32 ---
    N = 32
    k = np.fft.fftfreq(N, d=1.0 / N)
    k_sorted = np.sort(k)
    eigs = 1j * k_sorted  # eigenvalues of d/dx

    ax1.plot(eigs.real, eigs.imag, 'o', color=NAVY, markersize=6,
             markerfacecolor=NAVY, markeredgecolor='white',
             markeredgewidth=0.5, zorder=5)
    ax1.axhline(0, color='gray', linewidth=0.5, zorder=0)
    ax1.axvline(0, color='gray', linewidth=0.5, zorder=0)

    # Mark k_max
    ax1.annotate(r'$i k_{\max} = i \cdot 16$',
                 xy=(0, 16), xytext=(0.5, 14),
                 fontsize=10, color=CORAL,
                 arrowprops=dict(arrowstyle='->', color=CORAL, lw=1.2))
    ax1.annotate(r'$-i k_{\max}$',
                 xy=(0, -16), xytext=(0.5, -14),
                 fontsize=10, color=CORAL,
                 arrowprops=dict(arrowstyle='->', color=CORAL, lw=1.2))

    ax1.set_xlabel(r'Re$(\lambda)$')
    ax1.set_ylabel(r'Im$(\lambda)$')
    ax1.set_title(r'(a) Fourier eigenvalues, $N = 32$')
    ax1.set_xlim(-2, 2)
    ax1.set_ylim(-20, 20)
    ax1.grid(True, alpha=0.3, linewidth=0.5)

    # --- Panel (b): dt_crit vs N ---
    N_values = [16, 32, 64, 128, 256]
    dt_crits = []
    for Nv in N_values:
        dtc = find_dt_crit(Nv)
        dt_crits.append(dtc)
        print(f"N = {Nv:4d}, dt_crit = {dtc:.6f}, "
              f"dt_crit * (N/2) = {dtc * Nv / 2:.4f}")

    dt_crits = np.array(dt_crits)
    N_arr = np.array(N_values, dtype=float)

    ax2.loglog(N_arr, dt_crits, 'o-', color=NAVY, markersize=7, lw=1.5,
               label=r'$\Delta t_{\rm crit}$ (binary search)')

    # Reference line: 2.83 / (N/2)
    N_ref = np.linspace(12, 300, 100)
    ax2.loglog(N_ref, 2.83 / (N_ref / 2), '--', color=CORAL, lw=1.2,
               label=r'$2.83\, /\, (N/2)$')

    ax2.set_xlabel(r'$N$')
    ax2.set_ylabel(r'$\Delta t_{\rm crit}$')
    ax2.set_title(r'(b) RK4 CFL scaling for Fourier advection')
    ax2.legend(loc='upper right', fontsize=9)
    ax2.grid(True, which='both', alpha=0.3, linewidth=0.5)

    # Show slope
    ax2.text(80, 0.15, r'slope $= -1$', fontsize=11, color=CORAL,
             rotation=-28, ha='center')

    plt.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'fourier_cfl_scaling.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'fourier_cfl_scaling.png', bbox_inches='tight')
    print(f"\nFigure saved to {OUTPUT_DIR / 'fourier_cfl_scaling.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
