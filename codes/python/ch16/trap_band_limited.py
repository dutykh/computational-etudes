#!/usr/bin/env python3
"""
trap_band_limited.py
Chapter 16: Integration of Periodic Functions

Computational Etude 16.2: Band-Limited Exactness and One-Mode Aliasing.

The aliasing identity for the periodic trapezoidal rule states that
exp(i*j*theta) integrates to 2*pi*delta_{j,0} on the continuous level
and to 2*pi*sum_{l in Z} delta_{j, l*N} on the N-point grid.  Therefore:

  (a) If f is a trigonometric polynomial of degree m, then for any N > m
      the trapezoidal rule is EXACT (error at machine precision).

  (b) For a single mode exp(i*k*theta), the trapezoidal sum still gives
      the correct value of zero whenever k is not an integer multiple
      of N -- even when k is between N/2 and N (so the mode is aliased
      and cannot be reconstructed by interpolation).

This second observation resolves the apparent paradox of "one point per
wavelength versus two points per wavelength": reconstruction needs two,
but integration only needs one, because aliased nonzero modes still
integrate to zero.

This script demonstrates both phenomena:
  - Top panel: build a random trigonometric polynomial of degree m = 10
    and plot the trapezoidal error vs N from 4 to 30.  The error is
    machine-precision zero for all N > 10.
  - Bottom panel: take the single mode cos(15*theta) and the trapezoidal
    rule with N = 16 (so 15 is between N/2 = 8 and N = 16, an aliased
    frequency).  Plot how the error behaves as we sweep k from 0 to 30
    while N is held fixed at 16.

Generates Figure 16.2: band_limited.pdf

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
rcParams['legend.fontsize'] = 10
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
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch16' / 'python'


def trapezoidal_periodic(f, N):
    """N-point periodic trapezoidal rule on [0, 2*pi)."""
    theta = 2.0 * np.pi * np.arange(N) / N
    return (2.0 * np.pi / N) * np.sum(f(theta))


def main():
    # ------------------------------------------------------------------
    # Part (a): random trigonometric polynomial of degree m = 10
    # ------------------------------------------------------------------
    rng = np.random.default_rng(seed=42)
    m = 10
    # Real coefficients of cos(k*theta) and sin(k*theta), k = 1..m
    a_coeffs = rng.standard_normal(m + 1)  # a_0..a_m
    b_coeffs = rng.standard_normal(m + 1)  # b_1..b_m (b_0 unused)

    def trig_poly(theta):
        result = np.full_like(theta, a_coeffs[0])  # a_0
        for k in range(1, m + 1):
            result += a_coeffs[k] * np.cos(k * theta)
            result += b_coeffs[k] * np.sin(k * theta)
        return result

    # Exact integral of a trigonometric polynomial: only the constant
    # term contributes.  int_0^{2*pi} a_0 dtheta = 2*pi * a_0.
    I_exact = 2.0 * np.pi * a_coeffs[0]

    N_values_a = np.arange(4, 31)
    errors_a = np.zeros(len(N_values_a))
    for i, N in enumerate(N_values_a):
        errors_a[i] = abs(trapezoidal_periodic(trig_poly, N) - I_exact)
    errors_a = np.maximum(errors_a, 1e-17)

    print(f"Part (a): random trigonometric polynomial of degree m = {m}")
    print(f"Exact integral: 2*pi*a_0 = {I_exact:.10f}")
    print(f"Error at N = {m}:    {errors_a[m - 4]:.4e}")
    print(f"Error at N = {m+1}:    {errors_a[m + 1 - 4]:.4e}")
    print(f"Error at N = {m+5}:   {errors_a[m + 5 - 4]:.4e}")

    # ------------------------------------------------------------------
    # Part (b): single mode cos(k*theta), N = 16 fixed, sweep k
    # ------------------------------------------------------------------
    N_fixed = 16
    k_values = np.arange(0, 33)
    errors_b = np.zeros(len(k_values))
    for i, k in enumerate(k_values):
        f_k = lambda t, kk=k: np.cos(kk * t)
        # Exact integral: 2*pi if k = 0, else 0
        I_exact_k = 2.0 * np.pi if k == 0 else 0.0
        errors_b[i] = abs(trapezoidal_periodic(f_k, N_fixed) - I_exact_k)
    errors_b = np.maximum(errors_b, 1e-17)

    print(f"\nPart (b): single mode cos(k*theta), N = {N_fixed} fixed")
    print(f"Aliased frequencies (k = N, 2N, ...) where the rule fails:")
    for i, k in enumerate(k_values):
        if errors_b[i] > 1e-10:
            print(f"  k = {k}:  error = {errors_b[i]:.4e}")

    # ------------------------------------------------------------------
    # Two-panel plot
    # ------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # Panel (a)
    ax1.semilogy(N_values_a, errors_a, 'o-', color=CORAL,
                 markersize=5, linewidth=1.1,
                 markerfacecolor=CORAL, markeredgecolor=CORAL)
    ax1.axvline(x=m, color='gray', linestyle='--', linewidth=0.9, zorder=0)
    ax1.text(m + 0.3, 1e-3, r'exactness threshold $N = m$',
             rotation=90, va='center', ha='left', fontsize=9, color='gray')
    ax1.set_xlabel(r'Number of nodes $N$')
    ax1.set_ylabel(r'Absolute error $|I_N - I|$')
    ax1.set_title(r'(a) trigonometric polynomial of degree $m = 10$')
    ax1.grid(True, which='both', alpha=0.3, linewidth=0.5)
    ax1.set_xlim(3, 31)
    ax1.set_ylim(1e-18, 1e2)

    # Panel (b)
    colors_b = [CORAL if e > 1e-10 else NAVY for e in errors_b]
    ax2.bar(k_values, errors_b, color=colors_b, edgecolor='none')
    ax2.axvline(x=N_fixed, color='gray', linestyle='--', linewidth=0.9,
                zorder=0)
    ax2.axvline(x=2 * N_fixed, color='gray', linestyle='--', linewidth=0.9,
                zorder=0)
    ax2.text(N_fixed + 0.3, 4.0, r'$k = N$',
             rotation=90, va='center', ha='left', fontsize=9, color='gray')
    ax2.text(2 * N_fixed + 0.3, 4.0, r'$k = 2N$',
             rotation=90, va='center', ha='left', fontsize=9, color='gray')
    ax2.set_xlabel(r'Frequency $k$')
    ax2.set_ylabel(r'$|I_N(\cos k\theta) - I(\cos k\theta)|$')
    ax2.set_title(r'(b) single modes, fixed $N = 16$')
    ax2.set_yscale('symlog', linthresh=1e-15)
    ax2.set_ylim(0, 1e1)
    ax2.set_xlim(-0.5, 32.5)
    ax2.grid(True, which='major', alpha=0.3, linewidth=0.5)

    plt.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'band_limited.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'band_limited.png', bbox_inches='tight')
    print(f"\nFigure saved to {OUTPUT_DIR / 'band_limited.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
