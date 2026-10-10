#!/usr/bin/env python3
"""
trap_subgeometric.py
Chapter 16: Integration of Periodic Functions

Computational Etude 16.6: Subgeometric Decay (C^infty but not Analytic).

Weideman 2002 (Example 6) constructs the function

    f_6(x) = exp((cos x - 1) / (cos x + 1))    for x != pi,
             0                                  at x = pi.

It is 2*pi-periodic and infinitely differentiable, yet has an essential
singularity at x = pi.  Therefore it lies in NO strip of analyticity
around the real axis, so the strip-analyticity theorem (Theorem 3.1 of
Trefethen-Weideman 2014) cannot give a geometric error bound.

The Euler-Maclaurin formula nevertheless predicts a convergence rate
faster than any algebraic power 1/N^k, since all odd derivatives vanish
at the endpoints of the period.  The actual rate, derived from a
careful asymptotic analysis of the Meijer G-function representation of
the Fourier coefficients of f_6, is

    |I(f_6) - T_N(f_6)| ~ 8 sqrt(pi/3) * exp(-(3/2) * N^{2/3})

(up to oscillatory factors), the canonical example of "subgeometric"
convergence -- faster than any polynomial but slower than r^N for any
0 < r < 1.

The diagnostic for subgeometric decay is to plot the log error against
N^{2/3}: it becomes a straight line.  We do that here.

Exact integral: I(f_6) = 2 * e * pi * (1 - erf(1)) ~ 2.68658684...

Generates Figure 16.6: subgeometric.pdf  (two-panel figure)

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
from scipy.special import erf

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


def f6(x):
    """f_6 from Weideman (2002), Example 6 -- with safe handling at x = pi."""
    c = np.cos(x)
    # Safe evaluation: at x = pi we have cos x = -1, and the function
    # extends to zero by the analytic limit.
    out = np.zeros_like(x)
    mask = (c + 1.0) > 1e-15
    out[mask] = np.exp((c[mask] - 1.0) / (c[mask] + 1.0))
    return out


def main():
    # Reference value: I = 2 * e * pi * (1 - erf(1))
    I_exact = 2.0 * np.exp(1.0) * np.pi * (1.0 - erf(1.0))
    print(f"Exact integral: 2*e*pi*(1 - erf(1)) = {I_exact:.16f}")
    print("Reference (Weideman 2002):           2.68658684...\n")

    N_values = np.array([4, 6, 8, 10, 12, 16, 20, 24, 30, 40, 50, 60, 80, 100, 120, 160, 200])
    errors = np.zeros(len(N_values))
    for i, N in enumerate(N_values):
        errors[i] = abs(trapezoidal_periodic(f6, N) - I_exact)
    errors = np.maximum(errors, 1e-17)

    # Theoretical rate: 8 sqrt(pi/3) * exp(-(3/2) * N^{2/3})
    theory = 8.0 * np.sqrt(np.pi / 3.0) * np.exp(-1.5 * N_values**(2.0 / 3.0))

    print(f"{'N':>4s}  {'error':>14s}  {'envelope':>14s}")
    for i, N in enumerate(N_values):
        print(f"{N:4d}  {errors[i]:14.4e}  {theory[i]:14.4e}")

    # ------------------------------------------------------------------
    # Two-panel plot
    # ------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # Panel 1: error vs N (semilog)
    ax1.semilogy(N_values, errors, 'o-', color=CORAL, markersize=5,
                 linewidth=1.2, label=r'Trapezoidal error')
    ax1.semilogy(N_values, theory, '--', color=NAVY, linewidth=1.0,
                 label=r'$8\sqrt{\pi/3}\,e^{-3 N^{2/3}/2}$')
    ax1.set_xlabel(r'Number of nodes $N$')
    ax1.set_ylabel(r'Absolute error $|I_N - I|$')
    ax1.set_title(r'(a) error vs $N$: looks slower than geometric')
    ax1.legend(loc='upper right', frameon=True, fancybox=False,
               edgecolor='gray', framealpha=0.9)
    ax1.grid(True, which='both', alpha=0.3, linewidth=0.5)
    ax1.set_ylim(1e-17, 1e0)

    # Panel 2: error vs N^{2/3} (semilog) -- should be a straight line
    N23 = N_values.astype(float)**(2.0 / 3.0)
    ax2.semilogy(N23, errors, 'o-', color=CORAL, markersize=5,
                 linewidth=1.2, label=r'Trapezoidal error')
    ax2.semilogy(N23, theory, '--', color=NAVY, linewidth=1.0,
                 label=r'$8\sqrt{\pi/3}\,e^{-3 N^{2/3}/2}$')
    ax2.set_xlabel(r'$N^{2/3}$')
    ax2.set_ylabel(r'Absolute error $|I_N - I|$')
    ax2.set_title(r'(b) error vs $N^{2/3}$: linear on semilog')
    ax2.legend(loc='upper right', frameon=True, fancybox=False,
               edgecolor='gray', framealpha=0.9)
    ax2.grid(True, which='both', alpha=0.3, linewidth=0.5)
    ax2.set_ylim(1e-17, 1e0)

    plt.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'subgeometric.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'subgeometric.png', bbox_inches='tight')
    print(f"\nFigure saved to {OUTPUT_DIR / 'subgeometric.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
