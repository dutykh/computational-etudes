#!/usr/bin/env python3
"""
trap_supergeometric.py
Chapter 16: Integration of Periodic Functions

Computational Etude 16.5: Supergeometric Decay on an Entire Function.

The function f_5(x) = exp(cos x) is 2*pi-periodic and entire (analytic
in the whole complex plane).  Its exact integral over [0, 2*pi] is

    int_0^{2*pi} exp(cos x) dx = 2*pi * I_0(1) = 7.95492652101284...,

where I_0 is the modified Bessel function of order zero.  Because f_5
is entire, no fixed strip of analyticity bounds its convergence rate;
instead, Weideman 2002 derives the explicit asymptotic

    I(f_5) - T_N(f_5) ~ -4*pi*(-1)^N * I_N(1)
                       ~ 2*sqrt(2*pi/N) * (e/(2N))^N,

which is "faster than any geometric rate" -- a so-called supergeometric
convergence.  Each new sample point gives roughly one extra digit, and
in fact slightly more than that.

We reproduce the convergence table from Trefethen-Weideman SIAM Rev 2014
Section 4 and verify it against the asymptotic estimate.

Generates Figure 16.5: supergeometric.pdf

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
from scipy.special import iv

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
    f5 = lambda x: np.exp(np.cos(x))
    I_exact = 2.0 * np.pi * iv(0, 1.0)
    print(f"Exact integral: 2*pi*I_0(1) = {I_exact:.16f}")

    N_values = np.arange(1, 17)
    errors = np.zeros(len(N_values))
    for i, N in enumerate(N_values):
        errors[i] = abs(trapezoidal_periodic(f5, N) - I_exact)
    errors = np.maximum(errors, 1e-17)

    # Asymptotic from Weideman: 2 * sqrt(2*pi/N) * (e/(2N))^N
    e = np.exp(1.0)
    theory = 2.0 * np.sqrt(2.0 * np.pi / N_values) * (e / (2.0 * N_values))**N_values
    theory = np.maximum(theory, 1e-17)

    # Print the table
    print(f"\n{'N':>3s}  {'I_N':>22s}  {'|I_N - I|':>14s}  {'theory':>14s}")
    for i, N in enumerate(N_values):
        I_N = trapezoidal_periodic(f5, N)
        print(f"{N:3d}  {I_N:22.16f}  {errors[i]:14.4e}  {theory[i]:14.4e}")

    # ------------------------------------------------------------------
    # Plot: error vs N on a semilog scale, with the asymptotic curve
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 5))

    ax.semilogy(N_values, errors, 'o-', color=CORAL, markersize=6,
                linewidth=1.2,
                markerfacecolor=CORAL, markeredgecolor=CORAL,
                label=r'Trapezoidal error')
    ax.semilogy(N_values, theory, '--', color=NAVY, linewidth=1.0,
                label=r'$2\sqrt{2\pi/N}\,(e/2N)^N$')

    ax.set_xlabel(r'Number of nodes $N$')
    ax.set_ylabel(r'Absolute error $|I_N - I|$')
    ax.set_title(r'Supergeometric decay for $f(x) = e^{\cos x}$')
    ax.legend(loc='upper right', frameon=True, fancybox=False,
              edgecolor='gray', framealpha=0.9)
    ax.grid(True, which='both', alpha=0.3, linewidth=0.5)
    ax.set_xlim(0, 17)
    ax.set_ylim(1e-18, 1e2)

    plt.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'supergeometric.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'supergeometric.png', bbox_inches='tight')
    print(f"\nFigure saved to {OUTPUT_DIR / 'supergeometric.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
