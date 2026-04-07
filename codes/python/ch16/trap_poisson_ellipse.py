#!/usr/bin/env python3
"""
trap_poisson_ellipse.py
Chapter 16: Integration of Periodic Functions

Computational Etude 16.1: Poisson's Ellipse, the Original Paradox.

In the 1820s Poisson computed the perimeter of an ellipse with semi-axes
1/(2*pi) and 0.8/(2*pi) using the trapezoidal rule and obtained ten
correct decimal digits from just a handful of function evaluations.  The
integral is

    I = (1/(2*pi)) * int_0^{2*pi} sqrt(1 - 0.36 sin^2(theta)) d(theta)
      = (2/pi) * E(0.36)
      = 0.902779927772...

where E is the complete elliptic integral of the second kind.  The
integrand has branch points in the complex plane at theta = +/- i log(3),
so by the strip-analyticity theorem of Trefethen-Weideman the trapezoidal
error is O(3^{-N}).  Each new sample point adds about log10(3) ~ 0.48
correct digits, so the convergence is geometric.

This script reproduces Figure 1.1 of Trefethen and Weideman, SIAM Rev.
2014, exploiting the four-fold symmetry of the integrand to halve the
number of samples and using a high-precision reference value computed
from scipy's elliptic integral.

Generates Figure 16.1: poisson_ellipse.pdf

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
from scipy.special import ellipe

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
rcParams['xtick.major.width'] = 0.8
rcParams['ytick.major.width'] = 0.8

# Book color scheme
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'

# Output paths
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch16' / 'python'


def trapezoidal_periodic(f, N):
    """N-point periodic trapezoidal rule on [0, 2*pi)."""
    theta = 2.0 * np.pi * np.arange(N) / N
    return (2.0 * np.pi / N) * np.sum(f(theta))


def main():
    # The integrand of Poisson's example
    integrand = lambda t: np.sqrt(1.0 - 0.36 * np.sin(t)**2)

    # Reference value: I = (2/pi) * E(0.36)
    # scipy convention: ellipe takes m = k^2, here m = 0.36
    I_exact = (2.0 / np.pi) * ellipe(0.36)
    print(f"Exact value:  I = {I_exact:.16f}")
    print(f"Reference:        0.9027799277721857...\n")

    # Trapezoidal rule: divide the integral by 2*pi to match the
    # normalisation of the exact value above.
    # Use N = 4, 8, 12, ..., 200 (multiples of 4 to exploit symmetry)
    N_values = np.arange(4, 201, 4)
    errors = np.zeros(len(N_values))

    for i, N in enumerate(N_values):
        I_N = trapezoidal_periodic(integrand, N) / (2.0 * np.pi)
        errors[i] = abs(I_N - I_exact)

    # Clamp at machine epsilon for log plotting
    errors = np.maximum(errors, 1e-17)

    # Theoretical decay: 3^{-N} (i.e. 3^{-N/4} per quarter step)
    theory = 3.0 ** (-N_values)
    theory = np.maximum(theory, 1e-17)

    # Print a small table matching Figure 1.1 of Trefethen-Weideman
    print(f"{'N/4':>4s}  {'I_N':>22s}  {'|I_N - I|':>14s}")
    for i, N in enumerate(N_values[:10]):
        I_N = trapezoidal_periodic(integrand, N) / (2.0 * np.pi)
        print(f"{N//4:4d}  {I_N:22.18f}  {errors[i]:14.4e}")

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 4.8))

    ax.semilogy(N_values // 4, errors, 'o-', color=CORAL,
                markersize=4, linewidth=1.1,
                markerfacecolor=CORAL, markeredgecolor=CORAL,
                label=r'Trapezoidal rule')
    ax.semilogy(N_values // 4, theory, '--', color=NAVY,
                linewidth=1.0, label=r'$3^{-N}$')

    ax.set_xlabel(r'$N/4$ (number of independent samples)')
    ax.set_ylabel(r'Absolute error $|I_N - I|$')
    ax.set_title(r"Poisson's ellipse: trapezoidal convergence is geometric")
    ax.legend(loc='upper right', frameon=True, fancybox=False,
              edgecolor='gray', framealpha=0.9)
    ax.grid(True, which='both', alpha=0.3, linewidth=0.5)
    ax.set_xlim(0, 50)
    ax.set_ylim(1e-18, 1e0)

    plt.tight_layout()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'poisson_ellipse.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'poisson_ellipse.png', bbox_inches='tight')
    print(f"\nFigure saved to {OUTPUT_DIR / 'poisson_ellipse.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
