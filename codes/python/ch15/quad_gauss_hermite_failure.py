#!/usr/bin/env python3
"""
quad_gauss_hermite_failure.py
Chapter 15: Quadrature in Spectral Methods

Computational Etude 15.9: Gauss--Hermite vs Truncated Trapezoidal Rule.

Compares two quadrature strategies for the oscillatory integral

    I = int_{-inf}^{inf} e^{-x^2} cos(x^3) dx

which is a prototypical example of a Gaussian-weighted integral with a
non-polynomial integrand.

Gauss--Hermite quadrature:  The weight function e^{-x^2} is built into
the rule, so I_n = sum_k w_k cos(x_k^3). Theory predicts root-exponential
convergence: errors decay like O(exp(-C sqrt(n))), which appears as a
parabolic curve on a semilog plot of error vs n.

Truncated trapezoidal rule: The full integrand g(x) = e^{-x^2} cos(x^3)
is integrated on a finite interval [-L, L] using the composite trapezoidal
rule. Because g(x) decays super-algebraically as |x| -> inf, the
truncation error is negligible for moderate L. The trapezoidal rule then
converges exponentially in n (straight line on semilog), dramatically
outperforming Gauss--Hermite.

This etude demonstrates that Gauss rules matched to the weight function
are not always the best choice; a simpler rule on a truncated domain
can be far more effective when the integrand decays rapidly.

Generates Figure 15.9: gauss_hermite_failure.pdf

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
from scipy.special import roots_hermite
from scipy.integrate import quad as scipy_quad

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

# -----------------------------------------------------------------------------
# Book color scheme
# -----------------------------------------------------------------------------
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'
PURPLE = '#8E44AD'
ORANGE = '#E67E22'

# Output paths
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch15' / 'python'


def main():
    """Generate the Gauss--Hermite vs trapezoidal comparison figure."""

    # -------------------------------------------------------------------------
    # Reference value via high-accuracy adaptive quadrature
    # -------------------------------------------------------------------------
    integrand = lambda x: np.exp(-x**2) * np.cos(x**3)
    I_exact, _ = scipy_quad(integrand, -np.inf, np.inf,
                            limit=500, epsabs=1e-15, epsrel=1e-15)
    print(f"Reference integral value: {I_exact:.16e}")

    # -------------------------------------------------------------------------
    # Gauss--Hermite quadrature: I_n = sum w_k cos(x_k^3)
    # The weight function e^{-x^2} is absorbed into the weights.
    # -------------------------------------------------------------------------
    n_gh = np.arange(10, 1001, 10)
    err_gh = np.zeros(len(n_gh))

    for i, n in enumerate(n_gh):
        x_h, w_h = roots_hermite(n)
        I_gh = np.dot(w_h, np.cos(x_h**3))
        err_gh[i] = abs(I_gh - I_exact)

    # -------------------------------------------------------------------------
    # Truncated trapezoidal rule on [-L, L]
    # g(x) = e^{-x^2} cos(x^3), with L = 6 (truncation error < 10^{-15})
    # -------------------------------------------------------------------------
    L = 6.0
    n_trap = np.arange(10, 201, 2)
    err_trap = np.zeros(len(n_trap))

    for i, n in enumerate(n_trap):
        h = 2.0 * L / n
        x = np.linspace(-L, L, n + 1)
        g = np.exp(-x**2) * np.cos(x**3)
        # Composite trapezoidal rule
        I_trap = h * (0.5 * g[0] + np.sum(g[1:-1]) + 0.5 * g[-1])
        err_trap[i] = abs(I_trap - I_exact)

    # -------------------------------------------------------------------------
    # Create figure
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 5))

    floor = 1e-16

    # Gauss--Hermite errors
    ax.semilogy(n_gh, np.maximum(err_gh, floor), 'o', color=NAVY,
                markersize=3, markeredgecolor=NAVY, markeredgewidth=0.5,
                label=r'Gauss$\mathrm{-}$Hermite')

    # Trapezoidal errors
    ax.semilogy(n_trap, np.maximum(err_trap, floor), 's', color=CORAL,
                markersize=3.5, markerfacecolor='white',
                markeredgecolor=CORAL, markeredgewidth=0.7,
                label=r'Trapezoidal on $[-6, 6]$')

    # Machine epsilon reference
    ax.axhline(y=np.finfo(float).eps, color='gray', linestyle=':',
               linewidth=0.8, alpha=0.7)
    ax.text(950, 4e-16, r'$\epsilon_{\mathrm{mach}}$', fontsize=9,
            color='gray', ha='right')

    ax.set_xlabel(r'Number of quadrature points $n$')
    ax.set_ylabel(r'Absolute error $|I - I_n|$')
    ax.set_title(
        r'Quadrature for $\int_{-\infty}^{\infty} e^{-x^2} \cos(x^3)\, dx$'
    )
    ax.legend(loc='upper right', frameon=True, fancybox=False,
              edgecolor='gray', fontsize=10)
    ax.grid(True, alpha=0.3, linewidth=0.5)
    ax.set_xlim(0, 1050)
    ax.set_ylim(1e-16, 1e1)

    plt.tight_layout()

    # Save output
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'gauss_hermite_failure.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'gauss_hermite_failure.png', bbox_inches='tight')
    print(f"Figure saved to {OUTPUT_DIR / 'gauss_hermite_failure.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
