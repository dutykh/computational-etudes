#!/usr/bin/env python3
"""
quad_gauss_hermite_weights.py
Chapter 15: Quadrature in Spectral Methods

Computational Etude 15.8: Gauss--Hermite Weight Distribution.

Examines the pathological behaviour of Gauss--Hermite quadrature weights
for large n, which limits the practical utility of this classical rule
for high-accuracy integration.

Left panel: For n = 100, the Gauss--Hermite weights w_k are plotted
on a semilogarithmic scale against the corresponding nodes x_k. The
weights span a vast dynamic range: the central weights are of O(1),
while the outermost weights fall far below machine epsilon (approx.
10^{-16}). This means that when the weights are represented in
double-precision floating-point arithmetic, a significant fraction
carry no useful information.

Right panel: The fraction of weights below machine epsilon is plotted
as a function of n for n = 10, 20, ..., 500. This fraction increases
monotonically and approaches 1, showing that Gauss--Hermite quadrature
becomes numerically useless in standard double-precision for large n.

This motivates the search for alternative quadrature strategies on
unbounded domains, such as the truncated Gauss--Legendre approach
discussed in subsequent etudes.

Generates Figure 15.8: gauss_hermite_weights.pdf

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
    """Generate the two-panel Gauss--Hermite weight figure."""
    eps = np.finfo(float).eps  # Machine epsilon (~2.22e-16)

    # -------------------------------------------------------------------------
    # Left panel: weights for n = 100
    # -------------------------------------------------------------------------
    n_demo = 100
    x_h, w_h = roots_hermite(n_demo)

    # -------------------------------------------------------------------------
    # Right panel: fraction of weights below eps vs n
    # -------------------------------------------------------------------------
    n_values = np.arange(10, 501, 10)
    fraction_below = np.zeros(len(n_values))

    for i, n in enumerate(n_values):
        _, w = roots_hermite(n)
        fraction_below[i] = np.sum(np.abs(w) < eps) / n

    # -------------------------------------------------------------------------
    # Create figure
    # -------------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # --- Left panel: weight distribution ---
    ax1.semilogy(x_h, np.abs(w_h), 'o', color=NAVY, markersize=2.5,
                 markeredgecolor=NAVY, markeredgewidth=0.3,
                 label=r'$|w_k|$')

    # Machine epsilon reference line
    ax1.axhline(y=eps, color=CORAL, linestyle='--', linewidth=1.0,
                label=r'Machine $\epsilon \approx 2.2 \times 10^{-16}$')

    ax1.set_xlabel(r'Node $x_k$')
    ax1.set_ylabel(r'Weight $|w_k|$')
    ax1.set_title(r'(a) Gauss$\mathrm{-}$Hermite weights ($n = 100$)')
    ax1.legend(loc='lower center', frameon=True, fancybox=False,
               edgecolor='gray', fontsize=9)
    ax1.grid(True, alpha=0.3, linewidth=0.5)
    ax1.set_ylim(1e-300, 1e2)

    # --- Right panel: fraction below epsilon ---
    ax2.plot(n_values, fraction_below, 'o-', color=NAVY, markersize=3,
             linewidth=1.2, markeredgecolor=NAVY, markeredgewidth=0.5)

    ax2.axhline(y=1.0, color='gray', linestyle=':', linewidth=0.8)

    ax2.set_xlabel(r'Number of points $n$')
    ax2.set_ylabel(r'Fraction of $|w_k| < \epsilon$')
    ax2.set_title(r'(b) Fraction of underflowed weights')
    ax2.set_xlim(0, 510)
    ax2.set_ylim(-0.02, 1.05)
    ax2.grid(True, alpha=0.3, linewidth=0.5)

    plt.tight_layout()

    # Save output
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'gauss_hermite_weights.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'gauss_hermite_weights.png', bbox_inches='tight')
    print(f"Figure saved to {OUTPUT_DIR / 'gauss_hermite_weights.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
