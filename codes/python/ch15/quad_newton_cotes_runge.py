#!/usr/bin/env python3
"""
quad_newton_cotes_runge.py
Chapter 15: Quadrature in Spectral Methods

Computational Etude 15.2: The Runge Phenomenon in Quadrature.

Integrates the Runge function f(x) = 1/(1 + 25x^2) over [-1, 1] using
three quadrature rules of increasing order n:

    1. Newton--Cotes (equispaced nodes, Vandermonde-based weights).
    2. Gauss--Legendre (from numpy.polynomial.legendre.leggauss).
    3. Clenshaw--Curtis (Chebyshev nodes, FFT-based weights).

The exact integral is (2/5) arctan(5).

For Newton--Cotes, the weights are computed by solving the Vandermonde
system that enforces exactness on monomials 1, x, ..., x^n.  Because the
Vandermonde matrix becomes exponentially ill-conditioned for equispaced
nodes, the Newton--Cotes weights develop large alternating values and the
quadrature error grows without bound: this is the quadrature analogue of
the Runge phenomenon in polynomial interpolation.

In contrast, both Gauss--Legendre and Clenshaw--Curtis converge, though at
different rates for this particular integrand.

Reference: L. N. Trefethen, "Is Gauss Quadrature Better than
           Clenshaw--Curtis?", SIAM Review 50(1):67--87, 2008.

Generates Figure 15.2: newton_cotes_runge.pdf

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


def newton_cotes_weights(n):
    """
    Compute Newton--Cotes weights for n+1 equispaced nodes on [-1, 1].

    The weights satisfy sum_j w_j x_j^k = integral_{-1}^{1} x^k dx for
    k = 0, 1, ..., n.  This requires solving a Vandermonde system.  For
    large n the system is very ill-conditioned, which causes the weights
    to blow up: this is the source of the Runge phenomenon in quadrature.
    """
    x = np.linspace(-1, 1, n + 1)
    # Right-hand side: exact integrals of monomials
    k = np.arange(n + 1, dtype=float)
    rhs = np.where(k % 2 == 0, 2.0 / (k + 1), 0.0)
    # Vandermonde matrix: V[i, j] = x_j^i
    V = np.vander(x, increasing=True).T
    # Solve for weights
    w = np.linalg.solve(V, rhs)
    return x, w


def clenshaw_curtis_weights(n):
    """
    Compute Clenshaw--Curtis quadrature nodes and weights for n+1 points
    on [-1, 1] via DCT-I implemented through the FFT.

    Nodes: x_j = cos(j pi / n), j = 0, ..., n.

    The weights are computed from the Chebyshev moments
    mu_k = int_{-1}^{1} T_k(x) dx = 2/(1-k^2) for even k, 0 for odd k,
    using the DCT-I relation between weights and moments.
    """
    if n == 0:
        return np.array([0.0]), np.array([2.0])
    if n == 1:
        return np.array([1.0, -1.0]), np.array([1.0, 1.0])

    x = np.cos(np.pi * np.arange(n + 1) / n)

    # Chebyshev moments: mu_k = int_{-1}^{1} T_k(x) dx
    c = np.zeros(n + 1)
    c[0::2] = 2.0 / (1.0 - np.arange(0, n + 1, 2)**2)

    # DCT-I via FFT of length 2n (mirror the sequence)
    v = np.concatenate([c, c[-2:0:-1]])  # length 2n
    f = np.real(np.fft.fft(v))

    # Extract weights with endpoint halving
    w = f[:n + 1] / n
    w[0] /= 2
    w[-1] /= 2

    return x, w


def main():
    """Generate the Runge phenomenon in quadrature figure."""
    # The Runge function and its exact integral
    f = lambda x: 1.0 / (1.0 + 25.0 * x**2)
    I_exact = (2.0 / 5.0) * np.arctan(5.0)

    # Range of n values
    n_values = np.arange(4, 52, 2)

    err_nc = np.zeros(len(n_values))
    err_gl = np.zeros(len(n_values))
    err_cc = np.zeros(len(n_values))

    for i, n in enumerate(n_values):
        # Newton--Cotes
        x_nc, w_nc = newton_cotes_weights(n)
        I_nc = np.dot(w_nc, f(x_nc))
        err_nc[i] = abs(I_nc - I_exact)

        # Gauss--Legendre
        x_gl, w_gl = np.polynomial.legendre.leggauss(n + 1)
        I_gl = np.dot(w_gl, f(x_gl))
        err_gl[i] = abs(I_gl - I_exact)

        # Clenshaw--Curtis
        x_cc, w_cc = clenshaw_curtis_weights(n)
        I_cc = np.dot(w_cc, f(x_cc))
        err_cc[i] = abs(I_cc - I_exact)

    # Clamp zeros for log plot
    err_nc = np.maximum(err_nc, 1e-16)
    err_gl = np.maximum(err_gl, 1e-16)
    err_cc = np.maximum(err_cc, 1e-16)

    # Create figure
    fig, ax = plt.subplots(figsize=(7, 5))

    ax.semilogy(n_values, err_nc, 'o-', color=CORAL, markersize=4,
                linewidth=1.2, label=r'Newton$\mathrm{-}$Cotes')
    ax.semilogy(n_values, err_gl, 's-', color=NAVY, markersize=4,
                linewidth=1.2, label=r'Gauss$\mathrm{-}$Legendre')
    ax.semilogy(n_values, err_cc, 'D-', color=TEAL, markersize=4,
                linewidth=1.2, label=r'Clenshaw$\mathrm{-}$Curtis')

    ax.set_xlabel(r'Number of nodes $n$')
    ax.set_ylabel(r'Absolute error $|I - I_n|$')
    ax.set_title(r'Quadrature of $f(x) = 1/(1 + 25x^2)$ on $[-1, 1]$')
    ax.legend(loc='best', frameon=True, fancybox=False, edgecolor='gray')
    ax.grid(True, alpha=0.3, linewidth=0.5)

    # Set axis limits
    ax.set_xlim(n_values[0] - 1, n_values[-1] + 1)

    plt.tight_layout()

    # Save output
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'newton_cotes_runge.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'newton_cotes_runge.png', bbox_inches='tight')
    print(f"Figure saved to {OUTPUT_DIR / 'newton_cotes_runge.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
