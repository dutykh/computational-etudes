#!/usr/bin/env python3
"""
quad_polynomial_exactness.py
Chapter 15: Quadrature in Spectral Methods

Computational Etude 15.1: Polynomial Exactness Test.

For a fixed number of nodes n + 1 = 33, computes the absolute monomial
quadrature error

    |E_n(x^k)| = | sum_j w_j x_j^k - integral_{-1}^{1} x^k dx |

for k = 0, 1, ..., 2n = 64, using three classical interpolatory rules:

    1. Newton--Cotes (equispaced nodes, Vandermonde-based weights).
    2. Gauss--Legendre (numpy.polynomial.legendre.leggauss).
    3. Clenshaw--Curtis (Chebyshev nodes, FFT-based weights).

The exact integral is

    int_{-1}^{1} x^k dx = 2/(k+1)   for even k,
                        = 0          for odd k.

Theoretical degrees of precision:

    Newton--Cotes:    n   (here, k <= 32)
    Clenshaw--Curtis: n   (here, k <= 32)
    Gauss--Legendre:  2n + 1   (here, k <= 65, i.e. the entire test range)

The figure makes these degrees visible: NC and CC sit at machine precision
for k <= n and lift off slowly for larger k, while Gauss stays at the
machine-precision floor across the whole range. Odd k are exactly zero for
all three rules by symmetry and are clamped to the floor for plotting.

Generates Figure 15.1b: polynomial_exactness.pdf

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
    """Compute monomial-exactness errors and produce the figure."""
    n = 32  # n + 1 = 33 quadrature points
    k_vals = np.arange(0, 2 * n + 1)  # 0, 1, ..., 64

    # Build the three quadrature rules once
    x_nc, w_nc = newton_cotes_weights(n)
    x_gl, w_gl = np.polynomial.legendre.leggauss(n + 1)
    x_cc, w_cc = clenshaw_curtis_weights(n)

    err_nc = np.empty(len(k_vals), dtype=float)
    err_gl = np.empty(len(k_vals), dtype=float)
    err_cc = np.empty(len(k_vals), dtype=float)

    for i, k in enumerate(k_vals):
        # Exact integral: 2/(k+1) for even k, 0 for odd k
        I_exact = 2.0 / (k + 1) if k % 2 == 0 else 0.0
        err_nc[i] = abs(np.dot(w_nc, x_nc ** k) - I_exact)
        err_gl[i] = abs(np.dot(w_gl, x_gl ** k) - I_exact)
        err_cc[i] = abs(np.dot(w_cc, x_cc ** k) - I_exact)

    # Clamp errors at the machine-epsilon floor for log plotting
    eps_floor = 1e-17
    err_nc = np.maximum(err_nc, eps_floor)
    err_gl = np.maximum(err_gl, eps_floor)
    err_cc = np.maximum(err_cc, eps_floor)

    # Quick sanity prints
    print(f"n + 1 = {n + 1} quadrature points; testing k = 0 .. {2 * n}")
    print(f"max NC error for k <= n: {err_nc[k_vals <= n].max():.2e}")
    print(f"max CC error for k <= n: {err_cc[k_vals <= n].max():.2e}")
    print(f"max GL error for entire range: {err_gl.max():.2e}")
    print(f"max NC error for k > n: {err_nc[k_vals > n].max():.2e}")
    print(f"max CC error for k > n: {err_cc[k_vals > n].max():.2e}")

    # -------------------------------------------------------------------------
    # Create figure
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 4.8))

    ax.semilogy(k_vals, err_nc, 'o-', color=CORAL, markersize=4,
                linewidth=1.1, label=r'Newton$\mathrm{-}$Cotes')
    ax.semilogy(k_vals, err_gl, 's-', color=NAVY, markersize=4,
                linewidth=1.1, label=r'Gauss$\mathrm{-}$Legendre')
    ax.semilogy(k_vals, err_cc, 'D-', color=TEAL, markersize=4,
                linewidth=1.1, label=r'Clenshaw$\mathrm{-}$Curtis')

    # Vertical guide at k = n: NC/CC exactness boundary.  The label sits in
    # the empty band between the NC errors above and the GL/CC floor below.
    ax.axvline(x=n, color='gray', linestyle='--', linewidth=0.8, zorder=0)
    ax.text(n + 0.6, 1e-13, r'NC/CC exactness ($k = n$)',
            rotation=90, va='center', ha='left',
            fontsize=9, color='gray')

    # Vertical guide at k = 2n + 1: Gauss exactness boundary.  Placed to the
    # LEFT of the dashed line so the rotated text does not overrun the
    # right-hand axis.
    ax.axvline(x=2 * n + 1, color='gray', linestyle='--',
               linewidth=0.8, zorder=0)
    ax.text(2 * n + 1 - 0.6, 1e-13, r'Gauss exactness ($k = 2n + 1$)',
            rotation=90, va='center', ha='right',
            fontsize=9, color='gray')

    ax.set_xlabel(r'Monomial degree $k$')
    ax.set_ylabel(r'Absolute error $|E_n(x^k)|$')
    ax.set_title(
        r'Polynomial exactness of three quadrature rules '
        r'($n + 1 = 33$ points)')
    # The upper-left corner of the panel is empty (NC errors sit near
    # 1e-10 there, GL/CC are at the floor), so the legend fits cleanly.
    ax.legend(loc='upper left', frameon=True, fancybox=False,
              edgecolor='gray', framealpha=0.9)
    ax.grid(True, which='both', alpha=0.3, linewidth=0.5)

    ax.set_xlim(-1, 2 * n + 2)
    ax.set_ylim(1e-18, 1e1)
    ax.set_xticks(np.arange(0, 2 * n + 1, 8))

    plt.tight_layout()

    # Save output
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'polynomial_exactness.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'polynomial_exactness.png', bbox_inches='tight')
    print(f"Figure saved to {OUTPUT_DIR / 'polynomial_exactness.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
