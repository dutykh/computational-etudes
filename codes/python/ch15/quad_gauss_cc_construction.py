#!/usr/bin/env python3
"""
quad_gauss_cc_construction.py
Chapter 15: Quadrature in Spectral Methods

Computational Etude 15.4: Constructing Gauss and Clenshaw--Curtis Rules.

Implements two fundamental quadrature construction algorithms from scratch
and validates them against known results.

Gauss--Legendre via the Golub--Welsch Algorithm
------------------------------------------------
The n-point Gauss--Legendre nodes are eigenvalues of the n x n symmetric
tridiagonal Jacobi matrix J with:

    J[j, j]   = 0                          (diagonal)
    J[j, j+1] = beta_j = j / sqrt(4j^2 - 1)  (off-diagonal, j = 1, ..., n-1)

The weights are w_j = 2 * (v_j[0])^2, where v_j is the eigenvector
corresponding to the j-th eigenvalue.

Clenshaw--Curtis via FFT
-------------------------
The n+1 Clenshaw--Curtis nodes are x_j = cos(j pi / n), j = 0, ..., n.
The weights are computed by:

    1. Form the Chebyshev moments c_k = 2/(1-k^2) for even k, 0 for odd k.
    2. Construct a length-2n vector and compute its real IFFT.
    3. Extract the weights (with endpoint halving).

Both methods are validated by integrating e^x (analytic: e - 1/e) and
1/(1+x^2) (analytic: pi/2) and comparing with NumPy's reference
Gauss--Legendre implementation.

Generates Figure 15.4: gauss_cc_construction.pdf (two-panel figure)

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


def golub_welsch(n):
    """
    Compute n-point Gauss--Legendre nodes and weights via the Golub--Welsch
    algorithm.

    The Jacobi matrix for Legendre polynomials has zero diagonal and
    off-diagonal entries beta_j = j / sqrt(4j^2 - 1) for j = 1, ..., n-1.

    Parameters
    ----------
    n : int
        Number of quadrature points.

    Returns
    -------
    x : ndarray of shape (n,)
        Gauss--Legendre nodes (sorted).
    w : ndarray of shape (n,)
        Corresponding weights.
    """
    if n == 1:
        return np.array([0.0]), np.array([2.0])

    # Off-diagonal entries of the Jacobi matrix
    j = np.arange(1, n, dtype=float)
    beta = j / np.sqrt(4.0 * j**2 - 1.0)

    # Build the symmetric tridiagonal matrix
    J = np.diag(beta, 1) + np.diag(beta, -1)

    # Eigendecomposition
    nodes, V = np.linalg.eigh(J)

    # Weights from the first component of each eigenvector
    weights = 2.0 * V[0, :]**2

    # Sort by node position
    idx = np.argsort(nodes)
    return nodes[idx], weights[idx]


def clenshaw_curtis_fft(n):
    """
    Compute (n+1)-point Clenshaw--Curtis nodes and weights via DCT-I / FFT.

    The weights are derived from the Chebyshev moments
    mu_k = int_{-1}^{1} T_k(x) dx = 2/(1-k^2) for even k, 0 for odd k,
    using the DCT-I relation implemented through an FFT of length 2n.

    Parameters
    ----------
    n : int
        Polynomial degree (n+1 nodes).

    Returns
    -------
    x : ndarray of shape (n+1,)
        Clenshaw--Curtis nodes (x_j = cos(j pi / n), j = 0, ..., n).
    w : ndarray of shape (n+1,)
        Corresponding weights.
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
    """Generate the two-panel validation figure."""
    # Test functions and their exact integrals over [-1, 1]
    f1 = lambda x: np.exp(x)
    I1_exact = np.e - 1.0 / np.e  # e - e^{-1}

    f2 = lambda x: 1.0 / (1.0 + x**2)
    I2_exact = np.pi / 2.0  # 2 * arctan(1) = pi/2

    n_values = np.arange(2, 52, 2)

    err_gw_f1 = np.zeros(len(n_values))
    err_cc_f1 = np.zeros(len(n_values))
    err_ref_f1 = np.zeros(len(n_values))

    err_gw_f2 = np.zeros(len(n_values))
    err_cc_f2 = np.zeros(len(n_values))
    err_ref_f2 = np.zeros(len(n_values))

    for i, n in enumerate(n_values):
        # Golub--Welsch Gauss--Legendre (n+1 points)
        x_gw, w_gw = golub_welsch(n + 1)
        err_gw_f1[i] = abs(np.dot(w_gw, f1(x_gw)) - I1_exact)
        err_gw_f2[i] = abs(np.dot(w_gw, f2(x_gw)) - I2_exact)

        # Clenshaw--Curtis (n+1 points)
        x_cc, w_cc = clenshaw_curtis_fft(n)
        err_cc_f1[i] = abs(np.dot(w_cc, f1(x_cc)) - I1_exact)
        err_cc_f2[i] = abs(np.dot(w_cc, f2(x_cc)) - I2_exact)

        # NumPy reference Gauss--Legendre (n+1 points)
        x_ref, w_ref = np.polynomial.legendre.leggauss(n + 1)
        err_ref_f1[i] = abs(np.dot(w_ref, f1(x_ref)) - I1_exact)
        err_ref_f2[i] = abs(np.dot(w_ref, f2(x_ref)) - I2_exact)

    # Clamp zeros
    floor = 1e-16
    err_gw_f1 = np.maximum(err_gw_f1, floor)
    err_cc_f1 = np.maximum(err_cc_f1, floor)
    err_ref_f1 = np.maximum(err_ref_f1, floor)
    err_gw_f2 = np.maximum(err_gw_f2, floor)
    err_cc_f2 = np.maximum(err_cc_f2, floor)
    err_ref_f2 = np.maximum(err_ref_f2, floor)

    # Create two-panel figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # --- Panel (a): e^x ---
    ax1.semilogy(n_values, err_gw_f1, 'o-', color=NAVY, markersize=4,
                 linewidth=1.2, label='Golub--Welsch (Gauss)')
    ax1.semilogy(n_values, err_cc_f1, 'D-', color=TEAL, markersize=4,
                 linewidth=1.2, label='CC (FFT)')
    ax1.semilogy(n_values, err_ref_f1, 'x--', color=CORAL, markersize=5,
                 linewidth=1.0, label='NumPy reference (Gauss)')

    ax1.set_xlabel(r'Number of points $n$')
    ax1.set_ylabel(r'Absolute error $|I - I_n|$')
    ax1.set_title(r'(a) $f(x) = e^x$')
    ax1.legend(loc='best', frameon=True, fancybox=False, edgecolor='gray',
               fontsize=9)
    ax1.grid(True, alpha=0.3, linewidth=0.5)
    ax1.set_xlim(n_values[0] - 1, n_values[-1] + 1)

    # --- Panel (b): 1/(1+x^2) ---
    ax2.semilogy(n_values, err_gw_f2, 'o-', color=NAVY, markersize=4,
                 linewidth=1.2, label='Golub--Welsch (Gauss)')
    ax2.semilogy(n_values, err_cc_f2, 'D-', color=TEAL, markersize=4,
                 linewidth=1.2, label='CC (FFT)')
    ax2.semilogy(n_values, err_ref_f2, 'x--', color=CORAL, markersize=5,
                 linewidth=1.0, label='NumPy reference (Gauss)')

    ax2.set_xlabel(r'Number of points $n$')
    ax2.set_ylabel(r'Absolute error $|I - I_n|$')
    ax2.set_title(r'(b) $f(x) = 1/(1 + x^2)$')
    ax2.legend(loc='best', frameon=True, fancybox=False, edgecolor='gray',
               fontsize=9)
    ax2.grid(True, alpha=0.3, linewidth=0.5)
    ax2.set_xlim(n_values[0] - 1, n_values[-1] + 1)

    plt.tight_layout()

    # Save output
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'gauss_cc_construction.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'gauss_cc_construction.png', bbox_inches='tight')
    print(f"Figure saved to {OUTPUT_DIR / 'gauss_cc_construction.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
