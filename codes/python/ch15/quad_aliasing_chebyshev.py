#!/usr/bin/env python3
"""
quad_aliasing_chebyshev.py
Chapter 15: Quadrature in Spectral Methods

Computational Etude 15.6: Aliasing and Chebyshev Quadrature.

Illustrates the aliasing phenomenon in Clenshaw--Curtis quadrature by
examining the quadrature error E_n(T_k) = |I(T_k) - I_n(T_k)| for
individual Chebyshev polynomials T_k, and its impact on the integration
of |x|^3.

Left panel: For n = 30, the CC quadrature error |E_n(T_k)| is plotted
against even k from 0 to 200. Three regimes are visible:
    (1) k <= n: errors are at machine zero (exact integration),
    (2) n < k <= 2n (shaded): errors are small, O(n^{-2}),
    (3) k > 2n: aliasing errors of O(1) at multiples of 2n.
This structure underlies Trefethen's (2008) explanation of why CC
quadrature is nearly as accurate as Gauss for most integrands.

Right panel: For f(x) = |x|^3, the Chebyshev expansion coefficients
|a_k| (which decay algebraically as k^{-4}) are plotted alongside the
aliased error contributions |a_k * E_n(T_k)|. The product decays so
rapidly that the total aliasing error is negligible.

Reference: L. N. Trefethen, "Is Gauss Quadrature Better than
           Clenshaw--Curtis?", SIAM Review 50(1):67--87, 2008, Figure 3.

Generates Figure 15.6: aliasing_chebyshev.pdf

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


def clenshaw_curtis_weights(n):
    """
    Compute (n+1)-point Clenshaw--Curtis nodes and weights via DCT-I / FFT.

    Parameters
    ----------
    n : int
        Polynomial degree (n+1 nodes: x_j = cos(j pi / n), j = 0,...,n).

    Returns
    -------
    x : ndarray of shape (n+1,)
        Clenshaw--Curtis nodes.
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


def exact_integral_Tk(k):
    """
    Exact integral of T_k(x) over [-1, 1].

    I(T_k) = 2 / (1 - k^2)   for even k,
    I(T_k) = 0                for odd k.
    """
    if k % 2 == 1:
        return 0.0
    return 2.0 / (1.0 - k**2)


def chebyshev_coefficients(f, N):
    """
    Compute Chebyshev expansion coefficients a_k of f(x) using
    the DCT via evaluation at Chebyshev points.

    Parameters
    ----------
    f : callable
        Function to expand.
    N : int
        Number of Chebyshev points (returns N coefficients).

    Returns
    -------
    a : ndarray of shape (N,)
        Chebyshev coefficients a_0, a_1, ..., a_{N-1}.
    """
    j = np.arange(N)
    x = np.cos(np.pi * j / (N - 1))  # Chebyshev extremal points
    fx = f(x)

    # DCT-I via FFT
    v = np.concatenate([fx, fx[-2:1:-1]])
    coeffs = np.real(np.fft.fft(v))[:N] / (N - 1)
    coeffs[0] /= 2
    coeffs[-1] /= 2

    return coeffs


def main():
    """Generate the two-panel aliasing figure."""
    n = 30  # CC quadrature with n+1 = 31 points
    k_max = 200  # Maximum Chebyshev degree to examine

    # Compute CC nodes and weights
    x_cc, w_cc = clenshaw_curtis_weights(n)

    # -------------------------------------------------------------------------
    # Left panel: |E_n(T_k)| for even k = 0, 2, 4, ..., k_max
    # -------------------------------------------------------------------------
    k_even = np.arange(0, k_max + 1, 2)
    errors = np.zeros(len(k_even))

    for i, k in enumerate(k_even):
        # Evaluate T_k at CC nodes using the definition T_k(x) = cos(k*arccos(x))
        Tk_at_nodes = np.cos(k * np.arccos(x_cc))

        # Numerical integral via CC quadrature
        I_n = np.dot(w_cc, Tk_at_nodes)

        # Exact integral
        I_exact = exact_integral_Tk(k)

        errors[i] = abs(I_n - I_exact)

    # -------------------------------------------------------------------------
    # Right panel: Chebyshev coefficients of |x|^3 and aliased errors
    # -------------------------------------------------------------------------
    f_abs3 = lambda x: np.abs(x)**3

    # Compute many Chebyshev coefficients using a large number of points
    N_coeffs = 512
    a_k = chebyshev_coefficients(f_abs3, N_coeffs)

    # Only even coefficients are nonzero for this even function
    k_coeffs = np.arange(0, min(k_max + 1, N_coeffs), 2)
    a_even = np.abs(a_k[k_coeffs])

    # Compute aliased error contributions: |a_k * E_n(T_k)|
    n_common = min(len(k_coeffs), len(k_even))
    k_common = k_even[:n_common]
    a_common = np.zeros(n_common)
    e_common = np.zeros(n_common)

    for i, k in enumerate(k_common):
        if k < N_coeffs:
            a_common[i] = abs(a_k[k])
        # Compute E_n(T_k)
        Tk_at_nodes = np.cos(k * np.arccos(x_cc))
        I_n = np.dot(w_cc, Tk_at_nodes)
        I_exact = exact_integral_Tk(k)
        e_common[i] = abs(I_n - I_exact)

    aliased_product = a_common * e_common

    # -------------------------------------------------------------------------
    # Create figure
    # -------------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # --- Left panel: CC quadrature errors for T_k ---
    floor = 1e-17
    errors_plot = np.maximum(errors, floor)

    ax1.semilogy(k_even, errors_plot, 'o', color=NAVY, markersize=3,
                 markeredgecolor=NAVY, markeredgewidth=0.5)

    # Shade the region n <= k <= 2n
    ax1.axvspan(n, 2 * n, alpha=0.12, color=SKY, label=r'$n \leq k \leq 2n$')

    # Vertical dashed lines at k = n and k = 2n
    ax1.axvline(x=n, color=CORAL, linestyle='--', linewidth=0.8, alpha=0.7)
    ax1.axvline(x=2 * n, color=CORAL, linestyle='--', linewidth=0.8, alpha=0.7)

    # Annotate the regions
    ax1.text(n / 2, 2e-16, 'exact', ha='center', fontsize=9, color=TEAL,
             fontstyle='italic')
    ax1.text(1.5 * n, 1e-3, 'small\nerrors', ha='center', fontsize=9,
             color=NAVY, fontstyle='italic')
    ax1.text(2.5 * n + 20, 0.5, 'aliasing', ha='center', fontsize=9,
             color=CORAL, fontstyle='italic')

    ax1.set_xlabel(r'Even degree $k$')
    ax1.set_ylabel(r'$|E_n(T_k)|$')
    ax1.set_title(r'(a) CC error for Chebyshev polynomials ($n = 30$)')
    ax1.set_xlim(-2, k_max + 5)
    ax1.set_ylim(1e-17, 1e1)
    ax1.legend(loc='lower right', frameon=True, fancybox=False,
               edgecolor='gray', fontsize=9)
    ax1.grid(True, alpha=0.3, linewidth=0.5)

    # --- Right panel: |a_k| and |a_k * E_n(T_k)| for |x|^3 ---
    a_plot = np.maximum(a_even, floor)
    alias_plot = np.maximum(aliased_product[:len(k_coeffs)], floor)

    ax2.semilogy(k_coeffs, a_plot, 'o', color=NAVY, markersize=3,
                 markeredgecolor=NAVY, markeredgewidth=0.5,
                 label=r'$|a_k|$ (Chebyshev coefficients)')

    ax2.semilogy(k_common[:len(k_coeffs)],
                 alias_plot, 's', color=CORAL, markersize=3,
                 markerfacecolor='white', markeredgecolor=CORAL,
                 markeredgewidth=0.7,
                 label=r'$|a_k \cdot E_n(T_k)|$ (aliased error)')

    # Shade the region n <= k <= 2n
    ax2.axvspan(n, 2 * n, alpha=0.12, color=SKY)

    ax2.set_xlabel(r'Even degree $k$')
    ax2.set_ylabel('Magnitude')
    ax2.set_title(r'(b) Aliasing for $f(x) = |x|^3$')
    ax2.set_xlim(-2, k_max + 5)
    ax2.set_ylim(1e-17, 1e1)
    ax2.legend(loc='upper right', frameon=True, fancybox=False,
               edgecolor='gray', fontsize=9)
    ax2.grid(True, alpha=0.3, linewidth=0.5)

    plt.tight_layout()

    # Save output
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'aliasing_chebyshev.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'aliasing_chebyshev.png', bbox_inches='tight')
    print(f"Figure saved to {OUTPUT_DIR / 'aliasing_chebyshev.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
