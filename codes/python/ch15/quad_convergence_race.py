#!/usr/bin/env python3
"""
quad_convergence_race.py
Chapter 15: Quadrature in Spectral Methods

Computational Etude 15.5: The Convergence Race.

Reproduces the essential content of Figure 2 from Trefethen (2008),
comparing the convergence of Gauss--Legendre and Clenshaw--Curtis
quadrature for six functions of varying regularity on [-1, 1]:

    (a) x^{20}        -- polynomial, finite exactness threshold
    (b) e^x            -- entire function, supergeometric convergence
    (c) e^{-x^2}       -- entire, supergeometric convergence
    (d) 1/(1+16x^2)    -- analytic with poles at x = +- i/4, geometric rate
    (e) e^{-1/x^2}     -- C^infinity but not analytic at x = 0
    (f) |x|^3          -- C^2 but not C^3, algebraic convergence

The six-panel (3 x 2) figure shows the absolute quadrature error
|I - I_n| vs the number of nodes n on a semilogarithmic scale for
n = 1, 2, ..., 40.

For smooth analytic functions, both rules converge at similar rates
(geometric or supergeometric). For functions with limited regularity,
both degrade to algebraic convergence. The key insight of Trefethen
(2008) is that Clenshaw--Curtis is competitive with Gauss for most
practical integrands.

Reference: L. N. Trefethen, "Is Gauss Quadrature Better than
           Clenshaw--Curtis?", SIAM Review 50(1):67--87, 2008, Figure 2.

Generates Figure 15.5: convergence_race.pdf

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


def clenshaw_curtis_weights(n):
    """
    Compute (n+1)-point Clenshaw--Curtis nodes and weights via DCT-I / FFT.

    The weights are derived from the Chebyshev moments
    mu_k = int_{-1}^{1} T_k(x) dx = 2/(1-k^2) for even k, 0 for odd k,
    using the DCT-I relation implemented through an FFT of length 2n.

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


# -----------------------------------------------------------------------------
# Test functions and their exact integrals
# -----------------------------------------------------------------------------

def f_poly(x):
    """f(x) = x^20."""
    return x**20


def f_exp(x):
    """f(x) = e^x."""
    return np.exp(x)


def f_gauss(x):
    """f(x) = e^{-x^2}."""
    return np.exp(-x**2)


def f_runge(x):
    """f(x) = 1/(1 + 16x^2)."""
    return 1.0 / (1.0 + 16.0 * x**2)


def f_flat(x):
    """
    f(x) = e^{-1/x^2}, C^infinity but not analytic at x = 0.
    Defined as 0 at x = 0.
    """
    x = np.asarray(x, dtype=float)
    result = np.zeros_like(x)
    mask = x != 0
    result[mask] = np.exp(-1.0 / x[mask]**2)
    return result


def f_abs3(x):
    """f(x) = |x|^3."""
    return np.abs(x)**3


# Exact integrals over [-1, 1]
I_poly = 2.0 / 21.0                           # integral x^20 dx
I_exp = np.e - 1.0 / np.e                     # e - e^{-1}
I_gauss = np.sqrt(np.pi) * erf(1.0)           # sqrt(pi) * erf(1)
I_runge = np.arctan(4.0) / 2.0                # (1/4)*2*arctan(4) = arctan(4)/2
I_flat, _ = scipy_quad(f_flat, -1, 1,         # numerical high-accuracy reference
                        limit=200, epsabs=1e-15, epsrel=1e-15)
I_abs3 = 0.5                                  # 2 * integral_0^1 x^3 dx = 2*(1/4) = 1/2

# Collect into lists for iteration
functions = [f_poly, f_exp, f_gauss, f_runge, f_flat, f_abs3]
integrals = [I_poly, I_exp, I_gauss, I_runge, I_flat, I_abs3]
titles = [
    r'$x^{20}$',
    r'$e^x$',
    r'$e^{-x^2}$',
    r'$1/(1 + 16x^2)$',
    r'$e^{-1/x^2}$',
    r'$|x|^3$',
]


def main():
    """Generate the six-panel convergence race figure."""
    n_max = 40
    n_values = np.arange(1, n_max + 1)

    # Preallocate error arrays: [6 functions] x [n_max values]
    err_gauss = np.zeros((6, n_max))
    err_cc = np.zeros((6, n_max))

    for idx in range(6):
        f = functions[idx]
        I_exact = integrals[idx]

        for i, n in enumerate(n_values):
            # Gauss--Legendre with n+1 points
            x_gl, w_gl = np.polynomial.legendre.leggauss(n + 1)
            I_gl = np.dot(w_gl, f(x_gl))
            err_gauss[idx, i] = abs(I_gl - I_exact)

            # Clenshaw--Curtis with n+1 points
            x_cc, w_cc = clenshaw_curtis_weights(n)
            I_cc = np.dot(w_cc, f(x_cc))
            err_cc[idx, i] = abs(I_cc - I_exact)

    # Clamp zeros for log plot
    floor = 1e-16
    err_gauss = np.maximum(err_gauss, floor)
    err_cc = np.maximum(err_cc, floor)

    # Create 3x2 figure
    fig, axes = plt.subplots(3, 2, figsize=(10, 10))

    for idx, ax in enumerate(axes.flat):
        # Gauss--Legendre: filled circles
        ax.semilogy(n_values, err_gauss[idx], 'o', color=NAVY,
                     markersize=3.5, markerfacecolor=NAVY,
                     markeredgecolor=NAVY, markeredgewidth=0.5,
                     label=r'Gauss$\mathrm{-}$Legendre')

        # Clenshaw--Curtis: open circles
        ax.semilogy(n_values, err_cc[idx], 'o', color=TEAL,
                     markersize=3.5, markerfacecolor='white',
                     markeredgecolor=TEAL, markeredgewidth=0.8,
                     label=r'Clenshaw$\mathrm{-}$Curtis')

        ax.set_title(titles[idx], fontsize=11)
        ax.set_xlim(0, n_max + 1)
        ax.set_ylim(bottom=floor * 0.5)
        ax.grid(True, alpha=0.3, linewidth=0.5)

        # x-label only on bottom row
        if idx >= 4:
            ax.set_xlabel(r'$n$')

        # y-label only on left column
        if idx % 2 == 0:
            ax.set_ylabel(r'$|I - I_n|$')

        # Legend only on the first panel
        if idx == 0:
            ax.legend(loc='upper right', frameon=True, fancybox=False,
                      edgecolor='gray', fontsize=8)

    plt.tight_layout()

    # Save output
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'convergence_race.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'convergence_race.png', bbox_inches='tight')
    print(f"Figure saved to {OUTPUT_DIR / 'convergence_race.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
