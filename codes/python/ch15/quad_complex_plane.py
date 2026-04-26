#!/usr/bin/env python3
"""
quad_complex_plane.py
Chapter 15: Quadrature in Spectral Methods

Computational Etude 15.7: Quadrature in the Complex Plane.

Reproduces the essential content of Figure 4 from Trefethen (2008),
comparing Newton--Cotes, Gauss--Legendre, and Clenshaw--Curtis
quadrature rules through the lens of rational approximation in the
complex plane.

The key idea is that a quadrature rule with nodes x_k and weights w_k
defines a rational approximant

    r_n(z) = sum_{k} w_k / (z - x_k)

to the function

    phi(z) = log((z + 1) / (z - 1))

which is the Cauchy integral representation of the constant function 1
over [-1, 1]. The contour plot of |phi(z) - r_n(z)| reveals the
accuracy of each quadrature rule throughout the complex plane.

For Newton--Cotes (equispaced nodes), the approximation is poor near
the endpoints of [-1, 1], reflecting the Runge phenomenon. For
Gauss--Legendre and Clenshaw--Curtis, the approximation is far more
uniform, with CC showing slightly larger errors near the real axis
but comparable accuracy overall.

The figure shows a 1 x 3 panel layout for n = 16 (17 points):
    (a) Newton--Cotes
    (b) Gauss--Legendre
    (c) Clenshaw--Curtis

Reference: L. N. Trefethen, "Is Gauss Quadrature Better than
           Clenshaw--Curtis?", SIAM Review 50(1):67--87, 2008, Figure 4.

Generates Figure 15.7: complex_plane.pdf

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
from matplotlib import ticker
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

    c = np.zeros(n + 1)
    c[0::2] = 2.0 / (1.0 - np.arange(0, n + 1, 2)**2)

    v = np.concatenate([c, c[-2:0:-1]])
    f = np.real(np.fft.fft(v))

    w = f[:n + 1] / n
    w[0] /= 2
    w[-1] /= 2

    return x, w


def newton_cotes_weights(n):
    """
    Compute (n+1)-point Newton--Cotes (equispaced) nodes and weights on [-1, 1]
    via Vandermonde solve.

    Parameters
    ----------
    n : int
        Number of intervals (n+1 nodes).

    Returns
    -------
    x : ndarray of shape (n+1,)
        Equispaced nodes on [-1, 1].
    w : ndarray of shape (n+1,)
        Corresponding weights.
    """
    x = np.linspace(-1, 1, n + 1)

    # Vandermonde system: V^T w = b, where b_k = int_{-1}^{1} x^k dx
    V = np.vander(x, increasing=True).T
    b = np.zeros(n + 1)
    for k in range(n + 1):
        b[k] = (1.0 - (-1.0)**(k + 1)) / (k + 1)

    w = np.linalg.solve(V, b)
    return x, w


def phi(z):
    """
    Compute phi(z) = log((z+1)/(z-1)) using the principal branch.

    This is the Cauchy integral of the constant function 1 over [-1, 1]:
    phi(z) = int_{-1}^{1} 1/(z - t) dt = log((z+1)/(z-1)).
    """
    return np.log((z + 1.0) / (z - 1.0))


def rational_approx(z, x, w):
    """
    Compute the rational approximant r_n(z) = sum_k w_k / (z - x_k).

    Parameters
    ----------
    z : ndarray (complex)
        Evaluation points.
    x : ndarray
        Quadrature nodes.
    w : ndarray
        Quadrature weights.

    Returns
    -------
    r : ndarray (complex)
        Values of the rational approximant at z.
    """
    r = np.zeros_like(z, dtype=complex)
    for k in range(len(x)):
        r += w[k] / (z - x[k])
    return r


def main():
    """Generate the three-panel complex plane figure."""
    n = 16  # 17-point rules

    # Compute nodes and weights for all three rules
    x_nc, w_nc = newton_cotes_weights(n)
    x_gl, w_gl = np.polynomial.legendre.leggauss(n + 1)
    x_cc, w_cc = clenshaw_curtis_weights(n)

    rules = [
        (x_nc, w_nc, r'(a) Newton$\mathrm{-}$Cotes'),
        (x_gl, w_gl, r'(b) Gauss$\mathrm{-}$Legendre'),
        (x_cc, w_cc, r'(c) Clenshaw$\mathrm{-}$Curtis'),
    ]

    # Complex grid
    nx, ny = 400, 300
    x_re = np.linspace(-2.5, 2.5, nx)
    y_im = np.linspace(-2.0, 2.0, ny)
    X, Y = np.meshgrid(x_re, y_im)
    Z = X + 1j * Y

    # Mask points too close to the branch cut [-1, 1] on the real axis
    mask_distance = 0.03
    on_cut = (np.abs(Y) < mask_distance) & (X >= -1.0) & (X <= 1.0)

    # Contour levels: 10^{-13}, 10^{-12}, ..., 10^0 (must be increasing)
    levels = np.logspace(-13, 0, 14)

    # Compute phi(z) once
    phi_vals = phi(Z)

    # NEW (panel d): convergence on a benchmark integral.
    # We use the Runge-like integrand f(x) = 1 / (1 + 16 x^2), whose
    # exact integral on [-1, 1] is (1/2) atan(4) ≈ 0.66291...
    # Newton-Cotes diverges (Runge phenomenon); Gauss-Legendre and
    # Clenshaw-Curtis converge geometrically. Newton-Cotes only "wins"
    # on entire-function integrands, hence the controversy that
    # motivates Trefethen 2008.
    f_bench = lambda x: 1.0 / (1.0 + 16.0 * x ** 2)
    I_exact = 0.5 * np.arctan(4.0)
    ns_conv = np.array([4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128])
    err_nc, err_gl, err_cc = [], [], []
    for nv in ns_conv:
        x_n, w_n = newton_cotes_weights(int(nv))
        x_g, w_g = np.polynomial.legendre.leggauss(int(nv) + 1)
        x_c, w_c = clenshaw_curtis_weights(int(nv))
        err_nc.append(abs(np.sum(w_n * f_bench(x_n)) - I_exact))
        err_gl.append(abs(np.sum(w_g * f_bench(x_g)) - I_exact))
        err_cc.append(abs(np.sum(w_c * f_bench(x_c)) - I_exact))
    err_nc = np.array(err_nc)
    err_gl = np.array(err_gl)
    err_cc = np.array(err_cc)

    # Create 2x2 figure
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 9.0))

    panel_axes = [axes[0, 0], axes[0, 1], axes[1, 0]]
    for ax, (x_q, w_q, title) in zip(panel_axes, rules):
        # Compute rational approximation error
        r_vals = rational_approx(Z, x_q, w_q)
        error = np.abs(phi_vals - r_vals)

        # Mask branch cut region
        error[on_cut] = np.nan

        # Contour plot
        cs = ax.contour(X, Y, error, levels=levels,
                        colors='k', linewidths=0.5,
                        norm=plt.matplotlib.colors.LogNorm(vmin=1e-13, vmax=1))

        # Add contour labels (selected levels)
        label_levels = [1e0, 1e-3, 1e-6, 1e-9, 1e-12]
        fmt = {}
        for lev in levels:
            exp = int(np.round(np.log10(lev)))
            fmt[lev] = r'$10^{%d}$' % exp
        ax.clabel(cs, label_levels, fontsize=7, fmt=fmt, inline=True,
                  inline_spacing=3)

        # Mark the interval [-1, 1]
        ax.plot([-1, 1], [0, 0], '-', color=NAVY, linewidth=2.0, zorder=5)

        # Mark quadrature nodes
        ax.plot(x_q, np.zeros_like(x_q), 'o', color=CORAL, markersize=3,
                markeredgecolor=CORAL, markeredgewidth=0.5, zorder=6)

        ax.set_title(title, fontsize=11)
        ax.set_xlabel(r'$\mathrm{Re}(z)$')
        ax.set_ylabel(r'$\mathrm{Im}(z)$')
        ax.set_xlim(-2.5, 2.5)
        ax.set_ylim(-2.0, 2.0)
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.15, linewidth=0.5)

    # ---- Panel (d): NEW convergence on benchmark integral ----------
    ax = axes[1, 1]
    ax.semilogy(ns_conv, err_nc + 1e-18, '-o', color=CORAL, lw=1.0,
                mfc='white', ms=5,
                label=r'(a) Newton$-$Cotes (Runge: diverges)')
    ax.semilogy(ns_conv, err_gl + 1e-18, '-s', color=NAVY, lw=1.0,
                mfc='white', ms=5,
                label=r'(b) Gauss$-$Legendre (geometric)')
    ax.semilogy(ns_conv, err_cc + 1e-18, '-^', color=TEAL, lw=1.0,
                mfc='white', ms=5,
                label=r'(c) Clenshaw$-$Curtis (geometric)')
    ax.set_xlabel(r'$n$ (rule order, $n+1$ points)')
    ax.set_ylabel(r'$|\hat{I}_n - I|$')
    ax.set_title(r'(d) convergence on $\int_{-1}^{1}\!\frac{dx}{1 + 16 x^2}$',
                 fontsize=11)
    ax.grid(True, which='both', alpha=0.3, linewidth=0.4)
    ax.legend(loc='upper right', fontsize=8, frameon=False)

    plt.tight_layout()

    # Save output
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'complex_plane.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'complex_plane.png', bbox_inches='tight')
    print(f"Figure saved to {OUTPUT_DIR / 'complex_plane.pdf'}")
    print(f"  exact integral = {I_exact:.10f}")
    print(f"  NC at n=64: error = {err_nc[8]:.3e} (Runge: large!)")
    print(f"  GL at n=64: error = {err_gl[8]:.3e}")
    print(f"  CC at n=64: error = {err_cc[8]:.3e}")
    plt.close(fig)


if __name__ == '__main__':
    main()
