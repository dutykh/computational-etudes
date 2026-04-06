#!/usr/bin/env python3
"""
quad_convergence_rates.py
Chapter 15: Quadrature in Spectral Methods

Computational Etude 15.10: Asymptotic Convergence Rate Verification.

Verifies the theoretical asymptotic convergence rates for two quadrature
strategies applied to the integral

    I = int_{-inf}^{inf} e^{-x^2} cos(x^3) dx

using the rescaling technique of Trefethen and Weideman (2014).

Left panel: Gauss--Hermite errors plotted against n^{1/2}. Theory
predicts that the error decays as O(exp(-C * n^{1/2})), which
corresponds to "root-exponential" or "sub-geometric" convergence.
On a semilog plot with the horizontal axis scaled as sqrt(n), the
data should fall approximately on a straight line. A linear fit
confirms the predicted rate.

Right panel: Truncated Gauss--Legendre errors plotted against n^{2/3}.
For a Gauss--Legendre rule on the interval [-L * n^{1/3}, L * n^{1/3}]
(where the truncation width grows with n to capture more of the
integrand), theory predicts convergence at rate O(exp(-C * n^{2/3})).
On a semilog plot with the horizontal axis scaled as n^{2/3}, the
data should be approximately linear. This 2/3-power rate is
intermediate between root-exponential and geometric convergence.

Reference: L. N. Trefethen and J. A. C. Weideman, "The exponentially
           convergent trapezoidal rule", SIAM Review 56(3):385--458, 2014.

Generates Figure 15.10: convergence_rates.pdf

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
    """Generate the two-panel convergence rate verification figure."""

    # -------------------------------------------------------------------------
    # Reference value via high-accuracy adaptive quadrature
    # -------------------------------------------------------------------------
    integrand = lambda x: np.exp(-x**2) * np.cos(x**3)
    I_exact, _ = scipy_quad(integrand, -np.inf, np.inf,
                            limit=500, epsabs=1e-15, epsrel=1e-15)
    print(f"Reference integral value: {I_exact:.16e}")

    # -------------------------------------------------------------------------
    # Gauss--Hermite errors for n = 10, 15, ..., 500
    # -------------------------------------------------------------------------
    n_gh = np.arange(10, 501, 5)
    err_gh = np.zeros(len(n_gh))

    for i, n in enumerate(n_gh):
        x_h, w_h = roots_hermite(n)
        I_gh = np.dot(w_h, np.cos(x_h**3))
        err_gh[i] = abs(I_gh - I_exact)

    # -------------------------------------------------------------------------
    # Truncated Gauss--Legendre on [-L * n^{1/3}, L * n^{1/3}]
    # The full integrand e^{-x^2} cos(x^3) is used (no weight separation).
    # -------------------------------------------------------------------------
    L = 3.0  # Scaling parameter
    n_gl = np.arange(10, 501, 5)
    err_gl = np.zeros(len(n_gl))

    for i, n in enumerate(n_gl):
        # Truncation half-width grows as n^{1/3}
        half_width = L * n**(1.0 / 3.0)

        # Gauss--Legendre nodes on [-1, 1], mapped to [-half_width, half_width]
        s, ws = np.polynomial.legendre.leggauss(n)
        x_mapped = half_width * s
        w_mapped = half_width * ws

        # Evaluate the full integrand
        f_vals = np.exp(-x_mapped**2) * np.cos(x_mapped**3)
        I_gl = np.dot(w_mapped, f_vals)
        err_gl[i] = abs(I_gl - I_exact)

    # -------------------------------------------------------------------------
    # Prepare data for plotting: filter out zero errors
    # -------------------------------------------------------------------------
    floor = 1e-16

    # Gauss--Hermite: use errors above floor for fitting
    mask_gh = err_gh > floor
    sqrt_n_gh = np.sqrt(n_gh)
    log_err_gh = np.log10(np.maximum(err_gh, floor))

    # Truncated GL: use errors above floor for fitting
    mask_gl = err_gl > floor
    n23_gl = n_gl**(2.0 / 3.0)
    log_err_gl = np.log10(np.maximum(err_gl, floor))

    # Linear fits on the valid (above-floor) data
    if np.sum(mask_gh) > 2:
        coeffs_gh = np.polyfit(sqrt_n_gh[mask_gh], log_err_gh[mask_gh], 1)
        fit_gh = np.polyval(coeffs_gh, sqrt_n_gh)
        slope_gh = coeffs_gh[0]
    else:
        fit_gh = None
        slope_gh = None

    if np.sum(mask_gl) > 2:
        coeffs_gl = np.polyfit(n23_gl[mask_gl], log_err_gl[mask_gl], 1)
        fit_gl = np.polyval(coeffs_gl, n23_gl)
        slope_gl = coeffs_gl[0]
    else:
        fit_gl = None
        slope_gl = None

    # -------------------------------------------------------------------------
    # Create figure
    # -------------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # --- Left panel: Gauss--Hermite vs sqrt(n) ---
    ax1.semilogy(sqrt_n_gh, np.maximum(err_gh, floor), 'o', color=NAVY,
                 markersize=3, markeredgecolor=NAVY, markeredgewidth=0.5,
                 label=r'Gauss$\mathrm{-}$Hermite error')

    if fit_gh is not None:
        ax1.semilogy(sqrt_n_gh, 10**fit_gh, '--', color=CORAL, linewidth=1.2,
                     label=rf'Linear fit (slope $= {slope_gh:.2f}$)')

    ax1.axhline(y=np.finfo(float).eps, color='gray', linestyle=':',
                linewidth=0.8, alpha=0.7)

    ax1.set_xlabel(r'$n^{1/2}$')
    ax1.set_ylabel(r'Absolute error $|I - I_n|$')
    ax1.set_title(r'(a) Gauss$\mathrm{-}$Hermite: $O(e^{-C\sqrt{n}})$')
    ax1.legend(loc='upper right', frameon=True, fancybox=False,
               edgecolor='gray', fontsize=9)
    ax1.grid(True, alpha=0.3, linewidth=0.5)
    ax1.set_ylim(1e-16, 1e1)

    # --- Right panel: Truncated GL vs n^{2/3} ---
    ax2.semilogy(n23_gl, np.maximum(err_gl, floor), 'o', color=TEAL,
                 markersize=3, markeredgecolor=TEAL, markeredgewidth=0.5,
                 label=r'Truncated GL error')

    if fit_gl is not None:
        ax2.semilogy(n23_gl, 10**fit_gl, '--', color=CORAL, linewidth=1.2,
                     label=rf'Linear fit (slope $= {slope_gl:.2f}$)')

    ax2.axhline(y=np.finfo(float).eps, color='gray', linestyle=':',
                linewidth=0.8, alpha=0.7)

    ax2.set_xlabel(r'$n^{2/3}$')
    ax2.set_ylabel(r'Absolute error $|I - I_n|$')
    ax2.set_title(r'(b) Truncated GL: $O(e^{-C\, n^{2/3}})$')
    ax2.legend(loc='upper right', frameon=True, fancybox=False,
               edgecolor='gray', fontsize=9)
    ax2.grid(True, alpha=0.3, linewidth=0.5)
    ax2.set_ylim(1e-16, 1e1)

    plt.tight_layout()

    # Save output
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'convergence_rates.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'convergence_rates.png', bbox_inches='tight')
    print(f"Figure saved to {OUTPUT_DIR / 'convergence_rates.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
