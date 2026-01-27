#!/usr/bin/env python3
"""
fourier_decay.py

Demonstration of the decay hierarchy for Fourier coefficients, illustrating
Theorem 1: the connection between function smoothness and spectral decay.

Three test functions on [0, 2pi]:
    1. |sin(x)|^3    - Finite regularity (C^2 but not C^3)
                       Decay rate: O(k^{-4})
    2. 1/(1+sin^2(x/2)) - Analytic in a strip in the complex plane
                       Decay rate: O(c^{-|k|}) geometric
    3. exp(sin(x))   - Entire function (analytic everywhere)
                       Decay rate: faster than any exponential

The plot shows |f_hat_k| vs k on a log scale, revealing:
    - Algebraic decay appears as a straight line on log-log scale
    - Geometric decay appears as a straight line on semi-log scale
    - Super-geometric decay curves downward on semi-log scale

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

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'
PURPLE = '#9B59B6'

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch06' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'decay_hierarchy.pdf'


# -----------------------------------------------------------------------------
# Test functions
# -----------------------------------------------------------------------------
def f_finite_regularity(x):
    """
    |sin(x)|^3 - has bounded third derivative, fourth derivative in BV.
    Decay rate: O(k^{-4})
    """
    return np.abs(np.sin(x))**3


def f_analytic_strip(x):
    """
    1/(1 + sin^2(x/2)) - analytic in a strip |Im(z)| < a for some a > 0.
    Has poles where sin^2(z/2) = -1, i.e., at complex values.
    Decay rate: O(c^{-|k|}) geometric
    """
    return 1.0 / (1.0 + np.sin(x / 2.0)**2)


def f_entire(x):
    """
    exp(sin(x)) - entire function (analytic throughout the complex plane).
    Decay rate: faster than any exponential (super-geometric)
    """
    return np.exp(np.sin(x))


# -----------------------------------------------------------------------------
# Fourier coefficient computation
# -----------------------------------------------------------------------------
def compute_fourier_coefficients(f, N):
    """
    Compute the Fourier coefficients of a 2pi-periodic function f
    sampled at N equispaced points.

    Parameters
    ----------
    f : callable
        The function to analyze
    N : int
        Number of sample points

    Returns
    -------
    k : ndarray
        Wavenumbers (0, 1, 2, ..., N/2)
    f_hat : ndarray
        Absolute values of Fourier coefficients |f_hat_k|
    """
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    f_vals = f(x)

    # Compute FFT and normalize
    f_hat = np.fft.fft(f_vals) / N

    # Take absolute values and keep only positive frequencies
    f_hat_abs = np.abs(f_hat[:N // 2 + 1])

    # Wavenumbers
    k = np.arange(N // 2 + 1)

    return k, f_hat_abs


# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
def main():
    N = 128  # Use enough points to see the decay clearly

    # Compute Fourier coefficients for each function
    k1, f1_hat = compute_fourier_coefficients(f_finite_regularity, N)
    k2, f2_hat = compute_fourier_coefficients(f_analytic_strip, N)
    k3, f3_hat = compute_fourier_coefficients(f_entire, N)

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 5.5))

    # Plot coefficients (skip k=0 for cleaner visualization on log scale)
    ax.semilogy(k1[1:], f1_hat[1:], 'o-', color=TEAL, linewidth=1.5,
                markersize=4, label=r'$|\sin(x)|^3$ (finite regularity)')
    ax.semilogy(k2[1:], f2_hat[1:], 's-', color=CORAL, linewidth=1.5,
                markersize=4, label=r'$1/(1+\sin^2(x/2))$ (analytic in strip)')
    ax.semilogy(k3[1:], f3_hat[1:], 'D-', color=NAVY, linewidth=1.5,
                markersize=4, label=r'$\exp(\sin(x))$ (entire)')

    # Add theoretical decay rate references
    k_ref = np.arange(5, 60)

    # O(k^{-4}) reference for finite regularity
    C1 = f1_hat[10] * 10**4
    ax.semilogy(k_ref, C1 / k_ref**4, '--', color=TEAL, alpha=0.5, linewidth=1)
    ax.text(55, C1 / 55**4 * 2, r'$O(k^{-4})$', fontsize=10, color=TEAL)

    # Geometric decay reference for analytic strip
    # Find decay rate by fitting
    k_fit = np.arange(10, 40)
    idx_fit = k_fit
    log_f2 = np.log(f2_hat[idx_fit] + 1e-20)
    # Linear fit: log(f_hat) = a - b*k
    coeffs = np.polyfit(k_fit, log_f2, 1)
    decay_rate = -coeffs[0]
    ax.semilogy(k_ref, np.exp(coeffs[1]) * np.exp(-decay_rate * k_ref),
                '--', color=CORAL, alpha=0.5, linewidth=1)
    ax.text(55, np.exp(coeffs[1]) * np.exp(-decay_rate * 55) * 3,
            r'$O(e^{-%.2f k})$' % decay_rate, fontsize=10, color=CORAL)

    # Machine epsilon line
    ax.axhline(y=1e-16, color='gray', linestyle=':', linewidth=1)
    ax.text(62, 2e-16, 'Machine precision', fontsize=9, color='gray',
            va='bottom', ha='right')

    # Labels and title
    ax.set_xlabel(r'Wavenumber $k$', fontsize=11)
    ax.set_ylabel(r'Fourier coefficient magnitude $|\hat{f}_k|$', fontsize=11)
    ax.set_title('The Decay Hierarchy: Smoothness Determines Spectral Decay',
                 fontsize=12)

    # Axis settings
    ax.set_xlim(0, 65)
    ax.set_ylim(1e-17, 1e1)
    ax.legend(loc='upper right', fontsize=10)

    # Clean styling
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#444444')
    ax.spines['bottom'].set_color('#444444')
    ax.tick_params(colors='#444444')
    ax.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
    ax.xaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
    ax.set_axisbelow(True)

    plt.tight_layout()

    # Save figure
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_FILE, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(OUTPUT_FILE.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {OUTPUT_FILE.resolve()}')

    # Print coefficient values for verification
    print("\nFourier Coefficient Magnitudes:")
    print("=" * 60)
    print(f"{'k':>4}  {'|sin(x)|^3':>14}  {'1/(1+sin^2)':>14}  {'exp(sin)':>14}")
    print("-" * 60)
    for k in [1, 2, 5, 10, 20, 30, 40, 50]:
        if k < len(f1_hat):
            print(f"{k:4d}  {f1_hat[k]:14.6e}  {f2_hat[k]:14.6e}  {f3_hat[k]:14.6e}")
    print("=" * 60)

    plt.close(fig)


if __name__ == '__main__':
    main()
