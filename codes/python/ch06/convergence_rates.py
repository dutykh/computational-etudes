#!/usr/bin/env python3
"""
convergence_rates.py

Demonstration of spectral differentiation convergence rates for functions
of varying smoothness, illustrating Theorems 3 and 4.

Three test functions on [0, 2pi]:
    1. |sin(x)|^3       - Finite regularity, algebraic convergence O(N^{-3})
    2. 1/(1+sin^2(x/2)) - Analytic in strip, geometric convergence O(c^{-N})
    3. exp(sin(x))      - Entire function, super-geometric convergence

The plot shows differentiation error vs N on a semi-log scale, demonstrating
how function smoothness determines the convergence rate of spectral methods.

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
OUTPUT_FILE = OUTPUT_DIR / 'convergence_rates.pdf'


# -----------------------------------------------------------------------------
# Test functions and their exact derivatives
# -----------------------------------------------------------------------------
def f1(x):
    """Finite regularity: |sin(x)|^3"""
    return np.abs(np.sin(x))**3


def f1_deriv(x):
    """Exact derivative of |sin(x)|^3"""
    # d/dx |sin(x)|^3 = 3|sin(x)|^2 * sign(sin(x)) * cos(x)
    # = 3 sin(x) |sin(x)| cos(x)
    return 3 * np.sin(x) * np.abs(np.sin(x)) * np.cos(x)


def f2(x):
    """Analytic in strip: 1/(1 + sin^2(x/2))"""
    return 1.0 / (1.0 + np.sin(x / 2.0)**2)


def f2_deriv(x):
    """Exact derivative of 1/(1 + sin^2(x/2))"""
    # d/dx 1/(1 + sin^2(x/2)) = -sin(x/2)*cos(x/2) / (1 + sin^2(x/2))^2
    return -np.sin(x / 2.0) * np.cos(x / 2.0) / (1.0 + np.sin(x / 2.0)**2)**2


def f3(x):
    """Entire function: exp(sin(x))"""
    return np.exp(np.sin(x))


def f3_deriv(x):
    """Exact derivative of exp(sin(x))"""
    return np.cos(x) * np.exp(np.sin(x))


# -----------------------------------------------------------------------------
# Spectral differentiation matrix (periodic)
# -----------------------------------------------------------------------------
def spectral_diff_periodic(N):
    """
    Construct the periodic spectral differentiation matrix for N points.

    Parameters
    ----------
    N : int
        Number of grid points (should be even)

    Returns
    -------
    D : ndarray (N, N)
        Differentiation matrix
    x : ndarray (N,)
        Grid points on [0, 2pi)
    """
    h = 2 * np.pi / N
    x = h * np.arange(N)

    # Build differentiation matrix using the cotangent formula
    D = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if i != j:
                D[i, j] = 0.5 * ((-1)**(i - j)) / np.tan((i - j) * np.pi / N)

    return D, x


# -----------------------------------------------------------------------------
# Convergence study
# -----------------------------------------------------------------------------
def compute_errors(N_values):
    """Compute differentiation errors for each function at various N."""
    errors1 = []
    errors2 = []
    errors3 = []

    for N in N_values:
        D, x = spectral_diff_periodic(N)

        # Function 1: |sin(x)|^3
        u1 = f1(x)
        du1_exact = f1_deriv(x)
        du1_spectral = D @ u1
        errors1.append(np.max(np.abs(du1_spectral - du1_exact)))

        # Function 2: 1/(1 + sin^2(x/2))
        u2 = f2(x)
        du2_exact = f2_deriv(x)
        du2_spectral = D @ u2
        errors2.append(np.max(np.abs(du2_spectral - du2_exact)))

        # Function 3: exp(sin(x))
        u3 = f3(x)
        du3_exact = f3_deriv(x)
        du3_spectral = D @ u3
        errors3.append(np.max(np.abs(du3_spectral - du3_exact)))

    return np.array(errors1), np.array(errors2), np.array(errors3)


# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
def main():
    # Grid sizes to test
    N_values = np.array([6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 64])

    # Compute errors
    errors1, errors2, errors3 = compute_errors(N_values)

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 5.5))

    # Plot errors
    ax.semilogy(N_values, errors1, 'o-', color=TEAL, linewidth=1.5,
                markersize=6, label=r'$|\sin(x)|^3$ (finite regularity)')
    ax.semilogy(N_values, errors2, 's-', color=CORAL, linewidth=1.5,
                markersize=6, label=r'$1/(1+\sin^2(x/2))$ (analytic in strip)')
    ax.semilogy(N_values, errors3, 'D-', color=NAVY, linewidth=1.5,
                markersize=6, label=r'$\exp(\sin(x))$ (entire)')

    # Add theoretical reference lines
    N_ref = np.linspace(8, 64, 100)

    # O(N^{-3}) reference for finite regularity
    C1 = errors1[3] * N_values[3]**3
    ax.semilogy(N_ref, C1 / N_ref**3, '--', color=TEAL, alpha=0.5, linewidth=1)
    ax.text(66, C1 / 66**3 * 1.5, r'$O(N^{-3})$', fontsize=10, color=TEAL,
            ha='left', va='bottom')

    # Geometric decay reference for analytic strip
    # Fit the exponential decay rate
    idx_fit = (N_values >= 10) & (N_values <= 40) & (errors2 > 1e-15)
    if np.sum(idx_fit) >= 2:
        coeffs = np.polyfit(N_values[idx_fit], np.log(errors2[idx_fit]), 1)
        decay_rate = -coeffs[0]
        C2 = np.exp(coeffs[1])
        ax.semilogy(N_ref, C2 * np.exp(-decay_rate * N_ref), '--',
                    color=CORAL, alpha=0.5, linewidth=1)
        ax.text(30, C2 * np.exp(-decay_rate * 30) * 3,
                r'$O(e^{-%.2f N})$' % decay_rate, fontsize=10, color=CORAL)

    # Machine epsilon line
    ax.axhline(y=2.2e-16, color='gray', linestyle=':', linewidth=1)
    ax.text(66, 4e-16, 'Machine precision', fontsize=9, color='gray',
            va='bottom', ha='right')

    # Add annotation for spectral convergence
    ax.annotate('Super-geometric\nconvergence',
                xy=(24, errors3[6]), xytext=(35, 1e-8),
                fontsize=10, color=NAVY,
                ha='center',
                arrowprops=dict(arrowstyle='->', color=NAVY, lw=1.5))

    # Labels and title
    ax.set_xlabel(r'Number of grid points $N$', fontsize=11)
    ax.set_ylabel(r'Differentiation error $\|f^\prime - Df\|_\infty$', fontsize=11)
    ax.set_title('Spectral Differentiation Convergence: Smoothness Matters',
                 fontsize=12)

    # Axis settings
    ax.set_xlim(0, 70)
    ax.set_ylim(1e-16, 1e1)
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

    # Print results table
    print("\nConvergence Study Results:")
    print("=" * 65)
    print(f"{'N':>4}  {'|sin(x)|^3':>14}  {'1/(1+sin^2)':>14}  {'exp(sin)':>14}")
    print("-" * 65)
    for i, N in enumerate(N_values):
        print(f"{N:4d}  {errors1[i]:14.6e}  {errors2[i]:14.6e}  {errors3[i]:14.6e}")
    print("=" * 65)

    plt.close(fig)


if __name__ == '__main__':
    main()
