#!/usr/bin/env python3
"""
chebyshev_success.py

Demonstrates successful interpolation of the Runge function using Chebyshev
nodes. In contrast to the disastrous equispaced interpolation, Chebyshev
points provide rapid and uniform convergence.

Chebyshev-Gauss-Lobatto points:
    x_j = cos(jπ/N),  j = 0, 1, ..., N

These points cluster near the boundaries x = ±1, counteracting the growth
of the Lebesgue constant and ensuring convergence for smooth functions.

This script shows interpolation for N = 6, 10, 14 (Chebyshev nodes)
demonstrating excellent accuracy across the entire interval.

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Études: A Spectral Approach"
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
rcParams['text.usetex'] = False  # Set to True if LaTeX is available
rcParams['axes.linewidth'] = 0.8
rcParams['xtick.major.width'] = 0.8
rcParams['ytick.major.width'] = 0.8

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
N_FINE = 1000  # Number of points for smooth plotting

# Book color scheme
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'

# Output path
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch04' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'chebyshev_success.pdf'


# -----------------------------------------------------------------------------
# The Runge function
# -----------------------------------------------------------------------------
def runge(x):
    """
    The Runge function: f(x) = 1 / (1 + 25x²)
    """
    return 1.0 / (1.0 + 25.0 * x**2)


# -----------------------------------------------------------------------------
# Chebyshev points
# -----------------------------------------------------------------------------
def chebyshev_nodes(N):
    """
    Chebyshev-Gauss-Lobatto points on [-1, 1].

    x_j = cos(jπ/N),  j = 0, 1, ..., N

    These are the extrema of the Chebyshev polynomial T_N(x), plus the
    endpoints. They include x = ±1 and cluster more densely near the
    boundaries.

    Parameters
    ----------
    N : int
        Number of intervals (N+1 points)

    Returns
    -------
    x : ndarray
        Chebyshev nodes in decreasing order (from 1 to -1)
    """
    j = np.arange(N + 1)
    return np.cos(j * np.pi / N)


# -----------------------------------------------------------------------------
# Lagrange interpolation
# -----------------------------------------------------------------------------
def lagrange_interpolate(x_nodes, f_nodes, x_eval):
    """
    Evaluate the Lagrange interpolating polynomial at points x_eval.

    Parameters
    ----------
    x_nodes : array_like
        Interpolation nodes
    f_nodes : array_like
        Function values at nodes
    x_eval : array_like
        Points at which to evaluate the interpolant

    Returns
    -------
    p_eval : ndarray
        Interpolant values at x_eval
    """
    n = len(x_nodes)
    x_eval = np.atleast_1d(x_eval)
    p_eval = np.zeros_like(x_eval)

    for k in range(n):
        # Compute Lagrange basis polynomial L_k(x)
        L_k = np.ones_like(x_eval)
        for j in range(n):
            if j != k:
                L_k *= (x_eval - x_nodes[j]) / (x_nodes[k] - x_nodes[j])
        p_eval += f_nodes[k] * L_k

    return p_eval


# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
def main():
    # Fine grid for plotting the exact function
    x_fine = np.linspace(-1, 1, N_FINE)
    f_exact = runge(x_fine)

    # Different polynomial degrees to demonstrate
    N_values = [6, 10, 14]
    colors = [SKY, TEAL, CORAL]

    # Create figure with two panels
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    # Left panel: Function and interpolants
    ax1.plot(x_fine, f_exact, color=NAVY, linewidth=2, label='Runge function')

    for N, color in zip(N_values, colors):
        # Chebyshev nodes
        x_nodes = chebyshev_nodes(N)
        f_nodes = runge(x_nodes)

        # Evaluate interpolant
        p_interp = lagrange_interpolate(x_nodes, f_nodes, x_fine)

        # Plot interpolant
        ax1.plot(x_fine, p_interp, color=color, linewidth=1.2,
                 linestyle='--', label=f'$p_{{{N}}}(x)$ (Chebyshev)')

        # Plot nodes
        ax1.plot(x_nodes, f_nodes, 'o', color=color, markersize=5)

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$f(x)$')
    ax1.set_xlim(-1, 1)
    ax1.set_ylim(-0.1, 1.1)
    ax1.set_title('Chebyshev Interpolation: Success', fontsize=11)
    ax1.legend(loc='upper right', frameon=True, fancybox=False,
               edgecolor='none', facecolor='white', framealpha=0.9)

    # Add formula
    textstr = r'$f(x) = \frac{1}{1 + 25x^2}$'
    props = dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='none')
    ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=11,
             verticalalignment='top', bbox=props)

    # Right panel: Interpolation error
    for N, color in zip(N_values, colors):
        x_nodes = chebyshev_nodes(N)
        f_nodes = runge(x_nodes)
        p_interp = lagrange_interpolate(x_nodes, f_nodes, x_fine)
        error = np.abs(runge(x_fine) - p_interp)

        ax2.semilogy(x_fine, error + 1e-16, color=color, linewidth=1.5,
                     label=f'$N = {N}$')

    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r'$|f(x) - p_N(x)|$')
    ax2.set_xlim(-1, 1)
    ax2.set_ylim(1e-10, 1)
    ax2.set_title('Interpolation Error (log scale)', fontsize=11)
    ax2.legend(loc='upper right', frameon=True, fancybox=False,
               edgecolor='none', facecolor='white', framealpha=0.9)

    # Add annotation
    ax2.annotate('Error decreases\nwith N',
                 xy=(0, 1e-5), xytext=(0.4, 1e-3),
                 fontsize=9, color='gray',
                 arrowprops=dict(arrowstyle='->', color='gray', lw=0.8))

    # Clean styling for both panels
    for ax in [ax1, ax2]:
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

    # Print error statistics
    print("\nMaximum interpolation error (Chebyshev vs Equispaced):")
    print("-" * 60)
    print(f"{'N':>4} {'Chebyshev':>15} {'Equispaced':>15} {'Ratio':>12}")
    print("-" * 60)
    for N in N_values:
        # Chebyshev
        x_cheb = chebyshev_nodes(N)
        f_cheb = runge(x_cheb)
        p_cheb = lagrange_interpolate(x_cheb, f_cheb, x_fine)
        err_cheb = np.max(np.abs(runge(x_fine) - p_cheb))

        # Equispaced
        x_equi = np.linspace(-1, 1, N + 1)
        f_equi = runge(x_equi)
        p_equi = lagrange_interpolate(x_equi, f_equi, x_fine)
        err_equi = np.max(np.abs(runge(x_fine) - p_equi))

        ratio = err_equi / err_cheb if err_cheb > 0 else np.inf
        print(f"{N:>4d} {err_cheb:>15.2e} {err_equi:>15.2e} {ratio:>12.1f}x")
    print("-" * 60)
    print("Chebyshev interpolation is dramatically more accurate!")

    plt.close(fig)


if __name__ == '__main__':
    main()
