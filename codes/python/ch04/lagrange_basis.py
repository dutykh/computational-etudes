#!/usr/bin/env python3
"""
lagrange_basis.py

Visualizes Lagrange basis polynomials for equispaced vs Chebyshev nodes.

The Lagrange basis polynomials are:
    L_k(x) = ∏_{j≠k} (x - x_j) / (x_k - x_j)

These polynomials satisfy L_k(x_j) = δ_{kj} (Kronecker delta), making them
the "cardinal functions" for interpolation: the interpolant is simply
    p_N(x) = Σ_{k=0}^{N} f_k L_k(x)

For equispaced nodes, the L_k(x) develop large oscillations near the boundaries,
contributing to the Runge phenomenon. For Chebyshev nodes, the basis functions
remain well-behaved across the entire interval.

The maximum magnitude of the sum Σ|L_k(x)| is the Lebesgue function, whose
maximum is the Lebesgue constant Λ_N.

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
N = 10  # Number of intervals (N+1 nodes)
N_FINE = 1000  # Number of points for smooth plotting

# Book color scheme
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'
PURPLE = '#8E44AD'
ORANGE = '#E67E22'

# Output path
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch04' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'lagrange_basis.pdf'


# -----------------------------------------------------------------------------
# Node generation
# -----------------------------------------------------------------------------
def equispaced_nodes(N):
    """Equispaced nodes on [-1, 1]."""
    return np.linspace(-1, 1, N + 1)


def chebyshev_nodes(N):
    """Chebyshev-Gauss-Lobatto nodes on [-1, 1]."""
    j = np.arange(N + 1)
    return np.cos(j * np.pi / N)


# -----------------------------------------------------------------------------
# Lagrange basis computation
# -----------------------------------------------------------------------------
def lagrange_basis(x_nodes, k, x_eval):
    """
    Compute the k-th Lagrange basis polynomial at points x_eval.

    L_k(x) = ∏_{j≠k} (x - x_j) / (x_k - x_j)

    Parameters
    ----------
    x_nodes : array_like
        Interpolation nodes
    k : int
        Index of the basis polynomial
    x_eval : array_like
        Points at which to evaluate

    Returns
    -------
    L_k : ndarray
        Values of L_k at x_eval
    """
    n = len(x_nodes)
    x_eval = np.atleast_1d(x_eval)
    L_k = np.ones_like(x_eval, dtype=float)

    for j in range(n):
        if j != k:
            L_k *= (x_eval - x_nodes[j]) / (x_nodes[k] - x_nodes[j])

    return L_k


# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
def main():
    # Fine grid for plotting
    x_fine = np.linspace(-1, 1, N_FINE)

    # Generate nodes
    x_equi = equispaced_nodes(N)
    x_cheb = chebyshev_nodes(N)

    # Create figure with two panels
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # Indices of basis functions to plot (endpoints and near-endpoints)
    indices_to_plot = [0, 1, N//2, N-1, N]
    colors = [CORAL, ORANGE, TEAL, SKY, PURPLE]

    # Left panel: Equispaced basis functions
    for k, color in zip(indices_to_plot, colors):
        L_k = lagrange_basis(x_equi, k, x_fine)
        ax1.plot(x_fine, L_k, color=color, linewidth=1.2,
                 label=f'$L_{{{k}}}(x)$')

    # Mark nodes
    ax1.plot(x_equi, np.zeros_like(x_equi), 'o', color=NAVY, markersize=5)

    ax1.axhline(y=0, color='#888888', linewidth=0.5, linestyle='-')
    ax1.axhline(y=1, color='#888888', linewidth=0.5, linestyle='--', alpha=0.5)

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$L_k(x)$')
    ax1.set_title(f'Equispaced Nodes ($N = {N}$)', fontsize=11)
    ax1.set_xlim(-1, 1)
    ax1.set_ylim(-2.5, 2.5)
    ax1.legend(loc='upper left', fontsize=9, ncol=2)

    # Add annotation about oscillations
    ax1.annotate('Large oscillations\nnear boundaries',
                 xy=(-0.85, -1.8), xytext=(-0.5, -2.2),
                 fontsize=9, color='gray',
                 arrowprops=dict(arrowstyle='->', color='gray', lw=0.8))

    # Right panel: Chebyshev basis functions
    for k, color in zip(indices_to_plot, colors):
        L_k = lagrange_basis(x_cheb, k, x_fine)
        ax2.plot(x_fine, L_k, color=color, linewidth=1.2,
                 label=f'$L_{{{k}}}(x)$')

    # Mark nodes
    ax2.plot(x_cheb, np.zeros_like(x_cheb), 'o', color=NAVY, markersize=5)

    ax2.axhline(y=0, color='#888888', linewidth=0.5, linestyle='-')
    ax2.axhline(y=1, color='#888888', linewidth=0.5, linestyle='--', alpha=0.5)

    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r'$L_k(x)$')
    ax2.set_title(f'Chebyshev Nodes ($N = {N}$)', fontsize=11)
    ax2.set_xlim(-1, 1)
    ax2.set_ylim(-2.5, 2.5)
    ax2.legend(loc='upper left', fontsize=9, ncol=2)

    # Add annotation about stability
    ax2.annotate('Bounded oscillations',
                 xy=(0.6, 0.8), xytext=(0.3, 1.5),
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

    # Print maximum values of basis functions
    print(f"\nMaximum |L_k(x)| for N = {N}:")
    print("-" * 50)
    print(f"{'k':>4} {'Equispaced':>15} {'Chebyshev':>15}")
    print("-" * 50)
    for k in range(N + 1):
        L_equi = lagrange_basis(x_equi, k, x_fine)
        L_cheb = lagrange_basis(x_cheb, k, x_fine)
        print(f"{k:>4d} {np.max(np.abs(L_equi)):>15.4f} {np.max(np.abs(L_cheb)):>15.4f}")
    print("-" * 50)

    # Compute Lebesgue functions
    lebesgue_equi = sum(np.abs(lagrange_basis(x_equi, k, x_fine))
                        for k in range(N + 1))
    lebesgue_cheb = sum(np.abs(lagrange_basis(x_cheb, k, x_fine))
                        for k in range(N + 1))

    print(f"\nLebesgue constants (Λ_N = max Σ|L_k|):")
    print(f"  Equispaced: Λ_{N} = {np.max(lebesgue_equi):.2f}")
    print(f"  Chebyshev:  Λ_{N} = {np.max(lebesgue_cheb):.2f}")
    print(f"  Ratio: {np.max(lebesgue_equi) / np.max(lebesgue_cheb):.1f}x")

    plt.close(fig)


if __name__ == '__main__':
    main()
