#!/usr/bin/env python3
"""
cheb_cardinal.py

Visualizes Chebyshev cardinal functions (Lagrange interpolation basis) and
shows how the differentiation matrix entries are the derivatives of these
cardinal functions evaluated at the grid points.

This script generates Figure 6.3 for Chapter 6.

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

# Import Chebyshev matrix function
import sys
sys.path.insert(0, str(Path(__file__).parent))
from cheb_matrix import cheb_matrix

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
# Configuration
# -----------------------------------------------------------------------------
N = 10  # Number of intervals (smaller for clearer visualization)
N_FINE = 500  # Fine grid for plotting

# Book color scheme
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'
PURPLE = '#8E44AD'
ORANGE = '#E67E22'

# Output paths
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch07' / 'python'


def chebyshev_cardinal(x_eval, x_nodes, j):
    """
    Evaluate the j-th Chebyshev cardinal function at points x_eval.

    The cardinal function L_j(x) is the unique polynomial of degree N
    that satisfies L_j(x_k) = delta_{jk} (Kronecker delta).

    Parameters
    ----------
    x_eval : ndarray
        Points at which to evaluate the cardinal function
    x_nodes : ndarray
        Chebyshev grid points
    j : int
        Index of the cardinal function (0, 1, ..., N)

    Returns
    -------
    L_j : ndarray
        Values of L_j at x_eval
    """
    N = len(x_nodes) - 1

    # Use barycentric interpolation formula
    # L_j(x) = (w_j / (x - x_j)) / sum_k (w_k / (x - x_k))
    # where w_k are barycentric weights

    # For Chebyshev points, weights are: w_k = (-1)^k * (1/2 if k=0,N else 1)
    w = np.ones(N + 1)
    w[0] = 0.5
    w[N] = 0.5
    w = w * ((-1) ** np.arange(N + 1))

    L_j = np.zeros_like(x_eval)

    for i, x in enumerate(x_eval):
        # Check if x is at a node
        diffs = x - x_nodes
        if np.any(np.abs(diffs) < 1e-14):
            # x is at a node
            k = np.argmin(np.abs(diffs))
            L_j[i] = 1.0 if k == j else 0.0
        else:
            # Barycentric formula
            terms = w / diffs
            L_j[i] = (w[j] / diffs[j]) / np.sum(terms)

    return L_j


def main():
    """Create Chebyshev cardinal function visualization."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Chebyshev grid points
    D, x_nodes = cheb_matrix(N)
    x_fine = np.linspace(-1, 1, N_FINE)

    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))

    # Panel 1: Several cardinal functions
    ax1 = axes[0]

    # Select a few cardinal functions to display
    indices = [0, 3, 5, N]  # boundary, interior, interior, boundary
    colors = [CORAL, TEAL, PURPLE, ORANGE]
    linestyles = ['-', '--', '-.', ':']

    for idx, color, ls in zip(indices, colors, linestyles):
        L_j = chebyshev_cardinal(x_fine, x_nodes, idx)
        ax1.plot(x_fine, L_j, ls, color=color, linewidth=1.5,
                 label=f'$L_{{{idx}}}(x)$')
        # Mark the cardinal point where L_idx = 1
        ax1.plot(x_nodes[idx], 1.0, 'o', color=color, markersize=7,
                 markeredgecolor='white', markeredgewidth=1)

    # Mark all nodes on x-axis
    ax1.plot(x_nodes, np.zeros(N + 1), 'o', color=NAVY, markersize=5,
             markerfacecolor=NAVY, markeredgecolor='white', markeredgewidth=0.5)

    # Reference lines
    ax1.axhline(0, color='#888888', linewidth=0.5)
    ax1.axhline(1, color='#888888', linewidth=0.5, linestyle='--', alpha=0.5)

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$L_j(x)$')
    ax1.set_title(f'Chebyshev Cardinal Functions ($N = {N}$)', fontsize=11)
    ax1.legend(loc='upper right', fontsize=9, ncol=2)
    ax1.set_xlim(-1.05, 1.05)
    ax1.set_ylim(-0.3, 1.2)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')

    # Add annotation
    ax1.text(0.5, 0.85, r'$L_j(x_k) = \delta_{jk}$',
             transform=ax1.transAxes, fontsize=10, color='#666666')

    # Panel 2: Cardinal function with derivative slopes at nodes
    ax2 = axes[1]

    # Plot cardinal function L_5 and its derivatives
    j = 5
    L_j = chebyshev_cardinal(x_fine, x_nodes, j)
    ax2.plot(x_fine, L_j, '-', color=TEAL, linewidth=2, label=f'$L_{{{j}}}(x)$')

    # Mark the cardinal point
    ax2.plot(x_nodes[j], 1.0, 'o', color=TEAL, markersize=10,
             markeredgecolor='white', markeredgewidth=2)

    # Draw tangent lines at grid points (slopes from D matrix)
    for k in range(0, N + 1, 2):  # Every other point for clarity
        slope = D[k, j]  # This is L_j'(x_k)
        x_k = x_nodes[k]
        y_k = chebyshev_cardinal(np.array([x_k]), x_nodes, j)[0]

        # Draw short tangent line
        dx = 0.15
        x_tan = np.array([x_k - dx, x_k + dx])
        y_tan = y_k + slope * (x_tan - x_k)

        ax2.plot(x_tan, y_tan, '-', color=CORAL, linewidth=1.5, alpha=0.7)
        ax2.plot(x_k, y_k, 's', color=NAVY, markersize=4)

    # Reference line
    ax2.axhline(0, color='#888888', linewidth=0.5)

    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r'$L_j(x)$')
    ax2.set_title(f'Cardinal Function $L_{{{j}}}$ with Derivative Slopes', fontsize=11)
    ax2.set_xlim(-1.05, 1.05)
    ax2.set_ylim(-0.5, 1.2)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')

    # Add annotation
    ax2.text(0.02, 0.15, f'Slopes = column {j} of $D_N$:\n' +
             r'$D_{kj} = L_j^\prime(x_k)$',
             transform=ax2.transAxes, fontsize=9, color='#666666',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

    plt.tight_layout()

    # Save figure
    output_file = OUTPUT_DIR / 'cheb_cardinal.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig)

    # Verify cardinal property
    print(f"\nVerification of cardinal property (N = {N}):")
    print("-" * 50)
    print(f"{'j':>4} {'L_j(x_j)':>12} {'max|L_j(x_k)|, k!=j':>20}")
    print("-" * 50)
    for j in range(N + 1):
        L_j_at_nodes = chebyshev_cardinal(x_nodes, x_nodes, j)
        val_at_own = L_j_at_nodes[j]
        val_at_others = np.delete(L_j_at_nodes, j)
        max_at_others = np.max(np.abs(val_at_others))
        print(f"{j:>4d} {val_at_own:>12.6f} {max_at_others:>20.2e}")
    print("-" * 50)


if __name__ == '__main__':
    main()
