#!/usr/bin/env python3
"""
cheb_matrix_structure.py

Visualizes the structure of the Chebyshev differentiation matrix: heatmap
showing the matrix entries and row profile showing the decay pattern.

This script generates Figure 6.2 for Chapter 6.

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
from matplotlib.colors import TwoSlopeNorm
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
N = 16  # Number of intervals

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


def redblue_cmap():
    """Create a red-white-blue diverging colormap."""
    from matplotlib.colors import LinearSegmentedColormap
    colors = [(0.0, 0.0, 1.0),    # Blue
              (1.0, 1.0, 1.0),    # White
              (1.0, 0.0, 0.0)]    # Red
    return LinearSegmentedColormap.from_list('redblue', colors)


def main():
    """Create matrix structure visualization figure."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Construct differentiation matrix
    D, x = cheb_matrix(N)

    # Create figure with two panels
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))

    # Panel 1: Heatmap of D matrix
    ax1 = axes[0]
    vmax = np.max(np.abs(D))

    # Use diverging colormap centered at zero
    im1 = ax1.imshow(D, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='equal')
    cbar = plt.colorbar(im1, ax=ax1, shrink=0.8)
    cbar.ax.set_ylabel('Matrix entry value', fontsize=10)

    ax1.set_title(f'Chebyshev Differentiation Matrix $D_N$ ($N = {N}$)', fontsize=11)
    ax1.set_xlabel('Column index $j$')
    ax1.set_ylabel('Row index $i$')

    # Set tick labels
    tick_positions = [0, N//4, N//2, 3*N//4, N]
    ax1.set_xticks(tick_positions)
    ax1.set_yticks(tick_positions)
    ax1.set_xticklabels([str(i) for i in tick_positions])
    ax1.set_yticklabels([str(i) for i in tick_positions])

    # Add annotation about large corner entries
    ax1.annotate(f'$D_{{00}} = {D[0,0]:.1f}$', xy=(0, 0), xytext=(4, 4),
                 fontsize=9, color='white',
                 bbox=dict(boxstyle='round', facecolor=NAVY, alpha=0.8))
    ax1.annotate(f'$D_{{N,N}} = {D[N,N]:.1f}$', xy=(N, N), xytext=(N-5, N-3),
                 fontsize=9, color='white',
                 bbox=dict(boxstyle='round', facecolor=NAVY, alpha=0.8))

    # Panel 2: Row profile (first row and middle row)
    ax2 = axes[1]

    j_indices = np.arange(N + 1)

    # First row (boundary)
    markerline, stemlines, baseline = ax2.stem(j_indices, D[0, :], basefmt=' ',
                                                label=f'Row $i=0$ (boundary)')
    markerline.set_color(CORAL)
    markerline.set_marker('o')
    stemlines.set_color(CORAL)

    # Middle row
    mid_row = N // 2
    markerline2, stemlines2, baseline2 = ax2.stem(j_indices + 0.15, D[mid_row, :],
                                                   basefmt=' ', label=f'Row $i={mid_row}$ (interior)')
    markerline2.set_color(TEAL)
    markerline2.set_marker('s')
    stemlines2.set_color(TEAL)

    ax2.axhline(0, color='#888888', linewidth=0.5)
    ax2.set_xlabel('Column index $j$')
    ax2.set_ylabel('Matrix entry $D_{ij}$')
    ax2.set_title('Row Profiles of $D_N$', fontsize=11)
    ax2.legend(loc='upper right', fontsize=9)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # Add grid
    ax2.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
    ax2.set_xlim(-0.5, N + 0.5)

    plt.tight_layout()

    # Save figure
    output_file = OUTPUT_DIR / 'cheb_matrix_structure.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig)

    # Print matrix properties
    print(f"\nChebyshev Matrix Properties (N = {N}):")
    print("-" * 50)
    print(f"  Matrix size: ({N+1}, {N+1})")
    print(f"  Corner entry D[0,0]: {D[0,0]:.6f}")
    print(f"  Corner entry D[N,N]: {D[N,N]:.6f}")
    print(f"  Theoretical D[0,0] = (2N^2+1)/6: {(2*N**2+1)/6:.6f}")
    print(f"  Max entry: {np.max(D):.6f}")
    print(f"  Min entry: {np.min(D):.6f}")
    print(f"  Negative sum trick error: {np.max(np.abs(D @ np.ones(N+1))):.2e}")
    print("-" * 50)


if __name__ == '__main__':
    main()
