#!/usr/bin/env python3
"""
cheb_grid_comparison.py

Visualizes the comparison between equispaced and Chebyshev-Gauss-Lobatto grids.
Shows why Chebyshev points cluster near the boundaries and how this relates
to the projection from a circle.

This script generates Figure 6.1 for Chapter 6.

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
# Configuration
# -----------------------------------------------------------------------------
N = 16  # Number of grid points

# Book color scheme
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'
PURPLE = '#8E44AD'
ORANGE = '#E67E22'

# Output paths
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch06' / 'python'


def main():
    """Create grid comparison figure."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Grid points
    x_equi = np.linspace(-1, 1, N + 1)  # Equispaced
    x_cheb = np.cos(np.pi * np.arange(N + 1) / N)  # Chebyshev-GL

    # Create figure with three panels
    fig = plt.figure(figsize=(10, 7))

    # Use gridspec for flexible layout
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.5], width_ratios=[1, 1],
                          hspace=0.6, wspace=0.3)

    # Panel 1: Equispaced grid (top-left)
    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(x_equi, np.zeros_like(x_equi), 'o', color=CORAL, markersize=8,
             markerfacecolor=CORAL, markeredgecolor='white', markeredgewidth=1)
    ax1.axhline(0, color='#888888', linewidth=0.5, zorder=0)
    ax1.set_xlim(-1.1, 1.1)
    ax1.set_ylim(-0.2, 0.2)
    ax1.set_xlabel(r'$x$')
    ax1.set_title(f'Equispaced Grid ($N = {N}$ intervals)', fontsize=11)
    ax1.set_yticks([])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)

    # Add tick labels for key points
    ax1.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax1.set_xticklabels([r'$-1$', r'$-0.5$', r'$0$', r'$0.5$', r'$1$'])

    # Panel 2: Chebyshev grid (middle)
    ax2 = fig.add_subplot(gs[1, :])
    ax2.plot(x_cheb, np.zeros_like(x_cheb), 'o', color=TEAL, markersize=8,
             markerfacecolor=TEAL, markeredgecolor='white', markeredgewidth=1)
    ax2.axhline(0, color='#888888', linewidth=0.5, zorder=0)
    ax2.set_xlim(-1.1, 1.1)
    ax2.set_ylim(-0.2, 0.2)
    ax2.set_xlabel(r'$x$')
    ax2.set_title(f'Chebyshev-Gauss-Lobatto Grid ($N = {N}$ intervals)', fontsize=11)
    ax2.set_yticks([])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax2.set_xticklabels([r'$-1$', r'$-0.5$', r'$0$', r'$0.5$', r'$1$'])

    # Panel 3: Circle projection (bottom-left)
    ax3 = fig.add_subplot(gs[2, 0])

    # Draw circle
    theta = np.linspace(0, np.pi, 200)
    ax3.plot(np.cos(theta), np.sin(theta), '-', color=SKY, linewidth=2, alpha=0.8)

    # Draw projection lines and points for Chebyshev
    theta_cheb = np.pi * np.arange(N + 1) / N
    for i, t in enumerate(theta_cheb):
        x = np.cos(t)
        y = np.sin(t)
        # Point on circle
        ax3.plot(x, y, 'o', color=TEAL, markersize=6)
        # Projection line
        ax3.plot([x, x], [y, 0], '--', color='#888888', linewidth=0.5, alpha=0.6)
        # Point on x-axis
        ax3.plot(x, 0, 's', color=NAVY, markersize=5)

    ax3.axhline(0, color='#444444', linewidth=0.8)
    ax3.axvline(0, color='#444444', linewidth=0.5, alpha=0.5)
    ax3.set_xlim(-1.3, 1.3)
    ax3.set_ylim(-0.2, 1.3)
    ax3.set_aspect('equal')
    ax3.set_xlabel(r'$x = \cos\theta$')
    ax3.set_ylabel(r'$\sin\theta$')
    ax3.set_title('Circle Projection: ' + r'$x_j = \cos(j\pi/N)$', fontsize=11)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    # Panel 4: Boundary spacing zoom (bottom-right)
    ax4 = fig.add_subplot(gs[2, 1])

    # Compute spacing for both grids
    spacing_equi = np.diff(np.sort(x_equi))
    spacing_cheb = np.diff(np.sort(x_cheb))

    # Bar plot of spacings near boundary
    n_show = 8  # Show first 8 spacings
    x_pos = np.arange(n_show)
    width = 0.35

    ax4.bar(x_pos - width/2, spacing_equi[:n_show], width, color=CORAL,
            alpha=0.8, label='Equispaced')
    ax4.bar(x_pos + width/2, spacing_cheb[:n_show], width, color=TEAL,
            alpha=0.8, label='Chebyshev')

    ax4.set_xlabel('Interval index (from left boundary)')
    ax4.set_ylabel(r'Spacing $\Delta x$')
    ax4.set_title('Grid Spacing Near Boundary', fontsize=11)
    ax4.legend(loc='upper right', fontsize=9)
    ax4.set_xticks(x_pos)
    ax4.set_xticklabels([str(i+1) for i in x_pos])
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)

    # Add annotation about O(N^-2) clustering
    cheb_boundary_spacing = spacing_cheb[0]
    equi_spacing = spacing_equi[0]
    ax4.annotate(f'Chebyshev: ' + r'$\Delta x_1 \approx$' + f' {cheb_boundary_spacing:.4f}\n' +
                 r'$(\sim \pi^2/2N^2)$',
                 xy=(0.175, cheb_boundary_spacing + 0.005), xytext=(1.5, 0.14),
                 fontsize=9, color=TEAL, ha='center',
                 arrowprops=dict(arrowstyle='->', color=TEAL, alpha=0.6))

    # Main title
    fig.suptitle('Equispaced vs. Chebyshev-Gauss-Lobatto Grids', fontsize=13, y=0.98)

    # Save figure
    output_file = OUTPUT_DIR / 'grid_comparison.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig)

    # Print spacing information
    print(f"\nGrid spacing analysis (N = {N}):")
    print("-" * 50)
    print(f"{'Spacing':<12} {'Equispaced':>12} {'Chebyshev':>12} {'Ratio':>12}")
    print("-" * 50)
    print(f"{'Minimum':<12} {np.min(spacing_equi):>12.6f} {np.min(spacing_cheb):>12.6f} "
          f"{np.min(spacing_equi)/np.min(spacing_cheb):>12.2f}")
    print(f"{'Maximum':<12} {np.max(spacing_equi):>12.6f} {np.max(spacing_cheb):>12.6f} "
          f"{np.max(spacing_equi)/np.max(spacing_cheb):>12.2f}")
    print(f"{'At boundary':<12} {spacing_equi[0]:>12.6f} {spacing_cheb[0]:>12.6f} "
          f"{spacing_equi[0]/spacing_cheb[0]:>12.2f}")
    print("-" * 50)
    print(f"Theoretical boundary spacing: pi^2/(2*N^2) = {np.pi**2/(2*N**2):.6f}")


if __name__ == '__main__':
    main()
