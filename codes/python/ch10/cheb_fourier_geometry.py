#!/usr/bin/env python3
"""
cheb_fourier_geometry.py
Chapter 10: Spectral PDE Solvers with Chebyshev and Fourier Grids

Creates Figure 10.1: The Chebyshev-Fourier geometric connection.
Shows the unit circle projection and the "wrapped cosine" interpretation.

Author: Dr. Denys Dutykh
Date: February 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pathlib import Path

# Publication-quality settings
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman', 'CMU Serif', 'DejaVu Serif']
rcParams['mathtext.fontset'] = 'cm'
rcParams['axes.labelsize'] = 11
rcParams['axes.titlesize'] = 12
rcParams['xtick.labelsize'] = 10
rcParams['ytick.labelsize'] = 10
rcParams['figure.dpi'] = 150
rcParams['savefig.dpi'] = 300
rcParams['axes.linewidth'] = 0.8

# Colors
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'

# Output directory
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch10' / 'python'


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # =========================================================================
    # Left panel: Circle projection to Chebyshev points
    # =========================================================================
    ax = axes[0]

    N = 8
    theta = np.pi * np.arange(N + 1) / N
    x_cheb = np.cos(theta)
    y_cheb = np.sin(theta)

    # Draw unit semicircle
    theta_fine = np.linspace(0, np.pi, 200)
    ax.plot(np.cos(theta_fine), np.sin(theta_fine), 'k-', linewidth=1.2)
    ax.axhline(0, color='k', linewidth=0.8)

    # Draw projection lines and points
    for j in range(N + 1):
        # Point on circle
        ax.plot(x_cheb[j], y_cheb[j], 'o', color=NAVY, markersize=8, zorder=5)
        # Projection line
        ax.plot([x_cheb[j], x_cheb[j]], [0, y_cheb[j]], '--', color=SKY, linewidth=0.8)
        # Point on x-axis
        ax.plot(x_cheb[j], 0, 's', color=CORAL, markersize=7, zorder=5)

    # Annotate some angles
    ax.annotate(r'$\theta_0 = 0$', xy=(1.0, 0), xytext=(1.1, 0.2),
                fontsize=10, color=NAVY)
    ax.annotate(r'$\theta_N = \pi$', xy=(-1.0, 0), xytext=(-1.4, 0.2),
                fontsize=10, color=NAVY)

    # Show angle arc
    theta_arc = np.linspace(0, np.pi/4, 50)
    ax.plot(0.3*np.cos(theta_arc), 0.3*np.sin(theta_arc), '-', color=NAVY, linewidth=1)
    ax.annotate(r'$\theta$', xy=(0.35, 0.15), fontsize=11, color=NAVY)

    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-0.3, 1.3)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$x = \cos\theta$', fontsize=12)
    ax.set_ylabel(r'$y = \sin\theta$', fontsize=12)
    ax.set_title(r'Chebyshev points as circle projection', fontsize=12)

    # Label Chebyshev points
    ax.annotate(r'$x_0 = 1$', xy=(1.0, -0.05), xytext=(1.0, -0.2),
                fontsize=9, ha='center', color=CORAL)
    ax.annotate(r'$x_N = -1$', xy=(-1.0, -0.05), xytext=(-1.0, -0.2),
                fontsize=9, ha='center', color=CORAL)

    ax.grid(True, alpha=0.3)

    # =========================================================================
    # Right panel: Wrapped cosine interpretation
    # =========================================================================
    ax = axes[1]

    # Draw T_n(x) = cos(n*theta) as a wrapped cosine
    n = 5
    theta_fine = np.linspace(0, np.pi, 500)
    x_fine = np.cos(theta_fine)

    # Show theta on bottom x-axis, x on top
    ax2 = ax.twiny()

    # Plot cos(n*theta) vs theta
    y_cosine = np.cos(n * theta_fine)
    ax.plot(theta_fine, y_cosine, '-', color=NAVY, linewidth=1.5,
            label=r'$\cos(n\theta)$')

    # Show corresponding T_n(x)
    ax.set_xlabel(r'$\theta$', fontsize=12)
    ax.set_ylabel(r'$\cos(n\theta) = T_n(x)$', fontsize=12)
    ax.set_title(rf'Chebyshev polynomial $T_{n}(x)$ as wrapped cosine', fontsize=12)

    # Mark theta values
    ax.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
    ax.set_xticklabels([r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$'])

    # Top axis shows x = cos(theta)
    ax2.set_xlim(ax.get_xlim())
    x_ticks_theta = [0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi]
    x_ticks_x = [np.cos(t) for t in x_ticks_theta]
    ax2.set_xticks(x_ticks_theta)
    ax2.set_xticklabels([f'{x:.2f}' for x in x_ticks_x])
    ax2.set_xlabel(r'$x = \cos\theta$', fontsize=12)

    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', linewidth=0.5)
    ax.set_ylim(-1.3, 1.3)

    plt.tight_layout()

    # Save
    output_path = OUTPUT_DIR / 'cheb_fourier_geometry.pdf'
    plt.savefig(output_path, bbox_inches='tight')
    print(f"Figure saved to {output_path}")
    plt.show()


if __name__ == '__main__':
    main()
