#!/usr/bin/env python3
"""
aliasing_demo.py
Chapter 9: Physical and Fourier Space on Grids

Demonstrates aliasing: sin(pi*x) and sin(9*pi*x) coincide on grid h=1/4.
This is Etude 1 from the chapter.

Author: Dr. Denys Dutykh
Date: February 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def main():
    # Grid spacing
    h = 0.25

    # Continuous plotting grid
    x_cont = np.linspace(0, 2, 500)

    # Discrete grid points
    x_grid = np.arange(0, 2 + h/2, h)

    # Two functions with different frequencies
    f_cont = np.sin(np.pi * x_cont)      # Wavenumber k = pi
    g_cont = np.sin(9 * np.pi * x_cont)  # Wavenumber k = 9*pi

    # Samples at grid points
    f_samples = np.sin(np.pi * x_grid)
    g_samples = np.sin(9 * np.pi * x_grid)

    # Verify aliasing
    max_diff = np.max(np.abs(f_samples - g_samples))
    print(f"Maximum difference at grid points: {max_diff:.2e}")
    print(f"Aliasing relation: 9π = π + 2π·4 = π + 2π/h")

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 4))

    # Plot continuous functions
    ax.plot(x_cont, f_cont, 'b-', linewidth=1.5, label=r'$\sin(\pi x)$')
    ax.plot(x_cont, g_cont, 'r--', linewidth=1.5, label=r'$\sin(9\pi x)$')

    # Plot grid samples
    ax.plot(x_grid, f_samples, 'ko', markersize=10, markerfacecolor='white',
            markeredgewidth=2, label='Grid samples', zorder=5)

    # Add vertical lines at grid points
    for xg in x_grid:
        ax.axvline(x=xg, color='gray', linewidth=0.5, alpha=0.3)

    ax.set_xlabel(r'$x$', fontsize=12)
    ax.set_ylabel('Value', fontsize=12)
    ax.set_title(r'Aliasing on grid $h = 1/4$: $\sin(\pi x)$ and $\sin(9\pi x)$ coincide', fontsize=12)
    ax.legend(loc='upper right', fontsize=10)
    ax.set_xlim(-0.05, 2.05)
    ax.set_ylim(-1.3, 1.3)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', linewidth=0.5)

    plt.tight_layout()

    # Save figure
    output_dir = Path(__file__).parent.parent.parent.parent / 'textbook' / 'figures' / 'ch09' / 'python'
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / 'aliasing_demo.pdf', bbox_inches='tight')
    print(f"Figure saved to {output_dir / 'aliasing_demo.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
