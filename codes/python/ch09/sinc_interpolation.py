#!/usr/bin/env python3
"""
sinc_interpolation.py
Chapter 9: Physical and Fourier Space on Grids

Demonstrates band-limited interpolation using the sinc function.
Shows interpolation of delta, square wave, and hat functions.
This is Etude 2 from the chapter.

Author: Dr. Denys Dutykh
Date: February 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def sinc_kernel(z):
    """
    Compute sinc(z) = sin(pi*z)/(pi*z), with sinc(0) = 1.

    Parameters:
        z : array - Input values

    Returns:
        out : array - sinc(z) values
    """
    out = np.ones_like(z, dtype=float)
    mask = z != 0
    out[mask] = np.sin(np.pi * z[mask]) / (np.pi * z[mask])
    return out


def sinc_interpolate(x_grid, v, x_fine, h=1.0):
    """
    Band-limited interpolation via sinc functions.

    Parameters:
        x_grid : array - Grid points
        v : array - Function values at grid points
        x_fine : array - Fine grid for interpolation
        h : float - Grid spacing

    Returns:
        p : array - Interpolated values at x_fine
    """
    p = np.zeros_like(x_fine)
    for xi, vi in zip(x_grid, v):
        p += vi * sinc_kernel((x_fine - xi) / h)
    return p


def main():
    # Grid parameters
    h = 1.0
    xmax = 10
    x_grid = np.arange(-xmax, xmax + h, h)  # Grid points
    x_fine = np.arange(-xmax - h/20, xmax + h/20, h/10)  # Fine plotting grid

    # Create figure with 3 panels
    fig, axes = plt.subplots(3, 1, figsize=(8, 8))

    # Define three grid functions
    functions = [
        ('Discrete delta function', (x_grid == 0).astype(float)),
        ('Discrete square wave', (np.abs(x_grid) <= 3).astype(float)),
        ('Discrete hat function', np.maximum(0, 1 - np.abs(x_grid) / 3.0))
    ]

    for ax, (title, v) in zip(axes, functions):
        # Plot discrete samples
        ax.plot(x_grid, v, 'ko', markersize=8, label='Grid values')

        # Compute and plot interpolant
        p = sinc_interpolate(x_grid, v, x_fine, h)
        ax.plot(x_fine, p, 'b-', linewidth=1.0, label='Band-limited interpolant')

        ax.set_xlim(-xmax, xmax)
        ax.set_ylim(-0.5, 1.5)
        ax.set_yticks([0, 0.5, 1])
        ax.set_title(title, fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='k', linewidth=0.5)
        ax.legend(loc='upper right', fontsize=9)

    axes[-1].set_xlabel(r'$x$', fontsize=12)

    plt.tight_layout()

    # Save figure
    output_dir = Path(__file__).parent.parent.parent.parent / 'textbook' / 'figures' / 'ch09' / 'python'
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / 'sinc_interpolation.pdf', bbox_inches='tight')
    print(f"Figure saved to {output_dir / 'sinc_interpolation.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
