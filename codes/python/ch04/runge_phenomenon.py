#!/usr/bin/env python3
"""
runge_phenomenon.py

Demonstrates the Runge phenomenon: the failure of polynomial interpolation
on equispaced grids for certain smooth functions. Named after Carl Runge
who discovered this counterintuitive behavior in 1901.

The Runge function:
    f(x) = 1 / (1 + 25x²)

is smooth and infinitely differentiable, yet polynomial interpolation on
equispaced nodes diverges as the degree increases. The interpolant oscillates
wildly near the boundaries x = ±1.

This script shows interpolation for N = 6, 10, 14 (equispaced nodes)
demonstrating the progressive deterioration at the boundaries.

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
from numpy.polynomial.polynomial import Polynomial

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
OUTPUT_FILE = OUTPUT_DIR / 'runge_phenomenon.pdf'


# -----------------------------------------------------------------------------
# The Runge function
# -----------------------------------------------------------------------------
def runge(x):
    """
    The Runge function: f(x) = 1 / (1 + 25x²)

    This function has poles in the complex plane at z = ±i/5 = ±0.2i,
    which lie very close to the real axis. These nearby singularities
    are responsible for the divergence of equispaced interpolation.
    """
    return 1.0 / (1.0 + 25.0 * x**2)


# -----------------------------------------------------------------------------
# Lagrange interpolation
# -----------------------------------------------------------------------------
def lagrange_interpolate(x_nodes, f_nodes, x_eval):
    """
    Evaluate the Lagrange interpolating polynomial at points x_eval.

    Given N+1 nodes (x_0, ..., x_N) and values (f_0, ..., f_N),
    the Lagrange interpolant is:

        p_N(x) = Σ_{k=0}^{N} f_k L_k(x)

    where L_k(x) are the Lagrange basis polynomials:

        L_k(x) = Π_{j≠k} (x - x_j) / (x_k - x_j)

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

    # Create figure
    fig, ax = plt.subplots(figsize=(7, 5))

    # Plot exact function
    ax.plot(x_fine, f_exact, color=NAVY, linewidth=2, label='Runge function')

    # Interpolate for each N
    for N, color in zip(N_values, colors):
        # Equispaced nodes: x_j = -1 + 2j/N, j = 0, 1, ..., N
        x_nodes = np.linspace(-1, 1, N + 1)
        f_nodes = runge(x_nodes)

        # Evaluate interpolant
        p_interp = lagrange_interpolate(x_nodes, f_nodes, x_fine)

        # Plot interpolant
        ax.plot(x_fine, p_interp, color=color, linewidth=1.2,
                linestyle='--', label=f'$p_{{{N}}}(x)$ (equispaced)')

        # Plot nodes
        ax.plot(x_nodes, f_nodes, 'o', color=color, markersize=4)

    # Styling
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$f(x)$')
    ax.set_xlim(-1, 1)
    ax.set_ylim(-0.5, 1.5)
    ax.set_title('Runge Phenomenon: Failure of Equispaced Interpolation', fontsize=11)

    ax.legend(loc='upper right', frameon=True, fancybox=False,
              edgecolor='none', facecolor='white', framealpha=0.9)

    # Add annotation about oscillations
    ax.annotate('Oscillations grow\nnear boundaries',
                xy=(-0.9, -0.2), xytext=(-0.6, -0.4),
                fontsize=9, color='gray',
                arrowprops=dict(arrowstyle='->', color='gray', lw=0.8))

    # Clean styling
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#444444')
    ax.spines['bottom'].set_color('#444444')
    ax.tick_params(colors='#444444')
    ax.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
    ax.xaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
    ax.set_axisbelow(True)

    # Add text box explaining the phenomenon
    textstr = r'$f(x) = \frac{1}{1 + 25x^2}$'
    props = dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='none')
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', bbox=props)

    plt.tight_layout()

    # Save figure
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_FILE, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(OUTPUT_FILE.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {OUTPUT_FILE.resolve()}')

    # Print error statistics
    print("\nMaximum interpolation error at boundaries:")
    print("-" * 50)
    for N in N_values:
        x_nodes = np.linspace(-1, 1, N + 1)
        f_nodes = runge(x_nodes)
        p_interp = lagrange_interpolate(x_nodes, f_nodes, x_fine)
        max_error = np.max(np.abs(runge(x_fine) - p_interp))
        print(f"N = {N:2d}: max |f - p_N| = {max_error:.4f}")
    print("-" * 50)
    print("Note: Error INCREASES with N for equispaced nodes!")

    plt.close(fig)


if __name__ == '__main__':
    main()
