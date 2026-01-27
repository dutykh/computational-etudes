#!/usr/bin/env python3
"""
cheb_diff_demo.py

Demonstrates Chebyshev spectral differentiation using the "Witch of Agnesi"
test function u(x) = 1/(1 + 4x^2). Shows the function, its derivative,
and the spectral accuracy for different values of N.

This script generates Figure 6.4 for Chapter 6.

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
N_FINE = 500  # Fine grid for plotting exact solutions

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


# -----------------------------------------------------------------------------
# Test function: Witch of Agnesi
# -----------------------------------------------------------------------------
def u_func(x):
    """u(x) = 1 / (1 + 4x^2) - the Witch of Agnesi."""
    return 1.0 / (1.0 + 4.0 * x**2)


def u_prime_exact(x):
    """u'(x) = -8x / (1 + 4x^2)^2."""
    return -8.0 * x / (1.0 + 4.0 * x**2)**2


def main():
    """Create Chebyshev differentiation demonstration figure."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Fine grid for exact solution
    x_fine = np.linspace(-1, 1, N_FINE)

    # Create 2x2 figure
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    # Two values of N to demonstrate
    N_values = [10, 20]

    for col, N in enumerate(N_values):
        # Construct differentiation matrix
        D, x = cheb_matrix(N)

        # Evaluate function at Chebyshev points
        v = u_func(x)

        # Compute numerical derivative
        w = D @ v

        # Exact derivative at Chebyshev points
        w_exact = u_prime_exact(x)

        # Error
        error = w - w_exact
        max_error = np.max(np.abs(error))

        # Top row: Function and numerical points
        ax_top = axes[0, col]
        ax_top.plot(x_fine, u_func(x_fine), '-', color=NAVY, linewidth=1.5,
                    label='$u(x) = 1/(1+4x^2)$')
        ax_top.plot(x, v, 'o', color=CORAL, markersize=6,
                    markeredgecolor='white', markeredgewidth=0.5,
                    label=f'Chebyshev points')

        ax_top.set_xlabel(r'$x$')
        ax_top.set_ylabel(r'$u(x)$')
        ax_top.set_title(f'Function with $N = {N}$ points', fontsize=11)
        ax_top.legend(loc='upper right', fontsize=9)
        ax_top.set_xlim(-1.05, 1.05)
        ax_top.spines['top'].set_visible(False)
        ax_top.spines['right'].set_visible(False)
        ax_top.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
        ax_top.xaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')

        # Bottom row: Error in derivative
        ax_bot = axes[1, col]

        # Plot exact derivative
        ax_bot.plot(x_fine, u_prime_exact(x_fine), '-', color=NAVY, linewidth=1.5,
                    label="Exact $u'(x)$")
        ax_bot.plot(x, w, 'o', color=TEAL, markersize=6,
                    markeredgecolor='white', markeredgewidth=0.5,
                    label='Spectral')

        ax_bot.set_xlabel(r'$x$')
        ax_bot.set_ylabel(r"$u'(x)$")
        ax_bot.set_title(f'Derivative with $N = {N}$ (max error: {max_error:.2e})',
                         fontsize=11)
        ax_bot.legend(loc='upper right', fontsize=9)
        ax_bot.set_xlim(-1.05, 1.05)
        ax_bot.spines['top'].set_visible(False)
        ax_bot.spines['right'].set_visible(False)
        ax_bot.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
        ax_bot.xaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')

    # Main title
    fig.suptitle(r'Chebyshev Differentiation of $u(x) = 1/(1+4x^2)$ (Witch of Agnesi)',
                 fontsize=13, y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # Save figure
    output_file = OUTPUT_DIR / 'cheb_diff_demo.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig)

    # Print convergence table
    print(f"\nChebyshev differentiation accuracy for u(x) = 1/(1+4x^2):")
    print("-" * 50)
    print(f"{'N':>6} {'Max Error':>14} {'Convergence':>14}")
    print("-" * 50)

    prev_error = None
    for N in [4, 8, 10, 12, 16, 20, 24, 32, 40, 48, 64]:
        D, x = cheb_matrix(N)
        v = u_func(x)
        w = D @ v
        w_exact = u_prime_exact(x)
        error = np.max(np.abs(w - w_exact))

        if prev_error is not None and error > 0:
            rate = np.log10(prev_error / error) / np.log10(2)
            print(f"{N:>6d} {error:>14.6e} {rate:>14.2f}")
        else:
            print(f"{N:>6d} {error:>14.6e} {'---':>14}")

        prev_error = error

    print("-" * 50)
    print("Note: The Witch of Agnesi is analytic in [-1,1] with poles at x = ±i/2,")
    print("so exponential convergence is expected until machine precision is reached.")


if __name__ == '__main__':
    main()
