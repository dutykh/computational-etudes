#!/usr/bin/env python3
"""
chebfft_accuracy.py
Chapter 10: Spectral PDE Solvers with Chebyshev and Fourier Grids

Creates Figure 10.2: Demonstration of spectral accuracy for Chebyshev
differentiation via FFT.

Author: Dr. Denys Dutykh
Date: February 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pathlib import Path

from chebfft import chebfft

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


def f(x):
    """Test function f(x) = exp(x) * cos(4x)."""
    return np.exp(x) * np.cos(4 * x)


def f_prime(x):
    """Exact derivative f'(x) = exp(x) * (cos(4x) - 4*sin(4x))."""
    return np.exp(x) * (np.cos(4 * x) - 4 * np.sin(4 * x))


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))

    # Fine grid for plotting
    x_fine = np.linspace(-1, 1, 500)

    # =========================================================================
    # Panel 1: Function with Chebyshev nodes
    # =========================================================================
    ax = axes[0]

    N = 16
    x_cheb = np.cos(np.pi * np.arange(N + 1) / N)
    f_cheb = f(x_cheb)

    ax.plot(x_fine, f(x_fine), '-', color=NAVY, linewidth=1.5, label=r'$f(x) = e^x \cos(4x)$')
    ax.plot(x_cheb, f_cheb, 'o', color=CORAL, markersize=6, label=f'Chebyshev points ($N={N}$)')

    ax.set_xlabel(r'$x$', fontsize=12)
    ax.set_ylabel(r'$f(x)$', fontsize=12)
    ax.set_title('Function and Chebyshev nodes', fontsize=12)
    ax.legend(loc='lower left', fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-1.1, 1.1)

    # =========================================================================
    # Panel 2: Derivative comparison for N=16
    # =========================================================================
    ax = axes[1]

    f_prime_exact = f_prime(x_cheb)
    f_prime_fft = chebfft(f_cheb)

    ax.plot(x_fine, f_prime(x_fine), '-', color=NAVY, linewidth=1.5, label='Exact')
    ax.plot(x_cheb, f_prime_fft, 'o', color=CORAL, markersize=6, label='FFT')

    ax.set_xlabel(r'$x$', fontsize=12)
    ax.set_ylabel(r"$f'(x)$", fontsize=12)
    ax.set_title(f"Derivative comparison ($N={N}$)", fontsize=12)
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-1.1, 1.1)

    # =========================================================================
    # Panel 3: Convergence with N
    # =========================================================================
    ax = axes[2]

    N_values = np.arange(4, 32, 2)
    errors = []

    for N in N_values:
        x_cheb = np.cos(np.pi * np.arange(N + 1) / N)
        f_cheb = f(x_cheb)
        f_prime_exact = f_prime(x_cheb)
        f_prime_fft = chebfft(f_cheb)
        err = np.max(np.abs(f_prime_fft - f_prime_exact))
        errors.append(err)

    ax.semilogy(N_values, errors, 'o-', color=NAVY, linewidth=1.5, markersize=5)

    ax.set_xlabel(r'$N$', fontsize=12)
    ax.set_ylabel(r'Maximum error', fontsize=12)
    ax.set_title('Spectral convergence', fontsize=12)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_ylim(1e-16, 1)

    # Add reference line
    ax.axhline(1e-14, color=CORAL, linestyle='--', linewidth=0.8, label='Machine precision')
    ax.legend(loc='upper right', fontsize=9)

    plt.tight_layout()

    # Save
    output_path = OUTPUT_DIR / 'chebfft_accuracy.pdf'
    plt.savefig(output_path, bbox_inches='tight')
    print(f"Figure saved to {output_path}")
    plt.show()


if __name__ == '__main__':
    main()
