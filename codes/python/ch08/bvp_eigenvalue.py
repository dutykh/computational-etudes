#!/usr/bin/env python3
"""
bvp_eigenvalue.py

Solves the eigenvalue problem:

    u_xx = λ u,  x ∈ (-1, 1),  u(±1) = 0

Exact eigenvalues: λ_n = -(nπ/2)² for n = 1, 2, 3, ...
Exact eigenfunctions: u_n(x) = sin(nπ(x+1)/2)

This demonstrates the resolution limits of spectral methods - higher modes
require more points per wavelength (ppw) for accuracy.

This script generates Figure 7.6 for Chapter 7.

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

# Import Chebyshev matrix function from Chapter 6
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / 'ch07'))
from cheb_matrix import cheb_matrix, cheb_second_derivative_matrix

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
# Book color scheme
# -----------------------------------------------------------------------------
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'
PURPLE = '#8E44AD'
ORANGE = '#E67E22'

# Output paths
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch08' / 'python'


def exact_eigenvalue(n):
    """Exact eigenvalue λ_n = -(nπ/2)²."""
    return -(n * np.pi / 2)**2


def exact_eigenfunction(x, n):
    """Exact eigenfunction u_n(x) = sin(nπ(x+1)/2)."""
    return np.sin(n * np.pi * (x + 1) / 2)


def solve_eigenvalue_problem(N):
    """
    Solve the eigenvalue problem u_xx = λu with u(±1) = 0.

    Parameters
    ----------
    N : int
        Number of intervals

    Returns
    -------
    eigenvalues : ndarray
        Computed eigenvalues (sorted from largest to smallest, i.e., least negative)
    eigenvectors : ndarray
        Corresponding eigenvectors (columns)
    x : ndarray
        Grid points
    """
    D2, D, x = cheb_second_derivative_matrix(N)

    # Extract interior system
    D2_int = D2[1:N, 1:N]

    # Solve eigenvalue problem
    eigenvalues, eigenvectors = np.linalg.eig(D2_int)

    # Sort by eigenvalue (largest = least negative first)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx].real
    eigenvectors = eigenvectors[:, idx].real

    return eigenvalues, eigenvectors, x


def main():
    """Create eigenvalue problem figure."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    N = 36  # Grid size

    # Solve eigenvalue problem
    eigenvalues, eigenvectors, x = solve_eigenvalue_problem(N)

    # Create figure with 6 panels (2x3)
    fig = plt.figure(figsize=(12, 7))

    # Modes to display (1-indexed)
    modes = [5, 10, 15, 20, 25, 30]
    colors = [CORAL, TEAL, PURPLE, ORANGE, SKY, NAVY]

    # Fine grid for exact solutions
    x_fine = np.linspace(-1, 1, 500)

    for idx, (mode, color) in enumerate(zip(modes, colors)):
        ax = fig.add_subplot(2, 3, idx + 1)

        # Numerical eigenvalue and eigenfunction
        lam_num = eigenvalues[mode - 1]
        u_num = eigenvectors[:, mode - 1]

        # Exact values
        lam_exact = exact_eigenvalue(mode)
        u_exact_fine = exact_eigenfunction(x_fine, mode)
        u_exact_grid = exact_eigenfunction(x[1:N], mode)

        # Normalize numerical eigenfunction to match sign of exact
        if np.sign(u_num[0]) != np.sign(u_exact_grid[0]):
            u_num = -u_num

        # Normalize amplitude
        u_num = u_num / np.max(np.abs(u_num))
        u_exact_fine = u_exact_fine / np.max(np.abs(u_exact_fine))
        u_exact_grid = u_exact_grid / np.max(np.abs(u_exact_grid))

        # Eigenvalue error
        lam_error = np.abs(lam_num - lam_exact) / np.abs(lam_exact)

        # Points per wavelength
        wavelength = 4.0 / mode  # wavelength = 2L/n = 4/n for L=2
        avg_spacing = 2.0 / N
        ppw = wavelength / avg_spacing

        # Plot
        ax.plot(x_fine, u_exact_fine, '-', color=NAVY, linewidth=1.5,
                alpha=0.7, label='Exact')
        ax.plot(x[1:N], u_num, 'o', color=color, markersize=4,
                markeredgecolor='white', markeredgewidth=0.3, label='Spectral')

        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$u(x)$')
        ax.set_title(f'Mode {mode}: ppw = {ppw:.1f}, ' +
                     r'$\lambda$ error = ' + f'{lam_error:.1e}', fontsize=10)
        ax.set_xlim(-1.05, 1.05)
        ax.set_ylim(-1.3, 1.3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        if idx == 0:
            ax.legend(loc='upper right', fontsize=8)

        # Add reference line
        ax.axhline(0, color='#888888', linewidth=0.5)

    # Main title
    fig.suptitle(f'Eigenvalue Problem: ' + r'$u_{xx} = \lambda u$, $u(\pm 1) = 0$' +
                 f' (N = {N})', fontsize=13, y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Save figure
    output_file = OUTPUT_DIR / 'eigenvalue_problem.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig)

    # Print eigenvalue table
    print(f"\nEigenvalue accuracy for N = {N}:")
    print("-" * 70)
    print(f"{'Mode':>6} {'λ (exact)':>14} {'λ (numerical)':>14} "
          f"{'Rel. Error':>12} {'ppw':>8}")
    print("-" * 70)

    for mode in range(1, min(N, 36)):
        lam_exact = exact_eigenvalue(mode)
        lam_num = eigenvalues[mode - 1]
        rel_error = np.abs(lam_num - lam_exact) / np.abs(lam_exact)
        ppw = (4.0 / mode) / (2.0 / N)

        print(f"{mode:>6d} {lam_exact:>14.4f} {lam_num:>14.4f} "
              f"{rel_error:>12.2e} {ppw:>8.2f}")

    print("-" * 70)
    print("\nNote: Accuracy degrades when ppw < π (≈ 3.14).")
    print("The rule of thumb is: need ≥ π points per wavelength for spectral accuracy.")


if __name__ == '__main__':
    main()
