#!/usr/bin/env python3
"""
bvp_variable_coeff.py

Solves an Airy-type equation with variable coefficients:

    u_xx - (1 + x²)u = 1,  x ∈ (-1, 1),  u(±1) = 0

This demonstrates that variable coefficient problems require no additional
complexity with spectral methods - the variable coefficient becomes a
diagonal matrix multiplication.

This script generates Figure 7.4 for Chapter 7.

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


def solve_variable_coeff_bvp(N):
    """
    Solve u_xx - (1 + x²)u = 1 with u(±1) = 0.

    The problem becomes:
        [D² - diag(1 + x²)] u = 1

    Parameters
    ----------
    N : int
        Number of intervals

    Returns
    -------
    x : ndarray
        Chebyshev grid points
    u : ndarray
        Numerical solution
    """
    D2, D, x = cheb_second_derivative_matrix(N)

    # Variable coefficient: a(x) = 1 + x²
    a = 1.0 + x**2

    # Build the operator: L = D² - diag(a)
    L = D2 - np.diag(a)

    # Extract interior system
    L_int = L[1:N, 1:N]
    rhs_int = np.ones(N - 1)  # f(x) = 1

    # Solve
    u_int = np.linalg.solve(L_int, rhs_int)

    # Assemble with boundary conditions
    u = np.zeros(N + 1)
    u[1:N] = u_int

    return x, u


def solve_constant_coeff_bvp(N):
    """
    Solve u_xx - u = 1 with u(±1) = 0 for comparison.

    This has exact solution:
        u(x) = (cosh(1) - cosh(x)) / cosh(1) - 1
             = -1 + (cosh(1) - cosh(x)) / cosh(1)
    """
    D2, D, x = cheb_second_derivative_matrix(N)

    # Constant coefficient: a = 1
    L = D2 - np.eye(N + 1)

    # Extract interior system
    L_int = L[1:N, 1:N]
    rhs_int = np.ones(N - 1)

    # Solve
    u_int = np.linalg.solve(L_int, rhs_int)

    # Assemble with boundary conditions
    u = np.zeros(N + 1)
    u[1:N] = u_int

    return x, u


def exact_constant_coeff(x):
    """
    Exact solution for u_xx - u = 1 with u(±1) = 0.

    General solution: u = C1 cosh(x) + C2 sinh(x) - 1
    BCs: u(1) = C1 cosh(1) + C2 sinh(1) - 1 = 0
         u(-1) = C1 cosh(1) - C2 sinh(1) - 1 = 0
    Adding: 2 C1 cosh(1) - 2 = 0 => C1 = 1/cosh(1)
    Subtracting: 2 C2 sinh(1) = 0 => C2 = 0

    Therefore: u(x) = cosh(x)/cosh(1) - 1
    """
    return np.cosh(x) / np.cosh(1) - 1.0


def main():
    """Create variable coefficient BVP figure."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    N = 32

    # Solve both problems
    x_var, u_var = solve_variable_coeff_bvp(N)
    x_const, u_const = solve_constant_coeff_bvp(N)
    u_const_exact = exact_constant_coeff(x_const)

    # Fine grid for constant coefficient exact solution
    x_fine = np.linspace(-1, 1, 500)
    u_const_exact_fine = exact_constant_coeff(x_fine)

    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Panel 1: Compare solutions
    ax1 = axes[0]
    ax1.plot(x_var, u_var, 'o-', color=TEAL, linewidth=1.5, markersize=4,
             markeredgecolor='white', markeredgewidth=0.3,
             label=r'Variable: $u_{xx} - (1+x^2)u = 1$')
    ax1.plot(x_const, u_const, 's-', color=CORAL, linewidth=1.5, markersize=4,
             markeredgecolor='white', markeredgewidth=0.3,
             label=r'Constant: $u_{xx} - u = 1$')

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$u(x)$')
    ax1.set_title(f'Comparison of Solutions ($N = {N}$)', fontsize=11)
    ax1.legend(loc='lower center', fontsize=9)
    ax1.set_xlim(-1.05, 1.05)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
    ax1.xaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')

    # Annotate the effect of variable coefficient
    ax1.annotate('Variable coefficient\nreduces amplitude',
                 xy=(0, u_var[N//2]), xytext=(0.4, -0.4),
                 fontsize=9, color='#666666',
                 arrowprops=dict(arrowstyle='->', color='#666666', alpha=0.6))

    # Panel 2: Verify constant coefficient case
    ax2 = axes[1]
    ax2.plot(x_fine, u_const_exact_fine, '-', color=NAVY, linewidth=1.5,
             label='Exact')
    ax2.plot(x_const, u_const, 'o', color=CORAL, markersize=6,
             markeredgecolor='white', markeredgewidth=0.5,
             label='Spectral')

    error = np.max(np.abs(u_const - u_const_exact))

    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r'$u(x)$')
    ax2.set_title(f'Constant Coeff. Verification (error: {error:.2e})', fontsize=11)
    ax2.legend(loc='lower center', fontsize=9)
    ax2.set_xlim(-1.05, 1.05)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
    ax2.xaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')

    # Main title
    fig.suptitle(r'Variable Coefficient BVP: $u_{xx} - (1+x^2)u = 1$, $u(\pm 1) = 0$',
                 fontsize=13, y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    # Save figure
    output_file = OUTPUT_DIR / 'variable_coeff.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig)

    # Print solution properties
    print("\nSolution comparison:")
    print("-" * 50)
    print(f"Variable coefficient u(0) = {u_var[N//2]:.6f}")
    print(f"Constant coefficient u(0) = {u_const[N//2]:.6f}")
    print(f"Constant coeff exact u(0) = {exact_constant_coeff(0):.6f}")
    print(f"Numerical error (const coeff): {error:.2e}")
    print("-" * 50)


if __name__ == '__main__':
    main()
