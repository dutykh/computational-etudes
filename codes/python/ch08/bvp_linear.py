#!/usr/bin/env python3
"""
bvp_linear.py

Solves the 1D Poisson equation using Chebyshev spectral collocation:

    u_xx = sin(πx) + 2cos(2πx),  x ∈ (-1, 1),  u(±1) = 0

The exact solution is determined by integrating twice and applying
boundary conditions.

This script generates Figure 7.3 for Chapter 7.

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


# -----------------------------------------------------------------------------
# Problem setup
# -----------------------------------------------------------------------------
def rhs_function(x):
    """Right-hand side: f(x) = sin(πx) + 2cos(2πx)."""
    return np.sin(np.pi * x) + 2.0 * np.cos(2.0 * np.pi * x)


def exact_solution(x):
    """
    Exact solution of u_xx = sin(πx) + 2cos(2πx) with u(±1) = 0.

    Integrating sin(πx) twice: -sin(πx)/π²
    Integrating 2cos(2πx) twice: -cos(2πx)/(2π²)

    General solution: u(x) = -sin(πx)/π² - cos(2πx)/(2π²) + Ax + B

    Apply BCs:
    u(1) = 0:  -sin(π)/π² - cos(2π)/(2π²) + A + B = 0
              0 - 1/(2π²) + A + B = 0
    u(-1) = 0: -sin(-π)/π² - cos(-2π)/(2π²) - A + B = 0
              0 - 1/(2π²) - A + B = 0

    From these: 2B = 1/π², so B = 1/(2π²), and A = 0

    Therefore: u(x) = -sin(πx)/π² - cos(2πx)/(2π²) + 1/(2π²)
              u(x) = -sin(πx)/π² + (1 - cos(2πx))/(2π²)
    """
    return -np.sin(np.pi * x) / np.pi**2 + (1 - np.cos(2 * np.pi * x)) / (2 * np.pi**2)


def solve_poisson_dirichlet(N):
    """
    Solve u_xx = f(x) with u(±1) = 0 using Chebyshev collocation.

    Parameters
    ----------
    N : int
        Number of intervals (N+1 grid points)

    Returns
    -------
    x : ndarray
        Chebyshev grid points
    u : ndarray
        Numerical solution at grid points
    """
    # Get second derivative matrix and grid
    D2, D, x = cheb_second_derivative_matrix(N)

    # Extract interior points (remove first and last rows/columns)
    # x[0] = 1 (right boundary), x[N] = -1 (left boundary)
    D2_int = D2[1:N, 1:N]  # Interior submatrix
    f_int = rhs_function(x[1:N])  # RHS at interior points

    # Solve the linear system
    u_int = np.linalg.solve(D2_int, f_int)

    # Assemble full solution with boundary conditions
    u = np.zeros(N + 1)
    u[1:N] = u_int
    u[0] = 0.0  # u(1) = 0
    u[N] = 0.0  # u(-1) = 0

    return x, u


def main():
    """Create linear BVP figure."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Fine grid for exact solution
    x_fine = np.linspace(-1, 1, 500)
    u_exact_fine = exact_solution(x_fine)

    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Solve for N = 16
    N = 16
    x, u_num = solve_poisson_dirichlet(N)
    u_exact_grid = exact_solution(x)
    error = np.max(np.abs(u_num - u_exact_grid))

    # Panel 1: Solution
    ax1 = axes[0]
    ax1.plot(x_fine, u_exact_fine, '-', color=NAVY, linewidth=1.5, label='Exact')
    ax1.plot(x, u_num, 'o', color=TEAL, markersize=6,
             markeredgecolor='white', markeredgewidth=0.5, label='Spectral')

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$u(x)$')
    ax1.set_title(f'Solution ($N = {N}$, max error: {error:.2e})', fontsize=11)
    ax1.legend(loc='upper right', fontsize=9)
    ax1.set_xlim(-1.05, 1.05)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
    ax1.xaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')

    # Panel 2: Convergence study
    ax2 = axes[1]
    N_values = [4, 6, 8, 10, 12, 16, 20, 24, 32]
    errors = []

    for N_val in N_values:
        x_n, u_n = solve_poisson_dirichlet(N_val)
        u_exact_n = exact_solution(x_n)
        err = np.max(np.abs(u_n - u_exact_n))
        errors.append(err)

    ax2.semilogy(N_values, errors, 'o-', color=CORAL, linewidth=1.5, markersize=6,
                 markeredgecolor='white', markeredgewidth=0.5)

    ax2.set_xlabel(r'$N$')
    ax2.set_ylabel('Max Error')
    ax2.set_title('Convergence', fontsize=11)
    ax2.set_xlim(0, 35)
    ax2.set_ylim(1e-15, 1)
    ax2.grid(True, which='both', alpha=0.3)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # Add machine epsilon line
    ax2.axhline(1e-14, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax2.text(30, 3e-14, r'$\approx$ machine $\epsilon$', fontsize=9, color='gray')

    # Main title
    fig.suptitle(r'1D Poisson Equation: $u_{xx} = \sin(\pi x) + 2\cos(2\pi x)$, $u(\pm 1) = 0$',
                 fontsize=13, y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    # Save figure
    output_file = OUTPUT_DIR / 'poisson_1d.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig)

    # Print convergence table
    print("\nConvergence Table: 1D Poisson Equation")
    print("-" * 40)
    print(f"{'N':>6} {'Max Error':>14}")
    print("-" * 40)
    for N_val, err in zip(N_values, errors):
        print(f"{N_val:>6d} {err:>14.2e}")
    print("-" * 40)


if __name__ == '__main__':
    main()
