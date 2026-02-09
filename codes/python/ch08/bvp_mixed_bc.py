#!/usr/bin/env python3
"""
bvp_mixed_bc.py

Solves a BVP with mixed boundary conditions using Chebyshev spectral collocation:

    u_xx = -cos(πx/2),  x ∈ (-1, 1)
    u(-1) = 0           (Dirichlet: fixed temperature)
    u'(1) = 0           (Neumann: insulated boundary)

Exact solution: u(x) = (4/π²)cos(πx/2) + (2/π)(x + 1)

This demonstrates the "boundary bordering" technique: instead of extracting
the interior system (matrix stripping), we keep the full (N+1)×(N+1) system
and replace boundary rows with the appropriate conditions.

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

# Import Chebyshev matrix function from Chapter 7
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
    """Right-hand side: f(x) = -cos(πx/2)."""
    return -np.cos(np.pi * x / 2)


def exact_solution(x):
    """Exact solution: u(x) = (4/π²)cos(πx/2) + (2/π)(x + 1)."""
    return (4 / np.pi**2) * np.cos(np.pi * x / 2) + (2 / np.pi) * (x + 1)


def solve_mixed_bc(N):
    """
    Solve u_xx = f(x) with u(-1) = 0, u'(1) = 0 using boundary bordering.

    Instead of extracting the interior system (matrix stripping), we keep
    the full (N+1)×(N+1) system and replace boundary rows:
      - Row 0 (x=1, Neumann):  L[0,:] = D[0,:]   →  u'(1) = 0
      - Row N (x=-1, Dirichlet): L[N,:] = e_N     →  u(-1) = 0

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
    D2, D, x = cheb_second_derivative_matrix(N)

    # Build full system: L u = rhs
    L = D2.copy()
    rhs = rhs_function(x)

    # Neumann condition at x[0] = 1: u'(1) = 0
    L[0, :] = D[0, :]
    rhs[0] = 0.0

    # Dirichlet condition at x[N] = -1: u(-1) = 0
    L[N, :] = 0.0
    L[N, N] = 1.0
    rhs[N] = 0.0

    # Solve
    u = np.linalg.solve(L, rhs)

    return x, u


def main():
    """Create mixed BC figure."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Fine grid for exact solution
    x_fine = np.linspace(-1, 1, 500)
    u_exact_fine = exact_solution(x_fine)

    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Solve for N = 16
    N = 16
    x, u_num = solve_mixed_bc(N)
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
    ax1.legend(loc='upper left', fontsize=9)
    ax1.set_xlim(-1.05, 1.05)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
    ax1.xaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')

    # Annotate boundary conditions
    ax1.annotate(r'$u(-1) = 0$', xy=(-1, 0), xytext=(-0.7, 0.2),
                 fontsize=10, color=PURPLE,
                 arrowprops=dict(arrowstyle='->', color=PURPLE, alpha=0.7))
    # Neumann BC: show horizontal tangent line at x=1
    u_right = exact_solution(1.0)
    ax1.plot([0.9, 1.05], [u_right, u_right], '-', color=CORAL, linewidth=2)
    ax1.annotate(r"$u'(1) = 0$", xy=(1.0, u_right), xytext=(0.45, 1.1),
                 fontsize=10, color=CORAL,
                 arrowprops=dict(arrowstyle='->', color=CORAL, alpha=0.7))

    # Panel 2: Convergence study
    ax2 = axes[1]
    N_values = [4, 6, 8, 10, 12, 16, 20, 24, 28, 32, 36, 40, 48]
    errors = []

    for N_val in N_values:
        x_n, u_n = solve_mixed_bc(N_val)
        u_exact_n = exact_solution(x_n)
        err = max(np.max(np.abs(u_n - u_exact_n)), 1e-16)
        errors.append(err)

    ax2.semilogy(N_values, errors, 'o-', color=CORAL, linewidth=1.5, markersize=6,
                 markeredgecolor='white', markeredgewidth=0.5)

    ax2.set_xlabel(r'$N$')
    ax2.set_ylabel('Max Error')
    ax2.set_title('Convergence', fontsize=11)
    ax2.set_xlim(0, 52)
    ax2.set_ylim(1e-16, 1e1)
    ax2.grid(True, which='both', alpha=0.3)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # Add machine precision line
    ax2.axhline(2.2e-16, color='gray', linestyle=':', linewidth=1, alpha=0.7)
    ax2.text(50, 4e-16, 'Machine precision', fontsize=9, color='gray',
             ha='right', va='bottom')

    # Main title
    fig.suptitle(r"Mixed BCs: $u''= -\cos(\pi x/2)$, $u(-1)=0$, $u'(1)=0$",
                 fontsize=13, y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    # Save figure
    output_file = OUTPUT_DIR / 'mixed_bc.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig)

    # Print convergence table
    print("\nConvergence Table: Mixed BC Problem")
    print("-" * 40)
    print(f"{'N':>6} {'Max Error':>14}")
    print("-" * 40)
    for N_val, err in zip(N_values, errors):
        print(f"{N_val:>6d} {err:>14.2e}")
    print("-" * 40)


if __name__ == '__main__':
    main()
