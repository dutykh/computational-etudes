#!/usr/bin/env python3
"""
bc_inhom_lifting.py
Chapter 13: Boundary Conditions in Spectral Methods

Inhomogeneous Poisson equation solved via the lifting technique.

    u_xx = exp(4x),  x in (-1, 1),  u(-1) = 0,  u(1) = 1

The lifting function w(x) = (x + 1)/2 satisfies the boundary conditions:
    w(-1) = 0,  w(1) = 1

Since w'' = 0, the homogeneous variable v = u - w satisfies:
    v_xx = exp(4x),  v(-1) = 0,  v(1) = 0

After solving for v with homogeneous Dirichlet BCs, we reconstruct u = v + w.

Generates Figure 13.1: Inhomogeneous Poisson with lifting.

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
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch13' / 'python'


# -----------------------------------------------------------------------------
# Problem setup
# -----------------------------------------------------------------------------
def lifting_function(x):
    """Lifting function w(x) = (x + 1)/2 satisfying w(-1)=0, w(1)=1."""
    return (x + 1.0) / 2.0


def rhs_function(x):
    """Right-hand side: f(x) = exp(4x)."""
    return np.exp(4.0 * x)


def solve_poisson_lifting(N):
    """
    Solve u_xx = exp(4x) with u(-1)=0, u(1)=1 via lifting.

    Steps:
        1. Lifting: w(x) = (x+1)/2, so v = u - w satisfies v(+-1) = 0
        2. Since w'' = 0: v_xx = exp(4x)
        3. Solve for v using homogeneous Dirichlet BCs (matrix stripping)
        4. Reconstruct u = v + w

    Parameters
    ----------
    N : int
        Number of Chebyshev intervals (N+1 grid points)

    Returns
    -------
    x : ndarray
        Chebyshev grid points
    u : ndarray
        Numerical solution at grid points
    """
    D2, D, x = cheb_second_derivative_matrix(N)

    # Interior system for homogeneous Dirichlet BCs
    D2_int = D2[1:N, 1:N]
    f_int = rhs_function(x[1:N])

    # Solve v_xx = f(x) with v(+-1) = 0
    v_int = np.linalg.solve(D2_int, f_int)

    # Assemble full v with BCs
    v = np.zeros(N + 1)
    v[1:N] = v_int

    # Reconstruct u = v + w
    u = v + lifting_function(x)

    return x, u


def compute_reference_solution(N_ref=100):
    """Compute high-resolution reference solution."""
    x_ref, u_ref = solve_poisson_lifting(N_ref)
    return x_ref, u_ref


def main():
    """Create inhomogeneous lifting BVP figure."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Compute reference solution on high-resolution grid
    N_ref = 100
    x_ref, u_ref = compute_reference_solution(N_ref)

    # Interpolate reference onto fine plotting grid
    from numpy.polynomial.chebyshev import chebfit, chebval
    # Use barycentric interpolation via reference grid
    x_fine = np.linspace(-1, 1, 500)
    # Interpolate reference solution onto fine grid using polynomial fit
    coeffs = np.polynomial.chebyshev.chebfit(x_ref, u_ref, N_ref)
    u_fine = np.polynomial.chebyshev.chebval(x_fine, coeffs)

    # Create figure with two panels
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # ---- Panel 1: Solution for N=16 ----
    N = 16
    x, u_num = solve_poisson_lifting(N)

    # Interpolate reference to comparison grid
    u_exact_grid = np.polynomial.chebyshev.chebval(x, coeffs)
    error = np.max(np.abs(u_num - u_exact_grid))

    ax1 = axes[0]
    ax1.plot(x_fine, u_fine, '-', color=NAVY, linewidth=1.5, label='Reference')
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

    # ---- Panel 2: Convergence study ----
    ax2 = axes[1]
    N_values = list(range(4, 42, 2))
    errors = []

    for N_val in N_values:
        x_n, u_n = solve_poisson_lifting(N_val)
        # Compare against reference on the coarser grid
        u_ref_n = np.polynomial.chebyshev.chebval(x_n, coeffs)
        err = max(np.max(np.abs(u_n - u_ref_n)), 1e-16)
        errors.append(err)

    ax2.semilogy(N_values, errors, 'o-', color=CORAL, linewidth=1.5, markersize=5,
                 markeredgecolor='white', markeredgewidth=0.5)

    ax2.set_xlabel(r'$N$')
    ax2.set_ylabel('Max Error')
    ax2.set_title('Convergence', fontsize=11)
    ax2.set_xlim(0, 44)
    ax2.set_ylim(1e-16, 1e1)
    ax2.grid(True, which='both', alpha=0.3)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # Machine precision line
    ax2.axhline(2.2e-16, color='gray', linestyle=':', linewidth=1, alpha=0.7)
    ax2.text(42, 4e-16, 'Machine precision', fontsize=9, color='gray',
             ha='right', va='bottom')

    # Main title
    fig.suptitle(r'Inhomogeneous Poisson: $u_{xx} = e^{4x}$, $u(-1)=0$, $u(1)=1$ (Lifting)',
                 fontsize=13, y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    # Save figure
    output_file = OUTPUT_DIR / 'inhom_lifting.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig)

    # Print convergence table
    print("\nConvergence Table: Inhomogeneous Poisson with Lifting")
    print("-" * 40)
    print(f"{'N':>6} {'Max Error':>14}")
    print("-" * 40)
    for N_val, err in zip(N_values, errors):
        print(f"{N_val:>6d} {err:>14.2e}")
    print("-" * 40)


if __name__ == '__main__':
    main()
