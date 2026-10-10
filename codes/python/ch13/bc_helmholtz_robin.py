#!/usr/bin/env python3
"""
bc_helmholtz_robin.py
Chapter 13: Boundary Conditions in Spectral Methods

Helmholtz equation with Robin boundary conditions (boundary bordering / tau method).

    -u_xx + k^2 u = f(x),  x in (-1, 1)
    u'(-1) - alpha u(-1) = 0   (Robin at left)
    u(1) = 0                    (Dirichlet at right)

Parameters: k = 2, f(x) = exp(-x^2)

The full (N+1) x (N+1) system is built and boundary rows are overwritten:
    Row 0 (x = 1):   identity row for Dirichlet u(1) = 0
    Row N (x = -1):   D[N,:] - alpha * e_N^T for Robin condition

Generates:
    Figure 13.2: Helmholtz solutions for various alpha and convergence
    Figure 13.3: Condition number vs alpha

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
def rhs_function(x):
    """Right-hand side: f(x) = exp(-x^2)."""
    return np.exp(-x**2)


def solve_helmholtz_robin(N, k=2.0, alpha=5.0):
    """
    Solve -u_xx + k^2 u = f(x) with Robin/Dirichlet BCs using boundary bordering.

    Row 0 (x=1, Dirichlet):  L[0,:] = e_0^T          rhs[0] = 0
    Row N (x=-1, Robin):     L[N,:] = D[N,:] - alpha * e_N^T   rhs[N] = 0

    Parameters
    ----------
    N : int
        Number of Chebyshev intervals
    k : float
        Helmholtz parameter
    alpha : float
        Robin coefficient

    Returns
    -------
    x : ndarray
        Chebyshev grid points
    u : ndarray
        Numerical solution
    """
    D2, D, x = cheb_second_derivative_matrix(N)

    # Build operator: L = -D^2 + k^2 I
    I_mat = np.eye(N + 1)
    L = -D2 + k**2 * I_mat
    rhs = rhs_function(x)

    # Row 0 (x[0] = 1): Dirichlet u(1) = 0
    L[0, :] = 0.0
    L[0, 0] = 1.0
    rhs[0] = 0.0

    # Row N (x[N] = -1): Robin u'(-1) - alpha * u(-1) = 0
    L[N, :] = D[N, :] - alpha * I_mat[N, :]
    rhs[N] = 0.0

    # Solve
    u = np.linalg.solve(L, rhs)

    return x, u


def main():
    """Create Helmholtz Robin BC figures."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    k = 2.0

    # =========================================================================
    # Figure 1: Solutions for various alpha + convergence
    # =========================================================================
    fig1, axes1 = plt.subplots(1, 2, figsize=(11, 4.5))

    # ---- Panel 1: Solutions for different alpha ----
    ax1 = axes1[0]
    alpha_values = [0, 1, 5, 20]
    colors_alpha = [NAVY, TEAL, CORAL, PURPLE]

    N_plot = 32
    for alpha, color in zip(alpha_values, colors_alpha):
        x, u = solve_helmholtz_robin(N_plot, k=k, alpha=alpha)
        ax1.plot(x, u, '-o', color=color, linewidth=1.5, markersize=3,
                 markeredgecolor='white', markeredgewidth=0.3,
                 label=rf'$\alpha = {alpha}$')

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$u(x)$')
    ax1.set_title(f'Solutions ($N = {N_plot}$, $k = {k:.0f}$)', fontsize=11)
    ax1.legend(loc='upper left', fontsize=9)
    ax1.set_xlim(-1.05, 1.05)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
    ax1.xaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')

    # ---- Panel 2: Convergence for alpha = 5 ----
    ax2 = axes1[1]
    alpha_conv = 5.0

    # Reference solution with high N
    N_ref = 100
    x_ref, u_ref = solve_helmholtz_robin(N_ref, k=k, alpha=alpha_conv)
    coeffs_ref = np.polynomial.chebyshev.chebfit(x_ref, u_ref, N_ref)

    N_values = list(range(4, 42, 2))
    errors = []

    for N_val in N_values:
        x_n, u_n = solve_helmholtz_robin(N_val, k=k, alpha=alpha_conv)
        u_ref_n = np.polynomial.chebyshev.chebval(x_n, coeffs_ref)
        err = max(np.max(np.abs(u_n - u_ref_n)), 1e-16)
        errors.append(err)

    ax2.semilogy(N_values, errors, 'o-', color=CORAL, linewidth=1.5, markersize=5,
                 markeredgecolor='white', markeredgewidth=0.5)

    ax2.set_xlabel(r'$N$')
    ax2.set_ylabel('Max Error')
    ax2.set_title(rf'Convergence ($\alpha = {alpha_conv:.0f}$)', fontsize=11)
    ax2.set_xlim(0, 44)
    ax2.set_ylim(1e-16, 1e1)
    ax2.grid(True, which='both', alpha=0.3)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # Machine precision line
    ax2.axhline(2.2e-16, color='gray', linestyle=':', linewidth=1, alpha=0.7)
    ax2.text(42, 4e-16, 'Machine precision', fontsize=9, color='gray',
             ha='right', va='bottom')

    fig1.suptitle(r'Helmholtz with Robin BC: $-u_{xx} + k^2 u = e^{-x^2}$',
                  fontsize=13, y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    # Save figure 1
    output_file1 = OUTPUT_DIR / 'helmholtz_robin.pdf'
    fig1.savefig(output_file1, bbox_inches='tight', pad_inches=0.05)
    fig1.savefig(output_file1.with_suffix('.png'), bbox_inches='tight',
                 pad_inches=0.05, dpi=300)

    print(f'Figure 1 saved to: {output_file1.resolve()}')
    plt.close(fig1)

    # =========================================================================
    # Figure 2: Condition number vs alpha
    # =========================================================================
    fig2, ax3 = plt.subplots(1, 1, figsize=(6, 4.5))

    alpha_range = np.logspace(-2, 2, 80)
    N_cond_values = [16, 32, 64]
    colors_cond = [TEAL, CORAL, PURPLE]

    for N_c, color in zip(N_cond_values, colors_cond):
        cond_numbers = []
        for alpha_val in alpha_range:
            D2, D, x = cheb_second_derivative_matrix(N_c)
            I_mat = np.eye(N_c + 1)
            L = -D2 + k**2 * I_mat

            # Apply BCs
            L[0, :] = 0.0
            L[0, 0] = 1.0
            L[N_c, :] = D[N_c, :] - alpha_val * I_mat[N_c, :]

            cond_numbers.append(np.linalg.cond(L))

        ax3.loglog(alpha_range, cond_numbers, '-', color=color, linewidth=1.5,
                   label=f'$N = {N_c}$')

    ax3.set_xlabel(r'$\alpha$')
    ax3.set_ylabel(r'$\mathrm{cond}(L)$')
    ax3.set_title(r'Condition Number vs Robin Parameter $\alpha$', fontsize=12)
    ax3.legend(loc='best', fontsize=10)
    ax3.grid(True, which='both', alpha=0.3)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    plt.tight_layout()

    # Save figure 2
    output_file2 = OUTPUT_DIR / 'robin_conditioning.pdf'
    fig2.savefig(output_file2, bbox_inches='tight', pad_inches=0.05)
    fig2.savefig(output_file2.with_suffix('.png'), bbox_inches='tight',
                 pad_inches=0.05, dpi=300)

    print(f'Figure 2 saved to: {output_file2.resolve()}')
    plt.close(fig2)

    # Print convergence table
    print(f"\nConvergence Table: Helmholtz with Robin BC (alpha = {alpha_conv})")
    print("-" * 40)
    print(f"{'N':>6} {'Max Error':>14}")
    print("-" * 40)
    for N_val, err in zip(N_values, errors):
        print(f"{N_val:>6d} {err:>14.2e}")
    print("-" * 40)


if __name__ == '__main__':
    main()
