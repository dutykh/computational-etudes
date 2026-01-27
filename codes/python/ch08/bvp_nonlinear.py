#!/usr/bin/env python3
"""
bvp_nonlinear.py

Solves the Bratu equation (nonlinear BVP) using Newton iteration:

    u_xx + λ e^u = 0,  x ∈ (-1, 1),  u(±1) = 0

with λ = 0.5 (subcritical case, unique solution exists).

IMPORTANT: The critical value for [-1, 1] domain is λ_c ≈ 0.878.
For λ > λ_c, no solution exists! We use λ = 0.5 for a stable solution.

This demonstrates that nonlinear BVPs are easily handled with spectral
methods using Newton iteration.

This script generates Figure 7.5 and Table 7.1 for Chapter 7.

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


def solve_bratu_newton(N, lam=0.5, tol=1e-12, max_iter=30):
    """
    Solve the Bratu equation u_xx + λ e^u = 0 with u(±1) = 0 using Newton.

    Newton iteration for F(u) = D²u + λ e^u = 0:
        J δu = -F(u)
        J = D² + λ diag(e^u)
        u_new = u + α δu  (with line search for robustness)

    Parameters
    ----------
    N : int
        Number of intervals
    lam : float
        Parameter λ (default 1.0)
    tol : float
        Convergence tolerance
    max_iter : int
        Maximum iterations

    Returns
    -------
    x : ndarray
        Chebyshev grid points
    u : ndarray
        Solution
    residuals : list
        Convergence history (residual norms)
    n_iter : int
        Number of iterations used
    """
    D2, D, x = cheb_second_derivative_matrix(N)

    # Extract interior system
    D2_int = D2[1:N, 1:N]

    # Initial guess: u = 0
    u = np.zeros(N + 1)

    residuals = []

    for it in range(max_iter):
        # Residual: F(u) = D²u + λ e^u (interior points)
        exp_u = np.exp(u[1:N])
        F = D2_int @ u[1:N] + lam * exp_u

        res_norm = np.max(np.abs(F))
        residuals.append(res_norm)

        if res_norm < tol:
            return x, u, residuals, it + 1

        # Jacobian: J = D² + λ diag(e^u)
        J = D2_int + lam * np.diag(exp_u)

        # Newton step
        delta_u = np.linalg.solve(J, -F)

        # Backtracking line search to ensure residual decreases
        alpha = 1.0
        u_trial = u.copy()
        for _ in range(10):  # Max 10 backtracking steps
            u_trial[1:N] = u[1:N] + alpha * delta_u
            exp_u_trial = np.exp(u_trial[1:N])
            F_trial = D2_int @ u_trial[1:N] + lam * exp_u_trial
            res_trial = np.max(np.abs(F_trial))

            if res_trial < res_norm:
                break
            alpha *= 0.5

        # Update interior points with damped step
        u[1:N] = u[1:N] + alpha * delta_u

    return x, u, residuals, max_iter


def solve_bratu_fixedpoint(N, lam=0.5, tol=1e-10, max_iter=100):
    """
    Solve Bratu equation using fixed-point iteration for comparison.

    Fixed-point iteration: D² u_new = -λ exp(u_old)
    """
    D2, D, x = cheb_second_derivative_matrix(N)
    D2_int = D2[1:N, 1:N]

    # Initial guess
    u = np.zeros(N + 1)

    residuals = []

    for it in range(max_iter):
        # RHS
        rhs = -lam * np.exp(u[1:N])

        # Solve
        u_new_int = np.linalg.solve(D2_int, rhs)

        # Compute change
        change = np.max(np.abs(u_new_int - u[1:N]))
        residuals.append(change)

        u[1:N] = u_new_int

        if change < tol:
            return x, u, residuals, it + 1

    return x, u, residuals, max_iter


def main():
    """Create nonlinear BVP figure."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Panel 1: Solution for various N
    ax1 = axes[0]

    N_values = [8, 16, 32]
    colors = [CORAL, TEAL, PURPLE]

    for N, color in zip(N_values, colors):
        x, u, residuals, n_iter = solve_bratu_newton(N)
        ax1.plot(x, u, 'o-', color=color, linewidth=1.2, markersize=4,
                 markeredgecolor='white', markeredgewidth=0.3,
                 label=f'$N = {N}$ ({n_iter} iter)')

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$u(x)$')
    ax1.set_title(r'Bratu Equation: $u_{xx} + e^u = 0$', fontsize=11)
    ax1.legend(loc='lower center', fontsize=9, ncol=3)
    ax1.set_xlim(-1.05, 1.05)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
    ax1.xaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')

    # Panel 2: Convergence history
    ax2 = axes[1]

    # Compare fixed-point vs Newton for N=16
    N = 16
    _, _, res_fp, n_fp = solve_bratu_fixedpoint(N, tol=1e-14, max_iter=50)
    _, _, res_newton, n_newton = solve_bratu_newton(N, tol=1e-14, max_iter=20)

    ax2.semilogy(range(1, len(res_fp) + 1), res_fp, 'o-', color=CORAL,
                 linewidth=1.5, markersize=5, label='Fixed-point')
    ax2.semilogy(range(1, len(res_newton) + 1), res_newton, 's-', color=TEAL,
                 linewidth=1.5, markersize=5, label='Newton')

    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('Residual norm')
    ax2.set_title(f'Convergence History ($N = {N}$)', fontsize=11)
    ax2.legend(loc='upper right', fontsize=9)
    ax2.grid(True, which='both', alpha=0.3)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_ylim(1e-16, 10)

    # Add tolerance line
    ax2.axhline(1e-12, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax2.text(2, 3e-12, 'tol', fontsize=9, color='gray')

    # Main title
    fig.suptitle(r'Nonlinear BVP: Bratu Equation with $\lambda = 0.5$',
                 fontsize=13, y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    # Save figure
    output_file = OUTPUT_DIR / 'bratu.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig)

    # Generate convergence table
    print("\n" + "=" * 60)
    print("Table 7.1: Bratu Equation Convergence (Newton)")
    print("=" * 60)
    print(f"{'N':>6} {'Iterations':>12} {'u(0)':>16}")
    print("-" * 60)

    for N in [4, 8, 16, 32, 64]:
        x, u, _, n_iter = solve_bratu_newton(N)
        # Find index closest to x=0
        idx_0 = np.argmin(np.abs(x))
        print(f"{N:>6d} {n_iter:>12d} {u[idx_0]:>16.10f}")

    print("=" * 60)
    print("\nNote: Newton iteration converges quadratically (in 5-8 iterations)")
    print("while fixed-point converges linearly (in ~18 iterations).")


if __name__ == '__main__':
    main()
