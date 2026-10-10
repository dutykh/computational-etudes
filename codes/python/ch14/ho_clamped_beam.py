#!/usr/bin/env python3
"""
ho_clamped_beam.py
Chapter 14: Higher-Order Boundary Value Problems

Clamped beam under exponential load (Computational Etude 14.1).

Solves the fourth-order BVP (biharmonic-type equation):

    u^(4)(x) = e^x,   -1 < x < 1,   u(+-1) = u'(+-1) = 0

using a Chebyshev spectral method with the polynomial trick from
Trefethen's "Spectral Methods in MATLAB" (Program 38).

The key idea is to write u(x) = (1 - x^2) * q(x), which automatically
satisfies both Dirichlet conditions u(+-1) = 0.  The factor (1 - x^2)
also enforces u'(+-1) = 0 because u'(x) = -2x*q(x) + (1 - x^2)*q'(x),
so u'(+-1) = -2*(+-1)*q(+-1) + 0 = 0 only if q(+-1) = 0.  Instead of
enforcing this directly, the fourth derivative operator is reformulated
via the product rule as:

    D4full = [diag(1 - x^2) * D^4 - 8 * diag(x) * D^3 - 12 * D^2] * S

where S = diag(0, 1/(1 - x_1^2), ..., 1/(1 - x_{N-1}^2), 0) maps v
back to q on interior points and annihilates boundary values.

The exact solution is u(x) = e^x + c_3*x^3 + c_2*x^2 + c_1*x + c_0,
where the constants are determined by the four clamped boundary conditions.

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
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch14' / 'python'


# -----------------------------------------------------------------------------
# Biharmonic solver using Trefethen's polynomial trick (Program 38)
# -----------------------------------------------------------------------------
def biharmonic_clamped_matrix(N):
    """
    Build the fourth-order operator for clamped BCs via the polynomial trick.

    The approach writes u(x) = (1 - x^2) * q(x) and reformulates the
    fourth derivative operator using the product rule.  The resulting
    system is restricted to interior Chebyshev points.

    Parameters
    ----------
    N : int
        Number of Chebyshev intervals (grid has N+1 points).

    Returns
    -------
    D4int : ndarray, shape (N-1, N-1)
        Interior block of the fourth-order operator matrix.
    x : ndarray, shape (N+1,)
        Full Chebyshev grid (including boundary points).
    """
    D, x = cheb_matrix(N)

    # Powers of the differentiation matrix
    D2 = D @ D
    D3 = D2 @ D
    D4 = D2 @ D2

    # Diagonal matrices
    X = np.diag(x)                          # diag(x)
    one_minus_x2 = np.diag(1.0 - x**2)     # diag(1 - x^2)

    # S maps v -> q on interior points, zeros at boundaries
    s = np.zeros(N + 1)
    s[1:N] = 1.0 / (1.0 - x[1:N]**2)
    S = np.diag(s)

    # Full fourth-order operator with the polynomial trick
    # Derived from the Leibniz rule for d^4/dx^4 [(1 - x^2) * q(x)]
    D4full = (one_minus_x2 @ D4 - 8.0 * X @ D3 - 12.0 * D2) @ S

    # Extract interior block (strip first and last rows/columns)
    D4int = D4full[1:N, 1:N]

    return D4int, x


# -----------------------------------------------------------------------------
# Exact solution
# -----------------------------------------------------------------------------
def exact_solution(x):
    """
    Exact solution of u^(4)(x) = e^x with clamped BCs u(+-1) = u'(+-1) = 0.

    The general solution is:

        u(x) = e^x + c_3 * x^3 + c_2 * x^2 + c_1 * x + c_0

    since e^x is a particular solution of u^(4) = e^x (its fourth derivative
    is itself) and x^3, x^2, x, 1 span the nullspace of d^4/dx^4.

    The four constants are determined by the boundary conditions:

        u(1)  = 0  =>  e   + c_3 + c_2 + c_1 + c_0 = 0
        u(-1) = 0  =>  1/e - c_3 + c_2 - c_1 + c_0 = 0
        u'(1) = 0  =>  e   + 3*c_3 + 2*c_2 + c_1    = 0
        u'(-1)= 0  =>  1/e + 3*c_3 - 2*c_2 + c_1    = 0

    Parameters
    ----------
    x : ndarray
        Evaluation points.

    Returns
    -------
    u : ndarray
        Exact solution values at x.
    """
    e_plus = np.exp(1.0)
    e_minus = np.exp(-1.0)

    # 4x4 linear system: A * [c3, c2, c1, c0]^T = b
    A = np.array([
        [ 1.0,  1.0,  1.0, 1.0],      # u(1) = 0
        [-1.0,  1.0, -1.0, 1.0],      # u(-1) = 0
        [ 3.0,  2.0,  1.0, 0.0],      # u'(1) = 0
        [ 3.0, -2.0,  1.0, 0.0],      # u'(-1) = 0
    ])
    b = np.array([-e_plus, -e_minus, -e_plus, -e_minus])

    c = np.linalg.solve(A, b)
    c3, c2, c1, c0 = c

    return np.exp(x) + c3 * x**3 + c2 * x**2 + c1 * x + c0


# -----------------------------------------------------------------------------
# Main routine
# -----------------------------------------------------------------------------
def main():
    """Solve the clamped beam problem, perform convergence study, and plot."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # =========================================================================
    # Solve for N = 15 (demonstration)
    # =========================================================================
    N_demo = 15
    D4int, x = biharmonic_clamped_matrix(N_demo)

    # Right-hand side at interior points
    f = np.exp(x[1:N_demo])

    # Solve the linear system
    v = np.linalg.solve(D4int, f)

    # Reconstruct full solution (boundary values are zero)
    u_num = np.zeros(N_demo + 1)
    u_num[1:N_demo] = v

    # Exact solution on the full grid
    u_ex = exact_solution(x)

    demo_error = np.max(np.abs(u_num - u_ex))
    print(f"N = {N_demo}: max |u_num - u_exact| = {demo_error:.4e}")

    # Fine grid for smooth exact solution curve
    x_fine = np.linspace(-1.0, 1.0, 500)
    u_fine = exact_solution(x_fine)

    # =========================================================================
    # Convergence study
    # =========================================================================
    N_values = list(range(4, 34, 2))       # N = 4, 6, 8, ..., 32
    errors = []
    cond_numbers = []

    for N in N_values:
        D4int_n, x_n = biharmonic_clamped_matrix(N)
        f_n = np.exp(x_n[1:N])

        v_n = np.linalg.solve(D4int_n, f_n)

        u_n = np.zeros(N + 1)
        u_n[1:N] = v_n

        u_ex_n = exact_solution(x_n)
        err = np.max(np.abs(u_n - u_ex_n))
        errors.append(max(err, 1e-16))       # floor at machine precision
        cond_numbers.append(np.linalg.cond(D4int_n))

    # Print convergence table
    print(f"\n{'N':>4s}  {'Max Error':>12s}  {'cond(D4int)':>14s}")
    print("-" * 36)
    for N, err, kappa in zip(N_values, errors, cond_numbers):
        print(f"{N:4d}  {err:12.4e}  {kappa:14.4e}")
    print("-" * 36)

    # =========================================================================
    # Figure: two-panel plot
    # =========================================================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # ---- Left panel: solution for N = 15 ----
    ax1.plot(x_fine, u_fine, '-', color=NAVY, linewidth=1.5,
             label='Exact solution')
    ax1.plot(x, u_num, 'o', color=CORAL, markersize=6,
             markeredgecolor='white', markeredgewidth=0.5,
             label=f'Spectral ($N = {N_demo}$)')

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$u(x)$')
    ax1.set_title(r'Clamped beam: $u^{(4)} = e^x$', fontsize=11)
    ax1.legend(loc='upper left', fontsize=9)
    ax1.set_xlim(-1.05, 1.05)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
    ax1.xaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')

    # ---- Right panel: convergence + condition number ----
    # Primary y-axis: max error (semilog)
    color_err = CORAL
    ax2.semilogy(N_values, errors, 'o-', color=color_err, linewidth=1.5,
                 markersize=5, markeredgecolor='white', markeredgewidth=0.5,
                 label='Max error')
    ax2.set_xlabel(r'$N$')
    ax2.set_ylabel('Max error', color=color_err)
    ax2.tick_params(axis='y', labelcolor=color_err)
    ax2.set_ylim(1e-16, 1e0)
    ax2.set_xlim(2, 34)

    # Machine precision line
    ax2.axhline(2.2e-16, color='gray', linestyle=':', linewidth=0.8,
                alpha=0.7)
    ax2.text(33, 4e-16, r'$\epsilon_{\mathrm{mach}}$', fontsize=9,
             color='gray', ha='right', va='bottom')

    ax2.set_title('Convergence and conditioning', fontsize=11)
    ax2.spines['top'].set_visible(False)
    ax2.grid(True, which='both', alpha=0.15, color='#888888')

    # Secondary y-axis: condition number (log-log via semilogy)
    ax2b = ax2.twinx()
    color_cond = TEAL
    ax2b.semilogy(N_values, cond_numbers, 's--', color=color_cond,
                  linewidth=1.2, markersize=4, markeredgecolor='white',
                  markeredgewidth=0.4, label=r'$\kappa(D_4)$')
    ax2b.set_ylabel(r'$\mathrm{cond}(D_4)$', color=color_cond)
    ax2b.tick_params(axis='y', labelcolor=color_cond)
    ax2b.spines['top'].set_visible(False)

    # Combined legend
    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2b.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc='center left',
               fontsize=9)

    fig.suptitle(
        r'Computational Etude 14.1: Clamped Beam Under Exponential Load',
        fontsize=13, y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    # ---- Save ----
    output_file = OUTPUT_DIR / 'clamped_beam.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f"\nFigure saved to: {output_file.resolve()}")
    plt.close(fig)


if __name__ == '__main__':
    main()
