#!/usr/bin/env python3
"""
ho_coupled_comparison.py
Chapter 14: Higher-Order Boundary Value Problems

Direct vs. coupled system comparison (Computational Etude 14.3).

Solves the clamped beam problem:

    u^(4)(x) = e^x,   -1 < x < 1,   u(+-1) = u'(+-1) = 0

using TWO approaches and compares their accuracy and conditioning.

Approach 1 (Direct D^4):
    The polynomial trick from Trefethen's "Spectral Methods in MATLAB"
    (Program 38). Write u(x) = (1 - x^2) * q(x), which automatically
    satisfies u(+-1) = 0 and (via the product rule) u'(+-1) = 0.
    The fourth-order operator is reformulated as:
        D4full = [diag(1 - x^2) * D^4 - 8*diag(x) * D^3 - 12*D^2] * S
    and the interior system D4int * v = f is solved directly.

Approach 2 (Coupled second-order system):
    Introduce w = u'' to split the fourth-order ODE into two coupled
    second-order equations:
        u'' = w,        with u(+-1) = 0
        w'' = e^x,      with u'(+-1) = 0
    Assembled as a 2(N+1) x 2(N+1) block system:

        [ A11  A12 ] [ u ]   [ b1 ]
        [ A21  A22 ] [ w ] = [ b2 ]

    where the boundary rows of the u-block enforce u(+-1) = 0
    (Dirichlet) and the boundary rows of the w-block enforce
    u'(+-1) = 0 via D[0,:]*u = 0 and D[N,:]*u = 0.

The direct approach yields a condition number scaling as O(N^8),
while the coupled system, involving only D^2, scales as O(N^4).

Generates Figure 14.3: convergence and conditioning comparison.

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
# Exact solution (same as Etude 14.1)
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
# Approach 1: Direct D^4 via polynomial trick (same as Etude 14.1)
# -----------------------------------------------------------------------------
def solve_direct_d4(N):
    """
    Solve u^(4) = e^x with clamped BCs using the polynomial trick.

    Writes u(x) = (1 - x^2) * q(x) so that all four BCs are built in.
    The operator is reformulated using the Leibniz rule and restricted
    to interior points.

    Parameters
    ----------
    N : int
        Number of Chebyshev intervals.

    Returns
    -------
    u : ndarray, shape (N+1,)
        Numerical solution at all grid points.
    x : ndarray, shape (N+1,)
        Chebyshev grid.
    D4int : ndarray, shape (N-1, N-1)
        Interior fourth-order operator (for condition number).
    """
    D, x = cheb_matrix(N)

    # Powers of the differentiation matrix
    D2 = D @ D
    D3 = D2 @ D
    D4 = D2 @ D2

    # Diagonal matrices
    X = np.diag(x)
    one_minus_x2 = np.diag(1.0 - x**2)

    # S maps v -> q on interior points, zeros at boundaries
    s = np.zeros(N + 1)
    s[1:N] = 1.0 / (1.0 - x[1:N]**2)
    S = np.diag(s)

    # Full fourth-order operator with the polynomial trick
    D4full = (one_minus_x2 @ D4 - 8.0 * X @ D3 - 12.0 * D2) @ S

    # Extract interior block
    D4int = D4full[1:N, 1:N]

    # Right-hand side at interior points
    f = np.exp(x[1:N])

    # Solve
    v = np.linalg.solve(D4int, f)

    # Reconstruct full solution (boundary values are zero)
    u = np.zeros(N + 1)
    u[1:N] = v

    return u, x, D4int


# -----------------------------------------------------------------------------
# Approach 2: Coupled second-order system
# -----------------------------------------------------------------------------
def solve_coupled_d2(N):
    """
    Solve u^(4) = e^x with clamped BCs via a coupled system of two D^2 equations.

    Introduce w = u'' to obtain:
        u'' = w        (with u(+-1) = 0)
        w'' = e^x      (with u'(+-1) = 0 encoded via the D matrix)

    The full 2(N+1) x 2(N+1) block system is:

        [ A11  A12 ] [ u ]   [ b1 ]
        [ A21  A22 ] [ w ] = [ b2 ]

    Top block (rows for u-equation):
        Interior rows i = 1, ..., N-1:  D^2[i,:] * u - I[i,:] * w = 0
        Boundary rows i = 0, N:         u[0] = 0,  u[N] = 0  (Dirichlet)

    Bottom block (rows for w-equation):
        Interior rows i = 1, ..., N-1:  D^2[i,:] * w = exp(x_i)
        Boundary rows i = 0:            D[0,:] * u = 0  (u'(1) = 0)
        Boundary rows i = N:            D[N,:] * u = 0  (u'(-1) = 0)

    Parameters
    ----------
    N : int
        Number of Chebyshev intervals.

    Returns
    -------
    u : ndarray, shape (N+1,)
        Numerical solution at all grid points.
    x : ndarray, shape (N+1,)
        Chebyshev grid.
    A : ndarray, shape (2*(N+1), 2*(N+1))
        Full block system matrix (for condition number).
    """
    D, x = cheb_matrix(N)
    D2 = D @ D
    n = N + 1                # Number of grid points
    I_mat = np.eye(n)
    Z = np.zeros((n, n))

    # Block matrices
    A11 = D2.copy()          # u'' part
    A12 = -I_mat.copy()      # -w part  (u'' - w = 0)
    A21 = Z.copy()           # coupling for w equation
    A22 = D2.copy()          # w'' part

    # Right-hand sides
    b1 = np.zeros(n)
    b2 = np.exp(x)

    # ----- Boundary conditions for the u-block (Dirichlet: u(+-1) = 0) -----
    # Row 0: u(x_0) = u(1) = 0
    A11[0, :] = 0.0
    A11[0, 0] = 1.0
    A12[0, :] = 0.0
    b1[0] = 0.0

    # Row N: u(x_N) = u(-1) = 0
    A11[N, :] = 0.0
    A11[N, N] = 1.0
    A12[N, :] = 0.0
    b1[N] = 0.0

    # ----- Boundary conditions for the w-block (Neumann: u'(+-1) = 0) -----
    # Row 0 of w-block: u'(x_0) = u'(1) = D[0,:] * u = 0
    A21[0, :] = D[0, :]
    A22[0, :] = 0.0
    b2[0] = 0.0

    # Row N of w-block: u'(x_N) = u'(-1) = D[N,:] * u = 0
    A21[N, :] = D[N, :]
    A22[N, :] = 0.0
    b2[N] = 0.0

    # Assemble the full block system
    A = np.block([[A11, A12],
                  [A21, A22]])
    rhs = np.concatenate([b1, b2])

    # Solve
    sol = np.linalg.solve(A, rhs)
    u = sol[:n]
    # w = sol[n:]  # w = u'' (not needed for output)

    return u, x, A


# -----------------------------------------------------------------------------
# Main routine
# -----------------------------------------------------------------------------
def main():
    """Compare direct D^4 and coupled D^2 approaches with convergence plots."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # =========================================================================
    # Convergence and conditioning study
    # =========================================================================
    N_values = list(range(6, 42, 2))

    errors_direct = []
    errors_coupled = []
    cond_direct = []
    cond_coupled = []

    for N in N_values:
        u_ex = exact_solution

        # --- Direct D^4 approach ---
        u_d, x_d, D4int = solve_direct_d4(N)
        u_exact_d = exact_solution(x_d)
        err_d = np.max(np.abs(u_d - u_exact_d))
        errors_direct.append(max(err_d, 1e-16))
        cond_direct.append(np.linalg.cond(D4int))

        # --- Coupled D^2 approach ---
        u_c, x_c, A_block = solve_coupled_d2(N)
        u_exact_c = exact_solution(x_c)
        err_c = np.max(np.abs(u_c - u_exact_c))
        errors_coupled.append(max(err_c, 1e-16))
        cond_coupled.append(np.linalg.cond(A_block))

    # Print convergence table
    print(f"{'N':>4s}  {'Err (Direct)':>13s}  {'Err (Coupled)':>13s}"
          f"  {'cond(D4int)':>14s}  {'cond(Block)':>14s}")
    print("-" * 66)
    for i, N in enumerate(N_values):
        print(f"{N:4d}  {errors_direct[i]:13.4e}  {errors_coupled[i]:13.4e}"
              f"  {cond_direct[i]:14.4e}  {cond_coupled[i]:14.4e}")
    print("-" * 66)

    # =========================================================================
    # Figure: two-panel comparison
    # =========================================================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # ---- Left panel: max error vs N (semilog) ----
    ax1.semilogy(N_values, errors_direct, 'o-', color=CORAL, linewidth=1.5,
                 markersize=5, markeredgecolor='white', markeredgewidth=0.5,
                 label=r'Direct $D^4$ (polynomial trick)')
    ax1.semilogy(N_values, errors_coupled, 's-', color=TEAL, linewidth=1.5,
                 markersize=5, markeredgecolor='white', markeredgewidth=0.5,
                 label=r'Coupled $D^2$ system')

    ax1.set_xlabel(r'$N$')
    ax1.set_ylabel('Max error')
    ax1.set_title('Convergence comparison', fontsize=11)
    ax1.legend(loc='upper right', fontsize=9)
    ax1.set_xlim(4, 44)
    ax1.set_ylim(1e-16, 1e2)
    ax1.grid(True, which='both', alpha=0.15, color='#888888')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Machine precision line
    ax1.axhline(2.2e-16, color='gray', linestyle=':', linewidth=0.8, alpha=0.7)
    ax1.text(43, 4e-16, r'$\epsilon_{\mathrm{mach}}$', fontsize=9,
             color='gray', ha='right', va='bottom')

    # ---- Right panel: condition number vs N (loglog) ----
    ax2.loglog(N_values, cond_direct, 'o-', color=CORAL, linewidth=1.5,
               markersize=5, markeredgecolor='white', markeredgewidth=0.5,
               label=r'Direct $D^4$: $\mathcal{O}(N^8)$')
    ax2.loglog(N_values, cond_coupled, 's-', color=TEAL, linewidth=1.5,
               markersize=5, markeredgecolor='white', markeredgewidth=0.5,
               label=r'Coupled $D^2$: $\mathcal{O}(N^4)$')

    # Reference slopes
    N_arr = np.array(N_values, dtype=float)
    # Fit in log-log to get the actual scaling
    # Use last half of data for fitting (avoid pre-asymptotic regime)
    n_fit = len(N_values) // 2
    log_N = np.log(N_arr[n_fit:])

    # Direct D^4 slope
    log_cond_d = np.log(np.array(cond_direct[n_fit:]))
    slope_d = np.polyfit(log_N, log_cond_d, 1)[0]

    # Coupled D^2 slope
    log_cond_c = np.log(np.array(cond_coupled[n_fit:]))
    slope_c = np.polyfit(log_N, log_cond_c, 1)[0]

    # Reference lines (visual guides)
    N_ref = np.array([10, 40], dtype=float)
    # Scale reference lines to pass through approximate median of data
    ref_N_mid = N_arr[len(N_arr) // 2]

    # Direct reference: N^8
    cond_d_mid = cond_direct[len(N_values) // 2]
    ax2.loglog(N_ref, cond_d_mid * (N_ref / ref_N_mid)**8,
               '--', color=CORAL, linewidth=0.8, alpha=0.5)

    # Coupled reference: N^4
    cond_c_mid = cond_coupled[len(N_values) // 2]
    ax2.loglog(N_ref, cond_c_mid * (N_ref / ref_N_mid)**4,
               '--', color=TEAL, linewidth=0.8, alpha=0.5)

    ax2.set_xlabel(r'$N$')
    ax2.set_ylabel('Condition number')
    ax2.set_title('Conditioning comparison', fontsize=11)
    ax2.legend(loc='upper left', fontsize=9)
    ax2.grid(True, which='both', alpha=0.15, color='#888888')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # Print fitted slopes
    print(f"\nFitted condition number slopes (log-log):")
    print(f"  Direct D^4:  {slope_d:.2f}  (expected ~8)")
    print(f"  Coupled D^2: {slope_c:.2f}  (expected ~4)")

    fig.suptitle(
        r'Computational Etude 14.3: Direct $D^4$ vs. Coupled $D^2$ System',
        fontsize=13, y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    # ---- Save ----
    output_file = OUTPUT_DIR / 'coupled_comparison.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f"\nFigure saved to: {output_file.resolve()}")
    plt.close(fig)


if __name__ == '__main__':
    main()
