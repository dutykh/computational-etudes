#!/usr/bin/env python3
"""
ho_biharmonic_poisson.py
Chapter 14: Higher-Order Spectral Methods

Biharmonic Poisson BVP on a clamped square plate with manufactured solution:

    Delta^2 u = f(x,y),   (x, y) in (-1, 1)^2
    u = 0,  du/dn = 0      on the boundary

Manufactured exact solution:
    u(x,y) = sin^2(pi*x) * sin^2(pi*y)

which satisfies the clamped BCs automatically.

Right-hand side (obtained by applying Delta^2 to the exact solution):
    f(x,y) = 4*pi^4 * [4*cos(2*pi*x)*cos(2*pi*y)
                        - cos(2*pi*x) - cos(2*pi*y)]

Discretisation follows the 2D polynomial trick (Trefethen, Program 39):
  1. Build 1D clamped biharmonic operator D4 via the polynomial trick.
  2. Assemble the 2D operator via Kronecker products:
     L = I x D4 + D4 x I + 2 * D2int x D2int
  3. Solve the linear system L * u = f.
  4. Compare with the exact solution and measure convergence.

Generates Figure: biharmonic_poisson.pdf
  Left panel:  contourf of numerical solution at N = 24
  Right panel: convergence (max error vs N, semilogy)

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
# 1D clamped biharmonic operator via polynomial trick
# -----------------------------------------------------------------------------
def build_clamped_biharmonic(N):
    """
    Build the 1D clamped biharmonic operator D4 and interior second
    derivative D2int on the Chebyshev grid with N+1 points.

    The clamped boundary conditions u(+-1) = u'(+-1) = 0 are enforced
    implicitly through the "polynomial trick" (Trefethen, Program 39).

    Returns
    -------
    D4 : ndarray, shape (N-1, N-1)
        Clamped biharmonic operator restricted to interior points.
    D2int : ndarray, shape (N-1, N-1)
        Standard second derivative matrix restricted to interior points.
    x : ndarray, shape (N+1,)
        Full Chebyshev grid.
    """
    D, x = cheb_matrix(N)

    D2 = D @ D
    D3 = D2 @ D
    D4full = D3 @ D

    x2 = x ** 2
    S = np.diag(np.concatenate(([0.0],
                                 1.0 / (1.0 - x2[1:-1]),
                                 [0.0])))

    D4mat = (np.diag(1.0 - x2) @ D4full
             - 8.0 * np.diag(x) @ D3
             - 12.0 * D2) @ S

    D4 = D4mat[1:-1, 1:-1]
    D2int = D2[1:-1, 1:-1]

    return D4, D2int, x


# -----------------------------------------------------------------------------
# Exact solution and RHS
# -----------------------------------------------------------------------------
def u_exact(x, y):
    """Manufactured exact solution: sin^2(pi*x) * sin^2(pi*y)."""
    return np.sin(np.pi * x)**2 * np.sin(np.pi * y)**2


def f_rhs(x, y):
    """RHS: Delta^2 u_exact = 4*pi^4*(4*cos(2*pi*x)*cos(2*pi*y)
                                       - cos(2*pi*x) - cos(2*pi*y))."""
    c2x = np.cos(2.0 * np.pi * x)
    c2y = np.cos(2.0 * np.pi * y)
    return 4.0 * np.pi**4 * (4.0 * c2x * c2y - c2x - c2y)


# -----------------------------------------------------------------------------
# Solver
# -----------------------------------------------------------------------------
def solve_biharmonic_poisson(N):
    """
    Solve the biharmonic Poisson BVP on [-1,1]^2 with clamped BCs.

    Returns
    -------
    U : ndarray, shape (N+1, N+1)
        Numerical solution on the full grid (with zero boundaries).
    x : ndarray, shape (N+1,)
        Full Chebyshev grid.
    err : float
        Maximum absolute error against the manufactured solution.
    """
    D4, D2int, x = build_clamped_biharmonic(N)

    M = N - 1  # number of interior points per direction
    I_M = np.eye(M)

    # 2D biharmonic via Kronecker products
    L = (np.kron(I_M, D4) + np.kron(D4, I_M)
         + 2.0 * np.kron(D2int, D2int))

    # Interior grid
    xi = x[1:-1]
    XX, YY = np.meshgrid(xi, xi)

    # RHS vector
    f_vec = f_rhs(XX, YY).ravel()

    # Solve
    u_vec = np.linalg.solve(L, f_vec)

    # Exact solution at interior points
    u_ex_vec = u_exact(XX, YY).ravel()
    err = np.max(np.abs(u_vec - u_ex_vec))

    # Reshape to full grid with zero boundaries
    U = np.zeros((N + 1, N + 1))
    U[1:-1, 1:-1] = u_vec.reshape((M, M))

    return U, x, err


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    """Solve the biharmonic Poisson BVP and plot solution + convergence."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ----- Convergence study -----
    N_values = list(range(6, 26, 2))
    errors = []

    print("Biharmonic Poisson BVP with Manufactured Solution")
    print("=" * 55)
    print("  u(x,y) = sin^2(pi*x) * sin^2(pi*y)")
    print("  f(x,y) = 4*pi^4*(4*cos(2*pi*x)*cos(2*pi*y)")
    print("                    - cos(2*pi*x) - cos(2*pi*y))")
    print()
    print("  Convergence study:")
    print("  " + "-" * 40)
    print(f"  {'N':>4s}   {'Interior pts':>12s}   {'Max error':>12s}")
    print("  " + "-" * 40)

    U_plot = None
    x_plot = None

    for N in N_values:
        U, x, err = solve_biharmonic_poisson(N)
        errors.append(err)
        print(f"  {N:4d}   {(N-1)**2:12d}   {err:12.4e}")

        # Save solution at largest N for plotting
        if N == N_values[-1]:
            U_plot = U
            x_plot = x

    print("  " + "-" * 40)
    print()

    # ----- Figure: two panels -----
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # Left panel: contourf of solution at largest N
    XX, YY = np.meshgrid(x_plot, x_plot)
    max_val = np.max(np.abs(U_plot))
    levels = np.linspace(-max_val, max_val, 21)
    cf = ax1.contourf(XX, YY, U_plot, levels=levels, cmap='RdBu_r', alpha=0.4)
    ax1.contour(XX, YY, U_plot, levels=8, colors=NAVY, linewidths=0.6)
    ax1.plot([-1, 1, 1, -1, -1], [-1, -1, 1, 1, -1],
             color='k', linewidth=1.0)
    fig.colorbar(cf, ax=ax1, shrink=0.85, pad=0.02)
    ax1.set_xlim(-1.05, 1.05)
    ax1.set_ylim(-1.05, 1.05)
    ax1.set_aspect('equal')
    ax1.set_xticks([-1, 0, 1])
    ax1.set_yticks([-1, 0, 1])
    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$y$')
    ax1.set_title(r'Numerical solution ($N = %d$)' % N_values[-1])

    # Right panel: convergence
    ax2.semilogy(N_values, errors, 'o-', color=CORAL, markersize=5,
                 linewidth=1.5, label=r'$\|u - u_{\mathrm{exact}}\|_\infty$')
    ax2.axhline(y=2.2e-16, color='gray', linestyle=':', linewidth=0.8,
                label=r'Machine $\varepsilon$')
    ax2.set_xlabel(r'$N$')
    ax2.set_ylabel(r'Max error')
    ax2.set_title('Spectral convergence')
    ax2.legend(loc='upper right', framealpha=0.9)
    ax2.grid(True, which='both', alpha=0.3)
    ax2.set_xticks(N_values)

    fig.suptitle(
        r'Biharmonic Poisson BVP: $\Delta^2 u = f$, '
        r'$u = \partial_n u = 0$ on $\partial\Omega$',
        fontsize=13, y=1.02
    )

    plt.tight_layout()

    # Save
    output_file = OUTPUT_DIR / 'biharmonic_poisson.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f"  Figure saved to: {output_file.resolve()}")
    plt.close(fig)

    print("=" * 55)


if __name__ == '__main__':
    main()
