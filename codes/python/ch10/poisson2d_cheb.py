#!/usr/bin/env python3
"""
poisson2d_cheb.py
Chapter 10: Spectral PDE Solvers with Chebyshev and Fourier Grids

2D Poisson equation solver on Chebyshev tensor grid.

PDE: -Lap(u) = f, -1 < x, y < 1
BCs: u = 0 on boundary (Dirichlet)

Uses manufactured solution:
    u(x, y) = (1 - x^2)(1 - y^2) cos(pi*x/2) cos(pi*y/2)

which automatically satisfies homogeneous Dirichlet BCs.

Generates Figure 10.7: Exact, numerical, and error comparison.

Author: Dr. Denys Dutykh
Date: February 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pathlib import Path

from chebfft import cheb_matrix

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

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
NAVY = '#142D6E'

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch10' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'poisson2d_solution.pdf'


def exact_solution(xx, yy):
    """Manufactured exact solution that satisfies u = 0 on boundary."""
    return (1 - xx**2) * (1 - yy**2) * np.cos(np.pi * xx / 2) * np.cos(np.pi * yy / 2)


def exact_laplacian(xx, yy):
    """Analytical Laplacian of the exact solution.

    u = g(x)*h(y) with g(x) = (1-x^2)*cos(pi*x/2), h(y) = (1-y^2)*cos(pi*y/2).
    Lap(u) = g''(x)*h(y) + g(x)*h''(y).
    """
    p = np.pi
    cx, cy = np.cos(p * xx / 2), np.cos(p * yy / 2)
    sx, sy = np.sin(p * xx / 2), np.sin(p * yy / 2)
    gx = (1 - xx**2) * cx
    hy = (1 - yy**2) * cy
    gpp = -(2 + p**2 / 4 * (1 - xx**2)) * cx + 2 * p * xx * sx
    hpp = -(2 + p**2 / 4 * (1 - yy**2)) * cy + 2 * p * yy * sy
    return gpp * hy + gx * hpp


def poisson2d_cheb(N=32):
    """
    Solve 2D Poisson equation with manufactured solution.

    The discrete system is:
        A @ u_vec = f_vec
    where A is the Kronecker sum of D2_int with itself.

    Parameters:
        N : int
            Number of Chebyshev intervals in each direction

    Returns:
        xx, yy : arrays of shape (N+1, N+1)
            Meshgrid of Chebyshev points
        U : array of shape (N+1, N+1)
            Numerical solution
        U_exact : array of shape (N+1, N+1)
            Exact (manufactured) solution
        error : float
            Maximum absolute error
    """
    # Build differentiation matrices
    D, x = cheb_matrix(N)
    D2 = D @ D
    xx, yy = np.meshgrid(x, x)

    # Exact solution
    U_exact = exact_solution(xx, yy)

    # Interior indices
    ii = slice(1, N)
    D2_int = D2[np.ix_(range(1, N), range(1, N))]
    I_int = np.eye(N - 1)

    # Kronecker sum Laplacian for interior (represents Δ)
    L = np.kron(D2_int, I_int) + np.kron(I_int, D2_int)

    # For -Δu = f, we have -L @ u = f, so L @ u = -f
    # Or equivalently, use A = -L so that A @ u = f
    A = -L

    # RHS from analytical Laplacian (not discrete, to test convergence)
    f_full = -exact_laplacian(xx, yy)
    f_int = f_full[np.ix_(range(1, N), range(1, N))]

    # Solve A @ u = f where A = -Δ
    u_vec = np.linalg.solve(A, f_int.flatten())

    # Reconstruct full solution
    U = np.zeros((N + 1, N + 1))
    U[np.ix_(range(1, N), range(1, N))] = u_vec.reshape(N - 1, N - 1)

    # Compute error
    error = np.max(np.abs(U - U_exact))

    return xx, yy, U, U_exact, error


def main():
    # Solve Poisson equation
    N = 32
    xx, yy, U, U_exact, error = poisson2d_cheb(N=N)

    print(f"N = {N}, Max error = {error:.2e}")

    # Create figure with 3 panels
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # Common levels for first two plots
    vmin, vmax = U_exact.min(), U_exact.max()
    levels = np.linspace(vmin, vmax, 30)

    # Panel 1: Exact solution
    ax = axes[0]
    cf = ax.contourf(xx, yy, U_exact, levels=levels, cmap='RdBu_r')
    ax.contour(xx, yy, U_exact, levels=10, colors='k', linewidths=0.3, alpha=0.5)
    ax.set_xlabel(r'$x$', fontsize=11)
    ax.set_ylabel(r'$y$', fontsize=11)
    ax.set_title('Exact Solution', fontsize=12)
    ax.set_aspect('equal')
    cbar = plt.colorbar(cf, ax=ax, shrink=0.8)

    # Panel 2: Numerical solution
    ax = axes[1]
    cf = ax.contourf(xx, yy, U, levels=levels, cmap='RdBu_r')
    ax.contour(xx, yy, U, levels=10, colors='k', linewidths=0.3, alpha=0.5)
    ax.set_xlabel(r'$x$', fontsize=11)
    ax.set_ylabel(r'$y$', fontsize=11)
    ax.set_title('Numerical Solution', fontsize=12)
    ax.set_aspect('equal')
    cbar = plt.colorbar(cf, ax=ax, shrink=0.8)

    # Panel 3: Error (log scale)
    ax = axes[2]
    err = np.abs(U - U_exact)
    # Use log scale for error, but handle zeros
    err_plot = np.log10(np.maximum(err, 1e-16))
    cf = ax.contourf(xx, yy, err_plot, levels=20, cmap='viridis')
    ax.set_xlabel(r'$x$', fontsize=11)
    ax.set_ylabel(r'$y$', fontsize=11)
    ax.set_title(f'$\\log_{{10}}$ Error (max = {error:.2e})', fontsize=12)
    ax.set_aspect('equal')
    cbar = plt.colorbar(cf, ax=ax, shrink=0.8)
    cbar.set_label(r'$\log_{10}|u - u_{\rm exact}|$', fontsize=10)

    plt.suptitle('2D Poisson Equation: Spectral Solution', fontsize=14, y=1.02)
    plt.tight_layout()

    # Save figure
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_FILE, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(OUTPUT_FILE.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {OUTPUT_FILE.resolve()}')

    # Convergence study
    print("\nConvergence Study:")
    print("-" * 40)
    print(f"{'N':>6}  {'Max Error':>14}  {'Ratio':>10}")
    print("-" * 40)

    prev_error = None
    for N in [8, 12, 16, 20, 24, 28, 32, 40, 48]:
        _, _, _, _, err = poisson2d_cheb(N=N)
        if prev_error is not None and err > 1e-15:
            ratio = prev_error / err
        else:
            ratio = float('nan')
        print(f"{N:6d}  {err:14.6e}  {ratio:10.2f}")
        prev_error = err

    plt.close(fig)


if __name__ == '__main__':
    main()
