#!/usr/bin/env python3
"""
bvp_2d_poisson.py

Solves the 2D Poisson equation using tensor product Chebyshev methods:

    u_xx + u_yy = f(x,y),  (x,y) ∈ (-1,1)²,  u = 0 on boundary

Test problem: f(x,y) = -2π² sin(πx) sin(πy)
Exact solution: u(x,y) = sin(πx) sin(πy)

This demonstrates tensor product grids and Kronecker product operators.

This script generates Figures 7.7, 7.8, and 7.9 for Chapter 7.

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
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path

# Import Chebyshev matrix function from Chapter 6
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / 'ch06_chebyshev_differentiation'))
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
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch07' / 'python'


def rhs_function(X, Y):
    """f(x,y) = -2π² sin(πx) sin(πy)."""
    return -2 * np.pi**2 * np.sin(np.pi * X) * np.sin(np.pi * Y)


def exact_solution(X, Y):
    """u(x,y) = sin(πx) sin(πy)."""
    return np.sin(np.pi * X) * np.sin(np.pi * Y)


def solve_2d_poisson(N):
    """
    Solve 2D Poisson equation u_xx + u_yy = f with u = 0 on boundary.

    Uses Kronecker product formulation:
        L = I ⊗ D² + D² ⊗ I

    Parameters
    ----------
    N : int
        Number of intervals in each direction

    Returns
    -------
    X, Y : ndarray
        2D meshgrid of Chebyshev points
    U : ndarray
        Solution on the grid
    """
    D2, D, x = cheb_second_derivative_matrix(N)

    # Interior points only
    n_int = N - 1  # Number of interior points
    D2_int = D2[1:N, 1:N]
    x_int = x[1:N]

    # Build 2D Laplacian using Kronecker products
    # L = I ⊗ D² + D² ⊗ I
    I = np.eye(n_int)
    L = np.kron(I, D2_int) + np.kron(D2_int, I)

    # Create meshgrid for interior points
    X_int, Y_int = np.meshgrid(x_int, x_int)

    # Right-hand side (flattened in column-major order for Kronecker convention)
    F = rhs_function(X_int, Y_int)
    f_vec = F.flatten(order='F')  # Column-major to match Kronecker ordering

    # Solve the linear system
    u_vec = np.linalg.solve(L, f_vec)

    # Reshape solution
    U_int = u_vec.reshape((n_int, n_int), order='F')

    # Embed in full grid with boundary conditions
    U = np.zeros((N + 1, N + 1))
    U[1:N, 1:N] = U_int

    # Create full meshgrid
    X, Y = np.meshgrid(x, x)

    return X, Y, U, L


def main():
    """Create 2D Poisson figures."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    N = 16

    # Solve problem
    X, Y, U, L = solve_2d_poisson(N)
    U_exact = exact_solution(X, Y)
    error = np.max(np.abs(U - U_exact))

    # ==========================================================================
    # Figure 7.7: Tensor product grid
    # ==========================================================================
    fig1, ax1 = plt.subplots(1, 1, figsize=(6, 6))

    # Plot grid points
    ax1.plot(X.flatten(), Y.flatten(), 'o', color=TEAL, markersize=5,
             markeredgecolor='white', markeredgewidth=0.5)

    # Add grid lines
    for i in range(N + 1):
        ax1.plot(X[i, :], Y[i, :], '-', color=SKY, linewidth=0.5, alpha=0.6)
        ax1.plot(X[:, i], Y[:, i], '-', color=SKY, linewidth=0.5, alpha=0.6)

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$y$')
    ax1.set_title(f'Chebyshev Tensor Product Grid ($N = {N}$)', fontsize=12)
    ax1.set_xlim(-1.1, 1.1)
    ax1.set_ylim(-1.1, 1.1)
    ax1.set_aspect('equal')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Add annotation about clustering
    ax1.annotate('Grid clusters\nnear boundaries',
                 xy=(0.9, 0.9), xytext=(0.3, 0.5),
                 fontsize=9, color='#666666',
                 arrowprops=dict(arrowstyle='->', color='#666666', alpha=0.6))

    plt.tight_layout()

    output_file = OUTPUT_DIR / 'tensor_grid.pdf'
    fig1.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig1.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                 pad_inches=0.05, dpi=300)
    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig1)

    # ==========================================================================
    # Figure 7.8: 2D Poisson solution
    # ==========================================================================
    fig2 = plt.figure(figsize=(11, 5))

    # 3D surface plot
    ax2a = fig2.add_subplot(1, 2, 1, projection='3d')

    # Create finer grid for plotting
    x_fine = np.linspace(-1, 1, 100)
    X_fine, Y_fine = np.meshgrid(x_fine, x_fine)
    U_exact_fine = exact_solution(X_fine, Y_fine)

    surf = ax2a.plot_surface(X_fine, Y_fine, U_exact_fine, cmap='viridis',
                              alpha=0.8, linewidth=0, antialiased=True)
    ax2a.scatter(X.flatten(), Y.flatten(), U.flatten(), c=CORAL, s=10,
                 depthshade=False, label='Spectral')

    ax2a.set_xlabel(r'$x$')
    ax2a.set_ylabel(r'$y$')
    ax2a.set_zlabel(r'$u(x,y)$')
    ax2a.set_title(f'Solution (max error: {error:.2e})', fontsize=11)
    ax2a.view_init(elev=25, azim=45)

    # Contour plot
    ax2b = fig2.add_subplot(1, 2, 2)
    levels = np.linspace(-1, 1, 21)
    cs = ax2b.contourf(X_fine, Y_fine, U_exact_fine, levels=levels, cmap='viridis')
    ax2b.contour(X_fine, Y_fine, U_exact_fine, levels=levels, colors='white',
                 linewidths=0.5, alpha=0.5)

    # Mark numerical points
    ax2b.plot(X.flatten(), Y.flatten(), 'o', color=CORAL, markersize=3,
              alpha=0.6, markeredgecolor='none')

    ax2b.set_xlabel(r'$x$')
    ax2b.set_ylabel(r'$y$')
    ax2b.set_title('Contour Plot with Grid Points', fontsize=11)
    ax2b.set_aspect('equal')
    plt.colorbar(cs, ax=ax2b, shrink=0.8, label=r'$u(x,y)$')

    fig2.suptitle(r'2D Poisson: $\nabla^2 u = -2\pi^2 \sin(\pi x)\sin(\pi y)$',
                  fontsize=13, y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    output_file = OUTPUT_DIR / 'poisson_2d.pdf'
    fig2.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig2.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                 pad_inches=0.05, dpi=300)
    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig2)

    # ==========================================================================
    # Figure 7.9: Laplacian sparsity pattern
    # ==========================================================================
    fig3, ax3 = plt.subplots(1, 1, figsize=(6, 6))

    # Show sparsity pattern (actually dense, but shows Kronecker structure)
    n_int = N - 1
    ax3.spy(L, markersize=1, color=NAVY)
    ax3.set_title(f'2D Laplacian Structure ($N = {N}$, size {n_int**2}×{n_int**2})',
                  fontsize=11)
    ax3.set_xlabel('Column index')
    ax3.set_ylabel('Row index')

    # Add annotation
    ax3.text(0.5, -0.08, r'$L = I \otimes D^2 + D^2 \otimes I$',
             transform=ax3.transAxes, fontsize=11, ha='center')

    plt.tight_layout()

    output_file = OUTPUT_DIR / 'laplacian_sparsity.pdf'
    fig3.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig3.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                 pad_inches=0.05, dpi=300)
    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig3)

    # Print convergence table
    print("\nConvergence Table: 2D Poisson Equation")
    print("-" * 50)
    print(f"{'N':>6} {'Grid Size':>12} {'Max Error':>14}")
    print("-" * 50)

    for N_val in [4, 8, 12, 16, 20, 24, 32]:
        X_n, Y_n, U_n, _ = solve_2d_poisson(N_val)
        U_exact_n = exact_solution(X_n, Y_n)
        err = np.max(np.abs(U_n - U_exact_n))
        grid_size = (N_val + 1)**2
        print(f"{N_val:>6d} {grid_size:>12d} {err:>14.2e}")

    print("-" * 50)


if __name__ == '__main__':
    main()
