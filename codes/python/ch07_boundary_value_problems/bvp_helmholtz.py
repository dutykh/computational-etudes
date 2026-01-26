#!/usr/bin/env python3
"""
bvp_helmholtz.py

Solves the Helmholtz equation:

    u_xx + u_yy + k² u = f(x,y),  (x,y) ∈ (-1,1)²,  u = 0 on boundary

with localized Gaussian forcing:
    f(x,y) = exp(-20[(x-0.3)² + (y+0.4)²])

Uses k = 7 to demonstrate near-resonance behavior with the (2,4) eigenmode.
(Eigenvalue for mode (m,n) is k² = (π/2)²(m² + n²), so (2,4) gives k ≈ 7.02)

This script generates Figure 7.10 for Chapter 7.

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


def forcing_function(X, Y, x0=0.3, y0=-0.4, sigma=20.0):
    """
    Localized Gaussian forcing:
    f(x,y) = exp(-σ[(x-x0)² + (y-y0)²])
    """
    return np.exp(-sigma * ((X - x0)**2 + (Y - y0)**2))


def solve_helmholtz(N, k):
    """
    Solve Helmholtz equation u_xx + u_yy + k²u = f with u = 0 on boundary.

    Parameters
    ----------
    N : int
        Number of intervals in each direction
    k : float
        Wave number

    Returns
    -------
    X, Y : ndarray
        2D meshgrid of Chebyshev points
    U : ndarray
        Solution on the grid
    """
    D2, D, x = cheb_second_derivative_matrix(N)

    # Interior points only
    n_int = N - 1
    D2_int = D2[1:N, 1:N]
    x_int = x[1:N]

    # Build 2D Laplacian + k²I using Kronecker products
    # L = I ⊗ D² + D² ⊗ I + k²(I ⊗ I)
    I = np.eye(n_int)
    L = np.kron(I, D2_int) + np.kron(D2_int, I) + k**2 * np.eye(n_int**2)

    # Create meshgrid for interior points
    X_int, Y_int = np.meshgrid(x_int, x_int)

    # Right-hand side (flattened in column-major order)
    F = forcing_function(X_int, Y_int)
    f_vec = F.flatten(order='F')

    # Solve the linear system
    u_vec = np.linalg.solve(L, f_vec)

    # Reshape solution
    U_int = u_vec.reshape((n_int, n_int), order='F')

    # Embed in full grid with boundary conditions
    U = np.zeros((N + 1, N + 1))
    U[1:N, 1:N] = U_int

    # Create full meshgrid
    X, Y = np.meshgrid(x, x)

    return X, Y, U


def main():
    """Create Helmholtz equation figure."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    N = 32
    k = 7.0  # Near resonance with (2,4) mode

    # Theoretical resonance wavenumbers
    print("Resonance wavenumbers for first few modes:")
    print("-" * 50)
    print(f"{'Mode (m,n)':>15} {'k_mn':>12}")
    print("-" * 50)
    for m in range(1, 6):
        for n in range(1, 6):
            k_mn = (np.pi / 2) * np.sqrt(m**2 + n**2)
            if k_mn < 10:
                print(f"{f'({m},{n})':>15} {k_mn:>12.4f}")
    print("-" * 50)

    # Solve Helmholtz equation
    X, Y, U = solve_helmholtz(N, k)

    # Create figure
    fig = plt.figure(figsize=(12, 5))

    # Panel 1: 3D surface
    ax1 = fig.add_subplot(1, 2, 1, projection='3d')

    surf = ax1.plot_surface(X, Y, U, cmap='RdBu_r',
                            alpha=0.9, linewidth=0, antialiased=True,
                            rstride=1, cstride=1)

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$y$')
    ax1.set_zlabel(r'$u(x,y)$')
    ax1.set_title(f'Solution (k = {k})', fontsize=11)
    ax1.view_init(elev=25, azim=45)

    # Panel 2: Contour plot
    ax2 = fig.add_subplot(1, 2, 2)

    # Fine interpolation for smoother contours
    from scipy import interpolate
    x_fine = np.linspace(-1, 1, 200)
    y_fine = np.linspace(-1, 1, 200)
    X_fine, Y_fine = np.meshgrid(x_fine, y_fine)

    # Sort grid for interpolation (Chebyshev points are in decreasing order)
    x_sorted = np.sort(X[0, :])
    y_sorted = np.sort(Y[:, 0])
    sort_x_idx = np.argsort(X[0, :])
    sort_y_idx = np.argsort(Y[:, 0])
    U_sorted = U[np.ix_(sort_y_idx, sort_x_idx)]

    # Interpolate solution
    interp = interpolate.RectBivariateSpline(x_sorted, y_sorted, U_sorted.T)
    U_fine = interp(x_fine, y_fine).T

    # Contour plot
    vmax = np.max(np.abs(U_fine))
    levels = np.linspace(-vmax, vmax, 31)
    cs = ax2.contourf(X_fine, Y_fine, U_fine, levels=levels, cmap='RdBu_r')
    ax2.contour(X_fine, Y_fine, U_fine, levels=levels[::2], colors='k',
                linewidths=0.3, alpha=0.5)

    # Mark forcing location
    ax2.plot(0.3, -0.4, '*', color=ORANGE, markersize=15,
             markeredgecolor='white', markeredgewidth=1, label='Forcing center')

    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r'$y$')
    ax2.set_title('Contour Plot', fontsize=11)
    ax2.set_aspect('equal')
    ax2.legend(loc='upper right', fontsize=9)
    plt.colorbar(cs, ax=ax2, shrink=0.8, label=r'$u(x,y)$')

    # Main title
    fig.suptitle(r'Helmholtz Equation: $\nabla^2 u + k^2 u = f(x,y)$' +
                 f' (k = {k}, near (2,4) resonance at k = 7.02)',
                 fontsize=12, y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    # Save figure
    output_file = OUTPUT_DIR / 'helmholtz.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'\nFigure saved to: {output_file.resolve()}')
    plt.close(fig)

    # Compare solutions at different k values
    print(f"\nSolution amplitude at forcing location for different k:")
    print("-" * 50)
    print(f"{'k':>8} {'|u| at center':>16} {'Near mode':>15}")
    print("-" * 50)

    # Find forcing point in grid
    idx_x = np.argmin(np.abs(X[0, :] - 0.3))
    idx_y = np.argmin(np.abs(Y[:, 0] - (-0.4)))

    for k_test in [3.0, 5.0, 7.0, 7.02, 8.0, 10.0]:
        _, _, U_test = solve_helmholtz(N, k_test)
        u_at_force = np.abs(U_test[idx_y, idx_x])

        # Find nearest eigenvalue
        nearest_mode = ""
        min_dist = float('inf')
        for m in range(1, 10):
            for n in range(1, 10):
                k_mn = (np.pi / 2) * np.sqrt(m**2 + n**2)
                if abs(k_mn - k_test) < min_dist:
                    min_dist = abs(k_mn - k_test)
                    nearest_mode = f"({m},{n})"

        print(f"{k_test:>8.2f} {u_at_force:>16.6f} {nearest_mode:>15}")

    print("-" * 50)
    print("\nNote: Amplitude increases dramatically near resonance values.")


if __name__ == '__main__':
    main()
