#!/usr/bin/env python3
"""
bc_laplace_2d.py
Chapter 13: Boundary Conditions in Spectral Methods

Laplace equation with piecewise boundary data on [-1,1]^2:

    u_xx + u_yy = 0,  (x,y) in (-1,1)^2

Boundary conditions:
    u(x, 1)  = sin^4(pi*x)  for -1 < x < 0  (top, left half)
    u(1, y)  = (1/5)*sin(3*pi*y)              (right boundary)
    u = 0    everywhere else on the boundary

Solved via boundary bordering (Method II):
    1. Build full 2D Laplacian L = I x D^2 + D^2 x I
    2. Replace boundary rows with identity rows and set boundary data
    3. Solve the global linear system

Generates Figure 13.5: 2D Laplace with piecewise BCs.

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
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D
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
# Boundary data functions
# -----------------------------------------------------------------------------
def bc_top(x):
    """Top boundary u(x, 1) = sin^4(pi*x) for x < 0, else 0."""
    val = np.zeros_like(x)
    mask = x < 0
    val[mask] = np.sin(np.pi * x[mask]) ** 4
    return val


def bc_right(y):
    """Right boundary u(1, y) = (1/5)*sin(3*pi*y)."""
    return (1.0 / 5.0) * np.sin(3.0 * np.pi * y)


# -----------------------------------------------------------------------------
# Solver
# -----------------------------------------------------------------------------
def solve_laplace_2d(N):
    """
    Solve the 2D Laplace equation with piecewise boundary data.

    Uses boundary bordering (Method II): build the full (N+1)^2 system,
    then replace each boundary row with an identity row and the
    corresponding boundary value.

    Parameters
    ----------
    N : int
        Number of Chebyshev intervals in each direction.

    Returns
    -------
    XX, YY : ndarray, shape (N+1, N+1)
        2D meshgrid of Chebyshev points.
    U : ndarray, shape (N+1, N+1)
        Solution on the tensor product grid.
    """
    D2, D, x = cheb_second_derivative_matrix(N)
    n = N + 1  # number of points per direction

    # Build 2D Laplacian: L = I x D2 + D2 x I
    I = np.eye(n)
    L = np.kron(I, D2) + np.kron(D2, I)

    # Create meshgrid
    XX, YY = np.meshgrid(x, x)
    xx = XX.flatten()
    yy = YY.flatten()

    # Right-hand side (zero for Laplace equation)
    rhs = np.zeros(n * n)

    # Identify boundary nodes and apply boundary conditions
    for k in range(n * n):
        is_boundary = (np.abs(np.abs(xx[k]) - 1.0) < 1e-14 or
                       np.abs(np.abs(yy[k]) - 1.0) < 1e-14)

        if is_boundary:
            # Replace row with identity row
            L[k, :] = 0.0
            L[k, k] = 1.0

            # Set boundary value
            xk, yk = xx[k], yy[k]

            # Top boundary: y = 1
            if np.abs(yk - 1.0) < 1e-14:
                rhs[k] = bc_top(np.array([xk]))[0]
            # Right boundary: x = 1
            elif np.abs(xk - 1.0) < 1e-14:
                rhs[k] = bc_right(np.array([yk]))[0]
            # All other boundaries: u = 0
            else:
                rhs[k] = 0.0

    # Solve the linear system
    u_vec = np.linalg.solve(L, rhs)

    # Reshape to 2D
    U = u_vec.reshape((n, n))

    return XX, YY, U


def main():
    """Create 2D Laplace equation figure with piecewise boundary data."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    N = 24
    XX, YY, U = solve_laplace_2d(N)

    # Custom colormap using book colors
    colors_cmap = [NAVY, TEAL, SKY, '#FFFFFF']
    book_cmap = LinearSegmentedColormap.from_list('book', colors_cmap, N=256)

    # =========================================================================
    # Figure: 3D surface plot
    # =========================================================================
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_surface(XX, YY, U, cmap=book_cmap,
                           alpha=0.9, linewidth=0.3,
                           edgecolor='grey', antialiased=True,
                           rstride=1, cstride=1)

    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_zlabel(r'$u(x,y)$')
    ax.set_title(
        r'Laplace Equation: $\nabla^2 u = 0$ with Piecewise BCs ($N = %d$)'
        % N, fontsize=12
    )
    ax.view_init(elev=25, azim=-55)

    fig.colorbar(surf, ax=ax, shrink=0.6, aspect=15, pad=0.1,
                 label=r'$u(x,y)$')

    plt.tight_layout()

    # Save figure
    output_file = OUTPUT_DIR / 'laplace_2d_mixed.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig)

    # =========================================================================
    # Print summary
    # =========================================================================
    print(f"\n2D Laplace Equation with Piecewise Boundary Data")
    print("=" * 55)
    print(f"  Grid size:          N = {N}  ({(N+1)**2} unknowns)")
    print(f"  Domain:             [-1, 1]^2")
    print(f"  Solution range:     [{U.min():.6f}, {U.max():.6f}]")
    print(f"  BCs applied:")
    print(f"    Top (y=1, x<0):   sin^4(pi*x)")
    print(f"    Right (x=1):      (1/5)*sin(3*pi*y)")
    print(f"    Elsewhere:        u = 0")
    print("=" * 55)

    # Verify boundary conditions
    x = XX[0, :]
    y = YY[:, 0]

    # Top boundary check
    u_top = U[0, :]  # y = x[0] = 1
    bc_top_exact = bc_top(x)
    top_err = np.max(np.abs(u_top - bc_top_exact))
    print(f"\n  BC verification:")
    print(f"    Top boundary error:    {top_err:.2e}")

    # Right boundary check
    u_right = U[:, 0]  # x = x[0] = 1
    bc_right_exact = bc_right(y)
    right_err = np.max(np.abs(u_right - bc_right_exact))
    print(f"    Right boundary error:  {right_err:.2e}")

    # Bottom boundary check
    u_bottom = U[-1, :]  # y = x[N] = -1
    bottom_err = np.max(np.abs(u_bottom))
    print(f"    Bottom boundary error: {bottom_err:.2e}")

    # Left boundary check
    u_left = U[:, -1]  # x = x[N] = -1
    left_err = np.max(np.abs(u_left))
    print(f"    Left boundary error:   {left_err:.2e}")


if __name__ == '__main__':
    main()
