#!/usr/bin/env python3
"""
ho_plate_eigenmodes.py
Chapter 14: Higher-Order Spectral Methods

Eigenmodes of a clamped square plate -- biharmonic eigenvalue problem:

    Delta^2 u = lambda u,   (x, y) in (-1, 1)^2
    u = 0,  du/dn = 0       on the boundary

The clamped plate problem models the vibrations of a thin elastic plate
with all four edges rigidly fixed.  The biharmonic operator Delta^2 is
the composition of two Laplacians.

Discretization follows Trefethen's Program 39 (Spectral Methods in MATLAB):
  1. Build 1D clamped biharmonic operator D4 via a polynomial trick
     that automatically enforces u = du/dn = 0 on the boundary.
  2. Assemble the 2D operator via Kronecker products:
     L = I x D4 + D4 x I + 2 * D2 x D2
  3. Solve the eigenvalue problem L v = lambda v.
  4. Display the first 25 eigenmodes as nodal-line contour plots.

Generates Figure 14.4: Eigenmodes of a clamped square plate.

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""

import numpy as np
import scipy.linalg
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
    implicitly through the "polynomial trick" (Trefethen, Program 39):
    represent u(x) = (1 - x^2) * p(x), so that u and u' vanish at x = +-1.

    The fourth derivative of u = (1 - x^2) p is:
        u'''' = (1 - x^2) p'''' - 8x p''' - 12 p''

    so the effective operator on p at the interior nodes is:
        D4 = [diag(1 - x^2) D^4 - 8 diag(x) D^3 - 12 D^2] S

    where S = diag(1/(1 - x_j^2)) inverts the weight (1 - x^2) at
    interior points (and has zeros at the boundary).

    Parameters
    ----------
    N : int
        Number of Chebyshev intervals (grid has N+1 points).

    Returns
    -------
    D4 : ndarray, shape (N-1, N-1)
        Clamped biharmonic operator restricted to interior points.
    D2int : ndarray, shape (N-1, N-1)
        Standard second derivative matrix restricted to interior points.
    x : ndarray, shape (N+1,)
        Full Chebyshev grid (for reference / plotting).
    """
    D, x = cheb_matrix(N)

    # Powers of D
    D2 = D @ D
    D3 = D2 @ D
    D4full = D3 @ D

    # Diagonal weight matrices
    x2 = x ** 2
    S = np.diag(np.concatenate(([0.0],
                                 1.0 / (1.0 - x2[1:-1]),
                                 [0.0])))

    # Build D4 on full grid, then restrict to interior
    D4mat = (np.diag(1.0 - x2) @ D4full
             - 8.0 * np.diag(x) @ D3
             - 12.0 * D2) @ S

    # Restrict to interior points (indices 1 to N-1)
    D4 = D4mat[1:-1, 1:-1]
    D2int = D2[1:-1, 1:-1]

    return D4, D2int, x


# -----------------------------------------------------------------------------
# 2D biharmonic eigenvalue solver
# -----------------------------------------------------------------------------
def solve_plate_eigenmodes(N, n_modes=25):
    """
    Solve the clamped plate eigenvalue problem on [-1,1]^2.

    Assembles the 2D biharmonic operator via Kronecker products:
        L = I x D4 + D4 x I + 2 * (D2int x D2int)

    and computes the eigenvalues/eigenvectors of L.

    Parameters
    ----------
    N : int
        Number of Chebyshev intervals per direction.
    n_modes : int
        Number of eigenmodes to return (smallest eigenvalues).

    Returns
    -------
    Lam : ndarray, shape (n_modes,)
        Normalized eigenvalues sqrt(lambda_k / lambda_1).
    modes : list of ndarray
        Eigenvectors reshaped to (N+1, N+1) grids (with zero boundaries).
    x : ndarray, shape (N+1,)
        Full Chebyshev grid.
    """
    D4, D2int, x = build_clamped_biharmonic(N)

    M = N - 1  # number of interior points per direction
    I = np.eye(M)

    # 2D biharmonic via Kronecker products
    L = (np.kron(I, D4) + np.kron(D4, I)
         + 2.0 * np.kron(D2int, D2int))

    # Solve eigenvalue problem
    eigenvalues, eigenvectors = scipy.linalg.eig(L)

    # Eigenvalues should be real and positive for the biharmonic
    eigenvalues = np.real(eigenvalues)

    # Sort by ascending eigenvalue (smallest first)
    idx = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Keep only the first n_modes positive eigenvalues
    pos_mask = eigenvalues > 0
    eigenvalues = eigenvalues[pos_mask]
    eigenvectors = eigenvectors[:, pos_mask]

    eigenvalues = eigenvalues[:n_modes]
    eigenvectors = eigenvectors[:, :n_modes]

    # Normalize: Lam_k = sqrt(lambda_k / lambda_1)
    Lam = np.sqrt(eigenvalues / eigenvalues[0])

    # Reshape eigenvectors to 2D grids with zero-padded boundaries
    modes = []
    for k in range(n_modes):
        v = np.real(eigenvectors[:, k])
        V = v.reshape((M, M))

        # Pad with zeros on all boundaries
        U = np.zeros((N + 1, N + 1))
        U[1:-1, 1:-1] = V
        modes.append(U)

    return Lam, modes, x


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
def plot_eigenmodes(Lam, modes, x, n_modes=25):
    """
    Plot the first 25 eigenmodes as nodal-line contour plots in a 5x5 grid.

    Parameters
    ----------
    Lam : ndarray
        Normalized eigenvalues sqrt(lambda_k / lambda_1).
    modes : list of ndarray
        Eigenvector grids, each of shape (N+1, N+1).
    x : ndarray
        Chebyshev grid points.
    n_modes : int
        Number of modes to plot (should be 25 for 5x5 grid).

    Returns
    -------
    fig : matplotlib Figure
    """
    nrows, ncols = 5, 5
    fig, axes = plt.subplots(nrows, ncols, figsize=(10, 10))

    XX, YY = np.meshgrid(x, x)

    for k in range(n_modes):
        row = k // ncols
        col = k % ncols
        ax = axes[row, col]

        U = modes[k]

        # Plot nodal lines (contour at level 0)
        ax.contour(XX, YY, U, levels=[0], colors=NAVY, linewidths=0.8)

        # Draw domain boundary
        ax.plot([-1, 1, 1, -1, -1], [-1, -1, 1, 1, -1],
                color='k', linewidth=1.0)

        # Label with normalized eigenvalue
        ax.set_title(r'$%.4f$' % Lam[k], fontsize=9, pad=3)

        # Clean up axes
        ax.set_xlim(-1.05, 1.05)
        ax.set_ylim(-1.05, 1.05)
        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle(
        r'Eigenmodes of a Clamped Square Plate: '
        r'$\Delta^2 u = \lambda\, u$, '
        r'$u = \partial_n u = 0$ on $\partial\Omega$',
        fontsize=13, y=0.98
    )

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    return fig


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    """Compute and plot eigenmodes of a clamped square plate."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Discretization parameter (following Trefethen, Program 39)
    N = 17
    n_modes = 25

    print("Eigenmodes of a Clamped Square Plate")
    print("=" * 55)
    print(f"  Chebyshev grid:     N = {N}  ({N + 1} points per direction)")
    print(f"  Interior points:    {N - 1} x {N - 1} = {(N - 1) ** 2}")
    print(f"  Eigenvalue problem: {(N - 1) ** 2} x {(N - 1) ** 2} matrix")
    print(f"  Modes to compute:   {n_modes}")
    print()

    # Solve
    Lam, modes, x = solve_plate_eigenmodes(N, n_modes=n_modes)

    # Print eigenvalue table
    print("  Normalized eigenvalues  sqrt(lambda_k / lambda_1):")
    print("  " + "-" * 50)
    for k in range(n_modes):
        print(f"    Mode {k + 1:2d}:  {Lam[k]:.6f}")
    print("  " + "-" * 50)
    print()

    # Plot
    fig = plot_eigenmodes(Lam, modes, x, n_modes=n_modes)

    # Save
    output_file = OUTPUT_DIR / 'plate_eigenmodes.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f"  Figure saved to: {output_file.resolve()}")
    plt.close(fig)

    print("=" * 55)


if __name__ == '__main__':
    main()
