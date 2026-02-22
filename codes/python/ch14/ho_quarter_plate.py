#!/usr/bin/env python3
"""
ho_quarter_plate.py
Chapter 14: Higher-Order Boundary Value Problems

Eigenmodes of a clamped plate on the quarter-square [0,1]^2 using
symmetry reduction. Neumann (symmetry) conditions at x=0 and y=0
select even-even modes, while clamped conditions are imposed at
x=1 and y=1.

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
# 2D biharmonic eigenvalue solver on the quarter plate
# -----------------------------------------------------------------------------
def solve_quarter_plate(N, n_modes=12):
    """
    Solve the biharmonic eigenvalue problem on the quarter plate [0,1]^2
    with symmetry (Neumann) BCs at x=0, y=0 and clamped BCs at x=1, y=1.

    The Chebyshev grid on [-1,1] is mapped to [0,1] via x_phys = (x+1)/2,
    scaling all derivative operators by a factor of 2.

    The 2D biharmonic operator is assembled via Kronecker products:

        Delta^2 = d^4/dx^4 (x) I  +  2 d^2/dx^2 (x) d^2/dy^2  +  I (x) d^4/dy^4

    Boundary conditions are enforced by boundary bordering (row replacement):
    the rows of the 2D operator corresponding to boundary degrees of freedom
    are replaced by the appropriate BC equations, and the corresponding rows
    in the mass matrix are zeroed.

    For a fourth-order problem, two BCs per boundary segment are required:
      - At x=0 (symmetry):  u_x = 0  and  u_xxx = 0  (odd derivatives vanish)
      - At x=1 (clamped):   u = 0    and  u_x = 0
      - At y=0 (symmetry):  u_y = 0  and  u_yyy = 0
      - At y=1 (clamped):   u = 0    and  u_y = 0

    The resulting generalized eigenvalue problem L v = lambda M v has
    infinite eigenvalues at BC rows (where M is zero), which are filtered
    out to recover the physical spectrum.

    Parameters
    ----------
    N : int
        Number of Chebyshev intervals per direction.
    n_modes : int
        Number of eigenmodes to return (smallest real eigenvalues).

    Returns
    -------
    Lam : ndarray, shape (n_modes,)
        Normalized eigenvalues sqrt(lambda_k / lambda_1).
    modes : list of ndarray
        Eigenvectors reshaped to (N+1, N+1) grids.
    x_phys : ndarray, shape (N+1,)
        Physical grid on [0,1].
    eigenvalues : ndarray, shape (n_modes,)
        Raw eigenvalues lambda_k.
    """
    # 1D Chebyshev operators scaled for [0,1]
    D, x_cheb = cheb_matrix(N)
    D_phys = 2.0 * D                  # scale for half-domain mapping
    D2 = D_phys @ D_phys
    D3 = D2 @ D_phys
    D4 = D2 @ D2

    x_phys = (x_cheb + 1.0) / 2.0     # x_phys[0]=1, x_phys[N]=0

    # Index ordering:
    #   i = 0  -> x_phys = 1  (clamped boundary)
    #   i = N  -> x_phys = 0  (symmetry boundary)

    n = N + 1   # number of grid points per direction
    I1d = np.eye(n)

    # Assemble 2D biharmonic on the full n^2 grid via Kronecker products
    L = (np.kron(I1d, D4)
         + 2.0 * np.kron(D2, D2)
         + np.kron(D4, I1d))

    # Mass matrix for the generalized eigenvalue problem
    M = np.eye(n * n)

    # Grid ordering: u[i + j*n] corresponds to point (x_phys[i], x_phys[j])

    # ----- Boundary bordering: replace rows with BC equations -----
    # We process BCs in four passes, skipping rows already assigned.
    bc_rows = set()

    # Pass 1: Dirichlet u = 0 at the clamped boundaries (x=1 and y=1)
    for j in range(n):
        for i in range(n):
            k = i + j * n

            if i == 0:
                # x_phys = 1 (clamped): u(1, y) = 0
                L[k, :] = 0.0
                L[k, k] = 1.0
                M[k, :] = 0.0
                bc_rows.add(k)

            elif j == 0:
                # y_phys = 1 (clamped): u(x, 1) = 0
                L[k, :] = 0.0
                L[k, k] = 1.0
                M[k, :] = 0.0
                bc_rows.add(k)

    # Pass 2: Normal derivative u_n = 0 at clamped boundaries
    for j in range(n):
        for i in range(n):
            k = i + j * n
            if k in bc_rows:
                continue

            if i == 1:
                # du/dx = 0 at x_phys = 1 (row for i=1, next to clamped edge)
                L[k, :] = 0.0
                for ii in range(n):
                    L[k, ii + j * n] = D_phys[0, ii]
                M[k, :] = 0.0
                bc_rows.add(k)

            elif j == 1:
                # du/dy = 0 at y_phys = 1 (row for j=1)
                L[k, :] = 0.0
                for jj in range(n):
                    L[k, i + jj * n] = D_phys[0, jj]
                M[k, :] = 0.0
                bc_rows.add(k)

    # Pass 3: Neumann u_n = 0 at symmetry boundaries (x=0 and y=0)
    for j in range(n):
        for i in range(n):
            k = i + j * n
            if k in bc_rows:
                continue

            if i == N:
                # du/dx = 0 at x_phys = 0 (symmetry)
                L[k, :] = 0.0
                for ii in range(n):
                    L[k, ii + j * n] = D_phys[N, ii]
                M[k, :] = 0.0
                bc_rows.add(k)

            elif j == N:
                # du/dy = 0 at y_phys = 0 (symmetry)
                L[k, :] = 0.0
                for jj in range(n):
                    L[k, i + jj * n] = D_phys[N, jj]
                M[k, :] = 0.0
                bc_rows.add(k)

    # Pass 4: Third derivative = 0 at symmetry boundaries
    # (for even functions, u_xxx(0) = 0 and u_yyy(0) = 0)
    for j in range(n):
        for i in range(n):
            k = i + j * n
            if k in bc_rows:
                continue

            if i == N - 1:
                # d^3u/dx^3 = 0 at x_phys = 0 (symmetry)
                L[k, :] = 0.0
                for ii in range(n):
                    L[k, ii + j * n] = D3[N, ii]
                M[k, :] = 0.0
                bc_rows.add(k)

            elif j == N - 1:
                # d^3u/dy^3 = 0 at y_phys = 0 (symmetry)
                L[k, :] = 0.0
                for jj in range(n):
                    L[k, i + jj * n] = D3[N, jj]
                M[k, :] = 0.0
                bc_rows.add(k)

    # ----- Solve the generalized eigenvalue problem L v = lambda M v -----
    eigenvalues, eigenvectors = scipy.linalg.eig(L, M)

    # Filter: keep only finite, real, positive eigenvalues
    eigenvalues = np.array(eigenvalues)
    finite_mask = np.isfinite(eigenvalues)
    eigenvalues = eigenvalues[finite_mask]
    eigenvectors = eigenvectors[:, finite_mask]

    # Take real parts (self-adjoint problem)
    eigenvalues_real = np.real(eigenvalues)

    # Discard spurious near-zero and negative eigenvalues
    pos_mask = eigenvalues_real > 1.0
    eigenvalues_real = eigenvalues_real[pos_mask]
    eigenvectors = eigenvectors[:, pos_mask]

    # Sort by ascending eigenvalue
    idx = np.argsort(eigenvalues_real)
    eigenvalues_real = eigenvalues_real[idx]
    eigenvectors = eigenvectors[:, idx]

    # Keep the requested number of modes
    n_avail = min(n_modes, len(eigenvalues_real))
    eigenvalues_real = eigenvalues_real[:n_avail]
    eigenvectors = eigenvectors[:, :n_avail]

    # Normalize: Lam_k = sqrt(lambda_k / lambda_1)
    Lam = np.sqrt(eigenvalues_real / eigenvalues_real[0])

    # Reshape eigenvectors to 2D grids
    modes = []
    for k in range(n_avail):
        v = np.real(eigenvectors[:, k])
        U = v.reshape((n, n), order='F')   # column-major: u[i,j] at (x[i], y[j])
        modes.append(U)

    return Lam, modes, x_phys, eigenvalues_real


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
def plot_eigenmodes(Lam, modes, x_phys, n_modes=12):
    """
    Plot the first 12 even-even eigenmodes on the quarter plate [0,1]^2
    in a 4x3 grid of subplots showing nodal lines (zero contours).

    Parameters
    ----------
    Lam : ndarray
        Normalized eigenvalues sqrt(lambda_k / lambda_1).
    modes : list of ndarray
        Eigenvector grids, each of shape (N+1, N+1).
    x_phys : ndarray
        Physical grid points on [0,1].
    n_modes : int
        Number of modes to plot (should be 12 for 4x3 grid).

    Returns
    -------
    fig : matplotlib Figure
    """
    nrows, ncols = 4, 3
    fig, axes = plt.subplots(nrows, ncols, figsize=(8, 10))

    XX, YY = np.meshgrid(x_phys, x_phys, indexing='ij')

    for k in range(n_modes):
        row = k // ncols
        col = k % ncols
        ax = axes[row, col]

        U = modes[k]

        # Filled contours for the mode shape
        max_val = np.max(np.abs(U))
        levels = np.linspace(-max_val, max_val, 21)
        ax.contourf(XX, YY, U, levels=levels, cmap='RdBu_r', alpha=0.4)

        # Nodal lines (contour at level 0)
        ax.contour(XX, YY, U, levels=[0], colors=NAVY, linewidths=0.9)

        # Draw domain boundary
        ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color='k', linewidth=1.0)

        # Mark symmetry edges with dashed lines
        ax.plot([0, 0], [0, 1], '--', color=TEAL, linewidth=1.0)
        ax.plot([0, 1], [0, 0], '--', color=TEAL, linewidth=1.0)

        # Label with normalized eigenvalue
        ax.set_title(r'$%.4f$' % Lam[k], fontsize=9, pad=3)

        # Clean up axes
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.set_aspect('equal')
        ax.set_xticks([0, 0.5, 1])
        ax.set_yticks([0, 0.5, 1])
        ax.tick_params(labelsize=7)

    fig.suptitle(
        r'Even-Even Eigenmodes of a Clamped Plate (Quarter Domain)'
        '\n'
        r'$\Delta^2 u = \lambda\, u$ on $[0,1]^2$, '
        r'$\partial_n u = 0$ (symmetry), '
        r'$u = \partial_n u = 0$ (clamped)',
        fontsize=11, y=0.99
    )

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    return fig


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    """Compute and plot even-even eigenmodes of a clamped plate via quarter-domain."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Discretization parameter
    N = 13
    n_modes = 12

    print("Quarter-Plate Symmetry Reduction: Even-Even Eigenmodes")
    print("=" * 60)
    print(f"  Domain:             [0, 1]^2  (quarter of [-1, 1]^2)")
    print(f"  Chebyshev grid:     N = {N}  ({N + 1} points per direction)")
    print(f"  System size:        {(N + 1) ** 2} x {(N + 1) ** 2}")
    print(f"  Modes to compute:   {n_modes}")
    print(f"  BCs at x=0, y=0:   Neumann (symmetry)")
    print(f"  BCs at x=1, y=1:   Clamped (u = du/dn = 0)")
    print()

    # Solve
    Lam, modes, x_phys, eigenvalues = solve_quarter_plate(N, n_modes=n_modes)

    # Print eigenvalue table
    print("  Even-even eigenvalues of the clamped square plate:")
    print("  " + "-" * 55)
    print(f"  {'Mode':>4s}  {'lambda_k':>14s}  {'sqrt(lam_k/lam_1)':>18s}")
    print("  " + "-" * 55)
    for k in range(len(Lam)):
        print(f"  {k + 1:4d}  {eigenvalues[k]:14.6f}  {Lam[k]:18.6f}")
    print("  " + "-" * 55)
    print()

    # Plot
    fig = plot_eigenmodes(Lam, modes, x_phys, n_modes=len(Lam))

    # Save
    output_file = OUTPUT_DIR / 'quarter_plate.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f"  Figure saved to: {output_file.resolve()}")
    plt.close(fig)

    print("=" * 60)


if __name__ == '__main__':
    main()
