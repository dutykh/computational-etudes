#!/usr/bin/env python3
"""
bc_vibrating_string.py
Chapter 13: Boundary Conditions in Spectral Methods

Vibrating string with a free end -- eigenvalue problem:

    -u_xx = lambda u,   x in [0, 1]
    u(0) = 0            (Dirichlet, fixed end)
    u'(1) = 0           (Neumann, free end)

Exact eigenvalues:    lambda_n = ((2n - 1) pi / 2)^2,  n = 1, 2, ...
Exact eigenfunctions: u_n(x) = sin((2n - 1) pi x / 2)

The Chebyshev grid on [-1, 1] is mapped to [0, 1] via x_phys = (1 + x_cheb)/2,
so that x_cheb = -1 maps to x_phys = 0 and x_cheb = 1 maps to x_phys = 1.
The derivative matrix is scaled accordingly: D_phys = 2 D_cheb.

Boundary conditions are imposed via boundary bordering:
    Row 0  (x_cheb = 1, x_phys = 1):  D_phys[0,:] (Neumann u'(1) = 0)
    Row N  (x_cheb = -1, x_phys = 0): identity row (Dirichlet u(0) = 0)

Generates Figure 13.7: eigenmodes and eigenvalue errors.

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
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch13' / 'python'


# -----------------------------------------------------------------------------
# Exact solutions
# -----------------------------------------------------------------------------
def exact_eigenvalue(n):
    """Exact eigenvalue lambda_n = ((2n-1) pi / 2)^2 for n = 1, 2, ..."""
    return ((2 * n - 1) * np.pi / 2) ** 2


def exact_eigenfunction(x, n):
    """Exact eigenfunction u_n(x) = sin((2n-1) pi x / 2)."""
    return np.sin((2 * n - 1) * np.pi * x / 2)


# -----------------------------------------------------------------------------
# Solver
# -----------------------------------------------------------------------------
def solve_vibrating_string(N):
    """
    Solve -u_xx = lambda u on [0,1] with u(0)=0, u'(1)=0.

    Uses Chebyshev collocation with boundary bordering.  The eigenvalue
    problem is formulated as a generalised eigenproblem A u = lambda B u
    so that the boundary rows (which do not involve lambda) can be
    encoded by setting the corresponding rows of B to zero.

    Layout (Chebyshev ordering: x_0 = 1 -> x_N = -1):
        - Row 0  (x_phys = 1): Neumann  u'(1) = 0    -> A[0,:] = D, B[0,:] = 0
        - Row N  (x_phys = 0): Dirichlet u(0) = 0     -> A[N,:] = e_N, B[N,:] = 0
        - Interior rows:       -u'' = lambda u         -> A = -D2, B = I

    Parameters
    ----------
    N : int
        Number of Chebyshev intervals.

    Returns
    -------
    eigenvalues : ndarray
        Computed eigenvalues (sorted, positive, real).
    eigenvectors : ndarray, shape (N+1, ?)
        Corresponding eigenvectors on the physical grid.
    x_phys : ndarray, shape (N+1,)
        Physical grid points on [0, 1].
    """
    n = N + 1

    # Chebyshev grid and derivative on [-1, 1]
    D_cheb, x_cheb = cheb_matrix(N)

    # Map to [0, 1]:  x_phys = (1 + x_cheb) / 2
    x_phys = (1.0 + x_cheb) / 2.0
    D_phys = 2.0 * D_cheb          # d/dx_phys = 2 * d/dx_cheb
    D2_phys = D_phys @ D_phys      # second derivative

    # Build generalised eigenproblem  A u = lambda B u
    A = -D2_phys.copy()
    B = np.eye(n)

    # Neumann BC at x_phys = 1 (row 0): u'(1) = 0
    A[0, :] = D_phys[0, :]
    B[0, :] = 0.0

    # Dirichlet BC at x_phys = 0 (row N): u(0) = 0
    A[N, :] = 0.0
    A[N, N] = 1.0
    B[N, :] = 0.0

    # Solve the generalised eigenvalue problem
    eigenvalues, eigenvectors = scipy.linalg.eig(A, B)

    # Extract real parts (eigenvalues should be real for this self-adjoint problem)
    eigenvalues = eigenvalues.real
    eigenvectors = eigenvectors.real

    # Select finite, positive eigenvalues (boundary rows yield inf/nan)
    good_mask = np.isfinite(eigenvalues) & (eigenvalues > 1e-6)
    eigenvalues = eigenvalues[good_mask]
    eigenvectors = eigenvectors[:, good_mask]

    # Sort by eigenvalue
    idx = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    return eigenvalues, eigenvectors, x_phys


def main():
    """Create vibrating string eigenvalue figure."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    N = 36

    # Solve
    eigenvalues, eigenvectors, x_phys = solve_vibrating_string(N)

    # Fine grid for exact eigenfunctions
    x_fine = np.linspace(0, 1, 500)

    # Number of modes to compare
    n_modes_plot = 6
    n_modes_err = min(len(eigenvalues), N // 2)

    # Colors for the 6 modes
    mode_colors = [NAVY, CORAL, TEAL, PURPLE, ORANGE, SKY]

    # =========================================================================
    # Figure: 2 panels
    # =========================================================================
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # ---- Left panel: first 6 eigenmodes ----
    ax1 = axes[0]

    for k in range(n_modes_plot):
        n_mode = k + 1  # mode number (1-indexed)
        lam_exact = exact_eigenvalue(n_mode)

        # Numerical eigenfunction
        u_num = eigenvectors[:, k]

        # Exact eigenfunction on the spectral grid (for sign alignment)
        u_exact_grid = exact_eigenfunction(x_phys, n_mode)

        # Align sign: if inner product is negative, flip
        if np.dot(u_num, u_exact_grid) < 0:
            u_num = -u_num

        # Normalise to unit amplitude
        u_num = u_num / np.max(np.abs(u_num))

        # Exact eigenfunction on fine grid
        u_exact_fine = exact_eigenfunction(x_fine, n_mode)
        u_exact_fine = u_exact_fine / np.max(np.abs(u_exact_fine))

        color = mode_colors[k]
        label_num = f'$n = {n_mode}$' if k < 3 else None
        label_exact = 'Exact' if k == 0 else None

        ax1.plot(x_fine, u_exact_fine, '--', color=color, linewidth=1.0,
                 alpha=0.6, label=label_exact)
        ax1.plot(x_phys, u_num, 'o', color=color, markersize=4,
                 markeredgecolor='white', markeredgewidth=0.3,
                 label=label_num)

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$u_n(x)$')
    ax1.set_title(f'First {n_modes_plot} Eigenmodes ($N = {N}$)', fontsize=11)
    ax1.legend(loc='lower left', fontsize=8, ncol=2, framealpha=0.9)
    ax1.set_xlim(-0.02, 1.02)
    ax1.set_ylim(-1.3, 1.3)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.axhline(0, color='#888888', linewidth=0.5)
    ax1.grid(True, alpha=0.2)

    # ---- Right panel: eigenvalue errors ----
    ax2 = axes[1]

    mode_numbers = np.arange(1, n_modes_err + 1)
    errors = []

    for k in range(n_modes_err):
        n_mode = k + 1
        lam_exact = exact_eigenvalue(n_mode)
        lam_num = eigenvalues[k]
        err = max(np.abs(lam_num - lam_exact), 1e-16)
        errors.append(err)

    errors = np.array(errors)

    ax2.semilogy(mode_numbers, errors, 'o', color=CORAL, markersize=5,
                 markeredgecolor='white', markeredgewidth=0.5)

    # Resolution limit line at ~ 2N/3
    n_res = int(2 * N / 3)
    ax2.axvline(n_res, color=PURPLE, linestyle=':', linewidth=1.0, alpha=0.7)
    ax2.text(n_res + 0.5, 1e-2, r'$\approx 2N/3$', fontsize=9,
             color=PURPLE, va='center')

    # Machine precision line
    ax2.axhline(2.2e-16, color='gray', linestyle=':', linewidth=1, alpha=0.7)
    ax2.text(n_modes_err - 0.5, 5e-16, 'Machine\nprecision', fontsize=8,
             color='gray', ha='right', va='bottom')

    ax2.set_xlabel(r'Mode number $n$')
    ax2.set_ylabel(r'$|\lambda_n^{\mathrm{num}} - \lambda_n^{\mathrm{exact}}|$')
    ax2.set_title('Eigenvalue Errors', fontsize=11)
    ax2.set_xlim(0, n_modes_err + 1)
    ax2.set_ylim(1e-16, 1e6)
    ax2.grid(True, which='both', alpha=0.3)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # Main title
    fig.suptitle(r"Vibrating String: $-u'' = \lambda u$, $u(0) = 0$, $u'(1) = 0$"
                 f" ($N = {N}$)",
                 fontsize=13, y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.93])

    # Save figure
    output_file = OUTPUT_DIR / 'vibrating_string.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig)

    # =========================================================================
    # Print eigenvalue table
    # =========================================================================
    print(f"\nEigenvalue Accuracy for N = {N}:")
    print("-" * 60)
    print(f"{'n':>4}  {'lambda_exact':>16}  {'lambda_num':>16}  {'Error':>12}")
    print("-" * 60)
    for k in range(min(20, len(eigenvalues))):
        n_mode = k + 1
        lam_exact = exact_eigenvalue(n_mode)
        lam_num = eigenvalues[k]
        err = np.abs(lam_num - lam_exact)
        print(f"{n_mode:>4d}  {lam_exact:>16.8f}  {lam_num:>16.8f}  {err:>12.4e}")
    print("-" * 60)
    print(f"\nResolution limit: ~2N/3 = {2*N//3} well-resolved modes")
    print(f"Total positive eigenvalues found: {len(eigenvalues)}")


if __name__ == '__main__':
    main()
