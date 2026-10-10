#!/usr/bin/env python3
"""
ho_beam_eigenmodes.py
Chapter 14: Higher-Order Differential Equations

Computational Etude 14.2: Vibration Modes of a Clamped Beam.

Solves the eigenvalue problem for a clamped-clamped Euler--Bernoulli beam:

    u''''(x) = lambda * u(x),   -1 < x < 1,
    u(+-1) = u'(+-1) = 0.

The four boundary conditions (two Dirichlet and two Neumann) are imposed
using the boundary-bordering technique:

    1. Build the Chebyshev differentiation matrix D on N+1 points.
    2. Compute D4 = D^4 (fourth-derivative matrix).
    3. Extract the interior block D4_int = D4[1:N, 1:N] (Dirichlet u=0).
    4. Replace the first and last rows of D4_int with the derivative
       conditions D[0, 1:N] and D[N, 1:N] to enforce u'(+-1) = 0.
    5. Solve the resulting generalised eigenvalue problem A v = lambda B v
       where B has zero rows at the replaced positions.

Reference eigenvalues are computed from the transcendental equations that
arise from the exact solution of the clamped beam problem on [-1, 1].
Separating into symmetric and antisymmetric modes gives:

    Symmetric modes:       tan(beta) + tanh(beta) = 0
    Antisymmetric modes:   tan(beta) - tanh(beta) = 0

In both cases lambda = beta^4.

Parameters: N = 30 (default for figure generation).

Generates Figures 14.2a-b:
    (a) First 6 eigenmodes of the clamped beam.
    (b) Relative eigenvalue error vs mode number (semilog scale).

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Last modified: 2026-02-22

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""

import numpy as np
import scipy.linalg
import scipy.optimize
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
# Reference eigenvalues from the transcendental equation
# -----------------------------------------------------------------------------
def compute_reference_eigenvalues(n_modes=20):
    """
    Compute reference eigenvalues of u'''' = lambda u on [-1, 1]
    with clamped BCs u(+-1) = u'(+-1) = 0.

    The general solution u(x) = A cos(beta x) + B sin(beta x)
    + C cosh(beta x) + D sinh(beta x) with beta = lambda^{1/4}.
    Applying the four BCs and separating into symmetric/antisymmetric
    modes yields:

        Symmetric modes:       tan(beta) + tanh(beta) = 0
        Antisymmetric modes:   tan(beta) - tanh(beta) = 0

    The eigenvalue is lambda = beta^4.

    Parameters
    ----------
    n_modes : int
        Number of eigenvalues to compute.

    Returns
    -------
    lambdas : ndarray, shape (n_modes,)
        Exact eigenvalues sorted in ascending order.
    """
    # Symmetric mode equation: tan(beta) + tanh(beta) = 0
    def sym_eq(beta):
        return np.tan(beta) + np.tanh(beta)

    # Antisymmetric mode equation: tan(beta) - tanh(beta) = 0
    def antisym_eq(beta):
        return np.tan(beta) - np.tanh(beta)

    betas = []

    # Both equations have roots in intervals between the poles of tan,
    # which occur at (2k+1)*pi/2.  We search each continuity interval
    # (k*pi - pi/2 + eps, k*pi + pi/2 - eps) for roots of both equations.
    eps = 0.01
    n_search = n_modes + 10  # search enough intervals

    for k in range(n_search):
        # Continuity interval for tan: (k*pi - pi/2, k*pi + pi/2)
        lo = k * np.pi - np.pi / 2.0 + eps
        hi = k * np.pi + np.pi / 2.0 - eps
        if lo < eps:
            lo = eps  # beta must be positive

        # Check for symmetric mode root
        try:
            flo = sym_eq(lo)
            fhi = sym_eq(hi)
            if flo * fhi < 0:
                root = scipy.optimize.brentq(sym_eq, lo, hi, xtol=1e-15)
                if root > 0.5 and all(abs(root - b) > 0.01 for b in betas):
                    betas.append(root)
        except (ValueError, RuntimeError):
            pass

        # Check for antisymmetric mode root
        try:
            flo = antisym_eq(lo)
            fhi = antisym_eq(hi)
            if flo * fhi < 0:
                root = scipy.optimize.brentq(antisym_eq, lo, hi, xtol=1e-15)
                if root > 0.5 and all(abs(root - b) > 0.01 for b in betas):
                    betas.append(root)
        except (ValueError, RuntimeError):
            pass

    # Sort betas and compute eigenvalues lambda = beta^4
    betas = np.sort(betas)
    lambdas = betas[:n_modes] ** 4

    return lambdas


# -----------------------------------------------------------------------------
# Biharmonic clamped matrix (boundary-bordering technique)
# -----------------------------------------------------------------------------
def biharmonic_clamped_matrix(N):
    """
    Build the interior fourth-derivative operator with clamped BCs.

    For the problem u''''(x) = f(x) on [-1,1] with u(+-1) = u'(+-1) = 0,
    we use the boundary-bordering technique:

        1. Build D^4 on the (N+1)-point Chebyshev grid.
        2. Extract interior rows/columns [1:N, 1:N] for Dirichlet BCs.
        3. Replace the first row with D[0, 1:N]  (enforces u'(1) = 0).
        4. Replace the last row with D[N, 1:N]   (enforces u'(-1) = 0).

    This yields a generalised eigenvalue problem A v = lambda B v,
    where B has zero rows at positions 0 and N-2 (the replaced rows).

    Parameters
    ----------
    N : int
        Number of Chebyshev intervals.

    Returns
    -------
    A : ndarray, shape (N-1, N-1)
        Modified fourth-derivative matrix with boundary conditions.
    B : ndarray, shape (N-1, N-1)
        Mass matrix (identity with zeroed rows for BC equations).
    x_int : ndarray, shape (N-1,)
        Interior Chebyshev grid points.
    """
    D, x = cheb_matrix(N)

    # Fourth-derivative matrix
    D2 = D @ D
    D4 = D2 @ D2

    # Extract interior block (Dirichlet u(+-1) = 0)
    A = D4[1:N, 1:N].copy()

    # Replace first row with derivative condition at x = +1 (index 0)
    # u'(+1) = D[0, :] @ u = D[0, 1:N] @ u_int = 0
    A[0, :] = D[0, 1:N]

    # Replace last row with derivative condition at x = -1 (index N)
    # u'(-1) = D[N, :] @ u = D[N, 1:N] @ u_int = 0
    A[-1, :] = D[N, 1:N]

    # Mass matrix: identity except zero rows where BCs were imposed
    B = np.eye(N - 1)
    B[0, :] = 0.0
    B[-1, :] = 0.0

    x_int = x[1:N]

    return A, B, x_int


# -----------------------------------------------------------------------------
# Eigenvalue solver
# -----------------------------------------------------------------------------
def solve_beam_eigenmodes(N):
    """
    Solve the clamped beam eigenvalue problem u'''' = lambda u.

    Parameters
    ----------
    N : int
        Number of Chebyshev intervals.

    Returns
    -------
    eigenvalues : ndarray
        Computed eigenvalues sorted by real part (ascending).
    eigenvectors : ndarray
        Corresponding eigenvectors (columns), defined at interior
        Chebyshev points x_1, ..., x_{N-1}.
    x_int : ndarray
        Interior grid points.
    """
    A, B, x_int = biharmonic_clamped_matrix(N)

    # Solve generalised eigenvalue problem A v = lambda B v
    eigenvalues, eigenvectors = scipy.linalg.eig(A, B)

    # Filter out infinite/NaN eigenvalues (from the singular rows of B)
    finite_mask = np.isfinite(eigenvalues) & (np.abs(eigenvalues) < 1e15)
    eigenvalues = eigenvalues[finite_mask]
    eigenvectors = eigenvectors[:, finite_mask]

    # Take real parts (eigenvalues should be real and positive)
    eigenvalues = eigenvalues.real

    # Sort by ascending real part
    idx = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx].real

    # Keep only positive eigenvalues (physical modes)
    pos_mask = eigenvalues > 0
    eigenvalues = eigenvalues[pos_mask]
    eigenvectors = eigenvectors[:, pos_mask]

    return eigenvalues, eigenvectors, x_int


# -----------------------------------------------------------------------------
# Barycentric interpolation for smooth plotting
# -----------------------------------------------------------------------------
def cheb_interpolate(x_cheb, u_cheb, x_fine):
    """
    Interpolate from Chebyshev nodes to a fine grid via barycentric formula.

    Uses the barycentric weights for Chebyshev-Gauss-Lobatto points:
        w_j = (-1)^j * d_j,  where d_0 = d_N = 1/2, d_j = 1 otherwise.

    Parameters
    ----------
    x_cheb : ndarray, shape (M,)
        Chebyshev nodes (must be sorted ascending or descending).
    u_cheb : ndarray, shape (M,)
        Function values at Chebyshev nodes.
    x_fine : ndarray, shape (P,)
        Target evaluation points.

    Returns
    -------
    u_fine : ndarray, shape (P,)
        Interpolated values at x_fine.
    """
    M = len(x_cheb)
    # Barycentric weights for Chebyshev-Gauss-Lobatto points
    w = np.ones(M)
    w[0] = 0.5
    w[-1] = 0.5
    for j in range(M):
        if j % 2 == 1:
            w[j] = -w[j]

    # Evaluate using the second (true) form of the barycentric formula
    u_fine = np.zeros_like(x_fine)
    for i, xf in enumerate(x_fine):
        diff = xf - x_cheb
        # Check if xf coincides with a node
        exact = np.where(np.abs(diff) < 1e-15)[0]
        if len(exact) > 0:
            u_fine[i] = u_cheb[exact[0]]
        else:
            terms = w / diff
            u_fine[i] = np.dot(terms, u_cheb) / np.sum(terms)

    return u_fine


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
def plot_eigenmodes_and_convergence(output_dir):
    """
    Generate two-panel figure:
      Left:  First 6 eigenmodes of the clamped beam (N = 30).
      Right: Relative eigenvalue error vs mode number (semilog).
    """
    N = 30
    n_display = 6

    # Solve eigenvalue problem
    eigenvalues, eigenvectors, x_int = solve_beam_eigenmodes(N)

    # Full Chebyshev grid (for interpolation, including boundary zeros)
    _, x_full_cheb = cheb_matrix(N)

    # Reference eigenvalues
    n_ref = min(len(eigenvalues), 20)
    lam_exact = compute_reference_eigenvalues(n_ref)

    # Fine grid for smooth plotting
    x_fine = np.linspace(-1, 1, 500)

    # Colors for the 6 modes
    mode_colors = [NAVY, SKY, CORAL, TEAL, PURPLE, ORANGE]

    # ---- Create figure ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # ==== Left panel: Eigenmodes ====
    ax1 = axes[0]

    for k in range(n_display):
        if k >= len(eigenvalues):
            break

        # Extract eigenmode at interior points (x_int is in descending order)
        u_int = eigenvectors[:, k].copy()

        # Build full solution including boundary zeros
        u_full = np.zeros(N + 1)
        u_full[1:N] = u_int  # interior values (same order as x_full_cheb[1:N])

        # Normalise so max|u| = 1
        u_full = u_full / np.max(np.abs(u_full))

        # Fix sign convention: make the first peak (from the left) positive
        # x_full_cheb is descending (1 to -1), so reverse for left-to-right
        x_asc = x_full_cheb[::-1]
        u_asc = u_full[::-1]
        # Find the first local extremum from the left
        first_extremum = u_asc[np.argmax(np.abs(u_asc))]
        if first_extremum < 0:
            u_full = -u_full

        # Interpolate onto fine grid using barycentric formula
        u_fine = cheb_interpolate(x_full_cheb, u_full, x_fine)

        lam_k = eigenvalues[k]
        label = r'Mode %d ($\lambda_{%d} \approx %.1f$)' % (k + 1, k + 1, lam_k)

        ax1.plot(x_fine, u_fine, '-', color=mode_colors[k], linewidth=1.5,
                 label=label)

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$u(x)$')
    ax1.set_title('First %d eigenmodes ($N = %d$)' % (n_display, N))
    ax1.legend(loc='lower left', fontsize=7.5, framealpha=0.9, ncol=1)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.axhline(0, color='#888888', linewidth=0.5)
    ax1.grid(True, alpha=0.2)
    ax1.set_xlim(-1.05, 1.05)

    # ==== Right panel: Eigenvalue convergence ====
    ax2 = axes[1]

    n_compare = min(len(eigenvalues), len(lam_exact))
    modes = np.arange(1, n_compare + 1)
    rel_errors = np.abs(eigenvalues[:n_compare] - lam_exact[:n_compare]) / lam_exact[:n_compare]

    # Replace any zeros with machine epsilon for log plot
    rel_errors = np.maximum(rel_errors, 1e-16)

    ax2.semilogy(modes, rel_errors, 'o-', color=NAVY, markersize=6,
                 markeredgecolor='white', markeredgewidth=0.5, linewidth=1.2)

    ax2.set_xlabel('Mode number $n$')
    ax2.set_ylabel(r'$|\lambda_n^{\mathrm{num}} - \lambda_n^{\mathrm{exact}}| / \lambda_n^{\mathrm{exact}}$')
    ax2.set_title('Relative eigenvalue error ($N = %d$)' % N)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.grid(True, which='both', alpha=0.3)

    fig.suptitle(r"Clamped Beam Eigenmodes: $u''''= \lambda\, u$, "
                 r"$u(\pm 1) = u'(\pm 1) = 0$",
                 fontsize=13, y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.93])

    # Save figure
    output_file = output_dir / 'beam_eigenmodes.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)
    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig)

    return eigenvalues, lam_exact


def main():
    """Compute clamped beam eigenmodes and generate figures."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    N = 30

    print("Vibration Modes of a Clamped Beam")
    print("=" * 60)
    print(f"  Equation: u''''(x) = lambda u(x),  x in (-1, 1)")
    print(f"  BCs: u(+-1) = u'(+-1) = 0")
    print(f"  N = {N} Chebyshev points")
    print()

    # Compute reference eigenvalues
    lam_exact = compute_reference_eigenvalues(20)
    print("  Reference eigenvalues (from transcendental equation):")
    print("  " + "-" * 50)
    print(f"  {'Mode':>6} {'lambda (exact)':>20} {'beta':>14}")
    print("  " + "-" * 50)
    for k in range(min(10, len(lam_exact))):
        beta_k = lam_exact[k] ** 0.25
        print(f"  {k+1:>6d} {lam_exact[k]:>20.6f} {beta_k:>14.6f}")
    print("  " + "-" * 50)
    print()

    # Solve numerically
    eigenvalues, eigenvectors, x_int = solve_beam_eigenmodes(N)

    n_compare = min(len(eigenvalues), len(lam_exact))
    print(f"  Number of computed positive eigenvalues: {len(eigenvalues)}")
    print()

    print("  Comparison of numerical vs exact eigenvalues:")
    print("  " + "-" * 70)
    print(f"  {'Mode':>6} {'lambda (num)':>18} {'lambda (exact)':>18} {'Rel. Error':>14}")
    print("  " + "-" * 70)
    for k in range(n_compare):
        rel_err = abs(eigenvalues[k] - lam_exact[k]) / lam_exact[k]
        print(f"  {k+1:>6d} {eigenvalues[k]:>18.6f} {lam_exact[k]:>18.6f} {rel_err:>14.4e}")
    print("  " + "-" * 70)
    print()

    # Generate figures
    print("Generating figures...")
    plot_eigenmodes_and_convergence(OUTPUT_DIR)

    print("\nDone.")


if __name__ == '__main__':
    main()
