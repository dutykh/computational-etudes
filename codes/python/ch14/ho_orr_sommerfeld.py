#!/usr/bin/env python3
"""
ho_orr_sommerfeld.py
Chapter 14: Higher-Order Boundary Value Problems

Orr-Sommerfeld eigenvalue problem for plane Poiseuille flow at the
critical Reynolds number R = 5772, wavenumber alpha = 1.  The eigenvalues
form the famous Y-shaped spectrum in the complex plane.

The Orr-Sommerfeld equation governs linear stability of parallel shear
flows.  For plane Poiseuille flow U(x) = 1 - x^2 (parabolic profile
between walls at x = +-1), it reads:

    (1/R)(u'''' - 2u'' + u) - 2i u - i(1 - x^2)(u'' - u) = lambda (u'' - u)

with clamped (no-slip) boundary conditions u(+-1) = u'(+-1) = 0.

This is a generalized eigenvalue problem A v = lambda B v where:
    A = (1/R)(D4 - 2*D2 + I) - 2i*I - i*diag(1 - x^2)*(D2 - I)
    B = D2 - I

The fourth derivative matrix D4 is built via the polynomial trick
(Trefethen, Program 38/40) to enforce clamped boundary conditions.
The second derivative D2 uses Dirichlet conditions (boundary rows and
columns removed).

At R = 5772.22, the rightmost eigenvalue is approximately
lambda_max ~ -0.00007819..., placing the flow just barely stable --
this is the critical Reynolds number for plane Poiseuille flow.

Based on Program 40 from Trefethen's "Spectral Methods in MATLAB".

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
# Build clamped fourth-derivative and Dirichlet second-derivative matrices
# -----------------------------------------------------------------------------
def build_orr_sommerfeld_matrices(N):
    """
    Build the interior D4 (clamped) and D2 (Dirichlet) matrices needed
    for the Orr-Sommerfeld generalized eigenvalue problem.

    The clamped fourth-derivative operator is constructed via the polynomial
    trick from Trefethen (Program 38): writing u(x) = (1 - x^2)*q(x)
    automatically enforces u(+-1) = u'(+-1) = 0.  The fourth derivative
    then becomes:

        D4 = [diag(1 - x^2)*D^4 - 8*diag(x)*D^3 - 12*D^2] * S

    where S = diag(1/(1 - x_j^2)) on interior points (zero at boundaries).

    The second derivative D2 uses simple Dirichlet conditions: boundary
    rows and columns are stripped from D^2.

    Parameters
    ----------
    N : int
        Number of Chebyshev intervals (grid has N+1 points).

    Returns
    -------
    D4 : ndarray, shape (N-1, N-1)
        Clamped fourth-derivative operator on interior points.
    D2 : ndarray, shape (N-1, N-1)
        Dirichlet second-derivative operator on interior points.
    x_int : ndarray, shape (N-1,)
        Interior Chebyshev points (excluding boundaries).
    """
    D, x = cheb_matrix(N)

    # Powers of the differentiation matrix
    D2_full = D @ D
    D3_full = D2_full @ D
    D4_full = D3_full @ D

    # --- Clamped D4 via polynomial trick ---
    # S maps v -> q on interior points, zeros at boundaries
    s = np.zeros(N + 1)
    s[1:N] = 1.0 / (1.0 - x[1:N]**2)
    S = np.diag(s)

    # Full fourth-order operator with the polynomial trick
    D4mat = (np.diag(1.0 - x**2) @ D4_full
             - 8.0 * np.diag(x) @ D3_full
             - 12.0 * D2_full) @ S

    # Restrict to interior points
    D4 = D4mat[1:N, 1:N]

    # --- Dirichlet D2: simply strip boundary rows/columns ---
    D2 = D2_full[1:N, 1:N]

    # Interior grid points
    x_int = x[1:N]

    return D4, D2, x_int


# -----------------------------------------------------------------------------
# Solve the Orr-Sommerfeld eigenvalue problem
# -----------------------------------------------------------------------------
def orr_sommerfeld_eigenvalues(N, R=5772.0):
    """
    Compute the Orr-Sommerfeld eigenvalues for plane Poiseuille flow.

    Solves the generalized eigenvalue problem A v = lambda B v where:
        A = (1/R)(D4 - 2*D2 + I) - 2i*I - i*diag(1 - x^2)*(D2 - I)
        B = D2 - I

    Parameters
    ----------
    N : int
        Number of Chebyshev intervals.
    R : float
        Reynolds number (default: 5772, the critical value).

    Returns
    -------
    eigenvalues : ndarray, complex
        Orr-Sommerfeld eigenvalues.
    lam_max : complex
        The rightmost eigenvalue (largest real part).
    """
    D4, D2, x_int = build_orr_sommerfeld_matrices(N)

    M = N - 1  # number of interior points
    I_int = np.eye(M)

    # Assemble the A and B matrices
    # A = (1/R)(D4 - 2*D2 + I) - 2i*I - i*diag(1 - x^2)*(D2 - I)
    A = ((D4 - 2.0 * D2 + I_int) / R
         - 2j * I_int
         - 1j * np.diag(1.0 - x_int**2) @ (D2 - I_int))

    # B = D2 - I
    B = D2 - I_int

    # Solve generalized eigenvalue problem
    eigenvalues = scipy.linalg.eigvals(A, B)

    # Find the rightmost eigenvalue (largest real part)
    idx_max = np.argmax(eigenvalues.real)
    lam_max = eigenvalues[idx_max]

    return eigenvalues, lam_max


# -----------------------------------------------------------------------------
# Main routine
# -----------------------------------------------------------------------------
def main():
    """Compute and plot the Orr-Sommerfeld spectrum for several values of N."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    R = 5772.0  # Critical Reynolds number for plane Poiseuille flow
    N_values = [40, 60, 80, 100]

    print("Orr-Sommerfeld Eigenvalue Problem")
    print("=" * 60)
    print(f"  Flow:             Plane Poiseuille  U(x) = 1 - x^2")
    print(f"  Reynolds number:  R = {R:.0f}")
    print(f"  Wavenumber:       alpha = 1")
    print()

    # =========================================================================
    # Compute eigenvalues for each N
    # =========================================================================
    results = {}
    for N in N_values:
        eigenvalues, lam_max = orr_sommerfeld_eigenvalues(N, R)
        results[N] = (eigenvalues, lam_max)
        print(f"  N = {N:3d}:  lambda_max = {lam_max.real:+.10f} "
              f"{lam_max.imag:+.10f}i")

    print()
    print(f"  Expected:  lambda_max ~ -0.00007819...")
    print("=" * 60)

    # =========================================================================
    # Figure: 2x2 subplot grid
    # =========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    for k, N in enumerate(N_values):
        row = k // 2
        col = k % 2
        ax = axes[row, col]

        eigenvalues, lam_max = results[N]

        # Plot eigenvalues as dots in the complex plane
        ax.plot(eigenvalues.real, eigenvalues.imag, '.',
                color=NAVY, markersize=5)

        # Title with N and lambda_max
        ax.set_title(
            rf'$N = {N}$    $\lambda_{{\max}} = {lam_max.real:.8f}$',
            fontsize=11)

        # Axis labels
        ax.set_xlabel(r'$\mathrm{Re}\,c$')
        ax.set_ylabel(r'$\mathrm{Im}\,c$')

        # Axis limits
        ax.set_xlim(-0.8, 0.2)
        ax.set_ylim(-1.0, 0.0)

        # Grid
        ax.grid(True, linestyle='-', alpha=0.2, color='#888888')

        # Mark the imaginary axis (stability boundary)
        ax.axvline(0.0, color='gray', linestyle=':', linewidth=0.6,
                    alpha=0.5)

        # Spine styling
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle(
        r'Computational Etude 14.7: Orr--Sommerfeld Spectrum '
        rf'($R = {R:.0f}$, $\alpha = 1$)',
        fontsize=13, y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    # ---- Save ----
    output_file = OUTPUT_DIR / 'orr_sommerfeld.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f"\nFigure saved to: {output_file.resolve()}")
    plt.close(fig)


if __name__ == '__main__':
    main()
