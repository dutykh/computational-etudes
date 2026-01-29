#!/usr/bin/env python3
"""
harmonic_oscillator.py

Solving the quantum harmonic oscillator eigenvalue problem:
    -u'' + x^2 u = lambda u,  x in R

using the periodic spectral method on a truncated domain [-L, L].

Exact eigenvalues: lambda_n = 2n + 1 for n = 0, 1, 2, ...
Exact eigenfunctions: u_n(x) = H_n(x) * exp(-x^2/2) (Hermite functions)

This example demonstrates "spectral accuracy in action":
    - The eigenfunctions decay like exp(-x^2/2), so truncation is valid
    - The spectral method achieves machine precision for eigenvalues
    - Convergence is super-geometric (faster than any exponential)

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pathlib import Path
from scipy.special import hermite

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
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'
PURPLE = '#9B59B6'
ORANGE = '#E67E22'

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch08' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'harmonic_oscillator.pdf'


# -----------------------------------------------------------------------------
# Periodic spectral second derivative matrix
# -----------------------------------------------------------------------------
def periodic_d2_matrix(N):
    """
    Construct the periodic second derivative matrix for N points.

    Parameters
    ----------
    N : int
        Number of grid points (should be even)

    Returns
    -------
    D2 : ndarray (N, N)
        Second derivative matrix
    """
    h = 2 * np.pi / N

    # Build second derivative matrix using the formula from Trefethen
    # D2[j,j] = -pi^2/(3h^2) - 1/6
    # D2[j,k] = -0.5*(-1)^(j-k) / sin^2(h*(j-k)/2) for j != k
    column = np.zeros(N)
    column[0] = -np.pi**2 / (3 * h**2) - 1.0 / 6.0

    for k in range(1, N):
        column[k] = -0.5 * ((-1)**k) / np.sin(h * k / 2.0)**2

    # Build Toeplitz matrix
    from scipy.linalg import toeplitz
    D2 = toeplitz(column)

    return D2


# -----------------------------------------------------------------------------
# Solve harmonic oscillator
# -----------------------------------------------------------------------------
def solve_harmonic_oscillator(N, L):
    """
    Solve -u'' + x^2 u = lambda u on [-L, L] using periodic spectral method.

    Parameters
    ----------
    N : int
        Number of grid points
    L : float
        Half-width of domain

    Returns
    -------
    eigenvalues : ndarray
        Computed eigenvalues (sorted)
    eigenvectors : ndarray
        Corresponding eigenvectors (columns)
    x : ndarray
        Grid points
    """
    # Set up grid
    h = 2 * np.pi / N
    x_periodic = h * np.arange(N)
    x = L * (x_periodic - np.pi) / np.pi  # Map to [-L, L]

    # Second derivative matrix (with scaling)
    D2 = periodic_d2_matrix(N) * (np.pi / L)**2

    # Potential matrix
    V = np.diag(x**2)

    # Full operator: -D2 + V
    A = -D2 + V

    # Solve eigenvalue problem
    eigenvalues, eigenvectors = np.linalg.eigh(A)

    return eigenvalues, eigenvectors, x


# -----------------------------------------------------------------------------
# Exact Hermite functions (normalized)
# -----------------------------------------------------------------------------
def hermite_function(n, x):
    """
    Compute the n-th Hermite function (normalized eigenfunction of harmonic oscillator).

    psi_n(x) = (1/sqrt(2^n n! sqrt(pi))) * H_n(x) * exp(-x^2/2)

    where H_n is the physicist's Hermite polynomial.
    """
    # Normalization constant
    norm = 1.0 / np.sqrt(2**n * math.factorial(n) * np.sqrt(np.pi))

    # Hermite polynomial H_n(x)
    H_n = hermite(n)

    return norm * H_n(x) * np.exp(-x**2 / 2.0)


# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
def main():
    # Parameters
    L = 8.0  # Domain half-width (large enough for eigenfunction decay)

    # Create figure with subplots
    fig = plt.figure(figsize=(12, 5))

    # -------------------------------------------------------------------------
    # Panel 1: First few eigenfunctions
    # -------------------------------------------------------------------------
    ax1 = fig.add_subplot(1, 2, 1)

    N = 64
    eigenvalues, eigenvectors, x = solve_harmonic_oscillator(N, L)

    colors = [NAVY, CORAL, TEAL, PURPLE, ORANGE]

    for n in range(5):
        # Numerical eigenfunction
        u_num = eigenvectors[:, n]

        # Normalize to match exact eigenfunction
        # (eigenvectors from eigh are normalized, but sign may differ)
        u_exact = hermite_function(n, x)

        # Fix sign
        if np.dot(u_num, u_exact) < 0:
            u_num = -u_num

        # Scale to match (in case of different normalization)
        scale = np.max(np.abs(u_exact)) / np.max(np.abs(u_num))
        u_num = u_num * scale

        # Plot exact (line) and numerical (markers)
        ax1.plot(x, u_exact + 2 * n, '-', color=colors[n], linewidth=1.5,
                 label=f'$n={n}$, $\\lambda={2*n+1}$' if n == 0 else f'$n={n}$')
        ax1.plot(x[::4], (u_num + 2 * n)[::4], 'o', color=colors[n],
                 markersize=4, alpha=0.7)

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$u_n(x)$ (shifted)')
    ax1.set_title(f'Harmonic Oscillator Eigenfunctions ($N={N}$, $L={L}$)')
    ax1.set_xlim(-6, 6)
    ax1.legend(loc='upper right', fontsize=9)

    # Add annotation
    ax1.text(4.5, 7.5, 'Lines: exact\nDots: spectral', fontsize=9,
             ha='left', va='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # -------------------------------------------------------------------------
    # Panel 2: Eigenvalue convergence
    # -------------------------------------------------------------------------
    ax2 = fig.add_subplot(1, 2, 2)

    N_values = np.array([6, 12, 18, 24, 30, 36, 42, 48])
    exact_eigenvalues = np.array([1, 3, 5, 7])  # First four: 2n+1

    errors = np.zeros((len(N_values), 4))

    for i, N in enumerate(N_values):
        eigs, _, _ = solve_harmonic_oscillator(N, L)
        for j in range(4):
            errors[i, j] = np.abs(eigs[j] - exact_eigenvalues[j])

    # Plot convergence
    markers = ['o', 's', '^', 'D']
    for j in range(4):
        ax2.semilogy(N_values, errors[:, j], f'{markers[j]}-',
                     color=colors[j], linewidth=1.5, markersize=6,
                     label=rf'$\lambda_{j} = {exact_eigenvalues[j]}$')

    # Machine epsilon line
    ax2.axhline(y=2.2e-16, color='gray', linestyle=':', linewidth=1)
    ax2.text(50, 4e-16, 'Machine precision', fontsize=9, color='gray',
             va='bottom', ha='right')

    ax2.set_xlabel(r'Number of grid points $N$')
    ax2.set_ylabel(r'Eigenvalue error $|\lambda_\mathrm{computed} - \lambda_\mathrm{exact}|$')
    ax2.set_title(f'Eigenvalue Convergence ($L={L}$)')
    ax2.legend(loc='upper right', fontsize=9)
    ax2.set_xlim(0, 52)
    ax2.set_ylim(1e-16, 1e2)

    # Clean styling for both panels
    for ax in [ax1, ax2]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#444444')
        ax.spines['bottom'].set_color('#444444')
        ax.tick_params(colors='#444444')

    ax2.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
    ax2.xaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
    ax2.set_axisbelow(True)

    plt.tight_layout()

    # Save figure
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_FILE, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(OUTPUT_FILE.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {OUTPUT_FILE.resolve()}')

    # Print eigenvalue table (replicating Trefethen's Output 8)
    print("\n" + "=" * 70)
    print("Harmonic Oscillator Eigenvalue Convergence")
    print("Exact eigenvalues: lambda_n = 2n + 1")
    print("=" * 70)

    for N in [6, 12, 18, 24, 30, 36]:
        eigs, _, _ = solve_harmonic_oscillator(N, L)
        print(f"\nN = {N}")
        for j in range(4):
            print(f"  lambda_{j} = {eigs[j]:.14f}  (error: {np.abs(eigs[j] - (2*j+1)):.2e})")

    print("\n" + "=" * 70)

    plt.close(fig)


if __name__ == '__main__':
    main()
