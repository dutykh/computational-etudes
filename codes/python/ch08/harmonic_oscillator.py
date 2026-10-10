#!/usr/bin/env python3
"""
harmonic_oscillator.py

Solving the quantum harmonic oscillator eigenvalue problem:
    -u'' + x^2 u = lambda u,  x in R

using Chebyshev spectral methods on a truncated domain [-L, L],
with a comparison to the periodic (Fourier) spectral method.

Exact eigenvalues: lambda_n = 2n + 1 for n = 0, 1, 2, ...
Exact eigenfunctions: u_n(x) = H_n(x) * exp(-x^2/2) (Hermite functions)

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""

import math
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pathlib import Path
from scipy.special import hermite
from scipy.linalg import toeplitz

# Import Chebyshev matrix from Chapter 7
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'ch07'))
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
# Chebyshev method: solve harmonic oscillator
# -----------------------------------------------------------------------------
def solve_chebyshev(N, L):
    """
    Solve -u'' + x^2 u = lambda u on [-L, L] using Chebyshev spectral method.

    Parameters
    ----------
    N : int
        Number of Chebyshev intervals (N+1 points, N-1 interior)
    L : float
        Half-width of domain

    Returns
    -------
    eigenvalues : ndarray
        Computed eigenvalues (sorted)
    eigenvectors : ndarray
        Corresponding eigenvectors at interior points (columns)
    x_full : ndarray
        All grid points (including boundaries)
    """
    # Chebyshev points on [-1, 1]
    D, t = cheb_matrix(N)

    # Map to physical domain [-L, L]
    x_full = L * t

    # Second derivative matrix on physical domain
    D2 = (D @ D) / L**2

    # Extract interior system (strip boundary rows and columns)
    D2_int = D2[1:N, 1:N]
    x_int = x_full[1:N]

    # Build operator: A = -D2 + diag(x^2)
    A = -D2_int + np.diag(x_int**2)

    # Solve eigenvalue problem (A is not symmetric due to Chebyshev D^2)
    eigenvalues, eigenvectors = np.linalg.eig(A)
    idx = np.argsort(np.real(eigenvalues))
    eigenvalues = np.real(eigenvalues[idx])
    eigenvectors = np.real(eigenvectors[:, idx])

    return eigenvalues, eigenvectors, x_full


# -----------------------------------------------------------------------------
# Periodic (Fourier) method: solve harmonic oscillator
# -----------------------------------------------------------------------------
def periodic_d2_matrix(N):
    """Construct the periodic second derivative matrix for N equispaced points."""
    h = 2 * np.pi / N
    column = np.zeros(N)
    column[0] = -np.pi**2 / (3 * h**2) - 1.0 / 6.0
    for k in range(1, N):
        column[k] = -0.5 * ((-1)**k) / np.sin(h * k / 2.0)**2
    return toeplitz(column)


def solve_fourier(N, L):
    """
    Solve -u'' + x^2 u = lambda u on [-L, L] using periodic spectral method.

    Parameters
    ----------
    N : int
        Number of equispaced grid points
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
    h = 2 * np.pi / N
    x_periodic = h * np.arange(N)
    x = L * (x_periodic - np.pi) / np.pi

    D2 = periodic_d2_matrix(N) * (np.pi / L)**2
    V = np.diag(x**2)
    A = -D2 + V

    eigenvalues, eigenvectors = np.linalg.eigh(A)
    return eigenvalues, eigenvectors, x


# -----------------------------------------------------------------------------
# Exact Hermite functions (normalized)
# -----------------------------------------------------------------------------
def hermite_function(n, x):
    """Compute the n-th normalized Hermite function."""
    norm = 1.0 / np.sqrt(2**n * math.factorial(n) * np.sqrt(np.pi))
    H_n = hermite(n)
    return norm * H_n(x) * np.exp(-x**2 / 2.0)


# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
def main():
    L = 8.0

    fig = plt.figure(figsize=(12, 5))

    # -------------------------------------------------------------------------
    # Panel 1: Chebyshev eigenfunctions
    # -------------------------------------------------------------------------
    ax1 = fig.add_subplot(1, 2, 1)

    N = 64
    eigenvalues, eigenvectors, x_full = solve_chebyshev(N, L)
    x_int = x_full[1:N]

    colors = [NAVY, CORAL, TEAL, PURPLE, ORANGE]

    # Fine uniform grid for smooth exact curves
    x_fine = np.linspace(-6, 6, 500)

    for n in range(5):
        # Numerical eigenfunction (interior only)
        u_int = eigenvectors[:, n]

        # Embed in full grid with boundary zeros
        u_num = np.zeros(N + 1)
        u_num[1:N] = u_int

        # Exact eigenfunction on Chebyshev grid (for sign/scale matching)
        u_exact = hermite_function(n, x_full)

        # Exact eigenfunction on fine grid (for plotting)
        u_exact_fine = hermite_function(n, x_fine)

        # Fix sign
        if np.dot(u_num, u_exact) < 0:
            u_num = -u_num

        # Scale to match
        scale = np.max(np.abs(u_exact)) / np.max(np.abs(u_num))
        u_num = u_num * scale

        # Plot exact on fine grid (line) and numerical on Chebyshev grid (markers)
        ax1.plot(x_fine, u_exact_fine + 2 * n, '-', color=colors[n], linewidth=1.5)
        ax1.plot(x_full[::4], (u_num + 2 * n)[::4], 'o', color=colors[n],
                 markersize=4, alpha=0.7)

    # Generic legend entries (compact, 2 items only)
    ax1.plot([], [], '-', color='gray', linewidth=1.5, label='Exact')
    ax1.plot([], [], 'o', color='gray', markersize=4, label='Chebyshev')

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$u_n(x)$ (shifted)')
    ax1.set_title(f'Chebyshev Eigenfunctions ($N={N}$, $L={L:.0f}$)')
    ax1.set_xlim(-6, 6)
    ax1.legend(loc='upper right', fontsize=9)

    # -------------------------------------------------------------------------
    # Panel 2: Convergence comparison (Chebyshev vs Fourier)
    # -------------------------------------------------------------------------
    ax2 = fig.add_subplot(1, 2, 2)

    N_values = np.array([12, 18, 24, 30, 36, 42, 48])
    exact_eigenvalues = np.array([1, 3, 5, 7])

    errors_cheb = np.zeros((len(N_values), 4))
    errors_four = np.zeros((len(N_values), 4))

    for i, Nval in enumerate(N_values):
        # Chebyshev
        eigs_c, _, _ = solve_chebyshev(Nval, L)
        for j in range(4):
            errors_cheb[i, j] = max(np.abs(eigs_c[j] - exact_eigenvalues[j]), 1e-16)

        # Fourier
        eigs_f, _, _ = solve_fourier(Nval, L)
        for j in range(4):
            errors_four[i, j] = max(np.abs(eigs_f[j] - exact_eigenvalues[j]), 1e-16)

    # Plot lambda_0 and lambda_1 convergence for both methods
    ax2.semilogy(N_values, errors_four[:, 0], 'o-', color=NAVY, linewidth=1.5,
                 markersize=6, label=r'Fourier $\lambda_0$')
    ax2.semilogy(N_values, errors_cheb[:, 0], 's--', color=CORAL, linewidth=1.5,
                 markersize=6, label=r'Chebyshev $\lambda_0$')
    ax2.semilogy(N_values, errors_four[:, 1], '^-', color=TEAL, linewidth=1.5,
                 markersize=6, label=r'Fourier $\lambda_1$')
    ax2.semilogy(N_values, errors_cheb[:, 1], 'v--', color=PURPLE, linewidth=1.5,
                 markersize=6, label=r'Chebyshev $\lambda_1$')

    # Machine epsilon line
    ax2.axhline(y=2.2e-16, color='gray', linestyle=':', linewidth=1)
    ax2.text(50, 4e-16, 'Machine precision', fontsize=9, color='gray',
             va='bottom', ha='right')

    ax2.set_xlabel(r'Number of grid points $N$')
    ax2.set_ylabel(r'Eigenvalue error')
    ax2.set_title(f'Fourier vs Chebyshev ($L={L:.0f}$)')
    ax2.legend(loc='upper right', fontsize=9)
    ax2.set_xlim(8, 52)
    ax2.set_ylim(1e-16, 1e2)

    # Clean styling
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

    # -------------------------------------------------------------------------
    # Print comparison table
    # -------------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("Fourier vs Chebyshev Eigenvalue Comparison (L = 8)")
    print("Exact eigenvalues: lambda_n = 2n + 1")
    print("=" * 78)
    print(f"{'N':>4}  {'Fourier err(l0)':>16}  {'Cheb err(l0)':>16}  "
          f"{'Fourier err(l1)':>16}  {'Cheb err(l1)':>16}")
    print("-" * 78)

    for i, Nval in enumerate(N_values):
        print(f"{Nval:>4d}  {errors_four[i, 0]:>16.2e}  {errors_cheb[i, 0]:>16.2e}  "
              f"{errors_four[i, 1]:>16.2e}  {errors_cheb[i, 1]:>16.2e}")

    print("=" * 78)

    # Chebyshev eigenvalue table (N=36)
    print("\nChebyshev eigenvalues at N = 36, L = 8:")
    eigs_c, _, _ = solve_chebyshev(36, L)
    for j in range(4):
        print(f"  lambda_{j} = {eigs_c[j]:.14f}  (error: {np.abs(eigs_c[j] - (2*j+1)):.2e})")

    print("\nFourier eigenvalues at N = 36, L = 8:")
    eigs_f, _, _ = solve_fourier(36, L)
    for j in range(4):
        print(f"  lambda_{j} = {eigs_f[j]:.14f}  (error: {np.abs(eigs_f[j] - (2*j+1)):.2e})")

    plt.close(fig)


if __name__ == '__main__':
    main()
