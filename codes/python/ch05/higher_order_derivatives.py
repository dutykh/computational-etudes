#!/usr/bin/env python3
"""
higher_order_derivatives.py

Demonstrates spectral differentiation for higher-order derivatives (up to 4th order)
and compares D² construction methods: matrix squaring vs direct formula.

The test function is u(x) = exp(-sin(2x)) on [0, 2π), which is smooth, periodic,
and has higher frequency content than simpler test functions.

This script generates:
1. A 2×2 figure showing derivatives from 1st to 4th order
2. A comparison figure for D² construction methods
3. A convergence table showing accuracy vs derivative order

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
from pathlib import Path

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
# Configuration
# -----------------------------------------------------------------------------
N_DEMO = 32  # Grid points for demonstration figure
N_FINE = 500  # Fine grid for plotting exact solutions

# Book color scheme
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'
PURPLE = '#8E44AD'
ORANGE = '#E67E22'

# Output paths
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch05' / 'python'


# -----------------------------------------------------------------------------
# Test function and its exact derivatives
# -----------------------------------------------------------------------------
def test_function(x):
    """u(x) = exp(-sin(2x))"""
    return np.exp(-np.sin(2 * x))


def exact_derivative_1(x):
    """u'(x) = -2*cos(2x)*exp(-sin(2x))"""
    return -2 * np.cos(2 * x) * np.exp(-np.sin(2 * x))


def exact_derivative_2(x):
    """u''(x) = 4*(sin(2x) + cos²(2x))*exp(-sin(2x))"""
    s = np.sin(2 * x)
    c = np.cos(2 * x)
    return 4 * (s + c**2) * np.exp(-s)


def exact_derivative_3(x):
    """u'''(x) = 8*cos(2x)*(sin²(2x) - 3*sin(2x))*exp(-sin(2x))"""
    s = np.sin(2 * x)
    c = np.cos(2 * x)
    return 8 * c * (s**2 - 3 * s) * np.exp(-s)


def exact_derivative_4(x):
    """u''''(x) computed from the pattern."""
    s = np.sin(2 * x)
    c = np.cos(2 * x)
    # u'''' = 16*exp(-s)*[-s³ + 3s² + 5sc² - 3c² - s²c²]
    # Simplified: factor out and combine terms
    term = -s**3 + 3*s**2 + 5*s*c**2 - 3*c**2 - s**2 * c**2
    return 16 * term * np.exp(-s)


# -----------------------------------------------------------------------------
# Spectral differentiation matrix (periodic)
# -----------------------------------------------------------------------------
def spectral_diff_matrix(N):
    """
    Construct the periodic spectral differentiation matrix.

    Parameters:
        N : int - Number of grid points (should be even)

    Returns:
        D : ndarray (N, N) - Differentiation matrix
        x : ndarray (N,) - Grid points on [0, 2π)
    """
    h = 2 * np.pi / N
    x = h * np.arange(N)
    D = np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            if i != j:
                D[i, j] = 0.5 * ((-1) ** (i - j)) / np.tan((i - j) * np.pi / N)

    return D, x


# -----------------------------------------------------------------------------
# Main visualization functions
# -----------------------------------------------------------------------------
def create_derivatives_figure():
    """Create 2×2 figure showing derivatives from 1st to 4th order."""

    # Construct differentiation matrix
    D, x = spectral_diff_matrix(N_DEMO)

    # Sample test function at grid points
    u = test_function(x)

    # Compute numerical derivatives via matrix powers
    u1_num = D @ u
    u2_num = D @ D @ u
    u3_num = D @ D @ D @ u
    u4_num = D @ D @ D @ D @ u

    # Fine grid for exact solutions
    x_fine = np.linspace(0, 2 * np.pi, N_FINE)
    u1_exact = exact_derivative_1(x_fine)
    u2_exact = exact_derivative_2(x_fine)
    u3_exact = exact_derivative_3(x_fine)
    u4_exact = exact_derivative_4(x_fine)

    # Exact values at grid points for error computation
    u1_exact_grid = exact_derivative_1(x)
    u2_exact_grid = exact_derivative_2(x)
    u3_exact_grid = exact_derivative_3(x)
    u4_exact_grid = exact_derivative_4(x)

    # Compute errors
    err1 = np.max(np.abs(u1_num - u1_exact_grid))
    err2 = np.max(np.abs(u2_num - u2_exact_grid))
    err3 = np.max(np.abs(u3_num - u3_exact_grid))
    err4 = np.max(np.abs(u4_num - u4_exact_grid))

    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    derivatives = [
        (u1_exact, u1_num, err1, r"$u'(x)$", "First derivative"),
        (u2_exact, u2_num, err2, r"$u''(x)$", "Second derivative"),
        (u3_exact, u3_num, err3, r"$u'''(x)$", "Third derivative"),
        (u4_exact, u4_num, err4, r"$u''''(x)$", "Fourth derivative"),
    ]

    colors_exact = [NAVY, NAVY, NAVY, NAVY]
    colors_num = [CORAL, TEAL, PURPLE, ORANGE]

    for ax, (exact, num, err, ylabel, title), col_num in zip(
        axes.flat, derivatives, colors_num
    ):
        # Plot exact solution
        ax.plot(x_fine, exact, '-', color=NAVY, linewidth=1.5, label='Exact')
        # Plot numerical solution
        ax.plot(x, num, 'o', color=col_num, markersize=5, label='Spectral')

        ax.set_xlabel(r'$x$')
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=11)

        # Set x-axis ticks
        ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
        ax.set_xticklabels([r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
        ax.set_xlim(0, 2 * np.pi)

        ax.legend(loc='best', fontsize=9)

        # Add error annotation
        ax.text(0.95, 0.05, f'Max error: {err:.2e}',
                transform=ax.transAxes, fontsize=9,
                verticalalignment='bottom', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        # Clean styling
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#444444')
        ax.spines['bottom'].set_color('#444444')
        ax.tick_params(colors='#444444')
        ax.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
        ax.xaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
        ax.set_axisbelow(True)

    fig.suptitle(rf'Higher-Order Spectral Derivatives of $u(x) = e^{{-\sin(2x)}}$ ($N = {N_DEMO}$)',
                 fontsize=12, y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # Save figure
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_file = OUTPUT_DIR / 'higher_order_derivatives.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig)

    return err1, err2, err3, err4


def create_d2_comparison_figure():
    """Demonstrate D² via matrix squaring: structure and accuracy."""

    N = 16  # Grid size for visualization

    # Construct D and D²
    D, x = spectral_diff_matrix(N)
    D2 = D @ D

    # Test second derivative accuracy
    u = test_function(x)
    u2_num = D2 @ u
    u2_exact = exact_derivative_2(x)
    error = np.max(np.abs(u2_num - u2_exact))

    # Fine grid for plotting
    x_fine = np.linspace(0, 2 * np.pi, N_FINE)

    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    # Panel 1: Structure of D² matrix
    vmax = np.max(np.abs(D2))
    im1 = axes[0].imshow(D2, cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    axes[0].set_title(r'$D^2 = D \cdot D$ (matrix structure)', fontsize=11)
    axes[0].set_xlabel('Column $j$')
    axes[0].set_ylabel('Row $i$')
    plt.colorbar(im1, ax=axes[0], shrink=0.8)

    # Panel 2: Eigenvalue spectrum of D²
    eigvals = np.linalg.eigvals(D2)
    eigvals_sorted = np.sort(eigvals.real)
    # Theoretical eigenvalues for d²/dx²: -k² for k = -N/2+1, ..., N/2
    k_vals = np.arange(-N//2 + 1, N//2 + 1)
    theoretical_eigvals = -k_vals**2

    axes[1].plot(eigvals_sorted, 'o', color=NAVY, markersize=6, label='Numerical')
    axes[1].plot(np.sort(theoretical_eigvals), 's', color=CORAL, markersize=4,
                 markerfacecolor='none', label='Theoretical $-k^2$')
    axes[1].set_xlabel('Index')
    axes[1].set_ylabel('Eigenvalue')
    axes[1].set_title(r'Eigenvalues of $D^2$', fontsize=11)
    axes[1].legend(loc='best', fontsize=9)
    axes[1].grid(True, alpha=0.3)
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)

    # Panel 3: Second derivative accuracy
    axes[2].plot(x_fine, exact_derivative_2(x_fine), '-', color=NAVY,
                 linewidth=1.5, label='Exact')
    axes[2].plot(x, u2_num, 'o', color=TEAL, markersize=6, label='$D^2 u$')
    axes[2].set_xlabel(r'$x$')
    axes[2].set_ylabel(r"$u''(x)$")
    axes[2].set_title(rf"Second derivative (error = {error:.2e})", fontsize=11)
    axes[2].set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
    axes[2].set_xticklabels([r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
    axes[2].set_xlim(0, 2 * np.pi)
    axes[2].legend(loc='best', fontsize=9)
    axes[2].grid(True, alpha=0.3)
    axes[2].spines['top'].set_visible(False)
    axes[2].spines['right'].set_visible(False)

    fig.suptitle(rf'Matrix Squaring for Second Derivatives ($N = {N}$)',
                 fontsize=12, y=1.02)
    plt.tight_layout()

    # Save figure
    output_file = OUTPUT_DIR / 'd2_comparison.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig)

    return error


def generate_convergence_table():
    """Generate convergence table for different derivative orders."""

    N_values = [8, 16, 32, 64]
    results = []

    print("\n" + "=" * 70)
    print("Convergence Table: Higher-Order Derivatives")
    print(f"Test function: u(x) = exp(-sin(2x))")
    print("=" * 70)
    print(f"{'N':>6} {'Error u\'':>14} {'Error u\'\'':>14} "
          f"{'Error u\'\'\'':>14} {'Error u\'\'\'\'':>14}")
    print("-" * 70)

    for N in N_values:
        D, x = spectral_diff_matrix(N)
        u = test_function(x)

        # Numerical derivatives
        u1_num = D @ u
        u2_num = D @ D @ u
        u3_num = D @ D @ D @ u
        u4_num = D @ D @ D @ D @ u

        # Exact derivatives at grid points
        u1_exact = exact_derivative_1(x)
        u2_exact = exact_derivative_2(x)
        u3_exact = exact_derivative_3(x)
        u4_exact = exact_derivative_4(x)

        # Errors
        err1 = np.max(np.abs(u1_num - u1_exact))
        err2 = np.max(np.abs(u2_num - u2_exact))
        err3 = np.max(np.abs(u3_num - u3_exact))
        err4 = np.max(np.abs(u4_num - u4_exact))

        results.append((N, err1, err2, err3, err4))
        print(f"{N:>6d} {err1:>14.2e} {err2:>14.2e} {err3:>14.2e} {err4:>14.2e}")

    print("=" * 70)

    return results


def main():
    """Main function to generate all figures and tables."""

    print("=" * 70)
    print("Higher-Order Derivatives Demonstration")
    print("=" * 70)

    # Create main figure
    print("\n1. Creating higher-order derivatives figure...")
    err1, err2, err3, err4 = create_derivatives_figure()
    print(f"   Errors: u'={err1:.2e}, u''={err2:.2e}, "
          f"u'''={err3:.2e}, u''''={err4:.2e}")

    # Create D² comparison figure
    print("\n2. Creating D² matrix squaring demonstration...")
    d2_error = create_d2_comparison_figure()
    print(f"   Second derivative error via D²: {d2_error:.2e}")

    # Generate convergence table
    print("\n3. Generating convergence table...")
    results = generate_convergence_table()

    print("\n" + "=" * 70)
    print("All figures generated successfully!")
    print("=" * 70)


if __name__ == '__main__':
    main()
