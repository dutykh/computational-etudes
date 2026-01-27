#!/usr/bin/env python3
"""
collocation_vs_galerkin.py

Compares the Collocation (pseudospectral) and Galerkin methods for solving
a reaction-diffusion boundary value problem. This example illustrates two
fundamental approaches to the Method of Weighted Residuals.

Problem:
    u''(x) - 4u(x) = -1,   -1 ≤ x ≤ 1
    u(-1) = 0,  u(1) = 0

Exact solution:
    u_exact(x) = (1/4)(1 - cosh(2x)/cosh(2))

Basis functions (satisfy homogeneous BCs):
    φ₀(x) = (1 - x²)
    φ₁(x) = (1 - x²)x²

Trial function:
    u₁(x) = a₀φ₀(x) + a₁φ₁(x)

Methods compared:
    1. Collocation: Force residual to zero at x = 0 and x = 0.5
    2. Galerkin: Force residual orthogonal to basis functions

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Études: A Spectral Approach"
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
N_POINTS = 500

# Book color scheme
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
GREEN = '#27AE60'

# Output path
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch03' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'collocation_vs_galerkin.pdf'


# -----------------------------------------------------------------------------
# Exact solution
# -----------------------------------------------------------------------------
def u_exact(x):
    """
    Exact solution: u(x) = (1/4)(1 - cosh(2x)/cosh(2))

    Verification:
        u'(x) = -(1/4)(2sinh(2x)/cosh(2)) = -(1/2)sinh(2x)/cosh(2)
        u''(x) = -cosh(2x)/cosh(2)

    Substituting into ODE:
        u'' - 4u = -cosh(2x)/cosh(2) - 4·(1/4)(1 - cosh(2x)/cosh(2))
                 = -cosh(2x)/cosh(2) - 1 + cosh(2x)/cosh(2)
                 = -1 ✓

    Boundary conditions:
        u(±1) = (1/4)(1 - cosh(2)/cosh(2)) = (1/4)(1 - 1) = 0 ✓
    """
    return 0.25 * (1 - np.cosh(2*x) / np.cosh(2))


# -----------------------------------------------------------------------------
# Basis functions
# -----------------------------------------------------------------------------
def phi0(x):
    """φ₀(x) = 1 - x² (vanishes at x = ±1)"""
    return 1 - x**2


def phi1(x):
    """φ₁(x) = x² - x⁴ = (1 - x²)x² (vanishes at x = ±1)"""
    return x**2 - x**4


def phi0_xx(x):
    """Second derivative: φ₀''(x) = -2"""
    return -2 * np.ones_like(x) if hasattr(x, '__len__') else -2


def phi1_xx(x):
    """Second derivative: φ₁''(x) = 2 - 12x²"""
    return 2 - 12*x**2


# -----------------------------------------------------------------------------
# Collocation method
# -----------------------------------------------------------------------------
def solve_collocation():
    """
    Solve using collocation at x = 0 and x = 0.5.

    The residual is:
        R(x) = u₁''(x) - 4u₁(x) + 1
             = a₀(φ₀'' - 4φ₀) + a₁(φ₁'' - 4φ₁) + 1

    At x = 0:
        L[φ₀](0) = -2 - 4(1) = -6
        L[φ₁](0) = 2 - 4(0) = 2
        Equation: -6a₀ + 2a₁ + 1 = 0

    At x = 0.5:
        L[φ₀](0.5) = -2 - 4(0.75) = -5
        L[φ₁](0.5) = 2 - 12(0.25) - 4(0.1875) = 2 - 3 - 0.75 = -1.75
        Equation: -5a₀ - 1.75a₁ + 1 = 0

    System:
        6a₀ - 2a₁ = 1
        5a₀ + 1.75a₁ = 1
    """
    # Operator L = d²/dx² - 4
    def L_phi0(x):
        return phi0_xx(x) - 4*phi0(x)

    def L_phi1(x):
        return phi1_xx(x) - 4*phi1(x)

    # Collocation points
    x1, x2 = 0.0, 0.5

    # Build system matrix
    A = np.array([
        [L_phi0(x1), L_phi1(x1)],
        [L_phi0(x2), L_phi1(x2)]
    ])

    # Right-hand side: f = -1, so R = Lu - f = Lu + 1 = 0 means Lu = -1
    b = np.array([-1.0, -1.0])

    # Solve
    coeffs = np.linalg.solve(A, b)

    return coeffs[0], coeffs[1]


# -----------------------------------------------------------------------------
# Galerkin method
# -----------------------------------------------------------------------------
def solve_galerkin():
    """
    Solve using Galerkin method: ⟨R, φₖ⟩ = 0 for k = 0, 1.

    This requires computing integrals:
        A_{ij} = ∫₋₁¹ (φᵢ'' - 4φᵢ) φⱼ dx
        bᵢ = ∫₋₁¹ (-1) · φᵢ dx

    The matrix is symmetric for self-adjoint operators.
    """
    from scipy import integrate

    # Operator L = d²/dx² - 4
    def L_phi0(x):
        return phi0_xx(x) - 4*phi0(x)

    def L_phi1(x):
        return phi1_xx(x) - 4*phi1(x)

    # Compute matrix entries A_{ij} = ⟨L[φⱼ], φᵢ⟩
    A00, _ = integrate.quad(lambda x: L_phi0(x) * phi0(x), -1, 1)
    A01, _ = integrate.quad(lambda x: L_phi1(x) * phi0(x), -1, 1)
    A10, _ = integrate.quad(lambda x: L_phi0(x) * phi1(x), -1, 1)
    A11, _ = integrate.quad(lambda x: L_phi1(x) * phi1(x), -1, 1)

    A = np.array([[A00, A01], [A10, A11]])

    # Compute RHS: bᵢ = ⟨-1, φᵢ⟩ = -∫φᵢ dx
    b0, _ = integrate.quad(lambda x: -1 * phi0(x), -1, 1)
    b1, _ = integrate.quad(lambda x: -1 * phi1(x), -1, 1)

    b = np.array([b0, b1])

    # Solve
    coeffs = np.linalg.solve(A, b)

    return coeffs[0], coeffs[1]


# -----------------------------------------------------------------------------
# Approximate solutions
# -----------------------------------------------------------------------------
def u_collocation(x, a0, a1):
    """Collocation approximation."""
    return a0 * phi0(x) + a1 * phi1(x)


def u_galerkin(x, a0, a1):
    """Galerkin approximation."""
    return a0 * phi0(x) + a1 * phi1(x)


# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
def main():
    # Solve both systems
    a0_coll, a1_coll = solve_collocation()
    a0_gal, a1_gal = solve_galerkin()

    print("Collocation coefficients:")
    print(f"  a₀ = {a0_coll:.6f}")
    print(f"  a₁ = {a1_coll:.6f}")
    print("\nGalerkin coefficients:")
    print(f"  a₀ = {a0_gal:.6f}")
    print(f"  a₁ = {a1_gal:.6f}")

    # Spatial grid
    x = np.linspace(-1, 1, N_POINTS)

    # Compute solutions
    u_ex = u_exact(x)
    u_coll = u_collocation(x, a0_coll, a1_coll)
    u_gal = u_galerkin(x, a0_gal, a1_gal)

    error_coll = u_ex - u_coll
    error_gal = u_ex - u_gal

    # Collocation points for visualization
    x_coll_pts = np.array([0.0, 0.5])
    u_coll_pts = u_collocation(x_coll_pts, a0_coll, a1_coll)

    # Create two-panel figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.5))

    # Left panel: All solutions
    ax1.plot(x, u_ex, color=NAVY, linewidth=1.8, label='Exact')
    ax1.plot(x, u_coll, '--', color=CORAL, linewidth=1.5,
             label=f'Collocation')
    ax1.plot(x, u_gal, '-.', color=GREEN, linewidth=1.5,
             label=f'Galerkin')
    ax1.plot(x_coll_pts, u_coll_pts, 's', color=CORAL, markersize=8,
             markerfacecolor='none', markeredgewidth=1.5,
             label='Collocation pts')

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$u(x)$')
    ax1.set_xlim(-1, 1)
    ax1.set_title('Solutions Comparison', fontsize=11)
    ax1.legend(loc='lower center', frameon=True, fancybox=False,
               edgecolor='none', facecolor='white', framealpha=0.9,
               ncol=2)

    # Right panel: Errors
    ax2.plot(x, error_coll, color=CORAL, linewidth=1.5, label='Collocation')
    ax2.plot(x, error_gal, color=GREEN, linewidth=1.5, label='Galerkin')
    ax2.axhline(y=0, color='gray', linewidth=0.5, linestyle='--')

    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r'$u_{\mathrm{exact}} - u_{\mathrm{approx}}$')
    ax2.set_xlim(-1, 1)
    ax2.set_title('Error Comparison', fontsize=11)
    ax2.legend(loc='upper right', frameon=True, fancybox=False,
               edgecolor='none', facecolor='white', framealpha=0.9)

    # Clean styling for both panels
    for ax in [ax1, ax2]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#444444')
        ax.spines['bottom'].set_color('#444444')
        ax.tick_params(colors='#444444')
        ax.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
        ax.set_axisbelow(True)

    plt.tight_layout()

    # Save figure
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_FILE, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(OUTPUT_FILE.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'\nFigure saved to: {OUTPUT_FILE.resolve()}')

    # Print comparison table
    print("\n" + "=" * 65)
    print("Comparison at u(0) - the central maximum")
    print("=" * 65)
    print(f"{'Method':<20} {'u(0)':<15} {'Abs. Error':<15}")
    print("-" * 65)
    print(f"{'Exact':<20} {u_exact(0):<15.6f} {0:<15.6f}")
    print(f"{'Collocation (N=2)':<20} {u_collocation(0, a0_coll, a1_coll):<15.6f} "
          f"{abs(u_exact(0) - u_collocation(0, a0_coll, a1_coll)):<15.6f}")
    print(f"{'Galerkin (N=2)':<20} {u_galerkin(0, a0_gal, a1_gal):<15.6f} "
          f"{abs(u_exact(0) - u_galerkin(0, a0_gal, a1_gal)):<15.6f}")
    print("=" * 65)

    # RMS errors
    rms_coll = np.sqrt(np.mean(error_coll**2))
    rms_gal = np.sqrt(np.mean(error_gal**2))
    print(f"\nRMS Error (Collocation): {rms_coll:.6f}")
    print(f"RMS Error (Galerkin):    {rms_gal:.6f}")

    # Max errors
    max_coll = np.max(np.abs(error_coll))
    max_gal = np.max(np.abs(error_gal))
    print(f"\nMax Error (Collocation): {max_coll:.6f}")
    print(f"Max Error (Galerkin):    {max_gal:.6f}")

    plt.close(fig)


if __name__ == '__main__':
    main()
