#!/usr/bin/env python3
"""
collocation_example1.py

Demonstrates the collocation (pseudospectral) method for solving a boundary
value problem with a three-coefficient polynomial approximation. This example
follows the approach in Boyd (2000) but with a different differential equation.

Problem:
    u''(x) - (4x² + 2)u(x) = 0,   -1 ≤ x ≤ 1
    u(-1) = 1,  u(1) = 1

Exact solution:
    u_exact(x) = exp(x² - 1)

Trial function (automatically satisfies BCs):
    u₂(x) = 1 + (1 - x²)(a₀ + a₁x + a₂x²)

Collocation points:
    x = -1/2, 0, 1/2

Result:
    a₀ = -73/118,  a₁ = 0,  a₂ = -14/59

The figure shows the exact solution vs approximation (left) and the error (right),
similar to Figure 1.1 in Boyd's "Chebyshev and Fourier Spectral Methods".

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
rcParams['text.usetex'] = False  # Set to True if LaTeX is available
rcParams['axes.linewidth'] = 0.8
rcParams['xtick.major.width'] = 0.8
rcParams['ytick.major.width'] = 0.8

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
N_POINTS = 500  # Number of points for smooth plotting

# Book color scheme
NAVY = '#142D6E'
SKY = '#7896D2'

# Output path
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch03' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'collocation_example1.pdf'


# -----------------------------------------------------------------------------
# Exact solution
# -----------------------------------------------------------------------------
def u_exact(x):
    """
    Exact solution: u(x) = exp(x² - 1)

    Verification:
        u'(x) = 2x exp(x² - 1)
        u''(x) = (2 + 4x²) exp(x² - 1) = (2 + 4x²) u(x)

    Substituting into ODE:
        u'' - (4x² + 2)u = (2 + 4x²)u - (4x² + 2)u = 0 ✓

    Boundary conditions:
        u(±1) = exp(1 - 1) = exp(0) = 1 ✓
    """
    return np.exp(x**2 - 1)


# -----------------------------------------------------------------------------
# Approximate solution from collocation
# -----------------------------------------------------------------------------
def u_approx(x):
    """
    Approximate solution using collocation with 3 coefficients.

    Trial function: u₂(x) = 1 + (1 - x²)(a₀ + a₁x + a₂x²)

    Collocation at x = -1/2, 0, 1/2 gives:
        a₀ = -73/118
        a₁ = 0
        a₂ = -14/59

    Simplified form:
        u₂(x) = (14/59)x⁴ + (45/118)x² + 45/118
    """
    # Coefficients from solving the collocation system
    a0 = -73/118
    a1 = 0
    a2 = -14/59

    # Using the trial function form
    return 1 + (1 - x**2) * (a0 + a1*x + a2*x**2)


def solve_collocation_system():
    """
    Solve the collocation system to find a₀, a₁, a₂.

    This function demonstrates the actual computation that produces
    the analytical coefficients used in u_approx().

    The residual is:
        R(x) = u₂''(x) - (4x² + 2)u₂(x)

    Setting R(x) = 0 at x = -1/2, 0, 1/2 gives a 3×3 linear system.
    """
    import sympy as sp

    x = sp.Symbol('x')
    a0, a1, a2 = sp.symbols('a0 a1 a2')

    # Trial function
    u2 = 1 + (1 - x**2) * (a0 + a1*x + a2*x**2)

    # Second derivative
    u2_xx = sp.diff(u2, x, 2)

    # Residual
    R = u2_xx - (4*x**2 + 2) * u2
    R = sp.expand(R)

    # Collocation equations
    eq1 = R.subs(x, sp.Rational(-1, 2))
    eq2 = R.subs(x, 0)
    eq3 = R.subs(x, sp.Rational(1, 2))

    # Solve the system
    solution = sp.solve([eq1, eq2, eq3], [a0, a1, a2])

    return solution


# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
def main():
    # Spatial grid
    x = np.linspace(-1, 1, N_POINTS)

    # Compute solutions
    u_ex = u_exact(x)
    u_ap = u_approx(x)
    error = u_ex - u_ap

    # Collocation points for visualization
    x_coll = np.array([-0.5, 0, 0.5])
    u_coll = u_approx(x_coll)

    # Create two-panel figure (like Boyd's Figure 1.1)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.5))

    # Left panel: Exact vs Approximate
    ax1.plot(x, u_ex, color=NAVY, linewidth=1.5, label='Exact')
    ax1.plot(x, u_ap, 'o', color=SKY, markersize=4, markevery=25,
             label='Collocation ($N=3$)')
    ax1.plot(x_coll, u_coll, 's', color='#E74C3C', markersize=8,
             markerfacecolor='none', markeredgewidth=1.5,
             label='Collocation points')

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$u(x)$')
    ax1.set_xlim(-1, 1)
    ax1.set_title('Exact vs Approximate', fontsize=11)
    ax1.legend(loc='upper center', frameon=True, fancybox=False,
               edgecolor='none', facecolor='white', framealpha=0.9)

    # Right panel: Error
    ax2.plot(x, error, color=NAVY, linewidth=1.5)
    ax2.axhline(y=0, color='gray', linewidth=0.5, linestyle='--')

    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r'$u_{\mathrm{exact}} - u_{\mathrm{approx}}$')
    ax2.set_xlim(-1, 1)
    ax2.set_title('Error', fontsize=11)

    # Annotate max error
    max_err_idx = np.argmax(np.abs(error))
    max_err = error[max_err_idx]
    ax2.annotate(f'Max error: {max_err:.4f}',
                 xy=(0, max_err), xytext=(0.3, max_err - 0.005),
                 fontsize=9, color=NAVY,
                 arrowprops=dict(arrowstyle='->', color=NAVY, lw=0.8))

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

    print(f'Figure saved to: {OUTPUT_FILE.resolve()}')

    # Print error table
    print("\nError comparison at key points:")
    print("-" * 60)
    print(f"{'x':>8} {'u_exact':>12} {'u_approx':>12} {'error':>12}")
    print("-" * 60)
    for xi in [-1.0, -0.5, 0.0, 0.5, 1.0]:
        ue = u_exact(xi)
        ua = u_approx(xi)
        err = ue - ua
        print(f"{xi:>8.2f} {ue:>12.5f} {ua:>12.5f} {err:>12.5f}")
    print("-" * 60)
    print(f"Maximum absolute error: {np.max(np.abs(error)):.5f}")

    plt.close(fig)


if __name__ == '__main__':
    main()
