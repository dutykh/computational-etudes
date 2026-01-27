#!/usr/bin/env python3
"""
lebesgue_functions.py

Visualizes and compares Lebesgue functions for different node distributions.

The Lebesgue function is:
    Λ_N(x) = Σ_{k=0}^{N} |L_k(x)|

where L_k(x) are the Lagrange basis polynomials. The Lebesgue constant is:
    Λ_N = max_{x ∈ [-1,1]} Λ_N(x)

The Lebesgue constant bounds the interpolation error:
    ‖f - p_N‖_∞ ≤ (1 + Λ_N) E_N(f)

where E_N(f) is the best polynomial approximation error.

Asymptotic growth rates:
    - Chebyshev: Λ_N = (2/π) ln(N) + O(1) ≈ O(ln N)
    - Legendre:  Λ_N ≈ O(√N)
    - Equispaced: Λ_N ≈ 2^N / (N ln N) (exponential!)

The dramatically different growth rates explain why Chebyshev interpolation
converges while equispaced interpolation can diverge.

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
N_FINE = 1000  # Number of points for smooth plotting

# Book color scheme
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'

# Output path
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch04' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'lebesgue_functions.pdf'


# -----------------------------------------------------------------------------
# Node generation
# -----------------------------------------------------------------------------
def equispaced_nodes(N):
    """Equispaced nodes on [-1, 1]."""
    return np.linspace(-1, 1, N + 1)


def chebyshev_nodes(N):
    """Chebyshev-Gauss-Lobatto nodes on [-1, 1]."""
    j = np.arange(N + 1)
    return np.cos(j * np.pi / N)


def legendre_nodes(N):
    """Legendre-Gauss-Lobatto nodes on [-1, 1]."""
    # Use numpy's polynomial module to get Legendre roots
    # LGL nodes are roots of (1-x²)P'_N(x), plus endpoints ±1
    if N == 0:
        return np.array([0.0])

    # Get interior nodes as roots of P'_N
    from numpy.polynomial import legendre as L
    # P'_N has N-1 roots in (-1, 1)
    # These are the Legendre-Gauss nodes without endpoints

    # For LGL nodes, we need zeros of (1-x²)P'_N(x)
    # which are ±1 plus the N-1 zeros of P'_N(x)
    if N == 1:
        return np.array([-1.0, 1.0])

    # Derivative of P_N
    coeffs = np.zeros(N + 1)
    coeffs[N] = 1
    P_N = L.Legendre(coeffs)
    P_N_deriv = P_N.deriv()

    # Find roots of P'_N
    interior = P_N_deriv.roots()
    interior = np.real(interior)  # Remove tiny imaginary parts

    # Combine with endpoints
    return np.sort(np.concatenate([[-1], interior, [1]]))


# -----------------------------------------------------------------------------
# Lebesgue function computation
# -----------------------------------------------------------------------------
def lagrange_basis(x_nodes, k, x_eval):
    """Compute the k-th Lagrange basis polynomial at points x_eval."""
    n = len(x_nodes)
    x_eval = np.atleast_1d(x_eval)
    L_k = np.ones_like(x_eval, dtype=float)

    for j in range(n):
        if j != k:
            L_k *= (x_eval - x_nodes[j]) / (x_nodes[k] - x_nodes[j])

    return L_k


def lebesgue_function(x_nodes, x_eval):
    """
    Compute the Lebesgue function: Λ_N(x) = Σ_{k=0}^{N} |L_k(x)|
    """
    N = len(x_nodes) - 1
    Lambda = np.zeros_like(x_eval, dtype=float)

    for k in range(N + 1):
        Lambda += np.abs(lagrange_basis(x_nodes, k, x_eval))

    return Lambda


def lebesgue_constant(x_nodes, n_eval=2000):
    """Compute the Lebesgue constant: Λ_N = max_{x ∈ [-1,1]} Λ_N(x)"""
    x_eval = np.linspace(-1, 1, n_eval)
    Lambda = lebesgue_function(x_nodes, x_eval)
    return np.max(Lambda)


# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
def main():
    # Fine grid for plotting
    x_fine = np.linspace(-1, 1, N_FINE)

    # Create figure with two panels
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # ==========================================================================
    # Left panel: Lebesgue functions for N = 10
    # ==========================================================================
    N = 10

    # Compute Lebesgue functions
    x_equi = equispaced_nodes(N)
    x_cheb = chebyshev_nodes(N)
    x_leg = legendre_nodes(N)

    Lambda_equi = lebesgue_function(x_equi, x_fine)
    Lambda_cheb = lebesgue_function(x_cheb, x_fine)
    Lambda_leg = lebesgue_function(x_leg, x_fine)

    ax1.plot(x_fine, Lambda_equi, color=CORAL, linewidth=1.5,
             label=f'Equispaced ($\\Lambda_{{{N}}}$ = {np.max(Lambda_equi):.1f})')
    ax1.plot(x_fine, Lambda_leg, color=TEAL, linewidth=1.5,
             label=f'Legendre ($\\Lambda_{{{N}}}$ = {np.max(Lambda_leg):.2f})')
    ax1.plot(x_fine, Lambda_cheb, color=SKY, linewidth=1.5,
             label=f'Chebyshev ($\\Lambda_{{{N}}}$ = {np.max(Lambda_cheb):.2f})')

    ax1.axhline(y=1, color='#888888', linewidth=0.5, linestyle='--')

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$\Lambda_N(x)$')
    ax1.set_title(f'Lebesgue Functions ($N = {N}$)', fontsize=11)
    ax1.set_xlim(-1, 1)
    ax1.set_ylim(0, 35)
    ax1.legend(loc='upper center', fontsize=9)

    # Add annotation
    ax1.annotate('Peaks at boundaries\nfor equispaced nodes',
                 xy=(-0.95, 28), xytext=(-0.6, 25),
                 fontsize=9, color='gray',
                 arrowprops=dict(arrowstyle='->', color='gray', lw=0.8))

    # ==========================================================================
    # Right panel: Lebesgue constant growth with N
    # ==========================================================================
    N_values = np.arange(2, 31)

    Lambda_equi_values = []
    Lambda_cheb_values = []
    Lambda_leg_values = []

    for N in N_values:
        Lambda_equi_values.append(lebesgue_constant(equispaced_nodes(N)))
        Lambda_cheb_values.append(lebesgue_constant(chebyshev_nodes(N)))
        Lambda_leg_values.append(lebesgue_constant(legendre_nodes(N)))

    Lambda_equi_values = np.array(Lambda_equi_values)
    Lambda_cheb_values = np.array(Lambda_cheb_values)
    Lambda_leg_values = np.array(Lambda_leg_values)

    ax2.semilogy(N_values, Lambda_equi_values, 'o-', color=CORAL,
                 linewidth=1.5, markersize=4, label='Equispaced')
    ax2.semilogy(N_values, Lambda_leg_values, 's-', color=TEAL,
                 linewidth=1.5, markersize=4, label='Legendre')
    ax2.semilogy(N_values, Lambda_cheb_values, '^-', color=SKY,
                 linewidth=1.5, markersize=4, label='Chebyshev')

    # Plot asymptotic curves
    N_asy = np.linspace(5, 30, 100)
    # Chebyshev: (2/π) ln(N) + 0.6
    Lambda_cheb_asy = (2/np.pi) * np.log(N_asy) + 0.6
    ax2.semilogy(N_asy, Lambda_cheb_asy, '--', color=SKY, alpha=0.5,
                 linewidth=1, label=r'$(2/\pi)\ln N$')

    ax2.set_xlabel(r'$N$ (polynomial degree)')
    ax2.set_ylabel(r'Lebesgue constant $\Lambda_N$')
    ax2.set_title('Growth of Lebesgue Constants', fontsize=11)
    ax2.set_xlim(0, 32)
    ax2.legend(loc='upper left', fontsize=9)

    # Add text annotations for growth rates
    ax2.text(25, 1e6, r'$O(2^N)$', fontsize=10, color=CORAL)
    ax2.text(25, 4, r'$O(\sqrt{N})$', fontsize=10, color=TEAL)
    ax2.text(25, 2.5, r'$O(\ln N)$', fontsize=10, color=SKY)

    # Clean styling for both panels
    for ax in [ax1, ax2]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#444444')
        ax.spines['bottom'].set_color('#444444')
        ax.tick_params(colors='#444444')
        ax.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
        ax.xaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
        ax.set_axisbelow(True)

    plt.tight_layout()

    # Save figure
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_FILE, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(OUTPUT_FILE.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {OUTPUT_FILE.resolve()}')

    # Print summary table
    print("\nLebesgue Constants Comparison:")
    print("-" * 65)
    print(f"{'N':>4} {'Equispaced':>15} {'Legendre':>15} {'Chebyshev':>15}")
    print("-" * 65)
    for i, N in enumerate(N_values[::5]):  # Print every 5th value
        idx = i * 5
        print(f"{N:>4d} {Lambda_equi_values[idx]:>15.2f} "
              f"{Lambda_leg_values[idx]:>15.4f} {Lambda_cheb_values[idx]:>15.4f}")
    print("-" * 65)

    # Theoretical bounds
    print("\nAsymptotic Growth Rates:")
    print("  Chebyshev:  Λ_N ~ (2/π) ln(N) + O(1)")
    print("  Legendre:   Λ_N ~ O(√N)")
    print("  Equispaced: Λ_N ~ 2^(N+1) / (eN ln N)")

    plt.close(fig)


if __name__ == '__main__':
    main()
