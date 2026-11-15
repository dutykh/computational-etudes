#!/usr/bin/env python3
"""
lebesgue_constants_zoom.py

Zoomed comparison of Lebesgue constant growth for Legendre-Gauss-Lobatto and
Chebyshev-Gauss-Lobatto nodes.

This figure omits the exponentially growing equispaced curve, so that the two
logarithmic growths can be told apart:
    - Chebyshev-Gauss-Lobatto: Λ_N = (2/π)(ln N + γ + ln(8/π)) + o(1)
    - Legendre-Gauss-Lobatto:  Λ_N = O(ln N), close to (2/π) ln(N+1) + 0.685,
      the estimate of Hesthaven drawn from numerical experiments
The Legendre-Gauss-Lobatto nodes are the Fekete points of [-1, 1] (Fejér),
and their constant stays about 0.28 below the Chebyshev one.  Neither family
is optimal: every node set has Λ_N > (2/π) ln(N+1) + 0.5212 (Brutman, 1978),
and both curves stay a little above this bound.

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
# Book color scheme
NAVY = '#142D6E'
SKY = '#7896D2'
TEAL = '#16A085'

# Output path
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch04' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'lebesgue_constants_zoom.pdf'


# -----------------------------------------------------------------------------
# Node generation
# -----------------------------------------------------------------------------
def chebyshev_nodes(N):
    """Chebyshev-Gauss-Lobatto nodes on [-1, 1]."""
    j = np.arange(N + 1)
    return np.cos(j * np.pi / N)


def legendre_nodes(N):
    """Legendre-Gauss-Lobatto nodes on [-1, 1]."""
    from numpy.polynomial import legendre as L

    if N == 0:
        return np.array([0.0])
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
    # Create figure (single panel)
    fig, ax = plt.subplots(figsize=(7, 5))

    # Compute Lebesgue constants for range of N
    N_values = np.arange(2, 51)

    Lambda_cheb_values = []
    Lambda_leg_values = []

    for N in N_values:
        Lambda_cheb_values.append(lebesgue_constant(chebyshev_nodes(N)))
        Lambda_leg_values.append(lebesgue_constant(legendre_nodes(N)))

    Lambda_cheb_values = np.array(Lambda_cheb_values)
    Lambda_leg_values = np.array(Lambda_leg_values)

    # Plot computed values
    ax.plot(N_values, Lambda_leg_values, 's-', color=TEAL,
            linewidth=1.5, markersize=5, label='Legendre–Gauss–Lobatto')
    ax.plot(N_values, Lambda_cheb_values, '^-', color=SKY,
            linewidth=1.5, markersize=5, label='Chebyshev–Gauss–Lobatto')

    # Logarithmic asymptotics of the two families
    N_asy = np.linspace(2, 50, 200)

    # Chebyshev-Gauss-Lobatto: (2/π)(ln N + γ + ln(8/π)) + o(1)
    Lambda_cheb_asy = (2/np.pi) * (np.log(N_asy) + np.euler_gamma
                                   + np.log(8/np.pi))
    ax.plot(N_asy, Lambda_cheb_asy, '--', color=SKY, alpha=0.6,
            linewidth=1.5, label=r'$(2/\pi)(\ln N + \gamma + \ln(8/\pi))$')

    # Legendre-Gauss-Lobatto: the estimate (2/π) ln(N+1) + 0.685 of Hesthaven
    Lambda_leg_asy = (2/np.pi) * np.log(N_asy + 1) + 0.685
    ax.plot(N_asy, Lambda_leg_asy, '--', color=TEAL, alpha=0.6,
            linewidth=1.5, label=r'$(2/\pi)\ln(N+1) + 0.685$')

    # Lower bound for every node set, hence for the optimal nodes (Brutman)
    Lambda_low = (2/np.pi) * np.log(N_asy + 1) + 0.5212
    ax.plot(N_asy, Lambda_low, ':', color='#555555', linewidth=1.5,
            label=r'$(2/\pi)\ln(N+1) + 0.5212$, lower bound for all nodes')

    ax.set_xlabel(r'$N$ (polynomial degree)')
    ax.set_ylabel(r'Lebesgue constant $\Lambda_N$')
    ax.set_title('Lebesgue Constants: Legendre vs Chebyshev', fontsize=11)
    ax.set_xlim(0, 52)
    ax.set_ylim(1.0, 4.4)
    ax.legend(loc='upper left', fontsize=9.5)

    # Clean styling
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

    # Print comparison table
    print("\nLebesgue Constants Comparison (Legendre vs Chebyshev):")
    print("-" * 66)
    print(f"{'N':>4} {'Legendre':>12} {'Chebyshev':>12} {'Difference':>12}"
          f" {'Lower bound':>13}")
    print("-" * 66)
    for i, N in enumerate(N_values[::8]):  # Print every 8th value
        idx = i * 8
        low = (2/np.pi) * np.log(N + 1) + 0.5212
        print(f"{N:>4d} {Lambda_leg_values[idx]:>12.4f} "
              f"{Lambda_cheb_values[idx]:>12.4f} "
              f"{Lambda_cheb_values[idx] - Lambda_leg_values[idx]:>12.4f}"
              f" {low:>13.4f}")
    print("-" * 66)

    plt.close(fig)


if __name__ == '__main__':
    main()
