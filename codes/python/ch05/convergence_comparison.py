#!/usr/bin/env python3
"""
convergence_comparison.py

The "Computational Étude" of Chapter 5: Comparing the accuracy of finite
difference methods (orders 2, 4, 6) against the spectral method for
differentiating a smooth periodic function.

Test Function:
    u(x) = 1 / (2 + sin(x))   on [0, 2π)
    u'(x) = -cos(x) / (2 + sin(x))²

This function is:
    - Smooth (analytic) and periodic
    - Has a non-trivial Fourier spectrum
    - Complex-plane singularities at x = -π/2 ± i·arcsinh(2)

Expected Results:
    - FD2: Error ~ O(N^{-2})  [algebraic]
    - FD4: Error ~ O(N^{-4})  [algebraic]
    - FD6: Error ~ O(N^{-6})  [algebraic]
    - Spectral: Error ~ O(c^{-N})  [geometric/exponential]

The spectral method reaches machine precision around N ≈ 20, while
finite difference methods require much larger N for comparable accuracy.

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

# Import our modules
import sys
sys.path.insert(0, str(Path(__file__).parent))
from fdweights import fdweights
from spectral_matrix_periodic import spectral_diff_periodic

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

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch05' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'convergence_comparison.pdf'


# -----------------------------------------------------------------------------
# Test function and its derivative
# -----------------------------------------------------------------------------
def u_func(x):
    """Test function: u(x) = 1 / (2 + sin(x))"""
    return 1.0 / (2.0 + np.sin(x))


def u_deriv_exact(x):
    """Exact derivative: u'(x) = -cos(x) / (2 + sin(x))^2"""
    return -np.cos(x) / (2.0 + np.sin(x))**2


# -----------------------------------------------------------------------------
# Finite difference matrix construction
# -----------------------------------------------------------------------------
def fd_diff_periodic(N, order):
    """
    Construct a periodic finite difference differentiation matrix.
    """
    h = 2 * np.pi / N

    # Determine stencil size
    stencil_half = order // 2
    stencil_nodes = h * np.arange(-stencil_half, stencil_half + 1)

    # Compute FD weights for first derivative at center
    weights = fdweights(0, stencil_nodes, 1)

    # Build circulant matrix
    D = np.zeros((N, N))
    for i in range(N):
        for k, w in enumerate(weights):
            j = (i + k - stencil_half) % N
            D[i, j] += w

    return D


# -----------------------------------------------------------------------------
# Convergence study
# -----------------------------------------------------------------------------
def compute_errors(N_values, orders):
    """
    Compute maximum differentiation errors for various N and FD orders.
    """
    errors_fd = {order: [] for order in orders}
    errors_spectral = []

    for N in N_values:
        h = 2 * np.pi / N
        x = h * np.arange(N)

        u = u_func(x)
        du_exact = u_deriv_exact(x)

        # Finite difference methods
        for order in orders:
            if order < N:  # Need enough points for the stencil
                D_fd = fd_diff_periodic(N, order)
                du_fd = D_fd @ u
                error = np.max(np.abs(du_fd - du_exact))
            else:
                error = np.nan
            errors_fd[order].append(error)

        # Spectral method
        D_spec, _ = spectral_diff_periodic(N)
        du_spec = D_spec @ u
        error_spec = np.max(np.abs(du_spec - du_exact))
        errors_spectral.append(error_spec)

    return errors_fd, errors_spectral


# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
def main():
    # Grid sizes to test
    N_values = np.array([4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 64])
    orders = [2, 4, 6]

    # Compute errors
    errors_fd, errors_spectral = compute_errors(N_values, orders)

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 5.5))

    # Plot finite difference errors
    colors_fd = {2: CORAL, 4: TEAL, 6: SKY}
    markers_fd = {2: 'o', 4: 's', 6: '^'}

    for order in orders:
        ax.semilogy(N_values, errors_fd[order],
                    marker=markers_fd[order], color=colors_fd[order],
                    linewidth=1.5, markersize=5,
                    label=f'FD{order} (order {order})')

    # Plot spectral errors
    ax.semilogy(N_values, errors_spectral,
                marker='D', color=NAVY, linewidth=2, markersize=6,
                label='Spectral')

    # Add theoretical reference lines
    N_ref = np.linspace(8, 64, 100)

    # FD2: O(N^{-2}) = O(h^2)
    C2 = errors_fd[2][4] * N_values[4]**2  # Fit constant
    ax.semilogy(N_ref, C2 / N_ref**2, '--', color=CORAL, alpha=0.4, linewidth=1)

    # FD4: O(N^{-4})
    C4 = errors_fd[4][4] * N_values[4]**4
    ax.semilogy(N_ref, C4 / N_ref**4, '--', color=TEAL, alpha=0.4, linewidth=1)

    # FD6: O(N^{-6})
    C6 = errors_fd[6][4] * N_values[4]**6
    ax.semilogy(N_ref, C6 / N_ref**6, '--', color=SKY, alpha=0.4, linewidth=1)

    # Machine epsilon line
    ax.axhline(y=2.2e-16, color='gray', linestyle=':', linewidth=1)
    ax.text(66, 4e-16, 'Machine precision', fontsize=9, color='gray',
            va='bottom', ha='right')

    # Add slope annotations - positioned at end of reference lines, above curves
    ax.text(68, C2 / 68**2 * 1.8, r'$O(N^{-2})$', fontsize=10, color=CORAL,
            ha='right', va='bottom')
    ax.text(68, C4 / 68**4 * 2.5, r'$O(N^{-4})$', fontsize=10, color=TEAL,
            ha='right', va='bottom')
    ax.text(68, C6 / 68**6 * 3.0, r'$O(N^{-6})$', fontsize=10, color=SKY,
            ha='right', va='bottom')

    # Add "spectral accuracy" annotation with arrow
    ax.annotate('Spectral\naccuracy!',
                xy=(24, errors_spectral[8]), xytext=(12, 5e-10),
                fontsize=10, color=NAVY, fontweight='bold',
                ha='center',
                arrowprops=dict(arrowstyle='->', color=NAVY, lw=1.5))

    # Labels and title
    ax.set_xlabel(r'Number of grid points $N$', fontsize=11)
    ax.set_ylabel(r'Maximum error $\|u^\prime - Du\|_\infty$', fontsize=11)
    ax.set_title('Differentiation Error: Finite Differences vs Spectral Method',
                 fontsize=12)

    # Add test function info
    textbox = (r'Test function: $u(x) = \frac{1}{2 + \sin(x)}$' + '\n' +
               r'Domain: $[0, 2\pi)$ (periodic)')
    props = dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray')
    ax.text(0.98, 0.98, textbox, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='right', bbox=props)

    # Axis settings
    ax.set_xlim(0, 70)
    ax.set_ylim(1e-16, 1e1)
    ax.legend(loc='upper right', fontsize=9, bbox_to_anchor=(0.75, 0.85))

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

    # Print results table
    print("\nConvergence Study Results:")
    print("=" * 70)
    print(f"{'N':>4}  {'FD2':>12}  {'FD4':>12}  {'FD6':>12}  {'Spectral':>12}")
    print("-" * 70)
    for i, N in enumerate(N_values):
        print(f"{N:4d}  {errors_fd[2][i]:12.4e}  {errors_fd[4][i]:12.4e}  "
              f"{errors_fd[6][i]:12.4e}  {errors_spectral[i]:12.4e}")
    print("=" * 70)

    # Compute convergence rates
    print("\nConvergence Rates (estimated from last two points):")
    for order in orders:
        e1, e2 = errors_fd[order][-2], errors_fd[order][-1]
        n1, n2 = N_values[-2], N_values[-1]
        if e1 > 0 and e2 > 0:
            rate = np.log(e1/e2) / np.log(n2/n1)
            print(f"  FD{order}: {rate:.2f} (expected: {order})")

    # Spectral rate (geometric)
    e1, e2 = errors_spectral[6], errors_spectral[8]  # N=16, N=24
    if e1 > 1e-15 and e2 > 1e-15:
        # log(e) = log(C) - N * log(c)
        rate = (np.log(e1) - np.log(e2)) / (N_values[8] - N_values[6])
        print(f"  Spectral: geometric rate ≈ exp({-rate:.3f}·N)")

    plt.close(fig)


if __name__ == '__main__':
    main()
