#!/usr/bin/env python3
"""
spectral_derivatives_demo.py

Demonstrates the periodic spectral differentiation matrix by computing
first and second derivatives of a smooth periodic function.

Test function: u(x) = exp(sin^2(x))
- First derivative:  u'(x)  = sin(2x) * exp(sin^2(x))
- Second derivative: u''(x) = [sin^2(2x) + 2*cos(2x)] * exp(sin^2(x))

With just N = 64 grid points, the spectral method achieves near machine precision!

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

# Import our spectral differentiation module
import sys
sys.path.insert(0, str(Path(__file__).parent))
from spectral_matrix_periodic import spectral_diff_periodic

# -----------------------------------------------------------------------------
# Publication-quality matplotlib configuration
# -----------------------------------------------------------------------------
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman', 'CMU Serif', 'DejaVu Serif']
rcParams['mathtext.fontset'] = 'cm'
rcParams['axes.labelsize'] = 12
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
# Book color scheme
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'
PURPLE = '#8E44AD'

# Output path
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch05' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'spectral_derivatives_demo.pdf'


# -----------------------------------------------------------------------------
# Main computation and visualization
# -----------------------------------------------------------------------------
def main():
    # Number of grid points
    N = 64

    # Construct the spectral differentiation matrix
    D, x = spectral_diff_periodic(N)

    # Fine grid for plotting exact solutions
    x_fine = np.linspace(0, 2 * np.pi, 200)

    # Test function: u(x) = exp(sin^2(x))
    u = np.exp(np.sin(x)**2)
    u_fine = np.exp(np.sin(x_fine)**2)

    # Exact first derivative: u'(x) = sin(2x) * exp(sin^2(x))
    u1_exact = np.sin(2 * x) * np.exp(np.sin(x)**2)
    u1_fine = np.sin(2 * x_fine) * np.exp(np.sin(x_fine)**2)

    # Exact second derivative: u''(x) = [sin^2(2x) + 2*cos(2x)] * exp(sin^2(x))
    u2_exact = (np.sin(2 * x)**2 + 2 * np.cos(2 * x)) * np.exp(np.sin(x)**2)
    u2_fine = (np.sin(2 * x_fine)**2 + 2 * np.cos(2 * x_fine)) * np.exp(np.sin(x_fine)**2)

    # Numerical derivatives using spectral differentiation matrix
    u1_num = D @ u           # First derivative: D * u
    u2_num = D @ (D @ u)     # Second derivative: D^2 * u

    # Compute errors
    error_u1 = np.max(np.abs(u1_num - u1_exact))
    error_u2 = np.max(np.abs(u2_num - u2_exact))

    # -------------------------------------------------------------------------
    # Create figure
    # -------------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Left panel: u(x) and u'(x)
    ax1 = axes[0]
    ax1.plot(x_fine, u_fine, '-', color=NAVY, linewidth=2, label=r"$u(x) = e^{\sin^2 x}$")
    ax1.plot(x_fine, u1_fine, '-', color=CORAL, linewidth=2, label=r"$u'(x)$ exact")
    ax1.plot(x, u1_num, 'o', color=TEAL, markersize=7, markerfacecolor='white',
             markeredgewidth=2, label=r"$u'(x)$ spectral", zorder=5)

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$u(x)$, $u\'(x)$')
    ax1.set_title(rf'First Derivative ($N = {N}$)')
    ax1.set_xlim([0, 2 * np.pi])
    ax1.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
    ax1.set_xticklabels([r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)

    # Add error annotation
    ax1.text(0.05, 0.05, f'Max error: {error_u1:.2e}',
             transform=ax1.transAxes, fontsize=10, color=TEAL,
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor=TEAL))

    # Right panel: u''(x)
    ax2 = axes[1]
    ax2.plot(x_fine, u2_fine, '-', color=PURPLE, linewidth=2, label=r"$u''(x)$ exact")
    ax2.plot(x, u2_num, 's', color=TEAL, markersize=7, markerfacecolor='white',
             markeredgewidth=2, label=r"$u''(x)$ spectral", zorder=5)

    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r"$u''(x)$")
    ax2.set_title(rf'Second Derivative ($N = {N}$)')
    ax2.set_xlim([0, 2 * np.pi])
    ax2.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
    ax2.set_xticklabels([r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)

    # Add error annotation
    ax2.text(0.05, 0.05, f'Max error: {error_u2:.2e}',
             transform=ax2.transAxes, fontsize=10, color=TEAL,
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor=TEAL))

    plt.tight_layout()

    # -------------------------------------------------------------------------
    # Save figure
    # -------------------------------------------------------------------------
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_FILE, bbox_inches='tight', pad_inches=0.1)
    fig.savefig(OUTPUT_FILE.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.1, dpi=300)

    print(f'Figure saved to: {OUTPUT_FILE.resolve()}')

    # -------------------------------------------------------------------------
    # Print results
    # -------------------------------------------------------------------------
    print('\n' + '=' * 60)
    print('Spectral Differentiation Demo')
    print('=' * 60)
    print(f'Test function: u(x) = exp(sin²(x))')
    print(f'Number of grid points: N = {N}')
    print('-' * 60)
    print(f'First derivative  max error: {error_u1:.4e}')
    print(f'Second derivative max error: {error_u2:.4e}')
    print('-' * 60)
    print('With only 64 points, we achieve near machine precision!')
    print('=' * 60)

    # -------------------------------------------------------------------------
    # Convergence study for table
    # -------------------------------------------------------------------------
    print('\nConvergence Study (for table):')
    print('-' * 50)
    print(f'{"N":>6}  {"Error u\'":>14}  {"Error u\'\'":>14}')
    print('-' * 50)

    for N_test in [8, 16, 32, 64]:
        D_test, x_test = spectral_diff_periodic(N_test)
        u_test = np.exp(np.sin(x_test)**2)
        u1_exact_test = np.sin(2 * x_test) * np.exp(np.sin(x_test)**2)
        u2_exact_test = (np.sin(2 * x_test)**2 + 2 * np.cos(2 * x_test)) * np.exp(np.sin(x_test)**2)
        u1_num_test = D_test @ u_test
        u2_num_test = D_test @ (D_test @ u_test)
        err1 = np.max(np.abs(u1_num_test - u1_exact_test))
        err2 = np.max(np.abs(u2_num_test - u2_exact_test))
        print(f'{N_test:6d}  {err1:14.4e}  {err2:14.4e}')

    print('-' * 50)

    plt.close(fig)


if __name__ == '__main__':
    main()
