#!/usr/bin/env python3
"""
trap_poisson_kernel.py
Chapter 16: Integration of Periodic Functions

Computational Etude 16.4: Geometric Convergence on the Poisson Kernel.

For the function

    f_4(x) = 1 / (a - cos x),  a > 1

the exact integral is

    int_0^{2*pi} f_4(x) dx = 2*pi / sqrt(a^2 - 1).

Weideman (2002) derives the explicit error formula for the trapezoidal
rule by expanding f_4 as a Fourier series via partial fractions:

    I(f_4) - T_N(f_4) = -8*pi * r/(1 - r^2) * r^N / (1 - r^N),

where r = a - sqrt(a^2 - 1) is the smaller root of r^2 - 2*a*r + 1 = 0
and satisfies 0 < r < 1.  Therefore the error decays geometrically at
rate r^N.

The geometric rate has a complex-plane interpretation: f_4 has simple
poles in the complex theta-plane where cos(theta) = a, i.e. at

    theta_k = 2*k*pi +/- i*log(a + sqrt(a^2 - 1)) = 2*k*pi +/- i*alpha,

so the strip of analyticity has half-width alpha = log(a + sqrt(a^2-1))
= -log(r), and the strip-analyticity theorem gives O(e^{-alpha N}) =
O(r^N), in perfect agreement with Weideman's elementary derivation.

We take a = 2, so r = 2 - sqrt(3) ~ 0.2679 and the predicted decay
gives roughly log10(2 - sqrt(3)) ~ -0.572 digits per node.

Generates Figure 16.4: poisson_kernel.pdf  (two-panel figure)

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

NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch16' / 'python'


def trapezoidal_periodic(f, N):
    """N-point periodic trapezoidal rule on [0, 2*pi)."""
    theta = 2.0 * np.pi * np.arange(N) / N
    return (2.0 * np.pi / N) * np.sum(f(theta))


def main():
    a = 2.0
    f4 = lambda x: 1.0 / (a - np.cos(x))
    I_exact = 2.0 * np.pi / np.sqrt(a**2 - 1.0)

    # Geometric rate from the Weideman analysis
    r = a - np.sqrt(a**2 - 1.0)
    alpha = np.log(a + np.sqrt(a**2 - 1.0))  # half-width of analyticity strip

    print("Poisson kernel f_4 = 1 / (a - cos x)")
    print(f"a = {a}, exact integral I = 2*pi/sqrt(a^2 - 1) = {I_exact:.10f}")
    print(f"Geometric rate r = a - sqrt(a^2 - 1) = {r:.10f}")
    print(f"Strip half-width alpha = log(a + sqrt(a^2 - 1)) = {alpha:.10f}")
    print(f"log10(r) = {np.log10(r):.4f}  (expected slope of semilog plot)")
    print()

    N_values = np.arange(2, 51, 2)
    errors = np.zeros(len(N_values))
    for i, N in enumerate(N_values):
        errors[i] = abs(trapezoidal_periodic(f4, N) - I_exact)
    errors = np.maximum(errors, 1e-17)

    # Theoretical curve from Weideman (2002), eq. (21):
    #     err = 8*pi * r/(1 - r^2) * r^N / (1 - r^N)
    theory = 8.0 * np.pi * r / (1.0 - r**2) * r**N_values / (1.0 - r**N_values)
    theory = np.maximum(theory, 1e-17)

    print(f"{'N':>4s}  {'error':>14s}  {'theory':>14s}")
    for i, N in enumerate(N_values[:15]):
        print(f"{N:4d}  {errors[i]:14.4e}  {theory[i]:14.4e}")

    # ------------------------------------------------------------------
    # Two-panel plot: convergence + complex-plane pole picture
    # ------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # Panel 1: convergence
    ax1.semilogy(N_values, errors, 'o-', color=CORAL, markersize=5,
                 linewidth=1.2, label=r'Trapezoidal error')
    ax1.semilogy(N_values, theory, '--', color=NAVY, linewidth=1.0,
                 label=r'$\dfrac{8\pi r}{1-r^2}\dfrac{r^N}{1-r^N}$')
    ax1.set_xlabel(r'Number of nodes $N$')
    ax1.set_ylabel(r'Absolute error $|I_N - I|$')
    ax1.set_title(r'Geometric decay at rate $r^N$ for $a = 2$')
    ax1.legend(loc='upper right', frameon=True, fancybox=False,
               edgecolor='gray', framealpha=0.9)
    ax1.grid(True, which='both', alpha=0.3, linewidth=0.5)
    ax1.set_ylim(1e-18, 1e1)

    # Panel 2: complex-plane pole picture
    # Real axis: [0, 2*pi]
    ax2.axhline(y=0, color='black', linewidth=1.0, zorder=2)
    # Strip of analyticity: |Im theta| < alpha
    ax2.axhspan(-alpha, alpha, color=NAVY, alpha=0.15,
                label=r'strip of analyticity')
    ax2.axhline(y=alpha, color=NAVY, linewidth=0.8, linestyle='--')
    ax2.axhline(y=-alpha, color=NAVY, linewidth=0.8, linestyle='--')
    # Sample points (the trapezoidal nodes for N = 12)
    N_demo = 12
    theta_nodes = 2.0 * np.pi * np.arange(N_demo) / N_demo
    ax2.plot(theta_nodes, np.zeros_like(theta_nodes), 'o',
             color=CORAL, markersize=6, markerfacecolor=CORAL,
             label=f'$N = {N_demo}$ trapezoidal nodes', zorder=3)
    # Poles at theta = 0 +/- i*alpha and 2*pi +/- i*alpha
    pole_theta = [0.0, 2.0 * np.pi]
    for tp in pole_theta:
        ax2.plot([tp, tp], [alpha, -alpha], '', color='none')
        ax2.plot(tp, alpha, 'x', color=NAVY, markersize=10,
                 markeredgewidth=2.0, zorder=4)
        ax2.plot(tp, -alpha, 'x', color=NAVY, markersize=10,
                 markeredgewidth=2.0, zorder=4)
    ax2.text(np.pi, alpha + 0.15,
             r'pole at $\pm i\,\log(a + \sqrt{a^2 - 1})$',
             ha='center', fontsize=9, color=NAVY)
    ax2.set_xlim(-0.5, 2.0 * np.pi + 0.5)
    ax2.set_ylim(-2.0, 2.0)
    ax2.set_xlabel(r'$\mathrm{Re}\,\theta$')
    ax2.set_ylabel(r'$\mathrm{Im}\,\theta$')
    ax2.set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi])
    ax2.set_xticklabels(['$0$', r'$\pi/2$', r'$\pi$',
                         r'$3\pi/2$', r'$2\pi$'])
    ax2.set_title(r'Poles in the complex $\theta$-plane')
    ax2.legend(loc='lower right', frameon=True, fancybox=False,
               edgecolor='gray', framealpha=0.9, fontsize=9)
    ax2.grid(True, alpha=0.3, linewidth=0.5)

    plt.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'poisson_kernel.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'poisson_kernel.png', bbox_inches='tight')
    print(f"\nFigure saved to {OUTPUT_DIR / 'poisson_kernel.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
