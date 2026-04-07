#!/usr/bin/env python3
"""
trap_algebraic_decay.py
Chapter 16: Integration of Periodic Functions

Computational Etude 16.3: Algebraic Convergence and the Regularity
of the Periodic Extension.

Two functions on [0, 2*pi]:

  f_2(x) = |sin(x/2)|        with int = 4
  f_3(x) = |sin(x/2)|^3      with int = 8/3

Although both functions are continuous and periodic, neither is analytic:
their periodic extensions have jumps in higher derivatives at x = 0,
2*pi.  By Euler-Maclaurin, the trapezoidal error is governed by the
endpoint mismatch of odd derivatives:

  - f_2 has a jump in its first derivative at x = 0,
    so the trapezoidal error is O(N^{-2}).
  - f_3 has a jump in its third derivative at x = 0
    (the first and second derivatives match across the period),
    so the trapezoidal error is O(N^{-4}).

Weideman 2002 derives the explicit asymptotics:

  I(f_2) - T_N(f_2) ~ -pi^2/3 * N^{-2}
  I(f_3) - T_N(f_3) ~ -pi^4/30 * N^{-4}

We verify these rates on a log-log plot and report the linear-fit slopes.

Generates Figure 16.3: algebraic_decay.pdf

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
    # The two test functions
    f2 = lambda x: np.abs(np.sin(x / 2.0))
    f3 = lambda x: np.abs(np.sin(x / 2.0))**3

    I2_exact = 4.0
    I3_exact = 8.0 / 3.0

    # Geometric range of N values
    N_values = np.array([8, 16, 32, 64, 128, 256, 512, 1024, 2048])
    err2 = np.zeros(len(N_values))
    err3 = np.zeros(len(N_values))

    for i, N in enumerate(N_values):
        err2[i] = abs(trapezoidal_periodic(f2, N) - I2_exact)
        err3[i] = abs(trapezoidal_periodic(f3, N) - I3_exact)

    # Theoretical asymptotic constants from Weideman 2002, eqs. after (15)
    theory_f2 = (np.pi**2 / 3.0) * N_values.astype(float)**(-2)
    theory_f3 = (np.pi**4 / 30.0) * N_values.astype(float)**(-4)

    # Linear fit of log10 error vs log10 N (omit the first two points to
    # be safely in the asymptotic regime)
    slope_f2 = np.polyfit(np.log10(N_values[2:]), np.log10(err2[2:]), 1)[0]
    slope_f3 = np.polyfit(np.log10(N_values[2:]), np.log10(err3[2:]), 1)[0]

    print("Algebraic convergence test")
    print("-" * 50)
    print(f"f_2 = |sin(x/2)|     I = {I2_exact}")
    print(f"f_3 = |sin(x/2)|^3   I = {I3_exact:.10f}")
    print()
    print(f"{'N':>6s}  {'err(f_2)':>14s}  {'err(f_3)':>14s}")
    for i, N in enumerate(N_values):
        print(f"{N:6d}  {err2[i]:14.4e}  {err3[i]:14.4e}")
    print()
    print(f"Linear-fit slope for f_2: {slope_f2:.3f}  (theory: -2)")
    print(f"Linear-fit slope for f_3: {slope_f3:.3f}  (theory: -4)")

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 5))

    ax.loglog(N_values, err2, 'o-', color=CORAL, markersize=5,
              linewidth=1.2, label=r'$f_2 = |\sin(x/2)|$')
    ax.loglog(N_values, theory_f2, '--', color=CORAL, linewidth=0.8,
              alpha=0.6, label=r'$\pi^2 / (3 N^2)$')
    ax.loglog(N_values, err3, 's-', color=NAVY, markersize=5,
              linewidth=1.2, label=r'$f_3 = |\sin(x/2)|^3$')
    ax.loglog(N_values, theory_f3, '--', color=NAVY, linewidth=0.8,
              alpha=0.6, label=r'$\pi^4 / (30 N^4)$')

    ax.set_xlabel(r'Number of nodes $N$')
    ax.set_ylabel(r'Absolute error $|I_N - I|$')
    ax.set_title(r'Algebraic convergence: regularity of the periodic extension')
    ax.legend(loc='lower left', frameon=True, fancybox=False,
              edgecolor='gray', framealpha=0.9)
    ax.grid(True, which='both', alpha=0.3, linewidth=0.5)
    ax.set_xlim(6, 3000)

    plt.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'algebraic_decay.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'algebraic_decay.png', bbox_inches='tight')
    print(f"\nFigure saved to {OUTPUT_DIR / 'algebraic_decay.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
