#!/usr/bin/env python3
"""
trap_real_line.py
Chapter 16: Integration of Periodic Functions

Computational Etude 16.7: Trapezoidal Rule on the Real Line.

The trapezoidal rule applied with step h to a function w on the entire
real line takes the form

    I_h = h * sum_{k=-infty}^{infty} w(k*h).

Theorem 5.1 of Trefethen-Weideman (2014) says that if w is analytic
and bounded in a strip |Im z| < a and decays to zero at +/- infinity
in that strip, then

    |I_h - I| <= 2*M / (e^{2*pi*a/h} - 1),

i.e. the error decays as O(e^{-2*pi*a/h}) as h -> 0.  This is the same
exponential rate as for periodic integrands, with the natural
substitution N = 2*pi/h.

We illustrate this on the canonical example of Goodwin (1949):

    I = (1/sqrt(pi)) * int_{-infty}^{infty} exp(-x^2) dx = 1.

The integrand is entire, so a may be taken arbitrarily large; in
practice the truncation error of summing only |k| <= n becomes
relevant before the spectral error does.  We sum until the tail terms
fall below 10^{-300}.  This is a faithful reproduction of the
convergence table on page 400 of Trefethen-Weideman 2014.

Generates Figure 16.7: real_line_gaussian.pdf

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


def trapezoidal_real_line(w, h, n_max):
    """Approximate int_{-infty}^{infty} w(x) dx by h * sum_{|k| <= n_max} w(k*h)."""
    k = np.arange(-n_max, n_max + 1)
    return h * np.sum(w(k * h))


def main():
    # The integrand from Goodwin (1949) -- normalised to integrate to 1
    w = lambda x: np.exp(-x**2) / np.sqrt(np.pi)
    I_exact = 1.0

    # Following Trefethen-Weideman p. 400, table after eq. (5.3)
    N_values = np.arange(1, 13)  # so h = 2*pi / N
    h_values = 2.0 * np.pi / N_values

    errors = np.zeros(len(N_values))
    n_truncs = np.zeros(len(N_values), dtype=int)
    for i, h in enumerate(h_values):
        # Decide truncation: pick n_max such that exp(-(n_max*h)^2) < 1e-300
        # i.e. n_max*h > 26.4, so n_max > 26.4/h.  Add a generous buffer.
        n_max = max(int(np.ceil(28.0 / h)), 30)
        n_truncs[i] = n_max
        I_h = trapezoidal_real_line(w, h, n_max)
        errors[i] = abs(I_h - I_exact)
    errors = np.maximum(errors, 1e-17)

    # Theoretical bound: O(e^{-pi^2/h}) for the Gaussian (a -> infty for
    # an entire integrand; the optimal a depends on h).  We compare with
    # exp(-pi^2/h) (taking a = pi/2 as the natural choice that minimises
    # the bound for the entire Gaussian).
    theory = np.exp(-np.pi**2 / h_values)

    # Print the table (mimicking TW2014 p. 400)
    print(f"{'N':>3s}  {'h':>10s}  {'I_h':>22s}  {'|I_h - I|':>14s}")
    for i, N in enumerate(N_values):
        I_h = trapezoidal_real_line(w, h_values[i], n_truncs[i])
        print(f"{N:3d}  {h_values[i]:10.6f}  {I_h:22.16f}  {errors[i]:14.4e}")

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 5))

    ax.semilogy(N_values, errors, 'o-', color=CORAL, markersize=6,
                linewidth=1.2,
                markerfacecolor=CORAL, markeredgecolor=CORAL,
                label=r'$|I_h - I|$, $h = 2\pi/N$')
    ax.semilogy(N_values, theory, '--', color=NAVY, linewidth=1.0,
                label=r'$e^{-\pi^2/h}$')

    ax.set_xlabel(r'$N$ (with step $h = 2\pi/N$)')
    ax.set_ylabel(r'Absolute error')
    ax.set_title(r'Real-line trapezoidal rule on $e^{-x^2}/\sqrt{\pi}$')
    ax.legend(loc='upper right', frameon=True, fancybox=False,
              edgecolor='gray', framealpha=0.9)
    ax.grid(True, which='both', alpha=0.3, linewidth=0.5)
    ax.set_xlim(0, 13)
    ax.set_ylim(1e-18, 1e1)

    plt.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'real_line_gaussian.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'real_line_gaussian.png', bbox_inches='tight')
    print(f"\nFigure saved to {OUTPUT_DIR / 'real_line_gaussian.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
