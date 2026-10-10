#!/usr/bin/env python3
"""
trap_fft_coefficients.py
Chapter 16: Integration of Periodic Functions

Computational Etude 16.8: Computing Fourier Coefficients with the FFT.

The Fourier coefficients of a 2*pi-periodic function f are defined by

    c_k = (1/(2*pi)) * int_0^{2*pi} f(theta) * exp(-i*k*theta) d theta.

Each c_k is itself a periodic integral, so the natural numerical
estimator is the periodic trapezoidal rule:

    c_k_hat = (1/N) * sum_{j=0}^{N-1} f(theta_j) * exp(-2*pi*i*k*j/N).

This is exactly the FFT of the sample vector (f(theta_0), ..., f(theta_{N-1})),
divided by N.  All the convergence theorems we have proved for the
periodic trapezoidal rule therefore transfer directly to the FFT
approximation of Fourier coefficients: if f is analytic in a strip of
half-width a, then c_k_hat - c_k = O(e^{-aN}) for each fixed k, and
the FFT computes the entire vector of N coefficients in O(N log N)
operations.

This is the bridge from the periodic trapezoidal rule to all of Fourier
spectral methods (chapters 9 to 11): every spectral coefficient ever
computed by the FFT is, secretly, an instance of the trapezoidal rule
applied to a Fourier integral.

We test this on f(x) = 1/(2 - cos x), whose Fourier coefficients are
known in closed form:

    c_k = (1/sqrt(3)) * (2 - sqrt(3))^|k|.

Generates Figure 16.8: fft_coefficients.pdf  (two-panel figure)

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


def fourier_coeffs_trap(f, N):
    """Compute (c_0, c_1, ..., c_{N-1}) via the FFT of N samples of f."""
    theta = 2.0 * np.pi * np.arange(N) / N
    return np.fft.fft(f(theta)) / N


def main():
    a = 2.0
    f = lambda x: 1.0 / (a - np.cos(x))

    # Closed-form Fourier coefficients: c_k = (1/sqrt(3)) * (2 - sqrt(3))^|k|
    r = a - np.sqrt(a**2 - 1.0)
    print(f"a = {a}, r = a - sqrt(a^2 - 1) = {r:.10f}")
    print(f"Exact c_0 = 1/sqrt(a^2 - 1) = {1.0 / np.sqrt(a**2 - 1.0):.10f}\n")

    def c_exact(k):
        return (1.0 / np.sqrt(a**2 - 1.0)) * r**abs(k)

    # FFT-based coefficients with N = 32
    N = 32
    c_fft = fourier_coeffs_trap(f, N)

    # The FFT puts c_k at index k for k = 0, 1, ..., N/2,
    # and c_{-(N-k)} at index k for k = N/2 + 1, ..., N - 1.
    # Convert to the symmetric ordering k = -N/2, ..., N/2 - 1.
    c_fft_sym = np.fft.fftshift(c_fft)
    k_sym = np.arange(-N // 2, N // 2)

    # Exact reference, taking the real part only (symmetry => imag = 0)
    c_ref = np.array([c_exact(k) for k in k_sym])
    err = np.abs(c_fft_sym - c_ref)

    print(f"FFT vs exact Fourier coefficients (N = {N})")
    print(f"{'k':>5s}  {'|c_k| FFT':>14s}  {'|c_k| exact':>14s}  {'error':>14s}")
    # Print only k = 0..N/2 - 1 (the resolved band)
    for k in range(N // 2):
        idx = N // 2 + k
        print(f"{k:5d}  {np.abs(c_fft_sym[idx]):14.6e}  "
              f"{c_ref[idx]:14.6e}  {err[idx]:14.4e}")

    err = np.maximum(err, 1e-17)

    # ------------------------------------------------------------------
    # Two-panel plot
    # ------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # Panel 1: magnitudes |c_k| -- exact vs FFT
    ax1.semilogy(k_sym, np.abs(c_ref), '-', color=NAVY, linewidth=1.2,
                 label=r'exact: $\dfrac{1}{\sqrt{3}}(2 - \sqrt{3})^{|k|}$')
    ax1.semilogy(k_sym, np.abs(c_fft_sym), 'o', color=CORAL,
                 markersize=5, markerfacecolor=CORAL, markeredgecolor=CORAL,
                 label=f'FFT, $N = {N}$')
    ax1.set_xlabel(r'Mode index $k$')
    ax1.set_ylabel(r'$|\hat f_k|$')
    ax1.set_title(r'(a) Fourier-coefficient magnitudes for $1/(2 - \cos x)$')
    ax1.legend(loc='upper right', frameon=True, fancybox=False,
               edgecolor='gray', framealpha=0.9)
    ax1.grid(True, which='both', alpha=0.3, linewidth=0.5)
    ax1.set_ylim(1e-12, 1e1)

    # Panel 2: error |c_k_hat - c_k|
    ax2.semilogy(k_sym, err, 'o-', color=CORAL, markersize=5,
                 linewidth=1.0, markerfacecolor=CORAL,
                 markeredgecolor=CORAL,
                 label=r'$|\hat f_k^{\mathrm{FFT}} - \hat f_k|$')
    ax2.set_xlabel(r'Mode index $k$')
    ax2.set_ylabel(r'Absolute error in $\hat f_k$')
    ax2.set_title(r'(b) FFT error: machine precision in the resolved band')
    ax2.legend(loc='upper right', frameon=True, fancybox=False,
               edgecolor='gray', framealpha=0.9)
    ax2.grid(True, which='both', alpha=0.3, linewidth=0.5)
    ax2.set_ylim(1e-18, 1e0)

    plt.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'fft_coefficients.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'fft_coefficients.png', bbox_inches='tight')
    print(f"\nFigure saved to {OUTPUT_DIR / 'fft_coefficients.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
