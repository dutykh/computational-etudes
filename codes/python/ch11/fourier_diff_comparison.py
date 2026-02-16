#!/usr/bin/env python3
"""
fourier_diff_comparison.py
Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs

Compares Fourier differentiation via dense matrix multiplication (O(N^2))
with FFT-based differentiation (O(N log N)).

Test function: u(x) = exp(sin(x)) on [0, 2*pi)
Exact derivative: u'(x) = cos(x) * exp(sin(x))

Generates two figures:
  - fourier_diff_accuracy.pdf: spectral convergence of both methods
  - fourier_diff_timing.pdf: timing comparison (matrix vs FFT)

Author: Dr. Denys Dutykh
Date: February 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pathlib import Path
import time

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

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch11' / 'python'


def fourier_diff_matrix(N):
    """
    Build the N x N Fourier differentiation matrix for the periodic grid
    x_j = 2*pi*j/N, j = 0, ..., N-1.

    The (i, j) entry is the derivative of the j-th Lagrange basis function
    evaluated at x_i.
    """
    D = np.zeros((N, N))
    h = 2.0 * np.pi / N
    for i in range(N):
        for j in range(N):
            if i != j:
                D[i, j] = 0.5 * (-1.0)**(i + j) / np.tan((i - j) * h / 2.0)
    return D


def fourier_diff_fft(u):
    """
    Compute the derivative of u on a periodic grid using FFT.
    u : array of length N (values at x_j = 2*pi*j/N)
    """
    N = len(u)
    u_hat = np.fft.fft(u)
    # Wavenumbers in FFT order
    if N % 2 == 0:
        k = np.concatenate([np.arange(N // 2), [0], np.arange(1 - N // 2, 0)])
    else:
        k = np.concatenate([np.arange((N - 1) // 2 + 1),
                           np.arange(-(N - 1) // 2, 0)])
    du_hat = 1j * k * u_hat
    return np.fft.ifft(du_hat).real


def main():
    # Test function and exact derivative
    u_func = lambda x: np.exp(np.sin(x))
    du_exact_func = lambda x: np.cos(x) * np.exp(np.sin(x))

    # -------------------------------------------------------------------------
    # Figure 1: Accuracy comparison
    # -------------------------------------------------------------------------
    N_values = np.arange(4, 52, 2)
    err_matrix = np.zeros(len(N_values))
    err_fft = np.zeros(len(N_values))

    for idx, N in enumerate(N_values):
        x = 2.0 * np.pi * np.arange(N) / N
        u = u_func(x)
        du_exact = du_exact_func(x)

        # Matrix method
        D = fourier_diff_matrix(N)
        du_mat = D @ u
        err_matrix[idx] = np.max(np.abs(du_mat - du_exact))

        # FFT method
        du_fft = fourier_diff_fft(u)
        err_fft[idx] = np.max(np.abs(du_fft - du_exact))

    fig1, ax1 = plt.subplots(figsize=(6, 4.5))
    ax1.semilogy(N_values, err_matrix, 'o-', color=NAVY, markersize=5,
                 linewidth=1.5, label='Matrix $Du$')
    ax1.semilogy(N_values, err_fft, 's--', color=CORAL, markersize=5,
                 linewidth=1.5, label=r'FFT: $\mathcal{F}^{-1}(ik\,\hat{u})$')
    ax1.set_xlabel(r'$N$ (number of grid points)')
    ax1.set_ylabel(r'$\|u_x^{\mathrm{num}} - u_x^{\mathrm{exact}}\|_\infty$')
    ax1.set_title(r'Spectral Convergence: $u(x) = e^{\sin x}$')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(bottom=1e-16)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    fig1.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig1.savefig(OUTPUT_DIR / 'fourier_diff_accuracy.pdf',
                 bbox_inches='tight', pad_inches=0.05)
    print(f'Figure saved: {(OUTPUT_DIR / "fourier_diff_accuracy.pdf").resolve()}')
    plt.close(fig1)

    # -------------------------------------------------------------------------
    # Figure 2: Timing comparison
    # -------------------------------------------------------------------------
    N_timing = [2**p for p in range(3, 14)]  # 8 to 8192
    time_matrix = []
    time_fft = []
    n_repeats = 5

    for N in N_timing:
        x = 2.0 * np.pi * np.arange(N) / N
        u = u_func(x)

        # Time matrix method (only up to N=2048 to avoid long wait)
        if N <= 2048:
            D = fourier_diff_matrix(N)
            t0 = time.perf_counter()
            for _ in range(n_repeats):
                _ = D @ u
            t1 = time.perf_counter()
            time_matrix.append((t1 - t0) / n_repeats)
        else:
            time_matrix.append(np.nan)

        # Time FFT method
        t0 = time.perf_counter()
        for _ in range(n_repeats):
            _ = fourier_diff_fft(u)
        t1 = time.perf_counter()
        time_fft.append((t1 - t0) / n_repeats)

    N_timing = np.array(N_timing)
    time_matrix = np.array(time_matrix)
    time_fft = np.array(time_fft)

    fig2, ax2 = plt.subplots(figsize=(6, 4.5))

    # Matrix timing (only where valid)
    valid = ~np.isnan(time_matrix)
    ax2.loglog(N_timing[valid], time_matrix[valid], 'o-', color=NAVY,
               markersize=5, linewidth=1.5, label=r'Matrix $O(N^2)$')
    ax2.loglog(N_timing, time_fft, 's-', color=CORAL,
               markersize=5, linewidth=1.5, label=r'FFT $O(N \log N)$')

    # Reference slopes
    N_ref = np.array([64, 2048], dtype=float)
    scale_n2 = time_matrix[valid][-1] * (N_ref / N_timing[valid][-1])**2
    ax2.loglog(N_ref, scale_n2, ':', color=SKY, linewidth=1.0, label=r'$\sim N^2$')

    N_ref2 = np.array([64, 8192], dtype=float)
    scale_nlogn = time_fft[-1] * (N_ref2 * np.log(N_ref2)) / (N_timing[-1] * np.log(N_timing[-1]))
    ax2.loglog(N_ref2, scale_nlogn, ':', color='#E67E22', linewidth=1.0,
               label=r'$\sim N \log N$')

    ax2.set_xlabel(r'$N$')
    ax2.set_ylabel('Time (seconds)')
    ax2.set_title('Computational Cost: Matrix vs FFT Differentiation')
    ax2.legend(loc='upper left')
    ax2.grid(True, alpha=0.3, which='both')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    fig2.tight_layout()
    fig2.savefig(OUTPUT_DIR / 'fourier_diff_timing.pdf',
                 bbox_inches='tight', pad_inches=0.05)
    print(f'Figure saved: {(OUTPUT_DIR / "fourier_diff_timing.pdf").resolve()}')
    plt.close(fig2)


if __name__ == '__main__':
    main()
