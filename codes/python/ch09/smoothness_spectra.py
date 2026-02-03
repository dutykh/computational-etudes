#!/usr/bin/env python3
"""
smoothness_spectra.py
Chapter 9: Physical and Fourier Space on Grids

Demonstrates the relationship between function smoothness and spectral decay.
Shows spectra of: exp(sin(x)), periodic hat, and square wave.

Author: Dr. Denys Dutykh
Date: February 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def compute_spectrum(f_func, N=256):
    """
    Compute Fourier spectrum of a periodic function on [0, 2*pi).

    Parameters:
        f_func : callable - Function to analyze
        N : int - Number of grid points

    Returns:
        k : array - Wavenumbers (positive only)
        f_hat_mag : array - Magnitude of Fourier coefficients
    """
    x = 2 * np.pi * np.arange(N) / N
    f = f_func(x)
    f_hat = np.fft.fft(f) / N

    # Take positive wavenumbers only (excluding k=0)
    k = np.arange(1, N//2)
    f_hat_mag = np.abs(f_hat[1:N//2])

    return k, f_hat_mag


def main():
    N = 256

    # Define three functions with different smoothness

    # 1. Analytic: exp(sin(x))
    def f1(x):
        return np.exp(np.sin(x))

    # 2. C^0 but not C^1: periodic hat function centered at pi
    def f2(x):
        # Hat function: peak at x=pi, zero at x=0 and x=2*pi
        return np.maximum(0, 1 - np.abs(x - np.pi) / (np.pi/2))

    # 3. Discontinuous: square wave
    def f3(x):
        return np.sign(np.sin(x))

    # Compute spectra
    k1, spec1 = compute_spectrum(f1, N)
    k2, spec2 = compute_spectrum(f2, N)
    k3, spec3 = compute_spectrum(f3, N)

    # Create figure with 3 panels
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))

    # Panel 1: exp(sin(x)) - semilog plot (exponential decay)
    axes[0].semilogy(k1, spec1, 'b.-', markersize=4, linewidth=0.8)
    axes[0].set_xlabel(r'$k$', fontsize=11)
    axes[0].set_ylabel(r'$|\hat{f}_k|$', fontsize=11)
    axes[0].set_title(r'$\exp(\sin x)$ (analytic)', fontsize=11)
    axes[0].set_xlim(0, N//2)
    axes[0].set_ylim(1e-16, 10)
    axes[0].grid(True, alpha=0.3, which='both')

    # Panel 2: Hat function - log-log plot (algebraic decay ~k^{-2})
    axes[1].loglog(k2, spec2, 'b.-', markersize=4, linewidth=0.8)
    # Add reference line for k^{-2}
    k_ref = np.array([5.0, 100.0])
    axes[1].loglog(k_ref, 0.5 * k_ref**(-2), 'r--', linewidth=1.5, label=r'$k^{-2}$')
    axes[1].set_xlabel(r'$k$', fontsize=11)
    axes[1].set_ylabel(r'$|\hat{f}_k|$', fontsize=11)
    axes[1].set_title(r'Hat function ($C^0$)', fontsize=11)
    axes[1].set_xlim(1, N//2)
    axes[1].legend(loc='upper right', fontsize=9)
    axes[1].grid(True, alpha=0.3, which='both')

    # Panel 3: Square wave - log-log plot (algebraic decay ~k^{-1})
    axes[2].loglog(k3, spec3, 'b.-', markersize=4, linewidth=0.8)
    # Add reference line for k^{-1}
    k_ref = np.array([1.0, 100.0])
    axes[2].loglog(k_ref, 1.5 * k_ref**(-1), 'r--', linewidth=1.5, label=r'$k^{-1}$')
    axes[2].set_xlabel(r'$k$', fontsize=11)
    axes[2].set_ylabel(r'$|\hat{f}_k|$', fontsize=11)
    axes[2].set_title('Square wave (discontinuous)', fontsize=11)
    axes[2].set_xlim(1, N//2)
    axes[2].legend(loc='upper right', fontsize=9)
    axes[2].grid(True, alpha=0.3, which='both')

    plt.tight_layout()

    # Save figure
    output_dir = Path(__file__).parent.parent.parent.parent / 'textbook' / 'figures' / 'ch09' / 'python'
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / 'smoothness_spectra.pdf', bbox_inches='tight')
    print(f"Figure saved to {output_dir / 'smoothness_spectra.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
