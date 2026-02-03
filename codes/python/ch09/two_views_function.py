#!/usr/bin/env python3
"""
two_views_function.py
Chapter 9: Physical and Fourier Space on Grids

Demonstrates the two views of a function: physical space and Fourier space.
Shows a smooth function and its rapidly decaying Fourier transform.

Author: Dr. Denys Dutykh
Date: February 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def compute_spectrum(u_func, L=20.0, N=1024):
    """
    Approximate the Fourier transform of u_func via FFT.

    Parameters:
        u_func : callable - Function to transform
        L : float - Domain is [-L/2, L/2]
        N : int - Number of grid points

    Returns:
        k : array - Wavenumbers
        u_hat : array - Fourier coefficients
    """
    x = np.linspace(-L/2, L/2, N, endpoint=False)
    dx = x[1] - x[0]
    u = u_func(x)

    # Approximate continuous FT via FFT
    k = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    u_hat = dx * np.fft.fft(u)

    return x, u, k, u_hat


def main():
    # Define the test function: Gaussian-modulated cosine
    def u(x):
        return np.exp(-x**2) * np.cos(3*x)

    # Compute physical and Fourier representations
    x, u_vals, k, u_hat = compute_spectrum(u, L=20.0, N=2**10)

    # Shift for plotting (center k=0)
    k_shifted = np.fft.fftshift(k)
    u_hat_shifted = np.fft.fftshift(u_hat)

    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.5))

    # Left: Physical space
    axes[0].plot(x, u_vals, 'b-', linewidth=1.2)
    axes[0].set_xlabel(r'$x$', fontsize=12)
    axes[0].set_ylabel(r'$u(x)$', fontsize=12)
    axes[0].set_title('Physical space', fontsize=12)
    axes[0].set_xlim(-8, 8)
    axes[0].grid(True, alpha=0.3)
    axes[0].axhline(y=0, color='k', linewidth=0.5)

    # Right: Fourier space (log scale)
    axes[1].semilogy(k_shifted, np.abs(u_hat_shifted), 'b-', linewidth=1.2)
    axes[1].set_xlabel(r'$k$', fontsize=12)
    axes[1].set_ylabel(r'$|\hat{u}(k)|$', fontsize=12)
    axes[1].set_title('Fourier space (approximate spectrum)', fontsize=12)
    axes[1].set_xlim(-15, 15)
    axes[1].set_ylim(1e-16, 10)
    axes[1].grid(True, alpha=0.3, which='both')

    plt.tight_layout()

    # Save figure
    output_dir = Path(__file__).parent.parent.parent.parent / 'textbook' / 'figures' / 'ch09' / 'python'
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / 'two_views_function.pdf', bbox_inches='tight')
    print(f"Figure saved to {output_dir / 'two_views_function.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
