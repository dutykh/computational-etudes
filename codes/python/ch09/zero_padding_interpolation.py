#!/usr/bin/env python3
"""
zero_padding_interpolation.py
Chapter 9: Physical and Fourier Space on Grids

Demonstrates band-limited interpolation via zero-padding in Fourier space.
This is Etude 4 from the chapter.

Author: Dr. Denys Dutykh
Date: February 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def zero_pad_interpolate(v, q=4):
    """
    Interpolate a periodic function by zero-padding in Fourier space.

    Parameters:
        v : array - Function values on coarse grid (length N)
        q : int - Upsampling factor (fine grid has q*N points)

    Returns:
        v_fine : array - Interpolated values on fine grid
    """
    N = len(v)
    M = q * N

    # Forward FFT
    v_hat = np.fft.fft(v)

    # Zero-pad in Fourier space
    v_hat_padded = np.zeros(M, dtype=complex)

    # Copy low frequencies
    # Positive frequencies: 0, 1, ..., N/2-1
    v_hat_padded[:N//2] = v_hat[:N//2]
    # Negative frequencies: -N/2+1, ..., -1
    v_hat_padded[-(N//2-1):] = v_hat[-(N//2-1):]

    # Handle Nyquist mode: split between +N/2 and -N/2
    v_hat_padded[N//2] = v_hat[N//2] / 2
    v_hat_padded[M - N//2] = v_hat[N//2] / 2

    # Inverse FFT and scale
    v_fine = np.real(np.fft.ifft(v_hat_padded)) * q

    return v_fine


def main():
    # Coarse grid parameters
    N = 32
    x_coarse = 2 * np.pi * np.arange(N) / N

    # Test function: exp(sin(x))
    v = np.exp(np.sin(x_coarse))

    # Interpolate via zero-padding
    q = 4
    M = q * N
    x_fine = 2 * np.pi * np.arange(M) / M
    v_fine = zero_pad_interpolate(v, q)

    # Exact function for comparison
    v_exact = np.exp(np.sin(x_fine))

    # Compute interpolation error
    error = np.max(np.abs(v_fine - v_exact))
    print(f"N = {N} coarse points, M = {M} fine points")
    print(f"Maximum interpolation error: {error:.2e}")

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 4))

    # Plot interpolant
    ax.plot(x_fine, v_fine, 'b-', linewidth=1.2, label='Band-limited interpolant')

    # Plot coarse samples
    ax.plot(x_coarse, v, 'ko', markersize=8, markerfacecolor='white',
            markeredgewidth=2, label=f'Coarse samples (N = {N})')

    ax.set_xlabel(r'$x$', fontsize=12)
    ax.set_ylabel(r'$\exp(\sin x)$', fontsize=12)
    ax.set_title('Periodic band-limited interpolation via FFT zero-padding', fontsize=12)
    ax.legend(loc='upper right', fontsize=10)
    ax.set_xlim(0, 2*np.pi)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save figure
    output_dir = Path(__file__).parent.parent.parent.parent / 'textbook' / 'figures' / 'ch09' / 'python'
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / 'zero_padding_interpolation.pdf', bbox_inches='tight')
    print(f"Figure saved to {output_dir / 'zero_padding_interpolation.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
