#!/usr/bin/env python3
"""
zero_padding_interpolation.py
Chapter 9: Physical and Fourier Space on Grids

Demonstrates band-limited interpolation via zero-padding in Fourier space.
This is Etude 5 from the chapter.

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

    # Dense grid for the true function
    x_dense = np.linspace(0, 2 * np.pi, 500)
    v_true = np.exp(np.sin(x_dense))

    # Exact function at interpolated points for error
    v_exact = np.exp(np.sin(x_fine))

    # Compute interpolation error
    error = np.max(np.abs(v_fine - v_exact))
    print(f"N = {N} coarse points, M = {M} fine points")
    print(f"Maximum interpolation error: {error:.2e}")

    # Create two-panel figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 5.5),
                                    gridspec_kw={'height_ratios': [3, 1]})

    # --- Top panel: function, interpolant, samples ---
    ax1.plot(x_dense, v_true, '--', color='#E74C3C', linewidth=1.0,
             label=r'True $\exp(\sin x)$')
    ax1.plot(x_fine, v_fine, '-', color='#142D6E', linewidth=1.2,
             label=f'Band-limited interpolant ($M = {M}$)')
    ax1.plot(x_coarse, v, 'ko', markersize=7, markerfacecolor='white',
             markeredgewidth=1.5, label=f'Coarse samples ($N = {N}$)', zorder=5)

    ax1.set_ylabel(r'$\exp(\sin x)$', fontsize=12)
    ax1.set_title('Periodic band-limited interpolation via FFT zero-padding',
                  fontsize=12)
    ax1.legend(loc='upper right', fontsize=9)
    ax1.set_xlim(0, 2 * np.pi)
    ax1.grid(True, alpha=0.3)

    # --- Bottom panel: pointwise error ---
    ax2.semilogy(x_fine, np.abs(v_fine - v_exact) + 1e-16, '-',
                 color='#142D6E', linewidth=1.0)
    ax2.set_xlabel(r'$x$', fontsize=12)
    ax2.set_ylabel('Pointwise error', fontsize=11)
    ax2.set_xlim(0, 2 * np.pi)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save figure
    output_dir = Path(__file__).parent.parent.parent.parent / 'textbook' / 'figures' / 'ch09' / 'python'
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / 'zero_padding_interpolation.pdf', bbox_inches='tight')
    print(f"Figure saved to {output_dir / 'zero_padding_interpolation.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
