#!/usr/bin/env python3
"""
fft_aliasing.py
Chapter 9: Physical and Fourier Space on Grids

Demonstrates aliasing in FFT of high-frequency data.
sin(17x) sampled at N=32 points shows spike at k=-15 due to aliasing.
This is Etude 3 from the chapter.

Author: Dr. Denys Dutykh
Date: February 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def main():
    # Grid parameters
    N = 32
    x = 2 * np.pi * np.arange(N) / N

    # High-frequency function: k=17 > Nyquist limit N/2=16
    u = np.sin(17 * x)

    # Compute FFT
    u_hat = np.fft.fft(u) / N

    # Wavenumber vector
    k = np.arange(N)
    k[k > N//2] -= N  # Convert to -N/2+1, ..., N/2

    # Shift for plotting
    k_shifted = np.fft.fftshift(k)
    u_hat_shifted = np.fft.fftshift(u_hat)

    # Find the aliased location
    print(f"N = {N}, Nyquist limit = {N//2}")
    print(f"True wavenumber: 17")
    print(f"Alias: 17 = -15 + 32, so appears at k = -15")
    print(f"Maximum |û_k| at k = {k_shifted[np.argmax(np.abs(u_hat_shifted))]}")

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 4))

    # Stem plot of spectrum
    markerline, stemlines, baseline = ax.stem(k_shifted, np.abs(u_hat_shifted),
                                               linefmt='b-', markerfmt='bo',
                                               basefmt='k-')
    plt.setp(stemlines, linewidth=1.5)
    plt.setp(markerline, markersize=6)

    # Highlight the aliased peak
    aliased_k = -15
    aliased_idx = np.where(k_shifted == aliased_k)[0][0]
    ax.plot(aliased_k, np.abs(u_hat_shifted[aliased_idx]), 'ro', markersize=10,
            label=f'Aliased peak at k = {aliased_k}')

    # Add annotation
    ax.annotate(f'$k = 17$ aliases to\n$k = {aliased_k}$',
                xy=(aliased_k, np.abs(u_hat_shifted[aliased_idx])),
                xytext=(aliased_k + 5, 0.4),
                fontsize=10,
                arrowprops=dict(arrowstyle='->', color='red'))

    ax.set_xlabel(r'$k$ (wavenumber)', fontsize=12)
    ax.set_ylabel(r'$|\hat{u}_k|$', fontsize=12)
    ax.set_title(r'FFT of $\sin(17x)$ with $N = 32$ points: aliasing to $k = -15$', fontsize=12)
    ax.set_xlim(-N//2 - 1, N//2 + 1)
    ax.set_ylim(0, 0.6)
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save figure
    output_dir = Path(__file__).parent.parent.parent.parent / 'textbook' / 'figures' / 'ch09' / 'python'
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / 'fft_aliasing.pdf', bbox_inches='tight')
    print(f"Figure saved to {output_dir / 'fft_aliasing.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
