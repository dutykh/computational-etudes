#!/usr/bin/env python3
"""
spectral_derivative_accuracy.py
Chapter 9: Physical and Fourier Space on Grids

Demonstrates spectral accuracy of FFT-based differentiation.
Tests f(x) = exp(sin(x)) whose exact derivative is cos(x)*exp(sin(x)).
This is Etude 3 from the chapter.

Author: Dr. Denys Dutykh
Date: February 2026
"""

import numpy as np


def spectral_derivative(v):
    """Compute spectral derivative of periodic function."""
    N = len(v)
    v_hat = np.fft.fft(v)
    k = np.fft.fftfreq(N) * N  # Wavenumbers: 0,1,...,N/2,-N/2+1,...,-1
    k = 1j * k
    k[N//2] = 0  # Zero out Nyquist mode for odd derivatives
    w_hat = k * v_hat
    return np.real(np.fft.ifft(w_hat))


def main():
    # Test function: f(x) = exp(sin(x)), f'(x) = cos(x)*exp(sin(x))
    f = lambda x: np.exp(np.sin(x))
    f_prime = lambda x: np.cos(x) * np.exp(np.sin(x))

    # Grid sizes to test
    N_values = [4, 8, 12, 16, 20, 24, 28, 32]

    print(f"{'N':>4s}   {'Max error':>12s}")
    print("-" * 20)

    for N in N_values:
        x = 2 * np.pi * np.arange(N) / N
        v = f(x)
        w = spectral_derivative(v)
        w_exact = f_prime(x)
        error = np.max(np.abs(w - w_exact))
        print(f"{N:4d}   {error:12.2e}")


if __name__ == '__main__':
    main()
