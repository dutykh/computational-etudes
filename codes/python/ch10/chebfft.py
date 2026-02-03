#!/usr/bin/env python3
"""
chebfft.py
Chapter 10: Spectral PDE Solvers with Chebyshev and Fourier Grids

Chebyshev differentiation via FFT.
This module provides the core function for computing derivatives on
Chebyshev grids using the FFT.

Author: Dr. Denys Dutykh
Date: February 2026
"""

import numpy as np


def chebfft(v):
    """
    Compute the first derivative of a function sampled at Chebyshev points
    using the FFT.

    The algorithm transforms from the x-variable to the theta-variable
    (where x = cos(theta)), performs differentiation via FFT in theta-space,
    and transforms back.

    Parameters:
        v : array of shape (N+1,)
            Function values at Chebyshev points x_j = cos(j*pi/N), j=0,...,N

    Returns:
        w : array of shape (N+1,)
            Approximate derivative values at the same Chebyshev points
    """
    v = np.asarray(v, dtype=float)
    N = len(v) - 1

    if N == 0:
        return np.array([0.0])

    x = np.cos(np.pi * np.arange(N + 1) / N)

    # Extend to even periodic function in theta
    # V has length 2N: V = [v_0, v_1, ..., v_N, v_{N-1}, ..., v_1]
    V = np.concatenate([v, v[-2:0:-1]])

    # FFT to get Fourier coefficients
    U = np.fft.fft(V).real

    # Wavenumber vector for differentiation
    # k = [0, 1, 2, ..., N-1, 0, -(N-1), -(N-2), ..., -1]
    k = np.concatenate([np.arange(N), [0], np.arange(1-N, 0)])

    # Differentiate in theta: multiply by ik
    W = np.fft.ifft(1j * k * np.fft.fft(V)).real

    # Transform from theta derivative to x derivative
    # q'(x) = Q'(theta) / (dx/dtheta) = -Q'(theta) / sin(theta)
    #       = -Q'(theta) / sqrt(1 - x^2)
    w = np.zeros(N + 1)
    w[1:N] = -W[1:N] / np.sqrt(1 - x[1:N]**2)

    # Boundary values via summation formulas (l'Hopital's rule at x = +/- 1)
    ii = np.arange(N)
    w[0] = np.sum(ii**2 * U[ii]) / N + 0.5 * N * U[N]
    w[N] = np.sum((-1)**(ii+1) * ii**2 * U[ii]) / N + 0.5 * (-1)**(N+1) * N * U[N]

    return w


def cheb_matrix(N):
    """
    Construct the Chebyshev differentiation matrix D_N.

    Parameters:
        N : int
            Number of intervals (matrix is (N+1) x (N+1))

    Returns:
        D : array of shape (N+1, N+1)
            Chebyshev differentiation matrix
        x : array of shape (N+1,)
            Chebyshev points x_j = cos(j*pi/N), j=0,...,N
    """
    if N == 0:
        return np.array([[0.0]]), np.array([1.0])

    # Chebyshev-Gauss-Lobatto points
    x = np.cos(np.pi * np.arange(N + 1) / N)

    # Weights: c_0 = c_N = 2, others = 1
    c = np.ones(N + 1)
    c[0] = 2.0
    c[N] = 2.0

    # Build the matrix
    X = np.tile(x, (N + 1, 1))
    dX = X - X.T
    C = np.outer(c, 1.0 / c)
    sign = np.outer((-1.0)**np.arange(N + 1), (-1.0)**np.arange(N + 1))

    with np.errstate(divide='ignore', invalid='ignore'):
        D = C * sign / dX

    # Fix diagonal using negative sum trick
    np.fill_diagonal(D, 0.0)
    D[np.diag_indices(N + 1)] = -np.sum(D, axis=1)

    return D, x


if __name__ == '__main__':
    # Test with a smooth function
    N = 20
    x = np.cos(np.pi * np.arange(N + 1) / N)

    # Test function: f(x) = exp(x) * cos(4x)
    f = np.exp(x) * np.cos(4 * x)
    f_prime_exact = np.exp(x) * (np.cos(4 * x) - 4 * np.sin(4 * x))

    # Compute derivative via FFT
    f_prime_fft = chebfft(f)

    # Compute derivative via matrix
    D, _ = cheb_matrix(N)
    f_prime_mat = D @ f

    # Compare errors
    err_fft = np.max(np.abs(f_prime_fft - f_prime_exact))
    err_mat = np.max(np.abs(f_prime_mat - f_prime_exact))

    print(f"N = {N}")
    print(f"FFT method error:    {err_fft:.2e}")
    print(f"Matrix method error: {err_mat:.2e}")
