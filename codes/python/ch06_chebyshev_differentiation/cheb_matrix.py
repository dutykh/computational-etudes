#!/usr/bin/env python3
"""
cheb_matrix.py

Core Chebyshev differentiation matrix construction for spectral methods
on bounded (non-periodic) domains.

The Chebyshev-Gauss-Lobatto points are:
    x_j = cos(j*pi/N),  j = 0, 1, ..., N

These cluster near the boundaries, avoiding Runge phenomenon and enabling
spectral accuracy for smooth functions on [-1, 1].

The differentiation matrix D_N maps function values {v_0, ..., v_N} to
derivative values {w_0, ..., w_N} via polynomial interpolation and
differentiation at the collocation points.

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""

import numpy as np


def cheb_matrix(N):
    """
    Construct the Chebyshev differentiation matrix and grid points.

    Uses the explicit formulas from Trefethen's "Spectral Methods in MATLAB"
    with the "negative sum trick" for diagonal entries to ensure numerical
    stability (derivative of a constant is exactly zero).

    Parameters
    ----------
    N : int
        Number of intervals (matrix size is (N+1) x (N+1))

    Returns
    -------
    D : ndarray, shape (N+1, N+1)
        Chebyshev differentiation matrix
    x : ndarray, shape (N+1,)
        Chebyshev-Gauss-Lobatto points on [-1, 1], ordered from x_0=1 to x_N=-1

    Notes
    -----
    The matrix entries (before applying negative sum trick) are:

    Corner entries:
        D_{00} = (2N^2 + 1) / 6
        D_{NN} = -(2N^2 + 1) / 6

    Diagonal entries (j = 1, ..., N-1):
        D_{jj} = -x_j / (2(1 - x_j^2))

    Off-diagonal entries:
        D_{ij} = (c_i / c_j) * (-1)^{i+j} / (x_i - x_j)

    where c_0 = c_N = 2 and c_j = 1 for j = 1, ..., N-1.

    For numerical stability, we compute diagonals via:
        D_{jj} = -sum_{k != j} D_{jk}

    This ensures that D @ ones(N+1) = zeros(N+1) to machine precision.

    Examples
    --------
    >>> D, x = cheb_matrix(4)
    >>> x
    array([ 1.        ,  0.70710678,  0.        , -0.70710678, -1.        ])
    >>> np.allclose(D @ np.ones(5), np.zeros(5))
    True
    """
    if N == 0:
        return np.array([[0.0]]), np.array([1.0])

    # Chebyshev-Gauss-Lobatto points
    x = np.cos(np.pi * np.arange(N + 1) / N)

    # Coefficient vector c: c_0 = c_N = 2, others = 1
    c = np.ones(N + 1)
    c[0] = 2.0
    c[N] = 2.0

    # Construct the differentiation matrix
    # Use outer products for vectorized computation
    X = np.tile(x, (N + 1, 1))  # X[i,j] = x[j]
    dX = X - X.T                 # dX[i,j] = x[j] - x[i]

    # Off-diagonal entries: D[i,j] = c[i]/c[j] * (-1)^(i+j) / (x[i] - x[j])
    C = np.outer(c, 1.0 / c)     # C[i,j] = c[i] / c[j]

    # Sign matrix: (-1)^(i+j)
    i_idx = np.arange(N + 1)
    sign_mat = np.outer((-1.0) ** i_idx, (-1.0) ** i_idx)

    # Compute off-diagonal entries (avoid division by zero on diagonal)
    with np.errstate(divide='ignore', invalid='ignore'):
        D = C * sign_mat / (-dX)

    # Apply negative sum trick for diagonal entries
    # This ensures numerical stability: D @ ones = zeros
    np.fill_diagonal(D, 0.0)
    D[np.diag_indices(N + 1)] = -np.sum(D, axis=1)

    return D, x


def cheb_matrix_explicit(N):
    """
    Construct Chebyshev matrix using explicit formulas (for verification).

    This version uses the direct formulas for all entries including diagonals,
    without the negative sum trick. Useful for comparing with the stable version.

    Parameters
    ----------
    N : int
        Number of intervals

    Returns
    -------
    D : ndarray, shape (N+1, N+1)
        Chebyshev differentiation matrix
    x : ndarray, shape (N+1,)
        Chebyshev-Gauss-Lobatto grid points
    """
    if N == 0:
        return np.array([[0.0]]), np.array([1.0])

    x = np.cos(np.pi * np.arange(N + 1) / N)
    D = np.zeros((N + 1, N + 1))

    # Coefficient vector c
    c = np.ones(N + 1)
    c[0] = 2.0
    c[N] = 2.0

    for i in range(N + 1):
        for j in range(N + 1):
            if i != j:
                D[i, j] = (c[i] / c[j]) * ((-1) ** (i + j)) / (x[i] - x[j])

    # Explicit diagonal formulas
    D[0, 0] = (2 * N**2 + 1) / 6.0
    D[N, N] = -(2 * N**2 + 1) / 6.0

    for j in range(1, N):
        D[j, j] = -x[j] / (2 * (1 - x[j]**2))

    return D, x


def cheb_second_derivative_matrix(N):
    """
    Construct the second derivative matrix D^2 via matrix squaring.

    For second-order differential equations, we need D^2. While explicit
    formulas exist, matrix squaring D @ D is simpler and adequate for
    spectral N values (typically N < 200).

    Parameters
    ----------
    N : int
        Number of intervals

    Returns
    -------
    D2 : ndarray, shape (N+1, N+1)
        Second derivative matrix D^2
    D : ndarray, shape (N+1, N+1)
        First derivative matrix D
    x : ndarray, shape (N+1,)
        Chebyshev grid points
    """
    D, x = cheb_matrix(N)
    D2 = D @ D
    return D2, D, x


def verify_small_matrices():
    """
    Verify the matrix construction against known small cases.

    D_1 (2x2 matrix) and D_2 (3x3 matrix) have simple exact forms
    that can be computed by hand.

    Returns
    -------
    all_pass : bool
        True if all verification tests pass
    """
    all_pass = True

    # D_1: 2x2 matrix
    # x = [1, -1], linear interpolation derivative
    D1, x1 = cheb_matrix(1)
    D1_exact = np.array([[0.5, -0.5],
                         [0.5, -0.5]])
    if not np.allclose(D1, D1_exact):
        print("D_1 verification FAILED")
        print(f"Computed:\n{D1}")
        print(f"Expected:\n{D1_exact}")
        all_pass = False
    else:
        print("D_1 verification PASSED")

    # D_2: 3x3 matrix
    # x = [1, 0, -1], quadratic interpolation
    D2, x2 = cheb_matrix(2)
    D2_exact = np.array([[1.5, -2.0, 0.5],
                         [0.5, 0.0, -0.5],
                         [-0.5, 2.0, -1.5]])
    if not np.allclose(D2, D2_exact):
        print("D_2 verification FAILED")
        print(f"Computed:\n{D2}")
        print(f"Expected:\n{D2_exact}")
        all_pass = False
    else:
        print("D_2 verification PASSED")

    # Verify negative sum trick: D @ ones = zeros
    for N in [4, 8, 16, 32]:
        D, x = cheb_matrix(N)
        result = D @ np.ones(N + 1)
        max_error = np.max(np.abs(result))
        if max_error > 1e-12:
            print(f"Negative sum trick FAILED for N={N}: max error = {max_error:.2e}")
            all_pass = False
        else:
            print(f"Negative sum trick PASSED for N={N}: max error = {max_error:.2e}")

    return all_pass


def demo_differentiation():
    """
    Demonstrate spectral differentiation accuracy.

    Uses the "Witch of Agnesi" function u(x) = 1/(1 + 4x^2) as a test case.
    This is smooth on [-1, 1] and has a well-defined derivative.
    """
    # Test function: Witch of Agnesi
    u = lambda x: 1.0 / (1.0 + 4.0 * x**2)
    u_exact = lambda x: -8.0 * x / (1.0 + 4.0 * x**2)**2

    print("\nDifferentiation accuracy for u(x) = 1/(1 + 4x^2):")
    print("-" * 50)
    print(f"{'N':>6} {'Max Error':>14}")
    print("-" * 50)

    for N in [4, 8, 16, 32, 64]:
        D, x = cheb_matrix(N)
        v = u(x)
        w = D @ v
        w_exact = u_exact(x)
        error = np.max(np.abs(w - w_exact))
        print(f"{N:>6d} {error:>14.6e}")

    print("-" * 50)


if __name__ == '__main__':
    print("=" * 60)
    print("Chebyshev Differentiation Matrix Verification")
    print("=" * 60)

    # Verify small matrices
    print("\n1. Verifying small matrices (D_1, D_2)...")
    verify_small_matrices()

    # Demonstrate differentiation accuracy
    print("\n2. Demonstrating differentiation accuracy...")
    demo_differentiation()

    print("\n" + "=" * 60)
    print("Verification complete!")
    print("=" * 60)
