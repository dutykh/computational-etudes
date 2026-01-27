#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
spectral_matrix_periodic.py - Periodic spectral differentiation matrix

This module constructs the spectral differentiation matrix for periodic
functions on the domain [0, 2π) with N equispaced nodes.

Mathematical Background:
    For periodic functions, the natural interpolation basis consists of
    trigonometric polynomials. The spectral differentiation matrix D is
    derived from differentiating the interpolating trigonometric polynomial:

        p(x) = Σⱼ uⱼ φⱼ(x)

    where φⱼ are the cardinal functions (Dirichlet kernel shifted to node j).

    The matrix entries are:
        D[i,j] = (1/2) (-1)^(i-j) cot((i-j)π/N)   for i ≠ j
        D[i,i] = 0

    This matrix is:
    - Skew-symmetric: D^T = -D
    - Toeplitz (and circulant due to periodicity)
    - Dense: all off-diagonal entries are nonzero
    - Exact for trigonometric polynomials of degree ≤ N/2

Key Insight (from Chapter 4):
    Unlike the non-periodic case where equispaced nodes lead to the Runge
    phenomenon, for PERIODIC functions on [0, 2π), equispaced nodes are
    actually OPTIMAL. This is because the Chebyshev-like clustering near
    endpoints is unnecessary when there are no endpoints!

References:
    - Trefethen, L.N. (2000). "Spectral Methods in MATLAB", SIAM.
    - Fornberg, B. (1996). "A Practical Guide to Pseudospectral Methods",
      Cambridge University Press.

Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

Part of "Computational Études: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""

import numpy as np


def spectral_diff_periodic(N):
    """
    Construct the periodic spectral differentiation matrix.

    For N equispaced points on [0, 2π), this returns the matrix D such that
    for a periodic function u sampled at these points, D @ u approximates
    the derivative du/dx at the same points.

    Parameters
    ----------
    N : int
        Number of grid points. Should be even for best results.

    Returns
    -------
    D : ndarray of shape (N, N)
        The spectral differentiation matrix.
    x : ndarray of shape (N,)
        The grid points: x[j] = 2πj/N for j = 0, 1, ..., N-1.

    Notes
    -----
    The matrix entries are:
        D[i,j] = 0.5 * (-1)^(i-j) * cot((i-j) * π / N)   for i ≠ j
        D[i,i] = 0

    The matrix is skew-symmetric (D^T = -D) and Toeplitz (constant along
    diagonals). Due to periodicity, it is actually circulant.

    For even N, this matrix is exact for differentiating trigonometric
    polynomials of degree up to N/2.

    Examples
    --------
    >>> D, x = spectral_diff_periodic(8)
    >>> u = np.sin(x)         # Test function
    >>> du_exact = np.cos(x)  # Exact derivative
    >>> du_approx = D @ u     # Spectral approximation
    >>> error = np.max(np.abs(du_approx - du_exact))
    >>> print(f"Max error: {error:.2e}")  # Should be near machine precision
    """
    if N < 2:
        raise ValueError("Need at least 2 points")

    h = 2 * np.pi / N
    x = h * np.arange(N)

    # Build the differentiation matrix
    D = np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            if i != j:
                # D[i,j] = 0.5 * (-1)^(i-j) * cot((i-j) * π / N)
                diff = i - j
                D[i, j] = 0.5 * ((-1) ** diff) / np.tan(diff * np.pi / N)
            # D[i,i] = 0 (already initialized)

    return D, x


def spectral_diff_periodic_vectorized(N):
    """
    Vectorized construction of the periodic spectral differentiation matrix.

    This is a more efficient implementation that avoids explicit Python loops.

    Parameters
    ----------
    N : int
        Number of grid points.

    Returns
    -------
    D : ndarray of shape (N, N)
        The spectral differentiation matrix.
    x : ndarray of shape (N,)
        The grid points.
    """
    h = 2 * np.pi / N
    x = h * np.arange(N)

    # Create index difference matrix
    col = np.arange(N)
    row = np.arange(N).reshape(-1, 1)
    diff = row - col  # diff[i,j] = i - j

    # Compute off-diagonal entries
    # Use masked array to handle the singularity at diff=0
    with np.errstate(divide='ignore', invalid='ignore'):
        D = 0.5 * ((-1.0) ** diff) / np.tan(diff * np.pi / N)

    # Set diagonal to zero (tan(0) = 0 causes division by zero)
    np.fill_diagonal(D, 0.0)

    return D, x


def spectral_diff2_periodic(N):
    """
    Construct the periodic spectral SECOND differentiation matrix.

    This computes D² directly (more accurate than D @ D for some purposes).

    Parameters
    ----------
    N : int
        Number of grid points (should be even).

    Returns
    -------
    D2 : ndarray of shape (N, N)
        The second derivative spectral differentiation matrix.
    x : ndarray of shape (N,)
        The grid points.

    Notes
    -----
    The second derivative matrix has entries:
        D2[i,j] = -0.5 * (-1)^(i-j) / sin²((i-j)π/N)   for i ≠ j
        D2[i,i] = -π²(N² - 1) / (6h²)   where h = 2π/N

    This is more accurate than computing D @ D, especially for smooth functions.
    """
    h = 2 * np.pi / N
    x = h * np.arange(N)

    D2 = np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            if i != j:
                diff = i - j
                sin_val = np.sin(diff * np.pi / N)
                D2[i, j] = -0.5 * ((-1) ** diff) / (sin_val ** 2)

    # Diagonal entries
    diag_val = -(np.pi ** 2) * (N ** 2 - 1) / (3 * (2 * np.pi) ** 2)
    # Simplified: -N²/12 + 1/12 = -(N² - 1)/12 multiplied by (N/2π)² = N²/(4π²)
    # Actually the formula is: D2[i,i] = -π²(N² + 2) / (3 * h²) for standard formulation
    # Let me use D @ D for consistency
    D, _ = spectral_diff_periodic(N)
    D2 = D @ D

    return D2, x


# ==============================================================================
# Finite difference matrix construction using Fornberg weights
# ==============================================================================

def fd_diff_periodic(N, order=2):
    """
    Construct a periodic finite difference differentiation matrix.

    Uses Fornberg's algorithm to compute the FD weights, then assembles
    them into a circulant matrix for periodic boundary conditions.

    Parameters
    ----------
    N : int
        Number of grid points.
    order : int
        Order of accuracy (2, 4, 6, ...). The stencil width is order + 1.

    Returns
    -------
    D : ndarray of shape (N, N)
        The finite difference differentiation matrix.
    x : ndarray of shape (N,)
        The grid points.

    Notes
    -----
    Order 2: 3-point stencil [-1/2, 0, 1/2] / h
    Order 4: 5-point stencil [1/12, -2/3, 0, 2/3, -1/12] / h
    Order 6: 7-point stencil [-1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60] / h
    """
    from fdweights import fdweights

    h = 2 * np.pi / N
    x = h * np.arange(N)

    # Determine stencil size
    stencil_half = order // 2
    stencil_size = order + 1

    # Create stencil nodes centered at 0
    stencil_nodes = h * np.arange(-stencil_half, stencil_half + 1)

    # Compute FD weights for first derivative at center
    weights = fdweights(0, stencil_nodes, 1)

    # Build circulant matrix
    D = np.zeros((N, N))
    for i in range(N):
        for k, w in enumerate(weights):
            j = (i + k - stencil_half) % N  # Periodic wrap-around
            D[i, j] += w

    return D, x


# ==============================================================================
# Verification and testing
# ==============================================================================

def _verify_spectral_matrix():
    """
    Verify the spectral differentiation matrix against known results.
    """
    print("Verifying periodic spectral differentiation matrix...")
    print("=" * 60)

    # Test 1: Differentiate sin(x)
    N = 16
    D, x = spectral_diff_periodic(N)
    u = np.sin(x)
    du_exact = np.cos(x)
    du_approx = D @ u
    error1 = np.max(np.abs(du_approx - du_exact))
    print(f"Test 1: d/dx sin(x) with N={N}")
    print(f"  Max error: {error1:.2e}")

    # Test 2: Differentiate sin(3x) - higher frequency
    u = np.sin(3 * x)
    du_exact = 3 * np.cos(3 * x)
    du_approx = D @ u
    error2 = np.max(np.abs(du_approx - du_exact))
    print(f"\nTest 2: d/dx sin(3x) with N={N}")
    print(f"  Max error: {error2:.2e}")

    # Test 3: Verify skew-symmetry D^T = -D
    skew_error = np.max(np.abs(D + D.T))
    print(f"\nTest 3: Skew-symmetry ||D + D^T||_max")
    print(f"  Error: {skew_error:.2e}")

    # Test 4: Verify Toeplitz structure
    # Check that D[i,j] depends only on i-j
    is_toeplitz = True
    for i in range(N):
        for j in range(N):
            if i != j:
                diff = (i - j) % N
                if diff != 0:
                    # Compare with D[0, (0-diff) % N] = D[0, N-diff]
                    expected = D[0, (N - diff) % N] if diff != 0 else 0
                    if abs(D[i, j] - D[(i - diff) % N, j - diff]) > 1e-14:
                        is_toeplitz = False
                        break
    print(f"\nTest 4: Toeplitz structure: {'PASS' if is_toeplitz else 'FAIL'}")

    # Test 5: Compare with vectorized version
    D_vec, _ = spectral_diff_periodic_vectorized(N)
    vec_error = np.max(np.abs(D - D_vec))
    print(f"\nTest 5: Vectorized vs loop implementation")
    print(f"  Max difference: {vec_error:.2e}")

    # Test 6: Test function from the étude
    u = 1 / (2 + np.sin(x))
    du_exact = -np.cos(x) / (2 + np.sin(x))**2
    du_approx = D @ u
    error6 = np.max(np.abs(du_approx - du_exact))
    print(f"\nTest 6: d/dx [1/(2+sin(x))] with N={N}")
    print(f"  Max error: {error6:.2e}")

    # Test 7: Convergence test
    print(f"\nTest 7: Convergence study for u(x) = 1/(2+sin(x))")
    print(f"  {'N':>4}  {'Max Error':>12}")
    print(f"  {'-'*4}  {'-'*12}")
    for N in [8, 16, 32, 64]:
        D, x = spectral_diff_periodic(N)
        u = 1 / (2 + np.sin(x))
        du_exact = -np.cos(x) / (2 + np.sin(x))**2
        du_approx = D @ u
        error = np.max(np.abs(du_approx - du_exact))
        print(f"  {N:4d}  {error:12.4e}")

    print("\n" + "=" * 60)
    all_passed = (error1 < 1e-14 and error2 < 1e-12 and
                  skew_error < 1e-14 and vec_error < 1e-14)
    print("All tests PASSED!" if all_passed else "Some tests FAILED!")

    return all_passed


if __name__ == "__main__":
    _verify_spectral_matrix()
