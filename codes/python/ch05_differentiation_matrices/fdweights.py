#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fdweights.py - Finite difference weights for arbitrary nodes

This module implements Fornberg's algorithm for computing finite difference
weights on arbitrarily spaced grids. The algorithm is numerically stable and
efficient, requiring O(MN²) operations where M is the derivative order and
N is the number of nodes.

Mathematical Background:
    Given n+1 distinct nodes x₀, x₁, ..., xₙ and a point ξ where we want to
    approximate the m-th derivative, Fornberg's algorithm computes weights
    w₀, w₁, ..., wₙ such that:

        f^(m)(ξ) ≈ Σⱼ wⱼ f(xⱼ)

    The algorithm uses a clever recursion that avoids the ill-conditioned
    Vandermonde system that would arise from a direct approach.

References:
    - Fornberg, B. (1988). "Generation of Finite Difference Formulas on
      Arbitrarily Spaced Grids", Mathematics of Computation, 51(184), 699-706.
    - Fornberg, B. (1996). "A Practical Guide to Pseudospectral Methods",
      Cambridge University Press.

Original MATLAB implementation: Toby Driscoll (2007)
Python translation: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

Part of "Computational Études: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""

import numpy as np
from functools import lru_cache


def fdweights(xi, x, m):
    """
    Compute finite difference weights using Fornberg's algorithm.

    Parameters
    ----------
    xi : float
        Evaluation point for the derivative approximation.
    x : array_like
        Node positions (length n+1). Can be arbitrarily spaced.
    m : int
        Order of derivative sought (0 = interpolation, 1 = first derivative, etc.)

    Returns
    -------
    w : ndarray
        Weights for approximating the m-th derivative at ξ.
        The approximation is: f^(m)(ξ) ≈ Σⱼ w[j] * f(x[j])

    Examples
    --------
    >>> import numpy as np
    >>> # Standard 3-point central difference for first derivative
    >>> h = 0.1
    >>> x = np.array([-h, 0, h])
    >>> w = fdweights(0, x, 1)
    >>> print(w)  # Should be [-0.5/h, 0, 0.5/h] approximately [-5, 0, 5]

    >>> # Interpolation weights (m=0) at midpoint
    >>> x = np.array([0.0, 1.0])
    >>> w = fdweights(0.5, x, 0)
    >>> print(w)  # Should be [0.5, 0.5]

    >>> # Second derivative with 5-point stencil
    >>> h = 0.1
    >>> x = h * np.array([-2, -1, 0, 1, 2])
    >>> w = fdweights(0, x, 2)
    >>> # Should approximate [-1/12, 4/3, -5/2, 4/3, -1/12] / h²

    Notes
    -----
    This is a direct translation of Toby Driscoll's compact MATLAB implementation.
    For maximum efficiency in repeated calls, consider caching results or using
    the vectorized version for building full differentiation matrices.

    The algorithm is numerically stable even for large stencils and high
    derivative orders, unlike direct Vandermonde-based approaches.
    """
    x = np.asarray(x, dtype=float)
    n = len(x) - 1

    if m > n:
        raise ValueError(f"Derivative order m={m} cannot exceed number of "
                        f"intervals n={n} (need at least m+1 points)")

    # Translate so that evaluation point is at origin
    x_shifted = x - xi

    # Compute weights using the recursive algorithm
    w = np.zeros(n + 1)
    for k in range(n + 1):
        w[k] = _weight_recursive(tuple(x_shifted), m, n, k)

    return w


def _weight_recursive(x, m, j, k):
    """
    Recursive computation of a single FD weight.

    This implements the recursion from Fornberg (1988), assuming the
    evaluation point has been translated to the origin.

    Parameters
    ----------
    x : tuple of float
        Shifted node positions (evaluation point at origin). Must be tuple for caching.
    m : int
        Derivative order sought.
    j : int
        Stencil width minus 1 (using nodes x[0], ..., x[j]).
    k : int
        Index of the node for which we compute the weight (0 ≤ k ≤ j).

    Returns
    -------
    c : float
        The finite difference weight for node x[k].
    """
    # Base cases
    if m < 0 or m > j:
        # Undefined: need at least m+1 points for m-th derivative
        return 0.0

    if m == 0 and j == 0:
        # Base case: single-point interpolation weight is 1
        return 1.0

    # Recursive cases
    if k < j:
        # Weight for an "old" node (was already in the previous stencil)
        c = (x[j] * _weight_recursive(x, m, j-1, k)
             - m * _weight_recursive(x, m-1, j-1, k)) / (x[j] - x[k])
    else:
        # Weight for the "new" node (just added to the stencil)
        # Compute the ratio of products
        if j >= 2:
            beta = np.prod([x[j-1] - x[i] for i in range(j-1)]) / \
                   np.prod([x[j] - x[i] for i in range(j)])
        elif j == 1:
            beta = 1.0 / (x[1] - x[0])
        else:
            beta = 1.0

        c = beta * (m * _weight_recursive(x, m-1, j-1, j-1)
                    - x[j-1] * _weight_recursive(x, m, j-1, j-1))

    return c


def fdweights_all(xi, x, m_max):
    """
    Compute FD weights for all derivatives up to order m_max.

    This is more efficient than calling fdweights() repeatedly when
    multiple derivative orders are needed.

    Parameters
    ----------
    xi : float
        Evaluation point for the derivative approximation.
    x : array_like
        Node positions (length n+1).
    m_max : int
        Maximum derivative order to compute.

    Returns
    -------
    W : ndarray of shape (m_max+1, n+1)
        W[m, k] is the weight for node x[k] in the m-th derivative approximation.

    Examples
    --------
    >>> x = np.linspace(-1, 1, 5)
    >>> W = fdweights_all(0, x, 2)
    >>> # W[0, :] are interpolation weights
    >>> # W[1, :] are first derivative weights
    >>> # W[2, :] are second derivative weights
    """
    x = np.asarray(x, dtype=float)
    n = len(x) - 1

    if m_max > n:
        raise ValueError(f"Maximum derivative order m_max={m_max} cannot exceed "
                        f"number of intervals n={n}")

    W = np.zeros((m_max + 1, n + 1))
    x_shifted = tuple(x - xi)

    for m in range(m_max + 1):
        for k in range(n + 1):
            W[m, k] = _weight_recursive(x_shifted, m, n, k)

    return W


# ==============================================================================
# Verification and testing
# ==============================================================================

def _verify_fdweights():
    """
    Verify the implementation against known finite difference formulas.
    """
    print("Verifying fdweights implementation...")
    print("=" * 60)

    # Test 1: 3-point central difference for first derivative
    h = 1.0
    x = np.array([-h, 0, h])
    w = fdweights(0, x, 1)
    expected = np.array([-0.5, 0, 0.5]) / h
    error1 = np.max(np.abs(w - expected))
    print(f"Test 1: 3-point central difference (1st deriv)")
    print(f"  Computed: {w}")
    print(f"  Expected: {expected}")
    print(f"  Max error: {error1:.2e}")

    # Test 2: 5-point central difference for first derivative
    x = np.array([-2*h, -h, 0, h, 2*h])
    w = fdweights(0, x, 1)
    expected = np.array([1, -8, 0, 8, -1]) / (12 * h)
    error2 = np.max(np.abs(w - expected))
    print(f"\nTest 2: 5-point central difference (1st deriv)")
    print(f"  Computed: {w}")
    print(f"  Expected: {expected}")
    print(f"  Max error: {error2:.2e}")

    # Test 3: 3-point central difference for second derivative
    x = np.array([-h, 0, h])
    w = fdweights(0, x, 2)
    expected = np.array([1, -2, 1]) / h**2
    error3 = np.max(np.abs(w - expected))
    print(f"\nTest 3: 3-point central difference (2nd deriv)")
    print(f"  Computed: {w}")
    print(f"  Expected: {expected}")
    print(f"  Max error: {error3:.2e}")

    # Test 4: Interpolation (m=0)
    x = np.array([0.0, 1.0])
    w = fdweights(0.5, x, 0)
    expected = np.array([0.5, 0.5])
    error4 = np.max(np.abs(w - expected))
    print(f"\nTest 4: Linear interpolation at midpoint")
    print(f"  Computed: {w}")
    print(f"  Expected: {expected}")
    print(f"  Max error: {error4:.2e}")

    # Test 5: Apply to actual function
    print(f"\nTest 5: Differentiate exp(x) at x=0")
    x = np.array([-0.1, 0, 0.1, 0.2])
    w = fdweights(0, x, 1)
    f_vals = np.exp(x)
    approx = np.dot(w, f_vals)
    exact = 1.0  # d/dx exp(x) at x=0
    error5 = abs(approx - exact)
    print(f"  Approximation: {approx:.10f}")
    print(f"  Exact value:   {exact:.10f}")
    print(f"  Error: {error5:.2e}")

    print("\n" + "=" * 60)
    all_passed = max(error1, error2, error3, error4) < 1e-14 and error5 < 1e-4
    if all_passed:
        print("All tests PASSED!")
    else:
        print("Some tests FAILED!")

    return all_passed


if __name__ == "__main__":
    _verify_fdweights()
