#!/usr/bin/env python3
"""
quad_exactness_table.py
Chapter 15: Quadrature in Spectral Methods

Computational Etude 15.3: Exactness of Newton--Cotes Rules.

Demonstrates the stark difference in exactness between the monomial basis
{1, x, x^2, ...} and the Chebyshev basis {T_0, T_1, T_2, ...} for
Newton--Cotes quadrature with n = 30 equispaced nodes.

For even k from 26 to 38, this script computes and prints:

    |E_n(x^k)| = |integral_{-1}^{1} x^k dx - sum_j w_j x_j^k|
    |E_n(T_k)| = |integral_{-1}^{1} T_k(x) dx - sum_j w_j T_k(x_j)|

The exact integrals are:
    integral x^k dx  = 2/(k+1)  for even k,  0 for odd k,
    integral T_k dx   = 2/(1-k^2) for even k, 0 for odd k.

The key observation is that Newton--Cotes with n = 30 is exact for all
polynomials up to degree 30. Monomials x^k have most of their "energy"
concentrated in T_k, so for k > 30 the monomial error tracks the Chebyshev
error closely, and both blow up. However, the Chebyshev coefficients of
the monomials with k slightly above 30 are still small enough to keep the
monomial error moderate, while the Chebyshev polynomial T_k concentrates
its degree-k content entirely in a single basis function, yielding a larger
error immediately.

This table motivates the preference for Chebyshev-based quadrature rules.

This script prints a formatted table to stdout (no figure output).

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""

import numpy as np


def newton_cotes_weights(n):
    """
    Compute Newton--Cotes weights for n+1 equispaced nodes on [-1, 1].
    """
    x = np.linspace(-1, 1, n + 1)
    k = np.arange(n + 1, dtype=float)
    rhs = np.where(k % 2 == 0, 2.0 / (k + 1), 0.0)
    V = np.vander(x, increasing=True).T
    w = np.linalg.solve(V, rhs)
    return x, w


def chebyshev_T(k, x):
    """Evaluate Chebyshev polynomial T_k(x)."""
    return np.cos(k * np.arccos(x))


def main():
    """Compute and print the exactness comparison table."""
    n = 30
    x, w = newton_cotes_weights(n)

    # Header
    print(f"\nNewton--Cotes exactness table (n = {n})")
    print("=" * 58)
    print(f"{'k':>4s}   {'|E_n(x^k)|':>16s}   {'|E_n(T_k)|':>16s}   {'Ratio':>10s}")
    print("-" * 58)

    for k in range(26, 40, 2):
        # Monomial x^k
        # Exact integral of x^k over [-1, 1] for even k
        I_mono_exact = 2.0 / (k + 1)
        I_mono_quad = np.dot(w, x**k)
        err_mono = abs(I_mono_quad - I_mono_exact)

        # Chebyshev polynomial T_k
        # Exact integral of T_k(x) over [-1, 1] for even k: 2/(1-k^2)
        I_cheb_exact = 2.0 / (1.0 - k**2)
        T_vals = chebyshev_T(k, x)
        I_cheb_quad = np.dot(w, T_vals)
        err_cheb = abs(I_cheb_quad - I_cheb_exact)

        # Ratio
        if err_mono > 0 and err_cheb > 0:
            ratio = err_cheb / err_mono
        else:
            ratio = float('inf')

        print(f"{k:4d}   {err_mono:16.6e}   {err_cheb:16.6e}   {ratio:10.2e}")

    print("-" * 58)
    print("\nFor k <= n = 30, Newton--Cotes is exact for both bases.")
    print("For k > n, errors in T_k are typically much larger than in x^k,")
    print("because T_k concentrates all its degree-k content in one mode.\n")


if __name__ == '__main__':
    main()
