#!/usr/bin/env julia
#=
quad_exactness_table.jl
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
=#

using LinearAlgebra
using Printf

"""
    newton_cotes_weights(n)

Compute Newton--Cotes weights for n+1 equispaced nodes on [-1, 1].
"""
function newton_cotes_weights(n)
    x = range(-1, 1, length=n+1) |> collect
    V = [x[j]^k for k in 0:n, j in 1:n+1]
    rhs = [(1 - (-1)^(k+1)) / (k + 1) for k in 0:n]
    w = V \ rhs
    return x, w
end

"""
    chebyshev_T(k, x)

Evaluate Chebyshev polynomial T_k(x).
"""
chebyshev_T(k, x) = cos(k * acos(x))

function main()
    n = 30
    x, w = newton_cotes_weights(n)

    # Header
    println()
    @printf("Newton--Cotes exactness table (n = %d)\n", n)
    println("=" ^ 58)
    @printf("%4s   %16s   %16s   %10s\n", "k", "|E_n(x^k)|", "|E_n(T_k)|", "Ratio")
    println("-" ^ 58)

    for k in 26:2:38
        # Monomial x^k
        # Exact integral of x^k over [-1, 1] for even k
        I_mono_exact = 2.0 / (k + 1)
        I_mono_quad = dot(w, x .^ k)
        err_mono = abs(I_mono_quad - I_mono_exact)

        # Chebyshev polynomial T_k
        # Exact integral of T_k(x) over [-1, 1] for even k: 2/(1-k^2)
        I_cheb_exact = 2.0 / (1.0 - k^2)
        T_vals = chebyshev_T.(k, x)
        I_cheb_quad = dot(w, T_vals)
        err_cheb = abs(I_cheb_quad - I_cheb_exact)

        # Ratio
        if err_mono > 0 && err_cheb > 0
            ratio = err_cheb / err_mono
        else
            ratio = Inf
        end

        @printf("%4d   %16.6e   %16.6e   %10.2e\n", k, err_mono, err_cheb, ratio)
    end

    println("-" ^ 58)
    println()
    println("For k <= n = 30, Newton--Cotes is exact for both bases.")
    println("For k > n, errors in T_k are typically much larger than in x^k,")
    println("because T_k concentrates all its degree-k content in one mode.")
    println()
end

main()
