# rational_chebyshev.jl
#
# Chapter 18: Linear Spectral Eigenproblems --- shared utility.
#
# Rational Chebyshev basis on (-inf, +inf) via the algebraic map
#
#     x = ell * t / sqrt(1 - t^2),
#
# with t on the Chebyshev-Gauss-Lobatto grid in [-1, 1]. The endpoints
# t = +-1 map to x = +-inf, and the interior nodes provide (N-1) real-line
# collocation points.
#
# Chain-rule factors:
#     dt/dx         = (1 - t^2)^(3/2) / ell
#     d^2 t / d x^2 = -3 t (1 - t^2)^2 / ell^2
#
# First and second derivative matrices in x follow from the chain rule.
#
# The map parameter is named `ell` (Greek script ell) to avoid collision
# with the truncation half-width `L` used in chapter 20.
#
# Used by Etude 18.4 (harmonic oscillator) and Etude 18.6 (bound states
# plus continuum).
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
# Part of "Computational Etudes: A Spectral Approach"
# https://github.com/dutykh/computational-etudes

module RationalChebyshev

using LinearAlgebra

export rational_chebyshev_grid, rational_chebyshev_derivative_matrices

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

"""
    rational_chebyshev_grid(N, ell = 4.0)

Return `(x_int, t, D_t)` where `x_int` is the vector of (N-1) interior
collocation points on the real line, `t` is the full Chebyshev-Lobatto
grid on [-1, 1] with N+1 nodes, and `D_t` is the full (N+1) x (N+1)
Chebyshev differentiation matrix.

`ell` is the textbook's map parameter (Greek script ell), distinct from
the chapter-20 truncation half-width `L`.
"""
function rational_chebyshev_grid(N::Int, ell::Real = 4.0)
    D_t, t = cheb_matrix(N)
    t_int = t[2:N]
    x_int = ell .* t_int ./ sqrt.(1.0 .- t_int .^ 2)
    return x_int, t, D_t
end

"""
    rational_chebyshev_derivative_matrices(N, ell = 4.0)

Return `(D1_x, D2_x, x_int)` where `D1_x`, `D2_x` are the first and
second derivative matrices in x on the interior rational Chebyshev
grid, built by chain-rule composition of the Chebyshev t-derivatives.
"""
function rational_chebyshev_derivative_matrices(N::Int, ell::Real = 4.0)
    x_int, t, D_t = rational_chebyshev_grid(N, ell)

    t_int = t[2:N]
    one_minus_t2 = 1.0 .- t_int .^ 2
    dt_dx   = (one_minus_t2 .^ 1.5) ./ ell
    d2t_dx2 = -3.0 .* t_int .* (one_minus_t2 .^ 2) ./ (ell ^ 2)

    D_t_int  = D_t[2:N, 2:N]
    D_t2_int = (D_t * D_t)[2:N, 2:N]

    D1_x = Diagonal(dt_dx)       * D_t_int
    D2_x = Diagonal(dt_dx .^ 2)  * D_t2_int + Diagonal(d2t_dx2) * D_t_int

    return D1_x, D2_x, x_int
end

end # module

using Printf

# Self-test when run as a script
if abspath(PROGRAM_FILE) == @__FILE__
    using .RationalChebyshev: rational_chebyshev_derivative_matrices
    for N in (16, 32, 64)
        _, D2, x = rational_chebyshev_derivative_matrices(N, 4.0)
        u = exp.(-x .^ 2 ./ 2.0)
        u_xx_num = D2 * u
        u_xx_ex  = (x .^ 2 .- 1.0) .* u
        err = maximum(abs.(u_xx_num .- u_xx_ex))
        @printf("N = %3d  ell = 4   max error in u_xx (Gaussian) = %.3e\n", N, err)
    end
end
