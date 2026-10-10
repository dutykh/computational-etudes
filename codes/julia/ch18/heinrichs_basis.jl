# heinrichs_basis.jl
#
# Chapter 18: Linear Spectral Eigenproblems --- shared utility.
#
# Heinrichs boundary-adapted bases (1-x^2) T_j and (1-x^2)^2 T_j for
# second- and fourth-order problems on [-1, 1]. See Python docstring
# for full derivation.
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
# Part of "Computational Etudes: A Spectral Approach"

module HeinrichsBasis

using LinearAlgebra

export heinrichs_dirichlet_matrix, heinrichs_clamped_matrix,
       naive_dirichlet_matrix, naive_clamped_operator

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

function heinrichs_dirichlet_matrix(N::Int)
    D, x = cheb_matrix(N)
    D2 = D * D
    Di   = D[2:N, 2:N]
    D2i  = D2[2:N, 2:N]
    xi   = x[2:N]
    S = Diagonal(1.0 .- xi .^ 2)
    X = Diagonal(xi)
    L = -S * D2i + 4.0 * X * Di + 2.0 * Matrix{Float64}(I, N - 1, N - 1)
    A = S \ L
    return A, xi
end

function heinrichs_clamped_matrix(N::Int)
    D, x = cheb_matrix(N)
    D2 = D * D
    D3 = D2 * D
    D4 = D3 * D
    Di  = D[2:N, 2:N]
    D2i = D2[2:N, 2:N]
    D3i = D3[2:N, 2:N]
    D4i = D4[2:N, 2:N]
    xi  = x[2:N]
    I_int = Matrix{Float64}(I, N - 1, N - 1)
    S  = Diagonal(1.0 .- xi .^ 2)
    S2 = Diagonal((1.0 .- xi .^ 2) .^ 2)
    X  = Diagonal(xi)
    X2 = Diagonal(xi .^ 2)
    A = S2 * D4i -
        16.0 * X * S * D3i +
        (48.0 * X2 - 24.0 * I_int) * D2i +
        96.0 * X * Di +
        24.0 * I_int
    M = S2
    return A, M, xi
end

function naive_dirichlet_matrix(N::Int)
    D, x = cheb_matrix(N)
    D2 = D * D
    return -D2[2:N, 2:N], x[2:N]
end

function naive_clamped_operator(N::Int)
    D, x = cheb_matrix(N)
    D4 = (D * D) * (D * D)
    A = copy(D4)
    B = Matrix{Float64}(I, N + 1, N + 1)
    A[1, :] .= 0.0; A[1, 1] = 1.0; B[1, :] .= 0.0
    A[2, :] .= D[1, :];             B[2, :] .= 0.0
    A[N, :] .= D[N+1, :];           B[N, :] .= 0.0
    A[N+1, :] .= 0.0; A[N+1, N+1] = 1.0; B[N+1, :] .= 0.0
    return A, B
end

end # module

using Printf

if abspath(PROGRAM_FILE) == @__FILE__
    using .HeinrichsBasis
    using LinearAlgebra
    println("Condition-number comparison for second-order operator:")
    @printf("%6s %16s %20s\n", "N", "kappa(naive)", "kappa(Heinrichs)")
    for N in (16, 32, 64, 96)
        A_naive, _ = naive_dirichlet_matrix(N)
        A_hein,  _ = heinrichs_dirichlet_matrix(N)
        @printf("%6d %16.3e %20.3e\n", N, cond(A_naive), cond(A_hein))
    end
end
