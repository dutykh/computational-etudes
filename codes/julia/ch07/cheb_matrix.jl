# cheb_matrix.jl
#
# Core Chebyshev differentiation matrix construction for spectral methods
# on bounded (non-periodic) domains.
#
# The Chebyshev-Gauss-Lobatto points are:
#     x_j = cos(j*pi/N),  j = 0, 1, ..., N
#
# These cluster near the boundaries, avoiding the Runge phenomenon and
# enabling spectral accuracy for smooth functions on [-1, 1].
#
# Author: Dr. Denys Dutykh
#         Mathematics Department
#         Khalifa University of Science and Technology
#         Abu Dhabi, UAE
#
# Part of the book "Computational Etudes: A Spectral Approach"
# https://github.com/dutykh/computational-etudes

using LinearAlgebra

"""
    cheb_matrix(N)

Construct the Chebyshev differentiation matrix and grid points.

Uses the explicit formulas from Trefethen's "Spectral Methods in MATLAB"
with the "negative sum trick" for diagonal entries to ensure numerical
stability (derivative of a constant is exactly zero).

# Arguments
- `N::Int`: number of intervals (matrix size is (N+1) x (N+1))

# Returns
- `D::Matrix{Float64}`: Chebyshev differentiation matrix
- `x::Vector{Float64}`: Chebyshev-Gauss-Lobatto points on [-1, 1]
"""
function cheb_matrix(N::Int)
    N == 0 && return (zeros(1, 1), [1.0])

    # Chebyshev-Gauss-Lobatto points
    x = [cos(π * j / N) for j in 0:N]

    # Coefficient vector c: c_0 = c_N = 2, others = 1
    c = ones(N + 1)
    c[1] = 2.0
    c[N+1] = 2.0

    # Construct the differentiation matrix using outer products
    X = repeat(x', N + 1, 1)    # X[i,j] = x[j]
    dX = X .- X'                 # dX[i,j] = x[j] - x[i]

    # Off-diagonal entries: D[i,j] = c[i]/c[j] * (-1)^(i+j) / (x[i] - x[j])
    C = c * (1.0 ./ c')          # C[i,j] = c[i] / c[j]

    # Sign matrix: (-1)^(i+j), using 0-based indices
    idx = 0:N
    sign_mat = ((-1.0) .^ idx) * ((-1.0) .^ idx)'

    # Compute off-diagonal entries
    D = C .* sign_mat ./ (-dX .+ I(N + 1))  # add I to avoid 0/0 on diagonal
    D = D .- Diagonal(diag(D))               # zero out diagonal

    # Apply negative sum trick for diagonal entries
    for i in 1:N+1
        D[i, i] = -sum(D[i, :])
    end

    return D, x
end

"""
    cheb_second_derivative_matrix(N)

Construct the second derivative matrix D^2 via matrix squaring.

# Returns
- `D2::Matrix{Float64}`: second derivative matrix
- `D::Matrix{Float64}`: first derivative matrix
- `x::Vector{Float64}`: Chebyshev grid points
"""
function cheb_second_derivative_matrix(N::Int)
    D, x = cheb_matrix(N)
    D2 = D * D
    return D2, D, x
end

"""
    verify_small_matrices()

Verify the matrix construction against known small cases.
"""
function verify_small_matrices()
    all_pass = true

    # D_1: 2x2 matrix
    D1, _ = cheb_matrix(1)
    D1_exact = [0.5 -0.5; 0.5 -0.5]
    if !isapprox(D1, D1_exact, atol=1e-14)
        println("D_1 verification FAILED")
        all_pass = false
    else
        println("D_1 verification PASSED")
    end

    # D_2: 3x3 matrix
    D2, _ = cheb_matrix(2)
    D2_exact = [1.5 -2.0 0.5; 0.5 0.0 -0.5; -0.5 2.0 -1.5]
    if !isapprox(D2, D2_exact, atol=1e-14)
        println("D_2 verification FAILED")
        all_pass = false
    else
        println("D_2 verification PASSED")
    end

    # Verify negative sum trick: D * ones = zeros
    for N in [4, 8, 16, 32]
        D, _ = cheb_matrix(N)
        result = D * ones(N + 1)
        max_error = maximum(abs.(result))
        status = max_error < 1e-12 ? "PASSED" : "FAILED"
        println("Negative sum trick $status for N=$N: max error = $(round(max_error, sigdigits=2))")
        max_error >= 1e-12 && (all_pass = false)
    end

    return all_pass
end

"""
    demo_differentiation()

Demonstrate spectral differentiation accuracy using the Witch of Agnesi.
"""
function demo_differentiation()
    u(x) = 1.0 / (1.0 + 4.0 * x^2)
    u_exact(x) = -8.0 * x / (1.0 + 4.0 * x^2)^2

    println("\nDifferentiation accuracy for u(x) = 1/(1 + 4x^2):")
    println("-" ^ 50)
    println("     N     Max Error")
    println("-" ^ 50)

    for N in [4, 8, 16, 32, 64]
        D, x = cheb_matrix(N)
        v = u.(x)
        w = D * v
        w_ex = u_exact.(x)
        err = maximum(abs.(w .- w_ex))
        println("  $(lpad(N, 4))   $(err)")
    end
    println("-" ^ 50)
end

# Run verification if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    println("=" ^ 60)
    println("Chebyshev Differentiation Matrix Verification")
    println("=" ^ 60)
    println("\n1. Verifying small matrices (D_1, D_2)...")
    verify_small_matrices()
    println("\n2. Demonstrating differentiation accuracy...")
    demo_differentiation()
    println("\n" * "=" ^ 60)
    println("Verification complete!")
    println("=" ^ 60)
end
