#=
spectral_matrix_periodic.jl - Periodic spectral differentiation matrix

This module constructs the spectral differentiation matrix for periodic
functions on the domain [0, 2π) with N equispaced nodes.

Mathematical Background:
    For periodic functions, the natural interpolation basis consists of
    trigonometric polynomials. The spectral differentiation matrix D is
    derived from differentiating the interpolating trigonometric polynomial.

    The matrix entries are:
        D[i,j] = 0.5 * (-1)^(i-j) * cot((i-j)*π/N)   for i ≠ j
        D[i,i] = 0

    Properties:
    - Skew-symmetric: Dᵀ = -D
    - Toeplitz (and circulant due to periodicity)
    - Dense: all off-diagonal entries are nonzero
    - Exact for trigonometric polynomials of degree ≤ N/2

Reference:
    Trefethen, L.N. (2000). "Spectral Methods in MATLAB", SIAM.

Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
Part of "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
=#

"""
    spectral_diff_periodic(N)

Construct the N×N periodic spectral differentiation matrix for equispaced
points on [0, 2π).

Returns `(D, x)` where `D` is the differentiation matrix and `x` is the
grid vector.
"""
function spectral_diff_periodic(N)
    h = 2π / N
    x = h .* collect(0:N-1)

    D = zeros(N, N)
    for i in 1:N
        for j in 1:N
            if i != j
                d = i - j
                D[i, j] = 0.5 * ((-1)^d) / tan(d * π / N)
            end
        end
    end

    return D, x
end

# Verification when run as a script
using Printf, LinearAlgebra

if abspath(PROGRAM_FILE) == @__FILE__

    println("Verifying periodic spectral differentiation matrix...")
    println("=" ^ 60)

    N = 16
    D, x = spectral_diff_periodic(N)

    # Test 1: Differentiate sin(x)
    u = sin.(x)
    du_exact = cos.(x)
    du_approx = D * u
    err1 = maximum(abs.(du_approx .- du_exact))
    @printf("Test 1: d/dx sin(x) with N=%d, error = %.2e\n", N, err1)

    # Test 2: Differentiate sin(3x)
    u = sin.(3 .* x)
    du_exact = 3 .* cos.(3 .* x)
    du_approx = D * u
    err2 = maximum(abs.(du_approx .- du_exact))
    @printf("Test 2: d/dx sin(3x) with N=%d, error = %.2e\n", N, err2)

    # Test 3: Skew-symmetry
    skew_err = maximum(abs.(D .+ D'))
    @printf("Test 3: Skew-symmetry ||D + Dᵀ||_max = %.2e\n", skew_err)

    # Test 4: Convergence
    println("\nTest 4: Convergence for u(x) = 1/(2+sin(x))")
    @printf("  %4s  %12s\n", "N", "Max Error")
    for Nv in [8, 16, 32, 64]
        D, x = spectral_diff_periodic(Nv)
        u = 1.0 ./ (2.0 .+ sin.(x))
        du_exact = -cos.(x) ./ (2.0 .+ sin.(x)).^2
        du_approx = D * u
        err = maximum(abs.(du_approx .- du_exact))
        @printf("  %4d  %12.4e\n", Nv, err)
    end

    println("\n" * "=" ^ 60)
    all_passed = err1 < 1e-14 && err2 < 1e-12 && skew_err < 1e-14
    println(all_passed ? "All tests PASSED!" : "Some tests FAILED!")
end
