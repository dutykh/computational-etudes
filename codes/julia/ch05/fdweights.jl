# fdweights.jl - Finite difference weights for arbitrary nodes
#
# This module implements Fornberg's algorithm for computing finite difference
# weights on arbitrarily spaced grids. The algorithm is numerically stable and
# efficient, requiring O(MN^2) operations where M is the derivative order and
# N is the number of nodes.
#
# Mathematical Background:
#     Given n+1 distinct nodes x_0, x_1, ..., x_n and a point xi where we want
#     to approximate the m-th derivative, Fornberg's algorithm computes weights
#     w_0, w_1, ..., w_n such that:
#
#         f^(m)(xi) ≈ sum_j w_j f(x_j)
#
# References:
#     - Fornberg, B. (1988). "Generation of Finite Difference Formulas on
#       Arbitrarily Spaced Grids", Math. Comp., 51(184), 699-706.
#     - Fornberg, B. (1996). "A Practical Guide to Pseudospectral Methods",
#       Cambridge University Press.
#
# Original MATLAB implementation: Toby Driscoll (2007)
# Julia translation: Dr. Denys Dutykh
#                    Mathematics Department
#                    Khalifa University of Science and Technology
#                    Abu Dhabi, UAE
#
# Part of the book "Computational Etudes: A Spectral Approach"
# https://github.com/dutykh/computational-etudes

"""
    fdweights(xi, x, m)

Compute finite difference weights using Fornberg's algorithm.

# Arguments
- `xi::Real`: evaluation point for the derivative approximation.
- `x::AbstractVector`: node positions (length n+1), can be arbitrarily spaced.
- `m::Int`: order of derivative sought (0 = interpolation, 1 = first derivative, etc.)

# Returns
- `w::Vector{Float64}`: weights for approximating the m-th derivative at xi.
"""
function fdweights(xi::Real, x::AbstractVector, m::Int)
    n = length(x) - 1
    m > n && error("Derivative order m=$m cannot exceed number of intervals n=$n")

    x_shifted = Float64.(x) .- xi

    w = zeros(n + 1)
    for k in 0:n
        w[k+1] = _weight_recursive(x_shifted, m, n, k)
    end
    return w
end

"""
    _weight_recursive(x, m, j, k)

Recursive computation of a single FD weight (Fornberg 1988).
Assumes the evaluation point has been translated to the origin.
"""
function _weight_recursive(x::AbstractVector, m::Int, j::Int, k::Int)
    # Base cases
    (m < 0 || m > j) && return 0.0
    (m == 0 && j == 0) && return 1.0

    if k < j
        # Weight for an "old" node
        c = (x[j+1] * _weight_recursive(x, m, j-1, k) -
             m * _weight_recursive(x, m-1, j-1, k)) / (x[j+1] - x[k+1])
    else
        # Weight for the "new" node
        if j >= 2
            beta = prod(x[j] - x[i+1] for i in 0:j-2) /
                   prod(x[j+1] - x[i+1] for i in 0:j-1)
        elseif j == 1
            beta = 1.0 / (x[2] - x[1])
        else
            beta = 1.0
        end
        c = beta * (m * _weight_recursive(x, m-1, j-1, j-1) -
                    x[j] * _weight_recursive(x, m, j-1, j-1))
    end
    return c
end

"""
    fdweights_all(xi, x, m_max)

Compute FD weights for all derivatives up to order `m_max`.

# Returns
- `W::Matrix{Float64}`: `W[m+1, k+1]` is the weight for node `x[k]` in the
  m-th derivative approximation.
"""
function fdweights_all(xi::Real, x::AbstractVector, m_max::Int)
    n = length(x) - 1
    m_max > n && error("Maximum derivative order m_max=$m_max cannot exceed n=$n")

    x_shifted = Float64.(x) .- xi
    W = zeros(m_max + 1, n + 1)

    for m in 0:m_max
        for k in 0:n
            W[m+1, k+1] = _weight_recursive(x_shifted, m, n, k)
        end
    end
    return W
end

# ==============================================================================
# Verification
# ==============================================================================
function _verify_fdweights()
    println("Verifying fdweights implementation...")
    println("=" ^ 60)

    h = 1.0

    # Test 1: 3-point central difference for first derivative
    x = [-h, 0.0, h]
    w = fdweights(0.0, x, 1)
    expected = [-0.5, 0.0, 0.5] ./ h
    error1 = maximum(abs.(w .- expected))
    println("Test 1: 3-point central difference (1st deriv)")
    println("  Computed: $w")
    println("  Expected: $expected")
    println("  Max error: $(round(error1, sigdigits=2))")

    # Test 2: 5-point central difference for first derivative
    x = [-2h, -h, 0.0, h, 2h]
    w = fdweights(0.0, x, 1)
    expected = [1.0, -8.0, 0.0, 8.0, -1.0] ./ (12h)
    error2 = maximum(abs.(w .- expected))
    println("\nTest 2: 5-point central difference (1st deriv)")
    println("  Max error: $(round(error2, sigdigits=2))")

    # Test 3: 3-point central difference for second derivative
    x = [-h, 0.0, h]
    w = fdweights(0.0, x, 2)
    expected = [1.0, -2.0, 1.0] ./ h^2
    error3 = maximum(abs.(w .- expected))
    println("\nTest 3: 3-point central difference (2nd deriv)")
    println("  Max error: $(round(error3, sigdigits=2))")

    # Test 4: Interpolation (m=0)
    x = [0.0, 1.0]
    w = fdweights(0.5, x, 0)
    expected = [0.5, 0.5]
    error4 = maximum(abs.(w .- expected))
    println("\nTest 4: Linear interpolation at midpoint")
    println("  Max error: $(round(error4, sigdigits=2))")

    # Test 5: Apply to actual function
    x = [-0.1, 0.0, 0.1, 0.2]
    w = fdweights(0.0, x, 1)
    approx = dot(w, exp.(x))
    exact = 1.0
    error5 = abs(approx - exact)
    println("\nTest 5: Differentiate exp(x) at x=0")
    println("  Error: $(round(error5, sigdigits=2))")

    println("\n" * "=" ^ 60)
    all_passed = max(error1, error2, error3, error4) < 1e-14 && error5 < 1e-4
    println(all_passed ? "All tests PASSED!" : "Some tests FAILED!")
    return all_passed
end

# Run verification if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    using LinearAlgebra: dot
    _verify_fdweights()
end
