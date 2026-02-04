# chebfft.jl
#
# Chapter 10: Spectral PDE Solvers with Chebyshev and Fourier Grids
#
# Chebyshev differentiation via FFT.
# This module provides the core function for computing derivatives on
# Chebyshev grids using the FFT.
#
# Author: Dr. Denys Dutykh
#         Mathematics Department
#         Khalifa University of Science and Technology
#         Abu Dhabi, UAE
#
# Part of the book "Computational Etudes: A Spectral Approach"
# https://github.com/dutykh/computational-etudes

using FFTW
using LinearAlgebra

"""
    chebfft(v)

Compute the first derivative of a function sampled at Chebyshev points
using the FFT.

The algorithm transforms from the x-variable to the theta-variable
(where x = cos(theta)), performs differentiation via FFT in theta-space,
and transforms back.

# Arguments
- `v::AbstractVector`: function values at Chebyshev points
  x_j = cos(j*pi/N), j=0,...,N

# Returns
- `w::Vector{Float64}`: approximate derivative values at the same points
"""
function chebfft(v::AbstractVector)
    v = Float64.(v)
    N = length(v) - 1

    N == 0 && return [0.0]

    x = [cos(π * j / N) for j in 0:N]

    # Extend to even periodic function in theta
    V = [v; v[N:-1:2]]

    # FFT to get Fourier coefficients
    U = real.(fft(V))

    # Wavenumber vector for differentiation
    k = [collect(0:N-1); 0; collect(1-N:-1)]

    # Differentiate in theta: multiply by ik
    W = real.(ifft(1im .* k .* fft(V)))

    # Transform from theta derivative to x derivative
    w = zeros(N + 1)
    w[2:N] = -W[2:N] ./ sqrt.(1.0 .- x[2:N] .^ 2)

    # Boundary values via summation formulas (l'Hopital's rule at x = +/- 1)
    ii = 0:N-1
    w[1]   = sum(ii .^ 2 .* U[ii .+ 1]) / N + 0.5 * N * U[N+1]
    w[N+1] = sum((-1.0) .^ (ii .+ 1) .* ii .^ 2 .* U[ii .+ 1]) / N +
             0.5 * (-1.0)^(N + 1) * N * U[N+1]

    return w
end

# Run test if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    # Also include cheb_matrix for comparison
    include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

    N = 20
    x = [cos(π * j / N) for j in 0:N]

    # Test function: f(x) = exp(x) * cos(4x)
    f = exp.(x) .* cos.(4 .* x)
    f_prime_exact = exp.(x) .* (cos.(4 .* x) .- 4 .* sin.(4 .* x))

    # Compute derivative via FFT
    f_prime_fft = chebfft(f)

    # Compute derivative via matrix
    D, _ = cheb_matrix(N)
    f_prime_mat = D * f

    # Compare errors
    err_fft = maximum(abs.(f_prime_fft .- f_prime_exact))
    err_mat = maximum(abs.(f_prime_mat .- f_prime_exact))

    println("N = $N")
    println("FFT method error:    $(round(err_fft, sigdigits=4))")
    println("Matrix method error: $(round(err_mat, sigdigits=4))")
end
