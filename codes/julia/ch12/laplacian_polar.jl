#=
laplacian_polar.jl - Discrete polar Laplacian on the unit disk

Constructs the spectral Laplacian operator using Chebyshev discretisation
in the radial direction (doubling trick, r ∈ [-1,1]) and Fourier
discretisation in the angular direction.

Chapter 12: Spectral Methods in Polar Coordinates
Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
Last modified: February 2026
=#

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

using LinearAlgebra

"""
    laplacian_polar(Nr, M)

Build the discrete polar Laplacian on the unit disk.

Uses Chebyshev discretisation in r ∈ [-1,1] (doubling trick, odd Nr)
and Fourier discretisation in θ ∈ [0, 2π). Imposes u = 0 at r = 1.

# Arguments
- `Nr::Int`: number of Chebyshev radial points (must be odd, ≥ 3)
- `M::Int`: number of Fourier angular points (must be even, ≥ 4)

# Returns
- `L::Matrix{Float64}`: discrete Laplacian, size (N2*M) × (N2*M)
- `r_pos::Vector{Float64}`: positive radial grid points (decreasing)
- `theta::Vector{Float64}`: angular grid points
"""
function laplacian_polar(Nr::Int, M::Int)
    @assert isodd(Nr) "Nr must be odd to avoid a grid point at r = 0."
    @assert iseven(M) "M must be even for the Fourier angular grid."

    N2 = (Nr - 1) ÷ 2

    # --- Radial direction: Chebyshev on [-1, 1] ---
    D, r = cheb_matrix(Nr)
    D2full = D * D

    # Block decomposition with reversed negative-r columns
    D1 = D2full[2:N2+1, 2:N2+1]
    D2b = D2full[2:N2+1, Nr:-1:N2+2]
    E1 = D[2:N2+1, 2:N2+1]
    E2 = D[2:N2+1, Nr:-1:N2+2]

    r_pos = r[2:N2+1]
    R = diagm(1.0 ./ r_pos)

    # --- Angular direction: Fourier on [0, 2π) ---
    dt = 2π / M
    theta = [dt * (m - 1) for m in 1:M]
    M2 = M ÷ 2

    # Fourier second derivative Toeplitz matrix
    col = zeros(M)
    col[1] = -π^2 / (3 * dt^2) - 1 / 6
    for j in 1:M-1
        col[j+1] = 0.5 * (-1)^(j + 1) / sin(dt * j / 2)^2
    end
    # Build symmetric Toeplitz matrix
    D2t = [col[abs(i - j) + 1] for i in 1:M, j in 1:M]

    # --- Kronecker product assembly ---
    Z = zeros(M2, M2)
    I_half = Matrix{Float64}(I, M2, M2)
    S = [Z I_half; I_half Z]
    I_M = Matrix{Float64}(I, M, M)

    L = kron(D1 + R * E1, I_M) + kron(D2b + R * E2, S) + kron(R^2, D2t)

    return L, r_pos, theta
end
