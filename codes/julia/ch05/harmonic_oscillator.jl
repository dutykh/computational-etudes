# harmonic_oscillator.jl
#
# Solving the quantum harmonic oscillator eigenvalue problem:
#     -u'' + x^2 u = lambda u,  x in R
#
# using the periodic (Fourier) spectral method on a truncated domain [-L, L].
#
# Exact eigenvalues: lambda_n = 2n + 1 for n = 0, 1, 2, ...
# Exact eigenfunctions: u_n(x) = H_n(x) * exp(-x^2/2) (Hermite functions)
#
# This example demonstrates "spectral accuracy in action":
#     - The eigenfunctions decay like exp(-x^2/2), so truncation is valid
#     - The periodic spectral method achieves machine precision for eigenvalues
#     - Convergence is super-geometric (faster than any exponential)
#
# Author: Dr. Denys Dutykh
#         Mathematics Department
#         Khalifa University of Science and Technology
#         Abu Dhabi, UAE
#
# Part of the book "Computational Etudes: A Spectral Approach"
# https://github.com/dutykh/computational-etudes

using CairoMakie
using LaTeXStrings
using Colors
using LinearAlgebra
using Printf

# -----------------------------------------------------------------------------
# Publication-quality CairoMakie configuration
# -----------------------------------------------------------------------------
set_theme!(Theme(
    fontsize = 11,
    fonts = (regular = "CMU Serif", bold = "CMU Serif Bold", italic = "CMU Serif Italic"),
    Axis = (xlabelsize = 11, ylabelsize = 11, titlesize = 12,
            xticklabelsize = 10, yticklabelsize = 10,
            spinewidth = 0.8, xtickwidth = 0.8, ytickwidth = 0.8),
    Legend = (labelsize = 10,),
))

# -----------------------------------------------------------------------------
# Book color scheme
# -----------------------------------------------------------------------------
const NAVY   = colorant"#142D6E"
const SKY    = colorant"#7896D2"
const CORAL  = colorant"#E74C3C"
const TEAL   = colorant"#16A085"
const PURPLE = colorant"#9B59B6"
const ORANGE = colorant"#E67E22"

# Output paths
const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch05", "julia")

# -----------------------------------------------------------------------------
# Periodic spectral second derivative matrix
# -----------------------------------------------------------------------------
"""
    periodic_d2_matrix(N)

Construct the periodic second derivative matrix for N points.

Uses the formula from Trefethen's "Spectral Methods in MATLAB":
- D2[j,j] = -pi^2/(3h^2) - 1/6
- D2[j,k] = -0.5*(-1)^(j-k) / sin^2(h*(j-k)/2)  for j != k
"""
function periodic_d2_matrix(N)
    h = 2π / N

    # Build first column of Toeplitz matrix
    column = zeros(N)
    column[1] = -π^2 / (3 * h^2) - 1.0 / 6.0

    for k in 2:N
        column[k] = -0.5 * ((-1)^(k - 1)) / sin(h * (k - 1) / 2.0)^2
    end

    # Build symmetric Toeplitz matrix
    D2 = zeros(N, N)
    for i in 1:N
        for j in 1:N
            idx = mod(j - i, N)
            D2[i, j] = column[idx + 1]
        end
    end

    return D2
end

# -----------------------------------------------------------------------------
# Solve harmonic oscillator using periodic spectral method
# -----------------------------------------------------------------------------
"""
    solve_harmonic_oscillator(N, L)

Solve -u'' + x^2 u = lambda u on [-L, L] using periodic spectral method.

# Arguments
- `N::Int`: number of grid points
- `L::Float64`: half-width of domain

# Returns
- `eigenvalues::Vector{Float64}`: computed eigenvalues (sorted)
- `eigenvectors::Matrix{Float64}`: corresponding eigenvectors (columns)
- `x::Vector{Float64}`: grid points
"""
function solve_harmonic_oscillator(N, L)
    # Set up grid
    h = 2π / N
    x_periodic = h * collect(0:N-1)
    x = L * (x_periodic .- π) / π    # Map to [-L, L]

    # Second derivative matrix (with scaling)
    D2 = periodic_d2_matrix(N) * (π / L)^2

    # Potential matrix
    V = Diagonal(x .^ 2)

    # Full operator: -D2 + V
    A = -D2 + V

    # Solve eigenvalue problem (symmetric matrix)
    eig_result = eigen(Symmetric(A))

    return eig_result.values, eig_result.vectors, x
end

# -----------------------------------------------------------------------------
# Exact Hermite functions (normalized)
# -----------------------------------------------------------------------------
"""
    hermite_poly(n, x)

Evaluate the physicist's Hermite polynomial H_n(x) using the three-term recurrence:
    H_0(x) = 1
    H_1(x) = 2x
    H_{n+1}(x) = 2x H_n(x) - 2n H_{n-1}(x)
"""
function hermite_poly(n, x)
    n == 0 && return 1.0
    n == 1 && return 2.0 * x

    H_prev = 1.0
    H_curr = 2.0 * x
    for k in 2:n
        H_next = 2.0 * x * H_curr - 2.0 * (k - 1) * H_prev
        H_prev = H_curr
        H_curr = H_next
    end
    return H_curr
end

"""
    hermite_function(n, x_vals)

Compute the n-th Hermite function (normalized eigenfunction of harmonic oscillator).

psi_n(x) = (1/sqrt(2^n n! sqrt(pi))) * H_n(x) * exp(-x^2/2)

where H_n is the physicist's Hermite polynomial.
"""
function hermite_function(n, x_vals)
    # Normalization constant
    norm_const = 1.0 / sqrt(2.0^n * factorial(n) * sqrt(π))

    return [norm_const * hermite_poly(n, xi) * exp(-xi^2 / 2.0) for xi in x_vals]
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    mkpath(OUTPUT_DIR)

    # Parameters
    L = 8.0    # Domain half-width (large enough for eigenfunction decay)

    # Create figure with subplots
    fig = Figure(size=(1000, 430))

    # -------------------------------------------------------------------------
    # Panel 1: First few eigenfunctions
    # -------------------------------------------------------------------------
    ax1 = Axis(fig[1, 1],
        xlabel = L"x",
        ylabel = L"u_n(x) \;\mathrm{(shifted)}",
        title = @sprintf("Harmonic Oscillator Eigenfunctions (N=%d, L=%.0f)", 64, L),
        titlesize = 11)

    N = 64
    _, eigenvectors, x = solve_harmonic_oscillator(N, L)

    colors = [NAVY, CORAL, TEAL, PURPLE, ORANGE]

    # Fine uniform grid for smooth exact curves
    x_fine = collect(range(-6, 6, length=500))

    for n in 0:4
        # Numerical eigenfunction
        u_num = eigenvectors[:, n + 1]

        # Exact eigenfunction on spectral grid (for sign/scale matching)
        u_exact = hermite_function(n, x)

        # Exact eigenfunction on fine grid (for plotting)
        u_exact_fine = hermite_function(n, x_fine)

        # Fix sign
        if dot(u_num, u_exact) < 0
            u_num = -u_num
        end

        # Scale to match (in case of different normalization)
        scale = maximum(abs.(u_exact)) / maximum(abs.(u_num))
        u_num = u_num * scale

        # Plot exact on fine grid (line) and numerical on spectral grid (markers)
        lines!(ax1, x_fine, u_exact_fine .+ 2n, color=colors[n + 1], linewidth=1.5)
        scatter!(ax1, x[1:4:end], (u_num .+ 2n)[1:4:end], color=colors[n + 1],
                 markersize=5, alpha=0.7)
    end

    # Generic legend entries (compact, 2 items only)
    lines!(ax1, [NaN], [NaN], color=RGBf(0.4, 0.4, 0.4), linewidth=1.5, label="Exact")
    scatter!(ax1, [NaN], [NaN], color=RGBf(0.4, 0.4, 0.4), markersize=5, label="Spectral")

    xlims!(ax1, -6, 6)
    axislegend(ax1, position=:rt, labelsize=9)
    hidespines!(ax1, :t, :r)

    # -------------------------------------------------------------------------
    # Panel 2: Eigenvalue convergence
    # -------------------------------------------------------------------------
    ax2 = Axis(fig[1, 2],
        xlabel = L"Number of grid points $N$",
        ylabel = L"Eigenvalue error $|\lambda_\mathrm{computed} - \lambda_\mathrm{exact}|$",
        title = @sprintf("Eigenvalue Convergence (L=%.0f)", L),
        titlesize = 11,
        yscale = log10)

    N_values = [6, 12, 18, 24, 30, 36, 42, 48]
    exact_eigenvalues = [1, 3, 5, 7]    # First four: 2n+1

    errors = zeros(length(N_values), 4)

    for (i, Nval) in enumerate(N_values)
        eigs, _, _ = solve_harmonic_oscillator(Nval, L)
        for j in 1:4
            errors[i, j] = abs(eigs[j] - exact_eigenvalues[j])
        end
    end

    # Plot convergence
    markers = [:circle, :rect, :utriangle, :diamond]
    for j in 1:4
        scatterlines!(ax2, N_values, errors[:, j],
                      color=colors[j], linewidth=1.5, markersize=8,
                      marker=markers[j],
                      label=latexstring("\\lambda_$(j-1) = $(exact_eigenvalues[j])"))
    end

    # Machine epsilon line
    hlines!(ax2, [2.2e-16], color=:gray, linestyle=:dot, linewidth=1)
    text!(ax2, 50, 4e-16, text="Machine precision", fontsize=9, color=:gray,
          align=(:right, :bottom))

    xlims!(ax2, 0, 52)
    ylims!(ax2, 1e-16, 1e2)
    axislegend(ax2, position=:rt, labelsize=9)
    hidespines!(ax2, :t, :r)

    # Save figure
    save(joinpath(OUTPUT_DIR, "harmonic_oscillator.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "harmonic_oscillator.png"), fig, px_per_unit=4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "harmonic_oscillator.pdf"))

    # Print eigenvalue table
    println("\n" * "=" ^ 70)
    println("Harmonic Oscillator Eigenvalue Convergence (Fourier method)")
    println("Exact eigenvalues: lambda_n = 2n + 1")
    println("=" ^ 70)

    for Nval in [6, 12, 18, 24, 30, 36]
        eigs, _, _ = solve_harmonic_oscillator(Nval, L)
        @printf("\nN = %d\n", Nval)
        for j in 0:3
            @printf("  lambda_%d = %.14f  (error: %.2e)\n", j, eigs[j + 1],
                    abs(eigs[j + 1] - (2j + 1)))
        end
    end

    println("\n" * "=" ^ 70)
end

main()
