# harmonic_oscillator.jl
#
# Solving the quantum harmonic oscillator eigenvalue problem:
#     -u'' + x^2 u = lambda u,  x in R
#
# using Chebyshev spectral methods on a truncated domain [-L, L],
# with a comparison to the periodic (Fourier) spectral method.
#
# Exact eigenvalues: lambda_n = 2n + 1 for n = 0, 1, 2, ...
# Exact eigenfunctions: u_n(x) = H_n(x) * exp(-x^2/2) (Hermite functions)
#
# Author: Dr. Denys Dutykh
#         Mathematics Department
#         Khalifa University of Science and Technology
#         Abu Dhabi, UAE
#
# Part of the book "Computational Etudes: A Spectral Approach"
# https://github.com/dutykh/computational-etudes

using CairoMakie
using Colors
using LinearAlgebra
using Printf

# Import Chebyshev matrix from Chapter 7
include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

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
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch08", "julia")

# -----------------------------------------------------------------------------
# Chebyshev method: solve harmonic oscillator
# -----------------------------------------------------------------------------
"""
    solve_chebyshev(N, L)

Solve -u'' + x^2 u = lambda u on [-L, L] using Chebyshev spectral method.

# Arguments
- `N::Int`: number of Chebyshev intervals (N+1 points, N-1 interior)
- `L::Float64`: half-width of domain

# Returns
- `eigenvalues::Vector{Float64}`: computed eigenvalues (sorted)
- `eigenvectors::Matrix{Float64}`: corresponding eigenvectors at interior points
- `x_full::Vector{Float64}`: all grid points (including boundaries)
"""
function solve_chebyshev(N, L)
    # Chebyshev points on [-1, 1]
    D, t = cheb_matrix(N)

    # Map to physical domain [-L, L]
    x_full = L * t

    # Second derivative matrix on physical domain
    D2 = (D * D) / L^2

    # Extract interior system (strip boundary rows and columns)
    D2_int = D2[2:N, 2:N]
    x_int = x_full[2:N]

    # Build operator: A = -D2 + diag(x^2)
    A = -D2_int + Diagonal(x_int .^ 2)

    # Solve eigenvalue problem (A is not symmetric due to Chebyshev D^2)
    eig_result = eigen(Matrix(A))
    idx = sortperm(real.(eig_result.values))
    eigenvalues = real.(eig_result.values[idx])
    eigenvectors = real.(eig_result.vectors[:, idx])

    return eigenvalues, eigenvectors, x_full
end

# -----------------------------------------------------------------------------
# Periodic (Fourier) method: solve harmonic oscillator
# -----------------------------------------------------------------------------
"""
    periodic_d2_matrix(N)

Construct the periodic second derivative matrix for N equispaced points.
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

"""
    solve_fourier(N, L)

Solve -u'' + x^2 u = lambda u on [-L, L] using periodic spectral method.

# Arguments
- `N::Int`: number of equispaced grid points
- `L::Float64`: half-width of domain

# Returns
- `eigenvalues::Vector{Float64}`: computed eigenvalues (sorted)
- `eigenvectors::Matrix{Float64}`: corresponding eigenvectors (columns)
- `x::Vector{Float64}`: grid points
"""
function solve_fourier(N, L)
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

Evaluate the physicist's Hermite polynomial H_n(x) using the three-term recurrence.
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
"""
function hermite_function(n, x_vals)
    norm_const = 1.0 / sqrt(2.0^n * factorial(n) * sqrt(π))
    return [norm_const * hermite_poly(n, xi) * exp(-xi^2 / 2.0) for xi in x_vals]
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    mkpath(OUTPUT_DIR)

    # Parameters
    L = 8.0    # Domain half-width

    # Create figure with subplots
    fig = Figure(size=(1000, 430))

    # -------------------------------------------------------------------------
    # Panel 1: Chebyshev eigenfunctions
    # -------------------------------------------------------------------------
    ax1 = Axis(fig[1, 1],
        xlabel = L"x",
        ylabel = L"u_n(x) \;\mathrm{(shifted)}",
        title = @sprintf("Chebyshev Eigenfunctions (N=%d, L=%.0f)", 64, L),
        titlesize = 11)

    N = 64
    _, eigenvectors, x_full = solve_chebyshev(N, L)

    colors = [NAVY, CORAL, TEAL, PURPLE, ORANGE]

    # Fine uniform grid for smooth exact curves
    x_fine = collect(range(-6, 6, length=500))

    for n in 0:4
        # Numerical eigenfunction (interior only)
        u_int = eigenvectors[:, n + 1]

        # Embed in full grid with boundary zeros
        u_num = zeros(N + 1)
        u_num[2:N] = u_int

        # Exact eigenfunction on Chebyshev grid (for sign/scale matching)
        u_exact = hermite_function(n, x_full)

        # Exact eigenfunction on fine grid (for plotting)
        u_exact_fine = hermite_function(n, x_fine)

        # Fix sign
        if dot(u_num, u_exact) < 0
            u_num = -u_num
        end

        # Scale to match
        scale = maximum(abs.(u_exact)) / maximum(abs.(u_num))
        u_num = u_num * scale

        # Plot exact on fine grid (line) and numerical on Chebyshev grid (markers)
        lines!(ax1, x_fine, u_exact_fine .+ 2n, color=colors[n + 1], linewidth=1.5)
        scatter!(ax1, x_full[1:4:end], (u_num .+ 2n)[1:4:end], color=colors[n + 1],
                 markersize=5, alpha=0.7)
    end

    # Generic legend entries (compact, 2 items only)
    lines!(ax1, [NaN], [NaN], color=RGBf(0.4, 0.4, 0.4), linewidth=1.5, label="Exact")
    scatter!(ax1, [NaN], [NaN], color=RGBf(0.4, 0.4, 0.4), markersize=5, label="Chebyshev")

    xlims!(ax1, -6, 6)
    axislegend(ax1, position=:rt, labelsize=9)
    hidespines!(ax1, :t, :r)

    # -------------------------------------------------------------------------
    # Panel 2: Convergence comparison (Chebyshev vs Fourier)
    # -------------------------------------------------------------------------
    ax2 = Axis(fig[1, 2],
        xlabel = L"Number of grid points $N$",
        ylabel = L"Eigenvalue error",
        title = @sprintf("Fourier vs Chebyshev (L=%.0f)", L),
        titlesize = 11,
        yscale = log10)

    N_values = [12, 18, 24, 30, 36, 42, 48]
    exact_eigenvalues = [1, 3, 5, 7]

    errors_cheb = zeros(length(N_values), 4)
    errors_four = zeros(length(N_values), 4)

    for (i, Nval) in enumerate(N_values)
        eigs_c, _, _ = solve_chebyshev(Nval, L)
        for j in 1:4
            errors_cheb[i, j] = max(abs(eigs_c[j] - exact_eigenvalues[j]), 1e-16)
        end

        eigs_f, _, _ = solve_fourier(Nval, L)
        for j in 1:4
            errors_four[i, j] = max(abs(eigs_f[j] - exact_eigenvalues[j]), 1e-16)
        end
    end

    # Plot lambda_0 and lambda_1 convergence for both methods
    scatterlines!(ax2, N_values, errors_four[:, 1],
                  color=NAVY, linewidth=1.5, markersize=8, marker=:circle,
                  label=L"Fourier $\lambda_0$")
    scatterlines!(ax2, N_values, errors_cheb[:, 1],
                  color=CORAL, linewidth=1.5, markersize=8, marker=:rect,
                  linestyle=:dash,
                  label=L"Chebyshev $\lambda_0$")
    scatterlines!(ax2, N_values, errors_four[:, 2],
                  color=TEAL, linewidth=1.5, markersize=8, marker=:utriangle,
                  label=L"Fourier $\lambda_1$")
    scatterlines!(ax2, N_values, errors_cheb[:, 2],
                  color=PURPLE, linewidth=1.5, markersize=8, marker=:dtriangle,
                  linestyle=:dash,
                  label=L"Chebyshev $\lambda_1$")

    # Machine epsilon line
    hlines!(ax2, [2.2e-16], color=:gray, linestyle=:dot, linewidth=1)
    text!(ax2, 50, 4e-16, text="Machine precision", fontsize=9, color=:gray,
          align=(:right, :bottom))

    xlims!(ax2, 8, 52)
    ylims!(ax2, 1e-16, 1e2)
    axislegend(ax2, position=:rt, labelsize=9)
    hidespines!(ax2, :t, :r)

    # Save figure
    save(joinpath(OUTPUT_DIR, "harmonic_oscillator.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "harmonic_oscillator.png"), fig, px_per_unit=4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "harmonic_oscillator.pdf"))

    # Print comparison table
    println("\n" * "=" ^ 78)
    println("Fourier vs Chebyshev Eigenvalue Comparison (L = 8)")
    println("Exact eigenvalues: lambda_n = 2n + 1")
    println("=" ^ 78)
    @printf("%4s  %16s  %16s  %16s  %16s\n",
            "N", "Fourier err(l0)", "Cheb err(l0)", "Fourier err(l1)", "Cheb err(l1)")
    println("-" ^ 78)

    for (i, Nval) in enumerate(N_values)
        @printf("%4d  %16.2e  %16.2e  %16.2e  %16.2e\n",
                Nval, errors_four[i, 1], errors_cheb[i, 1],
                errors_four[i, 2], errors_cheb[i, 2])
    end

    println("=" ^ 78)

    # Chebyshev eigenvalue table (N=36)
    println("\nChebyshev eigenvalues at N = 36, L = 8:")
    eigs_c, _, _ = solve_chebyshev(36, L)
    for j in 0:3
        @printf("  lambda_%d = %.14f  (error: %.2e)\n", j, eigs_c[j + 1],
                abs(eigs_c[j + 1] - (2j + 1)))
    end

    println("\nFourier eigenvalues at N = 36, L = 8:")
    eigs_f, _, _ = solve_fourier(36, L)
    for j in 0:3
        @printf("  lambda_%d = %.14f  (error: %.2e)\n", j, eigs_f[j + 1],
                abs(eigs_f[j + 1] - (2j + 1)))
    end
end

main()
