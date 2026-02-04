# bvp_eigenvalue.jl
#
# Solves the eigenvalue problem:
#
#     u_xx = lambda * u,  x in (-1, 1),  u(+/-1) = 0
#
# Exact eigenvalues: lambda_n = -(n*pi/2)^2 for n = 1, 2, 3, ...
# Exact eigenfunctions: u_n(x) = sin(n*pi*(x+1)/2)
#
# This demonstrates the resolution limits of spectral methods: higher modes
# require more points per wavelength (ppw) for accuracy.
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

# Import Chebyshev matrix function from Chapter 7
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
const PURPLE = colorant"#8E44AD"
const ORANGE = colorant"#E67E22"

# Output paths
const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch08", "julia")

# -----------------------------------------------------------------------------
# Exact solutions
# -----------------------------------------------------------------------------
"""Exact eigenvalue lambda_n = -(n*pi/2)^2."""
exact_eigenvalue(n) = -(n * π / 2)^2

"""Exact eigenfunction u_n(x) = sin(n*pi*(x+1)/2)."""
exact_eigenfunction(x, n) = sin(n * π * (x + 1) / 2)

# -----------------------------------------------------------------------------
# Solver
# -----------------------------------------------------------------------------
"""
    solve_eigenvalue_problem(N)

Solve the eigenvalue problem u_xx = lambda*u with u(+/-1) = 0.

# Returns
- `eigenvalues::Vector{Float64}`: computed eigenvalues (sorted, least negative first)
- `eigenvectors::Matrix{Float64}`: corresponding eigenvectors (columns)
- `x::Vector{Float64}`: grid points
"""
function solve_eigenvalue_problem(N)
    D2, _, x = cheb_second_derivative_matrix(N)

    # Extract interior system
    D2_int = D2[2:N, 2:N]

    # Solve eigenvalue problem
    eig_result = eigen(D2_int)
    eigenvalues = real.(eig_result.values)
    eigenvectors = real.(eig_result.vectors)

    # Sort by eigenvalue (largest = least negative first)
    idx = sortperm(eigenvalues, rev=true)
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    return eigenvalues, eigenvectors, x
end

function main()
    mkpath(OUTPUT_DIR)

    N = 36    # Grid size

    # Solve eigenvalue problem
    eigenvalues, eigenvectors, x = solve_eigenvalue_problem(N)

    # Modes to display (1-indexed)
    modes = [5, 10, 15, 20, 25, 30]
    colors = [CORAL, TEAL, PURPLE, ORANGE, SKY, NAVY]

    # Fine grid for exact solutions
    x_fine = range(-1, 1, length=500)

    # Create figure with 6 panels (2x3)
    fig = Figure(size=(1000, 600))

    for (idx, (mode, color)) in enumerate(zip(modes, colors))
        row = div(idx - 1, 3) + 1
        col = mod(idx - 1, 3) + 1

        # Numerical eigenvalue and eigenfunction
        lam_num = eigenvalues[mode]
        u_num = eigenvectors[:, mode]

        # Exact values
        lam_exact = exact_eigenvalue(mode)
        u_exact_fine = exact_eigenfunction.(collect(x_fine), mode)
        u_exact_grid = exact_eigenfunction.(x[2:N], mode)

        # Normalize numerical eigenfunction to match sign of exact
        if sign(u_num[1]) != sign(u_exact_grid[1])
            u_num = -u_num
        end

        # Normalize amplitude
        u_num = u_num / maximum(abs.(u_num))
        u_exact_fine = u_exact_fine / maximum(abs.(u_exact_fine))

        # Eigenvalue error
        lam_error = abs(lam_num - lam_exact) / abs(lam_exact)

        # Points per wavelength
        wavelength = 4.0 / mode    # wavelength = 2L/n = 4/n for L=2
        avg_spacing = 2.0 / N
        ppw = wavelength / avg_spacing

        ax = Axis(fig[row, col],
            xlabel = L"x",
            ylabel = L"u(x)",
            title = @sprintf("Mode %d: ppw = %.1f, λ error = %.1e", mode, ppw, lam_error),
            titlesize = 10)

        lines!(ax, collect(x_fine), u_exact_fine, color=NAVY, linewidth=1.5,
               alpha=0.7, label="Exact")
        scatter!(ax, x[2:N], u_num, color=color, markersize=5,
                 strokecolor=:white, strokewidth=0.3, label="Spectral")

        xlims!(ax, -1.05, 1.05)
        ylims!(ax, -1.3, 1.3)
        hidespines!(ax, :t, :r)

        if idx == 1
            axislegend(ax, position=:rt, labelsize=8)
        end

        # Add reference line
        hlines!(ax, [0], color=RGBf(0.53, 0.53, 0.53), linewidth=0.5)
    end

    # Main title
    Label(fig[0, :],
        @sprintf("Eigenvalue Problem: u_xx = λu, u(±1) = 0 (N = %d)", N),
        fontsize=13)

    # Save figure
    save(joinpath(OUTPUT_DIR, "eigenvalue_problem.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "eigenvalue_problem.png"), fig, px_per_unit=4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "eigenvalue_problem.pdf"))

    # Print eigenvalue table
    @printf("\nEigenvalue accuracy for N = %d:\n", N)
    println("-" ^ 70)
    @printf("%6s %14s %14s %12s %8s\n", "Mode", "λ (exact)", "λ (numerical)", "Rel. Error", "ppw")
    println("-" ^ 70)

    for mode in 1:min(N - 1, 35)
        lam_exact = exact_eigenvalue(mode)
        lam_num = eigenvalues[mode]
        rel_error = abs(lam_num - lam_exact) / abs(lam_exact)
        ppw = (4.0 / mode) / (2.0 / N)

        @printf("%6d %14.4f %14.4f %12.2e %8.2f\n", mode, lam_exact, lam_num, rel_error, ppw)
    end

    println("-" ^ 70)
    println("\nNote: Accuracy degrades when ppw < pi (~ 3.14).")
    println("The rule of thumb is: need >= pi points per wavelength for spectral accuracy.")
end

main()
