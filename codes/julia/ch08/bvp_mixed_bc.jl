# bvp_mixed_bc.jl
#
# Solves a BVP with mixed boundary conditions using Chebyshev spectral collocation:
#
#     u_xx = -cos(πx/2),  x ∈ (-1, 1)
#     u(-1) = 0           (Dirichlet: fixed temperature)
#     u'(1) = 0           (Neumann: insulated boundary)
#
# Exact solution: u(x) = (4/π²)cos(πx/2) + (2/π)(x + 1)
#
# This demonstrates the "boundary bordering" technique: instead of extracting
# the interior system (matrix stripping), we keep the full (N+1)×(N+1) system
# and replace boundary rows with the appropriate conditions.
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
# Problem setup
# -----------------------------------------------------------------------------
rhs_function(x) = -cos(π * x / 2)
exact_solution(x) = (4 / π^2) * cos(π * x / 2) + (2 / π) * (x + 1)

"""
    solve_mixed_bc(N)

Solve u_xx = f(x) with u(-1) = 0, u'(1) = 0 using boundary bordering.

Instead of extracting the interior system, we keep the full (N+1)×(N+1)
system and replace boundary rows with the appropriate conditions.
"""
function solve_mixed_bc(N)
    D2, D, x = cheb_second_derivative_matrix(N)

    # Build full system: L u = rhs
    L = copy(D2)
    rhs = rhs_function.(x)

    # Neumann condition at x[1] = 1: u'(1) = 0
    L[1, :] = D[1, :]
    rhs[1] = 0.0

    # Dirichlet condition at x[N+1] = -1: u(-1) = 0
    L[N+1, :] .= 0.0
    L[N+1, N+1] = 1.0
    rhs[N+1] = 0.0

    # Solve
    u = L \ rhs

    return x, u
end

function main()
    mkpath(OUTPUT_DIR)

    # Create figure
    fig = Figure(size=(950, 400))

    # Solve for N = 16
    N = 16
    x, u_num = solve_mixed_bc(N)
    u_exact_grid = exact_solution.(x)
    error_N16 = maximum(abs.(u_num .- u_exact_grid))

    # Fine grid for exact solution
    x_fine = range(-1, 1, length=500)
    u_exact_fine = exact_solution.(x_fine)

    # Panel 1: Solution
    ax1 = Axis(fig[1, 1],
        xlabel = L"x",
        ylabel = L"u(x)",
        title = @sprintf("Solution (N = %d, max error: %.2e)", N, error_N16),
        titlesize = 11)

    lines!(ax1, collect(x_fine), u_exact_fine, color=NAVY, linewidth=1.5, label="Exact")
    scatter!(ax1, x, u_num, color=TEAL, markersize=8,
             strokecolor=:white, strokewidth=0.5, label="Spectral")

    xlims!(ax1, -1.05, 1.05)
    axislegend(ax1, position=:lt, labelsize=9)
    hidespines!(ax1, :t, :r)

    # Annotate boundary conditions
    u_right = exact_solution(1.0)
    text!(ax1, -0.7, 0.2, text=L"u(-1) = 0", fontsize=10, color=PURPLE)
    lines!(ax1, [0.9, 1.05], [u_right, u_right], color=CORAL, linewidth=2)
    text!(ax1, 0.45, 1.1, text=L"u'(1) = 0", fontsize=10, color=CORAL)

    # Panel 2: Convergence study
    ax2 = Axis(fig[1, 2],
        xlabel = L"N",
        ylabel = "Max Error",
        title = "Convergence",
        titlesize = 11,
        yscale = log10)

    N_values = [4, 6, 8, 10, 12, 16, 20, 24, 28, 32, 36, 40, 48]
    errors = Float64[]

    for N_val in N_values
        x_n, u_n = solve_mixed_bc(N_val)
        u_exact_n = exact_solution.(x_n)
        err = max(maximum(abs.(u_n .- u_exact_n)), 1e-16)
        push!(errors, err)
    end

    scatterlines!(ax2, N_values, errors, color=CORAL, linewidth=1.5, markersize=8,
                  strokecolor=:white, strokewidth=0.5)

    xlims!(ax2, 0, 52)
    ylims!(ax2, 1e-16, 1e1)
    hidespines!(ax2, :t, :r)

    # Machine precision line
    hlines!(ax2, [2.2e-16], color=:gray, linestyle=:dot, linewidth=1)
    text!(ax2, 50, 4e-16, text="Machine precision", fontsize=9, color=:gray,
          align=(:right, :bottom))

    # Main title
    Label(fig[0, :],
        L"Mixed BCs: $u'' = -\cos(\pi x/2)$, $u(-1) = 0$, $u'(1) = 0$",
        fontsize=13)

    # Save figure
    save(joinpath(OUTPUT_DIR, "mixed_bc.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "mixed_bc.png"), fig, px_per_unit=4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "mixed_bc.pdf"))

    # Print convergence table
    println("\nConvergence Table: Mixed BC Problem")
    println("-" ^ 40)
    @printf("%6s %14s\n", "N", "Max Error")
    println("-" ^ 40)
    for (N_val, err) in zip(N_values, errors)
        @printf("%6d %14.2e\n", N_val, err)
    end
    println("-" ^ 40)
end

main()
