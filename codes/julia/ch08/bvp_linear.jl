# bvp_linear.jl
#
# Solves the 1D Poisson equation using Chebyshev spectral collocation:
#
#     u_xx = sin(pi*x) + 2*cos(2*pi*x),  x in (-1, 1),  u(+/-1) = 0
#
# The exact solution is determined by integrating twice and applying
# boundary conditions.
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
"""
    rhs_function(x)

Right-hand side: f(x) = sin(pi*x) + 2*cos(2*pi*x).
"""
rhs_function(x) = sin(π * x) + 2.0 * cos(2.0 * π * x)

"""
    exact_solution(x)

Exact solution of u_xx = sin(pi*x) + 2*cos(2*pi*x) with u(+/-1) = 0.

Integrating sin(pi*x) twice: -sin(pi*x)/pi^2
Integrating 2*cos(2*pi*x) twice: -cos(2*pi*x)/(2*pi^2)

General solution: u(x) = -sin(pi*x)/pi^2 - cos(2*pi*x)/(2*pi^2) + Ax + B

Apply BCs:
u(1) = 0:  -1/(2*pi^2) + A + B = 0
u(-1) = 0: -1/(2*pi^2) - A + B = 0
=> B = 1/(2*pi^2), A = 0

Therefore: u(x) = -sin(pi*x)/pi^2 + (1 - cos(2*pi*x))/(2*pi^2)
"""
exact_solution(x) = -sin(π * x) / π^2 + (1 - cos(2π * x)) / (2π^2)

"""
    solve_poisson_dirichlet(N)

Solve u_xx = f(x) with u(+/-1) = 0 using Chebyshev collocation.

# Arguments
- `N::Int`: number of intervals (N+1 grid points)

# Returns
- `x::Vector{Float64}`: Chebyshev grid points
- `u::Vector{Float64}`: numerical solution at grid points
"""
function solve_poisson_dirichlet(N)
    # Get second derivative matrix and grid
    D2, _, x = cheb_second_derivative_matrix(N)

    # Extract interior points (remove first and last rows/columns)
    # x[1] = 1 (right boundary), x[N+1] = -1 (left boundary)
    D2_int = D2[2:N, 2:N]          # Interior submatrix
    f_int = rhs_function.(x[2:N])   # RHS at interior points

    # Solve the linear system
    u_int = D2_int \ f_int

    # Assemble full solution with boundary conditions
    u = zeros(N + 1)
    u[2:N] = u_int
    u[1] = 0.0      # u(1) = 0
    u[N+1] = 0.0    # u(-1) = 0

    return x, u
end

function main()
    mkpath(OUTPUT_DIR)

    # Fine grid for exact solution
    x_fine = range(-1, 1, length=500)
    u_exact_fine = exact_solution.(x_fine)

    # Solve for N = 16
    N = 16
    x, u_num = solve_poisson_dirichlet(N)
    u_exact_grid = exact_solution.(x)
    error = maximum(abs.(u_num .- u_exact_grid))

    # Create figure
    fig = Figure(size=(950, 400))

    # Panel 1: Solution
    ax1 = Axis(fig[1, 1],
        xlabel = L"x",
        ylabel = L"u(x)",
        title = @sprintf("Solution (N = %d, max error: %.2e)", N, error),
        titlesize = 11,
        xlabelsize = 11, ylabelsize = 11)

    lines!(ax1, collect(x_fine), u_exact_fine, color=NAVY, linewidth=1.5, label="Exact")
    scatter!(ax1, x, u_num, color=TEAL, markersize=8,
             strokecolor=:white, strokewidth=0.5, label="Spectral")

    xlims!(ax1, -1.05, 1.05)
    axislegend(ax1, position=:rt, labelsize=9)
    hidespines!(ax1, :t, :r)

    # Panel 2: Convergence study
    ax2 = Axis(fig[1, 2],
        xlabel = L"N",
        ylabel = "Max Error",
        title = "Convergence",
        titlesize = 11,
        yscale = log10,
        xlabelsize = 11, ylabelsize = 11)

    N_values = [4, 6, 8, 10, 12, 16, 20, 24, 32]
    errors = Float64[]

    for N_val in N_values
        x_n, u_n = solve_poisson_dirichlet(N_val)
        u_exact_n = exact_solution.(x_n)
        err = maximum(abs.(u_n .- u_exact_n))
        push!(errors, err)
    end

    scatterlines!(ax2, N_values, errors, color=CORAL, linewidth=1.5, markersize=8,
                  strokecolor=:white, strokewidth=0.5)

    xlims!(ax2, 0, 35)
    ylims!(ax2, 1e-15, 1)
    hidespines!(ax2, :t, :r)

    # Add machine epsilon line
    hlines!(ax2, [1e-14], color=:gray, linestyle=:dash, linewidth=1)
    text!(ax2, 30, 3e-14, text=L"\approx \mathrm{machine}\; \epsilon", fontsize=9, color=:gray)

    # Main title
    Label(fig[0, :],
        L"1D Poisson Equation: $u_{xx} = \sin(\pi x) + 2\cos(2\pi x)$, $u(\pm 1) = 0$",
        fontsize=13)

    # Save figure
    save(joinpath(OUTPUT_DIR, "poisson_1d.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "poisson_1d.png"), fig, px_per_unit=4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "poisson_1d.pdf"))

    # Print convergence table
    println("\nConvergence Table: 1D Poisson Equation")
    println("-" ^ 40)
    @printf("%6s %14s\n", "N", "Max Error")
    println("-" ^ 40)
    for (N_val, err) in zip(N_values, errors)
        @printf("%6d %14.2e\n", N_val, err)
    end
    println("-" ^ 40)
end

main()
