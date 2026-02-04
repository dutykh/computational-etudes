# bvp_variable_coeff.jl
#
# Solves an Airy-type equation with variable coefficients:
#
#     u_xx - (1 + x^2)*u = 1,  x in (-1, 1),  u(+/-1) = 0
#
# This demonstrates that variable coefficient problems require no additional
# complexity with spectral methods: the variable coefficient becomes a
# diagonal matrix multiplication.
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
# Solvers
# -----------------------------------------------------------------------------
"""
    solve_variable_coeff_bvp(N)

Solve u_xx - (1 + x^2)*u = 1 with u(+/-1) = 0.

The problem becomes: [D^2 - diag(1 + x^2)] u = 1
"""
function solve_variable_coeff_bvp(N)
    D2, _, x = cheb_second_derivative_matrix(N)

    # Variable coefficient: a(x) = 1 + x^2
    a = 1.0 .+ x .^ 2

    # Build the operator: L = D^2 - diag(a)
    L = D2 - Diagonal(a)

    # Extract interior system
    L_int = L[2:N, 2:N]
    rhs_int = ones(N - 1)    # f(x) = 1

    # Solve
    u_int = L_int \ rhs_int

    # Assemble with boundary conditions
    u = zeros(N + 1)
    u[2:N] = u_int

    return x, u
end

"""
    solve_constant_coeff_bvp(N)

Solve u_xx - u = 1 with u(+/-1) = 0 for comparison.
"""
function solve_constant_coeff_bvp(N)
    D2, _, x = cheb_second_derivative_matrix(N)

    # Constant coefficient: a = 1
    L = D2 - I(N + 1)

    # Extract interior system
    L_int = L[2:N, 2:N]
    rhs_int = ones(N - 1)

    # Solve
    u_int = L_int \ rhs_int

    # Assemble with boundary conditions
    u = zeros(N + 1)
    u[2:N] = u_int

    return x, u
end

"""
    exact_constant_coeff(x)

Exact solution for u_xx - u = 1 with u(+/-1) = 0.

General solution: u = C1*cosh(x) + C2*sinh(x) - 1
BCs give: C1 = 1/cosh(1), C2 = 0

Therefore: u(x) = cosh(x)/cosh(1) - 1
"""
exact_constant_coeff(x) = cosh(x) / cosh(1) - 1.0

function main()
    mkpath(OUTPUT_DIR)

    N = 32

    # Solve both problems
    x_var, u_var = solve_variable_coeff_bvp(N)
    x_const, u_const = solve_constant_coeff_bvp(N)
    u_const_exact = exact_constant_coeff.(x_const)

    # Fine grid for constant coefficient exact solution
    x_fine = range(-1, 1, length=500)
    u_const_exact_fine = exact_constant_coeff.(collect(x_fine))

    # Create figure
    fig = Figure(size=(950, 400))

    # Panel 1: Compare solutions
    ax1 = Axis(fig[1, 1],
        xlabel = L"x",
        ylabel = L"u(x)",
        title = @sprintf("Comparison of Solutions (N = %d)", N),
        titlesize = 11)

    scatterlines!(ax1, x_var, u_var, color=TEAL, linewidth=1.5, markersize=5,
                  strokecolor=:white, strokewidth=0.3,
                  label=L"Variable: $u_{xx} - (1+x^2)u = 1$")
    scatterlines!(ax1, x_const, u_const, color=CORAL, linewidth=1.5, markersize=5,
                  marker=:rect, strokecolor=:white, strokewidth=0.3,
                  label=L"Constant: $u_{xx} - u = 1$")

    xlims!(ax1, -1.05, 1.05)
    axislegend(ax1, position=:cb, labelsize=9)
    hidespines!(ax1, :t, :r)

    # Annotate the effect of variable coefficient
    mid_idx = div(N, 2) + 1
    text!(ax1, 0.4, -0.4,
        text="Variable coefficient\nreduces amplitude",
        fontsize=9, color=RGBf(0.4, 0.4, 0.4))

    # Panel 2: Verify constant coefficient case
    error = maximum(abs.(u_const .- u_const_exact))

    ax2 = Axis(fig[1, 2],
        xlabel = L"x",
        ylabel = L"u(x)",
        title = @sprintf("Constant Coeff. Verification (error: %.2e)", error),
        titlesize = 11)

    lines!(ax2, collect(x_fine), u_const_exact_fine, color=NAVY, linewidth=1.5, label="Exact")
    scatter!(ax2, x_const, u_const, color=CORAL, markersize=8,
             strokecolor=:white, strokewidth=0.5, label="Spectral")

    xlims!(ax2, -1.05, 1.05)
    axislegend(ax2, position=:cb, labelsize=9)
    hidespines!(ax2, :t, :r)

    # Main title
    Label(fig[0, :],
        L"Variable Coefficient BVP: $u_{xx} - (1+x^2)u = 1$, $u(\pm 1) = 0$",
        fontsize=13)

    # Save figure
    save(joinpath(OUTPUT_DIR, "variable_coeff.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "variable_coeff.png"), fig, px_per_unit=4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "variable_coeff.pdf"))

    # Print solution properties
    println("\nSolution comparison:")
    println("-" ^ 50)
    @printf("Variable coefficient u(0) = %.6f\n", u_var[mid_idx])
    @printf("Constant coefficient u(0) = %.6f\n", u_const[mid_idx])
    @printf("Constant coeff exact u(0) = %.6f\n", exact_constant_coeff(0))
    @printf("Numerical error (const coeff): %.2e\n", error)
    println("-" ^ 50)
end

main()
