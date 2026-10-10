# cheb_diff_demo.jl
#
# Demonstrates Chebyshev spectral differentiation using the "Witch of Agnesi"
# test function u(x) = 1/(1 + 4x^2). Shows the function, its derivative,
# and the spectral accuracy for different values of N.
#
# This script generates Figure 6.4 for Chapter 6.
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
using Printf

# Import Chebyshev matrix function
include(joinpath(@__DIR__, "cheb_matrix.jl"))

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
# Configuration
# -----------------------------------------------------------------------------
const N_FINE = 500  # Fine grid for plotting exact solutions

# Book color scheme
const NAVY   = colorant"#142D6E"
const SKY    = colorant"#7896D2"
const CORAL  = colorant"#E74C3C"
const TEAL   = colorant"#16A085"
const PURPLE = colorant"#9B59B6"
const ORANGE = colorant"#E67E22"

# Output paths
SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch07", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Test function: Witch of Agnesi
# -----------------------------------------------------------------------------
u_func(x) = 1.0 / (1.0 + 4.0 * x^2)
u_prime_exact(x) = -8.0 * x / (1.0 + 4.0 * x^2)^2

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Fine grid for exact solution
    x_fine = collect(range(-1, 1, length=N_FINE))

    # Create 2x2 figure
    fig = Figure(size = (700, 550))

    # Main title
    Label(fig[0, :], L"Chebyshev Differentiation of $u(x) = 1/(1+4x^2)$ (Witch of Agnesi)",
          fontsize = 13, font = "CMU Serif Bold")

    # Two values of N to demonstrate
    N_values = [10, 20]

    for (col, Nv) in enumerate(N_values)
        # Construct differentiation matrix
        D, x = cheb_matrix(Nv)

        # Evaluate function at Chebyshev points
        v = u_func.(x)

        # Compute numerical derivative
        w = D * v

        # Exact derivative at Chebyshev points
        w_exact = u_prime_exact.(x)

        # Error
        error_vec = w .- w_exact
        max_error = maximum(abs.(error_vec))

        # Top row: Function and numerical points
        ax_top = Axis(fig[1, col],
                      xlabel = L"x",
                      ylabel = L"u(x)",
                      title = "Function with \$N = $Nv\$ points")
        hidespines!(ax_top, :t, :r)
        xlims!(ax_top, -1.05, 1.05)

        # Grid
        ax_top.ygridvisible = true
        ax_top.ygridcolor = RGBAf(0.533, 0.533, 0.533, 0.2)
        ax_top.ygridstyle = :solid
        ax_top.xgridvisible = true
        ax_top.xgridcolor = RGBAf(0.533, 0.533, 0.533, 0.2)
        ax_top.xgridstyle = :solid

        lines!(ax_top, x_fine, u_func.(x_fine), color = NAVY, linewidth = 1.5,
               label = L"u(x) = 1/(1+4x^2)")
        scatter!(ax_top, x, v, color = CORAL, markersize = 8,
                 strokecolor = :white, strokewidth = 0.5,
                 label = "Chebyshev points")
        axislegend(ax_top, position = :rt, labelsize = 9, framevisible = false)

        # Bottom row: Derivative
        ax_bot = Axis(fig[2, col],
                      xlabel = L"x",
                      ylabel = L"u'(x)",
                      title = "Derivative with \$N = $Nv\$ (max error: $(@sprintf("%.2e", max_error)))")
        hidespines!(ax_bot, :t, :r)
        xlims!(ax_bot, -1.05, 1.05)

        # Grid
        ax_bot.ygridvisible = true
        ax_bot.ygridcolor = RGBAf(0.533, 0.533, 0.533, 0.2)
        ax_bot.ygridstyle = :solid
        ax_bot.xgridvisible = true
        ax_bot.xgridcolor = RGBAf(0.533, 0.533, 0.533, 0.2)
        ax_bot.xgridstyle = :solid

        lines!(ax_bot, x_fine, u_prime_exact.(x_fine), color = NAVY, linewidth = 1.5,
               label = "Exact \$u'(x)\$")
        scatter!(ax_bot, x, w, color = TEAL, markersize = 8,
                 strokecolor = :white, strokewidth = 0.5,
                 label = "Spectral")
        axislegend(ax_bot, position = :rt, labelsize = 9, framevisible = false)
    end

    rowgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 2, 15)

    # Save figure
    save(joinpath(OUTPUT_DIR, "cheb_diff_demo.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "cheb_diff_demo.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "cheb_diff_demo.pdf"))")

    # Print convergence table
    println("\nChebyshev differentiation accuracy for u(x) = 1/(1+4x^2):")
    println("-" ^ 50)
    @printf("%6s %14s %14s\n", "N", "Max Error", "Convergence")
    println("-" ^ 50)

    prev_error = nothing
    for Nv in [4, 8, 10, 12, 16, 20, 24, 32, 40, 48, 64]
        D, x = cheb_matrix(Nv)
        v = u_func.(x)
        w = D * v
        w_exact = u_prime_exact.(x)
        err = maximum(abs.(w .- w_exact))

        if prev_error !== nothing && err > 0
            rate = log10(prev_error / err) / log10(2)
            @printf("%6d %14.6e %14.2f\n", Nv, err, rate)
        else
            @printf("%6d %14.6e %14s\n", Nv, err, "---")
        end

        prev_error = err
    end

    println("-" ^ 50)
    println("Note: The Witch of Agnesi is analytic in [-1,1] with poles at x = +-i/2,")
    println("so exponential convergence is expected until machine precision is reached.")
end

main()
