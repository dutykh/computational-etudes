# collocation_example1.jl
#
# Demonstrates the collocation (pseudospectral) method for solving a boundary
# value problem with a three-coefficient polynomial approximation. This example
# follows the approach in Boyd (2000) but with a different differential equation.
#
# Problem:
#     u''(x) - (4x^2 + 2)u(x) = 0,   -1 <= x <= 1
#     u(-1) = 1,  u(1) = 1
#
# Exact solution:
#     u_exact(x) = exp(x^2 - 1)
#
# Trial function (automatically satisfies BCs):
#     u_2(x) = 1 + (1 - x^2)(a_0 + a_1*x + a_2*x^2)
#
# Collocation points:
#     x = -1/2, 0, 1/2
#
# Result:
#     a_0 = -73/118,  a_1 = 0,  a_2 = -14/59
#
# The figure shows the exact solution vs approximation (left) and the error
# (right), similar to Figure 1.1 in Boyd's "Chebyshev and Fourier Spectral
# Methods".
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

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
const N_POINTS = 500   # Number of points for smooth plotting

# Book color scheme
const NAVY  = colorant"#142D6E"
const SKY   = colorant"#7896D2"
const CORAL = colorant"#E74C3C"

# Output paths
const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch03", "julia")

# Publication-quality theme
set_theme!(Theme(
    fontsize = 11,
    fonts = (regular = "CMU Serif", bold = "CMU Serif Bold", italic = "CMU Serif Italic"),
    Axis = (xlabelsize = 11, ylabelsize = 11, titlesize = 12,
            xticklabelsize = 10, yticklabelsize = 10,
            spinewidth = 0.8, xtickwidth = 0.8, ytickwidth = 0.8),
    Legend = (labelsize = 10,),
))

# -----------------------------------------------------------------------------
# Exact solution
# -----------------------------------------------------------------------------
"""
    u_exact(x)

Exact solution: u(x) = exp(x^2 - 1).

Verification:
    u'(x)  = 2x exp(x^2 - 1)
    u''(x) = (2 + 4x^2) exp(x^2 - 1) = (2 + 4x^2) u(x)

Substituting into ODE:
    u'' - (4x^2 + 2)u = (2 + 4x^2)u - (4x^2 + 2)u = 0

Boundary conditions:
    u(+/-1) = exp(1 - 1) = exp(0) = 1
"""
u_exact(x) = exp(x^2 - 1)

# -----------------------------------------------------------------------------
# Approximate solution from collocation
# -----------------------------------------------------------------------------
"""
    u_approx(x)

Approximate solution using collocation with 3 coefficients.

Trial function: u_2(x) = 1 + (1 - x^2)(a_0 + a_1*x + a_2*x^2)

Collocation at x = -1/2, 0, 1/2 gives:
    a_0 = -73/118
    a_1 = 0
    a_2 = -14/59
"""
function u_approx(x)
    # Coefficients from solving the collocation system
    a0 = -73 / 118
    a1 = 0.0
    a2 = -14 / 59

    # Using the trial function form
    return 1 + (1 - x^2) * (a0 + a1 * x + a2 * x^2)
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Spatial grid
    x = range(-1, 1, length=N_POINTS)

    # Compute solutions
    u_ex = u_exact.(x)
    u_ap = u_approx.(x)
    error = u_ex .- u_ap

    # Collocation points for visualization
    x_coll = [-0.5, 0.0, 0.5]
    u_coll = u_approx.(x_coll)

    # Create two-panel figure (like Boyd's Figure 1.1)
    fig = Figure(size=(720, 280))

    # --- Left panel: Exact vs Approximate ---
    ax1 = Axis(fig[1, 1],
               xlabel = L"x",
               ylabel = L"u(x)",
               title = "Exact vs Approximate",
               titlesize = 11)
    xlims!(ax1, -1, 1)

    # Exact solution (solid line)
    lines!(ax1, collect(x), u_ex, color=NAVY, linewidth=1.5,
           label="Exact")

    # Collocation approximation (markers on the curve, every 25th point)
    idx = 1:25:N_POINTS
    scatter!(ax1, collect(x)[idx], u_ap[idx], color=SKY, markersize=6,
             marker=:circle, label=L"Collocation ($N=3$)")

    # Collocation points (open squares)
    scatter!(ax1, x_coll, u_coll, color=CORAL, markersize=10,
             marker=:rect, strokewidth=1.5, strokecolor=CORAL,
             label="Collocation points")

    axislegend(ax1, position=:ct, framevisible=false, bgcolor=(:white, 0.9))

    # Style: hide top/right spines, add subtle grid
    hidespines!(ax1, :t, :r)
    ax1.yminorgridvisible = false
    ax1.ygridvisible = true
    ax1.ygridstyle = :solid
    ax1.ygridcolor = (:gray, 0.2)

    # --- Right panel: Error ---
    ax2 = Axis(fig[1, 2],
               xlabel = L"x",
               ylabel = L"u_{\mathrm{exact}} - u_{\mathrm{approx}}",
               title = "Error",
               titlesize = 11)
    xlims!(ax2, -1, 1)

    lines!(ax2, collect(x), error, color=NAVY, linewidth=1.5)
    hlines!(ax2, [0.0], color=:gray, linewidth=0.5, linestyle=:dash)

    # Annotate max error
    max_err_idx = argmax(abs.(error))
    max_err = error[max_err_idx]
    x_max = collect(x)[max_err_idx]

    # Arrow annotation for max error
    text!(ax2, 0.3, max_err - 0.005,
          text = "Max error: $(round(max_err, digits=4))",
          fontsize = 9, color = NAVY)
    arrows!(ax2, [0.3], [max_err - 0.004], [x_max - 0.3], [max_err - (max_err - 0.004)],
            color = NAVY, linewidth = 0.8, arrowsize = 8)

    # Style: hide top/right spines, add subtle grid
    hidespines!(ax2, :t, :r)
    ax2.yminorgridvisible = false
    ax2.ygridvisible = true
    ax2.ygridstyle = :solid
    ax2.ygridcolor = (:gray, 0.2)

    # Save figure
    mkpath(OUTPUT_DIR)
    save(joinpath(OUTPUT_DIR, "collocation_example1.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "collocation_example1.png"), fig, px_per_unit=4)

    println("Figure saved to: $(joinpath(OUTPUT_DIR, "collocation_example1.pdf"))")

    # Print error table
    println("\nError comparison at key points:")
    println("-" ^ 60)
    println("       x    u_exact     u_approx        error")
    println("-" ^ 60)
    for xi in [-1.0, -0.5, 0.0, 0.5, 1.0]
        ue = u_exact(xi)
        ua = u_approx(xi)
        err = ue - ua
        @printf("%8.2f %12.5f %12.5f %12.5f\n", xi, ue, ua, err)
    end
    println("-" ^ 60)
    println("Maximum absolute error: $(round(maximum(abs.(error)), digits=5))")

    return fig
end

# Run if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
