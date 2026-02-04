# runge_phenomenon.jl
#
# Demonstrates the Runge phenomenon: the failure of polynomial interpolation
# on equispaced grids for certain smooth functions. Named after Carl Runge
# who discovered this counterintuitive behavior in 1901.
#
# The Runge function:
#     f(x) = 1 / (1 + 25x^2)
#
# is smooth and infinitely differentiable, yet polynomial interpolation on
# equispaced nodes diverges as the degree increases. The interpolant oscillates
# wildly near the boundaries x = +/-1.
#
# This script shows interpolation for N = 6, 10, 14 (equispaced nodes)
# demonstrating the progressive deterioration at the boundaries.
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
using LaTeXStrings

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
const N_FINE = 1000  # Number of points for smooth plotting

# Book color scheme
const NAVY  = colorant"#142D6E"
const SKY   = colorant"#7896D2"
const CORAL = colorant"#E74C3C"
const TEAL  = colorant"#16A085"

# Output path
const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch04", "julia")

# -----------------------------------------------------------------------------
# The Runge function
# -----------------------------------------------------------------------------
"""
    runge(x)

The Runge function: f(x) = 1 / (1 + 25x^2).

This function has poles in the complex plane at z = +/-i/5 = +/-0.2i,
which lie very close to the real axis. These nearby singularities
are responsible for the divergence of equispaced interpolation.
"""
runge(x) = 1.0 / (1.0 + 25.0 * x^2)

# -----------------------------------------------------------------------------
# Lagrange interpolation
# -----------------------------------------------------------------------------
"""
    lagrange_interpolate(x_nodes, f_nodes, x_eval)

Evaluate the Lagrange interpolating polynomial at points `x_eval`.

Given N+1 nodes and values, the Lagrange interpolant is:
    p_N(x) = sum_{k=0}^{N} f_k L_k(x)
where L_k(x) are the Lagrange basis polynomials.
"""
function lagrange_interpolate(x_nodes, f_nodes, x_eval)
    n = length(x_nodes)
    p_eval = zeros(length(x_eval))

    for k in 1:n
        L_k = ones(length(x_eval))
        for j in 1:n
            if j != k
                L_k .*= (x_eval .- x_nodes[j]) ./ (x_nodes[k] - x_nodes[j])
            end
        end
        p_eval .+= f_nodes[k] .* L_k
    end

    return p_eval
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Fine grid for plotting the exact function
    x_fine = range(-1, 1, length=N_FINE)
    f_exact = runge.(x_fine)

    # Different polynomial degrees to demonstrate
    N_values = [6, 10, 14]
    colors = [SKY, TEAL, CORAL]

    # Create figure
    fig = Figure(size = (700, 500))
    ax = Axis(fig[1, 1],
              xlabel = L"x",
              ylabel = L"f(x)",
              title = "Runge Phenomenon: Failure of Equispaced Interpolation",
              limits = ((-1, 1), (-0.5, 1.5)))

    # Plot exact function
    lines!(ax, collect(x_fine), f_exact, color = NAVY, linewidth = 2,
           label = "Runge function")

    # Interpolate for each N
    for (N, col) in zip(N_values, colors)
        # Equispaced nodes
        x_nodes = collect(range(-1, 1, length=N+1))
        f_nodes = runge.(x_nodes)

        # Evaluate interpolant
        p_interp = lagrange_interpolate(x_nodes, f_nodes, collect(x_fine))

        # Plot interpolant
        lines!(ax, collect(x_fine), p_interp, color = col, linewidth = 1.2,
               linestyle = :dash, label = latexstring("p_{$N}(x)") * " (equispaced)")

        # Plot nodes
        scatter!(ax, x_nodes, f_nodes, color = col, markersize = 8)
    end

    # Annotation about oscillations
    text!(ax, -0.6, -0.4, text = "Oscillations grow\nnear boundaries",
          fontsize = 9, color = :gray)

    # Formula box
    text!(ax, 0.05, 0.95, text = L"f(x) = \frac{1}{1 + 25x^2}",
          fontsize = 11, space = :relative, align = (:left, :top))

    # Legend
    axislegend(ax, position = :rt, framevisible = false, bgcolor = (:white, 0.9))

    # Save figure
    mkpath(OUTPUT_DIR)
    save(joinpath(OUTPUT_DIR, "runge_phenomenon.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "runge_phenomenon.png"), fig, px_per_unit = 4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "runge_phenomenon.pdf"))

    # Print error statistics
    println("\nMaximum interpolation error at boundaries:")
    println("-" ^ 50)
    for N in N_values
        x_nodes = collect(range(-1, 1, length=N+1))
        f_nodes = runge.(x_nodes)
        p_interp = lagrange_interpolate(x_nodes, f_nodes, collect(x_fine))
        max_error = maximum(abs.(runge.(x_fine) .- p_interp))
        println("N = $(lpad(N, 2)): max |f - p_N| = $(round(max_error, digits=4))")
    end
    println("-" ^ 50)
    println("Note: Error INCREASES with N for equispaced nodes!")

    return fig
end

# Run if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
