# chebyshev_success.jl
#
# Demonstrates successful interpolation of the Runge function using Chebyshev
# nodes. In contrast to the disastrous equispaced interpolation, Chebyshev
# points provide rapid and uniform convergence.
#
# Chebyshev-Gauss-Lobatto points:
#     x_j = cos(j*pi/N),  j = 0, 1, ..., N
#
# These points cluster near the boundaries x = +/-1, counteracting the growth
# of the Lebesgue constant and ensuring convergence for smooth functions.
#
# This script shows interpolation for N = 6, 10, 14 (Chebyshev nodes)
# demonstrating excellent accuracy across the entire interval.
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
const N_FINE = 1000

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
runge(x) = 1.0 / (1.0 + 25.0 * x^2)

# -----------------------------------------------------------------------------
# Chebyshev points
# -----------------------------------------------------------------------------
"""
    chebyshev_nodes(N)

Chebyshev-Gauss-Lobatto points on [-1, 1]:
    x_j = cos(j*pi/N),  j = 0, 1, ..., N
"""
chebyshev_nodes(N) = [cos(j * π / N) for j in 0:N]

# -----------------------------------------------------------------------------
# Lagrange interpolation
# -----------------------------------------------------------------------------
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
    x_fine = collect(range(-1, 1, length=N_FINE))
    f_exact = runge.(x_fine)

    # Different polynomial degrees to demonstrate
    N_values = [6, 10, 14]
    colors = [SKY, TEAL, CORAL]

    # Create figure with two panels
    fig = Figure(size = (1000, 400))

    # Left panel: Function and interpolants
    ax1 = Axis(fig[1, 1],
               xlabel = L"x", ylabel = L"f(x)",
               title = "Chebyshev Interpolation: Success",
               limits = ((-1, 1), (-0.1, 1.1)))

    lines!(ax1, x_fine, f_exact, color = NAVY, linewidth = 2,
           label = "Runge function")

    for (N, col) in zip(N_values, colors)
        x_nodes = chebyshev_nodes(N)
        f_nodes = runge.(x_nodes)
        p_interp = lagrange_interpolate(x_nodes, f_nodes, x_fine)

        lines!(ax1, x_fine, p_interp, color = col, linewidth = 1.2,
               linestyle = :dash, label = latexstring("p_{$N}(x)") * " (Chebyshev)")
        scatter!(ax1, x_nodes, f_nodes, color = col, markersize = 10)
    end

    # Formula
    text!(ax1, 0.05, 0.95, text = L"f(x) = \frac{1}{1 + 25x^2}",
          fontsize = 11, space = :relative, align = (:left, :top))

    axislegend(ax1, position = :rt, framevisible = false, bgcolor = (:white, 0.9))

    # Right panel: Interpolation error
    ax2 = Axis(fig[1, 2],
               xlabel = L"x", ylabel = L"|f(x) - p_N(x)|",
               title = "Interpolation Error (log scale)",
               limits = ((-1, 1), (1e-10, 1)),
               yscale = log10)

    for (N, col) in zip(N_values, colors)
        x_nodes = chebyshev_nodes(N)
        f_nodes = runge.(x_nodes)
        p_interp = lagrange_interpolate(x_nodes, f_nodes, x_fine)
        error = abs.(runge.(x_fine) .- p_interp)

        lines!(ax2, x_fine, error .+ 1e-16, color = col, linewidth = 1.5,
               label = latexstring("N = $N"))
    end

    # Annotation
    text!(ax2, 0.4, 1e-3, text = "Error decreases\nwith N",
          fontsize = 9, color = :gray)

    axislegend(ax2, position = :rt, framevisible = false, bgcolor = (:white, 0.9))

    # Save figure
    mkpath(OUTPUT_DIR)
    save(joinpath(OUTPUT_DIR, "chebyshev_success.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "chebyshev_success.png"), fig, px_per_unit = 4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "chebyshev_success.pdf"))

    # Print error statistics
    println("\nMaximum interpolation error (Chebyshev vs Equispaced):")
    println("-" ^ 60)
    println("$(rpad("N", 4)) $(lpad("Chebyshev", 15)) $(lpad("Equispaced", 15)) $(lpad("Ratio", 12))")
    println("-" ^ 60)
    for N in N_values
        # Chebyshev
        x_cheb = chebyshev_nodes(N)
        f_cheb = runge.(x_cheb)
        p_cheb = lagrange_interpolate(x_cheb, f_cheb, x_fine)
        err_cheb = maximum(abs.(runge.(x_fine) .- p_cheb))

        # Equispaced
        x_equi = collect(range(-1, 1, length=N+1))
        f_equi = runge.(x_equi)
        p_equi = lagrange_interpolate(x_equi, f_equi, x_fine)
        err_equi = maximum(abs.(runge.(x_fine) .- p_equi))

        ratio = err_cheb > 0 ? err_equi / err_cheb : Inf
        println("$(lpad(N, 4)) $(lpad(string(round(err_cheb, sigdigits=2)), 15)) $(lpad(string(round(err_equi, sigdigits=2)), 15)) $(lpad(string(round(ratio, digits=1)), 12))x")
    end
    println("-" ^ 60)
    println("Chebyshev interpolation is dramatically more accurate!")

    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
