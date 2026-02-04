# equipotential_curves.jl
#
# Visualizes the potential theory explanation for polynomial interpolation
# convergence. Based on Trefethen's Spectral Methods in MATLAB, Program 10.
#
# The potential function associated with a set of interpolation nodes is:
#     phi(z) = -integral_{-1}^{1} mu(x) ln|z - x| dx
#
# where mu(x) is the node density. The equipotential curves phi(z) = const
# determine the region of convergence for polynomial interpolation.
#
# For equispaced nodes:
#     mu(x) = 1/2 (uniform density)
#
# For Chebyshev density:
#     mu(x) = 1 / (pi * sqrt(1 - x^2))
#     phi(z) = log|z + sqrt(z^2 - 1)| - log 2
#     Equipotentials are Bernstein ellipses with foci at +/-1
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
const N_GRID = 500

# Book color scheme
const NAVY  = colorant"#142D6E"
const SKY   = colorant"#7896D2"
const CORAL = colorant"#E74C3C"
const TEAL  = colorant"#16A085"
const LIGHT_GRAY = colorant"#CCCCCC"

# Output path
const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch04", "julia")

# -----------------------------------------------------------------------------
# Potential functions
# -----------------------------------------------------------------------------
"""
    potential_uniform(Z; N_nodes=64)

Potential for uniform (equispaced) node distribution.
Approximates: phi(z) = -integral_{-1}^{1} (1/2) ln|z - x| dx
"""
function potential_uniform(Z; N_nodes=64)
    x_nodes = range(-1, 1, length=N_nodes)
    phi = zeros(size(Z))

    for x_k in x_nodes
        phi .+= log.(abs.(Z .- x_k) .+ 1e-15)
    end

    return phi ./ N_nodes
end

"""
    potential_chebyshev(Z)

Potential for Chebyshev node distribution (closed form):
    phi(z) = log|z + sqrt(z^2 - 1)| - log 2
"""
function potential_chebyshev(Z)
    # Compute z + sqrt(z^2 - 1) carefully
    sqrt_term = sqrt.(Complex.(Z .^ 2 .- 1))
    result = Z .+ sqrt_term

    # Choose the branch that gives |result| > 1
    mask = abs.(result) .< 1
    result[mask] .= Z[mask] .- sqrt_term[mask]

    return real.(log.(abs.(result))) .- log(2)
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Create complex grid
    x = range(-1.5, 1.5, length=N_GRID)
    y = range(-1.0, 1.0, length=N_GRID)
    X = [xi for yi in y, xi in x]
    Y = [yi for yi in y, xi in x]
    Z = X .+ im .* Y

    # Compute potentials
    phi_uniform = potential_uniform(Z)
    phi_chebyshev = potential_chebyshev(Z)

    # Contour levels
    levels_uniform = range(-0.6, 0.6, length=13)
    levels_chebyshev = log.([1.1, 1.2, 1.4, 1.6, 2.0, 2.5, 3.0]) .- log(2)

    # Runge function singularities
    pole_y = 0.2

    # Create two-panel figure
    fig = Figure(size = (1100, 450))

    # Left panel: Uniform (equispaced) potential
    ax1 = Axis(fig[1, 1],
               xlabel = L"\mathrm{Re}(z)", ylabel = L"\mathrm{Im}(z)",
               title = "Equispaced Nodes (Uniform Density)",
               limits = ((-1.5, 1.5), (-1, 1)),
               aspect = DataAspect())

    contour!(ax1, collect(x), collect(y), phi_uniform,
             levels = collect(levels_uniform), color = SKY, linewidth = 0.8)

    # Draw the interval [-1,1]
    lines!(ax1, [-1.0, 1.0], [0.0, 0.0], color = NAVY, linewidth = 2)
    scatter!(ax1, [-1.0, 1.0], [0.0, 0.0], color = NAVY, markersize = 10)

    # Mark singularities
    scatter!(ax1, [0.0, 0.0], [pole_y, -pole_y], color = CORAL, markersize = 15,
             marker = :xcross, strokewidth = 2,
             label = latexstring("Poles at \\pm $(pole_y)i"))

    # Cross-hair lines
    hlines!(ax1, 0, color = LIGHT_GRAY, linewidth = 0.5)
    vlines!(ax1, 0, color = LIGHT_GRAY, linewidth = 0.5)

    axislegend(ax1, position = :rt, framevisible = false, labelsize = 9)

    # Right panel: Chebyshev potential
    ax2 = Axis(fig[1, 2],
               xlabel = L"\mathrm{Re}(z)", ylabel = L"\mathrm{Im}(z)",
               title = "Chebyshev Nodes (Bernstein Ellipses)",
               limits = ((-1.5, 1.5), (-1, 1)),
               aspect = DataAspect())

    contour!(ax2, collect(x), collect(y), phi_chebyshev,
             levels = collect(levels_chebyshev), color = TEAL, linewidth = 0.8)

    # Draw the interval [-1,1]
    lines!(ax2, [-1.0, 1.0], [0.0, 0.0], color = NAVY, linewidth = 2)
    scatter!(ax2, [-1.0, 1.0], [0.0, 0.0], color = NAVY, markersize = 10)

    # Mark singularities
    scatter!(ax2, [0.0, 0.0], [pole_y, -pole_y], color = CORAL, markersize = 15,
             marker = :xcross, strokewidth = 2,
             label = latexstring("Poles at \\pm $(pole_y)i"))

    # Draw the critical Bernstein ellipse
    rho_critical = abs(0.2im + sqrt(complex(-1.04)))
    theta = range(0, 2π, length=200)
    a_crit = (rho_critical + 1/rho_critical) / 2
    b_crit = (rho_critical - 1/rho_critical) / 2
    lines!(ax2, a_crit .* cos.(theta), b_crit .* sin.(theta),
           linestyle = :dash, color = CORAL, linewidth = 1.5,
           label = latexstring("Critical ellipse (\\rho\\approx $(round(rho_critical, digits=2)))"))

    # Cross-hair lines
    hlines!(ax2, 0, color = LIGHT_GRAY, linewidth = 0.5)
    vlines!(ax2, 0, color = LIGHT_GRAY, linewidth = 0.5)

    axislegend(ax2, position = :rt, framevisible = false, labelsize = 9)

    # Save figure
    mkpath(OUTPUT_DIR)
    save(joinpath(OUTPUT_DIR, "equipotential_curves.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "equipotential_curves.png"), fig, px_per_unit = 4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "equipotential_curves.pdf"))

    # Print analysis
    println("\nPotential Theory Analysis:")
    println("-" ^ 60)
    println("Runge function poles at z = +/-$(pole_y)i")
    println("\nChebyshev (Bernstein ellipse) analysis:")
    println("  Critical rho = |$(pole_y)i + sqrt($(pole_y)i^2 - 1)| = $(round(rho_critical, digits=4))")
    println("  Convergence rate: O(rho^{-N}) = O($(round(1/rho_critical, digits=4))^N)")
    println("\nEquispaced analysis:")
    println("  The poles lie inside the smallest equipotential curve")
    println("  that encloses [-1,1], so interpolation diverges.")
    println("-" ^ 60)

    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
