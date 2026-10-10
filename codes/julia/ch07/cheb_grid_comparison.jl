# cheb_grid_comparison.jl
#
# Visualizes the comparison between equispaced and Chebyshev-Gauss-Lobatto grids.
# Shows why Chebyshev points cluster near the boundaries and how this relates
# to the projection from a circle.
#
# This script generates Figure 6.1 for Chapter 6.
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
const N = 16  # Number of grid points

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
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Grid points
    x_equi = range(-1, 1, length=N + 1) |> collect
    x_cheb = [cos(π * j / N) for j in 0:N]

    # Create figure with three rows
    fig = Figure(size = (700, 600))

    # Panel 1: Equispaced grid (top)
    ax1 = Axis(fig[1, 1:2],
               xlabel = L"x",
               title = "Equispaced Grid (\$N = $N\$ intervals)",
               yticksvisible = false, yticklabelsvisible = false,
               xticks = ([-1.0, -0.5, 0.0, 0.5, 1.0],
                         [L"-1", L"-0.5", L"0", L"0.5", L"1"]))
    xlims!(ax1, -1.1, 1.1)
    ylims!(ax1, -0.2, 0.2)
    hidespines!(ax1, :t, :r, :l)

    # Horizontal baseline
    hlines!(ax1, 0, color = RGBf(0.533, 0.533, 0.533), linewidth = 0.5)
    # Equispaced points
    scatter!(ax1, x_equi, zeros(N + 1), color = CORAL, markersize = 10,
             strokecolor = :white, strokewidth = 1)

    # Panel 2: Chebyshev grid (middle)
    ax2 = Axis(fig[2, 1:2],
               xlabel = L"x",
               title = "Chebyshev-Gauss-Lobatto Grid (\$N = $N\$ intervals)",
               yticksvisible = false, yticklabelsvisible = false,
               xticks = ([-1.0, -0.5, 0.0, 0.5, 1.0],
                         [L"-1", L"-0.5", L"0", L"0.5", L"1"]))
    xlims!(ax2, -1.1, 1.1)
    ylims!(ax2, -0.2, 0.2)
    hidespines!(ax2, :t, :r, :l)

    # Horizontal baseline
    hlines!(ax2, 0, color = RGBf(0.533, 0.533, 0.533), linewidth = 0.5)
    # Chebyshev points
    scatter!(ax2, x_cheb, zeros(N + 1), color = TEAL, markersize = 10,
             strokecolor = :white, strokewidth = 1)

    # Panel 3: Circle projection (bottom-left)
    ax3 = Axis(fig[3, 1],
               xlabel = L"x = \cos\theta",
               ylabel = L"\sin\theta",
               title = L"Circle Projection: $x_j = \cos(j\pi/N)$",
               aspect = DataAspect())
    xlims!(ax3, -1.3, 1.3)
    ylims!(ax3, -0.2, 1.3)
    hidespines!(ax3, :t, :r)

    # Draw semicircle
    theta_fine = range(0, π, length = 200)
    lines!(ax3, cos.(theta_fine), sin.(theta_fine), color = SKY, linewidth = 2)

    # Horizontal and vertical axes
    hlines!(ax3, 0, color = RGBf(0.267, 0.267, 0.267), linewidth = 0.8)
    vlines!(ax3, 0, color = RGBf(0.267, 0.267, 0.267), linewidth = 0.5, alpha = 0.5)

    # Draw projection lines and points for Chebyshev nodes
    theta_cheb = [π * j / N for j in 0:N]
    for t in theta_cheb
        xc = cos(t)
        yc = sin(t)
        # Point on circle
        scatter!(ax3, [xc], [yc], color = TEAL, markersize = 8)
        # Projection line
        lines!(ax3, [xc, xc], [yc, 0.0], color = RGBf(0.533, 0.533, 0.533),
               linewidth = 0.5, linestyle = :dash, alpha = 0.6)
        # Point on x-axis
        scatter!(ax3, [xc], [0.0], color = NAVY, markersize = 7, marker = :rect)
    end

    # Panel 4: Boundary spacing zoom (bottom-right)
    ax4 = Axis(fig[3, 2],
               xlabel = "Interval index (from left boundary)",
               ylabel = L"Spacing $\Delta x$",
               title = "Grid Spacing Near Boundary")
    hidespines!(ax4, :t, :r)

    # Compute spacing for both grids
    spacing_equi = diff(sort(x_equi))
    spacing_cheb = diff(sort(x_cheb))

    # Show first 8 spacings
    n_show = 8
    x_pos = 1:n_show
    width = 0.35

    # Bar plots
    barplot!(ax4, collect(x_pos) .- width/2, spacing_equi[1:n_show],
             width = width, color = CORAL, label = "Equispaced")
    barplot!(ax4, collect(x_pos) .+ width/2, spacing_cheb[1:n_show],
             width = width, color = TEAL, label = "Chebyshev")

    ax4.xticks = (collect(x_pos), [string(i) for i in 1:n_show])
    axislegend(ax4, position = :rt, labelsize = 9, framevisible = false)

    # Annotation about boundary spacing
    cheb_boundary_spacing = spacing_cheb[1]
    text!(ax4, 3.0, 0.14,
          text = "Chebyshev: \$\\Delta x_1 \\approx\$ $(round(cheb_boundary_spacing, digits=4))\n\$\\sim \\pi^2/2N^2\$",
          fontsize = 9, color = TEAL, align = (:center, :center))

    # Add arrows from annotation to bar
    arrows!(ax4, [2.5], [0.12], [-1.8], [-0.1],
            color = TEAL, linewidth = 0.8)

    # Main title
    Label(fig[0, :], "Equispaced vs. Chebyshev-Gauss-Lobatto Grids",
          fontsize = 13, font = "CMU Serif Bold")

    # Adjust row gaps
    rowgap!(fig.layout, 1, 15)
    rowgap!(fig.layout, 2, 15)
    rowgap!(fig.layout, 3, 20)

    # Save figure
    save(joinpath(OUTPUT_DIR, "grid_comparison.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "grid_comparison.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "grid_comparison.pdf"))")

    # Print spacing information
    println("\nGrid spacing analysis (N = $N):")
    println("-" ^ 50)
    println("$(rpad("Spacing", 12)) $(lpad("Equispaced", 12)) $(lpad("Chebyshev", 12)) $(lpad("Ratio", 12))")
    println("-" ^ 50)
    println("$(rpad("Minimum", 12)) $(lpad(round(minimum(spacing_equi), digits=6), 12)) $(lpad(round(minimum(spacing_cheb), digits=6), 12)) $(lpad(round(minimum(spacing_equi)/minimum(spacing_cheb), digits=2), 12))")
    println("$(rpad("Maximum", 12)) $(lpad(round(maximum(spacing_equi), digits=6), 12)) $(lpad(round(maximum(spacing_cheb), digits=6), 12)) $(lpad(round(maximum(spacing_equi)/maximum(spacing_cheb), digits=2), 12))")
    println("$(rpad("At boundary", 12)) $(lpad(round(spacing_equi[1], digits=6), 12)) $(lpad(round(spacing_cheb[1], digits=6), 12)) $(lpad(round(spacing_equi[1]/spacing_cheb[1], digits=2), 12))")
    println("-" ^ 50)
    println("Theoretical boundary spacing: pi^2/(2*N^2) = $(round(π^2/(2*N^2), digits=6))")
end

main()
