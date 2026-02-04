# aliasing_demo.jl
#
# Chapter 9: Physical and Fourier Space on Grids
#
# Demonstrates aliasing: sin(pi*x) and sin(9*pi*x) coincide on grid h=1/4.
# This is Etude 1 from the chapter.
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

# Book color scheme
const NAVY   = colorant"#142D6E"
const SKY    = colorant"#7896D2"
const CORAL  = colorant"#E74C3C"
const TEAL   = colorant"#16A085"
const PURPLE = colorant"#9B59B6"
const ORANGE = colorant"#E67E22"

# Output path
const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch09", "julia")

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Grid spacing
    h = 0.25

    # Continuous plotting grid
    x_cont = range(0, 2, length=500)

    # Discrete grid points
    x_grid = collect(0:h:2)

    # Two functions with different frequencies
    f_cont = sin.(π .* x_cont)       # Wavenumber k = pi
    g_cont = sin.(9π .* x_cont)      # Wavenumber k = 9*pi

    # Samples at grid points
    f_samples = sin.(π .* x_grid)
    g_samples = sin.(9π .* x_grid)

    # Verify aliasing
    max_diff = maximum(abs.(f_samples .- g_samples))
    println("Maximum difference at grid points: $(round(max_diff, sigdigits=2))")
    println("Aliasing relation: 9pi = pi + 2pi*4 = pi + 2pi/h")

    # Create figure
    fig = Figure(size = (600, 320))
    ax = Axis(fig[1, 1],
              xlabel = L"x",
              ylabel = "Value",
              title  = L"Aliasing on grid $h = 1/4$: $\sin(\pi x)$ and $\sin(9\pi x)$ coincide",
              limits = ((-0.05, 2.05), (-1.3, 1.3)))

    # Plot continuous functions
    lines!(ax, collect(x_cont), f_cont, color = NAVY, linewidth = 1.5,
           label = L"\sin(\pi x)")
    lines!(ax, collect(x_cont), g_cont, color = CORAL, linewidth = 1.5,
           linestyle = :dash, label = L"\sin(9\pi x)")

    # Plot grid samples
    scatter!(ax, x_grid, f_samples, color = :white, markersize = 10,
             strokecolor = :black, strokewidth = 2,
             label = "Grid samples")

    # Add vertical lines at grid points
    vlines!(ax, x_grid, color = (:gray, 0.3), linewidth = 0.5)

    # Horizontal zero line
    hlines!(ax, [0.0], color = :black, linewidth = 0.5)

    # Legend
    axislegend(ax, position = :rt, framevisible = false, bgcolor = (:white, 0.9))

    # Save figure
    mkpath(OUTPUT_DIR)
    save(joinpath(OUTPUT_DIR, "aliasing_demo.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "aliasing_demo.png"), fig, px_per_unit = 4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "aliasing_demo.pdf"))

    return fig
end

# Run if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
