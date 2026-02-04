# sinc_interpolation.jl
#
# Chapter 9: Physical and Fourier Space on Grids
#
# Demonstrates band-limited interpolation using the sinc function.
# Shows interpolation of delta, square wave, and hat functions.
# This is Etude 2 from the chapter.
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
# Sinc kernel
# -----------------------------------------------------------------------------
"""
    sinc_kernel(z)

Compute sinc(z) = sin(pi*z) / (pi*z), with sinc(0) = 1.
"""
function sinc_kernel(z::AbstractArray)
    out = ones(length(z))
    for i in eachindex(z)
        if z[i] != 0
            out[i] = sin(π * z[i]) / (π * z[i])
        end
    end
    return out
end

# -----------------------------------------------------------------------------
# Sinc interpolation
# -----------------------------------------------------------------------------
"""
    sinc_interpolate(x_grid, v, x_fine; h=1.0)

Band-limited interpolation via sinc functions.

# Arguments
- `x_grid`: grid points
- `v`: function values at grid points
- `x_fine`: fine grid for interpolation
- `h`: grid spacing
"""
function sinc_interpolate(x_grid, v, x_fine; h=1.0)
    p = zeros(length(x_fine))
    for (xi, vi) in zip(x_grid, v)
        p .+= vi .* sinc_kernel((x_fine .- xi) ./ h)
    end
    return p
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Grid parameters
    h    = 1.0
    xmax = 10
    x_grid = collect(-xmax:h:xmax)
    x_fine = collect(-xmax - h/20 : h/10 : xmax + h/20)

    # Define three grid functions
    functions = [
        ("Discrete delta function",  Float64.(x_grid .== 0)),
        ("Discrete square wave",     Float64.(abs.(x_grid) .<= 3)),
        ("Discrete hat function",    max.(0.0, 1.0 .- abs.(x_grid) ./ 3.0)),
    ]

    # Create figure with 3 panels
    fig = Figure(size = (600, 640))

    for (i, (title_str, v)) in enumerate(functions)
        ax = Axis(fig[i, 1],
                  title  = title_str,
                  limits = ((-xmax, xmax), (-0.5, 1.5)),
                  yticks = [0, 0.5, 1])

        if i == 3
            ax.xlabel = L"x"
        end

        # Compute interpolant
        p = sinc_interpolate(x_grid, v, x_fine; h=h)

        # Plot discrete samples
        scatter!(ax, x_grid, v, color = :black, markersize = 8,
                 label = "Grid values")

        # Plot interpolant
        lines!(ax, x_fine, p, color = NAVY, linewidth = 1.0,
               label = "Band-limited interpolant")

        # Zero line
        hlines!(ax, [0.0], color = :black, linewidth = 0.5)

        # Legend
        axislegend(ax, position = :rt, framevisible = false,
                   bgcolor = (:white, 0.9), labelsize = 9)
    end

    # Save figure
    mkpath(OUTPUT_DIR)
    save(joinpath(OUTPUT_DIR, "sinc_interpolation.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "sinc_interpolation.png"), fig, px_per_unit = 4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "sinc_interpolation.pdf"))

    return fig
end

# Run if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
