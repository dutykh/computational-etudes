# eye_aliasing.jl
#
# Chapter 9: Physical and Fourier Space on Grids
#
# Demonstrates visual aliasing: sin(n) for integer n appears to oscillate
# slowly when plotted, because the eye aliases the high frequency to a
# low-frequency pattern.
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

# Output path
const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch09", "julia")

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Integer arguments: sin(n) for n = 1, 2, ..., 3000
    n = collect(1:3000)
    y = sin.(n)

    # Create figure
    fig = Figure(size = (600, 280))
    ax = Axis(fig[1, 1],
              xlabel = L"n",
              ylabel = L"\sin(n)",
              title  = L"Visual aliasing: $\sin(n)$ for integer $n$ appears to oscillate slowly",
              limits = ((0, 3000), (-1.3, 1.3)))

    scatter!(ax, n, y, color = NAVY, markersize = 1)

    # Save figure
    mkpath(OUTPUT_DIR)
    save(joinpath(OUTPUT_DIR, "eye_aliasing.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "eye_aliasing.png"), fig, px_per_unit = 4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "eye_aliasing.pdf"))

    return fig
end

# Run if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
