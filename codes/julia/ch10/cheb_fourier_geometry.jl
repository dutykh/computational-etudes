# cheb_fourier_geometry.jl
#
# Chapter 10: Spectral PDE Solvers with Chebyshev and Fourier Grids
#
# Creates Figure 10.1: The Chebyshev-Fourier geometric connection.
# Shows the unit circle projection and the "wrapped cosine" interpretation.
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
const NAVY = colorant"#142D6E"
const SKY  = colorant"#7896D2"
const CORAL = colorant"#E74C3C"
const TEAL = colorant"#16A085"

# Output path
SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch10", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    fig = Figure(size = (880, 360))

    # =========================================================================
    # Left panel: Circle projection to Chebyshev points
    # =========================================================================
    ax1 = Axis(fig[1, 1],
               xlabel = L"x = \cos\theta",
               ylabel = L"y = \sin\theta",
               title  = "Chebyshev points as circle projection",
               aspect = DataAspect())

    N = 8
    theta = [pi * j / N for j in 0:N]
    x_cheb = cos.(theta)
    y_cheb = sin.(theta)

    # Draw unit semicircle
    theta_fine = range(0, pi, length=200)
    lines!(ax1, cos.(theta_fine), sin.(theta_fine), color = :black, linewidth = 1.2)
    hlines!(ax1, [0.0], color = :black, linewidth = 0.8)

    # Draw projection lines and points
    for j in 1:N+1
        # Projection line
        lines!(ax1, [x_cheb[j], x_cheb[j]], [0.0, y_cheb[j]],
               color = SKY, linewidth = 0.8, linestyle = :dash)
        # Point on circle
        scatter!(ax1, [x_cheb[j]], [y_cheb[j]], color = NAVY, markersize = 10)
        # Point on x-axis
        scatter!(ax1, [x_cheb[j]], [0.0], color = CORAL, markersize = 9, marker = :rect)
    end

    # Annotate some angles
    text!(ax1, 1.1, 0.2, text = L"\theta_0 = 0", fontsize = 10, color = NAVY)
    text!(ax1, -1.4, 0.2, text = L"\theta_N = \pi", fontsize = 10, color = NAVY)

    # Show angle arc
    theta_arc = range(0, pi/4, length=50)
    lines!(ax1, 0.3 .* cos.(theta_arc), 0.3 .* sin.(theta_arc),
           color = NAVY, linewidth = 1)
    text!(ax1, 0.35, 0.15, text = L"\theta", fontsize = 11, color = NAVY)

    # Label Chebyshev points
    text!(ax1, 1.0, -0.2, text = L"x_0 = 1", fontsize = 9, color = CORAL, align = (:center, :top))
    text!(ax1, -1.0, -0.2, text = L"x_N = -1", fontsize = 9, color = CORAL, align = (:center, :top))

    xlims!(ax1, -1.6, 1.6)
    ylims!(ax1, -0.3, 1.3)

    # =========================================================================
    # Right panel: Wrapped cosine interpretation
    # =========================================================================
    ax2 = Axis(fig[1, 2],
               xlabel = L"\theta",
               ylabel = L"\cos(n\theta) = T_n(x)",
               title  = L"Chebyshev polynomial $T_5(x)$ as wrapped cosine")

    # Draw T_n(x) = cos(n*theta) as a wrapped cosine
    n = 5
    theta_fine2 = range(0, pi, length=500)
    y_cosine = cos.(n .* theta_fine2)

    lines!(ax2, collect(theta_fine2), y_cosine, color = NAVY, linewidth = 1.5,
           label = L"\cos(n\theta)")

    # Mark theta values on bottom axis
    ax2.xticks = ([0, pi/4, pi/2, 3pi/4, pi],
                  [L"0", L"\pi/4", L"\pi/2", L"3\pi/4", L"\pi"])

    hlines!(ax2, [0.0], color = :black, linewidth = 0.5)
    ylims!(ax2, -1.3, 1.3)

    # Top axis showing x = cos(theta)
    ax2_top = Axis(fig[1, 2],
                   xaxisposition = :top,
                   xlabel = L"x = \cos\theta",
                   xlabelsize = 11,
                   xticklabelsize = 10)
    hideydecorations!(ax2_top)
    hidespines!(ax2_top, :b, :l, :r)
    ax2_top.xticks = ([0, pi/4, pi/2, 3pi/4, pi],
                      ["1.00", "0.71", "0.00", "-0.71", "-1.00"])

    # Link x-axes so they share the same limits
    linkxaxes!(ax2, ax2_top)

    # Save
    save(joinpath(OUTPUT_DIR, "cheb_fourier_geometry.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "cheb_fourier_geometry.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "cheb_fourier_geometry.pdf"))")
end

main()
