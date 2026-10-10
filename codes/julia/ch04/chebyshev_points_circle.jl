# chebyshev_points_circle.jl
#
# Visualizes the geometric construction of Chebyshev points from the unit circle.
#
# Chebyshev-Gauss-Lobatto points are obtained by:
# 1. Placing N+1 equally spaced points on the upper half of the unit circle
# 2. Projecting these points vertically onto the x-axis
#
# This geometric interpretation reveals why Chebyshev points cluster near +/-1:
# equal spacing on the circle maps to denser spacing near the boundaries when
# projected horizontally.
#
# The point density on [-1,1] follows:
#     rho(x) = N / (pi * sqrt(1 - x^2))
#
# which diverges at the boundaries, exactly compensating for the growth of
# interpolation error there.
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
const N = 8  # Number of intervals (N+1 points)

# Book color scheme
const NAVY       = colorant"#142D6E"
const SKY        = colorant"#7896D2"
const CORAL      = colorant"#E74C3C"
const TEAL       = colorant"#16A085"
const LIGHT_GRAY = colorant"#CCCCCC"

# Output path
const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch04", "julia")

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Create figure
    fig = Figure(size = (800, 500))
    ax = Axis(fig[1, 1],
              xlabel = L"x",
              title = "Geometric Construction of Chebyshev Points (N = $N)",
              limits = ((-1.4, 1.4), (-0.45, 1.15)),
              aspect = DataAspect())

    # Draw unit circle (upper half)
    theta_circle = range(0, π, length=200)
    x_circle = cos.(theta_circle)
    y_circle = sin.(theta_circle)
    lines!(ax, x_circle, y_circle, color = NAVY, linewidth = 1.5,
           label = "Unit circle")

    # Draw x-axis
    lines!(ax, [-1.3, 1.3], [0.0, 0.0], color = colorant"#666666", linewidth = 0.8)

    # Chebyshev points: angles equally spaced on [0, pi]
    theta_j = [j * π / N for j in 0:N]

    # Points on the circle
    x_circle_pts = cos.(theta_j)
    y_circle_pts = sin.(theta_j)

    # Points on the x-axis (Chebyshev nodes)
    x_cheb = cos.(theta_j)
    y_cheb = zeros(N + 1)

    # Draw projection lines
    for i in 1:N+1
        lines!(ax, [x_circle_pts[i], x_cheb[i]], [y_circle_pts[i], y_cheb[i]],
               color = LIGHT_GRAY, linewidth = 0.8, linestyle = :dash)
    end

    # Draw points on circle
    scatter!(ax, x_circle_pts, y_circle_pts, color = SKY, markersize = 12,
             strokecolor = NAVY, strokewidth = 1, label = "Points on circle")

    # Draw Chebyshev points on axis
    scatter!(ax, x_cheb, y_cheb, color = CORAL, markersize = 14,
             strokecolor = colorant"#B03A2E", strokewidth = 1.5, marker = :rect,
             label = "Chebyshev points")

    # Draw angle arc for middle point (j=4 for N=8)
    arc_idx = N ÷ 2 + 1  # Julia 1-based index
    arc_radius = 0.4
    arc_theta = range(0, theta_j[arc_idx], length=50)
    lines!(ax, arc_radius .* cos.(arc_theta), arc_radius .* sin.(arc_theta),
           color = TEAL, linewidth = 1.2)

    # Annotate the angle
    arc_label_angle = theta_j[arc_idx] / 2
    arc_label_x = (arc_radius + 0.15) * cos(arc_label_angle)
    arc_label_y = (arc_radius + 0.15) * sin(arc_label_angle)
    text!(ax, arc_label_x, arc_label_y,
          text = L"\theta_j = \frac{j\pi}{N}",
          fontsize = 10, color = TEAL, align = (:left, :center))

    # Label selected points
    label_indices = [1, 3, N÷2+1, N-1, N+1]  # 1-based: j=0,2,4,6,8
    for i in label_indices
        if i <= length(x_cheb)
            j = i - 1  # 0-based index for display
            if j == N ÷ 2  # Center point
                offset_y = -0.15
                offset_x = 0.08
            elseif j == 0  # rightmost
                offset_y = -0.18
                offset_x = 0.0
            elseif j == N  # leftmost
                offset_y = -0.18
                offset_x = 0.0
            else
                offset_y = 0.18
                offset_x = 0.0
            end
            text!(ax, x_cheb[i] + offset_x, offset_y,
                  text = latexstring("j = $j"),
                  fontsize = 9, align = (:center, :center), color = colorant"#444444")
        end
    end

    # Add formula box
    text!(ax, 0.98, 0.95,
          text = L"x_j = \cos\left(\frac{j\pi}{N}\right)" * "\n" * L"j = 0, 1, \ldots, N",
          fontsize = 11, space = :relative, align = (:right, :top))

    # Annotation: "Equal spacing on circle"
    text!(ax, -0.35, 0.95, text = "Equal spacing\non circle",
          fontsize = 9, color = SKY, align = (:center, :center), font = :italic)

    # Annotation: "Clustering near boundaries"
    text!(ax, 0.65, -0.28, text = "Clustering\nnear boundaries",
          fontsize = 9, color = CORAL, align = (:center, :center), font = :italic)

    # Remove spines
    hidespines!(ax)

    # Custom ticks
    ax.yticks = [0, 0.5, 1]
    ax.xticks = [-1, -0.5, 0, 0.5, 1]

    axislegend(ax, position = :lt, framevisible = false, bgcolor = (:white, 0.9))

    # Save figure
    mkpath(OUTPUT_DIR)
    save(joinpath(OUTPUT_DIR, "chebyshev_points_circle.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "chebyshev_points_circle.png"), fig, px_per_unit = 4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "chebyshev_points_circle.pdf"))

    # Print node positions
    println("\nChebyshev-Gauss-Lobatto points for N = $N:")
    println("-" ^ 40)
    println("$(rpad("j", 4)) $(lpad("theta_j", 12)) $(lpad("x_j", 12))")
    println("-" ^ 40)
    for i in 0:N
        theta_str = i == 0 ? "0" : (i == N ? "pi" : "$(i)pi/$N")
        println("$(rpad(i, 4)) $(lpad(theta_str, 12)) $(lpad(round(x_cheb[i+1], digits=6), 12))")
    end
    println("-" ^ 40)

    # Show clustering effect
    println("\nSpacing between consecutive points:")
    println("-" ^ 40)
    sorted_x = sort(x_cheb)
    spacings = diff(sorted_x)
    for (i, s) in enumerate(spacings)
        println("Interval $(i-1): dx = $(round(s, digits=6))")
    end
    println("-" ^ 40)
    println("Ratio (boundary/center): $(round(spacings[1] / spacings[N÷2], digits=2))")

    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
