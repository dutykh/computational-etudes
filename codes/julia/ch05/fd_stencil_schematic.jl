# fd_stencil_schematic.jl
#
# Creates a schematic illustration of finite difference stencils in 1D,
# showing how the stencil width increases with the order of accuracy,
# culminating in the spectral method that uses all nodes.
#
# Illustrates:
#     - 3-point stencil (2nd order): uses x_{i-1}, x_i, x_{i+1}
#     - 5-point stencil (4th order): uses x_{i-2}, ..., x_{i+2}
#     - 7-point stencil (6th order): uses x_{i-3}, ..., x_{i+3}
#     - Spectral stencil: uses ALL nodes
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
const NAVY       = colorant"#142D6E"
const SKY        = colorant"#7896D2"
const CORAL      = colorant"#E74C3C"
const TEAL       = colorant"#16A085"
const LIGHT_GRAY = colorant"#BBBBBB"
const VERY_LIGHT = colorant"#DDDDDD"

SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch05", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Grid setup
    n_points = 11   # Total points to show (odd number for symmetry)
    center_idx = 6  # Index of central point (1-indexed)

    # Stencil configurations: (half_width, label, is_spectral)
    stencils = [
        (1,  "2nd order FD (3 points)", false),
        (2,  "4th order FD (5 points)", false),
        (3,  "6th order FD (7 points)", false),
        (-1, "Spectral (all points)",   true),
    ]

    n_stencils = length(stencils)

    # Create figure
    fig = Figure(size = (880, 620))

    # Point positions in data coordinates
    x_left  = 0.08
    x_right = 0.72
    x_positions = range(x_left, x_right, length = n_points)
    h_spacing = x_positions[2] - x_positions[1]

    for (ax_idx, (half_width, label, is_spectral)) in enumerate(stencils)
        ax = Axis(fig[ax_idx, 1],
                  limits = (0, 1, -1, 1),
                  aspect = nothing)
        hidedecorations!(ax)
        hidespines!(ax)

        y_base = 0.0

        # Draw horizontal grid line
        lines!(ax, [0.05, 0.77], [y_base, y_base],
               color = VERY_LIGHT, linewidth = 1.2, linestyle = :dash)

        # Determine active indices
        if is_spectral
            active_indices = Set(1:n_points)
        else
            active_indices = Set((center_idx - half_width):(center_idx + half_width))
        end

        # Draw shaded region for stencil
        if is_spectral
            rect_left  = x_positions[1] - 0.012
            rect_right = x_positions[end] + 0.012
        else
            rect_left  = x_positions[center_idx - half_width] - 0.012
            rect_right = x_positions[center_idx + half_width] + 0.012
        end

        poly!(ax,
              Point2f[(rect_left, y_base - 0.45),
                       (rect_right, y_base - 0.45),
                       (rect_right, y_base + 0.45),
                       (rect_left, y_base + 0.45)],
              color = (SKY, 0.15),
              strokecolor = SKY,
              strokewidth = 1.2)

        # Draw all grid points
        for (i, xp) in enumerate(x_positions)
            if i == center_idx
                color = CORAL
                sz = 15
            elseif i in active_indices
                color = NAVY
                sz = 12
            else
                color = LIGHT_GRAY
                sz = 8
            end

            scatter!(ax, [xp], [y_base], color = color, markersize = sz,
                     strokecolor = :white, strokewidth = 2)
        end

        # Add continuation dots
        text!(ax, x_positions[1] - 0.025, y_base,
              text = "...", fontsize = 16, color = LIGHT_GRAY,
              align = (:right, :center))
        text!(ax, x_positions[end] + 0.025, y_base,
              text = "...", fontsize = 16, color = LIGHT_GRAY,
              align = (:left, :center))

        # Add labels below the points
        label_y = y_base - 0.7

        if is_spectral
            labels_spec = [
                (1,          L"x_0", NAVY),
                (center_idx, L"x_i", CORAL),
                (n_points,   L"x_N", NAVY),
            ]
            for (idx, lbl, col) in labels_spec
                text!(ax, x_positions[idx], label_y,
                      text = lbl, fontsize = 12, color = col,
                      align = (:center, :top))
            end
        else
            left_idx  = center_idx - half_width
            right_idx = center_idx + half_width

            left_lbl  = half_width == 1 ? L"x_{i-1}" :
                        half_width == 2 ? L"x_{i-2}" : L"x_{i-3}"
            right_lbl = half_width == 1 ? L"x_{i+1}" :
                        half_width == 2 ? L"x_{i+2}" : L"x_{i+3}"

            text!(ax, x_positions[left_idx], label_y,
                  text = left_lbl, fontsize = 12, color = NAVY,
                  align = (:center, :top))
            text!(ax, x_positions[center_idx], label_y,
                  text = L"x_i", fontsize = 12, color = CORAL,
                  align = (:center, :top))
            text!(ax, x_positions[right_idx], label_y,
                  text = right_lbl, fontsize = 12, color = NAVY,
                  align = (:center, :top))
        end

        # Add method label on the right
        text!(ax, 0.88, y_base,
              text = label, fontsize = 12, color = NAVY,
              font = :bold, align = (:left, :center))

        # Draw bracket above stencil
        bracket_y = y_base + 0.55
        tick_len = 0.12
        lines!(ax, [rect_left, rect_left], [bracket_y - tick_len, bracket_y],
               color = NAVY, linewidth = 1.5)
        lines!(ax, [rect_right, rect_right], [bracket_y - tick_len, bracket_y],
               color = NAVY, linewidth = 1.5)
        lines!(ax, [rect_left, rect_right], [bracket_y, bracket_y],
               color = NAVY, linewidth = 1.5)

        # Grid spacing annotation (first row only)
        if ax_idx == 1
            hx1 = x_positions[center_idx]
            hx2 = x_positions[center_idx + 1]
            hy  = y_base + 0.85

            # Double-arrow for h
            arrows!(ax, [hx1], [hy], [hx2 - hx1], [0.0],
                    color = TEAL, linewidth = 2, arrowsize = 10)
            arrows!(ax, [hx2], [hy], [hx1 - hx2], [0.0],
                    color = TEAL, linewidth = 2, arrowsize = 10)
            text!(ax, (hx1 + hx2) / 2, hy + 0.12,
                  text = L"h", fontsize = 14, color = TEAL,
                  font = :bold, align = (:center, :bottom))
        end
    end

    # Main title
    Label(fig[0, :],
          "Finite Difference Stencils: From Local to Global",
          fontsize = 15, font = :bold, color = NAVY)

    # Bottom annotation
    Label(fig[n_stencils + 1, :],
          "As stencil width increases, accuracy improves. The spectral method is the limit: all nodes contribute to u'(x_i).",
          fontsize = 11, color = :gray, font = :italic)

    rowgap!(fig.layout, 10)

    # Save figure
    save(joinpath(OUTPUT_DIR, "fd_stencil_schematic.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "fd_stencil_schematic.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "fd_stencil_schematic.pdf"))")
end

main()
