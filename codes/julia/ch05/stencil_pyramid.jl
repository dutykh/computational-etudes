# stencil_pyramid.jl
#
# Visualizes Fornberg's recursive algorithm for computing finite difference
# weights as a "stencil pyramid" showing how weights are computed progressively
# from 1-point to N-point stencils.
#
# The visualization shows:
#     - The recursive structure of the algorithm
#     - How weights for smaller stencils inform larger stencils
#     - The data dependencies in the computation
#
# Key Insight:
#     The algorithm builds up weights level by level:
#     - Level 0: Single node (weight = 1 for interpolation)
#     - Level 1: Two nodes (weights from level 0)
#     - Level 2: Three nodes (weights from level 1)
#     - ...
#     - Level n: n+1 nodes (weights from level n-1)
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
const NAVY        = colorant"#142D6E"
const SKY         = colorant"#7896D2"
const CORAL       = colorant"#E74C3C"
const TEAL        = colorant"#16A085"
const LIGHT_GRAY  = colorant"#E8E8E8"
const MEDIUM_GRAY = colorant"#CCCCCC"

SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch05", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    fig = Figure(size = (720, 500))
    ax = Axis(fig[1, 1],
              title = "Fornberg's Recursive Algorithm: The Stencil Pyramid",
              titlesize = 13,
              aspect = DataAspect())

    hidedecorations!(ax)
    hidespines!(ax)

    # Parameters
    n_levels  = 5    # Number of levels (m = 0, 1, 2, 3, 4)
    box_width  = 0.8
    box_height = 0.5
    h_spacing  = 1.0
    v_spacing  = 1.2

    # Store box positions for drawing arrows
    boxes = Dict{Tuple{Int,Int}, Tuple{Float64,Float64}}()

    for m in 0:n_levels-1  # m = stencil size - 1
        n_boxes = m + 1
        y = (n_levels - 1 - m) * v_spacing

        # Center the row
        total_width = n_boxes * box_width + (n_boxes - 1) * (h_spacing - box_width)
        x_start = -total_width / 2 + box_width / 2

        for k in 0:n_boxes-1
            x = x_start + k * h_spacing

            # Choose color based on position
            if m == n_levels - 1
                # Final level (result) - highlight
                color = NAVY
                text_color = :white
                alpha = 1.0
            elseif m == n_levels - 2
                # Previous level - medium highlight
                color = SKY
                text_color = :white
                alpha = 0.9
            else
                # Earlier levels - faded
                color = LIGHT_GRAY
                text_color = NAVY
                alpha = 0.8
            end

            # Draw box as a rounded rectangle
            poly!(ax,
                  Point2f[(x - box_width/2, y - box_height/2),
                           (x + box_width/2, y - box_height/2),
                           (x + box_width/2, y + box_height/2),
                           (x - box_width/2, y + box_height/2)],
                  color = (color, alpha),
                  strokecolor = NAVY,
                  strokewidth = 1.5)

            # Add text label
            if m <= 2
                label_text = rich("δ", subscript("$m"),
                                  superscript("($k,$m)"))
            else
                label_text = rich("δ", superscript("($k,$m)"))
            end

            text!(ax, x, y,
                  text = label_text,
                  fontsize = 9,
                  color = text_color,
                  font = :bold,
                  align = (:center, :center))

            boxes[(m, k)] = (x, y)
        end
    end

    # Draw arrows showing dependencies
    for m in 1:n_levels-1
        n_boxes = m + 1
        for k in 0:n_boxes-1
            x_to, y_to = boxes[(m, k)]

            # Arrows from previous level
            if k < m
                x_from, y_from = boxes[(m - 1, k)]
                arrows!(ax, [x_from], [y_from - box_height/2 - 0.05],
                        [x_to - x_from], [y_to + box_height/2 + 0.05 - (y_from - box_height/2 - 0.05)],
                        color = (TEAL, 0.6), linewidth = 1.2,
                        arrowsize = 8)
            end

            if k == m && m > 0
                x_from, y_from = boxes[(m - 1, m - 1)]
                # Arrow from rightmost of previous level to new node
                arrows!(ax, [x_from + box_width/4],
                        [y_from - box_height/2 - 0.05],
                        [x_to - box_width/4 - (x_from + box_width/4)],
                        [y_to + box_height/2 + 0.05 - (y_from - box_height/2 - 0.05)],
                        color = (CORAL, 0.7), linewidth = 1.5,
                        arrowsize = 8)
            end
        end
    end

    # Add level labels on the left
    for m in 0:n_levels-1
        y = (n_levels - 1 - m) * v_spacing
        x = -4.5

        text!(ax, x, y,
              text = "m = $m",
              fontsize = 11, color = NAVY,
              align = (:left, :center))
        text!(ax, x + 1.2, y,
              text = "($(m + 1) node$(m > 0 ? "s" : ""))",
              fontsize = 9, color = :gray,
              align = (:left, :center))
    end

    # Add explanation text at bottom
    text!(ax, 0.0, -0.8,
          text = "Each box represents a weight for a node in a stencil.\nArrows show data dependencies: weights at level m depend on weights from level m-1.",
          fontsize = 10, color = :gray, font = :italic,
          align = (:center, :top))

    # Add legend entries
    lines!(ax, [Float64(NaN)], [Float64(NaN)], color = (TEAL, 0.6), linewidth = 2,
           label = "Update existing node")
    lines!(ax, [Float64(NaN)], [Float64(NaN)], color = (CORAL, 0.7), linewidth = 2,
           label = "Compute new node (with beta factor)")
    axislegend(ax, position = :rt, labelsize = 9, framevisible = true,
               framecolor = :gray, bgcolor = RGBAf(1, 1, 1, 0.9))

    # Set axis limits
    xlims!(ax, -5.5, 5.5)
    ylims!(ax, -1.5, (n_levels - 1) * v_spacing + 1)

    # Save figure
    save(joinpath(OUTPUT_DIR, "stencil_pyramid.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "stencil_pyramid.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "stencil_pyramid.pdf"))")
end

main()
