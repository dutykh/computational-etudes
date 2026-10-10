# cheb_matrix_structure.jl
#
# Visualizes the structure of the Chebyshev differentiation matrix: heatmap
# showing the matrix entries and row profile showing the decay pattern.
#
# This script generates Figure 6.2 for Chapter 6.
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
using LinearAlgebra

# Import Chebyshev matrix function
include(joinpath(@__DIR__, "cheb_matrix.jl"))

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
const N = 16  # Number of intervals

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
    # Construct differentiation matrix
    D, x = cheb_matrix(N)

    # --- Figure 1: Heatmap + Row profiles ---
    fig = Figure(size = (800, 400))

    # Panel 1: Heatmap of D matrix
    ax1 = Axis(fig[1, 1],
               xlabel = L"Column index $j$",
               ylabel = L"Row index $i$",
               title = "Chebyshev Differentiation Matrix \$D_N\$ (\$N = $N\$)",
               yreversed = true,
               aspect = DataAspect())

    vmax = maximum(abs.(D))

    # Heatmap with diverging colormap centered at zero
    hm = heatmap!(ax1, 0:N, 0:N, D',
                  colormap = :RdBu, colorrange = (-vmax, vmax))
    Colorbar(fig[1, 1, Right()], hm, label = "Matrix entry value",
             labelsize = 10, width = 12)

    # Set ticks
    tick_positions = [0, N ÷ 4, N ÷ 2, 3 * N ÷ 4, N]
    ax1.xticks = (tick_positions, string.(tick_positions))
    ax1.yticks = (tick_positions, string.(tick_positions))

    # Annotations for corner entries
    text!(ax1, 4.0, 4.0,
          text = "\$D_{00} = $(round(D[1,1], digits=1))\$",
          fontsize = 9, color = :white, align = (:left, :center),
          font = "CMU Serif")
    text!(ax1, N - 5.0, N - 3.0,
          text = "\$D_{N,N} = $(round(D[N+1,N+1], digits=1))\$",
          fontsize = 9, color = :white, align = (:left, :center),
          font = "CMU Serif")

    # Panel 2: Row profiles (first row and middle row)
    ax2 = Axis(fig[1, 2],
               xlabel = L"Column index $j$",
               ylabel = L"Matrix entry $D_{ij}$",
               title = L"Row Profiles of $D_N$")
    hidespines!(ax2, :t, :r)

    j_indices = 0:N

    # First row (boundary)
    stemlines_x1 = Float64[]
    stemlines_y1 = Float64[]
    for j in j_indices
        push!(stemlines_x1, j); push!(stemlines_y1, 0.0)
        push!(stemlines_x1, j); push!(stemlines_y1, D[1, j + 1])
        push!(stemlines_x1, NaN); push!(stemlines_y1, NaN)
    end
    lines!(ax2, stemlines_x1, stemlines_y1, color = CORAL, linewidth = 1.0)
    scatter!(ax2, collect(Float64, j_indices), D[1, :],
             color = CORAL, markersize = 6, marker = :circle,
             label = L"Row $i=0$ (boundary)")

    # Middle row
    mid_row = N ÷ 2
    stemlines_x2 = Float64[]
    stemlines_y2 = Float64[]
    for j in j_indices
        push!(stemlines_x2, j + 0.15); push!(stemlines_y2, 0.0)
        push!(stemlines_x2, j + 0.15); push!(stemlines_y2, D[mid_row + 1, j + 1])
        push!(stemlines_x2, NaN); push!(stemlines_y2, NaN)
    end
    lines!(ax2, stemlines_x2, stemlines_y2, color = TEAL, linewidth = 1.0)
    scatter!(ax2, collect(Float64, j_indices) .+ 0.15, D[mid_row + 1, :],
             color = TEAL, markersize = 6, marker = :rect,
             label = "Row \$i=$mid_row\$ (interior)")

    # Reference line at zero
    hlines!(ax2, 0, color = RGBf(0.533, 0.533, 0.533), linewidth = 0.5)

    # Horizontal grid
    ax2.ygridvisible = true
    ax2.ygridcolor = RGBAf(0.533, 0.533, 0.533, 0.2)
    ax2.ygridstyle = :solid

    xlims!(ax2, -0.5, N + 0.5)
    axislegend(ax2, position = :rt, labelsize = 9, framevisible = false)

    # Save first figure
    save(joinpath(OUTPUT_DIR, "cheb_matrix_structure.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "cheb_matrix_structure.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "cheb_matrix_structure.pdf"))")

    # --- Figure 2: Interior row profile only ---
    fig2 = Figure(size = (500, 330))

    interior_row = 9  # 0-based index
    ax3 = Axis(fig2[1, 1],
               xlabel = L"Column index $j$",
               ylabel = L"Matrix entry $D_{9,j}$",
               title = L"Interior Row Profile ($i = 9$)")
    hidespines!(ax3, :t, :r)

    # Stem plot for row i=9
    stemlines_x3 = Float64[]
    stemlines_y3 = Float64[]
    for j in j_indices
        push!(stemlines_x3, j); push!(stemlines_y3, 0.0)
        push!(stemlines_x3, j); push!(stemlines_y3, D[interior_row + 1, j + 1])
        push!(stemlines_x3, NaN); push!(stemlines_y3, NaN)
    end
    lines!(ax3, stemlines_x3, stemlines_y3, color = TEAL, linewidth = 1.0)
    scatter!(ax3, collect(Float64, j_indices), D[interior_row + 1, :],
             color = TEAL, markersize = 8, marker = :rect)

    # Reference line at zero
    hlines!(ax3, 0, color = RGBf(0.533, 0.533, 0.533), linewidth = 0.5)

    # Horizontal grid
    ax3.ygridvisible = true
    ax3.ygridcolor = RGBAf(0.533, 0.533, 0.533, 0.2)
    ax3.ygridstyle = :solid

    xlims!(ax3, -0.5, N + 0.5)

    # Save second figure
    save(joinpath(OUTPUT_DIR, "cheb_matrix_interior_row.pdf"), fig2)
    save(joinpath(OUTPUT_DIR, "cheb_matrix_interior_row.png"), fig2, px_per_unit = 4)
    println("Interior row figure saved to: $(joinpath(OUTPUT_DIR, "cheb_matrix_interior_row.pdf"))")

    # Print matrix properties
    println("\nChebyshev Matrix Properties (N = $N):")
    println("-" ^ 50)
    println("  Matrix size: ($(N+1), $(N+1))")
    println("  Corner entry D[0,0]: $(round(D[1,1], digits=6))")
    println("  Corner entry D[N,N]: $(round(D[N+1,N+1], digits=6))")
    println("  Theoretical D[0,0] = (2N^2+1)/6: $(round((2*N^2+1)/6, digits=6))")
    println("  Max entry: $(round(maximum(D), digits=6))")
    println("  Min entry: $(round(minimum(D), digits=6))")
    println("  Negative sum trick error: $(round(maximum(abs.(D * ones(N+1))), sigdigits=2))")
    println("-" ^ 50)
end

main()
