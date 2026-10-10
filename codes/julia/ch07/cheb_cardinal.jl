# cheb_cardinal.jl
#
# Visualizes Chebyshev cardinal functions (Lagrange interpolation basis) and
# shows how the differentiation matrix entries are the derivatives of these
# cardinal functions evaluated at the grid points.
#
# This script generates Figure 6.3 for Chapter 6.
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
const N = 10       # Number of intervals (smaller for clearer visualization)
const N_FINE = 500 # Fine grid for plotting

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
# Chebyshev cardinal function evaluation via barycentric interpolation
# -----------------------------------------------------------------------------
"""
    chebyshev_cardinal(x_eval, x_nodes, j)

Evaluate the j-th Chebyshev cardinal function at points `x_eval`.

The cardinal function L_j(x) is the unique polynomial of degree N
that satisfies L_j(x_k) = delta_{jk} (Kronecker delta).

Uses the barycentric interpolation formula for Chebyshev points.
`j` is 1-based (Julia convention).
"""
function chebyshev_cardinal(x_eval, x_nodes, j)
    n = length(x_nodes) - 1

    # Barycentric weights for Chebyshev points:
    # w_k = (-1)^k * (1/2 if k=0,N else 1)
    w = ones(n + 1)
    w[1] = 0.5
    w[n + 1] = 0.5
    for k in 0:n
        w[k + 1] *= (-1)^k
    end

    L_j = zeros(length(x_eval))

    for (i, xv) in enumerate(x_eval)
        diffs = xv .- x_nodes
        # Check if xv is at a node
        close_idx = findfirst(d -> abs(d) < 1e-14, diffs)
        if close_idx !== nothing
            L_j[i] = (close_idx == j) ? 1.0 : 0.0
        else
            # Barycentric formula
            terms = w ./ diffs
            L_j[i] = (w[j] / diffs[j]) / sum(terms)
        end
    end

    return L_j
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Chebyshev grid points
    D, x_nodes = cheb_matrix(N)
    x_fine = collect(range(-1, 1, length=N_FINE))

    # Create figure with two panels
    fig = Figure(size = (800, 370))

    # Panel 1: Several cardinal functions
    ax1 = Axis(fig[1, 1],
               xlabel = L"x",
               ylabel = L"L_j(x)",
               title = "Chebyshev Cardinal Functions (\$N = $N\$)")
    hidespines!(ax1, :t, :r)
    xlims!(ax1, -1.05, 1.05)
    ylims!(ax1, -0.3, 1.2)

    # Horizontal grid
    ax1.ygridvisible = true
    ax1.ygridcolor = RGBAf(0.533, 0.533, 0.533, 0.2)
    ax1.ygridstyle = :solid

    # Select a few cardinal functions to display (1-based indices)
    indices = [1, 4, 6, N + 1]   # corresponding to 0-based: 0, 3, 5, N
    labels_0based = [0, 3, 5, N]
    colors_list = [CORAL, TEAL, PURPLE, ORANGE]
    linestyles = [:solid, :dash, :dashdot, :dot]

    for (idx, jl, c, ls) in zip(indices, labels_0based, colors_list, linestyles)
        L_j = chebyshev_cardinal(x_fine, x_nodes, idx)
        lines!(ax1, x_fine, L_j, color = c, linewidth = 1.5, linestyle = ls,
               label = "\$L_{$jl}(x)\$")
        # Mark the cardinal point where L_idx = 1
        scatter!(ax1, [x_nodes[idx]], [1.0], color = c, markersize = 9,
                 strokecolor = :white, strokewidth = 1)
    end

    # Mark all nodes on x-axis
    scatter!(ax1, x_nodes, zeros(N + 1), color = NAVY, markersize = 6,
             strokecolor = :white, strokewidth = 0.5)

    # Reference lines
    hlines!(ax1, 0, color = RGBf(0.533, 0.533, 0.533), linewidth = 0.5)
    hlines!(ax1, 1, color = RGBf(0.533, 0.533, 0.533), linewidth = 0.5,
            linestyle = :dash, alpha = 0.5)

    axislegend(ax1, position = :rt, labelsize = 9, nbanks = 2, framevisible = false)

    # Annotation
    text!(ax1, 0.5, 0.85 * 1.2,
          text = L"L_j(x_k) = \delta_{jk}",
          fontsize = 10, color = RGBf(0.4, 0.4, 0.4), align = (:center, :center),
          space = :relative)

    # Panel 2: Cardinal function with derivative slopes at nodes
    ax2 = Axis(fig[1, 2],
               xlabel = L"x",
               ylabel = L"L_j(x)",
               title = "Cardinal Function \$L_5\$ with Derivative Slopes")
    hidespines!(ax2, :t, :r)
    xlims!(ax2, -1.05, 1.05)
    ylims!(ax2, -0.5, 1.2)

    # Horizontal grid
    ax2.ygridvisible = true
    ax2.ygridcolor = RGBAf(0.533, 0.533, 0.533, 0.2)
    ax2.ygridstyle = :solid

    # Plot cardinal function L_5 (1-based index 6)
    j = 6  # 1-based Julia index, 0-based = 5
    L_j = chebyshev_cardinal(x_fine, x_nodes, j)
    lines!(ax2, x_fine, L_j, color = TEAL, linewidth = 2, label = L"L_5(x)")

    # Mark the cardinal point
    scatter!(ax2, [x_nodes[j]], [1.0], color = TEAL, markersize = 12,
             strokecolor = :white, strokewidth = 2)

    # Draw tangent lines at grid points (slopes from D matrix)
    for k in 1:2:N+1  # Every other point for clarity (1-based, step 2)
        slope = D[k, j]  # This is L_j'(x_k)
        x_k = x_nodes[k]
        y_k = chebyshev_cardinal([x_k], x_nodes, j)[1]

        # Draw short tangent line
        dx = 0.15
        x_tan = [x_k - dx, x_k + dx]
        y_tan = y_k .+ slope .* (x_tan .- x_k)

        lines!(ax2, x_tan, y_tan, color = CORAL, linewidth = 1.5, alpha = 0.7)
        scatter!(ax2, [x_k], [y_k], color = NAVY, markersize = 5, marker = :rect)
    end

    # Reference line
    hlines!(ax2, 0, color = RGBf(0.533, 0.533, 0.533), linewidth = 0.5)

    # Annotation
    text!(ax2, -0.95, -0.35,
          text = "Slopes = column 5 of \$D_N\$:\n\$D_{kj} = L_j'(x_k)\$",
          fontsize = 9, color = RGBf(0.4, 0.4, 0.4), align = (:left, :center))

    # Save figure
    save(joinpath(OUTPUT_DIR, "cheb_cardinal.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "cheb_cardinal.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "cheb_cardinal.pdf"))")

    # Verify cardinal property
    println("\nVerification of cardinal property (N = $N):")
    println("-" ^ 50)
    println("$(rpad("j", 4)) $(lpad("L_j(x_j)", 12)) $(lpad("max|L_j(x_k)|, k!=j", 20))")
    println("-" ^ 50)
    for jj in 1:N+1
        L_j_at_nodes = chebyshev_cardinal(x_nodes, x_nodes, jj)
        val_at_own = L_j_at_nodes[jj]
        val_at_others = [L_j_at_nodes[k] for k in 1:N+1 if k != jj]
        max_at_others = maximum(abs.(val_at_others))
        println("$(rpad(jj-1, 4)) $(lpad(round(val_at_own, digits=6), 12)) $(lpad(round(max_at_others, sigdigits=2), 20))")
    end
    println("-" ^ 50)
end

main()
