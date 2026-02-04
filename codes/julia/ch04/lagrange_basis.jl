# lagrange_basis.jl
#
# Visualizes Lagrange basis polynomials for equispaced vs Chebyshev nodes.
#
# The Lagrange basis polynomials are:
#     L_k(x) = prod_{j != k} (x - x_j) / (x_k - x_j)
#
# These polynomials satisfy L_k(x_j) = delta_{kj} (Kronecker delta), making them
# the "cardinal functions" for interpolation: the interpolant is simply
#     p_N(x) = sum_{k=0}^{N} f_k L_k(x)
#
# For equispaced nodes, the L_k(x) develop large oscillations near the boundaries,
# contributing to the Runge phenomenon. For Chebyshev nodes, the basis functions
# remain well-behaved across the entire interval.
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
const N = 10       # Number of intervals (N+1 nodes)
const N_FINE = 1000

# Book color scheme
const NAVY   = colorant"#142D6E"
const SKY    = colorant"#7896D2"
const CORAL  = colorant"#E74C3C"
const TEAL   = colorant"#16A085"
const PURPLE = colorant"#8E44AD"
const ORANGE = colorant"#E67E22"

# Output path
const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch04", "julia")

# -----------------------------------------------------------------------------
# Node generation
# -----------------------------------------------------------------------------
equispaced_nodes(N) = collect(range(-1, 1, length=N+1))
chebyshev_nodes(N) = [cos(j * π / N) for j in 0:N]

# -----------------------------------------------------------------------------
# Lagrange basis computation
# -----------------------------------------------------------------------------
"""
    lagrange_basis(x_nodes, k, x_eval)

Compute the k-th Lagrange basis polynomial at points `x_eval`.
    L_k(x) = prod_{j != k} (x - x_j) / (x_k - x_j)

Note: `k` uses 1-based indexing.
"""
function lagrange_basis(x_nodes, k, x_eval)
    n = length(x_nodes)
    L_k = ones(length(x_eval))

    for j in 1:n
        if j != k
            L_k .*= (x_eval .- x_nodes[j]) ./ (x_nodes[k] - x_nodes[j])
        end
    end

    return L_k
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Fine grid for plotting
    x_fine = collect(range(-1, 1, length=N_FINE))

    # Generate nodes
    x_equi = equispaced_nodes(N)
    x_cheb = chebyshev_nodes(N)

    # Create figure with two panels
    fig = Figure(size = (1100, 450))

    # Indices of basis functions to plot (1-based: corresponds to k=0,1,N/2,N-1,N)
    indices_to_plot = [1, 2, N÷2+1, N, N+1]
    colors = [CORAL, ORANGE, TEAL, SKY, PURPLE]
    # Labels use 0-based k for display
    labels_k = [0, 1, N÷2, N-1, N]

    # Left panel: Equispaced basis functions
    ax1 = Axis(fig[1, 1],
               xlabel = L"x", ylabel = L"L_k(x)",
               title = "Equispaced Nodes (N = $N)",
               limits = ((-1, 1), (-2.5, 2.5)))

    for (k, col, k_label) in zip(indices_to_plot, colors, labels_k)
        L_k = lagrange_basis(x_equi, k, x_fine)
        lines!(ax1, x_fine, L_k, color = col, linewidth = 1.2,
               label = latexstring("L_{$k_label}(x)"))
    end

    # Mark nodes
    scatter!(ax1, x_equi, zeros(N+1), color = NAVY, markersize = 8)

    hlines!(ax1, 0, color = colorant"#888888", linewidth = 0.5)
    hlines!(ax1, 1, color = colorant"#888888", linewidth = 0.5, linestyle = :dash)

    # Annotation
    text!(ax1, -0.5, -2.2, text = "Large oscillations\nnear boundaries",
          fontsize = 9, color = :gray)

    axislegend(ax1, position = :lt, framevisible = false, labelsize = 9, nbanks = 2)

    # Right panel: Chebyshev basis functions
    ax2 = Axis(fig[1, 2],
               xlabel = L"x", ylabel = L"L_k(x)",
               title = "Chebyshev Nodes (N = $N)",
               limits = ((-1, 1), (-2.5, 2.5)))

    for (k, col, k_label) in zip(indices_to_plot, colors, labels_k)
        L_k = lagrange_basis(x_cheb, k, x_fine)
        lines!(ax2, x_fine, L_k, color = col, linewidth = 1.2,
               label = latexstring("L_{$k_label}(x)"))
    end

    # Mark nodes
    scatter!(ax2, x_cheb, zeros(N+1), color = NAVY, markersize = 8)

    hlines!(ax2, 0, color = colorant"#888888", linewidth = 0.5)
    hlines!(ax2, 1, color = colorant"#888888", linewidth = 0.5, linestyle = :dash)

    # Annotation
    text!(ax2, 0.3, 1.5, text = "Bounded oscillations",
          fontsize = 9, color = :gray)

    axislegend(ax2, position = :lt, framevisible = false, labelsize = 9, nbanks = 2)

    # Save figure
    mkpath(OUTPUT_DIR)
    save(joinpath(OUTPUT_DIR, "lagrange_basis.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "lagrange_basis.png"), fig, px_per_unit = 4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "lagrange_basis.pdf"))

    # Print maximum values of basis functions
    println("\nMaximum |L_k(x)| for N = $N:")
    println("-" ^ 50)
    println("$(rpad("k", 4)) $(lpad("Equispaced", 15)) $(lpad("Chebyshev", 15))")
    println("-" ^ 50)
    for k in 1:N+1
        L_equi = lagrange_basis(x_equi, k, x_fine)
        L_cheb = lagrange_basis(x_cheb, k, x_fine)
        println("$(rpad(k-1, 4)) $(lpad(round(maximum(abs.(L_equi)), digits=4), 15)) $(lpad(round(maximum(abs.(L_cheb)), digits=4), 15))")
    end
    println("-" ^ 50)

    # Compute Lebesgue functions
    lebesgue_equi = sum(abs.(lagrange_basis(x_equi, k, x_fine)) for k in 1:N+1)
    lebesgue_cheb = sum(abs.(lagrange_basis(x_cheb, k, x_fine)) for k in 1:N+1)

    println("\nLebesgue constants (Lambda_N = max sum|L_k|):")
    println("  Equispaced: Lambda_$N = $(round(maximum(lebesgue_equi), digits=2))")
    println("  Chebyshev:  Lambda_$N = $(round(maximum(lebesgue_cheb), digits=2))")
    println("  Ratio: $(round(maximum(lebesgue_equi) / maximum(lebesgue_cheb), digits=1))x")

    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
