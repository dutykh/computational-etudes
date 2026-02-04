# spectral_matrix_structure.jl
#
# Visualizes the structure of the periodic spectral differentiation matrix,
# showing its Toeplitz/circulant structure and the characteristic cotangent
# decay of the off-diagonal entries.
#
# Key Properties Illustrated:
#     - Skew-symmetry: D^T = -D
#     - Toeplitz structure: D[i,j] depends only on (i-j)
#     - Circulant (due to periodicity)
#     - Dense: all off-diagonal entries nonzero
#     - Diagonal entries are zero
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
const NAVY  = colorant"#142D6E"
const SKY   = colorant"#7896D2"
const CORAL = colorant"#E74C3C"
const TEAL  = colorant"#16A085"

SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch05", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Periodic spectral differentiation matrix
# -----------------------------------------------------------------------------
"""
    spectral_diff_periodic(N)

Construct the periodic spectral differentiation matrix for N equispaced
points on [0, 2pi).
"""
function spectral_diff_periodic(N)
    N < 2 && error("Need at least 2 points")

    h = 2pi / N
    x = h .* (0:N-1)

    D = zeros(N, N)
    for i in 1:N
        for j in 1:N
            if i != j
                diff = i - j
                D[i, j] = 0.5 * ((-1)^diff) / tan(diff * pi / N)
            end
        end
    end

    return D, x
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    N = 16  # Number of grid points

    D, x = spectral_diff_periodic(N)

    # Create figure with two panels
    fig = Figure(size = (820, 380))

    # =========================================================================
    # Left panel: Heatmap of the matrix
    # =========================================================================
    ax1 = Axis(fig[1, 1],
               xlabel = L"Column index $j$",
               ylabel = L"Row index $i$",
               title = "Spectral Matrix D (N = $N)",
               yreversed = true,
               aspect = DataAspect())

    # Use diverging colormap centered at 0
    vmax = maximum(abs.(D))

    hm = heatmap!(ax1, 0:N-1, 0:N-1, D',
                  colormap = :RdBu,
                  colorrange = (-vmax, vmax))

    Colorbar(fig[1, 1], hm, label = "Matrix entry value",
             labelsize = 10, width = 12,
             tellwidth = true,
             vertical = true,
             halign = :right)

    # Add annotation about skew-symmetry
    text!(ax1, 0.3, 0.3,
          text = L"D^T = -D",
          fontsize = 11,
          font = :bold,
          color = NAVY,
          space = :relative)

    # Add ticks
    ax1.xticks = [0, N ÷ 2, N - 1]
    ax1.yticks = [0, N ÷ 2, N - 1]

    # =========================================================================
    # Right panel: First row profile (showing the cotangent decay)
    # =========================================================================
    ax2 = Axis(fig[1, 2],
               xlabel = L"Column offset $k = j - i$",
               ylabel = L"Matrix entry $D_{0,k}$",
               title = "First Row Profile (Toeplitz Structure)")

    row = D[1, :]  # First row (Julia 1-indexed)
    j_indices = 0:N-1

    # Shift indices to show the pattern centered around 0
    j_shifted = [j <= N ÷ 2 ? j : j - N for j in j_indices]

    # Sort for nicer plotting
    sort_idx = sortperm(j_shifted)
    j_plot = j_shifted[sort_idx]
    row_plot = row[sort_idx]

    # Stem plot
    for (jp, rp) in zip(j_plot, row_plot)
        lines!(ax2, [jp, jp], [0.0, rp], color = NAVY, linewidth = 1.2)
    end
    scatter!(ax2, j_plot, row_plot, color = NAVY, markersize = 6)

    # Add theoretical cotangent curve
    j_theory = range(-N / 2 + 0.01, N / 2 - 0.01, length = 500)
    j_theory_filt = filter(j -> abs(j) > 0.1, j_theory)
    D_theory = [0.5 * ((-1)^round(Int, j)) / tan(j * pi / N) for j in j_theory_filt]
    lines!(ax2, j_theory_filt, D_theory, color = CORAL, linewidth = 1.5,
           linestyle = :dash, alpha = 0.7,
           label = L"\frac{1}{2}(-1)^k \cot\!\left(\frac{k\pi}{N}\right)")

    hlines!(ax2, [0.0], color = :gray, linewidth = 0.5)
    vlines!(ax2, [0.0], color = :gray, linewidth = 0.5)

    xlims!(ax2, -N / 2 - 0.5, N / 2 + 0.5)

    axislegend(ax2, position = :rt, labelsize = 9)

    # Clean styling
    for ax in [ax1, ax2]
        hidespines!(ax, :t, :r)
        ax.leftspinecolor  = RGBf(0.267, 0.267, 0.267)
        ax.bottomspinecolor = RGBf(0.267, 0.267, 0.267)
    end

    ax2.ygridvisible = true
    ax2.ygridcolor   = RGBAf(0.533, 0.533, 0.533, 0.2)
    ax2.ygridstyle   = :solid
    ax2.xgridvisible = true
    ax2.xgridcolor   = RGBAf(0.533, 0.533, 0.533, 0.2)
    ax2.xgridstyle   = :solid

    # Save figure
    save(joinpath(OUTPUT_DIR, "spectral_matrix_structure.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "spectral_matrix_structure.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "spectral_matrix_structure.pdf"))")

    # Print matrix properties
    println("\nSpectral Matrix Properties (N = $N):")
    println("-" ^ 50)
    @printf("  Skew-symmetry error ||D + D^T||_max: %.2e\n", maximum(abs.(D + D')))
    @printf("  Diagonal sum: %.2e (should be 0)\n", sum(diag(D)))
    @printf("  Max entry: %.4f\n", maximum(D))
    @printf("  Min entry: %.4f\n", minimum(D))
    println("-" ^ 50)
end

using Printf
main()
