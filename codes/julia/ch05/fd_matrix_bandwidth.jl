# fd_matrix_bandwidth.jl
#
# Visualizes the sparsity patterns of finite difference and spectral
# differentiation matrices, showing how the bandwidth expands as the
# order of accuracy increases.
#
# Key Message:
#     - 2nd order FD: Tridiagonal (bandwidth 3)
#     - 4th order FD: Pentadiagonal (bandwidth 5)
#     - Spectral: Dense (bandwidth N)
#
#     This illustrates that spectral methods are the limit case of finite
#     differences as the stencil width extends to the full domain.
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
using Printf

# Include shared utilities
include(joinpath(@__DIR__, "fdweights.jl"))

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

# Output path
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

Returns `(D, x)` where D is the N x N differentiation matrix and x are
the grid points.
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
# Finite difference matrix construction
# -----------------------------------------------------------------------------
"""
    fd_diff_periodic(N, order)

Construct a periodic finite difference differentiation matrix.
"""
function fd_diff_periodic(N, order)
    h = 2pi / N

    # Determine stencil size
    stencil_half = order ÷ 2

    # Create stencil nodes centered at 0
    stencil_nodes = h .* (-stencil_half:stencil_half)

    # Compute FD weights for first derivative at center
    weights = fdweights(0.0, collect(stencil_nodes), 1)

    # Build circulant matrix
    D = zeros(N, N)
    for i in 1:N
        for (k, w) in enumerate(weights)
            j = mod(i - 1 + k - 1 - stencil_half, N) + 1  # Periodic wrap-around (1-indexed)
            D[i, j] += w
        end
    end

    return D
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    N = 20  # Number of grid points

    # Create the three matrices
    D_fd2 = fd_diff_periodic(N, 2)   # 2nd order: tridiagonal
    D_fd4 = fd_diff_periodic(N, 4)   # 4th order: pentadiagonal
    D_spec, _ = spectral_diff_periodic(N)  # Spectral: dense

    # Create figure with three panels
    fig = Figure(size = (880, 300))

    matrices = [D_fd2, D_fd4, D_spec]
    titles = ["2nd Order FD", "4th Order FD", "Spectral"]
    bandwidths = ["Bandwidth = 3", "Bandwidth = 5", "Dense (N x N)"]

    for (idx, (D, title, bw)) in enumerate(zip(matrices, titles, bandwidths))
        ax = Axis(fig[1, idx],
                  title = title,
                  titlesize = 12,
                  xlabel = "Column index j\n$bw",
                  xlabelsize = 10,
                  ylabel = idx == 1 ? "Row index i" : "",
                  ylabelsize = 10,
                  yreversed = true,
                  aspect = DataAspect())

        # Create binary sparsity pattern
        threshold = 1e-14
        sparsity = Float64.(abs.(D) .> threshold)

        # Plot using heatmap
        heatmap!(ax, 0:N-1, 0:N-1, sparsity',
                 colormap = [:white, NAVY],
                 colorrange = (0, 1))

        # Reduce tick density
        ax.xticks = [0, N ÷ 2, N - 1]
        ax.yticks = [0, N ÷ 2, N - 1]
    end

    # Add overall title
    Label(fig[0, :], "Finite Difference Matrix Sparsity: From Local to Global",
          fontsize = 13, font = :bold)

    # Save figure
    save(joinpath(OUTPUT_DIR, "fd_matrix_bandwidth.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "fd_matrix_bandwidth.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "fd_matrix_bandwidth.pdf"))")

    # Print matrix statistics
    println("\nMatrix Statistics (N = $N):")
    println("-" ^ 50)
    for (D, title) in zip(matrices, titles)
        nnz = count(x -> abs(x) > 1e-14, D)
        density = nnz / length(D) * 100
        println("$(rpad(title, 15)): $(lpad(string(nnz), 4)) nonzeros ($(lpad(@sprintf("%.1f", density), 5))% dense)")
    end
    println("-" ^ 50)
end

main()
