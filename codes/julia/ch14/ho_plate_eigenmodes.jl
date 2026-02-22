# ho_plate_eigenmodes.jl
#
# Chapter 14: Higher-Order Boundary Value Problems
#
# Eigenmodes of a clamped square plate -- biharmonic eigenvalue problem
# (Computational Etude 14.4).
#
#     Delta^2 u = lambda u,   (x, y) in (-1, 1)^2
#     u = 0,  du/dn = 0       on the boundary
#
# The clamped plate problem models the vibrations of a thin elastic plate
# with all four edges rigidly fixed.  The biharmonic operator Delta^2 is
# the composition of two Laplacians.
#
# Discretization follows Trefethen's Program 39 (Spectral Methods in MATLAB):
#   1. Build 1D clamped biharmonic operator D4 via a polynomial trick
#      that automatically enforces u = du/dn = 0 on the boundary.
#   2. Assemble the 2D operator via Kronecker products:
#      L = I x D4 + D4 x I + 2 * D2int x D2int
#   3. Solve the eigenvalue problem L v = lambda v.
#   4. Display the first 25 eigenmodes as nodal-line contour plots.
#
# Generates Figure 14.4: Eigenmodes of a clamped square plate.
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

# Import Chebyshev matrix function from Chapter 7
include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

# Publication-quality CairoMakie configuration
set_theme!(Theme(
    fontsize = 11,
    fonts = (regular = "CMU Serif", bold = "CMU Serif Bold", italic = "CMU Serif Italic"),
    Axis = (xlabelsize = 11, ylabelsize = 11, titlesize = 12,
            xticklabelsize = 10, yticklabelsize = 10,
            spinewidth = 0.8, xtickwidth = 0.8, ytickwidth = 0.8),
    Legend = (labelsize = 10,),
))

# Book color scheme
const NAVY   = colorant"#142D6E"
const SKY    = colorant"#7896D2"
const CORAL  = colorant"#E74C3C"
const TEAL   = colorant"#16A085"
const PURPLE = colorant"#8E44AD"
const ORANGE = colorant"#E67E22"

# Output paths
const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch14", "julia")

# -----------------------------------------------------------------------------
# 1D clamped biharmonic operator via polynomial trick
# -----------------------------------------------------------------------------

"""
    build_clamped_biharmonic(N)

Build the 1D clamped biharmonic operator D4 and interior second
derivative D2int on the Chebyshev grid with N+1 points.

The clamped boundary conditions u(+-1) = u'(+-1) = 0 are enforced
implicitly through the "polynomial trick" (Trefethen, Program 39):
represent u(x) = (1 - x^2) * p(x), so that u and u' vanish at x = +-1.

The fourth derivative of u = (1 - x^2) p is:
    u'''' = (1 - x^2) p'''' - 8x p''' - 12 p''

Returns (D4, D2int, x) where D4 is the (N-1)x(N-1) clamped biharmonic,
D2int is the (N-1)x(N-1) interior second derivative, and x is the full grid.
"""
function build_clamped_biharmonic(N)
    D, x = cheb_matrix(N)

    # Powers of D
    D2 = D * D
    D3 = D2 * D
    D4full = D3 * D

    # Diagonal weight matrices
    x2 = x.^2
    s = zeros(N + 1)
    s[2:N] = 1.0 ./ (1.0 .- x2[2:N])
    S = Diagonal(s)

    # Build D4 on full grid, then restrict to interior
    D4mat = (Diagonal(1.0 .- x2) * D4full - 8.0 * Diagonal(x) * D3 - 12.0 * D2) * S

    # Restrict to interior points (Julia 1-based: indices 2:N)
    D4 = D4mat[2:N, 2:N]
    D2int = D2[2:N, 2:N]

    return D4, D2int, x
end

# -----------------------------------------------------------------------------
# 2D biharmonic eigenvalue solver
# -----------------------------------------------------------------------------

"""
    solve_plate_eigenmodes(N; n_modes=25)

Solve the clamped plate eigenvalue problem on [-1,1]^2.

Assembles the 2D biharmonic operator via Kronecker products:
    L = I x D4 + D4 x I + 2 * (D2int x D2int)

and computes the eigenvalues/eigenvectors of L.

Returns (Lam, modes, x) where Lam contains normalized eigenvalues
sqrt(lambda_k / lambda_1), modes is a vector of (N+1)x(N+1) grids,
and x is the full Chebyshev grid.
"""
function solve_plate_eigenmodes(N; n_modes=25)
    D4, D2int, x = build_clamped_biharmonic(N)

    M = N - 1   # number of interior points per direction
    I_M = Matrix(1.0I, M, M)

    # 2D biharmonic via Kronecker products
    L = kron(I_M, D4) + kron(D4, I_M) + 2.0 * kron(D2int, D2int)

    # Solve eigenvalue problem
    F = eigen(L)
    eigenvalues = real.(F.values)
    eigenvectors = real.(F.vectors)

    # Sort by ascending eigenvalue (smallest first)
    idx = sortperm(eigenvalues)
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Keep only positive eigenvalues
    pos_mask = eigenvalues .> 0
    eigenvalues = eigenvalues[pos_mask]
    eigenvectors = eigenvectors[:, pos_mask]

    eigenvalues = eigenvalues[1:min(n_modes, length(eigenvalues))]
    eigenvectors = eigenvectors[:, 1:min(n_modes, size(eigenvectors, 2))]

    # Normalize: Lam_k = sqrt(lambda_k / lambda_1)
    Lam = sqrt.(eigenvalues ./ eigenvalues[1])

    # Reshape eigenvectors to 2D grids with zero-padded boundaries
    modes = Matrix{Float64}[]
    for k in 1:length(Lam)
        v = eigenvectors[:, k]
        V = reshape(v, (M, M))

        # Pad with zeros on all boundaries
        U = zeros(N + 1, N + 1)
        U[2:N, 2:N] = V
        push!(modes, U)
    end

    return Lam, modes, x
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    mkpath(OUTPUT_DIR)

    # Discretization parameter (following Trefethen, Program 39)
    N = 17
    n_modes = 25

    println("Eigenmodes of a Clamped Square Plate")
    println("="^55)
    @printf("  Chebyshev grid:     N = %d  (%d points per direction)\n", N, N + 1)
    @printf("  Interior points:    %d x %d = %d\n", N - 1, N - 1, (N - 1)^2)
    @printf("  Eigenvalue problem: %d x %d matrix\n", (N - 1)^2, (N - 1)^2)
    @printf("  Modes to compute:   %d\n\n", n_modes)

    # Solve
    Lam, modes, x = solve_plate_eigenmodes(N, n_modes=n_modes)

    # Print eigenvalue table
    println("  Normalized eigenvalues  sqrt(lambda_k / lambda_1):")
    println("  " * "-"^50)
    for k in 1:length(Lam)
        @printf("    Mode %2d:  %.6f\n", k, Lam[k])
    end
    println("  " * "-"^50)
    println()

    # =========================================================================
    # Figure: 5x5 contour grid of eigenmodes
    # =========================================================================
    nrows, ncols = 5, 5
    fig = Figure(size = (700, 750))

    XX = [x[j] for i in 1:N+1, j in 1:N+1]
    YY = [x[i] for i in 1:N+1, j in 1:N+1]

    for k in 1:min(n_modes, nrows * ncols)
        row = (k - 1) ÷ ncols + 1
        col = (k - 1) % ncols + 1

        ax = Axis(fig[row, col],
                  aspect = DataAspect(),
                  title = @sprintf("%.4f", Lam[k]),
                  titlesize = 9)

        U = modes[k]

        # Plot nodal lines (contour at level 0)
        contour!(ax, x, x, U', levels = [0.0], color = NAVY, linewidth = 0.8)

        # Draw domain boundary
        lines!(ax, [-1, 1, 1, -1, -1], [-1, -1, 1, 1, -1],
               color = :black, linewidth = 1.0)

        xlims!(ax, -1.05, 1.05)
        ylims!(ax, -1.05, 1.05)
        hidedecorations!(ax)
        hidespines!(ax)
    end

    # Main title
    Label(fig[0, :],
          L"Eigenmodes of a Clamped Plate: $\Delta^2 u = \lambda\, u$, $u = \partial_n u = 0$ on $\partial\Omega$",
          fontsize = 13)

    # Save figure
    save(joinpath(OUTPUT_DIR, "plate_eigenmodes.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "plate_eigenmodes.png"), fig, px_per_unit = 4)
    @printf("  Figure saved to: %s\n", joinpath(OUTPUT_DIR, "plate_eigenmodes.pdf"))

    println("="^55)
end

main()
