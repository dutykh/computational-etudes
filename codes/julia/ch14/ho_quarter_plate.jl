# ho_quarter_plate.jl
#
# Chapter 14: Higher-Order Boundary Value Problems
#
# Computational Etude 14.5: Quarter-plate symmetry reduction.
#
# Eigenmodes of a clamped plate on the quarter-square [0,1]^2 using
# symmetry reduction. Neumann (symmetry) conditions at x=0 and y=0
# select even-even modes, while clamped conditions are imposed at
# x=1 and y=1.
#
# The Chebyshev grid on [-1,1] is mapped to [0,1] via x_phys = (x+1)/2,
# scaling all derivative operators by a factor of 2.
#
# The 2D biharmonic operator is assembled via Kronecker products:
#
#     Delta^2 = d^4/dx^4 (x) I  +  2 d^2/dx^2 (x) d^2/dy^2  +  I (x) d^4/dy^4
#
# Boundary conditions are enforced by boundary bordering (row replacement):
# the rows of the 2D operator corresponding to boundary degrees of freedom
# are replaced by the appropriate BC equations, and the corresponding rows
# in the mass matrix are zeroed.
#
# For a fourth-order problem, two BCs per boundary segment are required:
#   - At x=0 (symmetry):  u_x = 0  and  u_xxx = 0  (odd derivatives vanish)
#   - At x=1 (clamped):   u = 0    and  u_x = 0
#   - At y=0 (symmetry):  u_y = 0  and  u_yyy = 0
#   - At y=1 (clamped):   u = 0    and  u_y = 0
#
# Generates Figure 14.5: First 12 even-even eigenmodes on the quarter plate.
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
# 2D biharmonic eigenvalue solver on the quarter plate
# -----------------------------------------------------------------------------
"""
    solve_quarter_plate(N, n_modes=12)

Solve the biharmonic eigenvalue problem on the quarter plate [0,1]^2
with symmetry (Neumann) BCs at x=0, y=0 and clamped BCs at x=1, y=1.

Returns normalized eigenvalues, eigenmodes, physical grid, and raw eigenvalues.
"""
function solve_quarter_plate(N::Int; n_modes::Int=12)
    # 1D Chebyshev operators scaled for [0,1]
    D, x_cheb = cheb_matrix(N)
    D_phys = 2.0 * D                  # scale for half-domain mapping
    D2 = D_phys * D_phys
    D3 = D2 * D_phys
    D4 = D2 * D2

    x_phys = (x_cheb .+ 1.0) ./ 2.0  # x_phys[1]=1, x_phys[N+1]=0

    # Index ordering (Julia 1-based):
    #   i = 1    -> x_phys = 1  (clamped boundary)
    #   i = N+1  -> x_phys = 0  (symmetry boundary)

    n = N + 1   # number of grid points per direction
    I1d = Matrix(1.0I, n, n)

    # Assemble 2D biharmonic on the full n^2 grid via Kronecker products
    L = kron(I1d, D4) + 2.0 * kron(D2, D2) + kron(D4, I1d)

    # Mass matrix for the generalized eigenvalue problem
    M = Matrix(1.0I, n * n, n * n)

    # Grid ordering: u[i + (j-1)*n] corresponds to point (x_phys[i], x_phys[j])

    # ----- Boundary bordering: replace rows with BC equations -----
    # We process BCs in four passes, skipping rows already assigned.
    bc_rows = Set{Int}()

    # Pass 1: Dirichlet u = 0 at the clamped boundaries (x=1 and y=1)
    for j in 1:n
        for i in 1:n
            k = i + (j - 1) * n

            if i == 1
                # x_phys = 1 (clamped): u(1, y) = 0
                L[k, :] .= 0.0
                L[k, k] = 1.0
                M[k, :] .= 0.0
                push!(bc_rows, k)
            elseif j == 1
                # y_phys = 1 (clamped): u(x, 1) = 0
                L[k, :] .= 0.0
                L[k, k] = 1.0
                M[k, :] .= 0.0
                push!(bc_rows, k)
            end
        end
    end

    # Pass 2: Normal derivative u_n = 0 at clamped boundaries
    for j in 1:n
        for i in 1:n
            k = i + (j - 1) * n
            k in bc_rows && continue

            if i == 2
                # du/dx = 0 at x_phys = 1 (row for i=2, next to clamped edge)
                L[k, :] .= 0.0
                for ii in 1:n
                    L[k, ii + (j - 1) * n] = D_phys[1, ii]
                end
                M[k, :] .= 0.0
                push!(bc_rows, k)
            elseif j == 2
                # du/dy = 0 at y_phys = 1 (row for j=2)
                L[k, :] .= 0.0
                for jj in 1:n
                    L[k, i + (jj - 1) * n] = D_phys[1, jj]
                end
                M[k, :] .= 0.0
                push!(bc_rows, k)
            end
        end
    end

    # Pass 3: Neumann u_n = 0 at symmetry boundaries (x=0 and y=0)
    for j in 1:n
        for i in 1:n
            k = i + (j - 1) * n
            k in bc_rows && continue

            if i == n
                # du/dx = 0 at x_phys = 0 (symmetry)
                L[k, :] .= 0.0
                for ii in 1:n
                    L[k, ii + (j - 1) * n] = D_phys[n, ii]
                end
                M[k, :] .= 0.0
                push!(bc_rows, k)
            elseif j == n
                # du/dy = 0 at y_phys = 0 (symmetry)
                L[k, :] .= 0.0
                for jj in 1:n
                    L[k, i + (jj - 1) * n] = D_phys[n, jj]
                end
                M[k, :] .= 0.0
                push!(bc_rows, k)
            end
        end
    end

    # Pass 4: Third derivative = 0 at symmetry boundaries
    # (for even functions, u_xxx(0) = 0 and u_yyy(0) = 0)
    for j in 1:n
        for i in 1:n
            k = i + (j - 1) * n
            k in bc_rows && continue

            if i == n - 1
                # d^3u/dx^3 = 0 at x_phys = 0 (symmetry)
                L[k, :] .= 0.0
                for ii in 1:n
                    L[k, ii + (j - 1) * n] = D3[n, ii]
                end
                M[k, :] .= 0.0
                push!(bc_rows, k)
            elseif j == n - 1
                # d^3u/dy^3 = 0 at y_phys = 0 (symmetry)
                L[k, :] .= 0.0
                for jj in 1:n
                    L[k, i + (jj - 1) * n] = D3[n, jj]
                end
                M[k, :] .= 0.0
                push!(bc_rows, k)
            end
        end
    end

    # ----- Solve the generalized eigenvalue problem L v = lambda M v -----
    eig_result = eigen(L, M)
    eigenvalues = eig_result.values
    eigenvectors = eig_result.vectors

    # Filter: keep only finite, real, positive eigenvalues
    finite_mask = isfinite.(eigenvalues)
    eigenvalues = eigenvalues[finite_mask]
    eigenvectors = eigenvectors[:, finite_mask]

    # Take real parts (self-adjoint problem)
    eigenvalues_real = real.(eigenvalues)

    # Discard spurious near-zero and negative eigenvalues
    pos_mask = eigenvalues_real .> 1.0
    eigenvalues_real = eigenvalues_real[pos_mask]
    eigenvectors = eigenvectors[:, pos_mask]

    # Sort by ascending eigenvalue
    idx = sortperm(eigenvalues_real)
    eigenvalues_real = eigenvalues_real[idx]
    eigenvectors = eigenvectors[:, idx]

    # Keep the requested number of modes
    n_avail = min(n_modes, length(eigenvalues_real))
    eigenvalues_real = eigenvalues_real[1:n_avail]
    eigenvectors = eigenvectors[:, 1:n_avail]

    # Normalize: Lam_k = sqrt(lambda_k / lambda_1)
    Lam = sqrt.(eigenvalues_real ./ eigenvalues_real[1])

    # Reshape eigenvectors to 2D grids
    modes = Matrix{Float64}[]
    for k in 1:n_avail
        v = real.(eigenvectors[:, k])
        U = reshape(v, (n, n))   # column-major: u[i,j] at (x_phys[i], x_phys[j])
        push!(modes, U)
    end

    return Lam, modes, x_phys, eigenvalues_real
end


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
"""
    plot_eigenmodes(Lam, modes, x_phys; n_modes=12)

Plot the first 12 even-even eigenmodes on the quarter plate [0,1]^2
in a 4x3 grid of subplots showing nodal lines (zero contours).
"""
function plot_eigenmodes(Lam, modes, x_phys; n_modes::Int=12)
    nrows, ncols = 4, 3
    fig = Figure(size = (700, 900))

    n_plot = min(n_modes, length(Lam))

    for k in 1:n_plot
        row = (k - 1) ÷ ncols + 1
        col = (k - 1) % ncols + 1

        ax = Axis(fig[row, col],
                  title = @sprintf("%.4f", Lam[k]),
                  titlesize = 9,
                  aspect = DataAspect(),
                  xticks = [0, 0.5, 1],
                  yticks = [0, 0.5, 1],
                  xticklabelsize = 7,
                  yticklabelsize = 7)

        U = modes[k]
        max_val = maximum(abs.(U))
        lvls = range(-max_val, max_val, length=21)

        # Filled contours for the mode shape
        contourf!(ax, x_phys, x_phys, U, levels = lvls, colormap = :RdBu, alpha = 0.4)

        # Nodal lines (contour at level 0)
        contour!(ax, x_phys, x_phys, U, levels = [0.0], color = NAVY, linewidth = 0.9)

        # Draw domain boundary
        lines!(ax, [0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color = :black, linewidth = 1.0)

        # Mark symmetry edges with dashed lines
        lines!(ax, [0, 0], [0, 1], color = TEAL, linewidth = 1.0, linestyle = :dash)
        lines!(ax, [0, 1], [0, 0], color = TEAL, linewidth = 1.0, linestyle = :dash)

        xlims!(ax, -0.05, 1.05)
        ylims!(ax, -0.05, 1.05)
    end

    Label(fig[0, :],
          L"Even-Even Eigenmodes of a Clamped Plate (Quarter Domain)",
          fontsize = 13)

    return fig
end


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
function main()
    mkpath(OUTPUT_DIR)

    # Discretization parameter
    N = 13
    n_modes = 12

    println("Quarter-Plate Symmetry Reduction: Even-Even Eigenmodes")
    println("=" ^ 60)
    @printf("  Domain:             [0, 1]^2  (quarter of [-1, 1]^2)\n")
    @printf("  Chebyshev grid:     N = %d  (%d points per direction)\n", N, N + 1)
    @printf("  System size:        %d x %d\n", (N + 1)^2, (N + 1)^2)
    @printf("  Modes to compute:   %d\n", n_modes)
    println("  BCs at x=0, y=0:   Neumann (symmetry)")
    println("  BCs at x=1, y=1:   Clamped (u = du/dn = 0)")
    println()

    # Solve
    Lam, modes, x_phys, eigenvalues = solve_quarter_plate(N, n_modes=n_modes)

    # Print eigenvalue table
    println("  Even-even eigenvalues of the clamped square plate:")
    println("  " * "-" ^ 55)
    @printf("  %4s  %14s  %18s\n", "Mode", "lambda_k", "sqrt(lam_k/lam_1)")
    println("  " * "-" ^ 55)
    for k in 1:length(Lam)
        @printf("  %4d  %14.6f  %18.6f\n", k, eigenvalues[k], Lam[k])
    end
    println("  " * "-" ^ 55)
    println()

    # Plot
    fig = plot_eigenmodes(Lam, modes, x_phys, n_modes=length(Lam))

    # Save
    save(joinpath(OUTPUT_DIR, "quarter_plate.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "quarter_plate.png"), fig, px_per_unit = 4)
    println("  Figure saved to: $(joinpath(OUTPUT_DIR, "quarter_plate.pdf"))")

    println("=" ^ 60)
end

main()
