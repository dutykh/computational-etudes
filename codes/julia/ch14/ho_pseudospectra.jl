# ho_pseudospectra.jl
#
# Chapter 14: Higher-Order Boundary Value Problems
#
# Computational Etude 14.8: Pseudospectra of the Orr-Sommerfeld Operator.
#
# Computes the epsilon-pseudospectra of the Orr-Sommerfeld operator for
# plane Poiseuille flow at the critical Reynolds number R = 5772,
# streamwise wavenumber alpha = 1, using N = 100 Chebyshev points.
#
# The Orr-Sommerfeld equation governs the linear stability of parallel
# shear flows. For the base flow U(y) = 1 - y^2 on [-1, 1], the
# generalized eigenvalue problem for the complex phase speed reads
# A v = c B v, where A and B encode the fourth-order ODE with
# no-slip boundary conditions v(+-1) = v'(+-1) = 0.
#
# The Dirichlet conditions are handled by restriction to interior
# Chebyshev points. The derivative (Neumann) conditions v'(+-1) = 0
# are imposed by boundary bordering. The resulting two algebraic
# constraints are then eliminated analytically via a projection that
# expresses the boundary unknowns in terms of the interior unknowns,
# yielding a non-singular reduced system whose matrix C = B_red^{-1} A_red
# is the one whose pseudospectra we compute.
#
# The epsilon-pseudospectrum is defined as
#
#     Lambda_eps(C) = { z in C : sigma_min(C - z I) < eps }
#
# where sigma_min denotes the smallest singular value. For normal
# operators the pseudospectra are simply eps-balls around the
# eigenvalues. The Orr-Sommerfeld operator is highly non-normal, so
# its pseudospectra extend far from the eigenvalues, revealing extreme
# sensitivity to perturbations.
#
# At R = 5772 and alpha = 1, just below the critical point of plane
# Poiseuille flow, all eigenvalues lie in the stable half-plane Im(c) < 0.
#
# Reference: Trefethen, "Spectral Methods in MATLAB", Program 40.
#            Trefethen & Embree, "Spectra and Pseudospectra" (2005).
#            Reddy, Schmid & Henningson, SIAM J. Appl. Math. 53 (1993).
#
# Generates Figure 14.7: Pseudospectra contour plot with eigenvalues.
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
# Build Orr-Sommerfeld matrix C = B^{-1} A via constraint elimination
# -----------------------------------------------------------------------------
"""
    build_orr_sommerfeld_matrix(N, R, alpha)

Build the matrix C = B^{-1} A for the Orr-Sommerfeld eigenvalue problem
for plane Poiseuille flow U(y) = 1 - y^2.

The Orr-Sommerfeld equation, written as a generalized eigenvalue problem
A v = c B v for the phase speed c, is

    A = (D^2 - alpha^2)^2 / (i alpha R) - diag(U) (D^2 - alpha^2) + diag(U'')
    B = -(D^2 - alpha^2)

restricted to interior Chebyshev points (Dirichlet conditions v(+-1) = 0).
The Neumann conditions v'(+-1) = 0 replace the first and last rows, which
makes B singular. These two algebraic constraints are then eliminated:
they express the unknowns at the interior points nearest +-1 as linear
combinations of the remaining unknowns, and substituting them into the
dynamical rows yields a reduced system with non-singular B_red, so that
C = B_red^{-1} A_red.

Returns C, of size (N-3) x (N-3), and its eigenvalues sorted by imaginary
part (descending).
"""
function build_orr_sommerfeld_matrix(N::Int, R::Float64, alpha::Float64)
    D, x = cheb_matrix(N)

    # Second and fourth derivative matrices
    D2 = D * D
    D4 = D2 * D2

    # Number of interior points
    n = N - 1

    # Interior grid (remove boundary points x[1] = 1, x[N+1] = -1)
    x_int = x[2:N]

    # Base flow: U(y) = 1 - y^2,  U''(y) = -2
    U = Diagonal(1.0 .- x_int.^2)
    U2 = -2.0 * Matrix(1.0I, n, n)

    # Interior blocks of D^2 and D^4
    D2_int = D2[2:N, 2:N]
    D4_int = D4[2:N, 2:N]

    # Shift operator: S = D^2 - alpha^2 I
    I_int = Matrix(1.0I, n, n)
    S = D2_int - alpha^2 * I_int

    # S^2 = D^4 - 2 alpha^2 D^2 + alpha^4 I
    S2 = D4_int - 2.0 * alpha^2 * D2_int + alpha^4 * I_int

    # Build A and B (boundary-bordering form)
    A = S2 / (1im * alpha * R) - U * S + U2
    B = complex.(-S)

    # Replace first/last rows with Neumann conditions v'(+-1) = 0
    A[1, :] = D[1, 2:N]
    B[1, :] .= 0.0
    A[n, :] = D[N+1, 2:N]
    B[n, :] .= 0.0

    # --- Constraint elimination ---
    # Rows 1 and n are algebraic constraints A[1,:] v = 0, A[n,:] v = 0
    # (B vanishes there, so the eigenvalue does not appear). Partition
    # v = [v_1; v_free; v_n] with v_free = v[2:n-1] and solve the two
    # constraints for [v_1; v_n] = P v_free, where P = -CC^{-1} F.
    CC = [A[1, 1] A[1, n];
          A[n, 1] A[n, n]]
    F = vcat(transpose(A[1, 2:n-1]), transpose(A[n, 2:n-1]))
    P = -(CC \ F)                    # size (2, n-2)

    # Lifting matrix T: v = T v_free, size (n, n-2)
    T = zeros(ComplexF64, n, n - 2)
    T[1, :] = P[1, :]                # v_1 row
    T[2:n-1, :] = Matrix(1.0I, n - 2, n - 2)   # v_free rows
    T[n, :] = P[2, :]                # v_n row

    # Reduced system (drop the two constraint rows)
    A_red = A[2:n-1, :] * T          # size (n-2, n-2)
    B_red = B[2:n-1, :] * T          # size (n-2, n-2)

    # C = B_red^{-1} A_red
    C = B_red \ A_red

    # Eigenvalues of C, sorted by imaginary part (most unstable first)
    eigenvalues = eigvals(C)
    eigenvalues = eigenvalues[sortperm(imag.(eigenvalues), rev=true)]

    return C, eigenvalues
end


# -----------------------------------------------------------------------------
# Pseudospectra computation
# -----------------------------------------------------------------------------
"""
    compute_pseudospectra(C, x_re, x_im)

Compute log10(sigma_min(C - z I)) on a grid of complex numbers z. The level
set where this equals log10(eps) is the boundary of the eps-pseudospectrum.
"""
function compute_pseudospectra(C, x_re, x_im)
    nx = length(x_re)
    ny = length(x_im)
    M = size(C, 1)
    Id = Matrix{ComplexF64}(I, M, M)

    sigma_grid = zeros(ny, nx)

    total = nx * ny
    count = 0

    for j in 1:ny
        for i in 1:nx
            z = x_re[i] + 1im * x_im[j]
            # Compute smallest singular value of (C - z I)
            sv = svdvals(C - z * Id)
            sigma_grid[j, i] = log10(max(minimum(sv), 1e-20))
            count += 1
        end

        # Progress report every 10 rows
        if j % 10 == 0
            pct = 100.0 * count / total
            @printf("    ... %.0f%% complete (%d/%d grid points)\n", pct, count, total)
        end
    end

    return sigma_grid
end


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
function main()
    mkpath(OUTPUT_DIR)

    # Physical parameters
    R = 5772.0         # Critical Reynolds number for plane Poiseuille flow
    alpha = 1.0        # Streamwise wavenumber
    N = 100            # Number of Chebyshev intervals

    println("Pseudospectra of the Orr-Sommerfeld Operator")
    println("=" ^ 60)
    @printf("  Reynolds number:      R = %d\n", Int(R))
    @printf("  Wavenumber:           alpha = %.1f\n", alpha)
    @printf("  Chebyshev points:     N = %d\n", N)
    @printf("  Interior points:      %d\n", N - 1)
    @printf("  Matrix C size:        %d x %d (after constraint elimination)\n", N - 3, N - 3)
    println()

    # ---- Build the matrix C and its eigenvalues ----
    println("  Building Orr-Sommerfeld matrix C = B_red^{-1} A_red ...")
    C, eigenvalues = build_orr_sommerfeld_matrix(N, R, alpha)
    @printf("  Matrix C shape: %s\n", string(size(C)))
    @printf("  Number of eigenvalues: %d\n", length(eigenvalues))

    if length(eigenvalues) > 0
        lam_crit = eigenvalues[1]
        println("  Most unstable eigenvalue:")
        @printf("    c = %.10f + %.10f i\n", real(lam_crit), imag(lam_crit))
        if imag(lam_crit) < 0
            println("    (stable: Im(c) < 0)")
        else
            println("    (UNSTABLE: Im(c) > 0)")
        end
    end

    # Print a few leading eigenvalues
    n_show = min(10, length(eigenvalues))
    @printf("\n  Leading %d eigenvalues (sorted by Im):\n", n_show)
    @printf("  %16s  %16s\n", "Re(c)", "Im(c)")
    println("  " * "-" ^ 36)
    for k in 1:n_show
        @printf("  %16.10f  %16.10f\n", real(eigenvalues[k]), imag(eigenvalues[k]))
    end
    println()

    # ---- Define grid for pseudospectra ----
    # The eigenvalues lie in Re ~ [0.2, 1.0], Im ~ [-1.0, 0],
    # so we center the grid on this region with some padding.
    nx, ny = 120, 120
    x_re = collect(range(-0.1, 1.1, length=nx))
    x_im = collect(range(-1.0, 0.1, length=ny))

    @printf("  Computing pseudospectra on %d x %d grid ...\n", nx, ny)
    @printf("  Grid: Re in [%.1f, %.1f], Im in [%.1f, %.1f]\n",
            x_re[1], x_re[end], x_im[1], x_im[end])
    @printf("  This requires %d SVD computations of %dx%d matrices.\n",
            nx * ny, size(C, 1), size(C, 1))
    println()

    t_start = time()
    sigma_grid = compute_pseudospectra(C, x_re, x_im)
    t_elapsed = time() - t_start
    @printf("\n  Pseudospectra computed in %.1f seconds.\n\n", t_elapsed)

    # ---- Plot ----
    println("  Generating figure ...")

    fig = Figure(size = (640, 520))

    ax = Axis(fig[1, 1],
              xlabel = L"\mathrm{Re}\,c",
              ylabel = L"\mathrm{Im}\,c",
              title = latexstring(@sprintf("\\varepsilon\\text{‐pseudospectra of Orr–Sommerfeld}{:}\\ R = %d,\\ \\alpha = %g,\\ N = %d",
                                           Int(R), alpha, N)))

    # Contour levels: log10(eps) = -8, -7, ..., -1. Small sigma_min (deep
    # inside the pseudospectra) is dark red, large sigma_min is pale yellow.
    levels = collect(-8.0:1.0:-1.0)
    cmap = Reverse(:YlOrRd)

    # Filled contours (the matrix is indexed [Im, Re], hence the transpose)
    cf = contourf!(ax, x_re, x_im, permutedims(sigma_grid),
                   levels = levels, colormap = cmap,
                   extendlow = :auto, extendhigh = :auto)

    # Contour lines for clarity
    contour!(ax, x_re, x_im, permutedims(sigma_grid),
             levels = levels, color = :black, linewidth = 0.5, alpha = 0.6)

    Colorbar(fig[1, 2], cf,
             label = L"\log_{10}\,\sigma_{\min}(C - z\,I)",
             ticks = levels)

    # Overlay eigenvalues (only those within the plot region)
    in_region = (real.(eigenvalues) .>= x_re[1]) .& (real.(eigenvalues) .<= x_re[end]) .&
                (imag.(eigenvalues) .>= x_im[1]) .& (imag.(eigenvalues) .<= x_im[end])
    eigs_plot = eigenvalues[in_region]

    scatter!(ax, real.(eigs_plot), imag.(eigs_plot),
             color = :white, markersize = 5,
             strokecolor = NAVY, strokewidth = 0.8, label = "Eigenvalues")

    # Stability boundary: Im(c) = 0
    hlines!(ax, [0.0], color = :gray, linestyle = :dash, linewidth = 0.8, alpha = 0.5)

    axislegend(ax, position = :lt, labelsize = 9, framevisible = true)

    xlims!(ax, x_re[1], x_re[end])
    ylims!(ax, x_im[1], x_im[end])

    # Save
    save(joinpath(OUTPUT_DIR, "pseudospectra.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "pseudospectra.png"), fig, px_per_unit = 4)
    println("  Figure saved to: $(joinpath(OUTPUT_DIR, "pseudospectra.pdf"))")

    println()
    println("=" ^ 60)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
