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
# The epsilon-pseudospectrum of the matrix pencil (A, B) is:
#
#     Lambda_eps = { z in C : sigma_min(A - z B) < eps }
#
# For normal operators the pseudospectra are simply eps-balls around the
# eigenvalues. The Orr-Sommerfeld operator is highly non-normal, so
# its pseudospectra extend far from the eigenvalues, revealing extreme
# sensitivity to perturbations.
#
# At R = 5772 (the critical Reynolds number for plane Poiseuille flow),
# all eigenvalues lie barely in the stable half-plane. The contour plot
# of log10(sigma_min) shows that eps-perturbations of size as small as
# 10^{-8} can push eigenvalues into the unstable half-plane, explaining
# the subcritical transition to turbulence.
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
# Orr-Sommerfeld matrices (boundary-bordering technique)
# -----------------------------------------------------------------------------
"""
    build_orr_sommerfeld_matrices(N, R, alpha)

Build the generalized eigenvalue matrices A, B for the Orr-Sommerfeld
equation for plane Poiseuille flow U(y) = 1 - y^2.

Boundary conditions: v(+-1) = v'(+-1) = 0 (no-slip).
Dirichlet conditions v(+-1) = 0 are enforced by restriction to interior
points (indices 2:N). Neumann conditions v'(+-1) = 0 are imposed by
replacing the first and last rows of both A and B.

Returns A, B matrices of size (N-1) x (N-1).
"""
function build_orr_sommerfeld_matrices(N::Int, R::Float64, alpha::Float64)
    D, x = cheb_matrix(N)

    # Second and fourth derivative matrices
    D2_full = D * D
    D4_full = D2_full * D2_full

    # Interior grid (remove boundary points x[1] = 1, x[N+1] = -1)
    x_int = x[2:N]

    # Base flow and its second derivative
    # U(y) = 1 - y^2,  U''(y) = -2
    U = Diagonal(1.0 .- x_int.^2)
    U2 = -2.0 * I(N - 1)    # diag(U'')

    # Interior blocks of D^2 and D^4
    D2_int = D2_full[2:N, 2:N]
    D4_int = D4_full[2:N, 2:N]

    # Shift operator: S = D^2 - alpha^2 I
    M_dim = N - 1
    I_int = Matrix(1.0I, M_dim, M_dim)
    S = D2_int - alpha^2 * I_int

    # S^2 = D^4 - 2 alpha^2 D^2 + alpha^4 I
    S2 = D4_int - 2.0 * alpha^2 * D2_int + alpha^4 * I_int

    # A = S^2 / (i alpha R) - U S + U''
    # B = -S
    A = S2 / (1im * alpha * R) - U * S + U2
    B = -copy(S)

    # Impose Neumann boundary conditions via row replacement:
    # v'(+1) = D[1, 2:N] @ v_int = 0  => row 1
    # v'(-1) = D[N+1, 2:N] @ v_int = 0  => row M_dim
    A[1, :] = D[1, 2:N]
    B[1, :] .= 0.0

    A[M_dim, :] = D[N+1, 2:N]
    B[M_dim, :] .= 0.0

    return A, B
end


# -----------------------------------------------------------------------------
# Eigenvalue computation
# -----------------------------------------------------------------------------
"""
    compute_eigenvalues(A, B)

Solve A v = lambda B v, filter spurious eigenvalues, sort by Im(lambda).
"""
function compute_eigenvalues(A, B)
    eigenvalues = eigvals(A, B)

    # Filter out infinite and NaN eigenvalues
    finite_mask = isfinite.(eigenvalues) .& (abs.(eigenvalues) .< 1e8)
    eigenvalues = eigenvalues[finite_mask]

    # Sort by imaginary part (most unstable first)
    idx = sortperm(imag.(eigenvalues), rev=true)
    eigenvalues = eigenvalues[idx]

    return eigenvalues
end


# -----------------------------------------------------------------------------
# Pseudospectra computation
# -----------------------------------------------------------------------------
"""
    compute_pseudospectra(A, B, x_re, x_im)

Compute log10(sigma_min(A - z B)) on a grid of complex numbers z.

For the generalized eigenvalue problem A v = lambda B v, the
epsilon-pseudospectrum is the set of z where sigma_min(A - z B) < eps.
"""
function compute_pseudospectra(A, B, x_re, x_im)
    nx = length(x_re)
    ny = length(x_im)

    sigma_grid = zeros(ny, nx)

    total = nx * ny
    count = 0

    for j in 1:ny
        for i in 1:nx
            z = x_re[i] + 1im * x_im[j]
            # Compute smallest singular value of (A - z B)
            sv = svdvals(A - z * B)
            sigma_grid[j, i] = log10(max(sv[end], 1e-20))
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
    @printf("  Wavenumber:           alpha = %.0f\n", alpha)
    @printf("  Chebyshev points:     N = %d\n", N)
    @printf("  Matrix size:          %d x %d\n", N - 1, N - 1)
    println()

    # ---- Build matrices and solve eigenvalue problem ----
    println("  Building Orr-Sommerfeld matrices ...")
    A, B = build_orr_sommerfeld_matrices(N, R, alpha)
    @printf("  Matrix shapes: A = %s, B = %s\n", string(size(A)), string(size(B)))

    println("  Computing eigenvalues ...")
    eigenvalues = compute_eigenvalues(A, B)
    @printf("  Number of finite eigenvalues: %d\n", length(eigenvalues))

    if length(eigenvalues) > 0
        # The most unstable eigenvalue has the largest Im(lambda)
        lam_crit = eigenvalues[1]
        @printf("  Most unstable eigenvalue:\n")
        @printf("    lambda = %.10f + %.10fi\n", real(lam_crit), imag(lam_crit))
        if imag(lam_crit) < 0
            println("    (stable: Im(lambda) < 0)")
        else
            println("    (UNSTABLE: Im(lambda) > 0)")
        end
    end

    # Print a few leading eigenvalues
    n_show = min(10, length(eigenvalues))
    @printf("\n  Leading %d eigenvalues (sorted by Im):\n", n_show)
    @printf("  %16s  %16s\n", "Re(lambda)", "Im(lambda)")
    println("  " * "-" ^ 36)
    for k in 1:n_show
        @printf("  %16.10f  %16.10f\n", real(eigenvalues[k]), imag(eigenvalues[k]))
    end
    println()

    # ---- Define grid for pseudospectra ----
    nx, ny = 120, 120
    x_re = range(-0.8, 0.2, length=nx)
    x_im = range(-1.0, 0.05, length=ny)

    @printf("  Computing pseudospectra on %d x %d grid ...\n", nx, ny)
    @printf("  Grid: Re in [%.1f, %.1f], Im in [%.2f, %.2f]\n",
            x_re[1], x_re[end], x_im[1], x_im[end])
    @printf("  This requires %d SVD computations of %dx%d matrices.\n",
            nx * ny, N - 1, N - 1)
    println()

    t_start = time()
    sigma_grid = compute_pseudospectra(A, B, x_re, x_im)
    t_elapsed = time() - t_start
    @printf("\n  Pseudospectra computed in %.1f seconds.\n\n", t_elapsed)

    # ---- Plot ----
    println("  Generating figure ...")

    fig = Figure(size = (700, 560))

    ax = Axis(fig[1, 1],
              xlabel = L"\mathrm{Re}(\lambda)",
              ylabel = L"\mathrm{Im}(\lambda)",
              title = @sprintf("\\varepsilon-Pseudospectra: R = %d, \\alpha = %g, N = %d",
                               Int(R), alpha, N))

    # Contour levels: log10(eps) = -1, -2, ..., -8
    levels = collect(-8.0:1.0:-1.0)

    # Filled contours
    cf = contourf!(ax, collect(x_re), collect(x_im), sigma_grid,
                   levels = levels, colormap = :YlOrRd, extendlow = :auto,
                   extendhigh = :auto)

    # Contour lines for clarity
    contour!(ax, collect(x_re), collect(x_im), sigma_grid,
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
             strokecolor = NAVY, strokewidth = 0.8)

    # Stability boundary: Im(lambda) = 0
    hlines!(ax, [0.0], color = :gray, linestyle = :dash, linewidth = 0.8, alpha = 0.5)

    xlims!(ax, x_re[1], x_re[end])
    ylims!(ax, x_im[1], x_im[end])

    # Save
    save(joinpath(OUTPUT_DIR, "pseudospectra.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "pseudospectra.png"), fig, px_per_unit = 4)
    println("  Figure saved to: $(joinpath(OUTPUT_DIR, "pseudospectra.pdf"))")

    println()
    println("=" ^ 60)
end

main()
