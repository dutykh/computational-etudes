# ho_orr_sommerfeld.jl
#
# Chapter 14: Higher-Order Boundary Value Problems
#
# Computational Etude 14.7: Orr-Sommerfeld eigenvalue problem for plane
# Poiseuille flow at the critical Reynolds number R = 5772, wavenumber
# alpha = 1. The eigenvalues form the famous Y-shaped spectrum in the
# complex plane.
#
# The Orr-Sommerfeld equation governs linear stability of parallel shear
# flows. For plane Poiseuille flow U(x) = 1 - x^2 (parabolic profile
# between walls at x = +-1) and alpha = 1, a normal mode proportional to
# exp(lambda t) satisfies
#
#     (1/R)(u'''' - 2u'' + u) - 2i u - i(1 - x^2)(u'' - u) = lambda (u'' - u)
#
# with clamped (no-slip) boundary conditions u(+-1) = u'(+-1) = 0. The
# eigenvalue is lambda = -i*alpha*c, where c is the complex phase speed, so
# that Re(lambda) = alpha*Im(c) is the growth rate of the mode.
#
# This is a generalized eigenvalue problem A v = lambda B v where:
#     A = (1/R)(D4 - 2*D2 + I) - 2i*I - i*diag(1 - x^2)*(D2 - I)
#     B = D2 - I
#
# The fourth derivative matrix D4 is built via the polynomial trick
# (Trefethen, Program 38/40) to enforce clamped boundary conditions.
# The second derivative D2 uses Dirichlet conditions (boundary rows and
# columns removed).
#
# At R = 5772, the rightmost eigenvalue is approximately
# lambda_max ~ -0.00007819..., placing the flow just barely stable.
#
# Based on Program 40 from Trefethen's "Spectral Methods in MATLAB".
#
# Generates Figure 14.6: 2x2 grid of Orr-Sommerfeld spectra for
# N = 40, 60, 80, 100.
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
# Build clamped fourth-derivative and Dirichlet second-derivative matrices
# -----------------------------------------------------------------------------
"""
    build_orr_sommerfeld_matrices(N)

Build the interior D4 (clamped) and D2 (Dirichlet) matrices needed for the
Orr-Sommerfeld generalized eigenvalue problem.

The clamped fourth-derivative operator is constructed via the polynomial
trick from Trefethen (Program 38): writing u(x) = (1 - x^2)*q(x)
automatically enforces u(+-1) = u'(+-1) = 0. The fourth derivative then
becomes

    D4 = [diag(1 - x^2)*D^4 - 8*diag(x)*D^3 - 12*D^2] * S

where S = diag(1/(1 - x_j^2)) on interior points (zero at boundaries).
The second derivative D2 uses simple Dirichlet conditions: boundary rows
and columns are stripped from D^2.

Returns D4 and D2, both (N-1) x (N-1), and the interior points x_int.
"""
function build_orr_sommerfeld_matrices(N::Int)
    D, x = cheb_matrix(N)

    # Powers of the differentiation matrix
    D2_full = D * D
    D3_full = D2_full * D
    D4_full = D3_full * D

    # --- Clamped D4 via polynomial trick ---
    # S maps v -> q on interior points, zeros at boundaries
    s = zeros(N + 1)
    s[2:N] .= 1.0 ./ (1.0 .- x[2:N].^2)
    S = Diagonal(s)

    # Full fourth-order operator with the polynomial trick
    D4mat = (Diagonal(1.0 .- x.^2) * D4_full
             - 8.0 * Diagonal(x) * D3_full
             - 12.0 * D2_full) * S

    # Restrict to interior points (indices 2:N in Julia 1-based)
    D4 = D4mat[2:N, 2:N]

    # --- Dirichlet D2: simply strip boundary rows/columns ---
    D2 = D2_full[2:N, 2:N]

    # Interior grid points
    x_int = x[2:N]

    return D4, D2, x_int
end


# -----------------------------------------------------------------------------
# Solve the Orr-Sommerfeld eigenvalue problem
# -----------------------------------------------------------------------------
"""
    orr_sommerfeld_eigenvalues(N; R=5772.0)

Compute the Orr-Sommerfeld eigenvalues for plane Poiseuille flow with
alpha = 1 by solving the generalized eigenvalue problem A v = lambda B v:

    A = (1/R)(D4 - 2*D2 + I) - 2i*I - i*diag(1 - x^2)*(D2 - I)
    B = D2 - I

Returns the eigenvalues and the rightmost eigenvalue (largest real part,
that is the largest growth rate).
"""
function orr_sommerfeld_eigenvalues(N::Int; R::Float64=5772.0)
    D4, D2, x_int = build_orr_sommerfeld_matrices(N)

    M = N - 1  # number of interior points
    I_int = Matrix(1.0I, M, M)

    # Assemble the A and B matrices
    A = ((D4 - 2.0 * D2 + I_int) / R
         - 2im * I_int
         - 1im * Diagonal(1.0 .- x_int.^2) * (D2 - I_int))
    B = D2 - I_int

    # Solve generalized eigenvalue problem
    eigenvalues = eigvals(A, B)

    # Find the rightmost eigenvalue (largest real part)
    idx_max = argmax(real.(eigenvalues))
    lam_max = eigenvalues[idx_max]

    return eigenvalues, lam_max
end


# -----------------------------------------------------------------------------
# Main routine
# -----------------------------------------------------------------------------
function main()
    mkpath(OUTPUT_DIR)

    R = 5772.0   # Critical Reynolds number for plane Poiseuille flow
    N_values = [40, 60, 80, 100]

    println("Orr-Sommerfeld Eigenvalue Problem")
    println("=" ^ 60)
    println("  Flow:             Plane Poiseuille  U(x) = 1 - x^2")
    @printf("  Reynolds number:  R = %.0f\n", R)
    println("  Wavenumber:       alpha = 1")
    println()

    # =========================================================================
    # Compute eigenvalues for each N
    # =========================================================================
    results = Dict{Int, Tuple{Vector{ComplexF64}, ComplexF64}}()
    for N in N_values
        eigenvalues, lam_max = orr_sommerfeld_eigenvalues(N, R=R)
        results[N] = (eigenvalues, lam_max)
        @printf("  N = %3d:  lambda_max = %+.10f %+.10fi\n",
                N, real(lam_max), imag(lam_max))
    end

    println()
    println("  Expected:  lambda_max ~ -0.00007819...")
    println("=" ^ 60)

    # =========================================================================
    # Figure: 2x2 subplot grid
    # =========================================================================
    fig = Figure(size = (800, 650))

    for (k, N) in enumerate(N_values)
        row = (k - 1) ÷ 2 + 1
        col = (k - 1) % 2 + 1

        eigenvalues, lam_max = results[N]

        ax = Axis(fig[row, col],
                  xlabel = L"\mathrm{Re}\,\lambda",
                  ylabel = L"\mathrm{Im}\,\lambda",
                  title = latexstring(@sprintf("N = %d,\\quad \\mathrm{Re}\\,\\lambda_{\\max} = %.11f", N, real(lam_max))),
                  titlesize = 11)

        # Plot eigenvalues as dots in the complex plane
        scatter!(ax, real.(eigenvalues), imag.(eigenvalues),
                 color = NAVY, markersize = 5)

        # Axis limits
        xlims!(ax, -0.8, 0.2)
        ylims!(ax, -1.0, 0.0)

        # Stability boundary: Re(lambda) = 0
        vlines!(ax, [0.0], color = :gray, linestyle = :dot, linewidth = 0.6,
                alpha = 0.5)

        hidespines!(ax, :t, :r)
    end

    Label(fig[0, :],
          L"Orr–Sommerfeld Spectrum ($R = 5772$, $\alpha = 1$)",
          fontsize = 13)

    # Save
    save(joinpath(OUTPUT_DIR, "orr_sommerfeld.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "orr_sommerfeld.png"), fig, px_per_unit = 4)
    println("\nFigure saved to: $(joinpath(OUTPUT_DIR, "orr_sommerfeld.pdf"))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
