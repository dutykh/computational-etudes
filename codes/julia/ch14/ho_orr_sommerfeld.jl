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
# between walls at x = +-1), it reads:
#
#     (1/(i*alpha*R))(D^2 - alpha^2)^2 v - U(D^2 - alpha^2)v + U''v
#         = lambda [-(D^2 - alpha^2)] v
#
# with clamped (no-slip) boundary conditions v(+-1) = v'(+-1) = 0.
#
# This is a generalized eigenvalue problem A v = lambda B v where:
#     A = (D4 - 2*alpha^2*D2 + alpha^4*I)/(i*alpha*R) - U*(D2 - alpha^2*I) + U''
#     B = -(D2 - alpha^2*I)
#
# The fourth derivative matrix D4 is built via the polynomial trick
# (Trefethen, Program 38/40) to enforce clamped boundary conditions.
# The second derivative D2 uses Dirichlet conditions (boundary rows and
# columns removed).
#
# At R = 5772.22, the rightmost eigenvalue is approximately
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
    build_orr_sommerfeld_matrices(N, R, alpha)

Build the generalized eigenvalue matrices A, B for the Orr-Sommerfeld
equation for plane Poiseuille flow U(y) = 1 - y^2.

The clamped D4 is constructed via the polynomial trick (Trefethen, Program 38):
writing u(x) = (1 - x^2)*q(x) automatically enforces u(+-1) = u'(+-1) = 0.

Returns A, B matrices (interior points only) for A v = lambda B v.
"""
function build_orr_sommerfeld_matrices(N::Int, R::Float64, alpha::Float64)
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
    D4_int = D4mat[2:N, 2:N]

    # --- Dirichlet D2: simply strip boundary rows/columns ---
    D2_int = D2_full[2:N, 2:N]

    # Interior grid points
    x_int = x[2:N]

    # Base flow and its second derivative
    # U(y) = 1 - y^2,  U''(y) = -2
    U_diag = Diagonal(1.0 .- x_int.^2)
    U2_diag = -2.0 * I(N - 1)

    M = N - 1  # number of interior points
    I_int = Matrix(1.0I, M, M)

    # Shift operator: S_op = D^2 - alpha^2 I
    S_op = D2_int - alpha^2 * I_int

    # S^2 = D^4 - 2 alpha^2 D^2 + alpha^4 I
    S2 = D4_int - 2.0 * alpha^2 * D2_int + alpha^4 * I_int

    # A = S^2 / (i alpha R) - U S_op + U''
    # B = -S_op
    A = S2 / (1im * alpha * R) - U_diag * S_op + U2_diag

    # B = -S_op
    B = -copy(S_op)

    # Impose Neumann boundary conditions via row replacement:
    # v'(+1) = D[1, :] @ v_full = D[1, 2:N] @ v_int = 0  => row 1
    # v'(-1) = D[N+1, :] @ v_full = D[N+1, 2:N] @ v_int = 0  => row M
    A[1, :] = D[1, 2:N]
    B[1, :] .= 0.0

    A[M, :] = D[N+1, 2:N]
    B[M, :] .= 0.0

    return A, B, x_int
end


# -----------------------------------------------------------------------------
# Compute eigenvalues
# -----------------------------------------------------------------------------
"""
    compute_eigenvalues(A, B)

Solve the generalized eigenvalue problem A v = lambda B v.
Filters out spurious infinite eigenvalues and sorts by imaginary part.
"""
function compute_eigenvalues(A::Matrix{ComplexF64}, B::AbstractMatrix)
    eigenvalues = eigvals(A, B)

    # Filter out infinite and NaN eigenvalues (from singular rows of B)
    finite_mask = isfinite.(eigenvalues) .& (abs.(eigenvalues) .< 1e8)
    eigenvalues = eigenvalues[finite_mask]

    # Sort by imaginary part (most unstable first)
    idx = sortperm(imag.(eigenvalues), rev=true)
    eigenvalues = eigenvalues[idx]

    return eigenvalues
end


# -----------------------------------------------------------------------------
# Solve the Orr-Sommerfeld eigenvalue problem
# -----------------------------------------------------------------------------
"""
    orr_sommerfeld_eigenvalues(N; R=5772.0, alpha=1.0)

Compute the Orr-Sommerfeld eigenvalues for plane Poiseuille flow.

Returns eigenvalues and the rightmost eigenvalue (largest real part).
"""
function orr_sommerfeld_eigenvalues(N::Int; R::Float64=5772.0, alpha::Float64=1.0)
    A, B, x_int = build_orr_sommerfeld_matrices(N, R, alpha)
    eigenvalues = compute_eigenvalues(A, B)

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
    alpha = 1.0  # Streamwise wavenumber
    N_values = [40, 60, 80, 100]

    println("Orr-Sommerfeld Eigenvalue Problem")
    println("=" ^ 60)
    println("  Flow:             Plane Poiseuille  U(x) = 1 - x^2")
    @printf("  Reynolds number:  R = %.0f\n", R)
    @printf("  Wavenumber:       alpha = %.0f\n", alpha)
    println()

    # =========================================================================
    # Compute eigenvalues for each N
    # =========================================================================
    results = Dict{Int, Tuple{Vector{ComplexF64}, ComplexF64}}()
    for N in N_values
        eigenvalues, lam_max = orr_sommerfeld_eigenvalues(N, R=R, alpha=alpha)
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
                  xlabel = L"\mathrm{Re}\,c",
                  ylabel = L"\mathrm{Im}\,c",
                  title = @sprintf("N = %d    \\lambda_{max} = %.8f", N, real(lam_max)),
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
          L"Orr--Sommerfeld Spectrum ($R = 5772$, $\alpha = 1$)",
          fontsize = 13)

    # Save
    save(joinpath(OUTPUT_DIR, "orr_sommerfeld.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "orr_sommerfeld.png"), fig, px_per_unit = 4)
    println("\nFigure saved to: $(joinpath(OUTPUT_DIR, "orr_sommerfeld.pdf"))")
end

main()
