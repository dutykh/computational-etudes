# ho_beam_eigenmodes.jl
#
# Chapter 14: Higher-Order Boundary Value Problems
#
# Computational Etude 14.2: Vibration Modes of a Clamped Beam.
#
# Solves the eigenvalue problem for a clamped-clamped Euler--Bernoulli beam:
#
#     u''''(x) = lambda * u(x),   -1 < x < 1,
#     u(+-1) = u'(+-1) = 0.
#
# The four boundary conditions (two Dirichlet and two Neumann) are imposed
# using the boundary-bordering technique:
#
#     1. Build the Chebyshev differentiation matrix D on N+1 points.
#     2. Compute D4 = D^4 (fourth-derivative matrix).
#     3. Extract the interior block D4_int = D4[2:N, 2:N] (Dirichlet u=0).
#     4. Replace the first and last rows of D4_int with the derivative
#        conditions D[1, 2:N] and D[N+1, 2:N] to enforce u'(+-1) = 0.
#     5. Solve the resulting generalised eigenvalue problem A v = lambda B v
#        where B has zero rows at the replaced positions.
#
# Reference eigenvalues are computed from the transcendental equations:
#     Symmetric modes:       tan(beta) + tanh(beta) = 0
#     Antisymmetric modes:   tan(beta) - tanh(beta) = 0
# In both cases lambda = beta^4.
#
# Parameters: N = 30 (default for figure generation).
#
# Generates Figure 14.2: First 6 eigenmodes and eigenvalue convergence.
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
# Reference eigenvalues from the transcendental equation
# -----------------------------------------------------------------------------

"""
    compute_reference_eigenvalues(n_modes=20)

Compute reference eigenvalues of u'''' = lambda u on [-1, 1]
with clamped BCs u(+-1) = u'(+-1) = 0.

The general solution u(x) = A cos(beta x) + B sin(beta x)
+ C cosh(beta x) + D sinh(beta x) with beta = lambda^{1/4}.
Applying the four BCs and separating into symmetric/antisymmetric
modes yields:

    Symmetric modes:       tan(beta) + tanh(beta) = 0
    Antisymmetric modes:   tan(beta) - tanh(beta) = 0

The eigenvalue is lambda = beta^4.
"""
function compute_reference_eigenvalues(n_modes=20)
    # Symmetric mode equation: tan(beta) + tanh(beta) = 0
    sym_eq(beta) = tan(beta) + tanh(beta)

    # Antisymmetric mode equation: tan(beta) - tanh(beta) = 0
    antisym_eq(beta) = tan(beta) - tanh(beta)

    betas = Float64[]
    eps_val = 0.01
    n_search = n_modes + 10

    for k in 0:n_search-1
        # Continuity interval for tan: (k*pi - pi/2, k*pi + pi/2)
        lo = k * pi - pi / 2.0 + eps_val
        hi = k * pi + pi / 2.0 - eps_val
        if lo < eps_val
            lo = eps_val
        end

        # Check for symmetric mode root via bisection
        try
            flo = sym_eq(lo)
            fhi = sym_eq(hi)
            if flo * fhi < 0
                root = bisection(sym_eq, lo, hi)
                if root > 0.5 && all(abs(root - b) > 0.01 for b in betas)
                    push!(betas, root)
                end
            end
        catch
        end

        # Check for antisymmetric mode root
        try
            flo = antisym_eq(lo)
            fhi = antisym_eq(hi)
            if flo * fhi < 0
                root = bisection(antisym_eq, lo, hi)
                if root > 0.5 && all(abs(root - b) > 0.01 for b in betas)
                    push!(betas, root)
                end
            end
        catch
        end
    end

    # Sort betas and compute eigenvalues lambda = beta^4
    sort!(betas)
    lambdas = betas[1:min(n_modes, length(betas))].^4

    return lambdas
end

"""
    bisection(f, a, b; tol=1e-15, maxiter=200)

Simple bisection root-finding method.
"""
function bisection(f, a, b; tol=1e-15, maxiter=200)
    fa = f(a)
    for _ in 1:maxiter
        c = (a + b) / 2.0
        fc = f(c)
        if abs(fc) < tol || (b - a) / 2.0 < tol
            return c
        end
        if fa * fc < 0
            b = c
        else
            a = c
            fa = fc
        end
    end
    return (a + b) / 2.0
end

# -----------------------------------------------------------------------------
# Biharmonic clamped matrix (boundary-bordering technique)
# -----------------------------------------------------------------------------

"""
    biharmonic_clamped_eigmatrix(N)

Build the interior fourth-derivative operator with clamped BCs
for the eigenvalue problem u''''(x) = lambda u(x).

Uses boundary-bordering:
    1. Extract interior block D4[2:N, 2:N] for Dirichlet BCs.
    2. Replace first row with D[1, 2:N] (enforces u'(1) = 0).
    3. Replace last row with D[N+1, 2:N] (enforces u'(-1) = 0).

Returns (A, B, x_int) for the generalised eigenvalue problem A v = lambda B v.
"""
function biharmonic_clamped_eigmatrix(N)
    D, x = cheb_matrix(N)

    # Fourth-derivative matrix
    D2 = D * D
    D4 = D2 * D2

    # Extract interior block (Dirichlet u(+-1) = 0)
    # Julia 1-based: x[1]=1, x[N+1]=-1, interior is 2:N
    A = copy(D4[2:N, 2:N])

    # Replace first row with derivative condition at x = +1 (index 1)
    # u'(+1) = D[1, :] * u = D[1, 2:N] * u_int = 0
    A[1, :] = D[1, 2:N]

    # Replace last row with derivative condition at x = -1 (index N+1)
    # u'(-1) = D[N+1, :] * u = D[N+1, 2:N] * u_int = 0
    A[end, :] = D[N+1, 2:N]

    # Mass matrix: identity except zero rows where BCs were imposed
    n_int = N - 1
    B = Matrix(1.0I, n_int, n_int)
    B[1, :] .= 0.0
    B[end, :] .= 0.0

    x_int = x[2:N]

    return A, B, x_int
end

# -----------------------------------------------------------------------------
# Eigenvalue solver
# -----------------------------------------------------------------------------

"""
    solve_beam_eigenmodes(N)

Solve the clamped beam eigenvalue problem u'''' = lambda u.

Returns sorted positive eigenvalues, corresponding eigenvectors
(at interior points), and the interior grid.
"""
function solve_beam_eigenmodes(N)
    A, B, x_int = biharmonic_clamped_eigmatrix(N)

    # Solve generalised eigenvalue problem A v = lambda B v
    F = eigen(A, B)
    eigenvalues = F.values
    eigenvectors = F.vectors

    # Filter out infinite/NaN eigenvalues (from the singular rows of B)
    finite_mask = isfinite.(eigenvalues) .& (abs.(eigenvalues) .< 1e15)
    eigenvalues = real.(eigenvalues[finite_mask])
    eigenvectors = real.(eigenvectors[:, finite_mask])

    # Sort by ascending real part
    idx = sortperm(eigenvalues)
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Keep only positive eigenvalues (physical modes)
    pos_mask = eigenvalues .> 0
    eigenvalues = eigenvalues[pos_mask]
    eigenvectors = eigenvectors[:, pos_mask]

    return eigenvalues, eigenvectors, x_int
end

# -----------------------------------------------------------------------------
# Barycentric interpolation for smooth plotting
# -----------------------------------------------------------------------------

"""
    cheb_interpolate(x_src, u_src, x_target)

Interpolate from Chebyshev nodes to a fine grid via barycentric formula.
"""
function cheb_interpolate(x_src, u_src, x_target)
    n = length(x_src)
    w = [(-1.0)^j for j in 0:n-1]
    w[1] *= 0.5
    w[end] *= 0.5

    u_target = similar(x_target)
    for (k, xt) in enumerate(x_target)
        idx = findfirst(i -> abs(x_src[i] - xt) < 1e-15, 1:n)
        if idx !== nothing
            u_target[k] = u_src[idx]
        else
            num = sum(w[j] / (xt - x_src[j]) * u_src[j] for j in 1:n)
            den = sum(w[j] / (xt - x_src[j]) for j in 1:n)
            u_target[k] = num / den
        end
    end
    return u_target
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    mkpath(OUTPUT_DIR)

    N = 30
    n_display = 6

    # Solve eigenvalue problem
    eigenvalues, eigenvectors, x_int = solve_beam_eigenmodes(N)

    # Full Chebyshev grid (for interpolation, including boundary zeros)
    _, x_full = cheb_matrix(N)

    # Reference eigenvalues
    n_ref = min(length(eigenvalues), 20)
    lam_exact = compute_reference_eigenvalues(n_ref)

    # Fine grid for smooth plotting
    x_fine = collect(range(-1, 1, length=500))

    # Colors for the 6 modes
    mode_colors = [NAVY, SKY, CORAL, TEAL, PURPLE, ORANGE]

    # =========================================================================
    # Figure: two-panel plot
    # =========================================================================
    fig = Figure(size = (900, 420))

    # ---- Panel 1: First 6 eigenmodes ----
    ax1 = Axis(fig[1, 1],
               xlabel = L"x",
               ylabel = L"u_n(x)",
               title = @sprintf("First %d eigenmodes (N = %d)", n_display, N))

    for k in 1:n_display
        if k > length(eigenvalues)
            break
        end

        # Extract eigenmode at interior points
        u_int = copy(eigenvectors[:, k])

        # Build full solution including boundary zeros
        u_full = zeros(N + 1)
        u_full[2:N] = u_int

        # Normalise so max|u| = 1
        u_full ./= maximum(abs.(u_full))

        # Fix sign convention: make the first peak (from the left) positive
        # x_full is descending (1 to -1), so reverse for left-to-right
        u_asc = reverse(u_full)
        first_extremum = u_asc[argmax(abs.(u_asc))]
        if first_extremum < 0
            u_full .*= -1
        end

        # Interpolate onto fine grid using barycentric formula
        u_fine = cheb_interpolate(x_full, u_full, x_fine)

        lam_k = eigenvalues[k]
        label = @sprintf("Mode %d (λ ≈ %.1f)", k, lam_k)

        lines!(ax1, x_fine, u_fine, color = mode_colors[k], linewidth = 1.5,
               label = label)
    end

    xlims!(ax1, -1.05, 1.05)
    axislegend(ax1, position = :lb, framevisible = true, labelsize = 7.5, nbanks = 1)
    hidespines!(ax1, :t, :r)
    hlines!(ax1, [0.0], color = :gray, linewidth = 0.5)

    # ---- Panel 2: Eigenvalue convergence ----
    ax2 = Axis(fig[1, 2],
               xlabel = "Mode number n",
               ylabel = L"|\lambda_n^{\mathrm{num}} - \lambda_n^{\mathrm{exact}}| / \lambda_n^{\mathrm{exact}}",
               title = @sprintf("Relative eigenvalue error (N = %d)", N),
               yscale = log10)

    n_compare = min(length(eigenvalues), length(lam_exact))
    modes = collect(1:n_compare)
    rel_errors = [max(abs(eigenvalues[k] - lam_exact[k]) / lam_exact[k], 1e-16)
                  for k in 1:n_compare]

    scatterlines!(ax2, modes, rel_errors, color = NAVY, linewidth = 1.2,
                  markersize = 6, markercolor = NAVY,
                  strokecolor = :white, strokewidth = 0.5)

    hidespines!(ax2, :t, :r)

    # Main title
    Label(fig[0, :],
          L"Clamped Beam Eigenmodes: $u'''' = \lambda\, u$, $u(\pm 1) = u'(\pm 1) = 0$",
          fontsize = 13)

    # Save figure
    save(joinpath(OUTPUT_DIR, "beam_eigenmodes.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "beam_eigenmodes.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "beam_eigenmodes.pdf"))")

    # =========================================================================
    # Print summary
    # =========================================================================
    println("\nVibration Modes of a Clamped Beam")
    println("="^60)
    @printf("  Equation: u''''(x) = lambda u(x),  x in (-1, 1)\n")
    @printf("  BCs: u(+-1) = u'(+-1) = 0\n")
    @printf("  N = %d Chebyshev points\n\n", N)

    println("  Reference eigenvalues (from transcendental equation):")
    println("  " * "-"^50)
    @printf("  %6s %20s %14s\n", "Mode", "lambda (exact)", "beta")
    println("  " * "-"^50)
    for k in 1:min(10, length(lam_exact))
        beta_k = lam_exact[k]^0.25
        @printf("  %6d %20.6f %14.6f\n", k, lam_exact[k], beta_k)
    end
    println("  " * "-"^50)
    println()

    println("  Comparison of numerical vs exact eigenvalues:")
    println("  " * "-"^70)
    @printf("  %6s %18s %18s %14s\n", "Mode", "lambda (num)", "lambda (exact)", "Rel. Error")
    println("  " * "-"^70)
    for k in 1:n_compare
        rel_err = abs(eigenvalues[k] - lam_exact[k]) / lam_exact[k]
        @printf("  %6d %18.6f %18.6f %14.4e\n", k, eigenvalues[k], lam_exact[k], rel_err)
    end
    println("  " * "-"^70)
end

main()
