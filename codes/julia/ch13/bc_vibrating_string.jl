# bc_vibrating_string.jl
#
# Chapter 13: Boundary Conditions in Spectral Methods
#
# Vibrating string eigenvalue problem with mixed boundary conditions.
#
#     -u_xx = lambda * u   on [0, 1]
#     u(0) = 0             (Dirichlet at left end)
#     u'(1) = 0            (Neumann / free end at right)
#
# The domain [0, 1] is mapped to the Chebyshev domain [-1, 1] via x = (xi + 1)/2,
# so D_x = 2 * D_xi.
#
# Chebyshev ordering: xi[1] = 1 -> physical x = 1 (free end, Neumann)
#                     xi[N+1] = -1 -> physical x = 0 (fixed end, Dirichlet)
#
# The operator A = -D_phys^2 with boundary rows replaced:
#     Row 1:   D_phys[1,:] (Neumann condition u'(1) = 0)
#     Row N+1: identity row (Dirichlet condition u(0) = 0)
#
# Exact eigenvalues: lambda_n = ((2n-1)*pi/2)^2,  n = 1, 2, ...
# Exact eigenfunctions: u_n(x) = sin((2n-1)*pi*x/2)
#
# Generates Figure 13.10: First 6 eigenmodes and eigenvalue errors.
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
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch13", "julia")

# -----------------------------------------------------------------------------
# Eigenvalue problem solver
# -----------------------------------------------------------------------------
"""
    solve_vibrating_string(N)

Solve the vibrating string eigenvalue problem:
    -u_xx = lambda * u on [0,1], u(0) = 0, u'(1) = 0

Returns sorted eigenvalues and eigenvectors on the physical grid x in [0, 1].
"""
function solve_vibrating_string(N)
    D, xi = cheb_matrix(N)
    n = N + 1

    # Map to physical domain [0, 1]: x = (xi + 1)/2
    # D_x = 2 * D_xi
    D_phys = 2.0 * D
    D2_phys = D_phys * D_phys
    x_phys = (xi .+ 1.0) ./ 2.0  # x_phys[1]=1 (free), x_phys[N+1]=0 (fixed)

    # Build operator A = -D2_phys
    A = -D2_phys

    # Row 1 (xi=1 -> x=1): Neumann u'(1) = 0
    A[1, :] = D_phys[1, :]

    # Row N+1 (xi=-1 -> x=0): Dirichlet u(0) = 0
    A[n, :] .= 0.0
    A[n, n] = 1.0

    # Solve eigenvalue problem
    F = eigen(A)
    lambdas = real.(F.values)
    vecs = real.(F.vectors)

    # Sort by eigenvalue magnitude
    idx = sortperm(lambdas)
    lambdas_sorted = lambdas[idx]
    vecs_sorted = vecs[:, idx]

    return lambdas_sorted, vecs_sorted, x_phys
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    mkpath(OUTPUT_DIR)

    N = 36
    lambdas, vecs, x_phys = solve_vibrating_string(N)

    # Exact eigenvalues: lambda_n = ((2n-1)*pi/2)^2
    n_exact = 20
    exact_lambdas = [((2*n - 1) * pi / 2)^2 for n in 1:n_exact]

    # Fine grid for exact eigenfunctions
    x_fine = collect(range(0, 1, length=200))

    # =========================================================================
    # Figure: Two panels
    # =========================================================================
    fig = Figure(size = (900, 420))

    # ---- Panel 1: First 6 eigenmodes ----
    ax1 = Axis(fig[1, 1],
               xlabel = L"x",
               ylabel = L"u_n(x)",
               title = "First 6 eigenmodes")

    colors_modes = [NAVY, TEAL, CORAL, PURPLE, ORANGE, SKY]
    n_modes = 6

    # Skip spurious eigenvalues near zero (from boundary rows)
    # Find the first positive eigenvalue that matches exact
    # The first few sorted eigenvalues include 0s from boundary conditions
    # Filter: keep eigenvalues > 1.0 (first exact is ~2.47)
    valid_mask = lambdas .> 1.0
    valid_idx = findall(valid_mask)

    for mode in 1:n_modes
        if mode > length(valid_idx)
            break
        end
        idx = valid_idx[mode]
        lam_computed = lambdas[idx]
        lam_exact = exact_lambdas[mode]

        # Computed eigenfunction
        psi = vecs[:, idx]
        # Normalize: max |psi| = 1
        psi ./= maximum(abs.(psi))
        # Ensure consistent sign (positive at first interior max)
        if psi[N ÷ 2] < 0
            psi .*= -1
        end

        # Exact eigenfunction
        u_exact = sin.((2 * mode - 1) .* pi .* x_fine ./ 2.0)
        u_exact ./= maximum(abs.(u_exact))

        # Sort physical grid for plotting
        idx_sort = sortperm(x_phys)

        # Plot exact (thin dashed) and computed (solid with markers)
        lines!(ax1, x_fine, u_exact, color = (colors_modes[mode], 0.4),
               linewidth = 1.0, linestyle = :dash)
        scatterlines!(ax1, x_phys[idx_sort], psi[idx_sort],
                       color = colors_modes[mode], linewidth = 1.5,
                       markersize = 3, markercolor = colors_modes[mode],
                       label = @sprintf("n=%d (λ=%.2f)", mode, lam_computed))
    end

    xlims!(ax1, -0.02, 1.02)
    axislegend(ax1, position = :lb, framevisible = true, labelsize = 8, nbanks = 2)
    hidespines!(ax1, :t, :r)

    # ---- Panel 2: Eigenvalue errors ----
    ax2 = Axis(fig[1, 2],
               xlabel = "Mode number n",
               ylabel = L"|\lambda_n - \lambda_n^{\mathrm{exact}}|",
               title = @sprintf("Eigenvalue errors (N = %d)", N),
               yscale = log10)

    n_compare = min(length(valid_idx), n_exact)
    errors = Float64[]
    mode_numbers = Int[]

    for mode in 1:n_compare
        idx = valid_idx[mode]
        err = abs(lambdas[idx] - exact_lambdas[mode])
        err = max(err, 1e-16)
        push!(errors, err)
        push!(mode_numbers, mode)
    end

    scatterlines!(ax2, mode_numbers, errors, color = NAVY, linewidth = 1.5,
                  markersize = 6, markercolor = NAVY,
                  strokecolor = :white, strokewidth = 0.5)

    hidespines!(ax2, :t, :r)
    hlines!(ax2, [2.2e-16], color = :gray, linestyle = :dot, linewidth = 1)
    text!(ax2, n_compare, 4e-16, text = "Machine precision", fontsize = 9,
          color = :gray, align = (:right, :bottom))

    # Main title
    Label(fig[0, :],
          L"Vibrating String: $-u_{xx} = \lambda u$, $u(0)=0$, $u'(1)=0$ ($N = %$(N)$)",
          fontsize = 13)

    # Save figure
    save(joinpath(OUTPUT_DIR, "vibrating_string.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "vibrating_string.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "vibrating_string.pdf"))")

    # =========================================================================
    # Print summary
    # =========================================================================
    println("\nVibrating String Eigenvalue Problem")
    println("="^60)
    @printf("  N = %d, domain [0, 1]\n", N)
    println("  BCs: u(0) = 0 (Dirichlet), u'(1) = 0 (Neumann)")
    println("  Exact eigenvalues: lambda_n = ((2n-1)*pi/2)^2")
    println()
    println("-"^60)
    @printf("%6s %16s %16s %12s\n", "Mode", "Computed", "Exact", "Error")
    println("-"^60)
    for mode in 1:min(10, n_compare)
        idx = valid_idx[mode]
        @printf("%6d %16.10f %16.10f %12.2e\n",
                mode, lambdas[idx], exact_lambdas[mode],
                abs(lambdas[idx] - exact_lambdas[mode]))
    end
    println("-"^60)
end

main()
