# ho_coupled_comparison.jl
#
# Chapter 14: Higher-Order Boundary Value Problems
#
# Direct vs. coupled system comparison (Computational Etude 14.3).
#
# Solves the clamped beam problem:
#
#     u^(4)(x) = e^x,   -1 < x < 1,   u(+-1) = u'(+-1) = 0
#
# using TWO approaches and compares their accuracy and conditioning.
#
# Approach 1 (Direct D^4):
#     The polynomial trick from Trefethen's "Spectral Methods in MATLAB"
#     (Program 38). Write u(x) = (1 - x^2) * q(x), which automatically
#     satisfies u(+-1) = 0 and (via the product rule) u'(+-1) = 0.
#     The fourth-order operator is reformulated as:
#         D4full = [diag(1 - x^2) * D^4 - 8*diag(x) * D^3 - 12*D^2] * S
#     and the interior system D4int * v = f is solved directly.
#
# Approach 2 (Coupled second-order system):
#     Introduce w = u'' to split the fourth-order ODE into two coupled
#     second-order equations:
#         u'' = w,        with u(+-1) = 0
#         w'' = e^x,      with u'(+-1) = 0
#     Assembled as a 2(N+1) x 2(N+1) block system.
#
# The direct approach yields a condition number scaling as O(N^8),
# while the coupled system, involving only D^2, scales as O(N^4).
#
# Generates Figure 14.3: convergence and conditioning comparison.
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
# Exact solution (same as Etude 14.1)
# -----------------------------------------------------------------------------

"""
    exact_solution(x)

Exact solution of u^(4)(x) = e^x with clamped BCs u(+-1) = u'(+-1) = 0.

The general solution is:
    u(x) = e^x + c_3*x^3 + c_2*x^2 + c_1*x + c_0

since e^x is a particular solution of u^(4) = e^x and x^3, x^2, x, 1
span the nullspace of d^4/dx^4.
"""
function exact_solution(x)
    e_plus = exp(1.0)
    e_minus = exp(-1.0)

    # 4x4 linear system: A * [c3, c2, c1, c0]^T = b
    A = [ 1.0  1.0  1.0  1.0;    # u(1) = 0
         -1.0  1.0 -1.0  1.0;    # u(-1) = 0
          3.0  2.0  1.0  0.0;    # u'(1) = 0
          3.0 -2.0  1.0  0.0]    # u'(-1) = 0
    b = [-e_plus, -e_minus, -e_plus, -e_minus]

    c = A \ b
    c3, c2, c1, c0 = c

    return exp.(x) .+ c3 .* x.^3 .+ c2 .* x.^2 .+ c1 .* x .+ c0
end

# -----------------------------------------------------------------------------
# Approach 1: Direct D^4 via polynomial trick
# -----------------------------------------------------------------------------

"""
    solve_direct_d4(N)

Solve u^(4) = e^x with clamped BCs using the polynomial trick.

Returns (u, x, D4int) where u is the full solution, x is the grid,
and D4int is the interior operator (for condition number).
"""
function solve_direct_d4(N)
    D, x = cheb_matrix(N)

    # Powers of the differentiation matrix
    D2 = D * D
    D3 = D2 * D
    D4 = D2 * D2

    # S maps v -> q on interior points, zeros at boundaries
    s = zeros(N + 1)
    s[2:N] = 1.0 ./ (1.0 .- x[2:N].^2)
    S = Diagonal(s)

    # Full fourth-order operator with the polynomial trick
    D4full = (Diagonal(1.0 .- x.^2) * D4 - 8.0 * Diagonal(x) * D3 - 12.0 * D2) * S

    # Extract interior block
    D4int = D4full[2:N, 2:N]

    # Right-hand side at interior points
    f = exp.(x[2:N])

    # Solve
    v = D4int \ f

    # Reconstruct full solution (boundary values are zero)
    u = zeros(N + 1)
    u[2:N] = v

    return u, x, D4int
end

# -----------------------------------------------------------------------------
# Approach 2: Coupled second-order system
# -----------------------------------------------------------------------------

"""
    solve_coupled_d2(N)

Solve u^(4) = e^x with clamped BCs via a coupled system of two D^2 equations.

Introduce w = u'' to obtain:
    u'' = w        (with u(+-1) = 0)
    w'' = e^x      (with u'(+-1) = 0 encoded via the D matrix)

Returns (u, x, A_block) where u is the full solution, x is the grid,
and A_block is the full block system (for condition number).
"""
function solve_coupled_d2(N)
    D, x = cheb_matrix(N)
    D2 = D * D
    n = N + 1   # Number of grid points
    I_mat = Matrix(1.0I, n, n)
    Z = zeros(n, n)

    # Block matrices
    A11 = copy(D2)        # u'' part
    A12 = -copy(I_mat)    # -w part  (u'' - w = 0)
    A21 = copy(Z)         # coupling for w equation
    A22 = copy(D2)        # w'' part

    # Right-hand sides
    b1 = zeros(n)
    b2 = exp.(x)

    # ----- Boundary conditions for the u-block (Dirichlet: u(+-1) = 0) -----
    # Row 1: u(x_1) = u(1) = 0  (Julia 1-based: index 1 corresponds to x=1)
    A11[1, :] .= 0.0
    A11[1, 1] = 1.0
    A12[1, :] .= 0.0
    b1[1] = 0.0

    # Row N+1: u(x_{N+1}) = u(-1) = 0
    A11[n, :] .= 0.0
    A11[n, n] = 1.0
    A12[n, :] .= 0.0
    b1[n] = 0.0

    # ----- Boundary conditions for the w-block (Neumann: u'(+-1) = 0) -----
    # Row 1 of w-block: u'(x_1) = u'(1) = D[1,:] * u = 0
    A21[1, :] = D[1, :]
    A22[1, :] .= 0.0
    b2[1] = 0.0

    # Row N+1 of w-block: u'(x_{N+1}) = u'(-1) = D[N+1,:] * u = 0
    A21[n, :] = D[n, :]
    A22[n, :] .= 0.0
    b2[n] = 0.0

    # Assemble the full block system
    A_block = [A11 A12; A21 A22]
    rhs = [b1; b2]

    # Solve
    sol = A_block \ rhs
    u = sol[1:n]

    return u, x, A_block
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    mkpath(OUTPUT_DIR)

    # =========================================================================
    # Convergence and conditioning study
    # =========================================================================
    N_values = collect(6:2:40)

    errors_direct = Float64[]
    errors_coupled = Float64[]
    cond_direct = Float64[]
    cond_coupled = Float64[]

    for N_val in N_values
        # --- Direct D^4 approach ---
        u_d, x_d, D4int = solve_direct_d4(N_val)
        u_exact_d = exact_solution(x_d)
        err_d = maximum(abs.(u_d .- u_exact_d))
        push!(errors_direct, max(err_d, 1e-16))
        push!(cond_direct, cond(D4int))

        # --- Coupled D^2 approach ---
        u_c, x_c, A_block = solve_coupled_d2(N_val)
        u_exact_c = exact_solution(x_c)
        err_c = maximum(abs.(u_c .- u_exact_c))
        push!(errors_coupled, max(err_c, 1e-16))
        push!(cond_coupled, cond(A_block))
    end

    # Print convergence table
    println("Convergence Table: Direct D^4 vs. Coupled D^2")
    println("-"^66)
    @printf("%4s  %13s  %13s  %14s  %14s\n",
            "N", "Err (Direct)", "Err (Coupled)", "cond(D4int)", "cond(Block)")
    println("-"^66)
    for (i, N_val) in enumerate(N_values)
        @printf("%4d  %13.4e  %13.4e  %14.4e  %14.4e\n",
                N_val, errors_direct[i], errors_coupled[i],
                cond_direct[i], cond_coupled[i])
    end
    println("-"^66)

    # Fit condition number slopes in log-log
    n_fit = length(N_values) ÷ 2
    log_N = log.(Float64.(N_values[n_fit:end]))

    # Direct D^4 slope
    log_cond_d = log.(cond_direct[n_fit:end])
    # Simple linear regression
    n_pts = length(log_N)
    slope_d = (n_pts * sum(log_N .* log_cond_d) - sum(log_N) * sum(log_cond_d)) /
              (n_pts * sum(log_N.^2) - sum(log_N)^2)

    # Coupled D^2 slope
    log_cond_c = log.(cond_coupled[n_fit:end])
    slope_c = (n_pts * sum(log_N .* log_cond_c) - sum(log_N) * sum(log_cond_c)) /
              (n_pts * sum(log_N.^2) - sum(log_N)^2)

    @printf("\nFitted condition number slopes (log-log):\n")
    @printf("  Direct D^4:  %.2f  (expected ~8)\n", slope_d)
    @printf("  Coupled D^2: %.2f  (expected ~4)\n", slope_c)

    # =========================================================================
    # Figure: two-panel comparison
    # =========================================================================
    fig = Figure(size = (900, 380))

    # ---- Left panel: max error vs N (semilog) ----
    ax1 = Axis(fig[1, 1],
               xlabel = L"N",
               ylabel = "Max error",
               title = "Convergence comparison",
               yscale = log10)

    scatterlines!(ax1, N_values, errors_direct, color = CORAL, linewidth = 1.5,
                  markersize = 5, markercolor = CORAL,
                  strokecolor = :white, strokewidth = 0.5,
                  label = L"Direct $D^4$ (polynomial trick)")
    scatterlines!(ax1, N_values, errors_coupled, color = TEAL, linewidth = 1.5,
                  markersize = 5, markercolor = TEAL, marker = :rect,
                  strokecolor = :white, strokewidth = 0.5,
                  label = L"Coupled $D^2$ system")

    xlims!(ax1, 4, 44)
    ylims!(ax1, 1e-16, 1e2)
    hidespines!(ax1, :t, :r)

    # Machine precision line
    hlines!(ax1, [2.2e-16], color = :gray, linestyle = :dot, linewidth = 0.8)
    text!(ax1, 43, 4e-16, text = L"\epsilon_{\mathrm{mach}}", fontsize = 9,
          color = :gray, align = (:right, :bottom))

    axislegend(ax1, position = :rt, framevisible = true, labelsize = 9)

    # ---- Right panel: condition number vs N (log-log) ----
    ax2 = Axis(fig[1, 2],
               xlabel = L"N",
               ylabel = "Condition number",
               title = "Conditioning comparison",
               xscale = log10,
               yscale = log10)

    scatterlines!(ax2, N_values, cond_direct, color = CORAL, linewidth = 1.5,
                  markersize = 5, markercolor = CORAL,
                  strokecolor = :white, strokewidth = 0.5,
                  label = L"Direct $D^4$: $\mathcal{O}(N^8)$")
    scatterlines!(ax2, N_values, cond_coupled, color = TEAL, linewidth = 1.5,
                  markersize = 5, markercolor = TEAL, marker = :rect,
                  strokecolor = :white, strokewidth = 0.5,
                  label = L"Coupled $D^2$: $\mathcal{O}(N^4)$")

    # Reference slopes
    N_ref = [10.0, 40.0]
    ref_N_mid = Float64(N_values[length(N_values) ÷ 2])

    # Direct reference: N^8
    cond_d_mid = cond_direct[length(N_values) ÷ 2]
    lines!(ax2, N_ref, cond_d_mid .* (N_ref ./ ref_N_mid).^8,
           color = (CORAL, 0.5), linewidth = 0.8, linestyle = :dash)

    # Coupled reference: N^4
    cond_c_mid = cond_coupled[length(N_values) ÷ 2]
    lines!(ax2, N_ref, cond_c_mid .* (N_ref ./ ref_N_mid).^4,
           color = (TEAL, 0.5), linewidth = 0.8, linestyle = :dash)

    hidespines!(ax2, :t, :r)
    axislegend(ax2, position = :lt, framevisible = true, labelsize = 9)

    # Main title
    Label(fig[0, :],
          L"Computational Etude 14.3: Direct $D^4$ vs. Coupled $D^2$ System",
          fontsize = 13)

    # Save figure
    save(joinpath(OUTPUT_DIR, "coupled_comparison.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "coupled_comparison.png"), fig, px_per_unit = 4)
    println("\nFigure saved to: $(joinpath(OUTPUT_DIR, "coupled_comparison.pdf"))")
end

main()
