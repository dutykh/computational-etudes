# ho_clamped_beam.jl
#
# Chapter 14: Higher-Order Boundary Value Problems
#
# Clamped beam under exponential load (Computational Etude 14.1).
#
# Solves the fourth-order BVP (biharmonic-type equation):
#
#     u^(4)(x) = e^x,   -1 < x < 1,   u(+-1) = u'(+-1) = 0
#
# using a Chebyshev spectral method with the polynomial trick from
# Trefethen's "Spectral Methods in MATLAB" (Program 38).
#
# The key idea is to write u(x) = (1 - x^2) * q(x), which automatically
# satisfies both Dirichlet conditions u(+-1) = 0.  The factor (1 - x^2)
# also enforces u'(+-1) = 0 because u'(x) = -2x*q(x) + (1 - x^2)*q'(x),
# so u'(+-1) = -2*(+-1)*q(+-1) + 0 = 0 only if q(+-1) = 0.  Instead of
# enforcing this directly, the fourth derivative operator is reformulated
# via the product rule as:
#
#     D4full = [diag(1 - x^2) * D^4 - 8 * diag(x) * D^3 - 12 * D^2] * S
#
# where S = diag(0, 1/(1 - x_1^2), ..., 1/(1 - x_{N-1}^2), 0) maps v
# back to q on interior points and annihilates boundary values.
#
# The exact solution is u(x) = e^x + c_3*x^3 + c_2*x^2 + c_1*x + c_0,
# where the constants are determined by the four clamped boundary conditions.
#
# Generates Figure 14.1: Clamped beam solution and convergence/conditioning.
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
# Biharmonic solver using Trefethen's polynomial trick (Program 38)
# -----------------------------------------------------------------------------

"""
    biharmonic_clamped_matrix(N)

Build the fourth-order operator for clamped BCs via the polynomial trick.

The approach writes u(x) = (1 - x^2) * q(x) and reformulates the
fourth derivative operator using the product rule.  The resulting
system is restricted to interior Chebyshev points.

Returns `(D4int, x)` where D4int is the (N-1)x(N-1) interior operator
and x is the full Chebyshev grid.
"""
function biharmonic_clamped_matrix(N)
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
    # Derived from the Leibniz rule for d^4/dx^4 [(1 - x^2) * q(x)]
    D4full = (Diagonal(1.0 .- x.^2) * D4 - 8.0 * Diagonal(x) * D3 - 12.0 * D2) * S

    # Extract interior block (strip first and last rows/columns)
    # Julia 1-based: interior indices are 2:N
    D4int = D4full[2:N, 2:N]

    return D4int, x
end

# -----------------------------------------------------------------------------
# Exact solution
# -----------------------------------------------------------------------------

"""
    exact_solution(x)

Exact solution of u^(4)(x) = e^x with clamped BCs u(+-1) = u'(+-1) = 0.

The general solution is:
    u(x) = e^x + c_3*x^3 + c_2*x^2 + c_1*x + c_0

since e^x is a particular solution of u^(4) = e^x and x^3, x^2, x, 1
span the nullspace of d^4/dx^4.

The four constants are determined by the boundary conditions:
    u(1) = 0,  u(-1) = 0,  u'(1) = 0,  u'(-1) = 0
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
# Main visualization
# -----------------------------------------------------------------------------
function main()
    mkpath(OUTPUT_DIR)

    # =========================================================================
    # Solve for N = 15 (demonstration)
    # =========================================================================
    N_demo = 15
    D4int, x = biharmonic_clamped_matrix(N_demo)

    # Right-hand side at interior points
    f = exp.(x[2:N_demo])

    # Solve the linear system
    v = D4int \ f

    # Reconstruct full solution (boundary values are zero)
    u_num = zeros(N_demo + 1)
    u_num[2:N_demo] = v

    # Exact solution on the full grid
    u_ex = exact_solution(x)

    demo_error = maximum(abs.(u_num .- u_ex))
    @printf("N = %d: max |u_num - u_exact| = %.4e\n", N_demo, demo_error)

    # Fine grid for smooth exact solution curve
    x_fine = collect(range(-1, 1, length=500))
    u_fine = exact_solution(x_fine)

    # =========================================================================
    # Convergence study
    # =========================================================================
    N_values = collect(4:2:32)
    errors = Float64[]
    cond_numbers = Float64[]

    for N_val in N_values
        D4int_n, x_n = biharmonic_clamped_matrix(N_val)
        f_n = exp.(x_n[2:N_val])

        v_n = D4int_n \ f_n

        u_n = zeros(N_val + 1)
        u_n[2:N_val] = v_n

        u_ex_n = exact_solution(x_n)
        err = maximum(abs.(u_n .- u_ex_n))
        push!(errors, max(err, 1e-16))
        push!(cond_numbers, cond(D4int_n))
    end

    # Print convergence table
    println("\nConvergence Table: Clamped Beam")
    println("-"^50)
    @printf("%6s %14s %14s\n", "N", "Max Error", "cond(D4int)")
    println("-"^50)
    for (N_val, err, kappa) in zip(N_values, errors, cond_numbers)
        @printf("%6d %14.4e %14.4e\n", N_val, err, kappa)
    end
    println("-"^50)

    # =========================================================================
    # Figure: two-panel plot
    # =========================================================================
    fig = Figure(size = (900, 380))

    # ---- Left panel: solution for N = 15 ----
    ax1 = Axis(fig[1, 1],
               xlabel = L"x",
               ylabel = L"u(x)",
               title = @sprintf("Clamped beam: N = %d", N_demo))

    lines!(ax1, x_fine, u_fine, color = NAVY, linewidth = 1.5,
           label = "Exact solution")
    scatter!(ax1, x, u_num, color = CORAL, markersize = 8,
             strokecolor = :white, strokewidth = 0.5,
             label = @sprintf("Spectral (N = %d)", N_demo))

    xlims!(ax1, -1.05, 1.05)
    axislegend(ax1, position = :lt, framevisible = true, labelsize = 9)
    hidespines!(ax1, :t, :r)

    # ---- Right panel: convergence + condition number ----
    ax2 = Axis(fig[1, 2],
               xlabel = L"N",
               ylabel = "Max error",
               title = "Convergence and conditioning",
               yscale = log10)

    scatterlines!(ax2, N_values, errors, color = CORAL, linewidth = 1.5,
                  markersize = 5, markercolor = CORAL,
                  strokecolor = :white, strokewidth = 0.5,
                  label = "Max error")

    # Machine precision line
    hlines!(ax2, [2.2e-16], color = :gray, linestyle = :dot, linewidth = 1)
    text!(ax2, 33, 4e-16, text = L"\epsilon_{\mathrm{mach}}", fontsize = 9,
          color = :gray, align = (:right, :bottom))

    xlims!(ax2, 2, 34)
    ylims!(ax2, 1e-16, 1e0)
    hidespines!(ax2, :t, :r)

    # Secondary y-axis: condition number
    ax2b = Axis(fig[1, 2],
                ylabel = L"\mathrm{cond}(D_4)",
                yscale = log10,
                yaxisposition = :right,
                ylabelcolor = TEAL,
                yticklabelcolor = TEAL,
                ytickcolor = TEAL)
    hidespines!(ax2b)
    hidexdecorations!(ax2b)

    scatterlines!(ax2b, N_values, cond_numbers, color = TEAL, linewidth = 1.2,
                  markersize = 4, markercolor = TEAL, marker = :rect,
                  linestyle = :dash,
                  strokecolor = :white, strokewidth = 0.4,
                  label = L"\kappa(D_4)")

    xlims!(ax2b, 2, 34)

    # Combined legend
    Legend(fig[1, 2],
           [LineElement(color = CORAL, linewidth = 1.5),
            LineElement(color = TEAL, linewidth = 1.2, linestyle = :dash)],
           ["Max error", L"\kappa(D_4)"],
           tellwidth = false, tellheight = false,
           halign = :left, valign = :center,
           labelsize = 9, framevisible = true,
           margin = (10, 10, 10, 10))

    # Main title
    Label(fig[0, :],
          L"Computational Etude 14.1: Clamped Beam $u^{(4)} = e^x$, $u(\pm 1) = u'(\pm 1) = 0$",
          fontsize = 13)

    # Save figure
    save(joinpath(OUTPUT_DIR, "clamped_beam.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "clamped_beam.png"), fig, px_per_unit = 4)
    println("\nFigure saved to: $(joinpath(OUTPUT_DIR, "clamped_beam.pdf"))")
end

main()
