# collocation_vs_galerkin.jl
#
# Compares the Collocation (pseudospectral) and Galerkin methods for solving
# a reaction-diffusion boundary value problem. This example illustrates two
# fundamental approaches to the Method of Weighted Residuals.
#
# Problem:
#     u''(x) - 4u(x) = -1,   -1 <= x <= 1
#     u(-1) = 0,  u(1) = 0
#
# Exact solution:
#     u_exact(x) = (1/4)(1 - cosh(2x)/cosh(2))
#
# Basis functions (satisfy homogeneous BCs):
#     phi_0(x) = (1 - x^2)
#     phi_1(x) = (1 - x^2)*x^2
#
# Trial function:
#     u_1(x) = a_0*phi_0(x) + a_1*phi_1(x)
#
# Methods compared:
#     1. Collocation: Force residual to zero at x = 0 and x = 0.5
#     2. Galerkin: Force residual orthogonal to basis functions
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
using QuadGK
using Printf

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
const N_POINTS = 500

# Book color scheme
const NAVY  = colorant"#142D6E"
const SKY   = colorant"#7896D2"
const CORAL = colorant"#E74C3C"
const GREEN = colorant"#27AE60"

# Output paths
const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch03", "julia")

# Publication-quality theme
set_theme!(Theme(
    fontsize = 11,
    fonts = (regular = "CMU Serif", bold = "CMU Serif Bold", italic = "CMU Serif Italic"),
    Axis = (xlabelsize = 11, ylabelsize = 11, titlesize = 12,
            xticklabelsize = 10, yticklabelsize = 10,
            spinewidth = 0.8, xtickwidth = 0.8, ytickwidth = 0.8),
    Legend = (labelsize = 10,),
))

# -----------------------------------------------------------------------------
# Exact solution
# -----------------------------------------------------------------------------
"""
    u_exact(x)

Exact solution: u(x) = (1/4)(1 - cosh(2x)/cosh(2)).

Verification:
    u'(x)  = -(1/4)(2sinh(2x)/cosh(2)) = -(1/2)sinh(2x)/cosh(2)
    u''(x) = -cosh(2x)/cosh(2)

Substituting into ODE:
    u'' - 4u = -cosh(2x)/cosh(2) - 4*(1/4)*(1 - cosh(2x)/cosh(2))
             = -cosh(2x)/cosh(2) - 1 + cosh(2x)/cosh(2)
             = -1

Boundary conditions:
    u(+/-1) = (1/4)(1 - cosh(2)/cosh(2)) = (1/4)(1 - 1) = 0
"""
u_exact(x) = 0.25 * (1 - cosh(2x) / cosh(2))

# -----------------------------------------------------------------------------
# Basis functions
# -----------------------------------------------------------------------------
"""phi_0(x) = 1 - x^2 (vanishes at x = +/-1)"""
phi0(x) = 1 - x^2

"""phi_1(x) = x^2 - x^4 = (1 - x^2)*x^2 (vanishes at x = +/-1)"""
phi1(x) = x^2 - x^4

"""Second derivative: phi_0''(x) = -2"""
phi0_xx(x) = -2.0

"""Second derivative: phi_1''(x) = 2 - 12x^2"""
phi1_xx(x) = 2 - 12x^2

# -----------------------------------------------------------------------------
# Collocation method
# -----------------------------------------------------------------------------
"""
    solve_collocation()

Solve using collocation at x = 0 and x = 0.5.

The residual is:
    R(x) = u_1''(x) - 4*u_1(x) + 1
         = a_0*(phi_0'' - 4*phi_0) + a_1*(phi_1'' - 4*phi_1) + 1

Setting R(x) = 0 at the collocation points yields a 2x2 linear system.
"""
function solve_collocation()
    # Operator L = d^2/dx^2 - 4
    L_phi0(x) = phi0_xx(x) - 4 * phi0(x)
    L_phi1(x) = phi1_xx(x) - 4 * phi1(x)

    # Collocation points
    x1, x2 = 0.0, 0.5

    # Build system matrix
    A = [L_phi0(x1) L_phi1(x1);
         L_phi0(x2) L_phi1(x2)]

    # Right-hand side: f = -1, so R = Lu - f = Lu + 1 = 0 means Lu = -1
    b = [-1.0, -1.0]

    # Solve
    coeffs = A \ b

    return coeffs[1], coeffs[2]
end

# -----------------------------------------------------------------------------
# Galerkin method
# -----------------------------------------------------------------------------
"""
    solve_galerkin()

Solve using Galerkin method: <R, phi_k> = 0 for k = 0, 1.

This requires computing integrals:
    A_{ij} = integral_{-1}^{1} (phi_i'' - 4*phi_i) * phi_j dx
    b_i    = integral_{-1}^{1} (-1) * phi_i dx

The matrix is symmetric for self-adjoint operators.
"""
function solve_galerkin()
    # Operator L = d^2/dx^2 - 4
    L_phi0(x) = phi0_xx(x) - 4 * phi0(x)
    L_phi1(x) = phi1_xx(x) - 4 * phi1(x)

    # Compute matrix entries A_{ij} = <L[phi_j], phi_i>
    A00, _ = quadgk(x -> L_phi0(x) * phi0(x), -1, 1)
    A01, _ = quadgk(x -> L_phi1(x) * phi0(x), -1, 1)
    A10, _ = quadgk(x -> L_phi0(x) * phi1(x), -1, 1)
    A11, _ = quadgk(x -> L_phi1(x) * phi1(x), -1, 1)

    A = [A00 A01; A10 A11]

    # Compute RHS: b_i = <-1, phi_i> = -integral(phi_i dx)
    b0, _ = quadgk(x -> -1 * phi0(x), -1, 1)
    b1, _ = quadgk(x -> -1 * phi1(x), -1, 1)

    b = [b0, b1]

    # Solve
    coeffs = A \ b

    return coeffs[1], coeffs[2]
end

# -----------------------------------------------------------------------------
# Approximate solutions
# -----------------------------------------------------------------------------
"""Collocation/Galerkin approximation with given coefficients."""
u_approx(x, a0, a1) = a0 * phi0(x) + a1 * phi1(x)

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Solve both systems
    a0_coll, a1_coll = solve_collocation()
    a0_gal, a1_gal = solve_galerkin()

    println("Collocation coefficients:")
    @printf("  a_0 = %.6f\n", a0_coll)
    @printf("  a_1 = %.6f\n", a1_coll)
    println("\nGalerkin coefficients:")
    @printf("  a_0 = %.6f\n", a0_gal)
    @printf("  a_1 = %.6f\n", a1_gal)

    # Spatial grid
    x = range(-1, 1, length=N_POINTS)

    # Compute solutions
    u_ex = u_exact.(x)
    u_coll = u_approx.(x, a0_coll, a1_coll)
    u_gal  = u_approx.(x, a0_gal, a1_gal)

    error_coll = u_ex .- u_coll
    error_gal  = u_ex .- u_gal

    # Collocation points for visualization
    x_coll_pts = [0.0, 0.5]
    u_coll_pts = u_approx.(x_coll_pts, a0_coll, a1_coll)

    # Create two-panel figure
    fig = Figure(size=(720, 280))

    # --- Left panel: All solutions ---
    ax1 = Axis(fig[1, 1],
               xlabel = L"x",
               ylabel = L"u(x)",
               title = "Solutions Comparison",
               titlesize = 11)
    xlims!(ax1, -1, 1)

    # Exact solution (solid line)
    lines!(ax1, collect(x), u_ex, color=NAVY, linewidth=1.8,
           label="Exact")

    # Collocation approximation (dashed line)
    lines!(ax1, collect(x), u_coll, color=CORAL, linewidth=1.5,
           linestyle=:dash, label="Collocation")

    # Galerkin approximation (dash-dot line)
    lines!(ax1, collect(x), u_gal, color=GREEN, linewidth=1.5,
           linestyle=:dashdot, label="Galerkin")

    # Collocation points (open squares)
    scatter!(ax1, x_coll_pts, u_coll_pts, color=:transparent, markersize=10,
             marker=:rect, strokewidth=1.5, strokecolor=CORAL,
             label="Collocation pts")

    axislegend(ax1, position=:cb, framevisible=false, bgcolor=(:white, 0.9),
               nbanks=2)

    # Style: hide top/right spines, add subtle grid
    hidespines!(ax1, :t, :r)
    ax1.ygridvisible = true
    ax1.ygridstyle = :solid
    ax1.ygridcolor = (:gray, 0.2)

    # --- Right panel: Errors ---
    ax2 = Axis(fig[1, 2],
               xlabel = L"x",
               ylabel = L"u_{\mathrm{exact}} - u_{\mathrm{approx}}",
               title = "Error Comparison",
               titlesize = 11)
    xlims!(ax2, -1, 1)

    lines!(ax2, collect(x), error_coll, color=CORAL, linewidth=1.5,
           label="Collocation")
    lines!(ax2, collect(x), error_gal, color=GREEN, linewidth=1.5,
           label="Galerkin")
    hlines!(ax2, [0.0], color=:gray, linewidth=0.5, linestyle=:dash)

    axislegend(ax2, position=:rt, framevisible=false, bgcolor=(:white, 0.9))

    # Style: hide top/right spines, add subtle grid
    hidespines!(ax2, :t, :r)
    ax2.ygridvisible = true
    ax2.ygridstyle = :solid
    ax2.ygridcolor = (:gray, 0.2)

    # Save figure
    mkpath(OUTPUT_DIR)
    save(joinpath(OUTPUT_DIR, "collocation_vs_galerkin.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "collocation_vs_galerkin.png"), fig, px_per_unit=4)

    println("\nFigure saved to: $(joinpath(OUTPUT_DIR, "collocation_vs_galerkin.pdf"))")

    # Print comparison table
    println("\n" * "=" ^ 65)
    println("Comparison at u(0) - the central maximum")
    println("=" ^ 65)
    @printf("%-20s %-15s %-15s\n", "Method", "u(0)", "Abs. Error")
    println("-" ^ 65)
    @printf("%-20s %-15.6f %-15.6f\n", "Exact", u_exact(0), 0.0)
    @printf("%-20s %-15.6f %-15.6f\n", "Collocation (N=2)",
            u_approx(0, a0_coll, a1_coll),
            abs(u_exact(0) - u_approx(0, a0_coll, a1_coll)))
    @printf("%-20s %-15.6f %-15.6f\n", "Galerkin (N=2)",
            u_approx(0, a0_gal, a1_gal),
            abs(u_exact(0) - u_approx(0, a0_gal, a1_gal)))
    println("=" ^ 65)

    # RMS errors
    rms_coll = sqrt(sum(error_coll .^ 2) / N_POINTS)
    rms_gal  = sqrt(sum(error_gal .^ 2) / N_POINTS)
    @printf("\nRMS Error (Collocation): %.6f\n", rms_coll)
    @printf("RMS Error (Galerkin):    %.6f\n", rms_gal)

    # Max errors
    max_coll = maximum(abs.(error_coll))
    max_gal  = maximum(abs.(error_gal))
    @printf("\nMax Error (Collocation): %.6f\n", max_coll)
    @printf("Max Error (Galerkin):    %.6f\n", max_gal)

    return fig
end

# Run if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
