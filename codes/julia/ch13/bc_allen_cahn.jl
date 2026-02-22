# bc_allen_cahn.jl
#
# Chapter 13: Boundary Conditions in Spectral Methods
#
# Allen-Cahn equation with fixed and time-dependent (driven) boundary conditions.
#
#     u_t = eps * u_xx + u - u^3,   x in (-1, 1),   eps = 0.01
#
# Initial condition: u(x, 0) = 0.53 x + 0.47 sin(-1.5 pi x)
#
# Part A -- Fixed BCs: u(-1) = -1, u(1) = 1
#     Lifting: w(x) = x satisfies w(-1) = -1, w(1) = 1
#     v = u - x satisfies homogeneous Dirichlet BCs
#     v_t = eps * v_xx + (v + x) - (v + x)^3
#     IMEX: diffusion implicit, reaction explicit
#
# Part B -- Driven BCs: u(-1) = -1, u(1) = 1 + sin^2(t/5)
#     Method II (boundary bordering): overwrite boundary values each step
#
# Generates:
#     Figure 13.4: Allen-Cahn with fixed BCs (space-time plot)
#     Figure 13.5: Allen-Cahn with driven BCs (space-time plot)
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
# Problem parameters
# -----------------------------------------------------------------------------
const EPS_AC = 0.01

"""Initial condition: u(x, 0) = 0.53 x + 0.47 sin(-1.5 pi x)."""
initial_condition(x) = 0.53 .* x .+ 0.47 .* sin.(-1.5 .* pi .* x)

# -----------------------------------------------------------------------------
# Part A: Fixed BCs via lifting
# -----------------------------------------------------------------------------
"""
    solve_allen_cahn_fixed(; N=80, tmax=60.0, n_snapshots=300)

Solve Allen-Cahn with fixed BCs u(-1)=-1, u(1)=1 via lifting w(x) = x.

The variable v = u - x satisfies v(+-1) = 0 and:
    v_t = eps * v_xx + (v + x) - (v + x)^3

IMEX time stepping:
    Implicit: eps * D2 v^{n+1}
    Explicit: reaction R(v^n) = (v^n + x) - (v^n + x)^3

Semi-discrete: (I - eps*dt*D2_int) v^{n+1} = v^n + dt * R(v^n)
"""
function solve_allen_cahn_fixed(; N=80, tmax=60.0, n_snapshots=300)
    D2, D, x = cheb_second_derivative_matrix(N)

    # Time step: stability-limited
    dt = min(0.5 / N^2, 0.01) * EPS_AC * 100
    dt = min(dt, 0.01)

    # Interior indices: 2:N (Julia 1-based)
    # x[1]=1 (right), x[N+1]=-1 (left), interior: x[2:N]
    D2_int = D2[2:N, 2:N]
    I_int = Matrix(1.0I, N - 1, N - 1)

    # Implicit matrix: (I - eps * dt * D2_int)
    A = I_int - EPS_AC * dt * D2_int
    A_lu = lu(A)

    # Initial condition
    u0 = initial_condition(x)
    v = u0 .- x  # v = u - w, where w(x) = x
    v_int = copy(v[2:N])

    # Interior grid points
    x_int = x[2:N]

    # Time stepping
    nsteps = Int(ceil(tmax / dt))
    save_every = max(1, nsteps ÷ n_snapshots)

    # Store full u = v + x
    u_full = zeros(N + 1)
    u_full[1] = 1.0     # u(1) = 1
    u_full[N+1] = -1.0   # u(-1) = -1
    u_full[2:N] = v_int .+ x_int

    U_save = [copy(u_full)]
    t_save = [0.0]

    t = 0.0
    for n in 1:nsteps
        # Reaction at interior points: R(v) = (v + x) - (v + x)^3
        u_int = v_int .+ x_int
        reaction = u_int .- u_int.^3

        # IMEX step: A v^{n+1} = v^n + dt * reaction
        rhs = v_int .+ dt .* reaction
        v_int = A_lu \ rhs

        t += dt

        if n % save_every == 0
            u_full = zeros(N + 1)
            u_full[1] = 1.0
            u_full[N+1] = -1.0
            u_full[2:N] = v_int .+ x_int
            push!(U_save, copy(u_full))
            push!(t_save, t)
        end
    end

    # Convert to matrix form (snapshots as rows)
    U_mat = reduce(vcat, [u' for u in U_save])
    return x, t_save, U_mat
end

# -----------------------------------------------------------------------------
# Part B: Driven BCs via Method II
# -----------------------------------------------------------------------------
"""
    solve_allen_cahn_driven(; N=80, tmax=60.0, n_snapshots=300)

Solve Allen-Cahn with driven BC: u(-1)=-1, u(1) = 1 + sin^2(t/5).

Method II (boundary bordering): at each time step, overwrite boundary
values of u after the implicit solve.
"""
function solve_allen_cahn_driven(; N=80, tmax=60.0, n_snapshots=300)
    D2, D, x = cheb_second_derivative_matrix(N)

    # Time step
    dt = min(0.5 / N^2, 0.01) * EPS_AC * 100
    dt = min(dt, 0.01)

    # Interior system
    D2_int = D2[2:N, 2:N]
    I_int = Matrix(1.0I, N - 1, N - 1)

    # Implicit matrix
    A = I_int - EPS_AC * dt * D2_int
    A_lu = lu(A)

    # Boundary-to-interior coupling columns
    # When BCs are non-zero, we need D2[2:N, 1] * u[1] + D2[2:N, N+1] * u[N+1]
    d2_col1 = D2[2:N, 1]     # coupling to x[1] = 1 (right boundary)
    d2_colN = D2[2:N, N+1]   # coupling to x[N+1] = -1 (left boundary)

    # Boundary functions
    g_right(t) = 1.0 + sin(t / 5.0)^2    # u(1, t) = 1 + sin^2(t/5)
    g_left(t) = -1.0                       # u(-1, t) = -1

    # Initial condition
    u = initial_condition(x)
    u[1] = g_right(0.0)
    u[N+1] = g_left(0.0)

    # Time stepping
    nsteps = Int(ceil(tmax / dt))
    save_every = max(1, nsteps ÷ n_snapshots)

    U_save = [copy(u)]
    t_save = [0.0]

    t = 0.0
    for n in 1:nsteps
        # Explicit reaction at interior points
        u_int = u[2:N]
        reaction = u_int .- u_int.^3

        # Time-advanced boundary values
        t_new = t + dt
        u_R = g_right(t_new)
        u_L = g_left(t_new)

        # RHS with boundary coupling at new time
        rhs = u_int .+ dt .* reaction .+ EPS_AC .* dt .* (d2_col1 .* u_R .+ d2_colN .* u_L)

        u_int_new = A_lu \ rhs

        # Update full solution
        u[2:N] = u_int_new
        u[1] = u_R
        u[N+1] = u_L

        t = t_new

        if n % save_every == 0
            push!(U_save, copy(u))
            push!(t_save, t)
        end
    end

    # Convert to matrix form
    U_mat = reduce(vcat, [u' for u in U_save])
    return x, t_save, U_mat
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    mkpath(OUTPUT_DIR)

    N = 80

    # =========================================================================
    # Figure 1: Fixed BCs
    # =========================================================================
    println("Solving Allen-Cahn with fixed BCs...")
    x, t_save, U_save = solve_allen_cahn_fixed(N=N, tmax=60.0, n_snapshots=400)

    fig1 = Figure(size = (580, 420))
    ax1 = Axis(fig1[1, 1],
               xlabel = L"x",
               ylabel = L"t",
               title = L"Allen-Cahn: Fixed BCs $u(\pm 1) = \pm 1$, $\varepsilon = 0.01$")

    # Sort x for proper heatmap display (Chebyshev points are cos-ordered)
    idx_sort = sortperm(x)
    x_sorted = x[idx_sort]
    U_sorted = U_save[:, idx_sort]

    hm1 = heatmap!(ax1, x_sorted, t_save, U_sorted,
                    colormap = :RdBu, colorrange = (-1.0, 1.0),
                    interpolate = true)
    Colorbar(fig1[1, 2], hm1, label = L"u(x, t)")

    # Save figure 1
    save(joinpath(OUTPUT_DIR, "allen_cahn_fixed.pdf"), fig1)
    save(joinpath(OUTPUT_DIR, "allen_cahn_fixed.png"), fig1, px_per_unit = 4)
    println("Figure 1 saved to: $(joinpath(OUTPUT_DIR, "allen_cahn_fixed.pdf"))")

    # =========================================================================
    # Figure 2: Driven BCs
    # =========================================================================
    println("Solving Allen-Cahn with driven BCs...")
    x2, t_save2, U_save2 = solve_allen_cahn_driven(N=N, tmax=60.0, n_snapshots=400)

    fig2 = Figure(size = (580, 420))
    ax2 = Axis(fig2[1, 1],
               xlabel = L"x",
               ylabel = L"t",
               title = L"Allen-Cahn: Driven BC $u(1) = 1 + \sin^2(t/5)$, $\varepsilon = 0.01$")

    idx_sort2 = sortperm(x2)
    x_sorted2 = x2[idx_sort2]
    U_sorted2 = U_save2[:, idx_sort2]

    hm2 = heatmap!(ax2, x_sorted2, t_save2, U_sorted2,
                    colormap = :RdBu, colorrange = (-1.0, 2.0),
                    interpolate = true)
    Colorbar(fig2[1, 2], hm2, label = L"u(x, t)")

    # Save figure 2
    save(joinpath(OUTPUT_DIR, "allen_cahn_driven.pdf"), fig2)
    save(joinpath(OUTPUT_DIR, "allen_cahn_driven.png"), fig2, px_per_unit = 4)
    println("Figure 2 saved to: $(joinpath(OUTPUT_DIR, "allen_cahn_driven.pdf"))")

    # Print summary
    println("\nAllen-Cahn simulation parameters:")
    @printf("  N = %d, eps = %.2f\n", N, EPS_AC)
    @printf("  Fixed BCs: %d snapshots, t in [%.1f, %.1f]\n",
            length(t_save), t_save[1], t_save[end])
    @printf("  Driven BCs: %d snapshots, t in [%.1f, %.1f]\n",
            length(t_save2), t_save2[1], t_save2[end])
end

main()
