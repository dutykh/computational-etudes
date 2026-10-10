# bc_radiative.jl
#
# Chapter 13: Boundary Conditions in Spectral Methods
#
# Radiative cooling problem: heat equation with nonlinear Robin BC.
#
#     u_t = u_xx   on [0, 1],  t > 0
#     u(0, t) = 1                    (hot wall, Dirichlet)
#     u_x(1, t) + sigma * u(1, t)^4 = 0   (radiative cooling, nonlinear Robin)
#
# Initial condition: u(x, 0) = 1 (uniform temperature)
#
# The domain [0, 1] is mapped to the Chebyshev domain [-1, 1] via x = (xi + 1)/2,
# so D_x = 2 * D_xi.
#
# Chebyshev ordering: xi[1] = 1 -> physical x = 1 (radiating end)
#                     xi[N+1] = -1 -> physical x = 0 (hot wall)
#
# Backward Euler with Newton iteration for the nonlinear boundary condition.
#
# Generates Figure 13.6: Radiative cooling profiles and boundary temperature.
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
# Solver
# -----------------------------------------------------------------------------
"""
    solve_radiative_cooling(; N=40, sigma=1.0, tmax=5.0, dt=0.005,
                              snapshot_times=[0.0, 0.1, 0.5, 1.0, 5.0])

Solve the radiative cooling problem using Backward Euler + Newton.

Returns `(x_phys, snapshots, t_history, u_boundary)` where
  - `x_phys`: physical grid points on [0, 1]
  - `snapshots`: dict of time => solution
  - `t_history`: vector of all time values
  - `u_boundary`: u(1, t) for each time step
"""
function solve_radiative_cooling(; N=40, sigma=1.0, tmax=5.0, dt=0.005,
                                   snapshot_times=[0.0, 0.1, 0.5, 1.0, 5.0])
    D, xi = cheb_matrix(N)
    D2 = D * D
    n = N + 1

    # Map to physical domain [0, 1]: x = (xi + 1)/2
    # D_x = 2 * D_xi, D2_x = 4 * D2_xi
    D_phys = 2.0 * D
    D2_phys = 4.0 * D2
    x_phys = (xi .+ 1.0) ./ 2.0  # x_phys[1]=1, x_phys[N+1]=0

    # Initial condition: u(x, 0) = 1
    u = ones(n)

    # Storage for boundary temperature history
    t_history = [0.0]
    u_boundary = [u[1]]  # u at x=1 (radiating end) = xi[1]=1

    # Storage for snapshots
    snapshots = Dict{Float64, Vector{Float64}}()
    # Store IC if 0.0 is in snapshot_times
    if 0.0 in snapshot_times
        snapshots[0.0] = copy(u)
    end

    # Identity matrix
    I_mat = Matrix(1.0I, n, n)

    t = 0.0
    nsteps = Int(ceil(tmax / dt))

    for step in 1:nsteps
        t += dt

        # Backward Euler: (I - dt * D2_phys) u^{n+1} = u^n
        # with BCs applied via Newton for the nonlinear Robin condition

        # Build the system matrix
        A = I_mat - dt * D2_phys

        # Newton iteration
        u_new = copy(u)
        for newton_iter in 1:20
            # Residual
            F = A * u_new - u

            # BC row N+1 (xi=-1 -> x=0): Dirichlet u(0) = 1
            F[n] = u_new[n] - 1.0

            # BC row 1 (xi=1 -> x=1): u_x(1) + sigma * u(1)^4 = 0
            # Using physical derivative: D_phys[1,:] * u + sigma * u[1]^4 = 0
            F[1] = dot(D_phys[1, :], u_new) + sigma * u_new[1]^4

            # Jacobian
            J = copy(A)

            # BC row N+1: Dirichlet
            J[n, :] .= 0.0
            J[n, n] = 1.0

            # BC row 1: Robin (nonlinear)
            J[1, :] = D_phys[1, :]
            J[1, 1] += 4.0 * sigma * u_new[1]^3

            # Newton update
            delta = J \ F
            u_new .-= delta

            if norm(delta, Inf) < 1e-12
                break
            end
        end

        u = u_new

        # Record boundary temperature
        push!(t_history, t)
        push!(u_boundary, u[1])

        # Save snapshot if time matches
        for t_snap in snapshot_times
            if abs(t - t_snap) < 0.5 * dt && !haskey(snapshots, t_snap)
                snapshots[t_snap] = copy(u)
            end
        end
    end

    return x_phys, snapshots, t_history, u_boundary
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    mkpath(OUTPUT_DIR)

    N = 40

    # =========================================================================
    # Figure: Two panels
    # =========================================================================
    fig = Figure(size = (900, 380))

    # ---- Panel 1: Temperature profiles at selected times (sigma = 1) ----
    ax1 = Axis(fig[1, 1],
               xlabel = L"x",
               ylabel = L"u(x, t)",
               title = L"Temperature profiles ($\sigma = 1$)")

    snapshot_times = [0.0, 0.1, 0.5, 1.0, 5.0]
    colors_snap = [NAVY, SKY, TEAL, CORAL, PURPLE]

    x_phys, snapshots, _, _ = solve_radiative_cooling(
        N=N, sigma=1.0, tmax=5.0, dt=0.005, snapshot_times=snapshot_times)

    # Physical x is ordered from x_phys[1]=1 to x_phys[N+1]=0
    # Sort for plotting
    idx_sort = sortperm(x_phys)
    x_plot = x_phys[idx_sort]

    for (i, t_snap) in enumerate(snapshot_times)
        if haskey(snapshots, t_snap)
            u_snap = snapshots[t_snap][idx_sort]
            lines!(ax1, x_plot, u_snap, color = colors_snap[i], linewidth = 1.5,
                   label = L"t = %$(t_snap)")
        end
    end

    xlims!(ax1, -0.02, 1.02)
    ylims!(ax1, 0.0, 1.1)
    axislegend(ax1, position = :lb, framevisible = true, labelsize = 9)
    hidespines!(ax1, :t, :r)

    # ---- Panel 2: Boundary temperature u(1,t) for various sigma ----
    ax2 = Axis(fig[1, 2],
               xlabel = L"t",
               ylabel = L"u(1, t)",
               title = "Boundary temperature at radiating end")

    sigma_values = [0.1, 1.0, 10.0]
    colors_sigma = [TEAL, CORAL, PURPLE]

    for (i, sigma) in enumerate(sigma_values)
        _, _, t_hist, u_bnd = solve_radiative_cooling(
            N=N, sigma=sigma, tmax=5.0, dt=0.005, snapshot_times=Float64[])

        lines!(ax2, t_hist, u_bnd, color = colors_sigma[i], linewidth = 1.5,
               label = L"\sigma = %$(sigma)")
    end

    xlims!(ax2, -0.1, 5.1)
    ylims!(ax2, 0.0, 1.05)
    axislegend(ax2, position = :rt, framevisible = true, labelsize = 9)
    hidespines!(ax2, :t, :r)

    # Main title
    Label(fig[0, :],
          L"Radiative Cooling: $u_t = u_{xx}$, $u(0)=1$, $u_x(1) + \sigma u(1)^4 = 0$",
          fontsize = 13)

    # Save figure
    save(joinpath(OUTPUT_DIR, "radiative_cooling.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "radiative_cooling.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "radiative_cooling.pdf"))")

    # Print summary
    println("\nRadiative cooling simulation completed.")
    println("  N = $N, dt = 0.005, tmax = 5.0")
    println("  Snapshot times: $snapshot_times")
    println("  Sigma values tested: $sigma_values")
end

main()
