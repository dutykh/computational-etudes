# wave1d_cheb.jl
#
# Chapter 10: Spectral PDE Solvers with Chebyshev and Fourier Grids
#
# 1D wave equation solver on Chebyshev grid with leapfrog time stepping.
#
# PDE: u_tt = c^2 u_xx, -1 < x < 1, t > 0
# BCs: u(-1, t) = u(1, t) = 0 (Dirichlet)
# ICs: u(x, 0) = sin(pi/2 * (1 + x)), u_t(x, 0) = 0
#
# Generates Figure 10.3: Waterfall plot of wave evolution.
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

# Import chebfft (provides chebfft() and cheb_matrix())
include(joinpath(@__DIR__, "chebfft.jl"))

# -----------------------------------------------------------------------------
# Publication-quality CairoMakie configuration
# -----------------------------------------------------------------------------
set_theme!(Theme(
    fontsize = 11,
    fonts = (regular = "CMU Serif", bold = "CMU Serif Bold", italic = "CMU Serif Italic"),
    Axis = (xlabelsize = 11, ylabelsize = 11, titlesize = 12,
            xticklabelsize = 10, yticklabelsize = 10,
            spinewidth = 0.8, xtickwidth = 0.8, ytickwidth = 0.8),
    Legend = (labelsize = 10,),
))

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
const NAVY   = colorant"#142D6E"
const SKY    = colorant"#7896D2"
const CORAL  = colorant"#E74C3C"
const TEAL   = colorant"#16A085"

# Output path
SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch10", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Wave equation solver
# -----------------------------------------------------------------------------
"""
    wave1d_cheb(; N=80, c=1.0, tmax=6.0, alpha=4.0, n_snapshots=100)

Solve 1D wave equation on Chebyshev grid using leapfrog time stepping.

Returns `(x, t_save, U_save)` where `U_save` is a matrix with snapshots as rows.
"""
function wave1d_cheb(; N=80, c=1.0, tmax=6.0, alpha=4.0, n_snapshots=100)
    # Chebyshev grid
    x = [cos(pi * j / N) for j in 0:N]

    # Build second derivative matrix
    D, _ = cheb_matrix(N)
    D2 = D * D

    # Initial data: half-sine
    u = sin.(0.5 * pi .* (1.0 .+ x))

    # Time stepping
    dt = alpha / N^2
    nsteps = Int(ceil(tmax / dt))

    # Store snapshots
    save_every = max(1, nsteps ÷ n_snapshots)
    U_save = [copy(u)]
    t_save = [0.0]

    # Taylor start for leapfrog
    # u_prev = u - dt * u_t + 0.5 * dt^2 * u_tt
    # With u_t = 0, u_tt = c^2 * D2 * u
    u_tt0 = c^2 .* (D2 * u)
    u_tt0[1] = 0.0        # Enforce BCs
    u_tt0[N+1] = 0.0
    u_prev = u .- 0.5 * dt^2 .* u_tt0

    t = 0.0
    for n in 1:nsteps
        # Compute Laplacian
        Lap_u = c^2 .* (D2 * u)
        Lap_u[1] = 0.0      # Enforce BCs
        Lap_u[N+1] = 0.0

        # Leapfrog step
        u_new = 2.0 .* u .- u_prev .+ dt^2 .* Lap_u

        # Boundary conditions
        u_new[1] = 0.0
        u_new[N+1] = 0.0

        u_prev = u
        u = u_new
        t += dt

        # Save snapshot
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
    # Solve wave equation
    N = 80
    tmax = 6.0
    x, t_save, U_save = wave1d_cheb(N=N, tmax=tmax, n_snapshots=80)

    # Create 3D surface plot
    fig = Figure(size = (800, 560))
    ax = Axis3(fig[1, 1],
               xlabel = L"x",
               ylabel = L"t",
               zlabel = L"u(x, t)",
               title  = "1D Wave Equation: Vibrating String",
               azimuth = -60 * pi / 180,
               elevation = 25 * pi / 180)

    surface!(ax, x, t_save, U_save, colormap = :coolwarm, alpha = 0.9)
    wireframe!(ax, x, t_save, U_save, color = (:black, 0.3), linewidth = 0.2,
               overdraw = false)

    # Save 3D figure
    save(joinpath(OUTPUT_DIR, "wave1d_waterfall.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "wave1d_waterfall.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "wave1d_waterfall.pdf"))")

    # Also create a 2D waterfall-style plot
    fig2 = Figure(size = (800, 480))
    ax2 = Axis(fig2[1, 1],
               xlabel = L"x",
               ylabel = L"u(x, t)" * " (offset by " * L"t" * ")",
               title  = "1D Wave Equation: Waterfall View")

    n_lines = 40
    n_total = length(t_save)
    indices = round.(Int, range(1, n_total, length=n_lines))

    for (i, idx) in enumerate(indices)
        offset = t_save[idx] * 0.15
        c_val = (i - 1) / (n_lines - 1)
        color = Makie.interpolated_getindex(Makie.to_colormap(:viridis), c_val)
        lines!(ax2, x, U_save[idx, :] .+ offset, color = color, linewidth = 0.8)
        band!(ax2, x, fill(offset, length(x)), U_save[idx, :] .+ offset,
              color = (color, 0.1))
    end

    hidespines!(ax2, :t, :r)

    # Save 2D figure
    save(joinpath(OUTPUT_DIR, "wave1d_waterfall_2d.pdf"), fig2)
    save(joinpath(OUTPUT_DIR, "wave1d_waterfall_2d.png"), fig2, px_per_unit = 4)
    println("Additional figure saved to: $(joinpath(OUTPUT_DIR, "wave1d_waterfall_2d.pdf"))")
end

main()
