# wave2d_cheb.jl
#
# Chapter 10: Spectral PDE Solvers with Chebyshev and Fourier Grids
#
# 2D wave equation solver on Chebyshev tensor grid with leapfrog time stepping.
#
# PDE: u_tt = c^2 (u_xx + u_yy), -1 < x, y < 1, t > 0
# BCs: u = 0 on boundary (Dirichlet)
# ICs: u(x, y, 0) = exp(-30*((x-0.2)^2 + (y+0.3)^2)), u_t(x, y, 0) = 0
#
# Generates Figure 10.4: Snapshots of wave evolution on 2D domain.
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

# Import cheb_matrix
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
const NAVY = colorant"#142D6E"

# Output path
SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch10", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Wave equation solver
# -----------------------------------------------------------------------------
"""
    laplacian_chebfft(U, N)

Compute the 2D Laplacian using chebfft along rows and columns.
"""
function laplacian_chebfft(U, N)
    Lap = zeros(N + 1, N + 1)
    for i in 1:N+1
        Lap[i, :] .+= chebfft(chebfft(U[i, :]))   # d²/dx²
    end
    for j in 1:N+1
        Lap[:, j] .+= chebfft(chebfft(U[:, j]))   # d²/dy²
    end
    return Lap
end

"""
    wave2d_cheb(; N=32, c=1.0, tmax=2.0, alpha=3.0)

Solve 2D wave equation on Chebyshev tensor grid using leapfrog time stepping.

Returns `(xx, yy, snapshots)` where `snapshots` is a list of `(time, U)` tuples.
"""
function wave2d_cheb(; N=32, c=1.0, tmax=2.0, alpha=3.0)
    # Chebyshev grid
    x = [cos(pi * j / N) for j in 0:N]
    xx = [x[i] for i in 1:N+1, j in 1:N+1]
    yy = [x[j] for i in 1:N+1, j in 1:N+1]

    # Initial data: offset Gaussian
    U = exp.(-30.0 .* ((xx .- 0.2).^2 .+ (yy .+ 0.3).^2))
    U_prev = copy(U)

    # Time stepping
    dt = alpha / N^2
    nsteps = Int(ceil(tmax / dt))

    # Times for snapshots
    t_snap = [0.0, 0.3, 0.7, 1.2]
    snapshots = [(0.0, copy(U))]
    current_snap_idx = 2  # next target index (1-based)

    # Taylor start
    Lap_U = laplacian_chebfft(U, N)
    U_prev = U .- 0.5 * (c * dt)^2 .* Lap_U

    t = 0.0
    for n in 1:nsteps
        # 2D Laplacian via chebfft along rows and columns
        Lap_U = laplacian_chebfft(U, N)

        # Leapfrog step
        U_new = 2.0 .* U .- U_prev .+ (c * dt)^2 .* Lap_U

        # Boundary conditions
        U_new[1, :]   .= 0.0
        U_new[end, :] .= 0.0
        U_new[:, 1]   .= 0.0
        U_new[:, end] .= 0.0

        U_prev = U
        U = U_new
        t += dt

        # Save snapshot if we are at a target time
        if current_snap_idx <= length(t_snap) && t >= t_snap[current_snap_idx]
            push!(snapshots, (t, copy(U)))
            current_snap_idx += 1
        end
    end

    return xx, yy, snapshots
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Solve wave equation
    N = 40
    xx, yy, snapshots = wave2d_cheb(N=N, tmax=1.5)

    # Create figure with 4 panels
    fig = Figure(size = (800, 720))

    vmin, vmax = -0.5, 1.0
    n_levels = 30

    for (i, (t, U)) in enumerate(snapshots[1:min(4, length(snapshots))])
        row = (i - 1) ÷ 2 + 1
        col = (i - 1) % 2 + 1
        ax = Axis(fig[row, col],
                  xlabel = L"x",
                  ylabel = L"y",
                  title  = latexstring("t = $(round(t, digits=2))"),
                  aspect = DataAspect())

        # Contourf plot
        contourf!(ax, xx[:, 1], yy[1, :], U,
                  levels = range(vmin, vmax, length=n_levels),
                  colormap = :RdBu, extendlow = :auto, extendhigh = :auto)

        # Add contour lines
        contour!(ax, xx[:, 1], yy[1, :], U,
                 levels = 10, color = (:black, 0.5), linewidth = 0.3)

        # Mark the initial center on the first panel
        if i == 1
            scatter!(ax, [0.2], [-0.3], marker = :cross, color = :black,
                     markersize = 12)
        end
    end

    Label(fig[0, :], "2D Wave Equation: Vibrating Membrane", fontsize = 14)

    # Save figure
    save(joinpath(OUTPUT_DIR, "wave2d_snapshots.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "wave2d_snapshots.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "wave2d_snapshots.pdf"))")
end

main()
