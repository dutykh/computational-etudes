# heat2d_cheb.jl
#
# Chapter 10: Spectral PDE Solvers with Chebyshev and Fourier Grids
#
# 2D heat equation solver on Chebyshev tensor grid with backward Euler.
#
# PDE: u_t = u_xx + u_yy, -1 < x, y < 1, t > 0
# BCs: u = 0 on boundary (Dirichlet)
# ICs: u(x, y, 0) = exp(-25*((x+0.3)^2 + (y-0.1)^2))
#
# Generates Figure 10.6: Heat diffusion evolution on 2D domain.
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
# Heat equation solver
# -----------------------------------------------------------------------------
"""
    heat2d_cheb(; N=24, tmax=0.5, dt=0.005)

Solve 2D heat equation on Chebyshev tensor grid using backward Euler.

The discrete Laplacian uses Kronecker sum:
    L = D2_int kron I + I kron D2_int

The system is:
    (I - dt * L) vec(U^{n+1}) = vec(U^n)

Returns `(xx, yy, snapshots)` where `snapshots` is a list of `(time, U)` tuples.
"""
function heat2d_cheb(; N=24, tmax=0.5, dt=0.005)
    # Build differentiation matrices
    D, x = cheb_matrix(N)
    D2 = D * D
    xx = [x[i] for i in 1:N+1, j in 1:N+1]
    yy = [x[j] for i in 1:N+1, j in 1:N+1]

    # Interior block
    D2_int = D2[2:N, 2:N]
    I_int = Matrix{Float64}(I, N - 1, N - 1)

    # Kronecker sum Laplacian for interior
    L = kron(D2_int, I_int) + kron(I_int, D2_int)
    A = Matrix{Float64}(I, (N - 1)^2, (N - 1)^2) - dt * L

    # Initial condition
    U = exp.(-25.0 .* ((xx .+ 0.3).^2 .+ (yy .- 0.1).^2))
    U[1, :] .= 0.0;   U[end, :] .= 0.0
    U[:, 1] .= 0.0;   U[:, end] .= 0.0

    # Times for snapshots
    t_snap = [0.0, 0.05, 0.15, 0.5]
    snapshots = [(0.0, copy(U))]
    current_snap_idx = 2

    # Time stepping
    nsteps = Int(ceil(tmax / dt))
    t = 0.0

    # Pre-factorize for efficiency
    A_lu = lu(A)

    for n in 1:nsteps
        # Extract interior and solve
        u_int = vec(U[2:N, 2:N])
        u_int_new = A_lu \ u_int

        # Reconstruct
        U_new = zeros(N + 1, N + 1)
        U_new[2:N, 2:N] = reshape(u_int_new, N - 1, N - 1)
        U = U_new

        t += dt

        # Save snapshot
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
    # Solve heat equation
    N = 32
    xx, yy, snapshots = heat2d_cheb(N=N, tmax=0.6, dt=0.002)

    # Create figure with 4 panels
    fig = Figure(size = (800, 720))

    # Find global vmax from initial condition
    vmax = maximum(snapshots[1][2])
    vmin = 0.0
    n_levels = 30

    for (i, (t, U)) in enumerate(snapshots[1:min(4, length(snapshots))])
        row = (i - 1) ÷ 2 + 1
        col = (i - 1) % 2 + 1
        ax = Axis(fig[row, col],
                  xlabel = L"x",
                  ylabel = L"y",
                  title  = latexstring("t = $(round(t, digits=3))"),
                  aspect = DataAspect())

        # Contourf plot
        contourf!(ax, xx[:, 1], yy[1, :], U,
                  levels = range(vmin, vmax, length=n_levels),
                  colormap = :hot, extendhigh = :auto)

        # Add contour lines
        contour!(ax, xx[:, 1], yy[1, :], U,
                 levels = 8, color = (:black, 0.5), linewidth = 0.3)

        # Mark the initial center on the first panel
        if i == 1
            scatter!(ax, [-0.3], [0.1], marker = :cross, color = :white,
                     markersize = 12)
        end
    end

    Label(fig[0, :], "2D Heat Equation: Diffusion on a Square", fontsize = 14)

    # Save figure
    save(joinpath(OUTPUT_DIR, "heat2d_snapshots.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "heat2d_snapshots.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "heat2d_snapshots.pdf"))")
end

main()
