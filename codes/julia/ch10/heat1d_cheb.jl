# heat1d_cheb.jl
#
# Chapter 10: Spectral PDE Solvers with Chebyshev and Fourier Grids
#
# 1D heat equation solver on Chebyshev grid with Crank--Nicolson time stepping.
#
# PDE: u_t = kappa * u_xx, -1 < x < 1, t > 0
# BCs: u(-1, t) = u(1, t) = 0 (Dirichlet)
# ICs: u(x, 0) = exp(-20*(x+0.5)^2) + 0.5*exp(-30*(x-0.4)^2) (two bumps)
#
# Generates Figure 10.5: Heat diffusion evolution.
#
# Author: Dr. Denys Dutykh
#         Mathematics Department
#         Khalifa University of Science and Technology
#         Abu Dhabi, UAE
#
# Part of the book "Computational Etudes: A Spectral Approach"
# https://github.com/dutykh/computational-etudes

using CairoMakie
using LaTeXStrings
using Colors
using LinearAlgebra

# Import cheb_matrix
include(joinpath(@__DIR__, "chebfft.jl"))
include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

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
const PURPLE = colorant"#9B59B6"

# Output path
SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch10", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Heat equation solver
# -----------------------------------------------------------------------------
"""
    heat1d_cheb(; N=64, kappa=1.0, tmax=1.0, dt=0.01, n_snapshots=100)

Solve 1D heat equation using Crank--Nicolson time stepping.

The semi-discrete system is:
    (I - kappa*dt/2 * D2_int) u^{n+1} = (I + kappa*dt/2 * D2_int) u^n

Returns `(x, t_save, U_save)` where `U_save` is a matrix with snapshots as rows.
"""
function heat1d_cheb(; N=64, kappa=1.0, tmax=1.0, dt=0.01, n_snapshots=100)
    # Build differentiation matrices
    D, x = cheb_matrix(N)
    D2 = D * D

    # Interior indices (Dirichlet BCs: u[1] = u[N+1] = 0, Julia 1-based)
    ii = 2:N
    D2_int = D2[ii, ii]
    I_int = Matrix{Float64}(I, N - 1, N - 1)

    # Crank--Nicolson matrices
    A = I_int - 0.5 * kappa * dt * D2_int
    B = I_int + 0.5 * kappa * dt * D2_int

    # Pre-factorize A (constant across all time steps)
    A_lu = lu(A)

    # Initial condition: two bumps
    u_full = exp.(-20.0 .* (x .+ 0.5).^2) .+ 0.5 .* exp.(-30.0 .* (x .- 0.4).^2)
    u = u_full[ii]

    # Time stepping
    nsteps = Int(ceil(tmax / dt))
    save_every = max(1, nsteps ÷ n_snapshots)

    U_save = [copy(u_full)]
    t_save = [0.0]

    t = 0.0
    for n in 1:nsteps
        # Solve the CN system using pre-factorized A
        u = A_lu \ (B * u)
        t += dt

        if n % save_every == 0
            # Reconstruct full solution
            u_full = zeros(N + 1)
            u_full[ii] = u
            push!(U_save, copy(u_full))
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
    # Solve heat equation
    N = 64
    tmax = 0.5
    x, t_save, U_save = heat1d_cheb(N=N, tmax=tmax, dt=0.005, n_snapshots=80)

    # Create figure with two panels
    fig = Figure(size = (960, 480))

    # Left panel: 3D surface
    ax_main = Axis3(fig[1, 1],
                    xlabel = L"x",
                    ylabel = L"t",
                    zlabel = L"u(x, t)",
                    title  = "Heat Diffusion: Space-Time View",
                    azimuth = -60 * pi / 180,
                    elevation = 25 * pi / 180)

    surface!(ax_main, x, t_save, U_save', colormap = :hot)

    # Right panel: snapshots at selected times
    ax_snap = Axis(fig[1, 2],
                   xlabel = L"x",
                   ylabel = L"u(x, t)",
                   title  = "Heat Diffusion: Time Snapshots")

    # Select specific times for snapshots
    t_targets = [0.0, 0.05, 0.1, 0.2, 0.5]
    colors = [NAVY, SKY, TEAL, CORAL, PURPLE]

    for (t_target, color) in zip(t_targets, colors)
        idx = argmin(abs.(t_save .- t_target))
        lines!(ax_snap, x, U_save[idx, :], color = color, linewidth = 1.5,
               label = latexstring("t = $(round(t_save[idx], digits=2))"))
    end

    xlims!(ax_snap, -1, 1)
    ylims!(ax_snap, 0, 1.1)
    hidespines!(ax_snap, :t, :r)
    axislegend(ax_snap, position = :rt, labelsize = 10)

    # Save figure
    save(joinpath(OUTPUT_DIR, "heat1d_evolution.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "heat1d_evolution.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "heat1d_evolution.pdf"))")
end

main()
