# allen_cahn.jl
#
# Chapter 10: Spectral PDE Solvers with Chebyshev and Fourier Grids
#
# Allen-Cahn equation with periodic BCs using Fourier spectral IMEX scheme.
#
# PDE: u_t = u_xx + u(1 - u^2), 0 < x < 2*pi, t > 0 (periodic)
#
# The diffusion term u_xx is treated implicitly and the nonlinear reaction
# u(1 - u^2) is treated explicitly (IMEX):
#   (I - dt * D_xx) u^{n+1} = u^n + dt * (u^n - (u^n)^3)
#
# In Fourier space this becomes diagonal:
#   u_hat_k^{n+1} = (u_hat_k^n + dt * FFT(u^n - (u^n)^3)) / (1 + k^2 * dt)
#
# Generates Figure 10.x: Allen-Cahn phase separation via IMEX spectral method.
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
using FFTW
using Random

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
# Allen-Cahn IMEX solver
# -----------------------------------------------------------------------------
"""
    allen_cahn_imex(; N=256, tmax=5.0, dt=0.01, n_snapshots=200)

Solve the Allen-Cahn equation using Fourier spectral IMEX scheme.

Diffusion (u_xx) is treated implicitly in Fourier space.
Reaction (u - u^3) is treated explicitly.

Returns `(x, t_save, U_save)` where `U_save` is a matrix with snapshots as rows.
"""
function allen_cahn_imex(; N=256, tmax=5.0, dt=0.01, n_snapshots=200)
    # Grid
    h = 2pi / N
    x = h .* collect(0:N-1)

    # Wavenumbers
    if N % 2 == 0
        k = [collect(0:N÷2-1); 0; collect(1-N÷2:-1)]
    else
        k = [collect(0:(N-1)÷2); collect(-(N-1)÷2:-1)]
    end

    # Initial condition: smoothed random perturbation
    Random.seed!(42)
    u = 0.1 .* randn(N)

    # Low-pass filter to smooth the initial data
    u_hat = fft(u)
    filter_mask = abs.(k) .<= N ÷ 8
    u_hat .*= filter_mask
    u = real.(ifft(u_hat))

    # IMEX multiplier: 1 / (1 + k^2 * dt)
    imex_denom = 1.0 .+ k.^2 .* dt

    # Time stepping
    nsteps = Int(ceil(tmax / dt))

    # Store snapshots
    save_every = max(1, nsteps ÷ n_snapshots)
    U_save = [copy(u)]
    t_save = [0.0]

    t = 0.0
    for n in 1:nsteps
        # Explicit reaction term: u - u^3 = u(1 - u^2)
        reaction = u .- u.^3

        # Transform to Fourier space
        u_hat = fft(u)
        reaction_hat = fft(reaction)

        # IMEX step in Fourier space:
        # u_hat^{n+1} = (u_hat^n + dt * reaction_hat) / (1 + k^2 * dt)
        u_hat = (u_hat .+ dt .* reaction_hat) ./ imex_denom

        # Transform back to physical space
        u = real.(ifft(u_hat))

        t += dt

        # Save snapshot
        if n % save_every == 0
            push!(U_save, copy(u))
            push!(t_save, t)
        end
    end

    # Convert to matrix form (snapshots as rows)
    U_mat = reduce(vcat, [u' for u in U_save])
    return x, t_save, U_mat
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Solve Allen-Cahn equation
    N = 256
    tmax = 5.0
    dt = 0.01
    x, t_save, U_save = allen_cahn_imex(N=N, tmax=tmax, dt=dt, n_snapshots=200)

    # Create figure with two panels
    fig = Figure(size = (960, 400))

    # Left panel: space-time heatmap
    ax_main = Axis(fig[1, 1],
                   xlabel = L"x",
                   ylabel = L"t",
                   title  = "Allen-Cahn: Phase Separation")

    hm = heatmap!(ax_main, x, t_save, U_save, colormap = :RdBu,
                  colorrange = (-1.0, 1.0), interpolate = true)
    Colorbar(fig[1, 1, Right()], hm, label = L"u(x, t)")

    # Right panel: snapshots at selected times
    ax_snap = Axis(fig[1, 2],
                   xlabel = L"x",
                   ylabel = L"u(x, t)",
                   title  = "Snapshots at Selected Times")

    snapshot_times = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0]
    colors_snap = [SKY, TEAL, colorant"#2ECC71", colorant"#F39C12", CORAL, NAVY]

    for (i, t_target) in enumerate(snapshot_times)
        idx = argmin(abs.(t_save .- t_target))
        t_actual = t_save[idx]
        lines!(ax_snap, x, U_save[idx, :], color = colors_snap[i], linewidth = 1.5,
               label = latexstring("t = $(round(t_actual, digits=1))"))
    end

    xlims!(ax_snap, 0, 2pi)
    ylims!(ax_snap, -1.3, 1.3)
    hlines!(ax_snap, [1.0, -1.0], color = :gray, linestyle = :dash,
            linewidth = 0.5, alpha = 0.5)
    axislegend(ax_snap, position = :rt, framevisible = true, labelsize = 9)
    hidespines!(ax_snap, :t, :r)

    # Save figure
    save(joinpath(OUTPUT_DIR, "allen_cahn.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "allen_cahn.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "allen_cahn.pdf"))")
end

main()
