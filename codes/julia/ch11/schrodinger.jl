# schrodinger.jl
#
# Chapter 10: Spectral PDE Solvers with Chebyshev and Fourier Grids
#
# Time-dependent Schrodinger equation with harmonic potential (Strang splitting).
#
# PDE: i u_t = -u_xx + x^2 u,  0 < x < 2*pi,  t > 0 (periodic)
# ICs: u(x, 0) = exp(-5*(x - pi)^2) * exp(3ix)  (Gaussian wavepacket)
#
# Method: Strang splitting (second-order, unitary)
#   1. Half potential step (physical space): u -> u * exp(-i x^2 dt/2)
#   2. Full kinetic step (Fourier space):   u_hat -> u_hat * exp(-i k^2 dt)
#   3. Half potential step (physical space): u -> u * exp(-i x^2 dt/2)
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
# Schrodinger equation solver (Strang splitting)
# -----------------------------------------------------------------------------
"""
    schrodinger_strang(; N=256, tmax=3.0, n_snapshots=200)

Solve the time-dependent Schrodinger equation with harmonic potential
using Strang splitting and FFT.

The Strang splitting preserves the L^2 norm (unitarity) and is
second-order accurate in time.

Returns `(x, t_save, U_save)` where `U_save` is a matrix of |u(x,t)|^2 snapshots.
"""
function schrodinger_strang(; N=256, tmax=3.0, n_snapshots=200)
    # Grid
    h = 2pi / N
    x = h .* collect(0:N-1)

    # Wavenumbers for FFT
    if N % 2 == 0
        k = [collect(0:N÷2-1); 0; collect(1-N÷2:-1)]
    else
        k = [collect(0:(N-1)÷2); collect(-(N-1)÷2:-1)]
    end

    # Potential: V(x) = x^2
    V = x.^2

    # Initial condition: Gaussian wavepacket with momentum
    u = exp.(-5.0 .* (x .- pi).^2) .* exp.(3im .* x)

    # Time stepping
    dt = 0.5 * h
    nsteps = Int(ceil(tmax / dt))
    dt = tmax / nsteps  # Adjust to hit tmax exactly

    # Precompute exponentials for splitting operators
    half_potential = exp.(-0.5im .* V .* dt)    # Half potential step
    full_kinetic  = exp.(-1im .* k.^2 .* dt)   # Full kinetic step

    # Store snapshots
    save_every = max(1, nsteps ÷ n_snapshots)
    U_save = [abs.(u).^2]
    t_save = [0.0]

    # Record initial L^2 norm for unitarity check
    norm0 = sqrt(h * sum(abs.(u).^2))

    t = 0.0
    for n in 1:nsteps
        # Strang splitting: potential/2 -> kinetic -> potential/2
        u .= half_potential .* u                   # Half potential step
        u_hat = fft(u)
        u_hat .= full_kinetic .* u_hat             # Full kinetic step
        u .= ifft(u_hat)
        u .= half_potential .* u                   # Half potential step

        t += dt

        # Save snapshot
        if n % save_every == 0
            push!(U_save, abs.(u).^2)
            push!(t_save, t)
        end
    end

    # Check unitarity preservation
    norm_final = sqrt(h * sum(abs.(u).^2))
    println("Initial L2 norm: $(round(norm0, digits=10))")
    println("Final   L2 norm: $(round(norm_final, digits=10))")
    println("Relative change: $(round(abs(norm_final - norm0) / norm0, sigdigits=2))")

    # Convert to matrix form (snapshots as rows)
    U_mat = reduce(vcat, [u' for u in U_save])
    return x, t_save, U_mat
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Solve Schrodinger equation
    N = 256
    tmax = 3.0
    x, t_save, U_save = schrodinger_strang(N=N, tmax=tmax, n_snapshots=200)

    # Create figure with two panels
    fig = Figure(size = (960, 400))

    # Left panel: space-time plot of |u|^2
    ax_main = Axis(fig[1, 1],
                   xlabel = L"x",
                   ylabel = L"t",
                   title  = L"Probability Density $|\psi(x,t)|^2$")

    hm = heatmap!(ax_main, x, t_save, U_save, colormap = :inferno,
                  interpolate = true)
    Colorbar(fig[1, 1, Right()], hm, label = L"|\psi|^2")

    # Right panel: waterfall view
    ax_water = Axis(fig[1, 2],
                    xlabel = L"x",
                    ylabel = L"|\psi(x,t)|^2" * " (offset by " * L"t" * ")",
                    title  = L"Waterfall View: $|\psi|^2$ Evolution")

    n_lines = 40
    n_total = length(t_save)
    indices = round.(Int, range(1, n_total, length=n_lines))

    for (i, idx) in enumerate(indices)
        offset = t_save[idx] * 0.3
        c_val = (i - 1) / (n_lines - 1)
        color = Makie.interpolated_getindex(Makie.to_colormap(:viridis), c_val)
        lines!(ax_water, x, U_save[idx, :] .+ offset, color = color, linewidth = 0.8)
    end

    # Add inset showing the potential V(x) = x^2
    ax_inset = Axis(fig[1, 2],
                    width = Relative(0.35), height = Relative(0.25),
                    halign = 0.95, valign = 0.95,
                    xlabel = L"x", ylabel = L"V(x)",
                    xlabelsize = 8, ylabelsize = 8,
                    xticklabelsize = 7, yticklabelsize = 7,
                    titlesize = 9, title = "Harmonic Potential",
                    backgroundcolor = :white)
    V = x.^2
    lines!(ax_inset, x, V, color = CORAL, linewidth = 1.5)

    hidespines!(ax_water, :t, :r)

    # Save figure
    save(joinpath(OUTPUT_DIR, "schrodinger.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "schrodinger.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "schrodinger.pdf"))")
end

main()
