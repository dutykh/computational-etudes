# nls_recurrence.jl
#
# Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs
#
# Split-step Fourier solver for the focusing nonlinear Schrodinger equation:
#     i u_t + u_xx + 2 |u|^2 u = 0 on [0, 2*pi), periodic.
#
# Demonstrates modulation instability and approximate recurrence using
# Strang splitting (linear-nonlinear-linear half-steps).
#
# Initial condition: u(x, 0) = 1 + 0.1 * cos(x)
#
# Generates figure: nls_recurrence.pdf
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
const NAVY  = colorant"#142D6E"
const SKY   = colorant"#7896D2"
const CORAL = colorant"#E74C3C"
const TEAL  = colorant"#16A085"

# Output path
SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch11", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# NLS split-step solver
# -----------------------------------------------------------------------------
"""
    nls_split_step(; N=256, t_final=5.0, dt=1e-3)

Solve i u_t + u_xx + 2|u|^2 u = 0 on [0, 2*pi) with periodic BCs
using the Strang split-step Fourier method.

Linear step:    u_hat_k <- exp(-i k^2 dt/2) * u_hat_k
Nonlinear step: u <- exp(2i |u|^2 dt) * u

Returns `(x, t_save, U_save, L2_save)` where `U_save` is a matrix of |u(x,t)|
snapshots (rows) and `L2_save` tracks the L^2 norm.
"""
function nls_split_step(; N::Int=256, t_final::Float64=5.0, dt::Float64=1e-3)
    x = 2.0 * pi .* collect(0:N-1) ./ N
    if N % 2 == 0
        k = Float64.([collect(0:N÷2-1); 0; collect(1-N÷2:-1)])
    else
        k = Float64.([collect(0:(N-1)÷2); collect(-(N-1)÷2:-1)])
    end

    # Initial condition: modulated plane wave
    u = ComplexF64.(1.0 .+ 0.1 .* cos.(x))

    # Time stepping
    n_steps = Int(ceil(t_final / dt))
    dt = t_final / n_steps

    # Precompute linear propagator (half-step)
    linear_half = exp.(-1im .* k.^2 .* dt ./ 2.0)

    n_save = 400
    save_every = max(1, n_steps ÷ n_save)
    U_save  = [abs.(u)]
    L2_save = [sqrt(mean_abs2(u))]
    t_save  = [0.0]

    t = 0.0
    for n in 1:n_steps
        # Strang splitting: L/2 -> N -> L/2
        u_hat = fft(u)
        u_hat .*= linear_half
        u .= ifft(u_hat)

        # Full nonlinear step
        u .*= exp.(2im .* abs.(u).^2 .* dt)

        # Second linear half-step
        u_hat = fft(u)
        u_hat .*= linear_half
        u .= ifft(u_hat)

        t += dt

        if n % save_every == 0
            push!(U_save, abs.(u))
            push!(L2_save, sqrt(mean_abs2(u)))
            push!(t_save, t)
        end
    end

    # Convert to matrix form
    U_mat = reduce(vcat, [u' for u in U_save])
    return x, t_save, U_mat, L2_save
end

"""Helper: mean of |u|^2."""
mean_abs2(u) = sum(abs.(u).^2) / length(u)

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
function main()
    N = 256
    t_final = 5.0
    dt = 5e-4

    x, t_save, U_save, L2_save = nls_split_step(N=N, t_final=t_final, dt=dt)

    # -------------------------------------------------------------------------
    # Figure: space-time diagram + L2 norm conservation
    # -------------------------------------------------------------------------
    fig = Figure(size = (960, 450),
                 figure_padding = (10, 10, 10, 10))

    # Left: space-time diagram of |u|
    ax1 = Axis(fig[1, 1],
               xlabel = L"x",
               ylabel = L"t",
               title  = L"NLS: $|u(x, t)|$ (Modulation Instability)")

    hm = heatmap!(ax1, x, t_save, U_save, colormap = :inferno, interpolate = true)
    Colorbar(fig[1, 1, Right()], hm, label = L"|u|")

    # Right: L2 norm conservation
    L2_0 = L2_save[1]
    rel_change = abs.(L2_save .- L2_0) ./ L2_0

    ax2 = Axis(fig[1, 2],
               xlabel = L"t",
               ylabel = L"|\|u\|_2 - \|u_0\|_2| / \|u_0\|_2",
               title  = L"$L^2$ Norm Conservation",
               yscale = log10)

    lines!(ax2, t_save, rel_change, color = NAVY, linewidth = 1.0)

    hidespines!(ax2, :t, :r)

    save(joinpath(OUTPUT_DIR, "nls_recurrence.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "nls_recurrence.png"), fig, px_per_unit = 4)
    println("Figure saved: $(joinpath(OUTPUT_DIR, "nls_recurrence.pdf"))")
end

main()
