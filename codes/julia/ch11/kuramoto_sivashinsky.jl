# kuramoto_sivashinsky.jl
#
# Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs
#
# Fourier pseudospectral solver for the Kuramoto-Sivashinsky equation:
#     u_t + u u_x + u_xx + u_xxxx = 0 on [0, L), periodic.
#
# Uses the per-step integrating factor RK4 (Trefethen SMIM style) to handle
# the stiff linear terms exactly. The linear symbol is L(k) = k^2 - k^4
# (energy production at low k, damping at high k).
#
# The per-step approach uses E = exp(L*dt/2) and E2 = E^2 applied at each step,
# avoiding the numerical overflow of cumulative integrating factors.
#
# Generates two figures:
#   - ks_spacetime.pdf: space-time "fingerprint" plot
#   - ks_spectrum.pdf: time-averaged energy spectrum
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
using Statistics

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
# Kuramoto-Sivashinsky solver
# -----------------------------------------------------------------------------
"""
    kuramoto_sivashinsky(; N=256, L=32.0*pi, t_final=200.0, dt=0.05)

Solve u_t + u u_x + u_xx + u_xxxx = 0 on [0, L), periodic,
using per-step integrating factor RK4 (Trefethen style) with 2/3 dealiasing.

Linear part: L(k) = k^2 - k^4 (treated exactly via IF)
Nonlinear part: N(u) = -u u_x = -(1/2)(u^2)_x

The per-step IF-RK4 uses E = exp(L*dt/2) and E2 = E^2 applied at each step.

Returns `(x, snap_times, snapshots)` where `snapshots` is a matrix with
snapshots as rows.
"""
function kuramoto_sivashinsky(; N::Int=256, L::Float64=32.0*pi,
                                t_final::Float64=200.0, dt::Float64=0.05)
    # Grid and wavenumbers
    x = L .* collect(0:N-1) ./ N
    k_base = collect(0:N-1)
    k_freq = [i <= N÷2 ? i : i - N for i in k_base]
    k = (2.0 * pi / L) .* Float64.(k_freq)

    # Linear symbol: lambda_k = k^2 - k^4
    lam = k.^2 .- k.^4

    # 2/3 dealiasing mask
    cutoff = N ÷ 3
    mask = trues(N)
    if cutoff > 0
        mask[cutoff+2:N-cutoff] .= false
    end

    # Initial condition: small random perturbation
    Random.seed!(42)
    u0 = 0.01 .* randn(N)
    # Low-pass filter
    u_hat_init = fft(u0)
    u_hat_init[abs.(k) .> 5.0 * (2.0 * pi / L)] .= 0.0
    u0 = real.(ifft(u_hat_init))

    # Per-step integrating factor (Trefethen's approach, avoids overflow)
    u_hat = fft(u0) |> copy

    n_steps = Int(ceil(t_final / dt))
    dt = t_final / n_steps

    E  = exp.(lam .* dt ./ 2.0)    # half-step propagator
    E2 = E.^2                       # full-step propagator

    function Nhat(u_hat_now)
        """Nonlinear term in Fourier space: -ik/2 * FFT(u^2)."""
        u = real.(ifft(u_hat_now))
        f_hat = -0.5im .* k .* fft(u.^2)
        f_hat[.!mask] .= 0.0
        return f_hat
    end

    n_save = 500
    save_every = max(1, n_steps ÷ n_save)
    snapshots  = [copy(u0)]
    snap_times = [0.0]

    t = 0.0
    for n in 1:n_steps
        # IF-RK4 per step (Trefethen SMIM style)
        Na = dt .* Nhat(u_hat)
        Nb = dt .* Nhat(E .* u_hat .+ Na ./ 2.0)
        Nc = dt .* Nhat(E .* u_hat .+ Nb ./ 2.0)
        Nd = dt .* Nhat(E2 .* u_hat .+ E .* Nc)
        u_hat .= E2 .* u_hat .+ (E2 .* Na .+ 2.0 .* E .* (Nb .+ Nc) .+ Nd) ./ 6.0
        t += dt

        if n % save_every == 0
            u = real.(ifft(u_hat))
            push!(snapshots, copy(u))
            push!(snap_times, t)
        end
    end

    # Convert to matrix form
    U_mat = reduce(vcat, [u' for u in snapshots])
    return x, snap_times, U_mat
end

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
function main()
    N = 256
    L = 32.0 * pi
    t_final = 200.0

    x, t_save, U_save = kuramoto_sivashinsky(N=N, L=L, t_final=t_final, dt=0.05)

    # -------------------------------------------------------------------------
    # Figure 1: Space-time "fingerprint"
    # -------------------------------------------------------------------------
    fig1 = Figure(size = (700, 520))

    ax1 = Axis(fig1[1, 1],
               xlabel = L"x",
               ylabel = L"t",
               title  = "Kuramoto--Sivashinsky: Spatiotemporal Chaos (\$L = $(round(Int, L / pi))\\pi\$)")

    # Determine color range from the developed part of the solution
    n_quarter = length(t_save) ÷ 4
    vmax = Float64(quantile(abs.(vec(U_save[n_quarter:end, :])), 0.98))

    hm = heatmap!(ax1, x, t_save, U_save, colormap = :RdBu,
                  colorrange = (-vmax, vmax), interpolate = true)
    Colorbar(fig1[1, 1, Right()], hm, label = L"u(x, t)")

    save(joinpath(OUTPUT_DIR, "ks_spacetime.pdf"), fig1)
    save(joinpath(OUTPUT_DIR, "ks_spacetime.png"), fig1, px_per_unit = 4)
    println("Figure saved: $(joinpath(OUTPUT_DIR, "ks_spacetime.pdf"))")

    # -------------------------------------------------------------------------
    # Figure 2: Time-averaged energy spectrum
    # -------------------------------------------------------------------------
    # Average over the second half (statistically steady state)
    n_start = length(t_save) ÷ 2
    n_total_snap = size(U_save, 1)
    spectra = zeros(N ÷ 2)
    for i in n_start:n_total_snap
        u_hat_snap = fft(U_save[i, :])
        spectra .+= abs.(u_hat_snap[1:N÷2]).^2 ./ N^2
    end
    spectra ./= (n_total_snap - n_start + 1)

    k_pos = (2.0 * pi / L) .* collect(0:N÷2-1)

    fig2 = Figure(size = (600, 420))

    ax2 = Axis(fig2[1, 1],
               xlabel = L"Wavenumber $k$",
               ylabel = L"\langle |\hat{u}_k|^2 \rangle",
               title  = "Time-Averaged Energy Spectrum (KS)",
               yscale = log10)

    lines!(ax2, k_pos[2:end], spectra[2:end], color = NAVY, linewidth = 1.5)

    hidespines!(ax2, :t, :r)

    save(joinpath(OUTPUT_DIR, "ks_spectrum.pdf"), fig2)
    save(joinpath(OUTPUT_DIR, "ks_spectrum.png"), fig2, px_per_unit = 4)
    println("Figure saved: $(joinpath(OUTPUT_DIR, "ks_spectrum.pdf"))")
end

main()
