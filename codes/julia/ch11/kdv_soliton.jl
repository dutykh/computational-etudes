# kdv_soliton.jl
#
# Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs
#
# Fourier pseudospectral solver for the KdV equation with integrating factor RK4:
#     u_t + 6 u u_x + u_xxx = 0 on [0, L), periodic.
#
# Simulates the elastic collision of two solitons of different amplitudes.
# Soliton solution: u(x,t) = 2*a^2 * sech^2(a*(x - x0 - 4*a^2*t))
#
# Uses the cumulative integrating factor approach (purely imaginary linear symbol).
#
# Generates figure: kdv_snapshots.pdf (snapshots before, during, after collision)
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
# Soliton profile
# -----------------------------------------------------------------------------
"""
    soliton(x, a, x0)

KdV soliton: 2*a^2 * sech^2(a*(x - x0)).
"""
function soliton(x::AbstractVector, a::Float64, x0::Float64)
    return 2.0 .* a^2 ./ cosh.(a .* (x .- x0)).^2
end

# -----------------------------------------------------------------------------
# KdV solver (two-soliton collision)
# -----------------------------------------------------------------------------
"""
    kdv_ifrk4(; N=512, L=80.0, a1=1.0, a2=0.5, x1=20.0, x2=40.0,
                t_final=25.0, dt=0.005, dealias=true)

Integrating factor RK4 solver for KdV:
    u_t + 6 u u_x + u_xxx = 0
on [0, L) with periodic BCs.

Uses cumulative integrating factor (purely imaginary linear symbol lam = i*k^3).

Returns `(x, snap_times, snapshots)` where `snapshots` is a matrix with
snapshots as rows.
"""
function kdv_ifrk4(; N::Int=512, L::Float64=80.0,
                     a1::Float64=1.0, a2::Float64=0.5,
                     x1::Float64=20.0, x2::Float64=40.0,
                     t_final::Float64=25.0, dt::Float64=0.005,
                     dealias::Bool=true)
    # Grid and wavenumbers
    x = L .* collect(0:N-1) ./ N
    k_base = collect(0:N-1)
    k_freq = [i <= N÷2 ? i : i - N for i in k_base]
    k = (2.0 * pi / L) .* Float64.(k_freq)

    # Linear symbol
    lam = 1im .* k.^3

    # Dealiasing mask
    cutoff = N ÷ 3
    mask = trues(N)
    if dealias && cutoff > 0
        mask[cutoff+2:N-cutoff] .= false
    end

    # Initial condition: two solitons
    u0 = soliton(x, a1, x1) .+ soliton(x, a2, x2)
    v_hat = fft(u0) |> copy

    function G(t_now, v_hat_now)
        """RHS in integrating factor variables."""
        E = exp.(lam .* t_now)
        u_hat = E .* v_hat_now
        u  = real.(ifft(u_hat))
        ux = real.(ifft(1im .* k .* u_hat))
        f_nl = -6.0 .* u .* ux
        f_nl_hat = fft(f_nl)
        f_nl_hat[.!mask] .= 0.0
        return exp.(-lam .* t_now) .* f_nl_hat
    end

    n_steps = Int(ceil(t_final / dt))
    dt = t_final / n_steps

    n_save = 200
    save_every = max(1, n_steps ÷ n_save)
    snapshots  = [copy(u0)]
    snap_times = [0.0]

    t = 0.0
    for n in 1:n_steps
        k1 = G(t, v_hat)
        k2 = G(t + 0.5 * dt, v_hat .+ 0.5 .* dt .* k1)
        k3 = G(t + 0.5 * dt, v_hat .+ 0.5 .* dt .* k2)
        k4 = G(t + dt, v_hat .+ dt .* k3)
        v_hat .= v_hat .+ (dt / 6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
        t += dt

        if n % save_every == 0
            u_hat = exp.(lam .* t) .* v_hat
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
    N = 512
    L = 80.0
    a1, a2 = 1.0, 0.5
    x1, x2 = 20.0, 40.0   # Taller soliton behind
    t_final = 25.0

    x, t_save, U_save = kdv_ifrk4(
        N=N, L=L, a1=a1, a2=a2, x1=x1, x2=x2, t_final=t_final, dt=0.002
    )

    # -------------------------------------------------------------------------
    # Figure: snapshots before, during, and after collision
    # -------------------------------------------------------------------------
    fig = Figure(size = (900, 700))

    t_targets = [0.0, 8.0, 12.0, t_final]
    titles = ["Initial condition (\$t = 0\$)",
              "Before collision",
              "During collision",
              "After collision (\$t = $(round(Int, t_final))\$)"]
    colors_snap = [NAVY, TEAL, CORAL, NAVY]

    for (i, (t_target, title_str, col)) in enumerate(zip(t_targets, titles, colors_snap))
        row = (i - 1) ÷ 2 + 1
        column = (i - 1) % 2 + 1
        ax = Axis(fig[row, column],
                  xlabel = L"x",
                  ylabel = L"u(x, t)",
                  title  = title_str)

        idx = argmin(abs.(t_save .- t_target))
        lines!(ax, x, U_save[idx, :], color = col, linewidth = 1.5)

        xlims!(ax, 0, L)
        ylims!(ax, -0.2, 2.5)
        hidespines!(ax, :t, :r)
    end

    Label(fig[0, :],
          "KdV Two-Soliton Collision (\$a_1 = $a1\$, \$a_2 = $a2\$, \$N = $N\$)",
          fontsize = 13)

    save(joinpath(OUTPUT_DIR, "kdv_snapshots.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "kdv_snapshots.png"), fig, px_per_unit = 4)
    println("Figure saved: $(joinpath(OUTPUT_DIR, "kdv_snapshots.pdf"))")
end

main()
