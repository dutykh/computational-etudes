# kdv_zabusky_kruskal.jl
#
# Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs
#
# Fourier pseudospectral solver for the KdV equation with integrating factor RK4:
#     u_t + u u_x + delta^2 u_xxx = 0 on [0, 2), periodic.
#
# Replicates the Zabusky-Kruskal (1965) experiment: a cosine initial condition
# breaks into a train of solitons, demonstrating the FPU recurrence phenomenon.
#
# Initial condition: u(x, 0) = cos(pi * x)
#
# Uses the cumulative integrating factor approach (purely imaginary linear symbol).
#
# Generates figure: kdv_zabusky_kruskal.pdf (space-time diagram)
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
# KdV solver (Zabusky-Kruskal form)
# -----------------------------------------------------------------------------
"""
    kdv_zabusky_kruskal(; N=256, L=2.0, t_final=3.6/pi, dt=1e-4)

Solve u_t + u u_x + delta^2 u_xxx = 0 on [0, L), periodic,
with u(x, 0) = cos(pi * x) and delta^2 = 0.022^2 (Zabusky-Kruskal form).

Uses cumulative integrating factor RK4 with 2/3 dealiasing.

Returns `(x, snap_times, snapshots)` where `snapshots` is a matrix with
snapshots as rows.
"""
function kdv_zabusky_kruskal(; N::Int=256, L::Float64=2.0,
                               t_final::Float64=3.6/pi, dt::Float64=1e-4)
    delta2 = 0.022^2

    # Grid and wavenumbers
    x = L .* collect(0:N-1) ./ N
    k_base = collect(0:N-1)
    k_freq = [i <= N÷2 ? i : i - N for i in k_base]
    k = (2.0 * pi / L) .* Float64.(k_freq)

    # Linear symbol: lambda_k = -i * delta^2 * k^3
    lam = -1im .* delta2 .* k.^3

    # 2/3 dealiasing mask
    cutoff = N ÷ 3
    mask = trues(N)
    if cutoff > 0
        mask[cutoff+2:N-cutoff] .= false
    end

    # Initial condition
    u0 = cos.(pi .* x)
    v_hat = fft(u0) |> copy

    function G(t_now, v_hat_now)
        """RHS in integrating factor variables."""
        E = exp.(lam .* t_now)
        u_hat = E .* v_hat_now
        u  = real.(ifft(u_hat))
        ux = real.(ifft(1im .* k .* u_hat))
        f_nl = -u .* ux
        f_nl_hat = fft(f_nl)
        f_nl_hat[.!mask] .= 0.0
        return exp.(-lam .* t_now) .* f_nl_hat
    end

    # Time stepping
    n_steps = Int(ceil(t_final / dt))
    dt = t_final / n_steps

    n_save = 300
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
    N = 256
    t_final = 3.6 / pi

    x, t_save, U_save = kdv_zabusky_kruskal(N=N, t_final=t_final, dt=5e-5)

    # -------------------------------------------------------------------------
    # Figure: Space-time diagram + snapshots
    # -------------------------------------------------------------------------
    fig = Figure(size = (960, 450))

    # Left: space-time plot
    ax1 = Axis(fig[1, 1],
               xlabel = L"x",
               ylabel = L"t",
               title  = "Zabusky--Kruskal Experiment")

    hm = heatmap!(ax1, x, t_save, U_save, colormap = :RdBu, interpolate = true)
    Colorbar(fig[1, 1, Right()], hm, label = L"u(x, t)")

    # Right: snapshots at selected times
    ax2 = Axis(fig[1, 2],
               xlabel = L"x",
               ylabel = L"u(x, t)",
               title  = "Soliton Train Formation")

    t_targets = [0.0, 0.3, 0.6, 0.9, t_final]
    colors_snap = [SKY, TEAL, colorant"#2ECC71", colorant"#F39C12", NAVY]

    for (i, t_target) in enumerate(t_targets)
        idx = argmin(abs.(t_save .- t_target))
        lines!(ax2, x, U_save[idx, :], color = colors_snap[i], linewidth = 1.3,
               label = latexstring("t = $(round(t_save[idx], digits=2))"))
    end

    xlims!(ax2, 0, 2.0)
    axislegend(ax2, position = :rt, framevisible = true, labelsize = 9)
    hidespines!(ax2, :t, :r)

    save(joinpath(OUTPUT_DIR, "kdv_zabusky_kruskal.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "kdv_zabusky_kruskal.png"), fig, px_per_unit = 4)
    println("Figure saved: $(joinpath(OUTPUT_DIR, "kdv_zabusky_kruskal.pdf"))")
end

main()
