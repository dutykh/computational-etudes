#!/usr/bin/env julia
#=
catastrophe_gallery.jl
Chapter 17: Stability of Spectral Methods

Computational Etude 17.1: Gallery of Deliberate Instabilities.

Three-panel figure showing what happens when CFL / stability limits are
deliberately violated in spectral discretisations:

  (a) Fourier advection u_t + u_x = 0 with forward Euler, dt 20% above CFL.
  (b) Chebyshev wave equation u_tt = u_xx (Dirichlet) with leapfrog,
      dt 50% above CFL.
  (c) Chebyshev heat equation u_t = u_xx (Dirichlet) with forward Euler,
      dt doubled above the stability limit.

Generates Figure 17.1: catastrophe_gallery.pdf

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
=#

using LinearAlgebra
using FFTW
using Plots
using LaTeXStrings

gr()

const NAVY  = RGB(20/255, 45/255, 110/255)
const SKY   = RGB(120/255, 150/255, 210/255)
const CORAL = RGB(231/255, 76/255, 60/255)
const TEAL  = RGB(22/255, 160/255, 133/255)
const AMBER = RGB(230/255, 126/255, 34/255)

const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch17", "julia")
mkpath(OUTPUT_DIR)

"""
Chebyshev differentiation matrix of size (N+1) x (N+1).
"""
function cheb(N)
    if N == 0; return [0.0;;], [1.0]; end
    x = [cos(pi * j / N) for j in 0:N]
    c = [j == 0 || j == N ? 2.0 : 1.0 for j in 0:N] .* [(-1.0)^j for j in 0:N]
    X = repeat(x, 1, N+1)
    dX = X - X' + I(N+1)
    D = (c * (1.0 ./ c)') ./ dX
    D .-= Diagonal(vec(sum(D, dims=2)))
    return D, x
end

"""
Forward Euler on Fourier advection, dt above CFL.
"""
function fourier_advection_instability()
    N = 64
    L = 2.0 * pi
    x = range(0, L, length=N+1)[1:N]  # endpoint=false

    # Wavenumbers
    k = fftfreq(N, N / 1.0)
    k .= 2.0 * pi .* k ./ L

    # Initial condition: Gaussian bump
    u0 = exp.(-10.0 .* (x .- pi) .^ 2)
    dt_stable = 0.5 * L / (N * 1.0)  # reference "CFL" ~ dx/c
    dt = 1.2 * dt_stable              # 20% above

    snap_times = [0, 50, 100, 150]

    # "Stable" run with RK4 for reference
    ik = 1im .* k
    u_hat = fft(u0)
    snapshots_stable = Vector{Vector{Float64}}()
    for n in 0:maximum(snap_times)
        if n in snap_times
            push!(snapshots_stable, real.(ifft(u_hat)))
        end
        # RK4 step
        k1 = -ik .* u_hat
        k2 = -ik .* (u_hat .+ 0.5 .* dt .* k1)
        k3 = -ik .* (u_hat .+ 0.5 .* dt .* k2)
        k4 = -ik .* (u_hat .+ dt .* k3)
        u_hat .= u_hat .+ dt / 6.0 .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
    end

    # Unstable run with forward Euler
    u_hat = fft(u0)
    snapshots_unstable = Vector{Vector{Float64}}()
    for n in 0:maximum(snap_times)
        if n in snap_times
            push!(snapshots_unstable, real.(ifft(u_hat)))
        end
        u_hat .= u_hat .+ dt .* (-ik .* u_hat)
    end

    return collect(x), snapshots_stable, snapshots_unstable, snap_times, dt
end

"""
Leapfrog on Chebyshev wave equation, dt above CFL.
"""
function chebyshev_wave_instability()
    N = 40
    D, x = cheb(N)
    D2 = D * D
    # Interior points (Dirichlet BCs): rows/cols 2:N in 1-based indexing
    D2_int = D2[2:N, 2:N]
    x_int = x[2:N]

    # Critical dt for leapfrog: dt < 2/sqrt(rho(D2_int))
    eigs = eigvals(D2_int)
    rho = maximum(abs.(eigs))
    dt_crit = 2.0 / sqrt(rho)
    dt = 1.5 * dt_crit  # 50% above

    # Initial condition: smooth bump
    u = exp.(-20.0 .* x_int .^ 2)
    u_prev = copy(u)  # zero velocity start

    n_steps = 60
    snap_full = zeros(N + 1)
    snap_full[2:N] .= u
    snapshots = [copy(snap_full)]

    snap_times = [0, 20, 40, 60]

    for n in 1:n_steps
        if n == 1
            # First step: use forward Euler for u^1
            u_new = u .+ 0.5 .* dt^2 .* (D2_int * u)
        else
            u_new = 2.0 .* u .- u_prev .+ dt^2 .* (D2_int * u)
        end
        u_prev .= u
        u .= u_new
        if n in snap_times
            snap_full = zeros(N + 1)
            snap_full[2:N] .= u
            push!(snapshots, copy(snap_full))
        end
    end

    return x, snapshots, snap_times, dt
end

"""
Forward Euler on Chebyshev heat equation, dt doubled.
"""
function chebyshev_heat_instability()
    N = 20
    D, x = cheb(N)
    D2 = D * D
    D2_int = D2[2:N, 2:N]
    x_int = x[2:N]

    eigs = eigvals(D2_int)
    rho = maximum(abs.(real.(eigs)))
    dt_crit = 2.0 / rho
    dt = 2.0 * dt_crit  # doubled

    u = sin.(pi .* x_int)
    snapshots = Vector{Vector{Float64}}()
    snap_times = [0, 1, 2, 3]

    for n in 0:maximum(snap_times)
        if n in snap_times
            snap_full = zeros(N + 1)
            snap_full[2:N] .= u
            push!(snapshots, copy(snap_full))
        end
        u .= u .+ dt .* (D2_int * u)
    end

    return x, snapshots, snap_times, dt
end

function main()
    x_f, snaps_s, snaps_u, times_f, dt_f = fourier_advection_instability()
    x_w, snaps_w, times_w, dt_w = chebyshev_wave_instability()
    x_h, snaps_h, times_h, dt_h = chebyshev_heat_instability()

    # --- Panel (a): Fourier advection ---
    p1 = plot(xlabel=L"$x$", ylabel=L"$u(x,t)$",
              title="(a) Fourier advection, forward Euler",
              ylims=(-3, 3), legend=:topright, legendfontsize=7,
              grid=true, gridalpha=0.3)
    alphas_val = [1.0, 0.5, 0.7, 1.0]
    for i in eachindex(times_f)
        if i == 1
            plot!(p1, x_f, snaps_s[i], color=NAVY, lw=1.5, label="initial")
        else
            plot!(p1, x_f, snaps_s[i], color=SKY, lw=1.0,
                  alpha=alphas_val[i], ls=:dash, label="")
            u_clip = clamp.(snaps_u[i], -3, 3)
            plot!(p1, x_f, u_clip, color=CORAL, lw=1.0,
                  alpha=alphas_val[i],
                  label="FE, n=$(times_f[i])")
        end
    end

    # --- Panel (b): Chebyshev wave ---
    colors_w = [NAVY, SKY, AMBER, CORAL]
    t_labels = vcat([0], times_w[2:end])
    p2 = plot(xlabel=L"$x$", ylabel=L"$u(x,t)$",
              title="(b) Chebyshev wave, leapfrog",
              ylims=(-5, 5), legend=:topleft, legendfontsize=7,
              grid=true, gridalpha=0.3)
    for (i, snap) in enumerate(snaps_w)
        snap_clip = clamp.(snap, -5, 5)
        plot!(p2, x_w, snap_clip, color=colors_w[i], lw=1.2,
              label="n = $(t_labels[i])")
    end

    # --- Panel (c): Chebyshev heat ---
    colors_h = [NAVY, SKY, AMBER, CORAL]
    p3 = plot(xlabel=L"$x$", ylabel=L"$u(x,t)$",
              title="(c) Chebyshev heat, forward Euler",
              ylims=(-5, 5), legend=:topright, legendfontsize=7,
              grid=true, gridalpha=0.3)
    for (i, snap) in enumerate(snaps_h)
        snap_clip = clamp.(snap, -5, 5)
        plot!(p3, x_h, snap_clip, color=colors_h[i], lw=1.2,
              label="n = $(times_h[i])")
    end

    fig = plot(p1, p2, p3, layout=(1, 3), size=(1400, 420))

    savefig(fig, joinpath(OUTPUT_DIR, "catastrophe_gallery.pdf"))
    savefig(fig, joinpath(OUTPUT_DIR, "catastrophe_gallery.png"))
    println("Figure saved to $(joinpath(OUTPUT_DIR, "catastrophe_gallery.pdf"))")

    display(fig)
end

main()
