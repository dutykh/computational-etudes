#!/usr/bin/env julia
#=
fair_comparison.jl
Chapter 17: Stability of Spectral Methods

Computational Etude 17.6: Fair Comparison of Time-Stepping Methods.

Two-panel figure comparing Forward Euler, Backward Euler, Crank-Nicolson,
and Integrating Factor (exact exponential) on the 1D Fourier heat equation
u_t = nu * u_xx on [0, 2*pi) with nu = 0.1.

  (a) Error vs dt (log-log).
  (b) Error vs wall-clock time (log-log).

Initial condition: u0 = sin(x) + 0.5*cos(3x).
Exact solution: each Fourier mode k decays as exp(-nu*k^2*T).

Generates Figure 17.6: fair_comparison.pdf

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
=#

using FFTW
using Plots
using LaTeXStrings
using Printf

gr()

const NAVY  = RGB(20/255, 45/255, 110/255)
const SKY   = RGB(120/255, 150/255, 210/255)
const CORAL = RGB(231/255, 76/255, 60/255)
const TEAL  = RGB(22/255, 160/255, 133/255)
const AMBER = RGB(230/255, 126/255, 34/255)

const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch17", "julia")
mkpath(OUTPUT_DIR)

"""Exact solution in Fourier space: each mode decays exponentially."""
function exact_solution_fourier(N, nu, T)
    x = range(0, 2pi, length=N+1)[1:N]
    u0 = sin.(x) .+ 0.5 .* cos.(3 .* x)
    u0_hat = fft(u0)
    k = fftfreq(N, N / 1.0)
    decay = exp.(-nu .* k.^2 .* T)
    u_exact_hat = u0_hat .* decay
    return real.(ifft(u_exact_hat))
end

"""Forward Euler in Fourier space."""
function forward_euler(N, nu, T, dt)
    x = range(0, 2pi, length=N+1)[1:N]
    u0 = sin.(x) .+ 0.5 .* cos.(3 .* x)
    u_hat = fft(u0)
    k = fftfreq(N, N / 1.0)
    lam = -nu .* k.^2

    n_steps = round(Int, T / dt)
    dt_actual = T / n_steps

    t0 = time_ns()
    for _ in 1:n_steps
        u_hat .= u_hat .+ dt_actual .* lam .* u_hat
        maximum(abs.(u_hat)) > 1e10 && return nothing, (time_ns() - t0) / 1e9
    end
    elapsed = (time_ns() - t0) / 1e9

    return real.(ifft(u_hat)), elapsed
end

"""Backward Euler in Fourier space."""
function backward_euler(N, nu, T, dt)
    x = range(0, 2pi, length=N+1)[1:N]
    u0 = sin.(x) .+ 0.5 .* cos.(3 .* x)
    u_hat = fft(u0)
    k = fftfreq(N, N / 1.0)
    lam = -nu .* k.^2

    n_steps = round(Int, T / dt)
    dt_actual = T / n_steps
    factor = 1.0 ./ (1.0 .- dt_actual .* lam)

    t0 = time_ns()
    for _ in 1:n_steps
        u_hat .= u_hat .* factor
    end
    elapsed = (time_ns() - t0) / 1e9

    return real.(ifft(u_hat)), elapsed
end

"""Crank-Nicolson in Fourier space."""
function crank_nicolson(N, nu, T, dt)
    x = range(0, 2pi, length=N+1)[1:N]
    u0 = sin.(x) .+ 0.5 .* cos.(3 .* x)
    u_hat = fft(u0)
    k = fftfreq(N, N / 1.0)
    lam = -nu .* k.^2

    n_steps = round(Int, T / dt)
    dt_actual = T / n_steps
    factor = (1.0 .+ 0.5 .* dt_actual .* lam) ./ (1.0 .- 0.5 .* dt_actual .* lam)

    t0 = time_ns()
    for _ in 1:n_steps
        u_hat .= u_hat .* factor
    end
    elapsed = (time_ns() - t0) / 1e9

    return real.(ifft(u_hat)), elapsed
end

"""Integrating Factor (exact exponential for linear problems)."""
function integrating_factor(N, nu, T, dt)
    x = range(0, 2pi, length=N+1)[1:N]
    u0 = sin.(x) .+ 0.5 .* cos.(3 .* x)
    u_hat = fft(u0)
    k = fftfreq(N, N / 1.0)
    lam = -nu .* k.^2

    n_steps = round(Int, T / dt)
    dt_actual = T / n_steps
    factor = exp.(lam .* dt_actual)

    t0 = time_ns()
    for _ in 1:n_steps
        u_hat .= u_hat .* factor
    end
    elapsed = (time_ns() - t0) / 1e9

    return real.(ifft(u_hat)), elapsed
end

function main()
    N  = 64
    nu = 0.1
    T  = 1.0

    u_exact = exact_solution_fourier(N, nu, T)

    # CFL for forward Euler: dt <= 2 / (nu * (N/2)^2)
    dt_cfl = 2.0 / (nu * (N / 2)^2)
    @printf("Forward Euler CFL limit: dt <= %.6e\n", dt_cfl)

    methods = [
        ("Forward Euler",      forward_euler,      NAVY,  :circle,   true),
        ("Backward Euler",     backward_euler,      TEAL,  :square,   false),
        ("Crank-Nicolson",     crank_nicolson,      CORAL, :utriangle, false),
        ("Integrating Factor", integrating_factor,  AMBER, :diamond,  false),
    ]

    dt_values_all = 10 .^ range(-5, -1, length=30)
    dt_values_fe  = dt_values_all[dt_values_all .< 0.9 * dt_cfl]

    p1 = plot(xlabel=L"\Delta t", ylabel=L"\|u_{\mathrm{num}} - u_{\mathrm{exact}}\|_\infty",
              title="(a) Error vs time step",
              xscale=:log10, yscale=:log10, ylims=(1e-15, 1e1),
              legend=:bottomright, legendfontsize=7,
              grid=true, gridalpha=0.3, minorgrid=true)

    p2 = plot(xlabel="Wall-clock time (s)",
              ylabel=L"\|u_{\mathrm{num}} - u_{\mathrm{exact}}\|_\infty",
              title="(b) Error vs computational cost",
              xscale=:log10, yscale=:log10, ylims=(1e-15, 1e1),
              legend=:topright, legendfontsize=7,
              grid=true, gridalpha=0.3, minorgrid=true)

    for (name, method, color, marker, is_fe) in methods
        dt_vals = is_fe ? dt_values_fe : dt_values_all
        errors  = Float64[]
        times   = Float64[]
        dts_used = Float64[]

        for dt in dt_vals
            result = method(N, nu, T, dt)
            u_num, elapsed = result
            u_num === nothing && continue
            err = maximum(abs.(u_num .- u_exact))
            err = max(err, 1e-14)
            push!(errors, err)
            push!(times, elapsed)
            push!(dts_used, dt)
        end

        isempty(errors) && continue

        plot!(p1, dts_used, errors, marker=marker, color=color,
              ms=4, lw=1.2, label=name, markevery=3)
        plot!(p2, times, errors, marker=marker, color=color,
              ms=4, lw=1.2, label=name, markevery=3)
    end

    # Reference slopes
    dt_ref = 10 .^ range(-5, -1, length=50)
    plot!(p1, dt_ref, 0.5 .* dt_ref, ls=:dash, color=:gray, lw=0.8, alpha=0.7, label="")
    annotate!(p1, 2e-3, 2e-3, text("slope 1", 7, :gray, rotation=35))
    plot!(p1, dt_ref, 2.0 .* dt_ref.^2, ls=:dot, color=:gray, lw=0.8, alpha=0.7, label="")
    annotate!(p1, 5e-3, 8e-5, text("slope 2", 7, :gray, rotation=42))

    fig = plot(p1, p2, layout=(1, 2), size=(1200, 500))

    savefig(fig, joinpath(OUTPUT_DIR, "fair_comparison.pdf"))
    savefig(fig, joinpath(OUTPUT_DIR, "fair_comparison.png"))
    println("\nFigure saved to $(joinpath(OUTPUT_DIR, "fair_comparison.pdf"))")

    display(fig)
end

main()
