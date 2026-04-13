#!/usr/bin/env julia
#=
kdv_rk_comparison.jl
Chapter 17: Stability of Spectral Methods

Computational Etude 17.7: Plain RK4 vs IF-RK4 on the KdV Equation.

Two-panel figure:
  (a) Maximum stable dt vs N (log-log) for plain RK4 and IF-RK4 on
      KdV u_t + 6*u*u_x + u_xxx = 0.
  (b) Solutions at t = 0.5 for N = 128 with both methods.

For plain RK4, the linear stability limit is set by the dispersive term
u_xxx whose eigenvalues are i*k^3.  RK4 stability on imaginary axis:
|y| <= 2*sqrt(2), giving dt_crit = 2*sqrt(2) / k_max^3 ~ N^{-3}.

For IF-RK4, the linear part is handled exactly, so the stability limit
is set by the nonlinear term only (much milder).

Generates Figure 17.7: kdv_rk_comparison.pdf

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

"""RK4 amplification factor."""
rk4_amplification(z) = 1.0 .+ z .+ z.^2 ./2 .+ z.^3 ./6 .+ z.^4 ./24

"""Analytical critical dt for RK4 on the linear part of KdV (u_xxx).
dt_crit = 2*sqrt(2) / k_max^3."""
function rk4_dt_crit_linear(N)
    L = 2.0
    k = 2.0 * pi .* fftfreq(N, N / L)
    k_max = maximum(abs.(k))
    return 2.0 * sqrt(2.0) / k_max^3
end

"""Set up wavenumbers and dealiasing mask for KdV on [0, 2)."""
function setup_kdv(N)
    L = 2.0
    x = range(0, L, length=N+1)[1:N]
    k = 2.0 * pi .* fftfreq(N, N / L)
    dealias = ones(N)
    kmax = N ÷ 3
    # Zero out modes above 2/3 rule: indices kmax+1 to N-kmax+1 (1-based)
    dealias[kmax+1 : N-kmax+1] .= 0.0
    return collect(x), k, dealias
end

"""IF-RK4 for KdV: integrating factor removes the linear u_xxx term."""
function kdv_ifrk4(N, dt, T)
    x, k, dealias = setup_kdv(N)
    u_hat = fft(cos.(pi .* x))

    ik  = 1im .* k
    lam = -1im .* k.^3

    n_steps = max(1, round(Int, T / dt))
    dt_actual = T / n_steps

    E  = exp.(lam .* dt_actual ./ 2.0)
    E2 = E.^2

    function NL(v_hat)
        v  = real.(ifft(v_hat))
        vx = real.(ifft(ik .* v_hat))
        return fft(-6.0 .* v .* vx) .* dealias
    end

    for _ in 1:n_steps
        Na = NL(u_hat)
        a  = E .* u_hat .+ 0.5 .* dt_actual .* E .* Na
        Nb = NL(a)
        b  = E .* u_hat .+ 0.5 .* dt_actual .* E .* Nb
        Nc = NL(b)
        c  = E2 .* u_hat .+ dt_actual .* E .* Nc
        Nd = NL(c)
        u_hat .= E2 .* u_hat .+
                 dt_actual / 6.0 .* (E2 .* Na .+ 2 .* E .* (Nb .+ Nc) .+ Nd)

        maximum(abs.(u_hat)) > 1e8 && return nothing, x
    end

    u = real.(ifft(u_hat))
    maximum(abs.(u)) > 100 && return nothing, x
    return u, x
end

"""Plain RK4 for KdV in Fourier space."""
function kdv_rk4(N, dt, T)
    x, k, dealias = setup_kdv(N)
    u_hat = fft(cos.(pi .* x))

    ik  = 1im .* k
    ik3 = 1im .* k.^3

    function rhs(v_hat)
        v  = real.(ifft(v_hat))
        vx = real.(ifft(ik .* v_hat))
        return fft(-6.0 .* v .* vx) .* dealias .- ik3 .* v_hat
    end

    n_steps = max(1, round(Int, T / dt))
    dt_actual = T / n_steps

    for _ in 1:n_steps
        k1 = rhs(u_hat)
        k2 = rhs(u_hat .+ 0.5 .* dt_actual .* k1)
        k3 = rhs(u_hat .+ 0.5 .* dt_actual .* k2)
        k4 = rhs(u_hat .+ dt_actual .* k3)
        u_hat .= u_hat .+ dt_actual / 6.0 .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)

        maximum(abs.(u_hat)) > 1e8 && return nothing, x
    end

    u = real.(ifft(u_hat))
    maximum(abs.(u)) > 100 && return nothing, x
    return u, x
end

"""Binary search for IF-RK4 critical dt (limited by nonlinear term)."""
function find_ifrk4_dt_crit(N; T=0.05, n_bisect=40)
    dt_lo = 1e-6
    dt_hi = 0.1

    # Find unstable upper bound
    for _ in 1:20
        u, _ = kdv_ifrk4(N, dt_hi, T)
        u === nothing && break
        dt_hi *= 2
        dt_hi > 10 && return dt_hi
    end

    # Verify dt_lo is stable
    u, _ = kdv_ifrk4(N, dt_lo, T)
    u === nothing && (dt_lo = 1e-8)

    for _ in 1:n_bisect
        dt_mid = sqrt(dt_lo * dt_hi)
        u, _ = kdv_ifrk4(N, dt_mid, T)
        if u !== nothing
            dt_lo = dt_mid
        else
            dt_hi = dt_mid
        end
        dt_hi / dt_lo < 1.02 && break
    end

    return sqrt(dt_lo * dt_hi)
end

function main()
    # --- Panel (a): dt_crit vs N ---
    N_values = [32, 64, 128, 256]
    dt_rk4_vals   = Float64[]
    dt_ifrk4_vals = Float64[]

    for Nv in N_values
        # RK4: analytical linear stability limit
        dtc_rk4 = rk4_dt_crit_linear(Nv)
        @printf("N = %4d: RK4 dt_crit = %.6e\n", Nv, dtc_rk4)
        push!(dt_rk4_vals, dtc_rk4)

        # IF-RK4: binary search
        dtc_if = find_ifrk4_dt_crit(Nv, T=0.05)
        @printf("         IF-RK4 dt_crit = %.6e\n", dtc_if)
        push!(dt_ifrk4_vals, dtc_if)
    end

    N_arr = Float64.(N_values)

    p1 = plot(N_arr, dt_rk4_vals, xscale=:log10, yscale=:log10,
              marker=:circle, ms=6, lw=1.5, color=NAVY,
              label="Plain RK4",
              xlabel=L"$N$", ylabel=L"\Delta t_{\mathrm{crit}}",
              title=L"(a) Maximum stable $\Delta t$ vs $N$",
              legend=:left, legendfontsize=7,
              grid=true, gridalpha=0.3, minorgrid=true)
    plot!(p1, N_arr, dt_ifrk4_vals, marker=:square, ms=6, lw=1.5, color=CORAL,
          label="IF-RK4")

    # Reference slopes via linear fit in log-log
    N_ref = range(20, 300, length=100)

    log_rk4 = log.(dt_rk4_vals)
    log_N   = log.(N_arr)
    # Simple linear fit: slope and intercept
    slope_rk4 = (length(log_N) * sum(log_N .* log_rk4) - sum(log_N) * sum(log_rk4)) /
                (length(log_N) * sum(log_N.^2) - sum(log_N)^2)
    c_rk4 = exp((sum(log_rk4) - slope_rk4 * sum(log_N)) / length(log_N))
    @printf("\nRK4 slope: %.2f\n", slope_rk4)
    plot!(p1, N_ref, c_rk4 .* N_ref .^ slope_rk4, ls=:dash, lw=1.0, color=SKY,
          alpha=0.7, label=latexstring("\\sim N^{$(round(slope_rk4, digits=1))}"))

    log_if = log.(dt_ifrk4_vals)
    slope_if = (length(log_N) * sum(log_N .* log_if) - sum(log_N) * sum(log_if)) /
               (length(log_N) * sum(log_N.^2) - sum(log_N)^2)
    c_if = exp((sum(log_if) - slope_if * sum(log_N)) / length(log_N))
    @printf("IF-RK4 slope: %.2f\n", slope_if)
    plot!(p1, N_ref, c_if .* N_ref .^ slope_if, ls=:dot, lw=1.0, color=AMBER,
          alpha=0.7, label=latexstring("\\sim N^{$(round(slope_if, digits=1))}"))

    # --- Panel (b): Solutions at t=0.5 for N=128 ---
    N_sol = 128
    T_sol = 0.5

    dt_rk4_safe = 0.8 * rk4_dt_crit_linear(N_sol)
    @printf("\nComputing RK4 solution at t=%.1f with N=%d, dt=%.6e...\n",
            T_sol, N_sol, dt_rk4_safe)
    u_rk4, x_rk4 = kdv_rk4(N_sol, dt_rk4_safe, T_sol)

    dt_if_use = 1e-3
    @printf("Computing IF-RK4 solution at t=%.1f with N=%d, dt=%.6e...\n",
            T_sol, N_sol, dt_if_use)
    u_if, x_if = kdv_ifrk4(N_sol, dt_if_use, T_sol)

    p2 = plot(xlabel=L"$x$", ylabel=L"$u(x,t)$",
              title=latexstring("(b) KdV solution at \$t = $T_sol\$, \$N = $N_sol\$"),
              legend=:topright, legendfontsize=9,
              grid=true, gridalpha=0.3)

    if u_rk4 !== nothing
        plot!(p2, x_rk4, u_rk4, color=NAVY, lw=1.8, label="RK4")
    else
        println("  WARNING: RK4 solution blew up!")
    end
    if u_if !== nothing
        plot!(p2, x_if, u_if, color=CORAL, lw=1.5, ls=:dash, label="IF-RK4")
    else
        println("  WARNING: IF-RK4 solution blew up!")
    end

    fig = plot(p1, p2, layout=(1, 2), size=(1200, 500))

    savefig(fig, joinpath(OUTPUT_DIR, "kdv_rk_comparison.pdf"))
    savefig(fig, joinpath(OUTPUT_DIR, "kdv_rk_comparison.png"))
    println("\nFigure saved to $(joinpath(OUTPUT_DIR, "kdv_rk_comparison.pdf"))")

    display(fig)
end

main()
