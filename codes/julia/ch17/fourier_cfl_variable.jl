#!/usr/bin/env julia
#=
fourier_cfl_variable.jl
Chapter 17: Stability of Spectral Methods

Computational Etude 17.3 (variable-coefficient extension):
Frozen-coefficient CFL prediction for variable-speed advection.

Two-panel figure:
  (a) Speed profile c(x) = 1 + 0.5*sin(x) on [0, 2*pi) with c_max
      and c_mean annotated.
  (b) Log-log plot of critical dt vs N for RK4 applied to Fourier
      advection u_t + c(x)*u_x = 0.  Experimental binary search
      compared with frozen-coefficient prediction 2.83/(c_max*N/2).

Generates Figure: fourier_cfl_variable.pdf

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

const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch17", "julia")
mkpath(OUTPUT_DIR)

"""RHS of u_t + c(x)*u_x = 0: returns -c(x) * F^{-1}(ik * F(u))."""
function variable_rhs(u, c, k)
    u_hat = fft(u)
    du = real.(ifft(1im .* k .* u_hat))
    return -c .* du
end

"""Run n_steps of RK4; return true if solution stays bounded."""
function is_stable(dt, u0, c, k; n_steps=500)
    u = copy(u0)
    for _ in 1:n_steps
        k1 = variable_rhs(u, c, k)
        k2 = variable_rhs(u .+ 0.5 .* dt .* k1, c, k)
        k3 = variable_rhs(u .+ 0.5 .* dt .* k2, c, k)
        k4 = variable_rhs(u .+ dt .* k3, c, k)
        u .+= dt / 6.0 .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
        if maximum(abs.(u)) > 1e6
            return false
        end
    end
    return true
end

"""Binary search for the critical dt of variable-coefficient advection."""
function find_dt_crit(N, c_func; tol=1e-6)
    x = 2π .* (0:N-1) ./ N
    c = c_func.(x)
    c_max = maximum(abs.(c))
    u0 = exp.(cos.(x))
    k = fftfreq(N, N / 1.0)

    dt_lo = 0.1 / (c_max * N)
    dt_hi = 10.0 / N

    # Ensure dt_hi is unstable
    while is_stable(dt_hi, u0, c, k)
        dt_hi *= 2.0
        dt_hi > 100 && return dt_hi
    end

    # Binary search
    for _ in 1:60
        dt_mid = 0.5 * (dt_lo + dt_hi)
        if is_stable(dt_mid, u0, c, k)
            dt_lo = dt_mid
        else
            dt_hi = dt_mid
        end
        (dt_hi - dt_lo < tol * dt_lo) && break
    end

    return 0.5 * (dt_lo + dt_hi)
end

function main()
    c_func(x) = 1.0 + 0.5 * sin(x)
    c_max = 1.5
    c_mean = 1.0

    # --- Panel (a): speed profile ---
    x_plot = range(0, 2π, length=500)
    c_vals = c_func.(x_plot)

    p1 = plot(x_plot, c_vals, color=NAVY, lw=2.0,
              label=L"c(x) = 1 + 0.5\sin(x)",
              xlabel=L"x", ylabel=L"c(x)",
              title="(a) Variable speed profile",
              xlims=(0, 2π), ylims=(0.3, 1.7),
              grid=true, gridalpha=0.3)
    hline!(p1, [c_max], ls=:dash, lw=1.2, color=CORAL,
           label=L"c_{\max} = 1.5")
    hline!(p1, [c_mean], ls=:dot, lw=1.2, color=SKY,
           label=L"\langle c \rangle = 1.0")

    # --- Panel (b): dt_crit vs N ---
    N_values = [16, 32, 64, 128, 256]
    dt_crits = Float64[]

    for Nv in N_values
        dtc = find_dt_crit(Nv, c_func)
        push!(dt_crits, dtc)
        @printf("N = %4d, dt_crit = %.6f, dt_crit * c_max * (N/2) = %.4f\n",
                Nv, dtc, dtc * c_max * Nv / 2)
    end

    N_arr = Float64.(N_values)

    p2 = plot(N_arr, dt_crits, xscale=:log10, yscale=:log10,
              marker=:circle, ms=6, lw=1.5, color=NAVY,
              label=L"\Delta t_{\mathrm{crit}}\ \mathrm{(experiment)}",
              xlabel=L"$N$", ylabel=L"\Delta t_{\mathrm{crit}}",
              title="(b) Frozen-coefficient prediction vs experiment",
              grid=true, gridalpha=0.3, minorgrid=true)

    # Frozen-coefficient prediction
    N_ref = range(12, 300, length=100)
    plot!(p2, N_ref, 2.83 ./ (c_max .* N_ref ./ 2), ls=:dash, lw=1.2,
          color=CORAL, label=L"2.83\, /\, (c_{\max} \cdot N/2)")

    # Constant-coefficient reference
    plot!(p2, N_ref, 2.83 ./ (N_ref ./ 2), ls=:dot, lw=1.2,
          color=SKY, label=L"2.83\, /\, (N/2)\ \mathrm{[const.\ coeff.]}")

    # Slope annotation
    annotate!(p2, 80, 0.04, text("slope = -1", 9, CORAL, rotation=-28))

    fig = plot(p1, p2, layout=(1, 2), size=(1100, 450))

    savefig(fig, joinpath(OUTPUT_DIR, "fourier_cfl_variable.pdf"))
    savefig(fig, joinpath(OUTPUT_DIR, "fourier_cfl_variable.png"))
    println("\nFigure saved to $(joinpath(OUTPUT_DIR, "fourier_cfl_variable.pdf"))")

    display(fig)
end

main()
