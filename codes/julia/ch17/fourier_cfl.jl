#!/usr/bin/env julia
#=
fourier_cfl.jl
Chapter 17: Stability of Spectral Methods

Computational Etude 17.3: Fourier CFL Scaling.

Two-panel figure:
  (a) Eigenvalues of the Fourier differentiation matrix for N=32,
      plotted on the imaginary axis from -16i to 16i.
  (b) Log-log plot of critical dt vs N for RK4 applied to Fourier
      advection u_t + u_x = 0.  Binary search for dt_crit.

Generates Figure 17.3: fourier_cfl_scaling.pdf

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

"""RK4 amplification factor g(z) = 1 + z + z^2/2 + z^3/6 + z^4/24."""
rk4_amplification(z) = 1.0 .+ z .+ z.^2 ./2 .+ z.^3 ./6 .+ z.^4 ./24

"""Maximum |g(dt * lambda)| over all Fourier eigenvalues for advection."""
function rk4_max_amplification(dt, N)
    k = fftfreq(N, N / 1.0)   # wavenumbers: 0,1,...,N/2-1,-N/2+1,...,-1
    z = -1im .* k .* dt       # dt * lambda_k
    return maximum(abs.(rk4_amplification(z)))
end

"""Binary search for critical dt where max|g| = 1."""
function find_dt_crit(N; tol=1e-8)
    dt_lo = 1e-10
    dt_hi = 1.0
    # Find unstable upper bound
    while rk4_max_amplification(dt_hi, N) <= 1.0
        dt_hi *= 2
        dt_hi > 100 && return dt_hi
    end

    for _ in 1:80
        dt_mid = 0.5 * (dt_lo + dt_hi)
        if rk4_max_amplification(dt_mid, N) <= 1.0
            dt_lo = dt_mid
        else
            dt_hi = dt_mid
        end
        (dt_hi - dt_lo < tol * dt_lo) && break
    end
    return 0.5 * (dt_lo + dt_hi)
end

function main()
    # --- Panel (a): Fourier eigenvalues for N=32 ---
    N = 32
    k = fftfreq(N, N / 1.0)
    k_sorted = sort(k)
    eigs = 1im .* k_sorted  # eigenvalues of d/dx

    p1 = plot(xlabel=L"\mathrm{Re}(\lambda)", ylabel=L"\mathrm{Im}(\lambda)",
              title=L"(a) Fourier eigenvalues, $N = 32$",
              xlims=(-2, 2), ylims=(-20, 20),
              grid=true, gridalpha=0.3)
    hline!(p1, [0], color=:gray, lw=0.5, label="")
    vline!(p1, [0], color=:gray, lw=0.5, label="")
    scatter!(p1, real.(eigs), imag.(eigs), color=NAVY, ms=5,
             markerstrokecolor=:white, markerstrokewidth=0.5, label="")
    # Annotate k_max
    annotate!(p1, 0.7, 14.5, text(L"ik_{\max} = i \cdot 16", 9, CORAL))
    annotate!(p1, 0.7, -14.5, text(L"-ik_{\max}", 9, CORAL))

    # --- Panel (b): dt_crit vs N ---
    N_values = [16, 32, 64, 128, 256]
    dt_crits = Float64[]

    for Nv in N_values
        dtc = find_dt_crit(Nv)
        push!(dt_crits, dtc)
        @printf("N = %4d, dt_crit = %.6f, dt_crit * (N/2) = %.4f\n",
                Nv, dtc, dtc * Nv / 2)
    end

    N_arr = Float64.(N_values)

    p2 = plot(N_arr, dt_crits, xscale=:log10, yscale=:log10,
              marker=:circle, ms=6, lw=1.5, color=NAVY,
              label=L"\Delta t_{\mathrm{crit}}\ \mathrm{(binary\ search)}",
              xlabel=L"$N$", ylabel=L"\Delta t_{\mathrm{crit}}",
              title="(b) RK4 CFL scaling for Fourier advection",
              grid=true, gridalpha=0.3, minorgrid=true)

    # Reference line: 2.83 / (N/2)
    N_ref = range(12, 300, length=100)
    plot!(p2, N_ref, 2.83 ./ (N_ref ./ 2), ls=:dash, lw=1.2, color=CORAL,
          label=L"2.83 / (N/2)")

    # Slope annotation
    annotate!(p2, 80, 0.06, text("slope = -1", 9, CORAL, rotation=-28))

    fig = plot(p1, p2, layout=(1, 2), size=(1100, 450))

    savefig(fig, joinpath(OUTPUT_DIR, "fourier_cfl_scaling.pdf"))
    savefig(fig, joinpath(OUTPUT_DIR, "fourier_cfl_scaling.png"))
    println("\nFigure saved to $(joinpath(OUTPUT_DIR, "fourier_cfl_scaling.pdf"))")

    display(fig)
end

main()
