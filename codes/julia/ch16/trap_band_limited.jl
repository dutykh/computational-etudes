#!/usr/bin/env julia
#=
trap_band_limited.jl
Chapter 16: Integration of Periodic Functions

Computational Etude 16.2: Band-Limited Exactness and One-Mode Aliasing.

(a) The N-point periodic trapezoidal rule is exact for any
    trigonometric polynomial of degree m < N.
(b) For a single mode cos(k*theta) with N fixed, the rule still
    gives the correct value of zero unless k is an integer multiple
    of N -- even when k > N/2 (the aliased Nyquist regime).

Generates Figure 16.2: band_limited.pdf  (two-panel figure)

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
=#

using Plots
using LaTeXStrings
using Random

gr()

const NAVY  = RGB(20/255, 45/255, 110/255)
const CORAL = RGB(231/255, 76/255, 60/255)

const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch16", "julia")
mkpath(OUTPUT_DIR)

trapezoidal_periodic(f, N) = (2π/N) * sum(f.(2π * (0:N-1) / N))

function main()
    # ----- Part (a): random trigonometric polynomial of degree m = 10 -----
    Random.seed!(42)
    m = 10
    a = randn(m + 1)
    b = randn(m + 1)
    function trig_poly(θ)
        s = a[1]
        for k in 1:m
            s += a[k+1] * cos(k * θ) + b[k+1] * sin(k * θ)
        end
        return s
    end
    I_exact = 2π * a[1]

    N_values_a = 4:30
    errors_a = [abs(trapezoidal_periodic(trig_poly, N) - I_exact) for N in N_values_a]
    errors_a = max.(errors_a, 1e-17)

    println("Part (a): random trig polynomial of degree m = $m")
    println("Exact integral 2π*a[1] = ", I_exact)
    println("Error at N = $(m): ", errors_a[m - 3])
    println("Error at N = $(m + 1): ", errors_a[m - 2])

    # ----- Part (b): single mode cos(k*theta), N = 16 fixed, sweep k -----
    N_fixed = 16
    k_values = 0:32
    errors_b = zeros(length(k_values))
    for (i, k) in enumerate(k_values)
        θ = 2π * (0:N_fixed-1) / N_fixed
        I_exact_k = (k == 0) ? 2π : 0.0
        errors_b[i] = abs((2π/N_fixed) * sum(cos.(k .* θ)) - I_exact_k)
    end
    errors_b = max.(errors_b, 1e-17)

    # ----- Two-panel plot -----
    p1 = plot(collect(N_values_a), errors_a,
              yscale = :log10,
              xlabel = L"Number of nodes $N$",
              ylabel = L"Absolute error $|I_N - I|$",
              title = L"(a) trig polynomial of degree $m = 10$",
              label = "", marker = :circle, markersize = 5, color = CORAL,
              linewidth = 1.1, grid = true, gridalpha = 0.3,
              framestyle = :box, xlims = (3, 31), ylims = (1e-18, 1e2))
    vline!(p1, [m], linestyle = :dash, color = RGB(0.5, 0.5, 0.5),
           linewidth = 0.9, label = "")

    bar_colors = [errors_b[i] > 1e-10 ? CORAL : NAVY for i in 1:length(k_values)]
    p2 = bar(collect(k_values), errors_b,
             yscale = :log10,
             xlabel = L"Frequency $k$",
             ylabel = L"$|I_N(\cos k\theta) - I(\cos k\theta)|$",
             title = L"(b) single modes, fixed $N = 16$",
             label = "", color = bar_colors, linecolor = bar_colors,
             grid = true, gridalpha = 0.3, framestyle = :box,
             xlims = (-0.5, 32.5), ylims = (1e-17, 1e1))
    vline!(p2, [N_fixed, 2*N_fixed], linestyle = :dash,
           color = RGB(0.5, 0.5, 0.5), linewidth = 0.9, label = "")

    p = plot(p1, p2, layout = (1, 2), size = (1100, 450))
    savefig(p, joinpath(OUTPUT_DIR, "band_limited.pdf"))
    savefig(p, joinpath(OUTPUT_DIR, "band_limited.png"))
    println("\nFigure saved to ", joinpath(OUTPUT_DIR, "band_limited.pdf"))
end

main()
