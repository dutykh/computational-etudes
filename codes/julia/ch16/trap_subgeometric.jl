#!/usr/bin/env julia
#=
trap_subgeometric.jl
Chapter 16: Integration of Periodic Functions

Computational Etude 16.6: Subgeometric Decay on Weideman's f_6.

f_6(x) = exp((cos x - 1)/(cos x + 1)) is C^∞-periodic but has an
essential singularity at x = π.  The trapezoidal-rule error decays at
the subgeometric rate exp(-(3/2) * N^{2/3}).

Diagnostic visual: replot vs N^{2/3}, the curve becomes a straight line.

Generates Figure 16.6: subgeometric.pdf  (two-panel figure)

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
=#

using Plots
using LaTeXStrings
using SpecialFunctions

gr()

const NAVY  = RGB(20/255, 45/255, 110/255)
const CORAL = RGB(231/255, 76/255, 60/255)

const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch16", "julia")
mkpath(OUTPUT_DIR)

trapezoidal_periodic(f, N) = (2π/N) * sum(f.(2π * (0:N-1) / N))

function f6(x)
    c = cos(x)
    (c + 1) > 1e-15 ? exp((c - 1) / (c + 1)) : 0.0
end

function main()
    I_exact = 2 * exp(1) * π * (1 - erf(1))
    println("Exact integral: 2*e*π*(1 - erf(1)) = ", I_exact)

    N_values = [4, 6, 8, 10, 12, 16, 20, 24, 30, 40, 50, 60, 80, 100, 120, 160, 200]
    errors = [abs(trapezoidal_periodic(f6, N) - I_exact) for N in N_values]
    errors = max.(errors, 1e-17)

    theory = [8 * sqrt(π/3) * exp(-1.5 * N^(2/3)) for N in N_values]

    p1 = plot(Float64.(N_values), errors,
              yscale = :log10,
              xlabel = L"Number of nodes $N$",
              ylabel = L"Absolute error $|I_N - I|$",
              title = L"(a) error vs $N$: looks slower than geometric",
              label = "Trapezoidal error",
              marker = :circle, markersize = 5, color = CORAL,
              linewidth = 1.2, legend = :topright,
              grid = true, gridalpha = 0.3, framestyle = :box,
              ylims = (1e-17, 1e0))
    plot!(p1, Float64.(N_values), theory,
          label = L"$8\sqrt{\pi/3}\,e^{-3 N^{2/3}/2}$",
          linestyle = :dash, color = NAVY, linewidth = 1.0)

    N23 = Float64.(N_values) .^ (2/3)
    p2 = plot(N23, errors,
              yscale = :log10,
              xlabel = L"$N^{2/3}$",
              ylabel = L"Absolute error $|I_N - I|$",
              title = L"(b) error vs $N^{2/3}$: linear on semilog",
              label = "Trapezoidal error",
              marker = :circle, markersize = 5, color = CORAL,
              linewidth = 1.2, legend = :topright,
              grid = true, gridalpha = 0.3, framestyle = :box,
              ylims = (1e-17, 1e0))
    plot!(p2, N23, theory,
          label = L"$8\sqrt{\pi/3}\,e^{-3 N^{2/3}/2}$",
          linestyle = :dash, color = NAVY, linewidth = 1.0)

    p = plot(p1, p2, layout = (1, 2), size = (1100, 450))
    savefig(p, joinpath(OUTPUT_DIR, "subgeometric.pdf"))
    savefig(p, joinpath(OUTPUT_DIR, "subgeometric.png"))
    println("\nFigure saved to ", joinpath(OUTPUT_DIR, "subgeometric.pdf"))
end

main()
