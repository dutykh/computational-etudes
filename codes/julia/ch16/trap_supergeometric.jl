#!/usr/bin/env julia
#=
trap_supergeometric.jl
Chapter 16: Integration of Periodic Functions

Computational Etude 16.5: Supergeometric Decay on e^cos(x).

f_5(x) = exp(cos x) is entire, so each new sample point gives more
than one extra digit.  The asymptotic envelope is

    |I - T_N| ~ 2*sqrt(2*pi/N) * (e/(2*N))^N.

Generates Figure 16.5: supergeometric.pdf

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

function main()
    f5(x) = exp(cos(x))
    I_exact = 2π * besseli(0, 1)
    println("Exact integral: 2π*I_0(1) = ", I_exact)

    N_values = 1:16
    errors = [abs(trapezoidal_periodic(f5, N) - I_exact) for N in N_values]
    errors = max.(errors, 1e-17)

    e_const = exp(1)
    theory = [2 * sqrt(2π / N) * (e_const / (2N))^N for N in N_values]
    theory = max.(theory, 1e-17)

    println("\n   N                I_N                |I_N - I|         theory")
    for (i, N) in enumerate(N_values)
        I_N = trapezoidal_periodic(f5, N)
        println(N, "  ", I_N, "  ", errors[i], "  ", theory[i])
    end

    p = plot(collect(N_values), errors,
             yscale = :log10,
             xlabel = L"Number of nodes $N$",
             ylabel = L"Absolute error $|I_N - I|$",
             title = L"Supergeometric decay for $f(x) = e^{\cos x}$",
             label = "Trapezoidal error",
             marker = :circle, markersize = 6, color = CORAL,
             linewidth = 1.2, legend = :topright,
             grid = true, gridalpha = 0.3, framestyle = :box,
             size = (700, 500),
             xlims = (0, 17), ylims = (1e-18, 1e2))
    plot!(p, collect(N_values), theory,
          label = L"$2\sqrt{2\pi/N}\,(e/2N)^N$",
          linestyle = :dash, color = NAVY, linewidth = 1.0)

    savefig(p, joinpath(OUTPUT_DIR, "supergeometric.pdf"))
    savefig(p, joinpath(OUTPUT_DIR, "supergeometric.png"))
    println("\nFigure saved to ", joinpath(OUTPUT_DIR, "supergeometric.pdf"))
end

main()
