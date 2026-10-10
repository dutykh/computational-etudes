#!/usr/bin/env julia
#=
trap_real_line.jl
Chapter 16: Integration of Periodic Functions

Computational Etude 16.7: Trapezoidal Rule on the Real Line.

Compute (1/sqrt(pi)) * int_{-infty}^{infty} exp(-x^2) dx = 1 by the
truncated trapezoidal rule with step h = 2*pi/N.  Theorem 5.1 of
Trefethen-Weideman (2014) predicts O(exp(-pi^2/h)) decay.

Generates Figure 16.7: real_line_gaussian.pdf

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
=#

using Plots
using LaTeXStrings

gr()

const NAVY  = RGB(20/255, 45/255, 110/255)
const CORAL = RGB(231/255, 76/255, 60/255)

const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch16", "julia")
mkpath(OUTPUT_DIR)

function trapezoidal_real_line(w, h, n_max)
    k = -n_max:n_max
    return h * sum(w.(k * h))
end

function main()
    w(x) = exp(-x^2) / sqrt(π)
    I_exact = 1.0

    N_values = 1:12
    h_values = [2π / N for N in N_values]
    errors = zeros(length(N_values))
    for (i, h) in enumerate(h_values)
        n_max = max(Int(ceil(28 / h)), 30)
        I_h = trapezoidal_real_line(w, h, n_max)
        errors[i] = abs(I_h - I_exact)
    end
    errors = max.(errors, 1e-17)

    theory = [exp(-π^2 / h) for h in h_values]

    println("   N        h                I_h                |I_h - I|")
    for (i, N) in enumerate(N_values)
        h = h_values[i]
        n_max = max(Int(ceil(28 / h)), 30)
        I_h = trapezoidal_real_line(w, h, n_max)
        println(N, "  ", h, "  ", I_h, "  ", errors[i])
    end

    p = plot(collect(N_values), errors,
             yscale = :log10,
             xlabel = L"$N$ (with step $h = 2\pi/N$)",
             ylabel = "Absolute error",
             title = L"Real-line trapezoidal rule on $e^{-x^2}/\sqrt{\pi}$",
             label = L"$|I_h - I|$, $h = 2\pi/N$",
             marker = :circle, markersize = 6, color = CORAL,
             linewidth = 1.2, legend = :topright,
             grid = true, gridalpha = 0.3, framestyle = :box,
             size = (700, 500),
             xlims = (0, 13), ylims = (1e-18, 1e1))
    plot!(p, collect(N_values), theory,
          label = L"$e^{-\pi^2/h}$",
          linestyle = :dash, color = NAVY, linewidth = 1.0)

    savefig(p, joinpath(OUTPUT_DIR, "real_line_gaussian.pdf"))
    savefig(p, joinpath(OUTPUT_DIR, "real_line_gaussian.png"))
    println("\nFigure saved to ", joinpath(OUTPUT_DIR, "real_line_gaussian.pdf"))
end

main()
