#!/usr/bin/env julia
#=
trap_algebraic_decay.jl
Chapter 16: Integration of Periodic Functions

Computational Etude 16.3: Algebraic Convergence on |sin(x/2)|^k.

f_2(x) = |sin(x/2)|     --> O(N^{-2})
f_3(x) = |sin(x/2)|^3   --> O(N^{-4})

Generates Figure 16.3: algebraic_decay.pdf

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

trapezoidal_periodic(f, N) = (2π/N) * sum(f.(2π * (0:N-1) / N))

function main()
    f2(x) = abs(sin(x / 2))
    f3(x) = abs(sin(x / 2))^3
    I2 = 4.0
    I3 = 8 / 3

    N_values = [8, 16, 32, 64, 128, 256, 512, 1024, 2048]
    err2 = [abs(trapezoidal_periodic(f2, N) - I2) for N in N_values]
    err3 = [abs(trapezoidal_periodic(f3, N) - I3) for N in N_values]

    theory_f2 = (π^2 / 3) ./ (Float64.(N_values) .^ 2)
    theory_f3 = (π^4 / 30) ./ (Float64.(N_values) .^ 4)

    println("f2 errors:  ", err2)
    println("f3 errors:  ", err3)

    p = plot(N_values, err2,
             xscale = :log10, yscale = :log10,
             xlabel = L"Number of nodes $N$",
             ylabel = L"Absolute error $|I_N - I|$",
             title = "Algebraic convergence: regularity of the periodic extension",
             label = L"$f_2 = |\sin(x/2)|$",
             marker = :circle, markersize = 5, color = CORAL,
             linewidth = 1.2, legend = :bottomleft,
             grid = true, gridalpha = 0.3,
             framestyle = :box, size = (700, 500),
             xlims = (6, 3000))
    plot!(p, N_values, theory_f2, label = L"$\pi^2/(3 N^2)$",
          linestyle = :dash, color = CORAL, linewidth = 0.8, alpha = 0.6)
    plot!(p, N_values, err3, label = L"$f_3 = |\sin(x/2)|^3$",
          marker = :square, markersize = 5, color = NAVY, linewidth = 1.2)
    plot!(p, N_values, theory_f3, label = L"$\pi^4/(30 N^4)$",
          linestyle = :dash, color = NAVY, linewidth = 0.8, alpha = 0.6)

    savefig(p, joinpath(OUTPUT_DIR, "algebraic_decay.pdf"))
    savefig(p, joinpath(OUTPUT_DIR, "algebraic_decay.png"))
    println("\nFigure saved to ", joinpath(OUTPUT_DIR, "algebraic_decay.pdf"))
end

main()
