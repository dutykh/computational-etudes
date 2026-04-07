#!/usr/bin/env julia
#=
trap_poisson_ellipse.jl
Chapter 16: Integration of Periodic Functions

Computational Etude 16.1: Poisson's Ellipse, the Original Paradox.

Reproduces Figure 1.1 of Trefethen and Weideman, SIAM Rev. 2014: the
trapezoidal-rule perimeter of an ellipse with semi-axes 1/(2*pi) and
0.6/(2*pi).  The integral is (2/pi) * E(0.36); the trapezoidal error
decays as 3^(-N) because the integrand has branch points at
theta = +/- i*log(3) in the complex plane.

Generates Figure 16.1: poisson_ellipse.pdf

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

# Book color scheme
const NAVY  = RGB(20/255, 45/255, 110/255)
const CORAL = RGB(231/255, 76/255, 60/255)

const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch16", "julia")
mkpath(OUTPUT_DIR)

trapezoidal_periodic(f, N) = (2π/N) * sum(f.(2π * (0:N-1) / N))

function main()
    integrand(t) = sqrt(1 - 0.36 * sin(t)^2)
    I_exact = (2 / π) * ellipe(0.36)
    println("Exact value: I = ", I_exact)

    N_values = collect(4:4:200)
    errors = [abs(trapezoidal_periodic(integrand, N) / (2π) - I_exact) for N in N_values]
    errors = max.(errors, 1e-17)

    theory = max.(3.0 .^ (-Float64.(N_values)), 1e-17)

    println("\nN/4    I_N                       |I_N - I|")
    for i in 1:10
        N = N_values[i]
        I_N = trapezoidal_periodic(integrand, N) / (2π)
        println(N ÷ 4, "   ", I_N, "   ", errors[i])
    end

    p = plot(N_values .÷ 4, errors,
             yscale = :log10,
             xlabel = L"$N/4$ (number of independent samples)",
             ylabel = L"Absolute error $|I_N - I|$",
             title = "Poisson's ellipse: trapezoidal convergence is geometric",
             label = "Trapezoidal rule",
             marker = :circle, markersize = 4, color = CORAL,
             linewidth = 1.1, legend = :topright,
             grid = true, gridalpha = 0.3,
             framestyle = :box, size = (700, 480),
             xlims = (0, 50), ylims = (1e-18, 1e0))
    plot!(p, N_values .÷ 4, theory,
          label = L"$3^{-N}$", linestyle = :dash, color = NAVY, linewidth = 1.0)

    savefig(p, joinpath(OUTPUT_DIR, "poisson_ellipse.pdf"))
    savefig(p, joinpath(OUTPUT_DIR, "poisson_ellipse.png"))
    println("\nFigure saved to ", joinpath(OUTPUT_DIR, "poisson_ellipse.pdf"))
end

main()
