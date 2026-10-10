#!/usr/bin/env julia
#=
trap_poisson_kernel.jl
Chapter 16: Integration of Periodic Functions

Computational Etude 16.4: Geometric Convergence on the Poisson Kernel.

f_4(x) = 1 / (a - cos x), a > 1.  Predicted geometric rate r^N with
r = a - sqrt(a^2 - 1), and an explicit closed-form trapezoidal error

  I(f_4) - T_N(f_4) = -8*pi * r/(1 - r^2) * r^N / (1 - r^N)

derived in Weideman (2002), eq. (21).

Generates Figure 16.4: poisson_kernel.pdf  (two-panel figure)

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
    a = 2.0
    f4(x) = 1 / (a - cos(x))
    I_exact = 2π / sqrt(a^2 - 1)
    r = a - sqrt(a^2 - 1)
    α = log(a + sqrt(a^2 - 1))

    println("Poisson kernel f_4 = 1/(a - cos x), a = $a")
    println("Exact integral: ", I_exact)
    println("Geometric rate r = ", r)

    N_values = 2:2:50
    errors = [abs(trapezoidal_periodic(f4, N) - I_exact) for N in N_values]
    errors = max.(errors, 1e-17)

    theory = [8π * r/(1 - r^2) * r^N / (1 - r^N) for N in N_values]
    theory = max.(theory, 1e-17)

    # ----- Two-panel plot -----
    p1 = plot(collect(N_values), errors,
              yscale = :log10,
              xlabel = L"Number of nodes $N$",
              ylabel = L"Absolute error $|I_N - I|$",
              title = L"Geometric decay at rate $r^N$ for $a = 2$",
              label = "Trapezoidal error",
              marker = :circle, markersize = 5, color = CORAL,
              linewidth = 1.2, legend = :topright,
              grid = true, gridalpha = 0.3, framestyle = :box,
              ylims = (1e-18, 1e1))
    plot!(p1, collect(N_values), theory,
          label = "closed-form (Weideman 2002)",
          linestyle = :dash, color = NAVY, linewidth = 1.0)

    # Pole picture
    p2 = plot([0, 2π], [0, 0], color = :black, linewidth = 1.0,
              label = "", framestyle = :box,
              xlims = (-0.5, 2π + 0.5), ylims = (-2, 2),
              xlabel = L"$\mathrm{Re}\,\theta$",
              ylabel = L"$\mathrm{Im}\,\theta$",
              title = L"Poles in the complex $\theta$-plane",
              xticks = ([0, π/2, π, 3π/2, 2π],
                        ["0", "π/2", "π", "3π/2", "2π"]),
              grid = true, gridalpha = 0.3)
    # Strip of analyticity (shaded)
    plot!(p2, [-0.5, 2π+0.5, 2π+0.5, -0.5], [-α, -α, α, α],
          seriestype = :shape, fillcolor = NAVY, fillalpha = 0.15,
          linecolor = :transparent, label = "strip of analyticity")
    hline!(p2, [α, -α], linestyle = :dash, color = NAVY, linewidth = 0.8,
           label = "")
    # Trapezoidal nodes
    N_demo = 12
    θ_nodes = 2π * (0:N_demo-1) / N_demo
    scatter!(p2, θ_nodes, zeros(length(θ_nodes)),
             color = CORAL, markersize = 6, markerstrokecolor = CORAL,
             label = "")
    # Poles
    scatter!(p2, [0, 2π, 0, 2π], [α, α, -α, -α],
             marker = :x, color = NAVY, markersize = 10,
             markerstrokewidth = 2, label = "")

    p = plot(p1, p2, layout = (1, 2), size = (1100, 450))
    savefig(p, joinpath(OUTPUT_DIR, "poisson_kernel.pdf"))
    savefig(p, joinpath(OUTPUT_DIR, "poisson_kernel.png"))
    println("\nFigure saved to ", joinpath(OUTPUT_DIR, "poisson_kernel.pdf"))
end

main()
