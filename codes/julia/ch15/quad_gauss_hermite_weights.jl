#!/usr/bin/env julia
#=
quad_gauss_hermite_weights.jl
Chapter 15: Quadrature in Spectral Methods

Computational Etude 15.8: Gauss--Hermite Weight Distribution.

Examines the pathological behaviour of Gauss--Hermite quadrature weights
for large n, which limits the practical utility of this classical rule
for high-accuracy integration.

Left panel: For n = 100, the Gauss--Hermite weights w_k are plotted
on a semilogarithmic scale against the corresponding nodes x_k. The
weights span a vast dynamic range: the central weights are of O(1),
while the outermost weights fall far below machine epsilon (approx.
10^{-16}). This means that when the weights are represented in
double-precision floating-point arithmetic, a significant fraction
carry no useful information.

Right panel: The fraction of weights below machine epsilon is plotted
as a function of n for n = 10, 20, ..., 500. This fraction increases
monotonically and approaches 1, showing that Gauss--Hermite quadrature
becomes numerically useless in standard double-precision for large n.

This motivates the search for alternative quadrature strategies on
unbounded domains, such as the truncated Gauss--Legendre approach
discussed in subsequent etudes.

Generates Figure 15.8: gauss_hermite_weights.pdf

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
=#

using LinearAlgebra
using Plots
using LaTeXStrings

gr()

# Book color scheme
const NAVY   = RGB(20/255, 45/255, 110/255)
const SKY    = RGB(120/255, 150/255, 210/255)
const CORAL  = RGB(231/255, 76/255, 60/255)
const TEAL   = RGB(22/255, 160/255, 133/255)
const PURPLE = RGB(142/255, 68/255, 173/255)
const ORANGE = RGB(230/255, 126/255, 34/255)

# Output path
const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch15", "julia")
mkpath(OUTPUT_DIR)

"""
    gauss_hermite(n)

Compute n-point Gauss--Hermite nodes and weights via the Golub--Welsch algorithm.
The Jacobi matrix for Hermite polynomials (physicists' convention) has zero
diagonal and off-diagonal entries beta_j = sqrt(j/2) for j = 1, ..., n-1.
"""
function gauss_hermite(n)
    if n == 1
        return [0.0], [sqrt(π)]
    end
    beta = [sqrt(j / 2) for j in 1:n-1]
    J = diagm(1 => beta, -1 => beta)
    F = eigen(Symmetric(J))
    idx = sortperm(F.values)
    x = F.values[idx]
    w = sqrt(π) .* F.vectors[1, idx] .^ 2
    return x, w
end

function main()
    eps_mach = eps(Float64)  # Machine epsilon (~2.22e-16)

    # -------------------------------------------------------------------------
    # Left panel: weights for n = 100
    # -------------------------------------------------------------------------
    n_demo = 100
    x_h, w_h = gauss_hermite(n_demo)

    # -------------------------------------------------------------------------
    # Right panel: fraction of weights below eps vs n
    # -------------------------------------------------------------------------
    n_values = 10:10:500
    fraction_below = zeros(length(n_values))

    for (i, n) in enumerate(n_values)
        _, w = gauss_hermite(n)
        fraction_below[i] = sum(abs.(w) .< eps_mach) / n
    end

    # -------------------------------------------------------------------------
    # Create figure
    # -------------------------------------------------------------------------

    # --- Left panel: weight distribution ---
    p1 = plot(yscale=:log10,
              xlabel=L"Node $x_k$",
              ylabel=L"Weight $|w_k|$",
              title=L"(a) Gauss-Hermite weights ($n = 100$)",
              titlefontsize=11, guidefontsize=11, tickfontsize=10,
              legendfontsize=9, legend=:bottom,
              grid=true, gridalpha=0.3, gridlinewidth=0.5,
              framestyle=:box,
              ylims=(1e-300, 1e2))

    scatter!(p1, x_h, abs.(w_h), markershape=:circle,
             markersize=2.5, color=NAVY, markerstrokecolor=NAVY,
             markerstrokewidth=0.3, label=L"|w_k|")

    # Machine epsilon reference line
    hline!(p1, [eps_mach], color=CORAL, linestyle=:dash, linewidth=1.0,
           label=L"Machine $\epsilon \approx 2.2 \times 10^{-16}$")

    # --- Right panel: fraction below epsilon ---
    p2 = plot(xlabel=L"Number of points $n$",
              ylabel=L"Fraction of $|w_k| < \epsilon$",
              title="(b) Fraction of underflowed weights",
              titlefontsize=11, guidefontsize=11, tickfontsize=10,
              legend=false,
              grid=true, gridalpha=0.3, gridlinewidth=0.5,
              framestyle=:box,
              xlims=(0, 510), ylims=(-0.02, 1.05))

    plot!(p2, collect(n_values), fraction_below, markershape=:circle,
          markersize=3, color=NAVY, markerstrokecolor=NAVY,
          markerstrokewidth=0.5, linewidth=1.2)

    hline!(p2, [1.0], color=:gray, linestyle=:dot, linewidth=0.8)

    # Combine into two-panel figure
    p = plot(p1, p2, layout=(1, 2), size=(1100, 450))

    # Save output
    savefig(p, joinpath(OUTPUT_DIR, "gauss_hermite_weights.pdf"))
    savefig(p, joinpath(OUTPUT_DIR, "gauss_hermite_weights.png"))
    println("Figure saved to $(joinpath(OUTPUT_DIR, "gauss_hermite_weights.pdf"))")

    return p
end

main()
