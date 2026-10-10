#!/usr/bin/env julia
#=
quad_gauss_hermite_failure.jl
Chapter 15: Quadrature in Spectral Methods

Computational Etude 15.9: Gauss--Hermite vs Truncated Trapezoidal Rule.

Compares two quadrature strategies for the oscillatory integral

    I = int_{-inf}^{inf} e^{-x^2} cos(x^3) dx

which is a prototypical example of a Gaussian-weighted integral with a
non-polynomial integrand.

Gauss--Hermite quadrature:  The weight function e^{-x^2} is built into
the rule, so I_n = sum_k w_k cos(x_k^3). Theory predicts root-exponential
convergence: errors decay like O(exp(-C sqrt(n))), which appears as a
parabolic curve on a semilog plot of error vs n.

Truncated trapezoidal rule: The full integrand g(x) = e^{-x^2} cos(x^3)
is integrated on a finite interval [-L, L] using the composite trapezoidal
rule. Because g(x) decays super-algebraically as |x| -> inf, the
truncation error is negligible for moderate L. The trapezoidal rule then
converges exponentially in n (straight line on semilog), dramatically
outperforming Gauss--Hermite.

This etude demonstrates that Gauss rules matched to the weight function
are not always the best choice; a simpler rule on a truncated domain
can be far more effective when the integrand decays rapidly.

Generates Figure 15.9: gauss_hermite_failure.pdf

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
=#

using LinearAlgebra
using QuadGK
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
    # -------------------------------------------------------------------------
    # Reference value via high-accuracy adaptive quadrature
    # -------------------------------------------------------------------------
    integrand(x) = exp(-x^2) * cos(x^3)
    I_exact, _ = quadgk(integrand, -Inf, Inf, rtol=1e-15, atol=1e-15)
    println("Reference integral value: $I_exact")

    # -------------------------------------------------------------------------
    # Gauss--Hermite quadrature: I_n = sum w_k cos(x_k^3)
    # The weight function e^{-x^2} is absorbed into the weights.
    # -------------------------------------------------------------------------
    n_gh = 10:10:1000
    err_gh = zeros(length(n_gh))

    for (i, n) in enumerate(n_gh)
        x_h, w_h = gauss_hermite(n)
        I_gh = dot(w_h, cos.(x_h .^ 3))
        err_gh[i] = abs(I_gh - I_exact)
    end

    # -------------------------------------------------------------------------
    # Truncated trapezoidal rule on [-L, L]
    # g(x) = e^{-x^2} cos(x^3), with L = 6 (truncation error < 10^{-15})
    # -------------------------------------------------------------------------
    L = 6.0
    n_trap = 10:2:200
    err_trap = zeros(length(n_trap))

    for (i, n) in enumerate(n_trap)
        h = 2.0 * L / n
        x = range(-L, L, length=n + 1) |> collect
        g = exp.(-x .^ 2) .* cos.(x .^ 3)
        # Composite trapezoidal rule
        I_trap = h * (0.5 * g[1] + sum(g[2:end-1]) + 0.5 * g[end])
        err_trap[i] = abs(I_trap - I_exact)
    end

    # -------------------------------------------------------------------------
    # Create figure
    # -------------------------------------------------------------------------
    floor_val = 1e-16

    p = plot(size=(700, 500), yscale=:log10,
             xlabel=L"Number of quadrature points $n$",
             ylabel=L"Absolute error $|I - I_n|$",
             title=L"Quadrature for $\int_{-\infty}^{\infty} e^{-x^2} \cos(x^3)\, dx$",
             titlefontsize=12, guidefontsize=11, tickfontsize=10,
             legendfontsize=10, legend=:topright,
             grid=true, gridalpha=0.3, gridlinewidth=0.5,
             framestyle=:box,
             xlims=(0, 1050), ylims=(1e-16, 1e1))

    # Gauss--Hermite errors
    scatter!(p, collect(n_gh), max.(err_gh, floor_val), markershape=:circle,
             markersize=3, color=NAVY, markerstrokecolor=NAVY,
             markerstrokewidth=0.5, label="Gauss-Hermite")

    # Trapezoidal errors
    scatter!(p, collect(n_trap), max.(err_trap, floor_val), markershape=:rect,
             markersize=3.5, color=:white, markerstrokecolor=CORAL,
             markerstrokewidth=0.7, label=L"Trapezoidal on $[-6, 6]$")

    # Machine epsilon reference
    hline!(p, [eps(Float64)], color=:gray, linestyle=:dot, linewidth=0.8,
           alpha=0.7, label="")
    annotate!(p, 950, 4e-16, text(L"\epsilon_{\mathrm{mach}}", 9, :gray, :right))

    # Save output
    savefig(p, joinpath(OUTPUT_DIR, "gauss_hermite_failure.pdf"))
    savefig(p, joinpath(OUTPUT_DIR, "gauss_hermite_failure.png"))
    println("Figure saved to $(joinpath(OUTPUT_DIR, "gauss_hermite_failure.pdf"))")

    return p
end

main()
