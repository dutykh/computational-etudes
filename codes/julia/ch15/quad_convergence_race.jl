#!/usr/bin/env julia
#=
quad_convergence_race.jl
Chapter 15: Quadrature in Spectral Methods

Computational Etude 15.5: The Convergence Race.

Reproduces the essential content of Figure 2 from Trefethen (2008),
comparing the convergence of Gauss--Legendre and Clenshaw--Curtis
quadrature for six functions of varying regularity on [-1, 1]:

    (a) x^{20}        -- polynomial, finite exactness threshold
    (b) e^x            -- entire function, supergeometric convergence
    (c) e^{-x^2}       -- entire, supergeometric convergence
    (d) 1/(1+16x^2)    -- analytic with poles at x = +- i/4, geometric rate
    (e) e^{-1/x^2}     -- C^infinity but not analytic at x = 0
    (f) |x|^3          -- C^2 but not C^3, algebraic convergence

The six-panel (3 x 2) figure shows the absolute quadrature error
|I - I_n| vs the number of nodes n on a semilogarithmic scale for
n = 1, 2, ..., 40.

For smooth analytic functions, both rules converge at similar rates
(geometric or supergeometric). For functions with limited regularity,
both degrade to algebraic convergence. The key insight of Trefethen
(2008) is that Clenshaw--Curtis is competitive with Gauss for most
practical integrands.

Reference: L. N. Trefethen, "Is Gauss Quadrature Better than
           Clenshaw--Curtis?", SIAM Review 50(1):67--87, 2008, Figure 2.

Generates Figure 15.5: convergence_race.pdf

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
=#

using LinearAlgebra
using FFTW
using SpecialFunctions
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
    gauss_legendre(n)

Compute n-point Gauss--Legendre nodes and weights via the Golub--Welsch algorithm.
"""
function gauss_legendre(n)
    if n == 1
        return [0.0], [2.0]
    end
    beta = [j / sqrt(4j^2 - 1) for j in 1:n-1]
    J = diagm(1 => beta, -1 => beta)
    F = eigen(Symmetric(J))
    idx = sortperm(F.values)
    x = F.values[idx]
    w = 2.0 .* F.vectors[1, idx] .^ 2
    return x, w
end

"""
    clenshaw_curtis_weights(n)

Compute (n+1)-point Clenshaw--Curtis nodes and weights via DCT-I / FFT.
"""
function clenshaw_curtis_weights(n)
    if n == 0
        return [0.0], [2.0]
    end
    if n == 1
        return [1.0, -1.0], [1.0, 1.0]
    end
    x = [cos(j * π / n) for j in 0:n]
    c = zeros(n + 1)
    for k in 0:2:n
        c[k+1] = 2.0 / (1.0 - k^2)
    end
    v = vcat(c, c[n:(-1):2])  # length 2n, mirror
    f = real.(fft(v))
    w = f[1:n+1] / n
    w[1] /= 2
    w[end] /= 2
    return x, w
end

# -------------------------------------------------------------------------
# Test functions and their exact integrals
# -------------------------------------------------------------------------
f_poly(x)  = x^20
f_exp(x)   = exp(x)
f_gauss(x) = exp(-x^2)
f_runge(x) = 1.0 / (1.0 + 16.0 * x^2)

function f_flat(x)
    # f(x) = e^{-1/x^2}, C^infinity but not analytic at x = 0
    x == 0.0 ? 0.0 : exp(-1.0 / x^2)
end

f_abs3(x) = abs(x)^3

# Exact integrals over [-1, 1]
const I_poly  = 2.0 / 21.0                          # integral x^20 dx
const I_exp   = exp(1) - exp(-1)                     # e - e^{-1}
const I_gauss_val = sqrt(π) * erf(1.0)              # sqrt(pi) * erf(1)
const I_runge = atan(4.0) / 2.0                      # arctan(4)/2
const I_flat, _ = quadgk(x -> x == 0.0 ? 0.0 : exp(-1.0 / x^2), -1, 1,
                          rtol=1e-15, atol=1e-15)
const I_abs3  = 0.5                                  # 2 * integral_0^1 x^3 dx = 1/2

# Collect into arrays for iteration
const functions = [f_poly, f_exp, f_gauss, f_runge, f_flat, f_abs3]
const integrals = [I_poly, I_exp, I_gauss_val, I_runge, I_flat, I_abs3]
const titles = [
    L"x^{20}",
    L"e^x",
    L"e^{-x^2}",
    L"1/(1 + 16x^2)",
    L"e^{-1/x^2}",
    L"|x|^3",
]

function main()
    n_max = 40
    n_values = 1:n_max
    floor_val = 1e-16

    # Preallocate error arrays: [6 functions] x [n_max values]
    err_gauss = zeros(6, n_max)
    err_cc    = zeros(6, n_max)

    for idx in 1:6
        fn = functions[idx]
        I_exact = integrals[idx]

        for (i, n) in enumerate(n_values)
            # Gauss--Legendre with n+1 points
            x_gl, w_gl = gauss_legendre(n + 1)
            I_gl = dot(w_gl, fn.(x_gl))
            err_gauss[idx, i] = abs(I_gl - I_exact)

            # Clenshaw--Curtis with n+1 points
            x_cc, w_cc = clenshaw_curtis_weights(n)
            I_cc = dot(w_cc, fn.(x_cc))
            err_cc[idx, i] = abs(I_cc - I_exact)
        end
    end

    # Clamp zeros for log plot
    err_gauss = max.(err_gauss, floor_val)
    err_cc    = max.(err_cc, floor_val)

    # Create 3x2 figure
    panels = []
    nv = collect(n_values)

    for idx in 1:6
        show_legend = (idx == 1)
        show_xlabel = (idx >= 5)
        show_ylabel = (idx % 2 == 1)

        pi = plot(yscale=:log10,
                  title=titles[idx], titlefontsize=11,
                  xlims=(0, n_max + 1), ylims=(floor_val * 0.5, 1e1),
                  grid=true, gridalpha=0.3, gridlinewidth=0.5,
                  framestyle=:box, tickfontsize=9,
                  guidefontsize=10)

        if show_xlabel
            xlabel!(pi, L"n")
        end
        if show_ylabel
            ylabel!(pi, L"|I - I_n|")
        end

        # Gauss--Legendre: filled circles
        scatter!(pi, nv, err_gauss[idx, :], markershape=:circle,
                 markersize=3, color=NAVY, markerstrokecolor=NAVY,
                 markerstrokewidth=0.5,
                 label=show_legend ? "Gauss-Legendre" : "")

        # Clenshaw--Curtis: open circles
        scatter!(pi, nv, err_cc[idx, :], markershape=:circle,
                 markersize=3, color=TEAL,
                 markerstrokecolor=TEAL, markerstrokewidth=0.8,
                 markerfillcolor=:white,
                 label=show_legend ? "Clenshaw-Curtis" : "")

        if show_legend
            plot!(pi, legend=:topright, legendfontsize=8)
        else
            plot!(pi, legend=false)
        end

        push!(panels, pi)
    end

    p = plot(panels..., layout=(3, 2), size=(1000, 1000))

    # Save output
    savefig(p, joinpath(OUTPUT_DIR, "convergence_race.pdf"))
    savefig(p, joinpath(OUTPUT_DIR, "convergence_race.png"))
    println("Figure saved to $(joinpath(OUTPUT_DIR, "convergence_race.pdf"))")

    return p
end

main()
