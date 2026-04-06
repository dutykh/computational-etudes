#!/usr/bin/env julia
#=
quad_newton_cotes_runge.jl
Chapter 15: Quadrature in Spectral Methods

Computational Etude 15.2: The Runge Phenomenon in Quadrature.

Integrates the Runge function f(x) = 1/(1 + 25x^2) over [-1, 1] using
three quadrature rules of increasing order n:

    1. Newton--Cotes (equispaced nodes, Vandermonde-based weights).
    2. Gauss--Legendre (Golub--Welsch algorithm).
    3. Clenshaw--Curtis (Chebyshev nodes, FFT-based weights).

The exact integral is (2/5) arctan(5).

For Newton--Cotes, the weights are computed by solving the Vandermonde
system that enforces exactness on monomials 1, x, ..., x^n.  Because the
Vandermonde matrix becomes exponentially ill-conditioned for equispaced
nodes, the Newton--Cotes weights develop large alternating values and the
quadrature error grows without bound: this is the quadrature analogue of
the Runge phenomenon in polynomial interpolation.

In contrast, both Gauss--Legendre and Clenshaw--Curtis converge, though at
different rates for this particular integrand.

Reference: L. N. Trefethen, "Is Gauss Quadrature Better than
           Clenshaw--Curtis?", SIAM Review 50(1):67--87, 2008.

Generates Figure 15.2: newton_cotes_runge.pdf

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
=#

using LinearAlgebra
using FFTW
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
    newton_cotes_weights(n)

Compute Newton--Cotes weights for n+1 equispaced nodes on [-1, 1].
"""
function newton_cotes_weights(n)
    x = range(-1, 1, length=n+1) |> collect
    V = [x[j]^k for k in 0:n, j in 1:n+1]
    rhs = [(1 - (-1)^(k+1)) / (k + 1) for k in 0:n]
    w = V \ rhs
    return x, w
end

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

function main()
    # The Runge function and its exact integral
    f(x) = 1.0 / (1.0 + 25.0 * x^2)
    I_exact = (2.0 / 5.0) * atan(5.0)

    # Range of n values
    n_values = 4:2:50

    err_nc = zeros(length(n_values))
    err_gl = zeros(length(n_values))
    err_cc = zeros(length(n_values))

    for (i, n) in enumerate(n_values)
        # Newton--Cotes
        x_nc, w_nc = newton_cotes_weights(n)
        I_nc = dot(w_nc, f.(x_nc))
        err_nc[i] = abs(I_nc - I_exact)

        # Gauss--Legendre (n+1 points)
        x_gl, w_gl = gauss_legendre(n + 1)
        I_gl = dot(w_gl, f.(x_gl))
        err_gl[i] = abs(I_gl - I_exact)

        # Clenshaw--Curtis
        x_cc, w_cc = clenshaw_curtis_weights(n)
        I_cc = dot(w_cc, f.(x_cc))
        err_cc[i] = abs(I_cc - I_exact)
    end

    # Clamp zeros for log plot
    err_nc = max.(err_nc, 1e-16)
    err_gl = max.(err_gl, 1e-16)
    err_cc = max.(err_cc, 1e-16)

    # Create figure
    p = plot(size=(700, 500), yscale=:log10,
             xlabel=L"Number of nodes $n$",
             ylabel=L"Absolute error $|I - I_n|$",
             title=L"Quadrature of $f(x) = 1/(1 + 25x^2)$ on $[-1, 1]$",
             titlefontsize=12, guidefontsize=11, tickfontsize=10,
             legendfontsize=10, legend=:best,
             grid=true, gridalpha=0.3, gridlinewidth=0.5,
             framestyle=:box,
             xlims=(first(n_values) - 1, last(n_values) + 1))

    plot!(p, collect(n_values), err_nc, seriestype=:scatter, markershape=:circle,
          markersize=4, color=CORAL, markerstrokecolor=CORAL, linewidth=1.2,
          label="Newton-Cotes")
    plot!(p, collect(n_values), err_nc, linewidth=1.2, color=CORAL, label="")

    plot!(p, collect(n_values), err_gl, seriestype=:scatter, markershape=:rect,
          markersize=4, color=NAVY, markerstrokecolor=NAVY, linewidth=1.2,
          label="Gauss-Legendre")
    plot!(p, collect(n_values), err_gl, linewidth=1.2, color=NAVY, label="")

    plot!(p, collect(n_values), err_cc, seriestype=:scatter, markershape=:diamond,
          markersize=4, color=TEAL, markerstrokecolor=TEAL, linewidth=1.2,
          label="Clenshaw-Curtis")
    plot!(p, collect(n_values), err_cc, linewidth=1.2, color=TEAL, label="")

    # Save output
    savefig(p, joinpath(OUTPUT_DIR, "newton_cotes_runge.pdf"))
    savefig(p, joinpath(OUTPUT_DIR, "newton_cotes_runge.png"))
    println("Figure saved to $(joinpath(OUTPUT_DIR, "newton_cotes_runge.pdf"))")

    return p
end

main()
