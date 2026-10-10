#!/usr/bin/env julia
#=
quad_polynomial_exactness.jl
Chapter 15: Quadrature in Spectral Methods

Computational Etude 15.1: Polynomial Exactness Test.

For a fixed number of nodes n + 1 = 33, computes the absolute monomial
quadrature error

    |E_n(x^k)| = | sum_j w_j x_j^k - integral_{-1}^{1} x^k dx |

for k = 0, 1, ..., 2n = 64, using three classical interpolatory rules:

    1. Newton--Cotes (equispaced nodes, Vandermonde-based weights).
    2. Gauss--Legendre (Golub--Welsch algorithm).
    3. Clenshaw--Curtis (Chebyshev nodes, FFT-based weights).

The exact integral is 2/(k+1) for even k and 0 for odd k.

Theoretical degrees of precision:

    Newton--Cotes:    n   (here, k <= 32)
    Clenshaw--Curtis: n   (here, k <= 32)
    Gauss--Legendre:  2n + 1   (here, k <= 65, the entire test range)

The figure shows that NC and CC sit at machine precision for k <= n and
lift off for larger k, while Gauss stays at the floor across the whole
range. Odd k are zero by symmetry and are clamped to the floor.

Generates Figure 15.1b: polynomial_exactness.pdf

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
    n = 32                        # n + 1 = 33 quadrature points
    k_vals = collect(0:2*n)       # 0, 1, ..., 64

    # Build the three rules once
    x_nc, w_nc = newton_cotes_weights(n)
    x_gl, w_gl = gauss_legendre(n + 1)
    x_cc, w_cc = clenshaw_curtis_weights(n)

    err_nc = zeros(length(k_vals))
    err_gl = zeros(length(k_vals))
    err_cc = zeros(length(k_vals))

    for (i, k) in enumerate(k_vals)
        I_exact = iseven(k) ? 2.0 / (k + 1) : 0.0
        err_nc[i] = abs(dot(w_nc, x_nc .^ k) - I_exact)
        err_gl[i] = abs(dot(w_gl, x_gl .^ k) - I_exact)
        err_cc[i] = abs(dot(w_cc, x_cc .^ k) - I_exact)
    end

    # Clamp errors at the machine-epsilon floor for log plotting
    eps_floor = 1e-17
    err_nc = max.(err_nc, eps_floor)
    err_gl = max.(err_gl, eps_floor)
    err_cc = max.(err_cc, eps_floor)

    # Quick sanity prints
    println("n + 1 = $(n + 1) quadrature points; testing k = 0 .. $(2 * n)")
    println("max NC error for k <= n: $(maximum(err_nc[k_vals .<= n]))")
    println("max CC error for k <= n: $(maximum(err_cc[k_vals .<= n]))")
    println("max GL error for entire range: $(maximum(err_gl))")
    println("max NC error for k > n:  $(maximum(err_nc[k_vals .> n]))")
    println("max CC error for k > n:  $(maximum(err_cc[k_vals .> n]))")

    # ---------------------------------------------------------------------
    # Create figure
    # ---------------------------------------------------------------------
    p = plot(size=(800, 480), yscale=:log10,
             xlabel=L"Monomial degree $k$",
             ylabel=L"Absolute error $|E_n(x^k)|$",
             title=L"Polynomial exactness of three quadrature rules ($n + 1 = 33$ points)",
             titlefontsize=12, guidefontsize=11, tickfontsize=10,
             legendfontsize=10, legend=:bottomright,
             grid=true, gridalpha=0.3, gridlinewidth=0.5,
             framestyle=:box,
             xlims=(-1, 2*n + 2), ylims=(1e-18, 1e1),
             xticks=0:8:2*n)

    plot!(p, k_vals, err_nc, seriestype=:scatter, markershape=:circle,
          markersize=4, color=CORAL, markerstrokecolor=CORAL,
          label="Newton-Cotes")
    plot!(p, k_vals, err_nc, linewidth=1.1, color=CORAL, label="")

    plot!(p, k_vals, err_gl, seriestype=:scatter, markershape=:rect,
          markersize=4, color=NAVY, markerstrokecolor=NAVY,
          label="Gauss-Legendre")
    plot!(p, k_vals, err_gl, linewidth=1.1, color=NAVY, label="")

    plot!(p, k_vals, err_cc, seriestype=:scatter, markershape=:diamond,
          markersize=4, color=TEAL, markerstrokecolor=TEAL,
          label="Clenshaw-Curtis")
    plot!(p, k_vals, err_cc, linewidth=1.1, color=TEAL, label="")

    # Vertical guides at k = n and k = 2n+1
    vline!(p, [n], linestyle=:dash, color=RGB(0.5, 0.5, 0.5),
           linewidth=0.8, label="")
    vline!(p, [2 * n + 1], linestyle=:dash, color=RGB(0.5, 0.5, 0.5),
           linewidth=0.8, label="")
    annotate!(p, n - 0.6, 1e-2,
              text(L"NC/CC exactness ($k = n$)", 9,
                   RGB(0.5, 0.5, 0.5), :right, rotation=90))
    annotate!(p, 2 * n + 1 - 0.6, 1e-2,
              text(L"Gauss exactness ($k = 2n + 1$)", 9,
                   RGB(0.5, 0.5, 0.5), :right, rotation=90))

    # Save output
    savefig(p, joinpath(OUTPUT_DIR, "polynomial_exactness.pdf"))
    savefig(p, joinpath(OUTPUT_DIR, "polynomial_exactness.png"))
    println("Figure saved to $(joinpath(OUTPUT_DIR, "polynomial_exactness.pdf"))")

    return p
end

main()
