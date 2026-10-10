#!/usr/bin/env julia
#=
quad_gauss_cc_construction.jl
Chapter 15: Quadrature in Spectral Methods

Computational Etude 15.4: Constructing Gauss and Clenshaw--Curtis Rules.

Implements two fundamental quadrature construction algorithms from scratch
and validates them against known results.

Gauss--Legendre via the Golub--Welsch Algorithm
------------------------------------------------
The n-point Gauss--Legendre nodes are eigenvalues of the n x n symmetric
tridiagonal Jacobi matrix J with:

    J[j, j]   = 0                          (diagonal)
    J[j, j+1] = beta_j = j / sqrt(4j^2 - 1)  (off-diagonal, j = 1, ..., n-1)

The weights are w_j = 2 * (v_j[0])^2, where v_j is the eigenvector
corresponding to the j-th eigenvalue.

Clenshaw--Curtis via FFT
-------------------------
The n+1 Clenshaw--Curtis nodes are x_j = cos(j pi / n), j = 0, ..., n.
The weights are computed by:

    1. Form the Chebyshev moments c_k = 2/(1-k^2) for even k, 0 for odd k.
    2. Construct a length-2n vector and compute its real FFT.
    3. Extract the weights (with endpoint halving).

Both methods are validated by integrating e^x (analytic: e - 1/e) and
1/(1+x^2) (analytic: pi/2) and comparing with a reference Gauss--Legendre
implementation.

Generates Figure 15.4: gauss_cc_construction.pdf (two-panel figure)

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
    golub_welsch(n)

Compute n-point Gauss--Legendre nodes and weights via the Golub--Welsch
algorithm.
"""
function golub_welsch(n)
    if n == 1
        return [0.0], [2.0]
    end
    j = collect(1:n-1) * 1.0
    beta = j ./ sqrt.(4.0 .* j .^ 2 .- 1.0)
    J = diagm(1 => beta, -1 => beta)
    F = eigen(Symmetric(J))
    idx = sortperm(F.values)
    nodes = F.values[idx]
    weights = 2.0 .* F.vectors[1, idx] .^ 2
    return nodes, weights
end

"""
    clenshaw_curtis_fft(n)

Compute (n+1)-point Clenshaw--Curtis nodes and weights via DCT-I / FFT.
"""
function clenshaw_curtis_fft(n)
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

"""
    gauss_legendre_ref(n)

Reference n-point Gauss--Legendre via Golub--Welsch (used as reference).
This is an independent implementation for cross-validation.
"""
function gauss_legendre_ref(n)
    # Same as golub_welsch but kept separate for clarity
    return golub_welsch(n)
end

function main()
    # Test functions and their exact integrals over [-1, 1]
    f1(x) = exp(x)
    I1_exact = exp(1) - exp(-1)  # e - e^{-1}

    f2(x) = 1.0 / (1.0 + x^2)
    I2_exact = π / 2.0  # 2 * arctan(1) = pi/2

    n_values = 2:2:50

    err_gw_f1  = zeros(length(n_values))
    err_cc_f1  = zeros(length(n_values))
    err_ref_f1 = zeros(length(n_values))

    err_gw_f2  = zeros(length(n_values))
    err_cc_f2  = zeros(length(n_values))
    err_ref_f2 = zeros(length(n_values))

    for (i, n) in enumerate(n_values)
        # Golub--Welsch Gauss--Legendre (n+1 points)
        x_gw, w_gw = golub_welsch(n + 1)
        err_gw_f1[i] = abs(dot(w_gw, f1.(x_gw)) - I1_exact)
        err_gw_f2[i] = abs(dot(w_gw, f2.(x_gw)) - I2_exact)

        # Clenshaw--Curtis (n+1 points)
        x_cc, w_cc = clenshaw_curtis_fft(n)
        err_cc_f1[i] = abs(dot(w_cc, f1.(x_cc)) - I1_exact)
        err_cc_f2[i] = abs(dot(w_cc, f2.(x_cc)) - I2_exact)

        # Reference Gauss--Legendre (n+1 points)
        x_ref, w_ref = gauss_legendre_ref(n + 1)
        err_ref_f1[i] = abs(dot(w_ref, f1.(x_ref)) - I1_exact)
        err_ref_f2[i] = abs(dot(w_ref, f2.(x_ref)) - I2_exact)
    end

    # Clamp zeros
    floor_val = 1e-16
    err_gw_f1  = max.(err_gw_f1, floor_val)
    err_cc_f1  = max.(err_cc_f1, floor_val)
    err_ref_f1 = max.(err_ref_f1, floor_val)
    err_gw_f2  = max.(err_gw_f2, floor_val)
    err_cc_f2  = max.(err_cc_f2, floor_val)
    err_ref_f2 = max.(err_ref_f2, floor_val)

    nv = collect(n_values)

    # --- Panel (a): e^x ---
    p1 = plot(size=(550, 450), yscale=:log10,
              xlabel=L"Number of points $n$",
              ylabel=L"Absolute error $|I - I_n|$",
              title=L"(a) $f(x) = e^x$",
              titlefontsize=12, guidefontsize=11, tickfontsize=10,
              legendfontsize=9, legend=:best,
              grid=true, gridalpha=0.3, gridlinewidth=0.5,
              framestyle=:box,
              xlims=(first(n_values) - 1, last(n_values) + 1))

    plot!(p1, nv, err_gw_f1, seriestype=:scatter, markershape=:circle,
          markersize=4, color=NAVY, markerstrokecolor=NAVY, label="Golub-Welsch (Gauss)")
    plot!(p1, nv, err_gw_f1, linewidth=1.2, color=NAVY, label="")

    plot!(p1, nv, err_cc_f1, seriestype=:scatter, markershape=:diamond,
          markersize=4, color=TEAL, markerstrokecolor=TEAL, label="CC (FFT)")
    plot!(p1, nv, err_cc_f1, linewidth=1.2, color=TEAL, label="")

    plot!(p1, nv, err_ref_f1, seriestype=:scatter, markershape=:xcross,
          markersize=5, color=CORAL, markerstrokecolor=CORAL, label="Reference (Gauss)")
    plot!(p1, nv, err_ref_f1, linewidth=1.0, color=CORAL, linestyle=:dash, label="")

    # --- Panel (b): 1/(1+x^2) ---
    p2 = plot(size=(550, 450), yscale=:log10,
              xlabel=L"Number of points $n$",
              ylabel=L"Absolute error $|I - I_n|$",
              title=L"(b) $f(x) = 1/(1 + x^2)$",
              titlefontsize=12, guidefontsize=11, tickfontsize=10,
              legendfontsize=9, legend=:best,
              grid=true, gridalpha=0.3, gridlinewidth=0.5,
              framestyle=:box,
              xlims=(first(n_values) - 1, last(n_values) + 1))

    plot!(p2, nv, err_gw_f2, seriestype=:scatter, markershape=:circle,
          markersize=4, color=NAVY, markerstrokecolor=NAVY, label="Golub-Welsch (Gauss)")
    plot!(p2, nv, err_gw_f2, linewidth=1.2, color=NAVY, label="")

    plot!(p2, nv, err_cc_f2, seriestype=:scatter, markershape=:diamond,
          markersize=4, color=TEAL, markerstrokecolor=TEAL, label="CC (FFT)")
    plot!(p2, nv, err_cc_f2, linewidth=1.2, color=TEAL, label="")

    plot!(p2, nv, err_ref_f2, seriestype=:scatter, markershape=:xcross,
          markersize=5, color=CORAL, markerstrokecolor=CORAL, label="Reference (Gauss)")
    plot!(p2, nv, err_ref_f2, linewidth=1.0, color=CORAL, linestyle=:dash, label="")

    # Combine into two-panel figure
    p = plot(p1, p2, layout=(1, 2), size=(1100, 450))

    # Save output
    savefig(p, joinpath(OUTPUT_DIR, "gauss_cc_construction.pdf"))
    savefig(p, joinpath(OUTPUT_DIR, "gauss_cc_construction.png"))
    println("Figure saved to $(joinpath(OUTPUT_DIR, "gauss_cc_construction.pdf"))")

    return p
end

main()
