#!/usr/bin/env julia
#=
quad_complex_plane.jl
Chapter 15: Quadrature in Spectral Methods

Computational Etude 15.7: Quadrature in the Complex Plane.

Reproduces the essential content of Figure 4 from Trefethen (2008),
comparing Newton--Cotes, Gauss--Legendre, and Clenshaw--Curtis
quadrature rules through the lens of rational approximation in the
complex plane.

The key idea is that a quadrature rule with nodes x_k and weights w_k
defines a rational approximant

    r_n(z) = sum_{k} w_k / (z - x_k)

to the function

    phi(z) = log((z + 1) / (z - 1))

which is the Cauchy integral representation of the constant function 1
over [-1, 1]. The contour plot of |phi(z) - r_n(z)| reveals the
accuracy of each quadrature rule throughout the complex plane.

For Newton--Cotes (equispaced nodes), the approximation is poor near
the endpoints of [-1, 1], reflecting the Runge phenomenon. For
Gauss--Legendre and Clenshaw--Curtis, the approximation is far more
uniform, with CC showing slightly larger errors near the real axis
but comparable accuracy overall.

The figure shows a 1 x 3 panel layout for n = 16 (17 points):
    (a) Newton--Cotes
    (b) Gauss--Legendre
    (c) Clenshaw--Curtis

Reference: L. N. Trefethen, "Is Gauss Quadrature Better than
           Clenshaw--Curtis?", SIAM Review 50(1):67--87, 2008, Figure 4.

Generates Figure 15.7: complex_plane.pdf

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

"""
    phi(z)

Compute phi(z) = log((z+1)/(z-1)) using the principal branch.
This is the Cauchy integral of the constant function 1 over [-1, 1].
"""
phi(z) = log((z + 1.0) / (z - 1.0))

"""
    rational_approx(z, x, w)

Compute the rational approximant r_n(z) = sum_k w_k / (z - x_k).
"""
function rational_approx(z, x_q, w_q)
    r = zeros(Complex{Float64}, size(z))
    for k in eachindex(x_q)
        r .+= w_q[k] ./ (z .- x_q[k])
    end
    return r
end

function main()
    n = 16  # 17-point rules

    # Compute nodes and weights for all three rules
    x_nc, w_nc = newton_cotes_weights(n)
    x_gl, w_gl = gauss_legendre(n + 1)
    x_cc, w_cc = clenshaw_curtis_weights(n)

    rules = [
        (x_nc, w_nc, "(a) Newton-Cotes"),
        (x_gl, w_gl, "(b) Gauss-Legendre"),
        (x_cc, w_cc, "(c) Clenshaw-Curtis"),
    ]

    # Complex grid
    nx, ny = 400, 300
    x_re = range(-2.5, 2.5, length=nx)
    y_im = range(-2.0, 2.0, length=ny)
    X = [xr for yi in y_im, xr in x_re]
    Y = [yi for yi in y_im, xr in x_re]
    Z = X .+ im .* Y

    # Mask points too close to the branch cut [-1, 1] on the real axis
    mask_distance = 0.03
    on_cut = (abs.(Y) .< mask_distance) .& (X .>= -1.0) .& (X .<= 1.0)

    # Contour levels: 10^{-13}, 10^{-12}, ..., 10^0
    levels = 10.0 .^ (-13:0)

    # Compute phi(z) once
    phi_vals = phi.(Z)

    # Create figure panels
    panels = []

    for (x_q, w_q, title) in rules
        # Compute rational approximation error
        r_vals = rational_approx(Z, x_q, w_q)
        error_vals = abs.(phi_vals .- r_vals)

        # Mask branch cut region
        error_masked = copy(error_vals)
        error_masked[on_cut] .= NaN

        # Clamp for log scale
        error_masked = max.(error_masked, 1e-15)

        pi = contour(collect(x_re), collect(y_im), log10.(error_masked),
                     levels=collect(-13.0:1.0:0.0),
                     color=:grays, linewidth=0.5,
                     title=title, titlefontsize=11,
                     xlabel=L"\mathrm{Re}(z)", ylabel=L"\mathrm{Im}(z)",
                     guidefontsize=10, tickfontsize=9,
                     xlims=(-2.5, 2.5), ylims=(-2.0, 2.0),
                     aspect_ratio=:equal, colorbar=false,
                     grid=true, gridalpha=0.15, gridlinewidth=0.5,
                     framestyle=:box, fill=false)

        # Mark the interval [-1, 1]
        plot!(pi, [-1, 1], [0, 0], color=NAVY, linewidth=2.0, label="")

        # Mark quadrature nodes
        scatter!(pi, x_q, zeros(length(x_q)), markershape=:circle,
                 markersize=3, color=CORAL, markerstrokecolor=CORAL,
                 markerstrokewidth=0.5, label="")

        push!(panels, pi)
    end

    p = plot(panels..., layout=(1, 3), size=(1300, 450))

    # Save output
    savefig(p, joinpath(OUTPUT_DIR, "complex_plane.pdf"))
    savefig(p, joinpath(OUTPUT_DIR, "complex_plane.png"))
    println("Figure saved to $(joinpath(OUTPUT_DIR, "complex_plane.pdf"))")

    return p
end

main()
