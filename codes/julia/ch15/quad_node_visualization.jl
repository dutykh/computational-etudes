#!/usr/bin/env julia
#=
quad_node_visualization.jl
Chapter 15: Quadrature in Spectral Methods

Computational Etude 15.1: Quadrature Node Distributions.

Visualises three classical sets of quadrature nodes on [-1, 1] for n = 32
points, displayed on parallel horizontal lines:

    1. Newton--Cotes (equispaced) nodes: x_j = -1 + 2j/n, j = 0, ..., n.
    2. Gauss--Legendre nodes: zeros of P_{n+1}(x).
    3. Clenshaw--Curtis (Chebyshev) nodes: x_j = cos(j pi / n), j = 0, ..., n.

The figure highlights the fundamental difference in node clustering near
the endpoints: equispaced nodes remain uniformly distributed, whereas both
Gauss--Legendre and Clenshaw--Curtis nodes cluster near x = +-1. This
clustering is essential for stable high-order polynomial interpolation and
for the convergence of spectral quadrature rules.

Reference: L. N. Trefethen, "Is Gauss Quadrature Better than
           Clenshaw--Curtis?", SIAM Review 50(1):67--87, 2008, Figure 1.

Generates Figure 15.1: node_visualization.pdf

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

function main()
    n = 32  # Number of intervals (n+1 points)

    # Newton--Cotes (equispaced) nodes
    nc_nodes = range(-1, 1, length=n+1) |> collect

    # Gauss--Legendre nodes (n+1 points)
    gl_nodes, _ = gauss_legendre(n + 1)

    # Clenshaw--Curtis (Chebyshev extremal) nodes
    cc_nodes = [cos(j * pi / n) for j in 0:n]

    # Vertical positions for the three rows
    y_nc = 3.0
    y_gl = 2.0
    y_cc = 1.0

    # Create figure
    p = plot(size=(800, 300), xlims=(-1.8, 1.15), ylims=(0.4, 3.6),
             xlabel=L"x", xticks=[-1, -0.5, 0, 0.5, 1],
             yticks=nothing, yshowaxis=false,
             title=L"Quadrature nodes on $[-1, 1]$ with $n = 32$",
             titlefontsize=12, guidefontsize=11, tickfontsize=10,
             legend=false, grid=false, framestyle=:box,
             left_margin=10Plots.mm, bottom_margin=5Plots.mm,
             top_margin=5Plots.mm, right_margin=5Plots.mm)

    # Plot horizontal base lines
    for y in [y_nc, y_gl, y_cc]
        plot!(p, [-1, 1], [y, y], color=:gray, linewidth=0.5)
    end

    # Plot nodes as filled circles
    ms = 4

    scatter!(p, nc_nodes, fill(y_nc, n+1), color=CORAL, markersize=ms,
             markerstrokecolor=CORAL, markerstrokewidth=0.5)

    scatter!(p, gl_nodes, fill(y_gl, n+1), color=NAVY, markersize=ms,
             markerstrokecolor=NAVY, markerstrokewidth=0.5)

    scatter!(p, cc_nodes, fill(y_cc, n+1), color=TEAL, markersize=ms,
             markerstrokecolor=TEAL, markerstrokewidth=0.5)

    # Labels on the left
    label_x = -1.25
    annotate!(p, label_x, y_nc, text("Newton-Cotes", 10, :right, CORAL))
    annotate!(p, label_x, y_gl, text("Gauss-Legendre", 10, :right, NAVY))
    annotate!(p, label_x, y_cc, text("Clenshaw-Curtis", 10, :right, TEAL))

    # Mark the interval endpoints
    for y in [y_nc, y_gl, y_cc]
        scatter!(p, [-1.0, 1.0], [y, y], markershape=:vline, color=:black,
                 markersize=6, markerstrokewidth=1.0)
    end

    # Save output
    savefig(p, joinpath(OUTPUT_DIR, "node_visualization.pdf"))
    savefig(p, joinpath(OUTPUT_DIR, "node_visualization.png"))
    println("Figure saved to $(joinpath(OUTPUT_DIR, "node_visualization.pdf"))")

    return p
end

main()
