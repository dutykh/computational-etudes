#!/usr/bin/env julia
#=
quad_aliasing_chebyshev.jl
Chapter 15: Quadrature in Spectral Methods

Computational Etude 15.6: Aliasing and Chebyshev Quadrature.

Illustrates the aliasing phenomenon in Clenshaw--Curtis quadrature by
examining the quadrature error E_n(T_k) = |I(T_k) - I_n(T_k)| for
individual Chebyshev polynomials T_k, and its impact on the integration
of |x|^3.

Left panel: For n = 30, the CC quadrature error |E_n(T_k)| is plotted
against even k from 0 to 200. Three regimes are visible:
    (1) k <= n: errors are at machine zero (exact integration),
    (2) n < k <= 2n (shaded): errors are small, O(n^{-2}),
    (3) k > 2n: aliasing errors of O(1) at multiples of 2n.
This structure underlies Trefethen's (2008) explanation of why CC
quadrature is nearly as accurate as Gauss for most integrands.

Right panel: For f(x) = |x|^3, the Chebyshev expansion coefficients
|a_k| (which decay algebraically as k^{-4}) are plotted alongside the
aliased error contributions |a_k * E_n(T_k)|. The product decays so
rapidly that the total aliasing error is negligible.

Reference: L. N. Trefethen, "Is Gauss Quadrature Better than
           Clenshaw--Curtis?", SIAM Review 50(1):67--87, 2008, Figure 3.

Generates Figure 15.6: aliasing_chebyshev.pdf

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
    exact_integral_Tk(k)

Exact integral of T_k(x) over [-1, 1].
I(T_k) = 2 / (1 - k^2) for even k, 0 for odd k.
"""
function exact_integral_Tk(k)
    if k % 2 == 1
        return 0.0
    end
    return 2.0 / (1.0 - k^2)
end

"""
    chebyshev_coefficients(f, N)

Compute Chebyshev expansion coefficients a_k of f(x) using the DCT
via evaluation at Chebyshev points.
"""
function chebyshev_coefficients(f, N)
    j = 0:N-1
    x = cos.(π .* j ./ (N - 1))  # Chebyshev extremal points
    fx = f.(x)

    # DCT-I via FFT
    v = vcat(fx, fx[end-1:-1:2])
    coeffs = real.(fft(v))[1:N] ./ (N - 1)
    coeffs[1] /= 2
    coeffs[end] /= 2

    return coeffs
end

function main()
    n = 30       # CC quadrature with n+1 = 31 points
    k_max = 200  # Maximum Chebyshev degree to examine

    # Compute CC nodes and weights
    x_cc, w_cc = clenshaw_curtis_weights(n)

    # -------------------------------------------------------------------------
    # Left panel: |E_n(T_k)| for even k = 0, 2, 4, ..., k_max
    # -------------------------------------------------------------------------
    k_even = 0:2:k_max
    errors = zeros(length(k_even))

    for (i, k) in enumerate(k_even)
        # Evaluate T_k at CC nodes: T_k(x) = cos(k * arccos(x))
        Tk_at_nodes = cos.(k .* acos.(x_cc))

        # Numerical integral via CC quadrature
        I_n = dot(w_cc, Tk_at_nodes)

        # Exact integral
        I_exact = exact_integral_Tk(k)

        errors[i] = abs(I_n - I_exact)
    end

    # -------------------------------------------------------------------------
    # Right panel: Chebyshev coefficients of |x|^3 and aliased errors
    # -------------------------------------------------------------------------
    f_abs3(x) = abs(x)^3

    # Compute many Chebyshev coefficients using a large number of points
    N_coeffs = 512
    a_k = chebyshev_coefficients(f_abs3, N_coeffs)

    # Only even coefficients are nonzero for this even function
    k_coeffs = 0:2:min(k_max, N_coeffs - 1)
    a_even = abs.(a_k[k_coeffs .+ 1])

    # Compute aliased error contributions: |a_k * E_n(T_k)|
    n_common = min(length(k_coeffs), length(k_even))
    k_common = collect(k_even)[1:n_common]
    a_common = zeros(n_common)
    e_common = zeros(n_common)

    for (i, k) in enumerate(k_common)
        if k + 1 <= N_coeffs
            a_common[i] = abs(a_k[k + 1])
        end
        # Compute E_n(T_k)
        Tk_at_nodes = cos.(k .* acos.(x_cc))
        I_n = dot(w_cc, Tk_at_nodes)
        I_exact = exact_integral_Tk(k)
        e_common[i] = abs(I_n - I_exact)
    end

    aliased_product = a_common .* e_common

    # -------------------------------------------------------------------------
    # Create figure
    # -------------------------------------------------------------------------
    floor_val = 1e-17

    # --- Left panel: CC quadrature errors for T_k ---
    errors_plot = max.(errors, floor_val)

    p1 = plot(yscale=:log10,
              xlabel=L"Even degree $k$",
              ylabel=L"|E_n(T_k)|",
              title=L"(a) CC error for Chebyshev polynomials ($n = 30$)",
              titlefontsize=11, guidefontsize=11, tickfontsize=10,
              legendfontsize=9, legend=:bottomright,
              grid=true, gridalpha=0.3, gridlinewidth=0.5,
              framestyle=:box,
              xlims=(-2, k_max + 5), ylims=(1e-17, 1e1))

    # Shade the region n <= k <= 2n
    vspan!(p1, [n, 2n], alpha=0.12, color=SKY, label=L"n \leq k \leq 2n")

    # Vertical dashed lines
    vline!(p1, [n], color=CORAL, linestyle=:dash, linewidth=0.8, alpha=0.7, label="")
    vline!(p1, [2n], color=CORAL, linestyle=:dash, linewidth=0.8, alpha=0.7, label="")

    scatter!(p1, collect(k_even), errors_plot, markershape=:circle,
             markersize=3, color=NAVY, markerstrokecolor=NAVY, markerstrokewidth=0.5,
             label="")

    # Annotate the regions
    annotate!(p1, n / 2, 2e-16, text("exact", 9, TEAL, :center))
    annotate!(p1, 1.5 * n, 1e-3, text("small\nerrors", 9, NAVY, :center))
    annotate!(p1, 2.5 * n + 20, 0.5, text("aliasing", 9, CORAL, :center))

    # --- Right panel: |a_k| and |a_k * E_n(T_k)| for |x|^3 ---
    a_plot = max.(a_even, floor_val)
    alias_plot = max.(aliased_product[1:length(k_coeffs)], floor_val)

    p2 = plot(yscale=:log10,
              xlabel=L"Even degree $k$",
              ylabel="Magnitude",
              title=L"(b) Aliasing for $f(x) = |x|^3$",
              titlefontsize=11, guidefontsize=11, tickfontsize=10,
              legendfontsize=9, legend=:topright,
              grid=true, gridalpha=0.3, gridlinewidth=0.5,
              framestyle=:box,
              xlims=(-2, k_max + 5), ylims=(1e-17, 1e1))

    # Shade the region n <= k <= 2n
    vspan!(p2, [n, 2n], alpha=0.12, color=SKY, label="")

    scatter!(p2, collect(k_coeffs), a_plot, markershape=:circle,
             markersize=3, color=NAVY, markerstrokecolor=NAVY, markerstrokewidth=0.5,
             label=L"|a_k|$ (Chebyshev coefficients)$")

    scatter!(p2, k_common[1:length(k_coeffs)], alias_plot,
             markershape=:rect, markersize=3, color=:white,
             markerstrokecolor=CORAL, markerstrokewidth=0.7,
             label=L"|a_k \cdot E_n(T_k)|$ (aliased error)$")

    # Combine into two-panel figure
    p = plot(p1, p2, layout=(1, 2), size=(1100, 450))

    # Save output
    savefig(p, joinpath(OUTPUT_DIR, "aliasing_chebyshev.pdf"))
    savefig(p, joinpath(OUTPUT_DIR, "aliasing_chebyshev.png"))
    println("Figure saved to $(joinpath(OUTPUT_DIR, "aliasing_chebyshev.pdf"))")

    return p
end

main()
