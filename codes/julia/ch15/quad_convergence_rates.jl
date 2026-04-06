#!/usr/bin/env julia
#=
quad_convergence_rates.jl
Chapter 15: Quadrature in Spectral Methods

Computational Etude 15.10: Asymptotic Convergence Rate Verification.

Verifies the theoretical asymptotic convergence rates for two quadrature
strategies applied to the integral

    I = int_{-inf}^{inf} e^{-x^2} cos(x^3) dx

using the rescaling technique of Trefethen and Weideman (2014).

Left panel: Gauss--Hermite errors plotted against n^{1/2}. Theory
predicts that the error decays as O(exp(-C * n^{1/2})), which
corresponds to "root-exponential" or "sub-geometric" convergence.
On a semilog plot with the horizontal axis scaled as sqrt(n), the
data should fall approximately on a straight line. A linear fit
confirms the predicted rate.

Right panel: Truncated Gauss--Legendre errors plotted against n^{2/3}.
For a Gauss--Legendre rule on the interval [-L * n^{1/3}, L * n^{1/3}]
(where the truncation width grows with n to capture more of the
integrand), theory predicts convergence at rate O(exp(-C * n^{2/3})).
On a semilog plot with the horizontal axis scaled as n^{2/3}, the
data should be approximately linear. This 2/3-power rate is
intermediate between root-exponential and geometric convergence.

Reference: L. N. Trefethen and J. A. C. Weideman, "The exponentially
           convergent trapezoidal rule", SIAM Review 56(3):385--458, 2014.

Generates Figure 15.10: convergence_rates.pdf

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
    # -------------------------------------------------------------------------
    # Reference value via high-accuracy adaptive quadrature
    # -------------------------------------------------------------------------
    integrand(x) = exp(-x^2) * cos(x^3)
    I_exact, _ = quadgk(integrand, -Inf, Inf, rtol=1e-15, atol=1e-15)
    println("Reference integral value: $I_exact")

    # -------------------------------------------------------------------------
    # Gauss--Hermite errors for n = 10, 15, ..., 500
    # -------------------------------------------------------------------------
    n_gh = 10:5:500
    err_gh = zeros(length(n_gh))

    for (i, n) in enumerate(n_gh)
        x_h, w_h = gauss_hermite(n)
        I_gh = dot(w_h, cos.(x_h .^ 3))
        err_gh[i] = abs(I_gh - I_exact)
    end

    # -------------------------------------------------------------------------
    # Truncated Gauss--Legendre on [-L * n^{1/3}, L * n^{1/3}]
    # The full integrand e^{-x^2} cos(x^3) is used (no weight separation).
    # -------------------------------------------------------------------------
    L = 3.0  # Scaling parameter
    n_gl = 10:5:500
    err_gl = zeros(length(n_gl))

    for (i, n) in enumerate(n_gl)
        # Truncation half-width grows as n^{1/3}
        half_width = L * n^(1.0 / 3.0)

        # Gauss--Legendre nodes on [-1, 1], mapped to [-half_width, half_width]
        s, ws = gauss_legendre(n)
        x_mapped = half_width .* s
        w_mapped = half_width .* ws

        # Evaluate the full integrand
        f_vals = exp.(-x_mapped .^ 2) .* cos.(x_mapped .^ 3)
        I_gl = dot(w_mapped, f_vals)
        err_gl[i] = abs(I_gl - I_exact)
    end

    # -------------------------------------------------------------------------
    # Prepare data for plotting: filter out zero errors
    # -------------------------------------------------------------------------
    floor_val = 1e-16

    # Gauss--Hermite: use errors above floor for fitting
    sqrt_n_gh = sqrt.(collect(n_gh))
    log_err_gh = log10.(max.(err_gh, floor_val))
    mask_gh = err_gh .> floor_val

    # Truncated GL: use errors above floor for fitting
    n23_gl = collect(n_gl) .^ (2.0 / 3.0)
    log_err_gl = log10.(max.(err_gl, floor_val))
    mask_gl = err_gl .> floor_val

    # Linear fits on the valid (above-floor) data
    fit_gh = nothing
    slope_gh = nothing
    if sum(mask_gh) > 2
        # polyfit: degree 1
        x_fit = sqrt_n_gh[mask_gh]
        y_fit = log_err_gh[mask_gh]
        # Least squares: y = a*x + b
        A = hcat(x_fit, ones(length(x_fit)))
        coeffs_gh = A \ y_fit
        slope_gh = coeffs_gh[1]
        fit_gh = coeffs_gh[1] .* sqrt_n_gh .+ coeffs_gh[2]
    end

    fit_gl = nothing
    slope_gl = nothing
    if sum(mask_gl) > 2
        x_fit = n23_gl[mask_gl]
        y_fit = log_err_gl[mask_gl]
        A = hcat(x_fit, ones(length(x_fit)))
        coeffs_gl = A \ y_fit
        slope_gl = coeffs_gl[1]
        fit_gl = coeffs_gl[1] .* n23_gl .+ coeffs_gl[2]
    end

    # -------------------------------------------------------------------------
    # Create figure
    # -------------------------------------------------------------------------

    # --- Left panel: Gauss--Hermite vs sqrt(n) ---
    p1 = plot(yscale=:log10,
              xlabel=L"n^{1/2}",
              ylabel=L"Absolute error $|I - I_n|$",
              title=L"(a) Gauss-Hermite: $O(e^{-C\sqrt{n}})$",
              titlefontsize=11, guidefontsize=11, tickfontsize=10,
              legendfontsize=9, legend=:topright,
              grid=true, gridalpha=0.3, gridlinewidth=0.5,
              framestyle=:box,
              ylims=(1e-16, 1e1))

    scatter!(p1, sqrt_n_gh, max.(err_gh, floor_val), markershape=:circle,
             markersize=3, color=NAVY, markerstrokecolor=NAVY,
             markerstrokewidth=0.5, label="Gauss-Hermite error")

    if fit_gh !== nothing
        plot!(p1, sqrt_n_gh, 10.0 .^ fit_gh, linestyle=:dash, color=CORAL,
              linewidth=1.2,
              label="Linear fit (slope = $(round(slope_gh, digits=2)))")
    end

    hline!(p1, [eps(Float64)], color=:gray, linestyle=:dot, linewidth=0.8,
           alpha=0.7, label="")

    # --- Right panel: Truncated GL vs n^{2/3} ---
    p2 = plot(yscale=:log10,
              xlabel=L"n^{2/3}",
              ylabel=L"Absolute error $|I - I_n|$",
              title=L"(b) Truncated GL: $O(e^{-C\, n^{2/3}})$",
              titlefontsize=11, guidefontsize=11, tickfontsize=10,
              legendfontsize=9, legend=:topright,
              grid=true, gridalpha=0.3, gridlinewidth=0.5,
              framestyle=:box,
              ylims=(1e-16, 1e1))

    scatter!(p2, n23_gl, max.(err_gl, floor_val), markershape=:circle,
             markersize=3, color=TEAL, markerstrokecolor=TEAL,
             markerstrokewidth=0.5, label="Truncated GL error")

    if fit_gl !== nothing
        plot!(p2, n23_gl, 10.0 .^ fit_gl, linestyle=:dash, color=CORAL,
              linewidth=1.2,
              label="Linear fit (slope = $(round(slope_gl, digits=2)))")
    end

    hline!(p2, [eps(Float64)], color=:gray, linestyle=:dot, linewidth=0.8,
           alpha=0.7, label="")

    # Combine into two-panel figure
    p = plot(p1, p2, layout=(1, 2), size=(1100, 450))

    # Save output
    savefig(p, joinpath(OUTPUT_DIR, "convergence_rates.pdf"))
    savefig(p, joinpath(OUTPUT_DIR, "convergence_rates.png"))
    println("Figure saved to $(joinpath(OUTPUT_DIR, "convergence_rates.pdf"))")

    return p
end

main()
