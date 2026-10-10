# eigen_complex_detour.jl
# Chapter 18, Etude 18.11: detour around the singularity.
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
# Part of "Computational Etudes: A Spectral Approach"

using CairoMakie
using Colors
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

_json_val(x::Bool) = x ? "true" : "false"
_json_val(x::Integer) = string(x)
_json_val(x::Real) = isfinite(x) ? string(Float64(x)) : (isnan(x) ? "NaN" : (x > 0 ? "Infinity" : "-Infinity"))
_json_val(x::Complex) = "\"" * string(x) * "\""
_json_val(x::AbstractString) = "\"" * replace(String(x), "\\" => "\\\\", "\"" => "\\\"") * "\""
_json_val(v::AbstractVector) = "[" * join(_json_val.(v), ", ") * "]"
_json_val(t::Tuple) = "[" * join(_json_val.(collect(t)), ", ") * "]"
_json_val(d::AbstractDict) = "{" * join([_json_val(string(k)) * ": " * _json_val(v) for (k, v) in d], ", ") * "}"
write_json(path, x) = open(io -> print(io, _json_val(x)), path, "w")

set_theme!(Theme(fontsize = 11,
    fonts = (regular = "CMU Serif", bold = "CMU Serif Bold", italic = "CMU Serif Italic"),
    Axis = (xlabelsize = 11, ylabelsize = 11, titlesize = 11,
            xticklabelsize = 10, yticklabelsize = 10, spinewidth = 0.8,
            xtickwidth = 0.8, ytickwidth = 0.8), Legend = (labelsize = 9,)))

const NAVY = colorant"#142D6E"
const CORAL = colorant"#E74C3C"
const TEAL = colorant"#16A085"
const PURPLE = colorant"#8E44AD"
const OUTPUT_DIR = joinpath(@__DIR__, "..", "..", "..", "textbook", "figures", "ch18", "julia")
mkpath(OUTPUT_DIR)

function solve_detoured(N::Int, Delta::Real)
    D, x = cheb_matrix(N)
    y = x .+ 1im .* Delta .* (x .^ 2 .- 1)
    dy_dx = 1.0 .+ 2im .* Delta .* x
    d2y_dx2 = fill(2im * Delta, length(x))
    idx = 2:N
    D_int = D[idx, idx]
    D2_int = (D * D)[idx, idx]
    dydx_i = dy_dx[idx]
    d2ydx2_i = d2y_dx2[idx]
    y_i = y[idx]
    L = Diagonal(dydx_i .^ (-2)) * D2_int -
        Diagonal(d2ydx2_i .* dydx_i .^ (-3)) * D_int +
        Diagonal(1.0 ./ y_i)
    lam = eigvals(L)
    order = sortperm(abs.(lam))
    return lam[order], y, y_i
end

function make_figure(; Delta_scan = (0.25, 0.50, 0.75), N_demo::Int = 40,
                     N_low::Int = 20, dump_path = nothing)
    lam_ref, _, _ = solve_detoured(N_demo, 0.5)
    first6_ref = lam_ref[1:6]

    fig = Figure(size = (1150, 450))

    ax1 = Axis(fig[1, 1]; xlabel = "Re y", ylabel = "Im y",
        title = "parabolic detour y(x) = x + iΔ(x² − 1)",
        limits = ((-1.2, 1.2), (-1.0, 0.3)))
    cols = [NAVY, CORAL, TEAL]
    for (Delta, c) in zip(Delta_scan, cols)
        xd = collect(range(-1, 1; length = 401))
        yd = xd .+ 1im .* Delta .* (xd .^ 2 .- 1)
        lines!(ax1, real.(yd), imag.(yd); color = c, linewidth = 1.3,
               label = "Δ = $Delta")
    end
    hlines!(ax1, [0]; color = :black, linewidth = 0.4, alpha = 0.5)
    vlines!(ax1, [0]; color = :black, linewidth = 0.4, alpha = 0.5)
    scatter!(ax1, [0.0], [0.0]; color = PURPLE, marker = :xcross,
        markersize = 14, strokewidth = 0, label = "singularity y = 0")
    scatter!(ax1, [-1.0, 1.0], [0.0, 0.0]; color = :white,
        strokecolor = :black, strokewidth = 1.5, markersize = 10)
    axislegend(ax1, position = :rb, labelsize = 9)

    ax2 = Axis(fig[1, 2]; xlabel = "Re λ", ylabel = "Im λ",
               title = "first eigenvalues in the complex plane")
    markers = [:circle, :rect, :utriangle]
    for (Delta, c, m) in zip(Delta_scan, cols, markers)
        lam, _, _ = solve_detoured(N_demo, Delta)
        first_few = lam[1:min(7, length(lam))]
        scatter!(ax2, real.(first_few), imag.(first_few);
            color = :white, strokecolor = c, strokewidth = 1.2,
            markersize = 11, marker = m,
            label = "Δ = $Delta, N = $N_demo")
    end
    hlines!(ax2, [0]; color = :black, linewidth = 0.4, alpha = 0.5)
    vlines!(ax2, [0]; color = :black, linewidth = 0.4, alpha = 0.5)
    axislegend(ax2, position = :lb, labelsize = 9)

    save(joinpath(OUTPUT_DIR, "eigen_complex_detour.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "eigen_complex_detour.png"), fig, px_per_unit = 4)

    lam_low, _, _ = solve_detoured(N_low, 0.5)
    r = Dict("N_demo" => N_demo, "N_low" => N_low,
        "Delta_scan" => collect(Delta_scan),
        "first6_ref_Delta0p5" => [(real(z), imag(z)) for z in first6_ref],
        "low_N_first6" => [(real(z), imag(z)) for z in lam_low[1:6]])
    dump_path !== nothing && write_json(dump_path, r)
    return r
end

function parse_args()
    dp = nothing; i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--dump" && i < length(ARGS); dp = ARGS[i+1]; i += 2
        else i += 1 end
    end
    return dp
end

if abspath(PROGRAM_FILE) == @__FILE__
    dp = parse_args()
    r = make_figure(; dump_path = dp)
    println("[Etude 18.11]  complex-plane detour")
    println("  first 6 at N = $(r["N_demo"]), Δ = 0.5:")
    for (zr, zi) in r["first6_ref_Delta0p5"]
        s = zi >= 0 ? "+" : "-"
        @printf("    %+.6f %s i %.6f\n", zr, s, abs(zi))
    end
    println("  low-N check (N=$(r["N_low"]), Δ=0.5):")
    for (zr, zi) in r["low_N_first6"]
        s = zi >= 0 ? "+" : "-"
        @printf("    %+.6f %s i %.6f\n", zr, s, abs(zi))
    end
    println("  figure: ", joinpath(OUTPUT_DIR, "eigen_complex_detour.pdf"))
    dp !== nothing && println("  dumped: ", dp)
end
