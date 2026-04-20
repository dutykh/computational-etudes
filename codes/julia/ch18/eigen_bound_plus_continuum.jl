# eigen_bound_plus_continuum.jl
# Chapter 18, Etude 18.6: bound states versus continuum (Pöschl-Teller).
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
# Part of "Computational Etudes: A Spectral Approach"

using CairoMakie
using Colors
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "rational_chebyshev.jl"))
include(joinpath(@__DIR__, "spectrum_verify.jl"))
using .RationalChebyshev: rational_chebyshev_derivative_matrices
using .SpectrumVerify: verify_spectrum

_json_val(x::Bool) = x ? "true" : "false"
_json_val(x::Integer) = string(x)
_json_val(x::Real) = isfinite(x) ? string(Float64(x)) : (isnan(x) ? "NaN" : (x > 0 ? "Infinity" : "-Infinity"))
_json_val(x::AbstractString) = "\"" * replace(String(x), "\\" => "\\\\", "\"" => "\\\"") * "\""
_json_val(v::AbstractVector) = "[" * join(_json_val.(v), ", ") * "]"
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
const OUTPUT_DIR = joinpath(@__DIR__, "..", "..", "..", "textbook", "figures", "ch18", "julia")
mkpath(OUTPUT_DIR)

function solve_poschl_teller(N::Int, nu::Real, L::Real = 6.0)
    _, D2, x = rational_chebyshev_derivative_matrices(N, L)
    V = -nu * (nu + 1.0) ./ cosh.(x) .^ 2
    H = -D2 + Diagonal(V)
    lam = eigvals(H)
    lam = real.(lam[abs.(imag.(lam)) .< 1e-6])
    return sort(lam)
end

function make_figure(; N1::Int = 60, N2::Int = 96, nu::Real = 4.0, tol::Real = 1e-3,
                     dump_path = nothing)
    lam1 = solve_poschl_teller(N1, nu)
    lam2 = solve_poschl_teller(N2, nu)
    expected = -collect(Int(nu):-1:1) .^ 2
    rep = verify_spectrum(lam1, lam2; tol = tol)

    fig = Figure(size = (1100, 450))

    ax1 = Axis(fig[1, 1];
        xlabel = "mode number j", ylabel = "eigenvalue E",
        title = "Pöschl-Teller V = -$(Int(nu))·$(Int(nu)+1) sech² x  (exact: $(length(expected)) bound states)",
        limits = (nothing, (min(-20.0, minimum(lam1[1:25])-1), 30.0)))
    j = collect(1:25)
    scatter!(ax1, j, lam1[1:25]; color = :white, strokecolor = NAVY,
             strokewidth = 1.1, markersize = 9, label = "N₁ = $N1")
    scatter!(ax1, j, lam2[1:25]; color = CORAL, marker = :xcross,
             markersize = 10, label = "N₂ = $N2")
    for E in expected
        hlines!(ax1, [Float64(E)]; color = TEAL, linestyle = :dash,
                linewidth = 0.6, alpha = 0.6)
    end
    hlines!(ax1, [0.0]; color = :black, linewidth = 0.6, alpha = 0.4)
    axislegend(ax1, position = :lt, labelsize = 9)

    ax2 = Axis(fig[1, 2];
        xlabel = "mode number j", ylabel = "1/δⱼ  (trusted = high)",
        yscale = log10,
        title = "drift diagnostic:  trusted = $(rep.n_trusted)",
        limits = ((0, 30), nothing))
    j = collect(1:length(rep.lam1))
    inv_ord = 1.0 ./ max.(rep.delta_ordinal, 1e-18)
    inv_nst = 1.0 ./ max.(rep.delta_nearest, 1e-18)
    scatter!(ax2, j, inv_ord; color = :white, strokecolor = NAVY,
             strokewidth = 1.0, markersize = 7, label = "1/δ_ordinal")
    scatter!(ax2, j, inv_nst; color = CORAL, marker = :xcross,
             markersize = 9, label = "1/δ_nearest")
    hlines!(ax2, [1.0 / tol]; color = :black, linestyle = :dash,
            linewidth = 0.6, alpha = 0.5)
    axislegend(ax2, position = :rt, labelsize = 9)

    save(joinpath(OUTPUT_DIR, "eigen_bound_plus_continuum.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "eigen_bound_plus_continuum.png"), fig, px_per_unit = 4)

    r = Dict{String,Any}("N1" => N1, "N2" => N2, "nu" => nu, "tol" => tol,
         "expected_bound" => Float64.(expected),
         "n_trusted" => rep.n_trusted,
         "lam1_low" => lam1[1:25], "lam2_low" => lam2[1:25],
         "trusted_first_k" => Int.(rep.trusted[1:10]))
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
    @printf("[Etude 18.6]  Pöschl-Teller with nu = %g\n", r["nu"])
    @printf("  expected bound states: %s\n", string(r["expected_bound"]))
    @printf("  trusted (drift < %g): %d\n", r["tol"], r["n_trusted"])
    println("  figure: ", joinpath(OUTPUT_DIR, "eigen_bound_plus_continuum.pdf"))
    dp !== nothing && println("  dumped: ", dp)
end
