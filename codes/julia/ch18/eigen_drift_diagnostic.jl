# eigen_drift_diagnostic.jl
#
# Chapter 18: Linear Spectral Eigenproblems
# Computational Etude 18.5: build a spectrum lie detector.
#
# Applies the drift-with-N diagnostic (spectrum_verify) to:
#   (A) Dirichlet Laplacian u_xx + lambda u = 0 on (-1, 1)
#   (B) Harmonic oscillator u_xx + (lambda - x^2) u = 0 on the real line.
#
# Reproduces Boyd (2000) Fig 7.7.
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
# Part of "Computational Etudes: A Spectral Approach"

using CairoMakie
using Colors
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))
include(joinpath(@__DIR__, "rational_chebyshev.jl"))
include(joinpath(@__DIR__, "spectrum_verify.jl"))
using .RationalChebyshev: rational_chebyshev_derivative_matrices
using .SpectrumVerify: verify_spectrum, SpectrumReport

_json_val(x::Bool)          = x ? "true" : "false"
_json_val(x::Integer)       = string(x)
_json_val(x::Real)          = isfinite(x) ? string(Float64(x)) : (isnan(x) ? "NaN" : (x > 0 ? "Infinity" : "-Infinity"))
_json_val(x::AbstractString)= "\"" * replace(String(x), "\\" => "\\\\", "\"" => "\\\"") * "\""
_json_val(v::AbstractVector)= "[" * join(_json_val.(v), ", ") * "]"
function _json_val(d::AbstractDict)
    pairs = String[]
    for (k, v) in d; push!(pairs, _json_val(string(k)) * ": " * _json_val(v)); end
    return "{" * join(pairs, ", ") * "}"
end
write_json(path, x) = open(io -> print(io, _json_val(x)), path, "w")

set_theme!(Theme(
    fontsize = 11,
    fonts = (regular = "CMU Serif", bold = "CMU Serif Bold", italic = "CMU Serif Italic"),
    Axis = (xlabelsize = 11, ylabelsize = 11, titlesize = 11,
            xticklabelsize = 10, yticklabelsize = 10,
            spinewidth = 0.8, xtickwidth = 0.8, ytickwidth = 0.8),
    Legend = (labelsize = 9,),
))

const NAVY  = colorant"#142D6E"
const CORAL = colorant"#E74C3C"

const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch18", "julia")
mkpath(OUTPUT_DIR)

function solve_laplacian(N::Int)
    D, _ = cheb_matrix(N)
    A = -(D * D)[2:N, 2:N]
    return sort(real.(eigvals(A)))
end

function solve_oscillator(N::Int, ell::Real = 4.0)
    _, D2, x = rational_chebyshev_derivative_matrices(N, ell)
    H = -D2 + Diagonal(x .^ 2)
    lam = eigvals(H)
    keep = abs.(imag.(lam)) .< 1e-6
    lam = real.(lam[keep])
    lam = lam[lam .> 0]
    return sort(lam)
end

function plot_panel!(ax, rep::SpectrumReport, title_str::String, tol::Real)
    j = collect(1:length(rep.lam1))
    inv_ord = 1.0 ./ max.(rep.delta_ordinal, 1e-18)
    inv_nst = 1.0 ./ max.(rep.delta_nearest, 1e-18)
    scatter!(ax, j, inv_ord; color = :white, strokecolor = NAVY,
             strokewidth = 1.0, markersize = 7,
             label = "1/δ_ordinal")
    scatter!(ax, j, inv_nst; color = CORAL, marker = :xcross,
             markersize = 9,
             label = "1/δ_nearest")
    hlines!(ax, [1.0 / tol]; color = :black, linestyle = :dash,
            linewidth = 0.6, alpha = 0.5)
    ax.title[] = title_str
end

function make_figure(; N1::Int = 32, N2::Int = 48, tol::Real = 1e-3,
                     dump_path = nothing)
    lam1_A = solve_laplacian(N1);  lam2_A = solve_laplacian(N2)
    rep_A  = verify_spectrum(lam1_A, lam2_A; tol = tol)

    lam1_B = solve_oscillator(N1); lam2_B = solve_oscillator(N2)
    rep_B  = verify_spectrum(lam1_B, lam2_B; tol = tol)

    fig = Figure(size = (1100, 450))

    ax1 = Axis(fig[1, 1];
        xlabel = "mode number j", ylabel = "1 / δ_j  (trusted = high)",
        yscale = log10,
        title = "Dirichlet Laplacian,  N₁ = $N1, N₂ = $N2")
    plot_panel!(ax1, rep_A, "Dirichlet Laplacian,  N₁ = $N1, N₂ = $N2", tol)
    axislegend(ax1, position = :rt, labelsize = 9, patchsize = (12, 12))

    ax2 = Axis(fig[1, 2];
        xlabel = "mode number j", ylabel = "1 / δ_j  (trusted = high)",
        yscale = log10,
        title = "Harmonic oscillator (real line),  N₁ = $N1, N₂ = $N2")
    plot_panel!(ax2, rep_B, "Harmonic oscillator,  N₁ = $N1, N₂ = $N2", tol)
    axislegend(ax2, position = :rt, labelsize = 9, patchsize = (12, 12))

    save(joinpath(OUTPUT_DIR, "eigen_drift_diagnostic.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "eigen_drift_diagnostic.png"), fig, px_per_unit = 4)

    r = Dict{String,Any}(
        "N1" => N1, "N2" => N2, "tol" => tol,
        "laplacian"  => Dict(
            "n_trusted"     => rep_A.n_trusted,
            "lam1"          => rep_A.lam1,
            "lam2"          => rep_A.lam2,
            "delta_ordinal" => rep_A.delta_ordinal,
            "delta_nearest" => rep_A.delta_nearest,
            "trusted"       => Int.(rep_A.trusted),
        ),
        "oscillator" => Dict(
            "n_trusted"     => rep_B.n_trusted,
            "lam1"          => rep_B.lam1,
            "lam2"          => rep_B.lam2,
            "delta_ordinal" => rep_B.delta_ordinal,
            "delta_nearest" => rep_B.delta_nearest,
            "trusted"       => Int.(rep_B.trusted),
        ),
    )
    dump_path !== nothing && write_json(dump_path, r)
    return r, rep_A, rep_B
end

function parse_args()
    N1 = 32; N2 = 48; tol = 1e-3; dump_path = nothing; i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--N1" && i < length(ARGS); N1 = parse(Int, ARGS[i+1]); i += 2
        elseif ARGS[i] == "--N2" && i < length(ARGS); N2 = parse(Int, ARGS[i+1]); i += 2
        elseif ARGS[i] == "--tol" && i < length(ARGS); tol = parse(Float64, ARGS[i+1]); i += 2
        elseif ARGS[i] == "--dump" && i < length(ARGS); dump_path = ARGS[i+1]; i += 2
        else i += 1 end
    end
    return N1, N2, tol, dump_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    N1, N2, tol, dp = parse_args()
    _, rA, rB = make_figure(; N1 = N1, N2 = N2, tol = tol, dump_path = dp)
    @printf("[Etude 18.5]  drift diagnostic at N1 = %d, N2 = %d, tol = %g\n", N1, N2, tol)
    @printf("  Dirichlet Laplacian : trusted = %d of %d\n", rA.n_trusted, length(rA.lam1))
    @printf("  Harmonic oscillator  : trusted = %d of %d\n", rB.n_trusted, length(rB.lam1))
    println("  figure: ", joinpath(OUTPUT_DIR, "eigen_drift_diagnostic.pdf"))
    dp !== nothing && println("  dumped: ", dp)
end
