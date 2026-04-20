# eigen_benchmark_oscillator.jl
#
# Chapter 18: Linear Spectral Eigenproblems
# Computational Etude 18.4: the infinite-interval tax.
#
# Solves u_xx + (lambda - x^2) u = 0, u -> 0 as |x| -> inf, using the
# rational Chebyshev algebraic map. Exact spectrum lambda_j = 2j+1,
# eigenfunctions H_j(x) exp(-x^2/2).
#
# Reproduces Boyd (2000) Figs 7.4 and 7.6.
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
# Part of "Computational Etudes: A Spectral Approach"
# https://github.com/dutykh/computational-etudes

using CairoMakie
using Colors
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "rational_chebyshev.jl"))
using .RationalChebyshev: rational_chebyshev_derivative_matrices

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

const NAVY   = colorant"#142D6E"
const SKY    = colorant"#7896D2"
const CORAL  = colorant"#E74C3C"
const PURPLE = colorant"#8E44AD"

const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch18", "julia")
mkpath(OUTPUT_DIR)

function solve_oscillator(N::Int, L::Real)
    _, D2_x, x = rational_chebyshev_derivative_matrices(N, L)
    H = -D2_x + Diagonal(x .^ 2)
    lam = eigvals(H)
    keep = abs.(imag.(lam)) .< 1e-6
    lam = real.(lam[keep])
    lam = lam[lam .> 0]
    return sort(lam)
end

function count_good(err, tol = 1e-2)
    for (k, e) in enumerate(err)
        e >= tol && return k - 1
    end
    return length(err)
end

function make_figure(; Ns = (16, 32), L_scan = (2.0, 4.0, 8.0), L_best::Real = 4.0,
                     dump_path = nothing)
    fig = Figure(size = (1100, 450))

    # Panel A: N-scan at fixed L_best
    ax1 = Axis(fig[1, 1];
        xlabel = "mode number j",
        ylabel = "|λ_num − λ_exact|",
        title  = "infinite-interval tax: oscillator spectrum",
        yscale = log10, limits = (nothing, (1e-17, 1e4)))
    markers = [:circle, :rect]
    cols = [NAVY, CORAL]
    counts_by_N = Dict{Int,Int}()
    for (N, m, c) in zip(Ns, markers, cols)
        lam = solve_oscillator(N, L_best)
        j = collect(0:length(lam)-1)
        lam_exact = 2 .* j .+ 1
        err = abs.(lam .- lam_exact)
        err_plot = max.(err, 1e-17)
        scatterlines!(ax1, j .+ 1, err_plot;
                      color = c, markercolor = :white,
                      strokecolor = c, strokewidth = 1.0,
                      linewidth = 0.6, marker = m, markersize = 8,
                      label = "N = $N, L = $L_best")
        counts_by_N[N] = count_good(err, 1e-2)
    end
    hlines!(ax1, [1e-2]; color = :black, linestyle = :dash,
            linewidth = 0.6, alpha = 0.5)
    axislegend(ax1, position = :rb, labelsize = 9, patchsize = (12, 12))

    # Panel B: L-scan at fixed N
    ax2 = Axis(fig[1, 2];
        xlabel = "mode number j",
        ylabel = "|λ_num − λ_exact|",
        title  = "L-scan at N = 32: map-parameter sensitivity",
        yscale = log10, limits = (nothing, (1e-17, 1e4)))
    markers2 = [:utriangle, :circle, :dtriangle]
    cols2 = [SKY, NAVY, PURPLE]
    counts_by_L = Dict{Float64,Int}()
    N_scan = 32
    for (L, m, c) in zip(L_scan, markers2, cols2)
        lam = solve_oscillator(N_scan, L)
        j = collect(0:length(lam)-1)
        lam_exact = 2 .* j .+ 1
        err = abs.(lam .- lam_exact)
        err_plot = max.(err, 1e-17)
        scatterlines!(ax2, j .+ 1, err_plot;
                      color = c, markercolor = :white,
                      strokecolor = c, strokewidth = 1.0,
                      linewidth = 0.6, marker = m, markersize = 8,
                      label = "L = $L")
        counts_by_L[L] = count_good(err, 1e-2)
    end
    hlines!(ax2, [1e-2]; color = :black, linestyle = :dash,
            linewidth = 0.6, alpha = 0.5)
    axislegend(ax2, position = :rb, labelsize = 9, patchsize = (12, 12))

    save(joinpath(OUTPUT_DIR, "eigen_benchmark_oscillator.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "eigen_benchmark_oscillator.png"), fig, px_per_unit = 4)

    spectra = Dict{String,Vector{Float64}}()
    for N in Ns
        spectra["N$(N)_L$(L_best)"] = solve_oscillator(N, L_best)
    end
    for L in L_scan
        spectra["N$(N_scan)_L$(L)"] = solve_oscillator(N_scan, L)
    end

    r = Dict{String,Any}(
        "Ns"            => collect(Ns),
        "L_best"        => L_best,
        "L_scan"        => collect(L_scan),
        "counts_by_N"   => Dict(string(k) => v for (k, v) in counts_by_N),
        "counts_by_L"   => Dict(string(k) => v for (k, v) in counts_by_L),
        "spectra"       => spectra,
    )
    dump_path !== nothing && write_json(dump_path, r)
    return r, counts_by_N, counts_by_L
end

function parse_args()
    dump_path = nothing; i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--dump" && i < length(ARGS); dump_path = ARGS[i+1]; i += 2
        else i += 1 end
    end
    return dump_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    dp = parse_args()
    _, cN, cL = make_figure(; dump_path = dp)
    println("[Etude 18.4]  harmonic oscillator on (-inf, +inf)")
    println("  L = 4.0 scan over N:")
    for N in sort(collect(keys(cN)))
        @printf("    N = %2d:  %d good eigenvalues (|err| < 1e-2)\n", N, cN[N])
    end
    println("  N = 32 scan over L:")
    for L in sort(collect(keys(cL)))
        @printf("    L = %.1f:  %d good eigenvalues\n", L, cL[L])
    end
    println("  figure: ", joinpath(OUTPUT_DIR, "eigen_benchmark_oscillator.pdf"))
    dp !== nothing && println("  dumped: ", dp)
end
