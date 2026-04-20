# eigen_benchmark_finite.jl
#
# Chapter 18: Linear Spectral Eigenproblems
# Computational Etude 18.3: Finite-interval benchmark and the N/2 rule.
#
# Repeats the Dirichlet Laplacian at N = 16, 32, 64 and plots the
# absolute eigenvalue error on a semilog scale vs. mode number j.
# Reproduces Boyd (2000) Figs. 7.1 and 7.3.
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
# Part of "Computational Etudes: A Spectral Approach"
# https://github.com/dutykh/computational-etudes

using CairoMakie
using Colors
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

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
const TEAL  = colorant"#16A085"

const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch18", "julia")
mkpath(OUTPUT_DIR)

function solve_dirichlet(N::Int)
    D, _ = cheb_matrix(N)
    D2 = D * D
    A = -D2[2:N, 2:N]
    lam = real.(eigvals(A))
    return sort(lam)
end

function count_good(err::Vector{<:Real}, tol::Real = 1e-2)
    for (k, e) in enumerate(err)
        if e >= tol
            return k - 1
        end
    end
    return length(err)
end

function make_figure(; Ns = (16, 32, 64), dump_path = nothing)
    colours = [NAVY, CORAL, TEAL]
    markers = [:circle, :rect, :utriangle]
    counts = Dict{Int,Int}()
    spectra = Dict{Int,Vector{Float64}}()

    fig = Figure(size = (760, 460))
    ax  = Axis(fig[1, 1];
        xlabel = "mode number j",
        ylabel = "|λ_num − λ_exact|",
        title  = "Dirichlet Laplacian: errors at three resolutions",
        yscale = log10,
        limits = ((0, maximum(Ns) + 2), (1e-17, 1e6)))

    for (N, c, m) in zip(Ns, colours, markers)
        lam = solve_dirichlet(N)
        j = collect(1:length(lam))
        lam_exact = (j .* π ./ 2) .^ 2
        err = abs.(lam .- lam_exact)
        spectra[N] = lam
        counts[N] = count_good(err, 1e-2)
        err_plot = max.(err, 1e-17)
        scatterlines!(ax, j, err_plot;
                      color = c, markercolor = :white,
                      strokecolor = c, strokewidth = 1.0,
                      linewidth = 0.6, marker = m, markersize = 8,
                      label = "N = $N")
        vlines!(ax, [N / 2]; color = c, linestyle = :dot,
                linewidth = 0.7, alpha = 0.5)
    end

    hlines!(ax, [1e-2]; color = :black, linestyle = :dash,
            linewidth = 0.6, alpha = 0.5)
    text!(ax, 1.5, 2e-2; text = "0.01 tolerance",
          fontsize = 8, color = (:black, 0.7))
    axislegend(ax, position = :rb, labelsize = 9, patchsize = (12, 12))

    save(joinpath(OUTPUT_DIR, "eigen_benchmark_finite.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "eigen_benchmark_finite.png"), fig, px_per_unit = 4)

    r = Dict(
        "Ns"           => collect(Ns),
        "spectra"      => Dict(string(k) => v for (k, v) in spectra),
        "good_counts"  => Dict(string(k) => v for (k, v) in counts),
    )
    dump_path !== nothing && write_json(dump_path, r)
    return r, counts, spectra
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
    _, counts, _ = make_figure(; dump_path = dp)
    println("[Etude 18.3]")
    for N in sort(collect(keys(counts)))
        @printf("  N = %2d:  %d eigenvalues with |err| < 1e-2   (heuristic N/2 = %d)\n",
                N, counts[N], N ÷ 2)
    end
    println("  figure: ", joinpath(OUTPUT_DIR, "eigen_benchmark_finite.pdf"))
    dp !== nothing && println("  dumped: ", dp)
end
