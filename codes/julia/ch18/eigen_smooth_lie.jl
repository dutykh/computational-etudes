# eigen_smooth_lie.jl
#
# Chapter 18: Linear Spectral Eigenproblems
# Computational Etude 18.1: The smooth lie.
#
# Solves the benchmark eigenvalue problem
#
#     u_xx + lambda u = 0,   u(+-1) = 0
#
# using Chebyshev pseudospectral collocation at N = 16. Exact eigenvalues
# are lambda_j = (j pi / 2)^2 with eigenfunctions cos(j pi x / 2) (j odd)
# or sin(j pi x / 2) (j even).
#
# This reproduces Boyd (2000) Fig 7.2: low modes recovered to near machine
# precision, high modes SMOOTH but quantitatively wrong.
#
# Output figure: textbook/figures/ch18/julia/eigen_smooth_lie.{pdf,png}
#
# Author: Dr. Denys Dutykh
#         Mathematics Department
#         Khalifa University of Science and Technology
#         Abu Dhabi, UAE
#
# Part of the book "Computational Etudes: A Spectral Approach"
# https://github.com/dutykh/computational-etudes

using CairoMakie
using Colors
using LinearAlgebra
using Printf

# Minimal inline JSON writer for the cross-language validation harness.
# The dump format uses only Dict{String, Union{Number, Vector{Number}}}
# so a full JSON library is unnecessary and would be a heavy dependency
# for scripts that only write a handful of bytes.
_json_val(x::Bool)    = x ? "true" : "false"
_json_val(x::Integer) = string(x)
_json_val(x::Real)    = isfinite(x) ? string(Float64(x)) : (isnan(x) ? "NaN" : (x > 0 ? "Infinity" : "-Infinity"))
_json_val(x::AbstractString) = "\"" * replace(String(x), "\\" => "\\\\", "\"" => "\\\"") * "\""
_json_val(v::AbstractVector) = "[" * join(_json_val.(v), ", ") * "]"
function _json_val(d::AbstractDict)
    pairs = String[]
    for (k, v) in d
        push!(pairs, _json_val(string(k)) * ": " * _json_val(v))
    end
    return "{" * join(pairs, ", ") * "}"
end
write_json(path::AbstractString, x) = open(io -> print(io, _json_val(x)), path, "w")

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

# publication-quality theme
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
const TEAL   = colorant"#16A085"
const PURPLE = colorant"#8E44AD"
const ORANGE = colorant"#E67E22"

const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch18", "julia")
mkpath(OUTPUT_DIR)

exact_eigenvalue(j) = (j * π / 2)^2
exact_eigenfunction(x, j) = isodd(j) ? cos(j * π * x / 2) : sin(j * π * x / 2)

"""
    solve_problem(N)

Assemble -D^2_int and return (sorted eigenvalues, sorted eigenvectors,
grid x).
"""
function solve_problem(N::Int)
    D, x = cheb_matrix(N)
    D2 = D * D
    A = -D2[2:N, 2:N]
    F = eigen(A)
    order = sortperm(real.(F.values))
    lam = real.(F.values[order])
    V = real.(F.vectors[:, order])
    # sign-fix: first interior component positive
    for k in 1:size(V, 2)
        if V[1, k] < 0
            V[:, k] .*= -1
        end
    end
    return lam, V, x
end

pad_boundary(v) = vcat(0.0, v, 0.0)

function normalise_max(u)
    m = maximum(abs, u)
    m == 0 && return u
    u = u ./ m
    idx = argmax(abs.(u))
    u[idx] < 0 && (u = -u)
    return u
end

function make_figure(; N::Int = 16, low_mode::Int = 3, high_mode::Int = 15,
                     dump_path::Union{Nothing,String} = nothing)
    lam, V, x = solve_problem(N)
    # indices into sorted spectrum
    idx_low, idx_high = low_mode, high_mode

    u_low_num  = normalise_max(pad_boundary(V[:, idx_low]))
    u_high_num = normalise_max(pad_boundary(V[:, idx_high]))

    xfine = range(-1.0, 1.0; length = 801)
    u_low_exact  = normalise_max(exact_eigenfunction.(xfine, low_mode))
    u_high_exact = normalise_max(exact_eigenfunction.(xfine, high_mode))

    lam_num_low   = lam[idx_low]
    lam_num_high  = lam[idx_high]
    lam_ex_low    = exact_eigenvalue(low_mode)
    lam_ex_high   = exact_eigenvalue(high_mode)

    fig = Figure(size = (1050, 420))

    ax1 = Axis(fig[1, 1];
        xlabel = "x", ylabel = "u(x)",
        title = @sprintf("mode j=%d:  λ_num=%.4f,  λ_exact=%.4f",
                         low_mode, lam_num_low, lam_ex_low),
        limits = ((-1.05, 1.05), (-1.15, 1.15)),
    )
    lines!(ax1, xfine, u_low_exact, color = NAVY, linewidth = 1.3,
           label = "exact mode $low_mode")
    scatter!(ax1, x, u_low_num, color = :white, strokecolor = CORAL,
             strokewidth = 1.2, markersize = 9, label = "numerical (N=$N)")
    axislegend(ax1, position = :cb, framevisible = true, patchsize = (12, 12))

    ax2 = Axis(fig[1, 2];
        xlabel = "x", ylabel = "u(x)",
        title = @sprintf("mode j=%d:  λ_num=%.4f,  λ_exact=%.4f",
                         high_mode, lam_num_high, lam_ex_high),
        limits = ((-1.05, 1.05), (-1.15, 1.15)),
    )
    lines!(ax2, xfine, u_high_exact, color = NAVY, linewidth = 1.3,
           label = "exact mode $high_mode")
    scatterlines!(ax2, x, u_high_num, color = CORAL, linewidth = 0.8,
                  markercolor = :white, strokecolor = CORAL,
                  strokewidth = 1.2, markersize = 9,
                  label = "numerical (N=$N)")
    axislegend(ax2, position = :cb, framevisible = true, patchsize = (12, 12))

    save(joinpath(OUTPUT_DIR, "eigen_smooth_lie.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "eigen_smooth_lie.png"), fig, px_per_unit = 4)

    results = Dict(
        "N"              => N,
        "low_mode"       => low_mode,
        "high_mode"      => high_mode,
        "eigvals_sorted" => lam,
        "lam_num_low"    => lam_num_low,
        "lam_num_high"   => lam_num_high,
        "lam_exact_low"  => lam_ex_low,
        "lam_exact_high" => lam_ex_high,
        "abs_err_low"    => abs(lam_num_low - lam_ex_low),
        "abs_err_high"   => abs(lam_num_high - lam_ex_high),
    )

    if dump_path !== nothing
        write_json(dump_path, results)
    end
    return results
end

function parse_args()
    N = 16
    dump_path = nothing
    i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--N" && i < length(ARGS)
            N = parse(Int, ARGS[i + 1]); i += 2
        elseif ARGS[i] == "--dump" && i < length(ARGS)
            dump_path = ARGS[i + 1]; i += 2
        else
            i += 1
        end
    end
    return N, dump_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    N, dump_path = parse_args()
    r = make_figure(; N = N, dump_path = dump_path)
    @printf("[Etude 18.1]  N = %d\n", N)
    @printf("  mode %2d:  lambda_num = %.6e   exact = %.6e   |err| = %.2e\n",
            r["low_mode"], r["lam_num_low"], r["lam_exact_low"], r["abs_err_low"])
    @printf("  mode %2d:  lambda_num = %.6e   exact = %.6e   |err| = %.2e\n",
            r["high_mode"], r["lam_num_high"], r["lam_exact_high"], r["abs_err_high"])
    println("  figure: ", joinpath(OUTPUT_DIR, "eigen_smooth_lie.pdf"))
    dump_path !== nothing && println("  dumped: ", dump_path)
end
