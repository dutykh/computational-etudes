# eigen_two_formulations.jl
#
# Chapter 18: Linear Spectral Eigenproblems
# Computational Etude 18.2: Assemble the pencil.
#
# Same problem as Etude 18.1 --- u_xx + lambda u = 0, u(+-1) = 0 ---
# solved two ways:
#
#   (A) Boundary bordering: full (N+1) x (N+1) pencil A v = lambda B v
#       with rows 1 and N+1 replaced by identity rows; B has zeros on
#       those rows. Two infinite eigenvalues are filtered out.
#
#   (B) Basis recombination: interior-block regular problem
#       -D^2_int u = lambda u.
#
# Both recover the same finite spectrum.
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

const NAVY   = colorant"#142D6E"
const CORAL  = colorant"#E74C3C"
const TEAL   = colorant"#16A085"
const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch18", "julia")
mkpath(OUTPUT_DIR)

function boundary_bordering(N::Int)
    D, x = cheb_matrix(N)
    D2 = D * D
    A = -D2
    B = Matrix{Float64}(I, N + 1, N + 1)
    A[1, :] .= 0.0; A[1, 1] = 1.0; B[1, :] .= 0.0
    A[N+1, :] .= 0.0; A[N+1, N+1] = 1.0; B[N+1, :] .= 0.0
    return A, B, x
end

function basis_recombination(N::Int)
    D, x = cheb_matrix(N)
    D2 = D * D
    A_int = -D2[2:N, 2:N]
    return A_int, x
end

function solve_pencil(A, B)
    F = eigen(A, B)
    lam = F.values
    keep = isfinite.(lam) .& (abs.(imag.(lam)) .< 1e-6)
    return sort(real.(lam[keep]))
end

solve_regular(A) = sort(real.(filter(z -> abs(imag(z)) < 1e-6, eigen(A).values)))

function make_figure(; N::Int = 16, dump_path = nothing)
    A_pen, B_pen, x = boundary_bordering(N)
    A_reg, _ = basis_recombination(N)

    lam_pen = solve_pencil(A_pen, B_pen)
    lam_reg = solve_regular(A_reg)
    m = min(length(lam_pen), length(lam_reg))
    lam_pen = lam_pen[1:m]; lam_reg = lam_reg[1:m]

    j = collect(1:m)
    lam_exact = (j .* π ./ 2) .^ 2
    err_pen = abs.(lam_pen .- lam_exact)
    err_reg = abs.(lam_reg .- lam_exact)
    diff    = abs.(lam_pen .- lam_reg)

    fig = Figure(size = (1050, 420))

    ax1 = Axis(fig[1, 1];
        xlabel = "mode number j", ylabel = "|λ_num − λ_exact|",
        title = @sprintf("both formulations,  N = %d", N),
        yscale = log10, limits = (nothing, (1e-17, 1e5)))
    scatter!(ax1, j, max.(err_pen, 1e-17); color = :white, strokecolor = NAVY,
             strokewidth = 1.2, markersize = 9, label = "pencil")
    scatter!(ax1, j, max.(err_reg, 1e-17); color = :white, strokecolor = CORAL,
             strokewidth = 1.2, marker = :rect, markersize = 8, label = "regular")
    scatter!(ax1, j, max.(diff, 1e-17); color = TEAL, marker = :xcross,
             markersize = 10, label = "|pencil − regular|")
    hlines!(ax1, [1e-2]; color = :black, linestyle = :dash, linewidth = 0.6)
    axislegend(ax1, position = :rb, labelsize = 8, patchsize = (12, 12))

    ax2 = Axis(fig[1, 2];
        xlabel = "column index", ylabel = "row index",
        title = "log₁₀|A| for boundary-bordered pencil",
        yreversed = true, limits = ((0.5, N + 1.5), (0.5, N + 1.5)))
    logA = log10.(max.(abs.(A_pen), 1e-20))
    heatmap!(ax2, 1:N+1, 1:N+1, logA'; colormap = :Blues, interpolate = false)

    save(joinpath(OUTPUT_DIR, "eigen_two_formulations.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "eigen_two_formulations.png"), fig, px_per_unit = 4)

    r = Dict(
        "N"                  => N,
        "lam_pencil"         => lam_pen,
        "lam_regular"        => lam_reg,
        "lam_exact"          => lam_exact,
        "max_abs_agreement"  => maximum(diff),
        "max_err_pencil"     => maximum(err_pen),
        "max_err_regular"    => maximum(err_reg),
    )
    dump_path !== nothing && write_json(dump_path, r)
    return r
end

function parse_args()
    N = 16; dump_path = nothing; i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--N" && i < length(ARGS); N = parse(Int, ARGS[i+1]); i += 2
        elseif ARGS[i] == "--dump" && i < length(ARGS); dump_path = ARGS[i+1]; i += 2
        else i += 1 end
    end
    return N, dump_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    N, dp = parse_args()
    r = make_figure(; N = N, dump_path = dp)
    @printf("[Etude 18.2]  N = %d\n", N)
    @printf("  max |lam_pencil - lam_regular| = %.3e\n", r["max_abs_agreement"])
    @printf("  max |lam_pencil  - exact|      = %.3e\n", r["max_err_pencil"])
    @printf("  max |lam_regular - exact|      = %.3e\n", r["max_err_regular"])
    println("  figure: ", joinpath(OUTPUT_DIR, "eigen_two_formulations.pdf"))
    dp !== nothing && println("  dumped: ", dp)
end
