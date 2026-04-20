# eigen_parameter_map.jl
# Chapter 18, Etude 18.10: map the neutral curve safely.
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
# Part of "Computational Etudes: A Spectral Approach"

using CairoMakie
using Colors
using LinearAlgebra
using Random
using Printf

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

_json_val(x::Bool) = x ? "true" : "false"
_json_val(x::Integer) = string(x)
_json_val(x::Real) = isfinite(x) ? string(Float64(x)) : (isnan(x) ? "NaN" : (x > 0 ? "Infinity" : "-Infinity"))
_json_val(x::AbstractString) = "\"" * replace(String(x), "\\" => "\\\\", "\"" => "\\\"") * "\""
_json_val(v::AbstractVector) = "[" * join(_json_val.(v), ", ") * "]"
_json_val(v::AbstractMatrix) = _json_val([collect(v[i, :]) for i in 1:size(v, 1)])
_json_val(d::AbstractDict) = "{" * join([_json_val(string(k)) * ": " * _json_val(v) for (k, v) in d], ", ") * "}"
write_json(path, x) = open(io -> print(io, _json_val(x)), path, "w")

set_theme!(Theme(fontsize = 11,
    fonts = (regular = "CMU Serif", bold = "CMU Serif Bold", italic = "CMU Serif Italic"),
    Axis = (xlabelsize = 11, ylabelsize = 11, titlesize = 11,
            xticklabelsize = 10, yticklabelsize = 10, spinewidth = 0.8,
            xtickwidth = 0.8, ytickwidth = 0.8), Legend = (labelsize = 9,)))

const NAVY = colorant"#142D6E"
const CORAL = colorant"#E74C3C"
const OUTPUT_DIR = joinpath(@__DIR__, "..", "..", "..", "textbook", "figures", "ch18", "julia")
mkpath(OUTPUT_DIR)

function build_op(N::Int, α::Real, β::Real)
    D, x = cheb_matrix(N)
    D2 = D * D
    return -D2[2:N, 2:N] + Diagonal(α .+ β .* x[2:N] .^ 2)
end

function qr_lowest_two(N, α, β)
    A = build_op(N, α, β)
    lam = sort(real.(eigvals(A)))
    return lam[1], lam[2]
end

function inv_iter_shift(N, α, β, shift; max_iter = 20)
    A = build_op(N, α, β)
    n = size(A, 1)
    Random.seed!(0)
    v = randn(n); v ./= norm(v)
    F = lu(A - shift * I)
    lam = 0.0
    for _ in 1:max_iter
        w = F \ v
        v = w ./ norm(w)
        lam = v ⋅ (A * v)
    end
    return lam
end

function make_figure(; N = 40, dump_path = nothing)
    alpha_coarse = collect(range(-2.0, 2.0; length = 9))
    beta_coarse = collect(range(0.0, 8.0; length = 9))
    lambda1 = zeros(length(alpha_coarse), length(beta_coarse))
    for (i, a) in enumerate(alpha_coarse), (j, b) in enumerate(beta_coarse)
        lambda1[i, j], _ = qr_lowest_two(N, a, b)
    end

    beta_fine = collect(range(0.0, 8.0; length = 40))
    lam1_safe = [qr_lowest_two(N, 0.0, b)[1] for b in beta_fine]
    l1_0, l2_0 = qr_lowest_two(N, 0.0, beta_fine[1])
    shift = l2_0
    lam1_naive = Float64[l2_0]
    for k in 2:length(beta_fine)
        lam = inv_iter_shift(N, 0.0, beta_fine[k], shift)
        push!(lam1_naive, lam); shift = lam
    end

    fig = Figure(size = (1150, 450))
    ax1 = Axis(fig[1, 1]; xlabel = "α", ylabel = "β",
               title = "coarse QR grid (9×9): lowest eigenvalue")
    cf = contourf!(ax1, alpha_coarse, beta_coarse, lambda1; levels = 18, colormap = :viridis)
    contour!(ax1, alpha_coarse, beta_coarse, lambda1; levels = 8, color = :white, linewidth = 0.5)
    scatter!(ax1, [(a, b) for a in alpha_coarse, b in beta_coarse] |> vec;
        color = :white, strokecolor = :black, strokewidth = 0.6, markersize = 6)
    Colorbar(fig[1, 1, Right()], cf; label = "λ₁")

    ax2 = Axis(fig[1, 2]; xlabel = "β (with α = 0)", ylabel = "eigenvalue",
               title = "eigenvalue tracking: safe vs. naive local iteration")
    scatterlines!(ax2, beta_fine, lam1_safe; color = NAVY, markercolor = :white,
        strokecolor = NAVY, strokewidth = 1.0, markersize = 6, linewidth = 0.8,
        label = "SAFE: QR re-anchor each step")
    scatterlines!(ax2, beta_fine, lam1_naive; color = CORAL, marker = :xcross,
        linewidth = 0.8, markersize = 7,
        label = "NAIVE: inverse iter seeded from wrong branch")
    axislegend(ax2, position = :rb, labelsize = 9)

    save(joinpath(OUTPUT_DIR, "eigen_parameter_map.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "eigen_parameter_map.png"), fig, px_per_unit = 4)

    r = Dict("N" => N, "alpha_grid" => alpha_coarse, "beta_grid" => beta_coarse,
        "lambda1_coarse" => lambda1, "beta_fine" => beta_fine,
        "lam1_fine_safe" => lam1_safe, "lam1_fine_naive" => lam1_naive)
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
    println("[Etude 18.10]  parameter map")
    @printf("  coarse grid: %d x %d = %d QR solves\n",
            length(r["alpha_grid"]), length(r["beta_grid"]),
            length(r["alpha_grid"]) * length(r["beta_grid"]))
    @printf("  safe tracker: starts at %.4f, ends at %.4f\n",
            r["lam1_fine_safe"][1], r["lam1_fine_safe"][end])
    @printf("  naive tracker: starts at %.4f, ends at %.4f\n",
            r["lam1_fine_naive"][1], r["lam1_fine_naive"][end])
    println("  figure: ", joinpath(OUTPUT_DIR, "eigen_parameter_map.pdf"))
    dp !== nothing && println("  dumped: ", dp)
end
