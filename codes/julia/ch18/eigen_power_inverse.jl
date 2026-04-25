# eigen_power_inverse.jl
# Chapter 18, Etude 18.9: one mode at a time (power and inverse iteration).
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

build_matrix(N::Int) = (D = cheb_matrix(N)[1]; -(D * D)[2:N, 2:N])

function power_method(A; max_iter = 80, seed = 0)
    Random.seed!(seed)
    n = size(A, 1)
    v = randn(n); v ./= norm(v)
    hist = Float64[]
    for _ in 1:max_iter
        w = A * v
        push!(hist, v ⋅ w)
        v = w ./ norm(w)
    end
    return hist, v
end

function inverse_iteration(A, shift; max_iter = 30, seed = 1)
    Random.seed!(seed)
    n = size(A, 1)
    v = randn(n); v ./= norm(v)
    hist = Float64[]
    F = lu(A - shift * I)
    for _ in 1:max_iter
        w = F \ v
        push!(hist, v ⋅ (A * v))
        v = w ./ norm(w)
    end
    push!(hist, v ⋅ (A * v))
    return hist, v
end

function make_figure(; N::Int = 32, dump_path = nothing)
    A = build_matrix(N)
    lam_ref = sort(real.(eigvals(A)))
    hist_pow, _ = power_method(A)
    lam_max = lam_ref[end]

    shifts = [5.0, 90.0, 250.0]
    histories = []; targets = Float64[]
    for mu in shifts
        h, _ = inverse_iteration(A, mu)
        final = h[end]
        t = lam_ref[argmin(abs.(lam_ref .- final))]
        push!(histories, h); push!(targets, t)
    end

    lam_a, lam_b = lam_ref[4], lam_ref[5]
    mu_bad = 0.5 * (lam_a + lam_b)
    hist_bad, _ = inverse_iteration(A, mu_bad)

    rate_gap = lam_ref[end-1] / lam_ref[end]
    fig = Figure(size = (1500, 480))
    ax1 = Axis(fig[1, 1]; xlabel = "iteration k", ylabel = "|μ_k − λ_max|",
        yscale = log10,
        title = "power method → λ_max = $(round(lam_max, digits=0))" *
                "  (stalls: |λ_{N-1}/λ_N| = $(round(rate_gap, digits=3)))")
    scatterlines!(ax1, 1:length(hist_pow), max.(abs.(hist_pow .- lam_max), 1e-18);
        color = NAVY, markercolor = :white, strokecolor = NAVY,
        strokewidth = 0.9, markersize = 6, linewidth = 0.8)

    ax2 = Axis(fig[1, 2]; xlabel = "iteration k", ylabel = "|μ_k − λ_target|",
        yscale = log10, title = "inverse iteration: three shifts, three modes")
    colors_B = [NAVY, CORAL, TEAL]
    for (mu, h, t, c) in zip(shifts, histories, targets, colors_B)
        err = max.(abs.(h .- t), 1e-18)
        scatterlines!(ax2, 1:length(h), err; color = c, markercolor = :white,
            strokecolor = c, strokewidth = 0.9, markersize = 6, linewidth = 0.8,
            label = "shift μ = $mu  →  λ = $(round(t, digits=3))")
    end
    axislegend(ax2, position = :rt, labelsize = 8)

    ax3 = Axis(fig[1, 3]; xlabel = "iteration k", ylabel = "|μ_k − λ|",
        yscale = log10, title = "cautionary: μ = $(round(mu_bad, digits=3)) between two modes")
    err_a = max.(abs.(hist_bad .- lam_a), 1e-18)
    err_b = max.(abs.(hist_bad .- lam_b), 1e-18)
    scatterlines!(ax3, 1:length(hist_bad), err_a; color = NAVY, markercolor = :white,
        strokecolor = NAVY, strokewidth = 0.9, markersize = 6, linewidth = 0.8,
        label = "dist to λ = $(round(lam_a, digits=3))")
    scatterlines!(ax3, 1:length(hist_bad), err_b; color = CORAL, markercolor = :white,
        strokecolor = CORAL, strokewidth = 0.9, markersize = 6, marker = :rect,
        linewidth = 0.8, label = "dist to λ = $(round(lam_b, digits=3))")
    axislegend(ax3, position = :rt, labelsize = 8)

    save(joinpath(OUTPUT_DIR, "eigen_power_inverse.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "eigen_power_inverse.png"), fig, px_per_unit = 4)

    r = Dict("N" => N, "lam_max" => lam_max, "final_power_method" => hist_pow[end],
        "shifts" => shifts, "targets" => targets, "bad_shift" => mu_bad,
        "bad_pair" => [lam_a, lam_b], "final_bad" => hist_bad[end])
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
    @printf("[Etude 18.9]  power and inverse iteration, N = %d\n", r["N"])
    @printf("  power method final -> %.6f (exact lam_max = %.6f)\n",
            r["final_power_method"], r["lam_max"])
    for (mu, t) in zip(r["shifts"], r["targets"])
        @printf("  inverse iter, shift %g  ->  lambda = %.6f\n", mu, t)
    end
    println("  figure: ", joinpath(OUTPUT_DIR, "eigen_power_inverse.pdf"))
    dp !== nothing && println("  dumped: ", dp)
end
