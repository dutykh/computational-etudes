# eigen_heinrichs_condition.jl
# Chapter 18, Etude 18.8: condition-number surgery.
# Fourth-order clamped eigenproblem u'''' = lambda u: compare naive D^4
# bordered pencil vs. Heinrichs (1-x^2)^2 T_j basis.
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
# Part of "Computational Etudes: A Spectral Approach"

using CairoMakie
using Colors
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "heinrichs_basis.jl"))
using .HeinrichsBasis: heinrichs_clamped_matrix, naive_clamped_operator

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

function exact_clamped_beta(n::Int)
    f(b) = cos(2b) * cosh(2b) - 1.0
    betas = Float64[]
    for j in 1:n
        b_asym = (2j + 1) * π / 4.0
        lo, hi = b_asym - 0.5, b_asym + 0.5
        flo, fhi = f(lo), f(hi)
        if flo * fhi < 0
            for _ in 1:200
                mid = 0.5 * (lo + hi); fmid = f(mid)
                if flo * fmid < 0; hi = mid; fhi = fmid
                else; lo = mid; flo = fmid end
                hi - lo < 1e-14 && break
            end
            push!(betas, 0.5 * (lo + hi))
        end
    end
    return betas
end

function cond_naive(N::Int)
    A, B = naive_clamped_operator(N)
    lam = eigvals(A, B)
    lam = lam[isfinite.(lam)]
    lam = sort(real.(lam[abs.(imag.(lam)) .< 1e-6]))
    lam = lam[lam .> 1.0]
    return cond(A), lam[1:min(4, length(lam))]
end

function cond_heinrichs(N::Int)
    A, M, _ = heinrichs_clamped_matrix(N)
    lam = eigvals(A, M)
    lam = lam[isfinite.(lam)]
    lam = sort(real.(lam[abs.(imag.(lam)) .< 1e-6]))
    lam = lam[lam .> 1.0]
    A_std = M \ A
    return cond(A_std), lam[1:min(4, length(lam))]
end

function make_figure(; Ns = (12, 16, 24, 32, 48, 64, 96), dump_path = nothing)
    betas = exact_clamped_beta(6)
    lam_exact = betas .^ 4

    kappa_naive = Float64[]; kappa_hein = Float64[]
    drift_naive = Float64[]; drift_hein = Float64[]
    for N in Ns
        kn, lam_n = cond_naive(N)
        kh, lam_h = cond_heinrichs(N)
        push!(kappa_naive, kn); push!(kappa_hein, kh)
        push!(drift_naive, length(lam_n) >= 1 ? abs(lam_n[1] - lam_exact[1]) : NaN)
        push!(drift_hein,  length(lam_h) >= 1 ? abs(lam_h[1] - lam_exact[1]) : NaN)
    end

    fig = Figure(size = (1100, 450))
    ax1 = Axis(fig[1, 1]; xlabel = "N", ylabel = "κ (condition number)",
        title = "conditioning of the fourth-derivative matrix",
        xscale = log10, yscale = log10)
    Ns_f = Float64.(collect(Ns))
    scatterlines!(ax1, Ns_f, kappa_naive; color = NAVY, markercolor = :white,
        strokecolor = NAVY, strokewidth = 1.1, markersize = 9, linewidth = 0.8,
        label = "naive (D⁴ bordered)")
    scatterlines!(ax1, Ns_f, kappa_hein; color = CORAL, markercolor = :white,
        strokecolor = CORAL, strokewidth = 1.1, marker = :rect, markersize = 8,
        linewidth = 0.8, label = "Heinrichs (1−x²)² Tⱼ")
    lines!(ax1, Ns_f, kappa_naive[1] .* (Ns_f ./ Ns_f[1]) .^ 8;
           color = NAVY, linestyle = :dash, linewidth = 0.6, label = "N⁸")
    lines!(ax1, Ns_f, kappa_hein[1] .* (Ns_f ./ Ns_f[1]) .^ 4;
           color = CORAL, linestyle = :dash, linewidth = 0.6, label = "N⁴")
    axislegend(ax1, position = :lt, labelsize = 9)

    ax2 = Axis(fig[1, 2]; xlabel = "N", ylabel = "error in first eigenvalue",
        title = "first-eigenvalue error vs. N",
        xscale = log10, yscale = log10)
    scatterlines!(ax2, Ns_f, max.(drift_naive, 1e-17); color = NAVY, markercolor = :white,
        strokecolor = NAVY, strokewidth = 1.1, markersize = 9, linewidth = 0.8,
        label = "naive: |λ₁ − exact|")
    scatterlines!(ax2, Ns_f, max.(drift_hein, 1e-17); color = CORAL, markercolor = :white,
        strokecolor = CORAL, strokewidth = 1.1, marker = :rect, markersize = 8,
        linewidth = 0.8, label = "Heinrichs: |λ₁ − exact|")
    axislegend(ax2, position = :rt, labelsize = 9)

    save(joinpath(OUTPUT_DIR, "eigen_heinrichs_condition.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "eigen_heinrichs_condition.png"), fig, px_per_unit = 4)

    r = Dict("Ns" => Ns_f, "kappa_naive" => kappa_naive,
        "kappa_heinrichs" => kappa_hein,
        "drift_naive" => drift_naive, "drift_heinrichs" => drift_hein,
        "lam_exact_first4" => lam_exact[1:4])
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
    println("[Etude 18.8]  Heinrichs condition-number surgery")
    @printf("  exact lambda_1..4 = %s\n", string(round.(r["lam_exact_first4"], digits=3)))
    println("  N, kappa(naive), kappa(Heinrichs):")
    for (N, kn, kh) in zip(r["Ns"], r["kappa_naive"], r["kappa_heinrichs"])
        @printf("    N=%3d   kappa_naive=%.2e   kappa_Heinrichs=%.2e\n", Int(N), kn, kh)
    end
    println("  figure: ", joinpath(OUTPUT_DIR, "eigen_heinrichs_condition.pdf"))
    dp !== nothing && println("  dumped: ", dp)
end
