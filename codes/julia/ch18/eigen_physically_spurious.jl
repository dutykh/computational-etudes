# eigen_physically_spurious.jl
# Chapter 18, Etude 18.7: manufacturing fake instability.
# Gottlieb-Orszag streamfunction: nu u_xxxx = lambda u_xx,
# u(+-1) = u_x(+-1) = 0. Naive formulation produces spurious large
# positive eigenvalues (Dawkins et al. 1998); cured formulation does not.
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

function naive_go(N::Int, nu::Real = 1.0)
    # Tau-style bordering: the four boundary conditions replace the
    # LAST four rows of the pencil.  This is the standard tau placement
    # and produces, with proper alpha/beta handling, four algebraically
    # infinite eigenvalues (one per BC row) and exactly ONE finite
    # spurious positive eigenvalue whose magnitude approaches O(N^4)
    # asymptotically (Dawkins-Dunbar-Douglass 1998).
    D, _ = cheb_matrix(N)
    D2 = D * D
    D4 = D2 * D2
    A = nu * copy(D4); B = copy(D2)
    ID = Matrix{Float64}(I, N + 1, N + 1)
    A[end-3, :] .= ID[1, :];   B[end-3, :] .= 0    # u(+1)  = 0
    A[end-2, :] .= ID[N+1, :]; B[end-2, :] .= 0    # u(-1)  = 0
    A[end-1, :] .= D[1, :];    B[end-1, :] .= 0    # u'(+1) = 0
    A[end,   :] .= D[N+1, :];  B[end,   :] .= 0    # u'(-1) = 0
    return A, B
end

function cured_go(N::Int, nu::Real = 1.0)
    D, _ = cheb_matrix(N)
    D2 = D * D
    D2i = D2[2:N, 2:N]
    A = nu * (D2i * D2i)
    B = copy(D2i)
    A[1, :]   .= D[1,   2:N]; B[1, :]   .= 0
    A[end, :] .= D[N+1, 2:N]; B[end, :] .= 0
    return A, B
end

function finite_real_eigvals(A, B; beta_tol = 1e-10)
    # Genuinely finite real eigenvalues of the generalised pencil
    # (A, B), via the homogeneous (alpha, beta) decomposition. An
    # eigenvalue is treated as algebraically infinite iff
    # |beta| < beta_tol * max(|alpha|, |beta|).  Forcing complex inputs
    # gives strictly 1x1 diagonal blocks in the generalised Schur
    # decomposition, so alpha and beta are unambiguous per eigenvalue.
    F = schur(complex(A), complex(B))
    alpha = F.alpha
    beta  = F.beta
    mag = max.(abs.(alpha), abs.(beta))
    finite = abs.(beta) .> beta_tol .* mag
    lam = alpha[finite] ./ beta[finite]
    return sort(real.(lam[abs.(imag.(lam)) .< 1e-6]))
end

function make_figure(; Ns_scan = (16, 24, 32, 48, 64, 96, 128, 192, 256),
                     nu::Real = 1.0, dump_path = nothing)
    A1, B1 = naive_go(32, nu); lam1 = finite_real_eigvals(A1, B1)
    A2, B2 = naive_go(48, nu); lam2 = finite_real_eigvals(A2, B2)
    pos1 = lam1[lam1 .> 1.0]; pos2 = lam2[lam2 .> 1.0]

    max_positive = Float64[]
    for N in Ns_scan
        A, B = naive_go(N, nu)
        lam = finite_real_eigvals(A, B)
        pos = lam[lam .> 1.0]
        push!(max_positive, isempty(pos) ? NaN : maximum(pos))
    end

    Ac1, Bc1 = cured_go(32, nu); lam_c1 = finite_real_eigvals(Ac1, Bc1)
    Ac2, Bc2 = cured_go(48, nu); lam_c2 = finite_real_eigvals(Ac2, Bc2)
    pc1 = lam_c1[lam_c1 .> 1.0]; pc2 = lam_c2[lam_c2 .> 1.0]

    fig = Figure(size = (1400, 420))

    ax1 = Axis(fig[1, 1]; xlabel = "mode number j", ylabel = "λ (naive)",
        title = "naive: $(length(pos1)) spurious λ > 0 at N=32, $(length(pos2)) at N=48")
    scatter!(ax1, 1:length(lam1), lam1; color = :white, strokecolor = NAVY,
             strokewidth = 1.0, markersize = 7, label = "N = 32")
    scatter!(ax1, 1:length(lam2), lam2; color = CORAL, marker = :xcross,
             markersize = 8, label = "N = 48")
    hlines!(ax1, [0.0]; color = :black, linewidth = 0.6, alpha = 0.5)
    axislegend(ax1, position = :rt, labelsize = 9)

    mask = isfinite.(max_positive)
    Nmask = Float64.(collect(Ns_scan))[mask]
    pmask = max_positive[mask]
    slope = log(pmask[end] / pmask[1]) / log(Nmask[end] / Nmask[1])
    ax2 = Axis(fig[1, 2]; xlabel = "N",
        ylabel = "spurious positive eigenvalue λ₊",
        title = "DDD scaling: empirical slope ≈ $(round(slope, digits=2)), asymptote N⁴",
        xscale = log10, yscale = log10)
    scatterlines!(ax2, Nmask, pmask; color = NAVY, markercolor = :white,
                  strokecolor = NAVY, strokewidth = 1.1, markersize = 9,
                  linewidth = 0.8, label = "spurious λ₊ (naive)")
    refy = pmask[1] .* (Nmask ./ Nmask[1]) .^ 4
    lines!(ax2, Nmask, refy; color = TEAL, linestyle = :dash, linewidth = 0.8,
           label = "N⁴ reference")
    axislegend(ax2, position = :lt, labelsize = 9)

    ax3 = Axis(fig[1, 3]; xlabel = "mode number j", ylabel = "λ (cured)",
        title = "cured: $(length(pc1)) spurious at N=32, $(length(pc2)) at N=48")
    scatter!(ax3, 1:length(lam_c1), lam_c1; color = :white, strokecolor = NAVY,
             strokewidth = 1.0, markersize = 7, label = "N = 32 (cured)")
    scatter!(ax3, 1:length(lam_c2), lam_c2; color = CORAL, marker = :xcross,
             markersize = 8, label = "N = 48 (cured)")
    hlines!(ax3, [0.0]; color = :black, linewidth = 0.6, alpha = 0.5)
    axislegend(ax3, position = :rt, labelsize = 9)

    save(joinpath(OUTPUT_DIR, "eigen_physically_spurious.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "eigen_physically_spurious.png"), fig, px_per_unit = 4)

    r = Dict{String,Any}(
        "nu" => nu,
        "Ns_scan" => collect(Float64, Ns_scan),
        "max_positive" => max_positive,
        "naive_n_positive" => Dict("32" => length(pos1), "48" => length(pos2)),
        "cured_n_positive" => Dict("32" => length(pc1), "48" => length(pc2)),
        "naive_positives_32" => pos1, "naive_positives_48" => pos2,
        "cured_spectrum_32_first5" => lam_c1[1:min(5, length(lam_c1))])
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
    @printf("[Etude 18.7]  Gottlieb-Orszag streamfunction (nu = %g)\n", r["nu"])
    @printf("  Naive positives at N=32: %d  at N=48: %d\n",
            r["naive_n_positive"]["32"], r["naive_n_positive"]["48"])
    @printf("  Cured positives at N=32: %d  at N=48: %d\n",
            r["cured_n_positive"]["32"], r["cured_n_positive"]["48"])
    println("  figure: ", joinpath(OUTPUT_DIR, "eigen_physically_spurious.pdf"))
    dp !== nothing && println("  dumped: ", dp)
end
