# ioakimidis_root.jl
# Chapter 21, Etude 21.7: Ioakimidis non-iterative root via Chebyshev quadrature.
#
# Test function f(x) = sin(x - pi/4) / sqrt(1 + 10 x^2) on [-1, 1] has the
# unique root x = pi/4.  Ioakimidis's quotient formula recovers it
# geometrically; bisection on the same number of evaluations is much slower.
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
# Part of "Computational Etudes: A Spectral Approach"

include(joinpath(@__DIR__, "tricks_common.jl"))
apply_theme!()
const OUT = output_dir(@__DIR__)

f_test(x) = sin(x - pi/4) / sqrt(1.0 + 10.0 * x^2)
const RHO = pi / 4.0


function ioakimidis_root(N::Int, a::Real, b::Real)
    j = collect(0:2*N)
    half = 0.5 .* ((a + b) .- (a - b) .* cos.(j .* pi ./ (2.0 * N)))
    fj = f_test.(half)
    sgn = [(-1.0) ^ jj for jj in j]
    weight = ones(length(j))
    weight[1] = 0.5
    weight[end] = 0.5
    num = sum(weight .* sgn .* half ./ fj)
    den = sum(weight .* sgn        ./ fj)
    return num / den
end


function bisection_root(n_evals::Int, a::Real, b::Real)
    fa, fb = f_test(a), f_test(b)
    @assert fa * fb < 0
    for _ in 1:(n_evals - 2)
        c = 0.5 * (a + b)
        fc = f_test(c)
        if fc * fa < 0
            b = c; fb = fc
        else
            a, fa = c, fc
        end
    end
    return 0.5 * (a + b)
end


function make_figure(; dump_path = nothing)
    a, b = -1.0, 1.0
    Ns = collect(2:30)
    err_io = [abs(ioakimidis_root(N, a, b) - RHO) for N in Ns]
    err_bi = [abs(bisection_root(2*N + 1, a, b) - RHO) for N in Ns]

    fig = Figure(size = (760, 440))
    ax = Axis(fig[1, 1]; xlabel = "number of f evaluations  (= 2N+1)",
                          ylabel = "|rho_N - pi/4|", yscale = log10,
                          title = "Ioakimidis: geometric, no iteration")
    n_evals = 2 .* Ns .+ 1
    scatterlines!(ax, n_evals, err_io .+ 1e-18; color = NAVY,
        markercolor = :white, strokecolor = NAVY, strokewidth = 1.0,
        markersize = 6, linewidth = 1.0,
        label = "Ioakimidis (non-iterative)")
    scatterlines!(ax, n_evals, err_bi .+ 1e-18; color = CORAL,
        marker = :rect, markercolor = :white, strokecolor = CORAL,
        strokewidth = 1.0, markersize = 5, linewidth = 0.9,
        label = "bisection (same # of evals)")
    hlines!(ax, [1e-15]; color = :gray, linewidth = 0.4, alpha = 0.5)
    axislegend(ax, position = :lb, labelsize = 10)
    ylims!(ax, 1e-16, 1.0)

    save(joinpath(OUT, "ioakimidis_root.pdf"), fig)
    save(joinpath(OUT, "ioakimidis_root.png"), fig, px_per_unit = 4)

    println("[Etude 21.7]  Ioakimidis non-iterative root")
    for N in (3, 6, 10, 20, 30)
        rho = ioakimidis_root(N, a, b)
        println("  N = $N: rho = $rho, err = $(abs(rho - RHO))")
    end
    println("  figure: ", joinpath(OUT, "ioakimidis_root.pdf"))

    if dump_path !== nothing
        write_json(dump_path, Dict("rho_exact" => RHO,
            "Ns" => Ns, "err_ioakimidis" => err_io, "err_bisection" => err_bi))
    end
end


function parse_args()
    dp = nothing; i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--dump" && i < length(ARGS)
            dp = ARGS[i+1]; i += 2
        else; i += 1; end
    end
    return dp
end


if abspath(PROGRAM_FILE) == @__FILE__
    make_figure(; dump_path = parse_args())
end
