# crossover_truncation.jl
# Chapter 21, Etude 21.4: cross-over truncation.
#
# Boyd's cartoon: a_n = 10 exp(-n/3) + 1e-6 / n^5.
# Asymptotically the algebraic tail wins, but only past n ~ 120.
# Real example: f(x) = exp(-((x-0.3)/0.1)^2) + 1e-7 (x+1)^(1/3) shows the
# same two-slope behaviour in actual Chebyshev coefficients.
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)

using FFTW
include(joinpath(@__DIR__, "tricks_common.jl"))
apply_theme!()
const OUT = output_dir(@__DIR__)


function cgl(N)
    return cos.(pi .* (0:N) ./ N)
end


function cheb_coeffs(v::AbstractVector)
    N = length(v) - 1
    ext = vcat(v, reverse(v[2:N]))
    A = real.(fft(ext)) ./ N
    A[1] *= 0.5; A[N+1] *= 0.5
    return A[1:N+1]
end


f_real(x) = exp(-((x - 0.3) / 0.1) ^ 2) + 1.0e-7 * (x + 1.0) ^ (1.0 / 3.0)


function make_figure(; dump_path = nothing)
    n_axis = collect(1:400)
    head_geom = 10.0 .* exp.(-n_axis ./ 3.0)
    tail_alg  = 1.0e-6 ./ n_axis .^ 5
    a_total = head_geom .+ tail_alg
    log_diff = log10.(tail_alg .+ 1e-300) .- log10.(head_geom .+ 1e-300)
    n_cross = float(n_axis[argmin(abs.(log_diff))])

    N = 400
    xg = cgl(N)
    a = cheb_coeffs(f_real.(xg))

    fig = Figure(size = (1150, 440))

    ax1 = Axis(fig[1, 1]; xlabel = "index n", ylabel = "|a_n|",
               yscale = log10,
               title = "Panel A.  Boyd's cartoon: slope change at n_cross",
               limits = ((0, 400), (1e-22, 50)))
    lines!(ax1, n_axis, head_geom; color = NAVY, linewidth = 1.2,
           label = "10 exp(-n/3) (geometric head)")
    lines!(ax1, n_axis, tail_alg; color = CORAL, linewidth = 1.2,
           label = "1e-6/n^5 (algebraic tail)")
    lines!(ax1, n_axis, a_total; color = PURPLE, linewidth = 1.6,
           label = "a_n = head + tail")
    vlines!(ax1, [n_cross]; color = :gray, linestyle = :dash,
            linewidth = 0.8, alpha = 0.6,
            label = "n_cross ~ $(Int(n_cross))")
    axislegend(ax1, position = :rt, labelsize = 9)

    ax2 = Axis(fig[1, 2]; xlabel = "Chebyshev degree n", ylabel = "|a_n|",
               yscale = log10,
               title = "Panel B.  f(x) = e^(-((x-0.3)/0.1)^2) + 1e-7 (x+1)^(1/3)",
               limits = ((0, N+10), (1e-18, 5)))
    n_axis2 = collect(0:N)
    scatterlines!(ax2, n_axis2, abs.(a) .+ 1e-300; color = NAVY,
        markersize = 2, linewidth = 1.0, label = "computed |a_n|")
    axislegend(ax2, position = :rt, labelsize = 9)

    save(joinpath(OUT, "crossover_truncation.pdf"), fig)
    save(joinpath(OUT, "crossover_truncation.png"), fig, px_per_unit = 4)

    println("[Etude 21.4]  cross-over truncation")
    println("  cartoon n_cross approx $(Int(n_cross))")
    println("  real example: |a_50| = $(abs(a[51])), |a_200| = $(abs(a[201])), "
            * "|a_400| = $(abs(a[401]))")
    println("  figure: ", joinpath(OUT, "crossover_truncation.pdf"))

    if dump_path !== nothing
        write_json(dump_path, Dict("n_cross_cartoon" => n_cross,
            "a_real" => abs.(a)))
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
