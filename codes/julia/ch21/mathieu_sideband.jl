# mathieu_sideband.jl
# Chapter 21, Etude 21.2: Mathieu's equation and sideband truncation.
# Reproduces Boyd Figs 19.1-19.2 for ce_15.
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)

using LinearAlgebra
include(joinpath(@__DIR__, "tricks_common.jl"))
apply_theme!()
const OUT = output_dir(@__DIR__)


odd_modes(M::Int) = 2 .* (1:M) .- 1


function galerkin_full(N_modes::Int, q::Real)
    modes = odd_modes(N_modes)
    M = Diagonal(Float64.(modes) .^ 2) |> Matrix
    for i in 1:N_modes, j in 1:N_modes
        if abs(modes[i] - modes[j]) == 2
            M[i, j] = q
        end
    end
    return M, modes
end


function reference_eigenpair(N_modes::Int, q::Real, n_carrier::Int)
    M, modes = galerkin_full(N_modes, q)
    F = eigen(Symmetric(M))
    lam, V = F.values, F.vectors
    idx_carrier = findfirst(==(n_carrier), modes)
    dominant = vec([argmax(abs.(V[:, k])) for k in 1:size(V, 2)])
    candidates = findall(==(idx_carrier), dominant)
    if isempty(candidates)
        k = argmin(abs.(lam .- n_carrier^2))
    else
        comps = abs.(V[idx_carrier, candidates])
        k = candidates[argmax(comps)]
    end
    return lam[k], V[:, k], modes
end


function sideband_eigenvalue(n_carrier::Int, q::Real, half_width::Int)
    modes = [n_carrier + 2 * k for k in -half_width:half_width]
    M = Diagonal(Float64.(modes) .^ 2) |> Matrix
    for i in 1:length(modes), j in 1:length(modes)
        if abs(modes[i] - modes[j]) == 2
            M[i, j] = q
        end
    end
    lam = eigvals(Symmetric(M))
    return lam[argmin(abs.(lam .- n_carrier^2))]
end


function make_figure(; dump_path = nothing)
    n_carrier = 15
    N_modes = 64
    q_demo = 10.0
    q_axis = collect(range(0.0, 25.0; length = 200))

    lam_demo, vec_demo, modes = reference_eigenpair(N_modes, q_demo, n_carrier)
    idx_carrier = findfirst(==(n_carrier), modes)
    vec_demo = vec_demo .* sign(vec_demo[idx_carrier])
    coeff_abs = abs.(vec_demo)

    delta_full = Float64[]
    delta_3 = Float64[]
    delta_5 = Float64[]
    for q in q_axis
        lam, _, _ = reference_eigenpair(N_modes, q, n_carrier)
        push!(delta_full, lam - n_carrier^2)
        push!(delta_3, sideband_eigenvalue(n_carrier, q, 1) - n_carrier^2)
        push!(delta_5, sideband_eigenvalue(n_carrier, q, 2) - n_carrier^2)
    end

    n_small = 3
    q_small_axis = collect(range(0.0, 10.0; length = 200))
    delta_full_small = Float64[]
    delta_3_small = Float64[]
    delta_5_small = Float64[]
    for q in q_small_axis
        lam, _, _ = reference_eigenpair(N_modes, q, n_small)
        push!(delta_full_small, lam - n_small^2)
        push!(delta_3_small, sideband_eigenvalue(n_small, q, 1) - n_small^2)
        push!(delta_5_small, sideband_eigenvalue(n_small, q, 2) - n_small^2)
    end

    # NEW (panel b in 3x2): ce_3 coefficient picture at q_demo
    _, vec_demo3, _ = reference_eigenpair(N_modes, q_demo, n_small)
    idx3 = findfirst(==(n_small), modes)
    vec_demo3 = vec_demo3 .* sign(vec_demo3[idx3])
    coeff_abs3 = abs.(vec_demo3)

    # NEW (panel c in 3x2): cluster width vs q at three carriers
    q_cluster = collect(range(0.5, 50.0, length=80))
    n_cs = (3, 7, 15)
    cluster_widths = Dict(n_c => Int[] for n_c in n_cs)
    for n_c in n_cs
        idx_c = findfirst(==(n_c), modes)
        for q in q_cluster
            _, vec, _ = reference_eigenpair(N_modes, q, n_c)
            ac = abs(vec[idx_c])
            if ac == 0
                push!(cluster_widths[n_c], 0); continue
            end
            mask = abs.(vec) .> 1e-3 * ac
            offsets = modes .- n_c
            if any(mask)
                push!(cluster_widths[n_c], div(maximum(abs.(offsets[mask])), 2))
            else
                push!(cluster_widths[n_c], 0)
            end
        end
    end

    # NEW (panel e in 3x2): convergence vs sideband size, n=15, q=q_demo
    half_widths = [1, 2, 3, 4, 5, 7, 10]
    k_demo = argmin(abs.(q_axis .- q_demo))
    delta_full_15 = delta_full[k_demo]
    err_vs_hw = Float64[]
    for hw in half_widths
        d_hw = sideband_eigenvalue(n_carrier, q_demo, hw) - n_carrier^2
        push!(err_vs_hw, abs(d_hw - delta_full_15) + 1e-18)
    end

    fig = Figure(size = (1100, 1200))

    # (a) |a_n| for ce_15
    ax1 = Axis(fig[1, 1]; xlabel = "Fourier degree n  (cos basis, odd n)",
               ylabel = "|a_n|",
               title = "(a) |a_n| for ce_15 at q = $(q_demo)",
               limits = ((0, 31), nothing))
    barplot!(ax1, modes, coeff_abs; width = 1.4, color = NAVY,
             strokecolor = NAVY, strokewidth = 0.5)
    vlines!(ax1, [n_carrier]; color = CORAL, linestyle = :dash,
            linewidth = 0.8, alpha = 0.7, label = "carrier n = 15")
    axislegend(ax1, position = :rt, labelsize = 9)

    # (b) cluster width vs q
    ax2 = Axis(fig[1, 2]; xlabel = "q",
               ylabel = "cluster half-width (modes)",
               title = "(b) cluster width vs q at three carriers")
    cluster_pal = Dict(3 => CORAL, 7 => TEAL, 15 => NAVY)
    for n_c in n_cs
        scatterlines!(ax2, q_cluster, cluster_widths[n_c];
                      color = cluster_pal[n_c],
                      markercolor = :white, strokecolor = cluster_pal[n_c],
                      strokewidth = 1.0, markersize = 4, linewidth = 1.0,
                      label = "n = $n_c")
    end
    hlines!(ax2, [2]; color = :gray, linewidth = 0.4, linestyle = :dot,
            label = "5x5 box capacity (hw=2)")
    axislegend(ax2, position = :lt, labelsize = 9)

    # (c) |a_n| for ce_3
    ax3 = Axis(fig[2, 1]; xlabel = "Fourier degree n  (cos basis, odd n)",
               ylabel = "|a_n|",
               title = "(c) |a_n| for ce_3 at q = $(q_demo) (cluster spread)",
               limits = ((0, 31), nothing))
    barplot!(ax3, modes, coeff_abs3; width = 1.4, color = NAVY,
             strokecolor = NAVY, strokewidth = 0.5)
    vlines!(ax3, [n_small]; color = CORAL, linestyle = :dash,
            linewidth = 0.8, alpha = 0.7, label = "carrier n = 3")
    axislegend(ax3, position = :rt, labelsize = 9)

    # (d) delta(q) for n=15
    ax4 = Axis(fig[2, 2]; xlabel = "q", ylabel = "delta(q) = lambda - 15^2",
               title = "(d) ce_15: 5x5 matches full to plot accuracy")
    lines!(ax4, q_axis, delta_full; color = NAVY, linewidth = 1.4,
           label = "delta_full (high-N)")
    lines!(ax4, q_axis, delta_5; color = TEAL, linewidth = 1.0,
           linestyle = :dash, label = "delta_5 (5x5)")
    lines!(ax4, q_axis, delta_3; color = CORAL, linewidth = 1.0,
           linestyle = :dot, label = "delta_3 (3x3)")
    axislegend(ax4, position = :lt, labelsize = 9)

    # (e) convergence vs sideband size
    ax5 = Axis(fig[3, 1]; xlabel = "sideband half-width hw",
               ylabel = "absolute error in delta",
               yscale = log10,
               title = "(e) convergence vs sideband size (n=15, q=$(q_demo))")
    scatterlines!(ax5, half_widths, err_vs_hw; color = NAVY,
                  markercolor = :white, strokecolor = NAVY,
                  strokewidth = 1.1, markersize = 6, linewidth = 1.1,
                  label = "|delta_hw - delta_full|")
    vlines!(ax5, [2]; color = CORAL, linestyle = :dot, linewidth = 0.8,
            label = "5x5 box (hw=2)")
    axislegend(ax5, position = :rt, labelsize = 9)

    # (f) delta(q) for n=3
    ax6 = Axis(fig[3, 2]; xlabel = "q", ylabel = "delta(q) = lambda - 9",
               title = "(f) ce_3: q/n^2 large, sideband breaks")
    lines!(ax6, q_small_axis, delta_full_small; color = NAVY, linewidth = 1.4,
           label = "delta_full")
    lines!(ax6, q_small_axis, delta_5_small; color = TEAL, linewidth = 1.0,
           linestyle = :dash, label = "delta_5")
    lines!(ax6, q_small_axis, delta_3_small; color = CORAL, linewidth = 1.0,
           linestyle = :dot, label = "delta_3")
    axislegend(ax6, position = :lb, labelsize = 9)

    save(joinpath(OUT, "mathieu_sideband.pdf"), fig)
    save(joinpath(OUT, "mathieu_sideband.png"), fig, px_per_unit = 4)

    println("[Etude 21.2]  Mathieu sideband truncation")
    println("  ce_$n_carrier, q = $q_demo")
    five_largest = sortperm(coeff_abs; rev = true)[1:5]
    for j in five_largest
        println("    n = $(modes[j])  |a_n| = $(coeff_abs[j])")
    end
    println("  figure: ", joinpath(OUT, "mathieu_sideband.pdf"))

    if dump_path !== nothing
        write_json(dump_path, Dict("n_carrier" => n_carrier,
            "modes_demo" => modes, "coeff_abs_demo" => coeff_abs,
            "q_axis" => q_axis, "delta_full" => delta_full,
            "delta_5" => delta_5, "delta_3" => delta_3))
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
