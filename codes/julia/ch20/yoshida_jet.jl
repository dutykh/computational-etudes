# yoshida_jet.jl
# Chapter 20: Spectral Methods on Unbounded Intervals
# Computational Etude 20.8: Pick the right half of the basis (odd SB for 1/y tail).
#
# Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

using CairoMakie
using LinearAlgebra
using Printf

sb_eval(n, y, ell) = sin.((n + 1) .* (pi / 2 .- atan.(y ./ ell)))

function sb_ddot(n, y, ell)
    t = pi / 2 .- atan.(y ./ ell)
    denom = ell^2 .+ y .^ 2
    dt = -ell ./ denom
    dt2 = 2 .* ell .* y ./ denom .^ 2
    c = cos.((n + 1) .* t); s = sin.((n + 1) .* t)
    return -((n + 1)^2) .* s .* dt .^ 2 .+ (n + 1) .* c .* dt2
end

function assemble(N, ell)
    idx = 2 .* (0:N - 1) .+ 1
    t_nodes = (1:N) .* pi ./ (2 * (N + 1))
    y_nodes = ell ./ tan.(t_nodes)
    A = zeros(N, N)
    for j in 1:N
        A[:, j] = sb_ddot(idx[j], y_nodes, ell) .- (y_nodes .^ 2) .* sb_eval(idx[j], y_nodes, ell)
    end
    c = A \ Vector(y_nodes)
    return c, idx
end

function evaluate(c, idx, y, ell)
    v = zeros(length(y))
    for k in 1:length(c)
        v .+= c[k] .* sb_eval(idx[k], y, ell)
    end
    return v
end

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch20", "julia")
    mkpath(outdir)
    ell = 3.0
    Ns = [2, 3, 4, 5, 6, 8, 10]
    y_fine = collect(range(0.1, 20.0, length=4001))
    c_ref, idx_ref = assemble(21, ell)
    v_ref = evaluate(c_ref, idx_ref, y_fine, ell)
    errs = [begin
        c, idx = assemble(N, ell)
        maximum(abs.(evaluate(c, idx, y_fine, ell) .- v_ref))
    end for N in Ns]

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"
    PURPLE = colorant"#8E44AD"
    fig = Figure(size=(1020, 340))

    ax1 = Axis(fig[1, 1], xlabel="y", ylabel="v(y)", title="(a) Yoshida jet")
    y_p = collect(range(0.05, 8.0, length=401))
    lines!(ax1, y_p, evaluate(c_ref, idx_ref, y_p, ell);
           color=NAVY, linewidth=1.4, label="ref (21 SB)")
    for (N, col, st) in [(2, CORAL, :dash), (3, TEAL, :dashdot), (4, PURPLE, :dot)]
        c, idx = assemble(N, ell)
        lines!(ax1, y_p, evaluate(c, idx, y_p, ell);
               color=col, linestyle=st, linewidth=1.0, label="N=$N")
    end
    axislegend(ax1; position=:rb)

    ax2 = Axis(fig[1, 2], xlabel="N (odd SB only)", ylabel="max error",
               yscale=log10, title="(b) Odd-SB resolves 1/y tail")
    scatterlines!(ax2, Ns, max.(errs, 1e-18); color=TEAL)

    save(joinpath(outdir, "yoshida_jet.pdf"), fig)
    save(joinpath(outdir, "yoshida_jet.png"), fig)
    @printf("[20.8-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
