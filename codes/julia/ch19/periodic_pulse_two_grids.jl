# periodic_pulse_two_grids.jl
# Chapter 19: Coordinate Transformations
# Computational Etude 19.1: A periodic pulse on two grids.
#
# Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
# Part of "Computational Etudes: A Spectral Approach"

using CairoMakie
using FFTW
using Printf

const KAPPA = 80.0
const L_MAP = 0.3

target(y) = exp(-KAPPA * (1.0 - cos(y)))

arctan_tan_map(x, ell) = 2.0 * atan(ell * tan(x / 2))
arctan_tan_inv(y, ell) = 2.0 * atan(tan(y / 2) / ell)

function fourier_interp(y_nodes, f_nodes, y_eval)
    N = length(y_nodes)
    coeffs = fft(f_nodes) ./ N
    ks = [0:div(N, 2)-1; -div(N, 2):-1]
    vals = zeros(length(y_eval))
    for m in 1:N
        vals .+= real.(coeffs[m] .* exp.(im .* ks[m] .* (y_eval .+ pi)))
    end
    return vals
end

function mapped_interp(x_nodes, f_nodes, y_eval, ell)
    x_eval = arctan_tan_inv.(y_eval, ell)
    return fourier_interp(x_nodes, f_nodes, x_eval)
end

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch19", "julia")
    mkpath(outdir)

    Ns = [8, 12, 16, 24, 32, 48, 64, 96, 128]
    y_eval = collect(range(-pi + 1e-9, pi - 1e-9, length=4097))
    truth = target.(y_eval)
    err_U = zeros(length(Ns)); err_M = zeros(length(Ns))
    for (i, N) in enumerate(Ns)
        yU = [-pi + 2pi * k / N for k in 0:N-1]
        fU = target.(yU)
        err_U[i] = maximum(abs.(fourier_interp(yU, fU, y_eval) .- truth))

        x = [-pi + 2pi * k / N for k in 0:N-1]
        yM = arctan_tan_map.(x, L_MAP)
        fM = target.(yM)
        err_M[i] = maximum(abs.(mapped_interp(x, fM, y_eval, L_MAP) .- truth))
    end

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"
    GREY = colorant"#A6A6A6"
    fig = Figure(size=(1100, 760))

    Nshow = 32
    yU = [-pi + 2pi * k / Nshow for k in 0:Nshow-1]
    xU = yU
    yM = arctan_tan_map.(xU, L_MAP)
    y_line = collect(range(-pi, pi, length=401))
    x_line = collect(range(-pi + 1e-9, pi - 1e-9, length=1001))

    # (a) Pulse and grids in physical y
    ax1 = Axis(fig[1, 1];
               xlabel = "y",
               ylabel = "f(y), grid points",
               title  = "(a) Pulse and grids in physical y",
               limits = ((-pi, pi), (-0.28, 1.18)))
    h_f  = lines!(ax1, y_line, target.(y_line); color=NAVY, linewidth=1.4)
    h_uU = scatter!(ax1, yU, fill(-0.08, length(yU));
                    color=CORAL, marker=:xcross, markersize=10)
    h_uM = scatter!(ax1, yM, fill(-0.18, length(yM));
                    color=TEAL, marker=:circle, markersize=7)
    axislegend(ax1, [h_f, h_uU, h_uM],
               ["f(y)", "uniform, N = 32", "arctan/tan, ell = 0.3"];
               position=:lt, framevisible=false)

    # (b) Pulse in computational coordinate x  -- the punchline
    ax2 = Axis(fig[1, 2];
               xlabel = "x",
               ylabel = "f(y(x)), grid points",
               title  = "(b) Pulse in computational x",
               limits = ((-pi, pi), (-0.28, 1.18)))
    f_of_x = target.(arctan_tan_map.(x_line, L_MAP))
    h_ref  = lines!(ax2, y_line, target.(y_line);
                    color=GREY, linestyle=:dash, linewidth=1.0)
    h_ftil = lines!(ax2, x_line, f_of_x; color=TEAL, linewidth=1.6)
    h_xg   = scatter!(ax2, xU, fill(-0.18, length(xU));
                      color=TEAL, marker=:circle, markersize=7)
    axislegend(ax2, [h_ref, h_ftil, h_xg],
               ["original f (physical y)",
                "f(y(x)) (computational x)",
                "uniform x-grid, N = 32"];
               position=:lt, framevisible=false)

    # (c) Convergence
    ax3 = Axis(fig[2, 1];
               xlabel = "N", ylabel = "max error",
               title  = "(c) Convergence", yscale=log10)
    scatterlines!(ax3, Ns, max.(err_U, 1e-18);
                  color=CORAL, label="uniform Fourier")
    scatterlines!(ax3, Ns, max.(err_M, 1e-18);
                  color=TEAL, label="arctan/tan, ell = 0.3")
    axislegend(ax3; position=:rt, framevisible=false)

    # (d) Coefficient decay
    N = 96
    yU96 = [-pi + 2pi * k / N for k in 0:N-1]
    cU = abs.(fftshift(fft(target.(yU96)) ./ N))
    xU96 = yU96
    yM96 = arctan_tan_map.(xU96, L_MAP)
    cM = abs.(fftshift(fft(target.(yM96)) ./ N))
    kfull = fftshift([0:div(N, 2)-1; -div(N, 2):-1])
    ax4 = Axis(fig[2, 2];
               xlabel = "wave-number k", ylabel = "|c_k|",
               title  = "(d) Coefficient decay, N = 96",
               yscale = log10,
               limits = (nothing, (1e-16, 1)))
    scatter!(ax4, kfull[kfull .>= 0], cU[kfull .>= 0] .+ 1e-18;
             color=CORAL, marker=:xcross, markersize=8, label="uniform Fourier")
    scatter!(ax4, kfull[kfull .>= 0], cM[kfull .>= 0] .+ 1e-18;
             color=TEAL, marker=:circle, markersize=6, label="arctan/tan")
    axislegend(ax4; position=:rt, framevisible=false)

    save(joinpath(outdir, "periodic_pulse_two_grids.pdf"), fig)
    save(joinpath(outdir, "periodic_pulse_two_grids.png"), fig)
    @printf("[19.1-julia] saved figure to %s\n",
            joinpath(outdir, "periodic_pulse_two_grids.pdf"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
