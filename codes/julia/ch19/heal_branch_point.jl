# heal_branch_point.jl
# Chapter 19: Coordinate Transformations
# Computational Etude 19.5: Healing the square-root branch point.
#
# Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

using CairoMakie
using FFTW
using Printf

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

function chebyshev_coeffs(v)
    N = length(v) - 1
    V = vcat(v, reverse(v[2:N]))
    A = real.(fft(V)) ./ N
    A[1] *= 0.5; A[N + 1] *= 0.5
    return A[1:N + 1]
end

function cheb_eval(a, xfine, N)
    T0 = ones(length(xfine)); T1 = copy(xfine)
    val = a[1] .* T0 .+ a[2] .* T1
    for n in 2:N
        Tk = 2 .* xfine .* T1 .- T0
        val .+= a[n + 1] .* Tk
        T0 = T1; T1 = Tk
    end
    return val
end

function chebyshev_error(N)
    _, x = cheb_matrix(N)
    v = sqrt.(max.(1 .- x .^ 2, 0.0))
    a = chebyshev_coeffs(v)
    xfine = collect(range(-1 + 1e-10, 1 - 1e-10, length=2001))
    val = cheb_eval(a, xfine, N)
    return maximum(abs.(val .- sqrt.(max.(1 .- xfine .^ 2, 0.0)))), abs.(a)
end

function mapped_error(N, Y)
    _, xi = cheb_matrix(N)
    y = Y .* xi
    f = 1.0 ./ cosh.(y)
    a = chebyshev_coeffs(f)
    xfine = collect(range(-1, 1, length=2001))
    yfine = Y .* xfine
    val = cheb_eval(a, xfine, N)
    return maximum(abs.(val .- 1.0 ./ cosh.(yfine))), abs.(a)
end

function chebyshev_pw_error(N, Xfine)
    _, x = cheb_matrix(N)
    a = chebyshev_coeffs(sqrt.(max.(1 .- x .^ 2, 0.0)))
    val = cheb_eval(a, Xfine, N)
    return abs.(val .- sqrt.(max.(1 .- Xfine .^ 2, 0.0)))
end

function mapped_pw_error(N, Y, Xfine)
    _, xi = cheb_matrix(N)
    a = chebyshev_coeffs(1.0 ./ cosh.(Y .* xi))
    Xc = clamp.(Xfine, -1.0 + 1e-12, 1.0 - 1e-12)
    yfine = atanh.(Xc)
    inside = abs.(yfine) .<= Y
    xi_fine = ifelse.(inside, yfine ./ Y, 0.0)
    val = cheb_eval(a, xi_fine, N)
    err = abs.(val .- 1.0 ./ cosh.(yfine))
    err[.!inside] .= NaN
    return err
end

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch19", "julia")
    mkpath(outdir)

    Ns = [8, 16, 24, 32, 48, 64, 96, 128]
    err_d = zeros(length(Ns)); err_m = zeros(length(Ns))
    for (i, N) in enumerate(Ns)
        err_d[i], _ = chebyshev_error(N)
        err_m[i], _ = mapped_error(N, 10.0)
    end
    _, a_d = chebyshev_error(64)
    _, a_m = mapped_error(64, 10.0)

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"
    fig = Figure(size=(1100, 760))

    # (a) Two views of the same function
    ax1 = Axis(fig[1, 1]; xlabel="coordinate", ylabel="value",
               title="(a) Two views of the same function")
    Xp = collect(range(-1, 1, length=401))
    yp = collect(range(-6, 6, length=401))
    h_g = lines!(ax1, Xp, sqrt.(max.(1 .- Xp .^ 2, 0));
                 color=CORAL, linewidth=1.4)
    h_s = lines!(ax1, yp, 1.0 ./ cosh.(yp);
                 color=TEAL, linewidth=1.4, linestyle=:dash)
    axislegend(ax1, [h_g, h_s], ["sqrt(1 - X^2)", "sech y"];
               position=:rt, framevisible=false)

    # (b) NEW: pointwise error at fixed N = 32
    Nshow = 32
    Xfine = collect(range(-1.0 + 1e-6, 1.0 - 1e-6, length=1001))
    err_dir_pw = chebyshev_pw_error(Nshow, Xfine)
    err_map_pw = mapped_pw_error(Nshow, 10.0, Xfine)
    ax2 = Axis(fig[1, 2]; xlabel="X",
               ylabel="|g(X) - g_N(X)|",
               title="(b) Pointwise error at N = $Nshow",
               yscale=log10,
               limits=((-1.0, 1.0), (1e-16, 1e0)))
    h_d = lines!(ax2, Xfine, err_dir_pw .+ 1e-18;
                 color=CORAL, linewidth=1.4)
    # mask NaN entries (outside truncation window) with finite max for plotting
    err_map_safe = ifelse.(isnan.(err_map_pw), 1e0, err_map_pw .+ 1e-18)
    h_m = lines!(ax2, Xfine, err_map_safe;
                 color=TEAL, linewidth=1.4)
    axislegend(ax2, [h_d, h_m],
               ["direct, sqrt(1-X^2)", "tanh-mapped, sech y"];
               position=:lb, framevisible=false)

    # (c) Convergence: algebraic vs geometric
    ax3 = Axis(fig[2, 1]; xlabel="N", ylabel="max error",
               xscale=log10, yscale=log10,
               title="(c) Convergence: algebraic vs geometric")
    scatterlines!(ax3, Ns, err_d; color=CORAL, label="direct")
    scatterlines!(ax3, Ns, max.(err_m, 1e-18); color=TEAL, label="tanh-mapped")
    lines!(ax3, Ns, 1.0 ./ Ns; color=NAVY, linestyle=:dot, label="1/N")
    axislegend(ax3; position=:rt, framevisible=false)

    # (d) Coefficient decay at N = 64
    ax4 = Axis(fig[2, 2]; xlabel="Chebyshev index n", ylabel="|a_n|",
               yscale=log10, title="(d) Coefficient decay, N = 64",
               limits=(nothing, (1e-17, 10.0)))
    scatter!(ax4, 0:length(a_d)-1, max.(a_d, 1e-18);
             color=CORAL, marker=:xcross, label="direct Chebyshev")
    scatter!(ax4, 0:length(a_m)-1, max.(a_m, 1e-18);
             color=TEAL, marker=:circle, label="tanh-mapped")
    axislegend(ax4; position=:rt, framevisible=false)

    save(joinpath(outdir, "heal_branch_point.pdf"), fig)
    save(joinpath(outdir, "heal_branch_point.png"), fig)
    @printf("[19.5-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
