# lebesgue_random_nodes.jl
#
# Computational experiment: Lebesgue constant for random nodes on [-1, 1].
#
# This script investigates the growth of the Lebesgue constant when interpolation
# nodes are chosen randomly (uniformly distributed) on [-1, 1]. Through Monte Carlo
# simulation, we discover the statistical behavior and derive an empirical
# asymptotic formula.
#
# The experiment:
#     1. For each N, generate N+1 random points uniformly on [-1, 1]
#     2. Compute the Lebesgue constant for these random nodes
#     3. Repeat M times to obtain statistics (mean, std, min, max)
#     4. Fit an empirical asymptotic formula to the averaged data
#
# Key finding: Random nodes exhibit exponential growth similar to equispaced nodes,
# demonstrating that the special structure of Chebyshev points is essential
# for stable interpolation.
#
# Author: Dr. Denys Dutykh
#         Mathematics Department
#         Khalifa University of Science and Technology
#         Abu Dhabi, UAE
#
# Part of the book "Computational Etudes: A Spectral Approach"
# https://github.com/dutykh/computational-etudes

using CairoMakie
using Colors
using LaTeXStrings
using Random
using Statistics
using LsqFit

# -----------------------------------------------------------------------------
# Publication-quality CairoMakie configuration
# -----------------------------------------------------------------------------
set_theme!(Theme(
    fontsize = 11,
    fonts = (regular = "CMU Serif", bold = "CMU Serif Bold", italic = "CMU Serif Italic"),
    Axis = (xlabelsize = 11, ylabelsize = 11, titlesize = 12,
            xticklabelsize = 10, yticklabelsize = 10,
            spinewidth = 0.8, xtickwidth = 0.8, ytickwidth = 0.8),
    Legend = (labelsize = 10,),
))

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
const N_FINE = 2000       # Number of points for Lebesgue function evaluation
const M_TRIALS = 200      # Number of Monte Carlo trials per N
const N_MAX = 30          # Maximum polynomial degree
const SEED = 42           # Random seed for reproducibility

# Book color scheme
const NAVY  = colorant"#142D6E"
const SKY   = colorant"#7896D2"
const CORAL = colorant"#E74C3C"
const TEAL  = colorant"#16A085"
const GOLD  = colorant"#F39C12"

# Output path
const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch04", "julia")

# -----------------------------------------------------------------------------
# Node generation
# -----------------------------------------------------------------------------
random_nodes(N, rng) = sort(2.0 .* rand(rng, N + 1) .- 1.0)
equispaced_nodes(N) = collect(range(-1, 1, length=N+1))
chebyshev_nodes(N) = [cos(j * π / N) for j in 0:N]

# -----------------------------------------------------------------------------
# Lebesgue function computation
# -----------------------------------------------------------------------------
function lagrange_basis(x_nodes, k, x_eval)
    n = length(x_nodes)
    L_k = ones(length(x_eval))
    for j in 1:n
        if j != k
            L_k .*= (x_eval .- x_nodes[j]) ./ (x_nodes[k] - x_nodes[j])
        end
    end
    return L_k
end

function lebesgue_function(x_nodes, x_eval)
    N = length(x_nodes) - 1
    Lambda = zeros(length(x_eval))
    for k in 1:N+1
        Lambda .+= abs.(lagrange_basis(x_nodes, k, x_eval))
    end
    return Lambda
end

function lebesgue_constant(x_nodes; n_eval=2000)
    x_eval = collect(range(-1, 1, length=n_eval))
    Lambda = lebesgue_function(x_nodes, x_eval)
    return maximum(Lambda)
end

# -----------------------------------------------------------------------------
# Monte Carlo experiment
# -----------------------------------------------------------------------------
function monte_carlo_lebesgue(N, M, rng; n_eval=2000)
    samples = zeros(M)
    for m in 1:M
        x_rand = random_nodes(N, rng)
        samples[m] = lebesgue_constant(x_rand, n_eval=n_eval)
    end
    return samples
end

# -----------------------------------------------------------------------------
# Main computation and visualization
# -----------------------------------------------------------------------------
function main()
    println("=" ^ 70)
    println("Computational Experiment: Lebesgue Constant for Random Nodes")
    println("=" ^ 70)

    rng = MersenneTwister(SEED)

    N_values = 2:N_MAX
    n_N = length(N_values)

    # Statistics arrays
    Lambda_mean = zeros(n_N)
    Lambda_std = zeros(n_N)
    Lambda_min = zeros(n_N)
    Lambda_max = zeros(n_N)
    Lambda_median = zeros(n_N)
    Lambda_p25 = zeros(n_N)
    Lambda_p75 = zeros(n_N)

    Lambda_equi = zeros(n_N)
    Lambda_cheb = zeros(n_N)

    println("\nRunning Monte Carlo simulation with M = $M_TRIALS trials per N...")
    println("N ranges from 2 to $N_MAX")
    println()

    for (i, N) in enumerate(N_values)
        samples = monte_carlo_lebesgue(N, M_TRIALS, rng, n_eval=N_FINE)

        Lambda_mean[i] = mean(samples)
        Lambda_std[i] = std(samples)
        Lambda_min[i] = minimum(samples)
        Lambda_max[i] = maximum(samples)
        Lambda_median[i] = median(samples)
        Lambda_p25[i] = quantile(samples, 0.25)
        Lambda_p75[i] = quantile(samples, 0.75)

        Lambda_equi[i] = lebesgue_constant(equispaced_nodes(N), n_eval=N_FINE)
        Lambda_cheb[i] = lebesgue_constant(chebyshev_nodes(N), n_eval=N_FINE)

        if (i % 5 == 0) || N == N_MAX
            println("  N = $(lpad(N, 2)): Lambda_mean = $(lpad(round(Lambda_mean[i], digits=4), 12)), std = $(lpad(round(Lambda_std[i], digits=4), 10)), range = [$(round(Lambda_min[i], digits=2)), $(round(Lambda_max[i], digits=2))]")
        end
    end

    # =========================================================================
    # Fit empirical formula: Lambda ~ a * exp(b * N) / N^c
    # =========================================================================
    println("\n" * "-" ^ 70)
    println("Fitting empirical asymptotic formula...")

    mask = collect(N_values) .>= 5
    N_fit_data = Float64.(collect(N_values)[mask])
    Lambda_fit_data = Lambda_mean[mask]

    # Exponential model
    model(N, p) = p[1] .* exp.(p[2] .* N) ./ N .^ p[3]
    p0 = [0.5, 0.5, 0.5]

    fit_result = try
        curve_fit(model, N_fit_data, Lambda_fit_data, p0)
    catch
        nothing
    end

    if fit_result !== nothing
        a_fit, b_fit, c_fit = fit_result.param
    else
        # Fallback: simple exponential fit in log space
        log_Lambda = log.(Lambda_fit_data)
        A = hcat(ones(length(N_fit_data)), N_fit_data)
        coeffs = A \ log_Lambda
        a_fit, b_fit, c_fit = exp(coeffs[1]), coeffs[2], 0.0
    end

    println("\nEmpirical formula: Lambda_N ~ $(round(a_fit, digits=4)) * exp($(round(b_fit, digits=4)) * N) / N^$(round(c_fit, digits=4))")

    # =========================================================================
    # Create figure
    # =========================================================================
    fig = Figure(size = (1100, 450))

    # -------------------------------------------------------------------------
    # Left panel: Growth of Lebesgue constant with confidence bands
    # -------------------------------------------------------------------------
    ax1 = Axis(fig[1, 1],
               xlabel = L"N" * " (polynomial degree)",
               ylabel = L"\Lambda_N" * " (Lebesgue constant)",
               title = "Growth of Lebesgue Constant",
               limits = ((0, N_MAX + 2), nothing),
               yscale = log10)

    # Shaded regions
    band!(ax1, collect(N_values), Lambda_min, Lambda_max,
          color = (GOLD, 0.2), label = "Random (min-max range)")
    band!(ax1, collect(N_values), Lambda_p25, Lambda_p75,
          color = (GOLD, 0.3), label = "Random (25th-75th percentile)")

    # Mean line
    scatterlines!(ax1, collect(N_values), Lambda_mean, color = GOLD,
                  linewidth = 1.5, markersize = 6, label = "Random (mean)")

    # Equispaced
    scatterlines!(ax1, collect(N_values), Lambda_equi, color = CORAL,
                  linewidth = 1.2, markersize = 5, marker = :rect,
                  linestyle = :dash, label = "Equispaced")

    # Chebyshev
    scatterlines!(ax1, collect(N_values), Lambda_cheb, color = SKY,
                  linewidth = 1.2, markersize = 5, marker = :utriangle,
                  linestyle = :dash, label = "Chebyshev")

    # Fitted curve
    N_fit_plot = collect(range(5, N_MAX, length=100))
    Lambda_fit_plot = a_fit .* exp.(b_fit .* N_fit_plot) ./ N_fit_plot .^ c_fit
    lines!(ax1, N_fit_plot, Lambda_fit_plot, color = NAVY, linewidth = 1.5,
           alpha = 0.8,
           label = latexstring("Fit: \\Lambda_N \\sim $(round(a_fit, digits=2)) e^{$(round(b_fit, digits=3))N}"))

    # Growth rate annotations
    text!(ax1, N_MAX - 5, Lambda_mean[end] * 2, text = L"O(e^{bN})",
          fontsize = 10, color = GOLD)
    text!(ax1, N_MAX - 5, Lambda_cheb[end] * 1.5, text = L"O(\ln N)",
          fontsize = 10, color = SKY)

    axislegend(ax1, position = :lt, framevisible = false, labelsize = 8)

    # -------------------------------------------------------------------------
    # Right panel: Distribution at fixed N
    # -------------------------------------------------------------------------
    N_hist = 15
    M_hist = 500
    rng2 = MersenneTwister(SEED + 1)
    samples_hist = monte_carlo_lebesgue(N_hist, M_hist, rng2, n_eval=N_FINE)

    mean_val = mean(samples_hist)
    median_val = median(samples_hist)
    equi_val = lebesgue_constant(equispaced_nodes(N_hist), n_eval=N_FINE)
    cheb_val = lebesgue_constant(chebyshev_nodes(N_hist), n_eval=N_FINE)

    log_samples = log10.(samples_hist)

    ax2 = Axis(fig[1, 2],
               xlabel = L"\log_{10}(\Lambda_{15})",
               ylabel = "Count",
               title = "Distribution of \$\\Lambda_N\$ for N = $N_hist ($M_hist trials)")

    hist!(ax2, log_samples, bins = 35, color = (GOLD, 0.7),
          strokecolor = :white, strokewidth = 0.5)

    # Vertical lines for statistics
    vlines!(ax2, log10(median_val), color = TEAL, linewidth = 2, linestyle = :dash,
            label = "Median = 10^$(round(log10(median_val), digits=1))")
    vlines!(ax2, log10(mean_val), color = NAVY, linewidth = 2,
            label = "Mean = 10^$(round(log10(mean_val), digits=1))")
    vlines!(ax2, log10(equi_val), color = CORAL, linewidth = 2, linestyle = :dot,
            label = "Equispaced = 10^$(round(log10(equi_val), digits=1))")
    vlines!(ax2, log10(cheb_val), color = SKY, linewidth = 2, linestyle = :dashdot,
            label = "Chebyshev = $(round(cheb_val, digits=1))")

    # Spread annotation
    spread = log10(maximum(samples_hist)) - log10(minimum(samples_hist))
    text!(ax2, 0.05, 0.95, text = "Spread: $(round(spread, digits=1)) orders\nof magnitude",
          fontsize = 9, space = :relative, align = (:left, :top))

    axislegend(ax2, position = :rt, framevisible = false, labelsize = 8)

    # Save figure
    mkpath(OUTPUT_DIR)
    save(joinpath(OUTPUT_DIR, "lebesgue_random_nodes.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "lebesgue_random_nodes.png"), fig, px_per_unit = 4)

    println("\nFigure saved to: ", joinpath(OUTPUT_DIR, "lebesgue_random_nodes.pdf"))

    # =========================================================================
    # Print summary
    # =========================================================================
    println("\n" * "=" ^ 70)
    println("Summary Statistics")
    println("=" ^ 70)
    println("$(rpad("N", 4)) $(lpad("Mean", 12)) $(lpad("Std", 12)) $(lpad("Min", 12)) $(lpad("Max", 12)) $(lpad("Equispaced", 12))")
    println("-" ^ 70)
    for i in 1:3:n_N
        N = N_values[i]
        println("$(lpad(N, 4)) $(lpad(round(Lambda_mean[i], digits=4), 12)) $(lpad(round(Lambda_std[i], digits=4), 12)) $(lpad(round(Lambda_min[i], digits=4), 12)) $(lpad(round(Lambda_max[i], digits=4), 12)) $(lpad(round(Lambda_equi[i], digits=4), 12))")
    end
    println("-" ^ 70)

    println("\n" * "=" ^ 70)
    println("EMPIRICAL ASYMPTOTIC FORMULA")
    println("=" ^ 70)
    println("\nFor random nodes uniformly distributed on [-1, 1]:")
    println("\n    Lambda_N ~ $(round(a_fit, digits=4)) * exp($(round(b_fit, digits=4)) * N) / N^$(round(c_fit, digits=4))")
    println("\nor approximately:")
    println("\n    Lambda_N ~ $(round(a_fit, digits=4)) * $(round(exp(b_fit), digits=4))^N")
    println("\nThis confirms exponential growth similar to equispaced nodes.")
    println("=" ^ 70)

    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
