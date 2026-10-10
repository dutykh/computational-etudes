#!/usr/bin/env julia
#=
trap_fft_coefficients.jl
Chapter 16: Integration of Periodic Functions

Computational Etude 16.8: FFT Computation of Fourier Coefficients.

The FFT of a sample vector, divided by N, is exactly the periodic
trapezoidal rule applied to the Fourier coefficient integrals.
Test on f(x) = 1/(2 - cos x), whose Fourier coefficients are known
in closed form: c_k = (1/sqrt(3)) * (2 - sqrt(3))^|k|.

Generates Figure 16.8: fft_coefficients.pdf  (two-panel figure)

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
=#

using FFTW
using Plots
using LaTeXStrings

gr()

const NAVY  = RGB(20/255, 45/255, 110/255)
const CORAL = RGB(231/255, 76/255, 60/255)

const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch16", "julia")
mkpath(OUTPUT_DIR)

function main()
    a = 2.0
    f(x) = 1 / (a - cos(x))
    r = a - sqrt(a^2 - 1)
    N = 32
    θ = 2π * (0:N-1) / N
    c_fft = fft(f.(θ)) / N
    c_fft_sym = fftshift(c_fft)
    k_sym = collect(-N÷2 : N÷2 - 1)
    c_exact = [(1 / sqrt(a^2 - 1)) * r ^ abs(k) for k in k_sym]
    err = abs.(c_fft_sym .- c_exact)
    err = max.(err, 1e-17)

    println("Max error in resolved band: ", maximum(err[2:end-1]))

    p1 = plot(k_sym, abs.(c_exact),
              yscale = :log10,
              xlabel = L"Mode index $k$",
              ylabel = L"$|\hat f_k|$",
              title = L"(a) Fourier-coefficient magnitudes for $1/(2 - \cos x)$",
              label = "exact",
              color = NAVY, linewidth = 1.2,
              grid = true, gridalpha = 0.3, framestyle = :box,
              ylims = (1e-12, 1e1), legend = :topright)
    scatter!(p1, k_sym, abs.(c_fft_sym),
             color = CORAL, markersize = 5, markerstrokecolor = CORAL,
             label = "FFT, N = $N")

    p2 = plot(k_sym, err,
              yscale = :log10,
              xlabel = L"Mode index $k$",
              ylabel = L"Absolute error in $\hat f_k$",
              title = L"(b) FFT error: machine precision in the resolved band",
              label = L"$|\hat f_k^{\mathrm{FFT}} - \hat f_k|$",
              marker = :circle, markersize = 5, color = CORAL,
              linewidth = 1.0, legend = :topright,
              grid = true, gridalpha = 0.3, framestyle = :box,
              ylims = (1e-18, 1e0))

    p = plot(p1, p2, layout = (1, 2), size = (1100, 450))
    savefig(p, joinpath(OUTPUT_DIR, "fft_coefficients.pdf"))
    savefig(p, joinpath(OUTPUT_DIR, "fft_coefficients.png"))
    println("\nFigure saved to ", joinpath(OUTPUT_DIR, "fft_coefficients.pdf"))
end

main()
