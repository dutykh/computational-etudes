# symbolic_quartic_oscillator.jl
# Chapter 21, Etude 21.10: quartic oscillator on the infinite line via the
# rational-Chebyshev map y = ell x / sqrt(1 - x^2), ell = 2.
# Symbolic 5-term Galerkin produces a 5th-degree secular polynomial in E
# (Boyd Eq 20.16); its roots are compared to Bender-Orszag references.
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)

using SymPy
using Polynomials
include(joinpath(@__DIR__, "tricks_common.jl"))
apply_theme!()
const OUT = output_dir(@__DIR__)


function build_secular_determinant()
    @syms x::real E::real
    ell = SymPy.Sym(2)
    a = [SymPy.symbols("a_$(k)") for k in 1:5]
    u = (1 - x^2) * sum(a[k] * x^(2(k-1)) for k in 1:5)
    u_x  = diff(u, x)
    u_xx = diff(u, x, 2)

    R = (1 - x^2)^4 * ((1 - x^2) * u_xx - 3 * x * u_x) / ell^2 +
        ((1 - x^2)^2 * E - ell^4 * x^4) * u

    eqs = [SymPy.integrate(x^(2j) * R, (x, -1, 1)) for j in 0:4]

    M = SymPy.Sym[simplify(eqs[i].coeff(a[j])) for i in 1:5, j in 1:5]
    Mmat = SymPy.Sym.(reshape(M, 5, 5))
    D_E = simplify(symbolic_det(Mmat))
    poly_E = sympy.Poly(D_E, E)
    return poly_E, x, E
end


# Hand-rolled cofactor-expansion determinant for symbolic matrices.
function symbolic_det(A::AbstractMatrix{<:SymPy.Sym})
    n = size(A, 1)
    n == 1 && return A[1, 1]
    s = SymPy.Sym(0)
    for j in 1:n
        sub = A[setdiff(1:n, 1), setdiff(1:n, j)]
        s += (-1)^(j + 1) * A[1, j] * symbolic_det(sub)
    end
    return s
end


function main_run(dump_path)
    poly_E, x_sym, E_sym = build_secular_determinant()
    coeffs_E = poly_E.all_coeffs()       # high-degree first
    println("Symbolic secular determinant (high degree -> low):")
    for (k, c) in enumerate(coeffs_E)
        println("   E^$(length(coeffs_E) - k): $(simplify(c))")
    end
    norm_const = Float64(coeffs_E[end])
    println("\nBoyd Eq 20.16 form  D(E) / D(0):")
    for (k, c) in enumerate(reverse(coeffs_E))
        println("   E^$(k - 1): $(Float64(c) / norm_const)")
    end

    coeffs_num = Float64.(coeffs_E)
    poly = Polynomial(reverse(coeffs_num))   # Polynomials.jl: low to high
    rts = Polynomials.roots(poly)
    rts_real = sort([real(r) for r in rts if abs(imag(r)) < 1e-6])

    bender_orszag = [1.060362, 7.45570, 16.2618, 26.528, 37.92]
    println("\nEigenvalues from secular determinant vs Bender-Orszag:")
    for (r, ref) in zip(rts_real, bender_orszag)
        rel = (r - ref) / ref
        println("   E_num = $(round(r; digits = 4))    E_exact = $ref    rel err = $(round(rel; sigdigits=3))")
    end

    fig = Figure(size = (1150, 440))

    ax1 = Axis(fig[1, 1]; xlabel = "eigenvalue E", ylabel = "D(E)",
               title = "secular determinant D(E) and its roots")
    Es = collect(range(0.0, 50.0; length = 1000))
    Ds = poly.(Es)
    lines!(ax1, Es, Ds; color = NAVY, linewidth = 1.2, label = "D(E)")
    hlines!(ax1, [0.0]; color = :gray, linewidth = 0.4, alpha = 0.5)
    for ref in bender_orszag
        vlines!(ax1, [ref]; color = TEAL, linewidth = 0.8, alpha = 0.5,
                linestyle = :dot)
    end
    scatter!(ax1, rts_real, zeros(length(rts_real));
             color = CORAL, marker = :xcross, markersize = 14,
             label = "numerical roots of D")
    scatter!(ax1, bender_orszag, zeros(length(bender_orszag));
             color = :white, strokecolor = TEAL, strokewidth = 1.0, markersize = 10,
             label = "Bender-Orszag E_n")
    axislegend(ax1, position = :lt, labelsize = 9)

    ns = [0, 2, 4, 6, 8]
    rels = abs.([(r - ref) / ref for (r, ref) in zip(rts_real, bender_orszag)])
    ax2 = Axis(fig[1, 2]; xlabel = "physical mode index n", ylabel = "relative error",
               yscale = log10,
               title = "trust the lower spectrum: bottom 3 good, top 2 unusable",
               limits = ((-0.5, 8.5), (1e-4, 1e2)))
    scatterlines!(ax2, ns, rels; color = NAVY, markercolor = :white,
        strokecolor = NAVY, strokewidth = 1.0, markersize = 8, linewidth = 1.0,
        label = "|E_n^sym - E_n| / E_n")
    hlines!(ax2, [0.01]; color = TEAL, linewidth = 0.4, alpha = 0.5,
            label = "1% threshold")
    axislegend(ax2, position = :lt, labelsize = 9)

    save(joinpath(OUT, "symbolic_quartic_oscillator.pdf"), fig)
    save(joinpath(OUT, "symbolic_quartic_oscillator.png"), fig, px_per_unit = 4)
    println("  figure: ", joinpath(OUT, "symbolic_quartic_oscillator.pdf"))

    if dump_path !== nothing
        write_json(dump_path, Dict("ell" => 2.0,
            "secular_coeffs" => coeffs_num, "roots" => rts_real,
            "bender_orszag" => bender_orszag, "rel_err" => rels))
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
    main_run(parse_args())
end
