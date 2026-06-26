// textbook/chapters/special_tricks.typ
// Chapter 21 (bonus): Special Tricks for the Spectral Researcher
// Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
// Email: denys.dutykh@ku.ac.ae
// Homepage: https://www.denys-dutykh.com/
// Last modified: April 2026

#import "../styles/template.typ": dropcap, etude-conclusion, idx, chapter-abstract, exercise, hint-for


= Special Tricks for the Spectral Researcher <ch-special-tricks>

#chapter-abstract(keywords: [Sideband truncation · Singularity subtraction · Radiation basis functions · Spectral root-finding · Hilbert transform · Tau method])[
This closing, bonus chapter is a farewell handbook rather than a lecture: a curated collection of special tricks, drawn from Boyd's chapters 19 and 20 and presented as computational études in the book's three languages. Each trick exists because some standard recipe of the preceding chapters fails for a specific structural reason, so the chapter teaches diagnosis by symptom --- slow coefficient decay, a wrong dominant mode, an untrustworthy upper spectrum, expensive function evaluations, an endpoint singularity, an oscillatory tail. The remedies follow in turn: sideband truncation when the dominant Fourier coefficient sits near a high carrier rather than the fundamental; singularity subtraction to peel a known algebraic singularity from a smooth remainder; radiation basis functions that make asymptotic plane waves explicit unknowns in scattering problems; a Chebyshev surrogate that replaces a costly determinant search with non-iterative root-finding; and the Hilbert transform computed as a one-line FFT operation. A set of six precepts shows how the symbolic, small-$N$ regime reverses the usual large-$N$ advice, and the Lanczos tau philosophy --- solve the exact solution of a nearby problem rather than an approximate solution of the original --- supplies the chapter's guiding maxim.
]

#dropcap[A spectral method becomes elegant when the discretisation mirrors the structure of the solution. Twenty chapters of this book have built up a vocabulary --- bases, maps, weightings, eigenproblems, BVPs --- in which 'mirror the structure' can be made operational. The chapter you are reading now is different from the others. It is not a lecture chapter. It is a farewell handbook: a small collection of *special tricks* drawn from John Boyd's chapters 19 and 20 of @Boyd2000, presented as computational études with the same three-language treatment as the rest of the book. Each étude here exists because some standard recipe of chapters 1--20 fails for some specific structural reason --- a carrier-frequency mode that hides far from the lowest harmonic, a corner singularity that defeats the algebra, an asymptotic plane wave that no rational basis can represent, an expensive determinant whose roots demand many evaluations, a low-$N$ symbolic question that asks for the *opposite* of what large-$N$ wisdom recommends. The chapter is bonus material; it does not appear in the lecture schedule, and you will not be examined on it. It exists for the moment, two or five or ten years from now, when you find yourself stuck at a research problem that does not respond to the standard methods, and you need to remember that a method is a marriage of basis and structure --- and that *the right trick is not a shortcut around understanding; it is understanding made operational*.]

By the end of this chapter, you should be able to:

1. Diagnose a failing spectral computation by *symptom* (slow coefficient decay, wrong dominant mode, untrustworthy upper spectrum, expensive evaluations, endpoint singularity, oscillatory tail) rather than by trial and error.
2. Apply *sideband truncation#idx("sideband truncation")* when the dominant Fourier coefficient sits near a high carrier rather than the fundamental.
3. Apply *singularity subtraction#idx("singularity subtraction")* to isolate a known algebraic singularity from a smooth remainder.
4. Build *radiation basis#idx("radiation basis") functions* into a scattering formulation so the asymptotic plane waves become explicit unknowns rather than impossible targets for a rational basis.
5. Replace a costly determinant search by a *Chebyshev surrogate* and find its roots non-iteratively.
6. Compute a *Hilbert transform#idx("Hilbert transform")* of an analytic periodic function as a one-line FFT operation with super-geometric accuracy.
7. Recognise the *six precepts* in which the symbolic small-$N$ regime reverses the usual numerical advice (Galerkin over collocation, Legendre or Gegenbauer over Chebyshev, polynomialise transcendental structure, rationalise irrationals, exploit symmetry aggressively, prefer basis recombination over boundary bordering).
8. Read the *Lanczos#idx("Lanczos") $tau$ philosophy* operationally: solve the *exact* solution of a *nearby* problem instead of the approximate solution of the original.

The chapter ends with a deliberately brief literature pointer (which the author plans to enrich in a later revision) and a one-page closing checklist for the student moving into research.

== How to Read This Chapter --- A Diagnostic Prelude <sec-st-diagnostic>

The chapter is organised by *failure mode*, not by basis or operator. Each section answers a question of the form *what went wrong with the naive spectral method, and what structural modification fixes it?* The diagnostic table below summarises the seven motifs you will meet, the symptom that should suggest each, and the section in which it appears. We urge you to use the table backwards: when you encounter one of these symptoms in your own research, scan the table for the structural clue it carries.

#figure(
  table(
    columns: (1fr, 1.1fr, 1.0fr, 0.55fr),
    align: (left + horizon, left + horizon, left + horizon, center + horizon),
    stroke: none,
    inset: (x: 6pt, y: 5pt),
    fill: (col, row) => if row == 0 { rgb(20, 45, 110) }
                       else if calc.rem(row, 2) == 0 { rgb(248, 250, 254) }
                       else { white },
    [#text(white, weight: "bold")[Symptom]],
    [#text(white, weight: "bold")[Structural clue]],
    [#text(white, weight: "bold")[Likely remedy]],
    [#text(white, weight: "bold")[Section]],
    [Coefficients peak at a high mode, not the fundamental], [a fast carrier modulated by a slow envelope], [sideband truncation around the carrier], [@sec-st-sideband],
    [Slow algebraic decay of $|a_n|$, root in $1/n^p$], [endpoint or corner singularity], [singularity subtraction; tailored singular basis], [@sec-st-tailored],
    [Slope change in the coefficient plot], [smooth bulk plus weak singularity (cross-over truncation)], [accept algebraic decay only beyond $n_(times)$; reformulate before refining], [@sec-st-tailored],
    [Oscillating asymptote, no decay at infinity], [radiation condition or scattering plane wave], [augment the basis with one-sided radiation functions $H(x) cos(k x), H(x) sin(k x)$], [@sec-st-radiation],
    [Each evaluation of the residual is expensive], [parameter-dependent determinant whose root we want], [Chebyshev surrogate and colleague-matrix root-finding], [@sec-st-rootfinding],
    [Need a Hilbert transform of an analytic periodic function], [Fourier multiplier $-i thin op("sgn")(k)$ is exact on the circle], [FFT --- multiply by $-i thin op("sgn")(k)$ --- inverse FFT], [@sec-st-transforms],
    [Working symbolically with $N$ small], [the usual numerical defaults invert: prefer Galerkin, Legendre, polynomialisation, rationalisation], [the six precepts of @sec-st-symbolic], [@sec-st-symbolic],
    [Polynomial cannot satisfy the original ODE], [solve a *modified* ODE *exactly*], [Lanczos $tau$ method (@sec-st-tau)], [@sec-st-tau],
  ),
  caption: [Diagnostic roadmap for the special tricks of this chapter. The seven rows correspond to the seven failure modes addressed in sections @sec-st-meta to @sec-st-tau.],
) <tab-st-diagnostic>

The chapter is not a substitute for the underlying analysis you have already done; it is an organisation of the moments at which that analysis becomes operational. We start with four meta-tricks that recur as a thread through the entire chapter.

== Meta-tricks Before the Main Tricks <sec-st-meta>

Boyd opens his chapter on special tricks by listing four meta-ideas which, although discussed elsewhere in his book, are too important to omit here. We adopt the same opening and place each in conversation with études from the present book:

+ *Nonlinear degrees of freedom.* The unknowns of a spectral approximation need not be only linear series coefficients. A *width* parameter, a *scale* parameter, or any other parameter that enters the ansatz nonlinearly is fair game. Localized pulses, boundary layer#idx("boundary layer")s, and solitary waves are natural beneficiaries; @etude-st-quartic in @sec-st-symbolic uses a map parameter $ell$ as exactly such a free width.

+ *Amplitude--phase splitting for slowly-decaying oscillations.* When oscillations persist with slowly-decaying amplitude as $|y| arrow.r infinity$, no rational basis can represent the infinitely many crests --- but representing *amplitude* and *phase* separately can. This is exactly what @etude-un-j0 of chapter 20 does for $sqrt(1+y) thin J_0 (y)$.

+ *Complex-plane detours around interior singularities.* When the differential operator has a singular coefficient on the *interior* of the integration interval, a *contour* in the complex plane that arcs around the singularity gives a numerically tractable problem on a curve where everything is finite. This is the trick of @etude-complex-detour in chapter 18.

+ *Coordinate maps and point clustering.* When the action is concentrated in a small region of the domain, a coordinate transformation that clusters grid points there is often more effective than refining $N$ uniformly. Chapter 19 in its entirety is devoted to this; see also Étude 18.4.1 (the rational Chebyshev $T B_n$ illustration).

These four ideas, taken together, license the rest of the chapter: every trick that follows is some specialisation of *change the discretisation, not the resolution*.

We open with the smallest possible computational rescue story.

=== Computational Étude 21.1: A Rescue Story in Miniature <etude-st-rescue>

Approximate $f(y) = sech(y)$ on a *truncated* finite interval $[-L, L]$ by Chebyshev collocation, *same* $N = 48$, *two* truncation lengths:

- $L = 30$: a conservative, defensive choice --- 'make the window large enough to contain everything'. The pulse $sech(y)$ has support of order $6$ around $y = 0$, so $88%$ of the Chebyshev grid lies outside the support of the function. That resolution is wasted, and it costs us decimals.
- $L = 8$: just enough room that $sech(L) approx 6.7 times 10^(-4)$ sits at the level of expected truncation noise. Now the Chebyshev points live in the support of $f$, and the polynomial sees a problem it can resolve.

The chapter's first lesson lives in the contrast: *the truncation length $L$ is part of the discretisation, not a cosmetic detail*.

In Python:

```python
def truncated_chebyshev_error(N, L, y_ref):
    t_grid = cgl(N)
    y_grid = L * t_grid                 # CGL points pulled back to [-L, L]
    a = cheb_coeffs(f(y_grid))
    t_ref = y_ref / L
    inside = np.abs(t_ref) <= 1.0
    approx = np.full_like(y_ref, np.nan)
    approx[inside] = reconstruct(a, t_ref[inside])
    err = np.abs(approx - f(y_ref))
    err[~inside] = np.nan
    return a, approx, err
```

In MATLAB:

```matlab
function [a, err] = trunc_cheb_err(N, L, y_ref)
    t_grid = cgl(N);
    y_grid = L * t_grid;
    a = cheb_coeffs(sech(y_grid));
    t_ref = y_ref / L;
    inside = abs(t_ref) <= 1.0;
    approx = nan(size(y_ref));
    approx(inside) = clenshaw(a, t_ref(inside));
    err = abs(approx - sech(y_ref));
    err(~inside) = NaN;
end
```

In Julia:

```julia
function trunc_cheb_err(N::Int, L::Real, y_ref::AbstractVector)
    t_grid = cgl(N)
    y_grid = L .* t_grid
    a = cheb_coeffs(sech.(y_grid))
    t_ref = y_ref ./ L
    inside = abs.(t_ref) .<= 1.0
    approx = fill(NaN, length(y_ref))
    approx[inside] .= clenshaw(a, t_ref[inside])
    err = abs.(approx .- sech.(y_ref))
    err[.!inside] .= NaN
    return a, err
end
```

#figure(
  image("../figures/ch21/python/rescue_naive_vs_tailored.pdf", width: 100%),
  caption: [Étude 21.1: same function, same $N=48$, two truncation lengths. *Left*: $sech(y)$ overlaid with the naive grid ($L=30$, coral circles) and the tailored grid ($L=8$, teal triangles). The naive grid lavishes resolution where the function is essentially zero. *Right*: pointwise error on the inner window of each truncation: naive max error $approx 9.6 times 10^(-2)$, tailored max error $approx 1.0 times 10^(-4)$ --- a factor of $approx 10^3$ recovered by editing one constant.],
) <fig-st-rescue>

The three scripts are available in:
- `codes/python/ch21/rescue_naive_vs_tailored.py`
- `codes/matlab/ch21/rescue_naive_vs_tailored.m`
- `codes/julia/ch21/rescue_naive_vs_tailored.jl`

#etude-conclusion[
  The moral is *deliberately small*. We have not changed the basis, the resolution, or the algorithm; we have changed only the *length scale* on which the basis lives. The error has dropped by a factor of $approx 10^3$ from a single edited constant. The remainder of the chapter is a sequence of similar interventions: each one shifts the discretisation toward the structure that the solution already has.
]

== Sideband Truncation: Mathieu's Equation <sec-st-sideband>

Truncate a Fourier expansion at the *first* $N$ harmonics. This is the standard procedure of chapters 9 and 11. For Mathieu's equation
$ u_(x x) + [lambda - 2 q cos(2 x)] u = 0, $ <eq-st-mathieu>
on the periodic interval $[0, 2 pi]$ this default is *catastrophically wasteful* whenever the eigenfunction we want is *high*. Boyd's archetypal example @Boyd2000 is $"ce"_(15)$, the eigenfunction whose unperturbed ($q=0$) form is $cos(15 x)$. For moderate $q$, perturbation theory shows that the largest Fourier coefficient sits not at $n=0$ or $n=1$ but at $n=15$ --- with two small *sidebands* at $n=13, 17$ and a few smaller ones at $n=11, 19$, falling off rapidly. The dominant cluster occupies a *band* of width $approx 9$ centred at the carrier, not the customary band of width $N$ at the bottom of the spectrum. Truncating at the first $N=10$ Fourier modes would discard the carrier itself; truncating at the first $N=20$ modes would carry a great many irrelevant low harmonics in order to reach $n=15$.

#block(stroke: 0.8pt + rgb(20, 45, 110), radius: 3pt, inset: 10pt,
  fill: rgb(248, 250, 254))[
  *Definition (sideband truncation, after Boyd).*  When the basis set is restricted to basis functions of the form $phi_(n plus.minus m)(x)$ where $m << n$ and the dominant coefficient sits at a *carrier* $n$, the truncation is called a *sideband truncation*. The eigenfunction $phi_n (x)$ is the *fundamental*; the few neighbours $phi_(n plus.minus 2)(x), phi_(n plus.minus 4)(x), thin dots.h.c$ are the *sidebands*. The cost of locating the perturbed eigenvalue is the cost of diagonalising a small matrix of size $2 m + 1$, instead of a large matrix of size $N$.
]

The trigonometric identity $cos(2 x) cos(n x) = (1 \/ 2) [cos((n+2) x) + cos((n-2) x)]$ shows that the perturbation $-2 q cos(2 x)$ couples mode $n$ only to modes $n plus.minus 2$. The eigenvalue equation (@eq-st-mathieu) therefore decouples into four parity classes (even-even, even-odd, odd-odd, odd-even), and within each class the coupling is tridiagonal in the parity-restricted index. For $"ce"_(15)$ (which lives in the *even-odd* class with basis ${cos(x), cos(3 x), cos(5 x), thin dots.h.c}$), the natural sideband matrix is $5 times 5$:
$ mat(
  lambda - (n-4)^2, -q, 0, 0, 0;
  -q, lambda - (n-2)^2, -q, 0, 0;
  0, -q, lambda - n^2, -q, 0;
  0, 0, -q, lambda - (n+2)^2, -q;
  0, 0, 0, -q, lambda - (n+4)^2
) mat(a_(n-4); a_(n-2); a_n; a_(n+2); a_(n+4)) = bold(0). $ <eq-st-mathieu-5x5>
The inner $3 times 3$ block, indicated in the étude as $delta_3$, gives a quintic in $lambda$; the full $5 times 5$ block gives a $5 times 5$ secular determinant which Boyd resolves by a single Newton step. The *eigenvalue correction* $delta(q) eq.def lambda(q) - n^2$ for $n = 15$ comes out of the Galerkin diagonalisation directly.

=== Computational Étude 21.2: Sideband Truncation for $"ce"_(15)$ <etude-st-mathieu>

We compute the first $30$ Fourier coefficients of $"ce"_(15)$ at $q = 10$ on the full basis (high-$N$ Galerkin, $N_("modes") = 64$), reproducing the bar chart of Boyd Fig 19.1, and then sweep $q in [0, 25]$ and compare the eigenvalue correction $delta(q)$ as obtained from the high-$N$ reference, the $5 times 5$ sideband, and the inner $3 times 3$ sideband. A third panel breaks the trick deliberately by repeating the experiment at the small carrier $n = 3$, where $q \/ n^2$ is large and the sidebands extend down to the lowest mode, so truncation fails.

In Python:

```python
def galerkin_full(N_modes, q):
    """Build the (N_modes x N_modes) Galerkin matrix on the basis
    {cos(1 x), cos(3 x), cos(5 x), ...} for Mathieu's equation."""
    modes = odd_modes(N_modes)
    M = np.diag(modes.astype(float) ** 2)
    for i in range(N_modes):
        for j in range(N_modes):
            if abs(modes[i] - modes[j]) == 2:
                M[i, j] = q
    return M, modes


def sideband_eigenvalue(n_carrier, q, half_width):
    """Sideband truncation: pick (2 hw + 1) modes around n_carrier,
    diagonalise, return the eigenvalue closest to n_carrier^2."""
    modes = np.array([n_carrier + 2 * k
                      for k in range(-half_width, half_width + 1)], dtype=int)
    M = np.diag(modes.astype(float) ** 2)
    for i in range(len(modes)):
        for j in range(len(modes)):
            if abs(modes[i] - modes[j]) == 2:
                M[i, j] = q
    lam = np.linalg.eigvalsh(M)
    return lam[np.argmin(np.abs(lam - n_carrier ** 2))]
```

In MATLAB:

```matlab
function lam = sideband_eigenvalue(n_carrier, q, half_width)
    modes = arrayfun(@(k) n_carrier + 2 * k, -half_width:half_width);
    M = diag(double(modes).^2);
    for i = 1:numel(modes)
        for j = 1:numel(modes)
            if abs(modes(i) - modes(j)) == 2
                M(i, j) = q;
            end
        end
    end
    lams = sort(eig((M + M.') / 2));
    [~, idx] = min(abs(lams - n_carrier^2));
    lam = lams(idx);
end
```

In Julia:

```julia
function sideband_eigenvalue(n_carrier::Int, q::Real, half_width::Int)
    modes = [n_carrier + 2 * k for k in -half_width:half_width]
    M = Diagonal(Float64.(modes) .^ 2) |> Matrix
    for i in 1:length(modes), j in 1:length(modes)
        abs(modes[i] - modes[j]) == 2 && (M[i, j] = q)
    end
    lam = eigvals(Symmetric(M))
    return lam[argmin(abs.(lam .- n_carrier^2))]
end
```

#figure(
  image("../figures/ch21/python/mathieu_sideband.pdf", width: 100%),
  caption: [Étude 21.2: Mathieu sideband truncation (six panels in a 3$times$2 grid). Top row = coefficient pictures, bottom row = eigenvalue corrections. Left column = the easy regime ($"ce"_(15)$, $q \/ n^2 = 10 \/ 225$), right column = the hard regime ($"ce"_3$, $q \/ n^2 = 10 \/ 9$), middle column = mechanism diagnostics. _(a)_: $|a_n|$ for $"ce"_(15)$ at $q = 10$ --- carrier $|a_(15)| approx 0.97$, sidebands $|a_(13)|, |a_(17)| approx 0.17, 0.15$, and $|a_(11)|, |a_(19)| approx 0.017, 0.011$ --- a tight cluster around the carrier. _(b)_: cluster half-width versus $q$ at three carriers $n in {3, 7, 15}$; the dotted grey line at half-width 2 marks the capacity of a $5 times 5$ sideband box, and the curves show at what $q$ each carrier first exceeds that capacity. _(c)_: $|a_n|$ for $"ce"_3$ at the _same_ $q = 10$; the cluster has spread to all small $n$ --- the visual reason a $5 times 5$ box no longer captures it. _(d)_: $delta(q) eq.def lambda(q) - 225$ for $"ce"_(15)$, computed three ways; the $5 times 5$ sideband (teal dashed) is indistinguishable from the high-$N$ reference (navy) at the resolution of the plot, the $3 times 3$ (coral dotted) lags slightly. _(e)_: convergence of the truncated $delta$ to the full value as the sideband half-width grows from 1 to 10, at fixed $n = 15$, $q = 10$; the $5 times 5$ box already lies near machine precision. _(f)_: the same three-way comparison for $"ce"_3$ where $q \/ n^2$ is large --- the $5 times 5$ sideband is wildly wrong because the relevant cluster has spread all the way down to the lowest mode.],
) <fig-st-mathieu>

#etude-conclusion[
  At $q = 10$ the carrier coefficient is $|a_(15)| approx 0.972$, the dominant sidebands $|a_(13)|, |a_(17)| approx 0.17, 0.15$, and the next pair $|a_(11)|, |a_(19)| approx 0.017, 0.011$. The $5 times 5$ secular determinant gives $delta_5(10) approx 0.22547$ against the high-$N$ value $delta_("full")(10) approx 0.22561$ --- *four matching decimals from a $5 times 5$ matrix*. The $3 times 3$ already gives $delta_3(10) approx 0.21334$, two-decimal accuracy. The breakdown panel at $n = 3$ shows the trick is *local in spectral space*: when $q \/ n^2$ becomes order one the relevant sideband band reaches all the way down to the lowest mode and the truncation loses its locality. The principle that follows generalises this observation across all carrier-plus-envelope problems.
]

The three scripts are available in:
- `codes/python/ch21/mathieu_sideband.py`
- `codes/matlab/ch21/mathieu_sideband.m`
- `codes/julia/ch21/mathieu_sideband.jl`

#block(stroke: 0.8pt + rgb(20, 45, 110), radius: 3pt, inset: 10pt,
  fill: rgb(248, 250, 254))[
  *Principle (sideband truncation, generalised).*  Whenever a problem has a *fast carrier* and a *slow envelope* --- envelope solitons of the nonlinear Schrödinger family, sidebands of a parametric oscillator, perturbed eigenmodes of a high-order quantum well --- the right truncation is a band around the carrier, not a band at the bottom. The reduction in matrix size is from $O(n)$ to $O(1)$ --- *if* the sideband locality holds.
]

== Tailored Basis Functions for Pathologies <sec-st-tailored>

The lesson of @sec-st-sideband generalises: *match the basis to the structure of the solution*. When the structure is a high carrier, the match is a sideband band; when the structure is a *singularity*, the match is to *include the singular function in the basis*. We open this section with the simplest possible illustration: singularity subtraction. We close it with a deliberately constructed *cross-over truncation* that demonstrates why a pure asymptotic statement is not always a numerical statement.

=== Computational Étude 21.3: Singularity Subtraction <etude-st-subtract>

Take the test function
$ f(x) = e^x + 0.1 (1 + x)^(2 \/ 3),  quad x in [-1, 1]. $ <eq-st-subtract-f>
The first term is entire; the second has a *corner-type* singularity at $x = -1$ --- continuous but not smooth (its first derivative is unbounded). Expanded directly in Chebyshev polynomials, $f$ has coefficients that decay only algebraically (driven by the singular term), and a maximum pointwise error that decays correspondingly slowly. The trick is to subtract the singular part *analytically*, expand the smooth remainder in Chebyshev, and re-add the singular function when reconstructing.

In Python:

```python
def cheb_coeffs(v):
    N = len(v) - 1
    extended = np.concatenate([v, v[N - 1:0:-1]])
    A = np.real(np.fft.fft(extended)) / N
    A[0] *= 0.5; A[N] *= 0.5
    return A[:N + 1]

# naive: expand f directly
a_naive = cheb_coeffs(f_full(cgl(N)))
# trick: expand f - 0.1 (1+x)^(2/3), reconstruct adds the singular back
a_trick = cheb_coeffs(f_smooth(cgl(N)))
approx_trick = reconstruct(a_trick, x_eval) + f_singular(x_eval)
```

In MATLAB:

```matlab
function a = cheb_coeffs(v)
    N = numel(v) - 1;
    ext = [v(:); flipud(v(2:N).')]';
    A = real(fft(ext)) / N;
    A(1) = 0.5 * A(1);
    A(N+1) = 0.5 * A(N+1);
    a = A(1:N+1);
end
% naive
a_naive = cheb_coeffs(f_full(cgl(N)));
% trick: expand f - 0.1 (1+x)^(2/3), then add singular back
a_trick = cheb_coeffs(f_smooth(cgl(N)));
approx_trick = clenshaw(a_trick, x_eval) + f_singular(x_eval);
```

In Julia:

```julia
function cheb_coeffs(v::AbstractVector)
    N = length(v) - 1
    ext = vcat(v, reverse(v[2:N]))
    A = real.(fft(ext)) ./ N
    A[1] *= 0.5; A[N+1] *= 0.5
    return A[1:N+1]
end
# naive
a_naive = cheb_coeffs(f_full.(cgl(N)))
# trick
a_trick = cheb_coeffs(f_smooth.(cgl(N)))
approx_trick = reconstruct(a_trick, x_eval) .+ f_singular.(x_eval)
```

#figure(
  image("../figures/ch21/python/singularity_subtraction.pdf", width: 100%),
  caption: [Étude 21.3: singularity subtraction. *Left*: Chebyshev coefficient magnitudes at $N = 80$. The naive expansion (coral) decays algebraically, $|a_n|$ tracking the decay of the singular term. The trick (teal) reaches machine $epsilon$ at modest $n$. *Right*: maximum pointwise error vs $N$. Naive error is $approx 6 times 10^(-6)$ at $N = 256$; the trick reaches $approx 1.8 times 10^(-15)$ already at $N = 16$. *Same problem, same machinery, $approx 10$ decimals recovered* by acknowledging the singularity analytically.],
) <fig-st-subtract>

The three scripts are available in:
- `codes/python/ch21/singularity_subtraction.py`
- `codes/matlab/ch21/singularity_subtraction.m`
- `codes/julia/ch21/singularity_subtraction.jl`

#etude-conclusion[
  The trick is *structural*: the basis ${T_n}$ is bad at representing $(1 + x)^(2 \/ 3)$ no matter how large $N$ is, but the basis ${T_n} union {(1 + x)^(2 \/ 3)}$ is *excellent* --- ten decimals recovered at $N = 16$ where the naive expansion plateaus at $approx 6 times 10^(-6)$. In two-dimensional problems with corner singularities (the Stokes streamfunction in a driven cavity is the classical example), the local asymptotic analysis at each corner provides the singular functions to add to the basis, and the same exponential recovery results. *Acknowledging the singularity analytically* is the secret; the basis just has to know what the singular function looks like.
]

=== Computational Étude 21.4: Cross-over Truncation <etude-st-crossover>

A pure asymptotic statement is not always a numerical statement. Boyd's cartoon @Boyd2000 makes the point in a single line:
$ a_n tilde.op 10 thin e^(-n \/ 3) + (10^(-6)) / n^5. $ <eq-st-crossover>
Asymptotically the algebraic $1 \/ n^5$ wins, but it does not win until $n$ is past a *cross-over* $n_(times) approx 120$. For all the $N$ that any practical computation will ever reach, the *exponential* head dominates, even though it is the algebraic tail that has the larger asymptotic order. The 'asymptotic' answer is irrelevant.

In Python:

```python
n_axis = np.arange(1, 401)
head_geom = 10.0 * np.exp(-n_axis / 3.0)
tail_alg  = 1.0e-6 / n_axis ** 5
log_diff  = np.log10(tail_alg + 1e-300) - np.log10(head_geom + 1e-300)
n_cross   = float(n_axis[np.argmin(np.abs(log_diff))])
```

In MATLAB:

```matlab
n_axis = 1:400;
head_geom = 10.0 * exp(-n_axis / 3.0);
tail_alg  = 1.0e-6 ./ n_axis.^5;
log_diff  = log10(tail_alg + 1e-300) - log10(head_geom + 1e-300);
[~, idx]  = min(abs(log_diff));
n_cross   = double(n_axis(idx));
```

In Julia:

```julia
n_axis = collect(1:400)
head_geom = 10.0 .* exp.(-n_axis ./ 3.0)
tail_alg  = 1.0e-6 ./ n_axis .^ 5
log_diff  = log10.(tail_alg .+ 1e-300) .- log10.(head_geom .+ 1e-300)
n_cross   = float(n_axis[argmin(abs.(log_diff))])
```

#figure(
  image("../figures/ch21/python/crossover_truncation.pdf", width: 100%),
  caption: [Étude 21.4: cross-over truncation. *Left*: Boyd's cartoon: the exponential head $10 thin e^(-n / 3)$ (navy) and the algebraic tail $10^(-6) / n^5$ (coral) both contribute to $a_n$; the sum (purple) shows a clean slope change at $n_(times) approx 120$. *Right*: real Chebyshev coefficients of $f(x) = e^(-((x-0.3)/0.1)^2) + 10^(-7) (x+1)^(1\/3)$ on $[-1, 1]$. The smooth bulk gives geometric decay near the origin; the weak endpoint singularity gives algebraic decay in the far tail. The picture is not a cartoon: the same two-slope structure appears in real spectra.],
) <fig-st-crossover>

The three scripts are available in:
- `codes/python/ch21/crossover_truncation.py`
- `codes/matlab/ch21/crossover_truncation.m`
- `codes/julia/ch21/crossover_truncation.jl`

#etude-conclusion[
  The lesson is *operational*: when you see a coefficient plot with a clean slope change, *do not refine $N$ to reach the asymptotic regime*. Most likely you will hit roundoff before you reach it. Either accept the head-region accuracy, or remove the singularity analytically (as in @etude-st-subtract) so that only the smooth head remains. This generalises the textbook warning that "asymptotic" answers are not always "numerical" answers --- a $1 \/ n^5$ tail with prefactor $10^(-6)$ has crossover near $n approx 120$, far beyond any practical resolution.
]

== Radiation Basis Functions for Wave Scattering <sec-st-radiation>

Of the special tricks in this chapter, the *radiation basis* is the most striking demonstration of what 'tailor the basis to the asymptotics' actually means. Standard rational Chebyshev functions on the line decay algebraically to constants at $plus.minus infinity$ --- they cannot represent the *non-decaying* sinusoidal plane waves of an asymptotic wavefunction at all, no matter how large $N$ becomes. Boyd's repair is to *augment* the rational basis with two *radiation functions*

$ H(x) cos(k x), quad H(x) sin(k x), quad H(x) eq.def (1 + tanh x) \/ 2, $ <eq-st-radiation>

so that the unknown asymptotic constants are carried by the coefficients of these special functions. The smoothed step $H(x)$ vanishes as $x arrow.r -infinity$ and saturates to $1$ as $x arrow.r +infinity$, so the radiation functions contribute $0$ on one side of the line and the desired plane wave on the other.

The pedagogical setup is the one-dimensional Schrödinger scattering problem
$ psi_(x x) + [k^2 - V(x)] psi = 0, quad V(x) = -v thin sech^2(x), quad v = 1, $
with the asymptotic ansatz $psi tilde.op e^(i k x) + alpha thin e^(-i k x)$ as $x arrow.r -infinity$ and $psi tilde.op beta thin e^(i k x)$ as $x arrow.r +infinity$. The reflection coefficient $cal(R) eq.def |alpha|^2$ admits a closed form for this Pöschl--Teller potential @Boyd2000:
$ cal(R) = (1 + cos(2 pi sqrt(v + 1\/4))) / (cosh(2 pi k) + cos(2 pi sqrt(v + 1\/4))). $ <eq-st-reflection>
Boyd's algorithm (his Eq 19.20--19.28) decomposes $psi = (1 + alpha) C(x) + i (1 - alpha) S(x)$ where $C, S$ each solve an inhomogeneous linear equation with $V(x) cos(k x)$ or $V(x) sin(k x)$ as forcing; the asymptotic constants $gamma_1, gamma_2, sigma_1, sigma_2$ then determine $alpha$ via a $2 times 2$ matching system (his Eq 19.29).

=== Computational Étude 21.5: Schrödinger Scattering with the Radiation Basis <etude-st-radiation>

We reproduce Boyd's Table 19.1: $cal(R)(k)$ for $v = 1$ at $k in {0.3, 0.6, 0.9, dots, 3.0}$, using the rational Chebyshev parameter $ell = 2$ (Boyd's choice) and $N = 48$ collocation points on the grid $x_i = ell cot(pi (2i-1) \/ (2 N))$.

In Python:

```python
def solve_scattering(N, ell, k):
    x, t = collocation_grid(N, ell)
    Phi   = np.zeros((N, N))
    Phipp = np.zeros((N, N))
    Phi[:, :N - 2]   = basis_TB(N, t)
    Phipp[:, :N - 2] = basis_TB_double_prime(N, t, ell)
    Phi[:, N - 2]   = H(x) * np.cos(k * x)
    Phipp[:, N - 2] = (Hpp(x) * np.cos(k * x)
                       - 2.0 * k * Hp(x) * np.sin(k * x)
                       - k ** 2 * H(x) * np.cos(k * x))
    Phi[:, N - 1]   = H(x) * np.sin(k * x)
    Phipp[:, N - 1] = (Hpp(x) * np.sin(k * x)
                       + 2.0 * k * Hp(x) * np.cos(k * x)
                       - k ** 2 * H(x) * np.sin(k * x))
    M = Phipp + np.outer(k ** 2 - V(x), np.ones(N)) * Phi
    a = np.linalg.solve(M, V(x) * np.cos(k * x))    # C-tilde
    b = np.linalg.solve(M, V(x) * np.sin(k * x))    # S-tilde
    gamma1, gamma2 = a[N - 2] + 1.0, a[N - 1]
    sigma1, sigma2 = b[N - 2],       b[N - 1] + 1.0
    A = np.array([[gamma1 + sigma2, sigma1 - gamma2],
                  [gamma2 - sigma1, gamma1 + sigma2]])
    re, im = np.linalg.solve(A, [sigma2 - gamma1, -sigma1 - gamma2])
    return complex(re, im)
```

In MATLAB:

```matlab
function [alpha, drift] = solve_scattering(N, ell, k, v)
    i = 1:N;  t = pi * (2 * i - 1) / (2 * N);  x = ell ./ tan(t);
    Vx = -v ./ cosh(x).^2;
    Phi = zeros(N, N);  Phipp = zeros(N, N);
    n_idx = 0:N-3;
    for j = 1:numel(n_idx)
        nn = n_idx(j);
        Phi(:, j) = cos(nn * t).';
        s = sin(t).'; c = cos(t).';
        Phipp(:, j) = -(nn / ell^2) * s.^3 .* (nn * cos(nn * t).' .* s + 2 * sin(nn * t).' .* c);
    end
    Hx = 0.5 * (1 + tanh(x));  Hpx = 0.5 ./ cosh(x).^2;  Hppx = -tanh(x) ./ cosh(x).^2;
    Phi(:, N-1) = (Hx .* cos(k * x)).';
    Phipp(:, N-1) = (Hppx .* cos(k * x) - 2 * k * Hpx .* sin(k * x) ...
                     - k^2 * Hx .* cos(k * x)).';
    Phi(:, N) = (Hx .* sin(k * x)).';
    Phipp(:, N) = (Hppx .* sin(k * x) + 2 * k * Hpx .* cos(k * x) ...
                   - k^2 * Hx .* sin(k * x)).';
    M = Phipp + (k^2 - Vx.') .* Phi;
    a = M \ (Vx .* cos(k * x)).';   b = M \ (Vx .* sin(k * x)).';
    gamma1 = a(N-1) + 1;  gamma2 = a(N);
    sigma1 = b(N-1);      sigma2 = b(N) + 1;
    A = [gamma1 + sigma2, sigma1 - gamma2; gamma2 - sigma1, gamma1 + sigma2];
    re_im = A \ [sigma2 - gamma1; -sigma1 - gamma2];
    alpha = complex(re_im(1), re_im(2));   drift = sum(a(1:N-2));
end
```

In Julia:

```julia
function solve_scattering(N::Int, ell::Real, k::Real)
    x, t = collocation_grid(N, ell)
    Phi   = zeros(N, N);  Phipp = zeros(N, N)
    Phi[:, 1:N-2]   = basis_TB(N, t)
    Phipp[:, 1:N-2] = basis_TB_double_prime(N, t, ell)
    Phi[:, N-1]   = H.(x) .* cos.(k .* x)
    Phipp[:, N-1] = Hpp.(x) .* cos.(k .* x) .- 2 .* k .* Hp.(x) .* sin.(k .* x) .-
                    k^2 .* H.(x) .* cos.(k .* x)
    Phi[:, N]     = H.(x) .* sin.(k .* x)
    Phipp[:, N]   = Hpp.(x) .* sin.(k .* x) .+ 2 .* k .* Hp.(x) .* cos.(k .* x) .-
                    k^2 .* H.(x) .* sin.(k .* x)
    M = Phipp .+ (k^2 .- V.(x)) .* Phi
    a = M \ (V.(x) .* cos.(k .* x))
    b = M \ (V.(x) .* sin.(k .* x))
    gamma1, gamma2 = a[N-1] + 1.0, a[N]
    sigma1, sigma2 = b[N-1],       b[N] + 1.0
    A = [gamma1 + sigma2  sigma1 - gamma2; gamma2 - sigma1  gamma1 + sigma2]
    re_im = A \ [sigma2 - gamma1, -sigma1 - gamma2]
    return complex(re_im[1], re_im[2]), sum(a[1:N-2])
end
```

#figure(
  image("../figures/ch21/python/radiation_scattering.pdf", width: 100%),
  caption: [Étude 21.5: Schrödinger $sech^2$ scattering with the radiation basis. *Left*: numerical $cal(R)(k)$ (coral circles) overlaid on the closed form (navy line, @eq-st-reflection); the spectral computation reproduces seven decades of magnitude. *Right*: absolute error (navy squares) stays at the $10^(-8)$ to $10^(-9)$ level uniformly over $k$, even where $cal(R)$ itself becomes as small as $10^(-8)$. The 'asymptotic-constant drift' $|sum_n a_n|$ (teal triangles) is the residual unrepresentable constant inherited from $T B_n (plus infinity) = 1$; it is small because $V(x)$ decays exponentially.],
) <fig-st-radiation>

#etude-conclusion[
  The reproduction of Boyd Table 19.1 is striking. At $k = 0.3$ the numerical and exact reflection coefficients agree to eight significant digits; at $k = 3.0$ the reflection coefficient itself has dropped to $cal(R) approx 2.2 times 10^(-8)$ but the absolute error is still only $9 times 10^(-10)$. *Spectral accuracy is essential here*: any linearly-convergent method would have lost the exponentially-small reflection coefficient long before $k = 3$. The same idea --- *augment the basis with the asymptotic non-decaying mode* --- generalises beyond scattering: *weakly nonlocal solitary wave#idx("weakly nonlocal solitary wave")s* (the KdV soliton with a small oscillatory tail, the AEW family, certain kink--antikink pairs in $phi^4$ theory) all consist of a localised core plus a small persistent oscillation. In each case, write the trial function as standard rational basis *plus* one or two known non-decaying modes, and let the linear system find the asymptotic amplitudes for you.
]

The three scripts are available in:
- `codes/python/ch21/radiation_scattering.py`
- `codes/matlab/ch21/radiation_scattering.m`
- `codes/julia/ch21/radiation_scattering.jl`

== Spectral Root-finding and Nonlinear Eigenparameters <sec-st-rootfinding>

Spectral ideas are useful not only for solving differential equations, but also for solving the *algebraic* equations that spectral discretisations generate. We pick two Chebyshev-based root-finders here. The first --- *Lanczos economization#idx("Lanczos economization")* --- replaces the original (expensive) function by a cheap polynomial surrogate, then root-finds the surrogate. The second --- *Ioakimidis' non-iterative formula* --- expresses a single root on an interval as the *quotient of two finite Chebyshev quadratures*, with no Newton iteration anywhere.

=== Computational Étude 21.6: Lanczos Economization <etude-st-lanczos>

When each evaluation of the determinant $D(lambda) = det A(lambda)$ of an $M$-by-$M$ matrix is expensive (one LU per evaluation, $O(M^3)$ flops), repeated Newton iteration costs *many* LUs --- one per Newton step per root. The Lanczos remedy is to evaluate $D$ at $K$ Chebyshev nodes, build a degree-$(K - 1)$ Chebyshev interpolant, and find its roots in closed form via the colleague-matrix eigenvalue problem.

In real applications $M$ may be several hundred to a few thousand; the toy below uses $M = 20$ to keep the demonstration fast, but the *cost ratio* (LUs avoided) is what matters.

In Python:

```python
def roots_chebyshev_surrogate(K, a, b):
    j = np.arange(K)
    t = np.cos(np.pi * j / (K - 1))
    lam = 0.5 * (a + b) + 0.5 * (b - a) * t
    D_samples = np.array([expensive_D(L) for L in lam])     # K LU calls
    # Chebyshev coefficients via DCT-I:
    extended = np.concatenate([D_samples, D_samples[K - 2:0:-1]])
    A = np.real(np.fft.fft(extended)) / (K - 1)
    A[0] *= 0.5;  A[K - 1] *= 0.5
    coeffs = A[:K]
    # Colleague-matrix root-finding (cheap):
    n = K - 1
    C = np.zeros((n, n))
    for i in range(n - 1):
        C[i, i + 1] = 0.5;  C[i + 1, i] = 0.5
    C[0, 1] = 1.0
    last = -coeffs[:n] / (2.0 * coeffs[n])
    last[n - 2] += 0.5
    C[n - 1, :] = last
    eigs = np.linalg.eigvals(C)
    real_in = np.array([float(np.real(z)) for z in eigs
                        if abs(z.imag) < 1e-6 and -1.0 <= z.real <= 1.0])
    return 0.5 * (a + b) + 0.5 * (b - a) * np.sort(real_in)
```

In MATLAB:

```matlab
function [roots_out, lam, D_samples] = roots_chebyshev(K, a, b, fn)
    j = 0:K-1;
    t = cos(pi * j / (K - 1));
    lam = 0.5 * (a + b) + 0.5 * (b - a) * t;
    D_samples = arrayfun(fn, lam);
    ext = [D_samples, fliplr(D_samples(2:K-1))];
    A = real(fft(ext)) / (K - 1);
    A(1) = 0.5 * A(1); A(K) = 0.5 * A(K);
    coeffs = A(1:K);
    n = K - 1;  C = zeros(n);
    for i = 1:n-1
        C(i, i + 1) = 0.5;  C(i + 1, i) = 0.5;
    end
    C(1, 2) = 1.0;
    last = -coeffs(1:n) / (2.0 * coeffs(n + 1));
    last(n - 1) = last(n - 1) + 0.5;
    C(n, :) = last;
    eigs_C = eig(C);
    real_in = sort(real(eigs_C(abs(imag(eigs_C)) < 1e-6 & ...
                              -1.0 <= real(eigs_C) & real(eigs_C) <= 1.0)));
    roots_out = 0.5 * (a + b) + 0.5 * (b - a) * real_in;
end
```

In Julia:

```julia
function roots_chebyshev_surrogate(K::Int, a::Real, b::Real)
    j = collect(0:K-1)
    t = cos.(pi .* j ./ (K - 1))
    lam = 0.5 * (a + b) .+ 0.5 * (b - a) .* t
    D_samples = [expensive_D(L) for L in lam]
    ext = vcat(D_samples, reverse(D_samples[2:K-1]))
    A = real.(fft(ext)) ./ (K - 1)
    A[1] *= 0.5; A[K] *= 0.5
    coeffs = A[1:K]
    n = K - 1
    C = zeros(n, n)
    for i in 1:n-1
        C[i, i+1] = 0.5;  C[i+1, i] = 0.5
    end
    C[1, 2] = 1.0
    last_row = -coeffs[1:n] ./ (2.0 * coeffs[n+1])
    last_row[n-1] += 0.5
    C[n, :] = last_row
    eigs = eigvals(C)
    real_in_window = sort([real(z) for z in eigs
                          if abs(imag(z)) < 1e-6 && -1.0 <= real(z) <= 1.0])
    return 0.5 * (a + b) .+ 0.5 * (b - a) .* real_in_window
end
```

#figure(
  image("../figures/ch21/python/lanczos_economization.pdf", width: 100%),
  caption: [Étude 21.6: Lanczos economization for an expensive determinant. *Left*: $D(lambda)$ on $[0.5, 30]$ with the five roots detected by both methods. The Lanczos sample points (teal dots) are 17 Chebyshev-extrema nodes; the colleague-matrix eigenvalue problem on the resulting interpolant returns the root locations directly. *Right*: the cost story. The naive scan-and-bisect strategy uses $60 + 30 times 5 = 210$ expensive LU factorisations to locate the same roots; Lanczos uses $17$. Both methods recover the (perturbed) eigenvalues to $approx 10^(-11)$ root-position accuracy.],
) <fig-st-lanczos>

The three scripts are available in:
- `codes/python/ch21/lanczos_economization.py`
- `codes/matlab/ch21/lanczos_economization.m`
- `codes/julia/ch21/lanczos_economization.jl`

#etude-conclusion[
  The cost ratio is $approx 12 : 1$ here ($210$ expensive factorisations versus $17$), and the lesson scales: when function-evaluation cost grows like $M^3$ for a matrix of size $M$, *the savings from a Chebyshev surrogate are decisive*. The strategy is general: *replace expensive functions by cheap polynomial surrogates, then root-find or optimise the surrogate*. Boyd's scaffolding generalises the trick to two unknowns (Ioakimidis), with the caveat that the polynomial degree grows quickly with dimensionality. Lanczos economization is a workhorse trick whenever each evaluation of $f$ involves a hidden eigensolve, integral, or simulation.
]

=== Computational Étude 21.7: Ioakimidis' Non-Iterative Root <etude-st-ioakimidis>

Suppose $f(x)$ has a single simple root $rho$ on $[a, b]$. Boyd writes (after Ioakimidis, unpublished preprint, cited in @Boyd2000 chapter 19.6) that one can compute $rho$ *non-iteratively* as the quotient of two Chebyshev quadratures of $1 / f$:

$ rho = (sum_(j=0)^(2 N) " " (-1)^j x_j / f(x_j)) / (sum_(j=0)^(2 N) " " (-1)^j 1 / f(x_j)), quad x_j = (1/2) [(a+b) - (a-b) cos(pi j \/ (2 N))], $ <eq-st-ioakimidis>

with the doubly-primed sums denoting the Clenshaw-Curtis half-weighting at the endpoints. Convergence is *geometric* in $N$; no Newton iteration is required.

In Python:

```python
def ioakimidis_root(N, a, b):
    j = np.arange(2 * N + 1)
    half = 0.5 * ((a + b) - (a - b) * np.cos(j * np.pi / (2.0 * N)))
    fj = f_test(half)
    sign = (-1.0) ** j
    weight = np.ones_like(j, dtype=float)
    weight[0] = 0.5; weight[-1] = 0.5
    num = np.sum(weight * sign * half / fj)
    den = np.sum(weight * sign        / fj)
    return num / den
```

In MATLAB:

```matlab
function rho = ioakimidis(N, a, b)
    j = 0:(2*N);
    half = 0.5 * ((a + b) - (a - b) * cos(j * pi / (2 * N)));
    fj = f_test(half);
    sign_j = (-1).^j;
    weight = ones(size(j));
    weight(1) = 0.5; weight(end) = 0.5;
    num = sum(weight .* sign_j .* half ./ fj);
    den = sum(weight .* sign_j        ./ fj);
    rho = num / den;
end
```

In Julia:

```julia
function ioakimidis_root(N::Int, a::Real, b::Real)
    j = collect(0:2*N)
    half = 0.5 .* ((a + b) .- (a - b) .* cos.(j .* pi ./ (2.0 * N)))
    fj = f_test.(half)
    sgn = [(-1.0) ^ jj for jj in j]
    weight = ones(length(j))
    weight[1] = 0.5; weight[end] = 0.5
    num = sum(weight .* sgn .* half ./ fj)
    den = sum(weight .* sgn        ./ fj)
    return num / den
end
```

#figure(
  image("../figures/ch21/python/ioakimidis_root.pdf", width: 80%),
  caption: [Étude 21.7: Ioakimidis non-iterative root for $f(x) = sin(x - pi \/ 4) \/ sqrt(1 + 10 x^2)$ on $[-1, 1]$ (the unique root is $rho = pi \/ 4 approx 0.7853981633974483$). Geometric convergence: the Ioakimidis quotient (navy) reaches $|"err"| approx 10^(-12)$ at $N = 30$ ($61$ function evaluations); plain bisection on the same number of evaluations (coral) only reaches $|"err"| approx 10^(-9)$. The Ioakimidis formula#idx("Ioakimidis formula") has *no iteration anywhere* --- it is a single quotient of two finite sums.],
) <fig-st-ioakimidis>

The three scripts are available in:
- `codes/python/ch21/ioakimidis_root.py`
- `codes/matlab/ch21/ioakimidis_root.m`
- `codes/julia/ch21/ioakimidis_root.jl`

#etude-conclusion[
  The Ioakimidis trick generalises modestly (to a pair of equations in two unknowns) but does not survive into high dimensions. Its great pedagogical value is that it *turns quadrature into root-finding*, a structurally satisfying unification: the unique root of $f$ on $[a, b]$ is the centre-of-mass of $1 \/ f$ measured by Clenshaw--Curtis. The geometric convergence comes for free from the spectral-quadrature accuracy of the chapter @ch-quadrature; the *absence of any iteration* is a reminder that some classical "iterative" problems are secretly direct calculations in disguise.
]

== A Spectral Transform as a Research Tool: the Hilbert Transform <sec-st-transforms>

Quadrature itself is the subject of chapters @ch-quadrature and @ch-periodic-trap; we will not duplicate that material here. We pick a single transform-based étude --- the *Hilbert transform via Fourier multiplier* --- because it does not appear elsewhere in the book and is one of the most useful *one-line* tools a researcher will ever write.

The Hilbert transform of a real function $f$ is the principal-value integral
$ H{f}(y) = (1 / pi) "PV" integral_(-infinity)^(infinity) f(x) / (x - y) thin d x. $ <eq-st-hilbert>
For a *periodic* analytic $f$ with period $2 pi$, the Hilbert transform on the circle is *exactly* a Fourier multiplier:
$ H{f}(y) = -i sum_(k != 0) op("sgn")(k) thin c_k thin e^(i k y), quad f(x) = sum_k c_k thin e^(i k x). $ <eq-st-hilbert-fourier>
In code this is three lines: FFT, multiply by $-i thin op("sgn")(k)$, inverse FFT. Errors fall geometrically (in fact super-geometrically for analytic periodic test functions whose Fourier coefficients decay factorially).

=== Computational Étude 21.8: Hilbert Transform via Fourier Multiplier <etude-st-hilbert>

We test on $f(x) = e^(cos x)$, whose Fourier series is
$ e^(cos x) = I_0 (1) + 2 sum_(k = 1)^(infinity) I_k (1) cos(k x), $
so that $H{f}(y) = 2 sum_(k=1)^(infinity) I_k (1) sin(k y)$. The modified Bessel coefficients $I_k (1)$ decay factorially, $I_k (1) tilde.op (1 \/ 2)^k \/ k!$ as $k arrow.r infinity$.

In Python:

```python
def hilbert_via_fft(f_values):
    N = len(f_values)
    F = np.fft.fft(f_values)
    k = np.fft.fftfreq(N, d=1.0 / N)        # integer wavenumbers
    multiplier = -1j * np.sign(k)
    multiplier[0] = 0.0                     # zero-mode contribution vanishes
    return np.real(np.fft.ifft(multiplier * F))
```

In MATLAB:

```matlab
function Hf = hilbert_via_fft(fv)
    N = numel(fv);
    F = fft(fv);
    k = [0:N/2-1, -N/2:-1];                  % integer wavenumbers (fftfreq order)
    multiplier = -1i * sign(k);
    multiplier(1) = 0;
    Hf = real(ifft(multiplier .* F));
end
```

In Julia:

```julia
function hilbert_via_fft(f_values::AbstractVector)
    N = length(f_values)
    F = fft(f_values)
    k = fftfreq(N, 1.0) .* N
    multiplier = -1im .* sign.(k)
    multiplier[1] = 0.0
    return real.(ifft(multiplier .* F))
end
```

#figure(
  image("../figures/ch21/python/hilbert_fourier.pdf", width: 100%),
  caption: [Étude 21.8: Hilbert transform via Fourier multiplier. *Left*: $f(x) = e^(cos x)$ (navy) and its exact Hilbert transform $2 sum I_k (1) sin(k y)$ (teal); coral circles show the FFT-based numerical $H{f}_N (y_j)$ at $N = 32$ Fourier modes per period. *Right*: max error vs $N$ on a semilog axis. The error reaches $10^(-15)$ at $N = 32$ and is bounded above by $2 thin I_(N/2)(1)$, the tail of the Bessel-coefficient series (coral dashed). Convergence is *super-geometric* on this analytic periodic test.],
) <fig-st-hilbert>

The three scripts are available in:
- `codes/python/ch21/hilbert_fourier.py`
- `codes/matlab/ch21/hilbert_fourier.m`
- `codes/julia/ch21/hilbert_fourier.jl`

#etude-conclusion[
  The *FFT-multiplier* implementation of the Hilbert transform is *one of the most useful one-line tools a researcher will ever write* --- on the periodic interval, the entire algorithm is `H = ifft(-i * sgn(k) * fft(f))`, and the error reaches $10^(-15)$ at $N = 32$ Fourier modes. The same three lines apply, with caveats about non-periodic boundary contamination, to functions decaying exponentially on the line: pick a domain $[-L, L]$ large enough that $|f(L)|$ is below the desired tolerance, sample on $N$ points, apply the FFT multiplier. Convergence becomes algebraic in $L$ (governed by the slow $1 \/ y$ decay of the Hilbert transform), but the algorithm itself is the same one line. Many "principal-value" calculations in physics and signal processing reduce to this idiom.
]

== Symbolic Spectral Methods: When Numerical Rules Reverse <sec-st-symbolic>

Up to this point, every étude in this book has lived in floating-point arithmetic on grids of size $N >> 1$, and the algorithmic choices have been driven by the asymptotic regime of large $N$. In this section we step into a *small-$N$ symbolic* regime in which the standard advice of chapters 1--20 *reverses*. The contrast is summarised in @tab-st-symbolic-vs-numerical.

#figure(
  table(
    columns: (1.0fr, 1.0fr),
    align: (left + horizon, left + horizon),
    stroke: none,
    inset: (x: 8pt, y: 6pt),
    fill: (col, row) => if row == 0 { rgb(20, 45, 110) }
                       else if calc.rem(row, 2) == 0 { rgb(248, 250, 254) }
                       else { white },
    [#text(white, weight: "bold")[Symbolic, $N$ small ($N tilde.op O(3)$)]],
    [#text(white, weight: "bold")[Numerical, $N$ large ($N >> 1$)]],
    [Galerkin or methods of moments], [collocation / pseudospectral],
    [Legendre, Gegenbauer, or even powers of $x$], [Chebyshev polynomials],
    [polynomialise transcendental structure (turn $cos, sin, sech$ into polynomials by trigonometric identities or substitutions like $z = tanh x$)], [unnecessary],
    [rationalise floating-point numbers and irrational collocation points (replace $1 / sqrt(5)$ by $1/sqrt(5)$ symbolically, or by a rational approximation $4/9$)], [unnecessary],
    [exploit symmetry aggressively (parity, antisymmetry, decoupled modes)], [important but less decisive],
    [basis recombination (build the BC into the ansatz, e.g. $u = (1 - x^2) p(x)$)], [boundary bordering also acceptable],
  ),
  caption: [The six precepts in which the symbolic small-$N$ regime reverses the usual numerical advice (after @Boyd2000 chapter 20.2).],
) <tab-st-symbolic-vs-numerical>

Two cautions are worth attaching to the table immediately. First, *common sense outranks rules*: the table is a guide, not a commandment. Second, *integer swell* and *nonlinearity* are the two main practical dangers: a symbolic boundary-layer expansion of a polynomial of degree $40$ in $epsilon$ may have integer coefficients hundreds of digits long, even when the final answer is a polynomial of modest degree. Knowing when to switch from exact rationals to high-precision floating-point is part of the trade.

We illustrate two of the precepts in detail.

=== Computational Étude 21.9: Carrier--Pearson Boundary-layer BVP <etude-st-bl>

Carrier and Pearson's textbook example (@Boyd2000 Eq 20.3):
$ epsilon^2 thin u_(x x) - u = -1, quad u(-1) = u(1) = 0, $ <eq-st-carrier-pearson>
has the exact solution $u(x) = 1 - cosh(x \/ epsilon) \/ cosh(1 \/ epsilon)$. For small $epsilon$ the solution has a steep boundary layer of width $approx epsilon$ at each endpoint. The Galerkin precepts apply in concert: incorporate the homogeneous BC by writing $u(x) = (1 - x^2) p(x)$ (basis recombination); restrict $p$ to even powers (parity / symmetry); and use *unweighted moments* $\{1, x^2, x^4, x^6\}$ as test functions, *not* Legendre polynomials weighted by $(1 - x^2)$, because we do not want to dilute the boundary-layer contribution.

The four-term ansatz $u_4 (x) = (1 - x^2) (a_0 + a_2 x^2 + a_4 x^4 + a_6 x^6)$ produces a $4 times 4$ symbolic linear system whose entries are polynomials in $epsilon^2$. The result, in fully simplified rational form (matching @Boyd2000 Eq 20.8), is
$ a_0 = (3(43243200 epsilon^6 + 2162160 epsilon^4 + 27720 epsilon^2 + 31)) / D, $
$ a_2 = (33 (327600 epsilon^4 + 41)) / D, $
$ a_4 = (429 (840 epsilon^2 - 13)) / D, $
$ a_6 = 6435 / D, $ <eq-st-bl-coeffs>
where the common denominator
$ D = 128 (2027025 epsilon^8 + 945945 epsilon^6 + 51975 epsilon^4 + 630 epsilon^2 + 1) $
is a polynomial of degree $8$ in $epsilon$. This is a *closed-form* approximate solution: a four-term Galerkin approximation expressed as rational functions of the small parameter $epsilon$.

In Python (using SymPy):

```python
def galerkin_4_term():
    x, eps = sp.symbols("x epsilon", positive=True)
    a = sp.symbols("a_0 a_2 a_4 a_6")
    p = a[0] + a[1] * x ** 2 + a[2] * x ** 4 + a[3] * x ** 6
    u4 = (1 - x ** 2) * p
    R4 = eps ** 2 * sp.diff(u4, x, 2) - u4 + 1
    eqs = [sp.integrate(x ** (2 * j) * R4, (x, -1, 1)) for j in range(4)]
    sol = sp.solve(eqs, a, dict=True)[0]
    return u4.subs({ai: sp.together(sp.simplify(sol[ai])) for ai in a})
```

In MATLAB (Symbolic Math Toolbox):

```matlab
function [u4_subs, sol, x, eps_sym] = galerkin_4_term()
    syms x eps_sym real
    a0 = sym('a_0'); a2 = sym('a_2'); a4 = sym('a_4'); a6 = sym('a_6');
    p = a0 + a2 * x^2 + a4 * x^4 + a6 * x^6;
    u4 = (1 - x^2) * p;
    R4 = eps_sym^2 * diff(u4, x, 2) - u4 + 1;
    eqs = sym(zeros(4, 1));
    for j = 0:3
        eqs(j + 1) = int(x^(2 * j) * R4, x, -1, 1);
    end
    s = solve(eqs, [a0 a2 a4 a6]);
    sol = struct('a_0', s.a_0, 'a_2', s.a_2, 'a_4', s.a_4, 'a_6', s.a_6);
    u4_subs = subs(u4, [a0, a2, a4, a6], [s.a_0, s.a_2, s.a_4, s.a_6]);
end
```

In Julia (using SymPy.jl, the same Python-SymPy backend):

```julia
function galerkin_4_term()
    @syms x::real epsilon::real
    a = [SymPy.symbols("a_$(2k)") for k in 0:3]
    p = sum(a[k+1] * x ^ (2k) for k in 0:3)
    u = (1 - x^2) * p
    R = epsilon^2 * diff(u, x, 2) - u + 1
    eqs = [integrate(x^(2j) * R, (x, -1, 1)) for j in 0:3]
    sol = solve(eqs, a)
    return u, sol, x, epsilon, a
end
```

#figure(
  image("../figures/ch21/python/symbolic_boundary_layer.pdf", width: 100%),
  caption: [Étude 21.9: Carrier--Pearson boundary-layer BVP solved by four-term symbolic Galerkin. *Left*: solutions at $epsilon = 1 \/ 10$ and $epsilon = 1 \/ 20$. The four-term Galerkin (coral circles) tracks the cosh-layer profile (navy) at modest cost. *Right*: $L^infinity$ error vs $epsilon$ on a log-log axis. The closed-form rational solution achieves $approx 10^(-9)$ error at $epsilon = 1$, and the error grows monotonically as $epsilon$ shrinks --- the boundary layer is harder to resolve with only four moments. The numerical results reproduce the values of @Boyd2000 Table 20.3 across the entire $epsilon$ range.],
) <fig-st-bl>

The three scripts are available in:
- `codes/python/ch21/symbolic_boundary_layer.py`
- `codes/matlab/ch21/symbolic_boundary_layer.m`
- `codes/julia/ch21/symbolic_boundary_layer.jl`

#etude-conclusion[
  The closed-form rational solution from four-term symbolic Galerkin is *highly uniform in $x$* but *highly non-uniform in $epsilon$* --- a hallmark of singular-perturbation problems. *Symbolic spectral methods reveal such structure cleanly* in a way that pure numerics cannot: the honest answer is a rational function in $epsilon$ that you can examine analytically, asymptotically expand around any chosen value, and use as input to further analysis. This is the central reason for the symbolic--numerical reversal of @sec-st-symbolic: a *small-$N$ symbolic* answer often gives more scientific insight than a *large-$N$ numerical* answer to the same problem.
]

=== Computational Étude 21.10: Symbolic Quartic Oscillator on the Line <etude-st-quartic>

The eigenvalue problem
$ u_(y y) + (E - y^4) u = 0, quad |u| arrow.r 0 thin "as" thin |y| arrow.r infinity, $ <eq-st-quartic>
is the classical *quartic oscillator#idx("quartic oscillator")* of quantum mechanics. We map the line to $[-1, 1]$ via the rational-Chebyshev coordinate $y = ell x \/ sqrt(1 - x^2)$ with $ell = 2$ (a free *width* parameter that we tune to roughly match the scale of the bound states), build in the BC by writing $u(x) = (1 - x^2) (a_1 + a_2 x^2 + a_3 x^4 + a_4 x^6 + a_5 x^8)$, and demand orthogonality of the residual against the test functions $\{1, x^2, x^4, x^6, x^8\}$. The result is a $5 times 5$ linear system whose entries depend *linearly* on $E$. The vanishing of its determinant is therefore a *degree-5 polynomial in $E$* --- the *secular determinant* $D(E)$ --- whose roots approximate the lowest five even-parity eigenvalues.

After symbolic simplification (with $ell = 2$), the secular polynomial reproduces @Boyd2000 Eq 20.16 *to all printed digits*:
$ D(E) \/ D(0) = 1 - 1.143704 thin E + 0.203243 thin E^2 - 0.010199 thin E^3 + 1.235 times 10^(-4) thin E^4 - 1.532 times 10^(-7) thin E^5. $ <eq-st-quartic-secular>

In Python (using SymPy):

```python
def build_secular_determinant():
    x, E = sp.symbols("x E", real=True)
    ell = sp.Integer(2)
    a = sp.symbols("a1:6")
    u = (1 - x ** 2) * sum(a[k] * x ** (2 * k) for k in range(5))
    R = (1 - x ** 2) ** 4 * ((1 - x ** 2) * sp.diff(u, x, 2)
                              - 3 * x * sp.diff(u, x)) / ell ** 2 \
        + ((1 - x ** 2) ** 2 * E - ell ** 4 * x ** 4) * u
    eqs = [sp.integrate(x ** (2 * j) * R, (x, -1, 1)) for j in range(5)]
    M = sp.Matrix([[sp.simplify(eqs[i].coeff(a[j])) for j in range(5)]
                   for i in range(5)])
    return sp.Poly(sp.simplify(M.det()), E)
```

In MATLAB (Symbolic Math Toolbox):

```matlab
function [coeffs_E, M] = build_secular()
    syms x E real
    ell = sym(2);
    a = sym('a_', [1 5]);
    p = sum(a .* x.^(2 * (0:4)));
    u = (1 - x^2) * p;
    u_x = diff(u, x);  u_xx = diff(u, x, 2);
    R = (1 - x^2)^4 * ((1 - x^2) * u_xx - 3 * x * u_x) / ell^2 ...
        + ((1 - x^2)^2 * E - ell^4 * x^4) * u;
    eqs = sym(zeros(5, 1));
    for j = 0:4
        eqs(j + 1) = int(x^(2 * j) * R, x, -1, 1);
    end
    M = sym(zeros(5));
    for i = 1:5
        for j = 1:5
            M(i, j) = simplify(diff(eqs(i), a(j)));
        end
    end
    DE = simplify(det(M));
    coeffs_E = sym2poly(DE);
end
```

In Julia (using SymPy.jl):

```julia
function build_secular_determinant()
    @syms x::real E::real
    ell = SymPy.Sym(2)
    a = [SymPy.symbols("a_$(k)") for k in 1:5]
    u = (1 - x^2) * sum(a[k] * x^(2(k-1)) for k in 1:5)
    R = (1 - x^2)^4 * ((1 - x^2) * diff(u, x, 2) - 3 * x * diff(u, x)) / ell^2 +
        ((1 - x^2)^2 * E - ell^4 * x^4) * u
    eqs = [SymPy.integrate(x^(2j) * R, (x, -1, 1)) for j in 0:4]
    M = SymPy.Sym[simplify(eqs[i].coeff(a[j])) for i in 1:5, j in 1:5]
    Mmat = SymPy.Sym.(reshape(M, 5, 5))
    return sympy.Poly(simplify(symbolic_det(Mmat)), E), x, E
end
```

#figure(
  image("../figures/ch21/python/symbolic_quartic_oscillator.pdf", width: 100%),
  caption: [Étude 21.10: symbolic quartic oscillator on the line. *Left*: $log_(10) |D(E)|$ on the lower window $E in [0, 50]$, where the deep V-shaped dips are the zeros of the secular determinant. A rug strip beneath the curve marks the five Bender--Orszag references (teal circles) and the numerical roots inside the window (coral crosses) for the five lowest even-parity eigenvalues $E_n, thin n = 0, 2, 4, 6, 8$. The fifth numerical root is a mirage at $E approx 720$, signalled by the annotated arrow rather than allowed to stretch the axis. *Right*: relative error per mode index. The lower three modes ($E_0, E_2, E_4$) are accurate to $0.5%, 2.7%, 6.9%$ --- usable for a five-term symbolic approximation. The upper two modes ($E_6, E_8$) are *numerical mirages* of the truncation: $E_6$ has $142%$ relative error, $E_8$ has $1790%$ relative error, both wholly worthless. *Trust the lower spectrum first.*],
) <fig-st-quartic>

The three scripts are available in:
- `codes/python/ch21/symbolic_quartic_oscillator.py`
- `codes/matlab/ch21/symbolic_quartic_oscillator.m`
- `codes/julia/ch21/symbolic_quartic_oscillator.jl`

#etude-conclusion[
  Five-term Galerkin gives $E_0 = 1.0654$, $E_2 = 7.6551$, $E_4 = 17.380$ against the Bender--Orszag reference $1.060, 7.456, 16.262$ --- the lowest three are $0.5 %$, $2.7 %$, $6.9 %$ off, *usable for back-of-envelope work*, and the algebraic structure of the secular determinant $D(E)$ (a polynomial in $E$ with rational coefficients) allows further symbolic analysis. The upper two modes ($E_6, E_8$) are *spectral mirages*: $142 %$ and $1790 %$ off the true values. *As a rule of thumb, a five-term symbolic approximation gives reliable answers for roughly $N \/ 2 = 2$--$3$ of its eigenvalues; the rest are decorative.* The same $N \/ 2$ rule recurs throughout numerical eigenanalysis (cf.\ Étude 18.3 for the high-$N$ numerical version).
]

== A Short Appendix Inside the Chapter: the Tau-Method <sec-st-tau>

Lanczos invented the *tau-method* in 1938 (the same paper that gave the world the *pseudospectral* method). The tau-method is both an algorithm and a philosophy. As an *algorithm*, the tau-method is essentially Galerkin with boundary bordering using Chebyshev test functions; the formal definition is in Chapter 21 of @Boyd2000. As a *philosophy*, it is the more interesting object: instead of approximating the solution of the *exact* problem, one finds the *exact* solution of a *modified* problem in which the modification is a small perturbation chosen in a spectrally good basis.

This second strategy --- *to solve approximate equations exactly* --- is the philosophy of the tau-method. We illustrate it on a one-line first-order ODE.

=== Computational Étude 21.11: Tau Method on a First-Order ODE <etude-st-tau>

Consider the elementary problem (Boyd Eq 21.13)
$ u'(x) + u(x) = 0, quad u(-1) = 1, $ <eq-st-tau-original>
whose exact solution is $u(x) = e^(-(x + 1))$. A polynomial cannot satisfy this equation identically (its derivative would need to satisfy the equation, but the derivative of a degree-$N$ polynomial is degree $N - 1$, with one fewer constraint). Lanczos's idea: solve the *modified* equation
$ v'(x) + v(x) = tau thin T_N (x), quad v(-1) = 1, $ <eq-st-tau-modified>
which *does* have a unique polynomial solution of degree $N$. The amplitude $|tau|$ of the perturbation is the natural measure of how far we have departed from the original problem; for a smooth solution, $|tau|$ should decay geometrically with $N$.

In Python:

```python
def tau_solve(N):
    """v'(x) + v(x) = tau T_N(x), v(-1) = 1, on the (N+1)-pt CGL grid."""
    D, x = cheb_matrix(N)
    L = D + np.eye(N + 1)
    TN = np.cos(N * np.arccos(x))
    A = np.zeros((N + 2, N + 2))
    b = np.zeros(N + 2)
    A[0:N+1, 0:N+1] = L
    A[0:N+1, N+1]   = -TN
    A[N+1, :]   = 0.0
    A[N+1, N]   = 1.0     # v at x = -1 (last CGL node)
    b[N+1] = 1.0
    sol = np.linalg.solve(A, b)
    return x, sol[0:N+1], sol[N+1]
```

In MATLAB:

```matlab
function [x, v, tau] = tau_solve(N)
    [D, x] = cheb_matrix(N);
    L = D + eye(N + 1);
    TN = cos(N * acos(x)).';
    A = zeros(N + 2, N + 2);  b = zeros(N + 2, 1);
    A(1:N+1, 1:N+1) = L;
    A(1:N+1, N+2)   = -TN;
    A(N+2, :)   = 0;  A(N+2, N+1) = 1;  b(N+2) = 1;
    sol = A \ b;
    v = sol(1:N+1);  tau = sol(N+2);
end
```

In Julia:

```julia
function tau_solve(N::Int)
    D, x = cheb_matrix(N)
    L = D + I
    TN = cos.(N .* acos.(x))
    A = zeros(N + 2, N + 2);  b = zeros(N + 2)
    A[1:N+1, 1:N+1] = L
    A[1:N+1, N+2]   = -TN
    A[N+2, :] .= 0.0;  A[N+2, N+1] = 1.0;  b[N+2] = 1.0
    sol = A \ b
    return x, sol[1:N+1], sol[N+2]
end
```

#figure(
  image("../figures/ch21/python/tau_first_order.pdf", width: 100%),
  caption: [Étude 21.11: tau method#idx("tau method") on $u' + u = 0$, $u(-1) = 1$. *Left*: amplitude $|tau|$ of the Lanczos perturbation (navy circles) decays geometrically with $N$, reaching machine zero by $N = 16$. The pointwise error of the tau approximation $v_N (x_j)$ (teal squares) and of a standard pseudospectral solution (coral triangles) decay at the same rate, both reaching machine $epsilon$ by $N approx 12$. *Right*: at $N = 16$ the tau solution and the exact $u(x) = e^(-(x+1))$ are visually indistinguishable; the residual perturbation has amplitude $approx 5 times 10^(-16)$.],
) <fig-st-tau>

The three scripts are available in:
- `codes/python/ch21/tau_first_order.py`
- `codes/matlab/ch21/tau_first_order.m`
- `codes/julia/ch21/tau_first_order.jl`

#etude-conclusion[
  The pedagogical value lies in the *picture*: $|tau|$ falls geometrically because the right-hand side $tau T_N (x)$ that we have added to the equation is itself spectrally small in the Chebyshev norm. *The "modification" Lanczos invokes is, quantitatively, exactly as small as the truncation order allows.* The same construction generalises to higher-order ODEs, eigenvalue problems, and (with care about boundary conditions) to PDEs. The closing trick of the chapter expresses its philosophy: *to solve approximate equations exactly* is sometimes more revealing --- and more honest --- than to solve exact equations approximately.
]

== A Non-Exhaustive Literature Review <sec-st-lit>

The tricks gathered in this chapter are a small sampler of a much larger and rapidly growing toolbox. The 2020--2026 literature has been particularly active: the techniques classically called 'special tricks' have matured from computational workarounds into formalised frameworks that now underpin contemporary work in computational physics, signal processing, scientific machine learning, and applied mathematics. The recurring lesson of that recent literature, as Dr. Dutykh's own report on advanced structural interventions makes plain, is that *the most significant breakthroughs in spectral computing arise not from increasing floating-point precision or hardware throughput, but from structurally modifying the basis, the domain, or the operator itself to capture the underlying physics analytically*. The following pointers organise that literature thematically; they are intended as starting points for further reading, not as a complete survey.

#heading(level: 3, numbering: none)[General references]

The principal source of the tricks compiled in this chapter is John P. Boyd's encyclopaedic _Chebyshev and Fourier Spectral Methods_ @Boyd2000, whose chapters 19 (special tricks), 20 (symbolic spectral method#idx("symbolic spectral method")s), and 21 (the tau method) deserve close reading. Boyd's own bibliography in those chapters is itself a significant resource. The companion textbook by Trefethen @Trefethen2000 supplies the underlying Chebyshev / Fourier machinery in concise computational form. Together they are the two volumes that the present text most directly extends.

#heading(level: 3, numbering: none)[Domain scaling and truncation length]

The étude that opens this chapter (@etude-st-rescue) restates the modern view, articulated repeatedly in the literature since 2020, that the truncation length $L$ is not a cosmetic parameter but an *active, nonlinear degree of freedom* coupled to the physics @Boyd2000. Whenever the action is concentrated in a microscopic region of a macroscopic domain --- boundary layers in fluid mechanics, localised potential wells in quantum mechanics --- static uniform grids inevitably exhaust their resolution before reaching convergence. Modern solvers integrate active length-scale parameterisations directly into their basis functions, treating widths and scales as free variables to be tuned alongside the linear coefficients.

#heading(level: 3, numbering: none)[Sideband truncation and frequency-domain locality]

Boyd's archetypal sideband-truncation example for $"ce"_(15)$, reproduced in @etude-st-mathieu, has gained renewed life in nonlinear optics and condition monitoring. Between 2020 and 2026 the ANR OPTIMAL project @Polya2024Sideband, coordinated by the FEMTO-ST institute, has used sideband truncation as the operational tool for managing the highly noise-driven dynamics of four-wave mixing in optical fibres; truncation-induced spectral broadening in Airy-pulse-driven cascaded four-wave mixing has been investigated in detail by @Anastassiou2025. In a different application area --- electromechanical fault diagnosis --- the same philosophy underlies the *sparse Bayesian learning* (SBL) frameworks of @Pradeep2023 for broken-rotor-bar fault identification: by *truncating the spectral grid* to the symmetric pair of fault-sideband frequencies adjacent to the supply harmonic and feeding only that sub-grid to the Bayesian estimator, robust diagnosis becomes possible in the noisy short-sampling-time regimes where classical FFT analysis fails. The carrier--envelope generalisation also underlies the literature on envelope solitons; see @Boyd2000 chapter 19.5 on weakly nonlocal solitary waves.

#heading(level: 3, numbering: none)[Singularity subtraction, swapping, and cross-over truncation]

The 1-D singularity-subtraction surrogate of @etude-st-subtract sits at the elementary end of a substantial recent literature on the *analytical isolation* of pathological behaviour. In computational electromagnetics, high-order singularity-subtraction methods for Nyström discretisations of boundary integral equations have been developed by @PerezAranciba2018, while @Polimeridis2016 exploit the technique with high-order polynomial vector basis functions on planar triangles. A genuinely new family, the *singularity swapping methods* of @Klinteberg2024, treats nearly-singular layer-potential evaluations: instead of subtracting the singularity, the integrand is reformulated in the parameterisation space of the boundary curve via complex-plane root-finding, so that the singularity is geometrically swapped out of the integration. This delivers true spectral accuracy for arbitrarily close field evaluations on closed and piecewise-analytic curves. The cross-over truncation behaviour visualised in @etude-st-crossover --- the operational moral that the asymptotically dominant tail may be numerically irrelevant for any practical $N$ --- is discussed in this connection by @Boyd2000 chapter 19.3 and in many of the cited modern works.

#heading(level: 3, numbering: none)[Radiation basis functions and complex source beams]

Boyd's radiation-basis trick for 1-D Schrödinger scattering, reproduced in @etude-st-radiation, is the prototype of a much broader strategy: *augment the basis with the asymptotic non-decaying mode* whenever standard rational or polynomial bases are structurally incapable of representing it. In computational electromagnetics, the Eibert group has extensively developed this idea for three-dimensional antenna characterisation. Their Multi-Level Spherical Wave Expansion @Eibert2019NFFF and Complex Source Beam reformulations @EibertCSB2017 of the Fast Irregular Antenna Field Transformation Algorithm use radiation-type basis functions launched from the antenna under test as exact wave-packet solutions of the wave equation, exploiting the Gaussian-like field taper of complex source beams (after @Felsen1976CSB) to reach far-field accuracy with far fewer modes than generic spherical harmonics. The recent work of @NorouzipourEibert2025 examines the truncation-region question that the radiation-basis idea raises in practice. A different application of the same philosophy of *embedding radiation conditions directly into the operator* appears in the divide-and-conquer waveguide-circuit solvers of @Cassaro2025, which formulate exact radiation boundary conditions through mathematical projectors onto outgoing modes.

#heading(level: 3, numbering: none)[Spectral root-finding: Lanczos economization and Ioakimidis]

Spectral surrogates for expensive determinants (@etude-st-lanczos) and non-iterative Chebyshev quadrature for root-finding (@etude-st-ioakimidis) sit at the algebraic end of the spectral toolbox. The colleague-matrix algorithm exploited in our Lanczos-economization étude has been formalised and analysed in modern terms by @BrubeckNakatsukasaTrefethen2021 (the *Vandermonde with Arnoldi* approach to polynomial fitting and evaluation), and @HalikiasTownsend2024 provide rigorous forward rounding-error analyses that justify Chebyshev surrogates as numerically stable alternatives to classical root-finding even on flat or noisy determinant landscapes.

The Ioakimidis non-iterative root formula @AustinKravanjaTrefethen2014 has had the wider afterlife: its core idea --- locating zeros via quadrature in the complex plane --- is the mathematical foundation of modern contour-integral eigensolvers (the FEAST and Sakurai--Sugiura families), which use Cauchy-integral formulations to isolate eigenvalues of massive matrices inside specified contours. The original quadrature-rules application to plane equipotential lines and other harmonic-function curves is @Ioakimidis1996; an analogous closed-form integral formula for roots of nonlinear equations appears in @Ioakimidis2019Closedform. Recent work on interval quadratic equations @Hladik2022 uses Ioakimidis-style integral enclosures to compute exact, sharp bounds for real roots of parametric quadratics in interval arithmetic --- a domain notoriously hostile to iterative methods.

#heading(level: 3, numbering: none)[Hilbert transform and the analytic signal in scientific machine learning]

The Fourier-multiplier Hilbert-transform algorithm of @etude-st-hilbert was made operationally accessible to spectral-method audiences by Weideman and James (1992); its historical connection to the Hilbert--Huang transform (HHT) underlies a long line of high-resolution spectral-denoising work in fluid dynamics, climatology, and signal processing @HuangShenWind2010, including very recent applications to X-ray and photoelectron spectra @Spectral2022HHTFtest.

The truly striking development of 2025 is that this classical algorithm has migrated into deep learning. The *Hilbert Neural Operator* (HNO) of @HilbertNeuralOperator2025 explicitly addresses the stationarity and phase-blindness limitations of the Fourier Neural Operator (FNO @LiKovachki2021) by mapping the input PDE state to its analytic representation via the continuous Hilbert transform: amplitude and phase become independent learnable features, and the network's spectral convolution is performed in the analytic-signal domain. Empirically, HNO outperforms FNO on causal and phase-sensitive systems by margins large enough to suggest that *the most advanced AI architectures for PDEs remain deeply reliant on rigorous analytical restructuring* to overcome the representational bottlenecks that purely data-driven approaches inherit. The transition from FNO to HNO is the cleanest contemporary demonstration of a structural-alignment principle that a researcher reading this chapter should internalise as a habit of mind.

#heading(level: 3, numbering: none)[Symbolic spectral methods: hybrid pipelines]

The Carrier--Pearson example of @etude-st-bl and the symbolic quartic oscillator of @etude-st-quartic are minimal illustrations of a larger machinery that has been industrialised over the last decade. Dr. Davide Stocco's PhD thesis @Stocco2023 and the work of the Bertolazzi group at Trento @BertolazziBiralDaLio2006 develop systematic *symbolic--numeric pipelines* for solving differential--algebraic equations and optimal control problems in large multibody systems: the boundary-value problem from the variational form is generated *symbolically* (with all integration-by-parts steps and Jacobian assembly performed in computer algebra), and floating-point arithmetic is reserved for the final solve. Their most recent IMECE 2024 contribution @StoccoLarcherBertolazzi2024 reports automatic symbolic index reduction for DAE systems, enabling the real-time simulation of high-fidelity tire--ground enveloping models and complex vehicle suspension kinematics --- problems that are simply unreachable in pure floating-point because of integer swell in intermediate Jacobians. The Carrier--Pearson reference values for the boundary-layer test problem and the Bender--Orszag reference for the quartic oscillator (Bender and Orszag, _Advanced Mathematical Methods for Scientists and Engineers_, McGraw-Hill 1978, p 523) are the classical benchmarks; the integer-swell warning is in @Boyd2000 chapter 20.2.

#heading(level: 3, numbering: none)[The Lanczos $tau$-philosophy and time-delay systems]

Cornelius Lanczos's 1938 paper @Lanczos1938 is the founding document of the $tau$-method; the modern theoretical renaissance is due to Provoost and Michiels at Leuven, whose work between 2024 and 2026 has elevated the method from a specialised numerical trick to a foundational tenet of modern systems control theory. Their _Lanczos $tau$ framework for time-delay systems_ @Provoost2024 establishes three groundbreaking equivalences. First, when a shifted Legendre basis is used, the $tau$-method discretisation of a delay differential equation in the time domain is *strictly mathematically equivalent* to a high-order Padé rational approximation evaluated in the frequency domain --- bridging the long-standing gap between time-domain spectral convergence and frequency-domain transfer-function methods. Second, the reformulated $tau$-method gives rise to discretisations whose matrices are *highly sparse and self-nesting*, in stark contrast to the dense matrices of pseudospectral collocation on delay systems. Third, and most significantly, the $cal(H)_2$-norm of retarded delay differential equations --- a quantity central to robust controller synthesis whose previous spectral approximations converged only algebraically --- is shown to admit *supergeometric convergence* via the $tau$-method, with precise equivalence to pseudospectral collocation when the nonzero collocation points align with zeros of orthogonal polynomials. The follow-up paper @ProvoostMichielsH2norm2026 develops the algorithmic side of computing and optimising this $cal(H)_2$ norm. Canonical polynomials, the recursive shortcut Lanczos invoked for one-off problems with a single $tau$ term, are treated in @Boyd2000 chapter 21 and in the substantial subsequent literature pioneered by E. L. Ortiz and his collaborators (Boyd's Table 21.1 lists references spanning 1938--1997).

#heading(level: 3, numbering: none)[Cosmological applications: spectral methods in quantum gravity]

A particularly vivid contemporary application of the structural-alignment principle of this chapter --- and the étude on complex-plane detour in chapter 18 (@etude-complex-detour) --- is the spectral analysis of *quasinormal modes* of black holes and wormholes in quantum-gravity-inspired geometries. The recent work of @BaticDutykhScardigli2026 on _Planck stars_ (hypothetical black-hole geometries modified by scale-dependent gravity to remove the central singularity) deploys a highly tailored Chebyshev-collocation spectral method matched precisely to the effective Newtonian potential of the modified geometry, resolving for the first time the structurally complex *Martini-glass morphology* of the QNM spectrum: nearly equally-spaced overdamped modes with characteristic anomalous gaps, alongside isolated gravitational modes separated from the main sequence by massive frequency intervals. Standard low-order WKB approximations are simply incapable of resolving these features. A companion line of work @DutykhBaticWormholes2025 applies the same spectral methodology to noncommutative geometry-inspired wormholes, exploiting compactified spatial domains and high-order Chebyshev mappings to definitively establish the absence of overdamped modes in that geometry. These are computational successes that depend, at every step, on tailoring the discretisation exactly to the structural topology of the problem domain --- the same lesson the present chapter has tried to teach in miniature.

#heading(level: 3, numbering: none)[Three macroscopic trends]

A glance across this literature reveals three trends that the researcher leaving this chapter should keep in mind. First, *the dominance of structural embeddings*: rather than integrating through a pathology by grid refinement, modern algorithms enforce the physics through explicit basis selection. Second, *the evolution of transform frameworks in deep learning*: the leap from the Fourier Neural Operator to the Hilbert Neural Operator demonstrates that classical spectral insights remain indispensable even at the frontier of operator-learning architectures. Third, *the rigorous unification of continuous and discrete paradigms*: the precise equivalences proved between the Lanczos $tau$-method, pseudospectral collocation, and Padé rational approximation point to a deeper structural unity in numerical analysis that is only now being mapped out. The efficacy of any spectral method, in the end, is dictated unconditionally by how intimately its discretisation respects the geometry, symmetries, and asymptotes of the underlying problem; *optimal spectral intervention is never a shortcut around analytical understanding --- it is, profoundly and exclusively, mathematical understanding made fully operational*.

== A Farewell to the Student-Researcher <sec-st-farewell>

We close with a six-point checklist condensed from the entire chapter:

+ *Identify the dominant structure of the solution* before you choose a basis. A boundary layer, a corner singularity, a fast carrier modulated by a slow envelope, a non-decaying oscillatory tail, a discrete spectrum of nearby eigenvalues: each of these calls for a different match between basis and problem.

+ *Choose the basis to respect that structure*. Standard Chebyshev for smooth functions on $[-1, 1]$; rational Chebyshev for algebraic decay on the line; Hermite for Gaussian decay; sideband truncation around a known carrier; basis augmentation with singular functions when corners or endpoint singularities are unavoidable; radiation functions when the asymptote is a non-decaying plane wave.

+ *Choose the truncation to respect where the coefficients actually live*. The dominant Fourier coefficient may not be at $n = 0$. The slope of the coefficient plot may change near a cross-over $n_(times)$. The truncation length $L$ on a finite-interval surrogate is part of the discretisation, not a cosmetic detail.

+ *Use plots of coefficients, residuals, and errors as diagnostic instruments*. A coefficient plot tells you whether the basis is right. A residual plot tells you whether the equation is being satisfied. A pointwise-error plot tells you where the worst error lives. Treat them as instruments, not as final outputs.

+ *When symbolic work is appropriate, reverse the usual numerical instincts*. Galerkin over collocation. Legendre or Gegenbauer over Chebyshev. Polynomialise transcendental structure. Rationalise irrational constants. Exploit symmetry aggressively. Build the boundary conditions into the basis. Trust the lower spectrum; distrust the upper.

+ *When the standard approach fails, reformulate boldly but verify carefully*. Each of the eleven études of this chapter rests on an analytical insight that *changes the formulation*: a contour deformation, a basis augmentation, a coordinate map, a low-$N$ symbolic Galerkin, a $tau$ perturbation. The verification is always the same: a coefficient plot, a residual norm, and a comparison against an independently computed reference. *Reformulate boldly but verify carefully* --- this is the working method of the spectral researcher.

The central lesson of spectral methods is not that one basis fits all problems. The central lesson is that *analysis and computation must speak to each other*. When a computation struggles, ask what feature of the solution has not yet been honoured: its symmetry, its scale, its singularity, its oscillatory tail, its geometry, or its asymptotics. The right trick is not a shortcut around understanding; it is *understanding made operational*.

We wish you good luck with your spectral methods, and all the best in applying them to the problems that you, in turn, will choose to investigate.

== Exercises <sec-st-exercises>

The exercises below revisit the eleven études of this bonus chapter. The conceptual problems prove the structural facts on which each trick rests; the computational problems extend the études with convergence studies, parameter sweeps, and deliberate breakdowns; and the projects open onto the recent literature gathered in @sec-st-lit. Every computational task may be carried out in any of the book's three languages, and the named scripts under `codes/python/ch21`, `codes/matlab/ch21`, and `codes/julia/ch21` give a starting point.

=== Conceptual Exercises

#exercise(title: [Tridiagonal Sideband Coupling])[
  Mathieu's equation @eq-st-mathieu is discretised in the Fourier cosine basis. (a) Using the product identity $cos(2 x) cos(n x) = (1\/2) [cos((n + 2) x) + cos((n - 2) x)]$, prove that the perturbation $-2 q cos(2 x)$ couples the mode $cos(n x)$ only to its neighbours $cos((n plus.minus 2) x)$, so that the Galerkin matrix is tridiagonal within each parity class. (b) Derive the entries of the secular matrix @eq-st-mathieu-5x5 directly from this coupling and identify the role of the carrier index $n$. (c) Argue that the four parity classes decouple, and explain in one sentence why $"ce"_(15)$ lives in the class generated by ${cos(x), cos(3 x), cos(5 x), dots}$.
] <ex-st-sideband-tridiag>

#hint-for(<ex-st-sideband-tridiag>)[The cosine product identity converts each multiplication by $cos(2 x)$ into a shift of $plus.minus 2$ in the mode index; collect the coefficient of $cos(n x)$ in the residual to read off the $n$-th row of the matrix.]

#exercise(title: [Algebraic Decay and Singularity Subtraction])[
  Consider the singular term $g(x) = (1 + x)^(2\/3)$ of @eq-st-subtract-f on $[-1, 1]$. (a) Explain why the Chebyshev coefficients $a_n$ of a function carrying an endpoint branch point $(1 + x)^gamma$, with $gamma$ not a non-negative integer, decay only algebraically, $|a_n| tilde.op C thin n^(-(2 gamma + 1))$; evaluate the exponent for $gamma = 2\/3$. (b) Deduce the algebraic rate at which the maximum pointwise error of the naive expansion falls with $N$. (c) Explain why expanding the smooth remainder $f - 0.1 thin g$ and re-adding $0.1 thin g$ on reconstruction restores geometric decay, and state precisely what analytic knowledge the trick requires.
] <ex-st-algebraic-decay>

#exercise(title: [Where the Cross-over Lives])[
  The cartoon spectrum @eq-st-crossover superposes a geometric head $10 thin e^(-n\/3)$ and an algebraic tail $10^(-6)\/n^5$. (a) Write the condition that equates the two contributions and reduce it to $n - 15 ln n = 21 ln 10 approx 48.3$. (b) Solve it by hand iteration or graphically to confirm the cross-over index $n_(times) approx 120$. (c) A typical computation is truncated near $N approx 40$; argue that on $1 lt.eq.slant n lt.eq.slant N$ the geometric head dominates everywhere, so the asymptotically larger algebraic tail is numerically invisible, the operational moral of @etude-st-crossover.
] <ex-st-crossover-index>

#exercise(title: [Radiation Functions and the Smoothed Step])[
  The radiation basis @eq-st-radiation augments the rational Chebyshev set with $H(x) cos(k x)$ and $H(x) sin(k x)$, where $H(x) eq.def (1 + tanh x)\/2$. (a) Show that $H(x) arrow.r 0$ as $x arrow.r -infinity$ and $H(x) arrow.r 1$ as $x arrow.r +infinity$, and compute $H' (x)$ and $H'' (x)$ in closed form, confirming that both decay like $sech^2 x$. (b) Apply the Helmholtz operator $partial_(x x) + k^2$ to $H(x) cos(k x)$ and show that the result equals $H'' (x) cos(k x) - 2 k thin H' (x) sin(k x)$, exponentially localised near the origin because the non-decaying plane wave is annihilated. (c) Explain why this localisation is exactly the structural property that lets the linear solve of @etude-st-radiation carry the unknown asymptotic amplitudes in the coefficients of the two radiation functions.
] <ex-st-radiation-step>

#exercise(title: [Reflectionless Potentials])[
  The closed-form reflection coefficient @eq-st-reflection of the Pöschl--Teller well $V(x) = -v thin sech^2 x$ is $cal(R) = (1 + cos(2 pi sqrt(v + 1\/4)))\/(cosh(2 pi k) + cos(2 pi sqrt(v + 1\/4)))$. (a) Show that $cal(R) arrow.r 0$ as $k arrow.r infinity$ and identify the exponential decay rate set by $cosh(2 pi k)$. (b) Prove that the well is reflectionless, $cal(R) equiv 0$ for every $k > 0$, exactly when $v = ell (ell + 1)$ for a non-negative integer $ell$, and list the three smallest positive such $v$. (c) Explain why @fig-st-radiation requires a logarithmic vertical axis to display $cal(R)(k)$ for $v = 1$ across the tested range of $k$.
] <ex-st-reflectionless>

#exercise(title: [Colleague-Matrix Root-finding])[
  The Lanczos surrogate of @etude-st-lanczos finds the roots of a Chebyshev interpolant $p(x) = sum_(n = 0)^(K - 1) c_n thin T_n (x)$ as eigenvalues of a colleague matrix. (a) Using the three-term recurrence $2 x thin T_n (x) = T_(n + 1)(x) + T_(n - 1)(x)$, show that multiplication by $x$ acts on the coefficient vector through an almost-tridiagonal matrix, and identify the rank-one correction carried by the leading coefficient $c_(K - 1)$. (b) Conclude that the roots of $p$ in $[-1, 1]$ are the real eigenvalues of the resulting $(K - 1) times (K - 1)$ colleague matrix that lie in $[-1, 1]$. (c) Explain why this is numerically preferable to forming the companion matrix of $p$ in the monomial basis, with reference to the conditioning of the two bases.
] <ex-st-colleague>

#exercise(title: [Ioakimidis as a Centre of Mass])[
  Suppose $f$ has a single simple root $rho$ in $(a, b)$ and no other zero there. The Ioakimidis formula @eq-st-ioakimidis expresses $rho$ as a quotient of two Clenshaw--Curtis sums of $1\/f$. (a) Show that near $rho$ the integrand $1\/f(x)$ has a simple pole with residue $1\/f' (rho)$, so both sums are dominated by the same pole. (b) Interpret the quotient as the centre of mass $rho = (integral_a^b x thin mu(x) thin d x)\/(integral_a^b mu(x) thin d x)$ of the density $mu = 1\/f$, and explain why the shared residue cancels to leave $rho$. (c) State the role of the half-weighting of the endpoint nodes (the doubly-primed sums) and explain why the convergence is geometric rather than algebraic, citing the spectral quadrature of @ch-quadrature.
] <ex-st-ioakimidis-com>

#hint-for(<ex-st-ioakimidis-com>)[Linearise $f(x) approx f' (rho)(x - rho)$ near the root; then $1\/f tilde.op 1\/(f' (rho)(x - rho))$, and the common factor $1\/f' (rho)$ appears in both the numerator and denominator sums, cancelling in the ratio to leave the pole location $rho$.]

#exercise(title: [Hilbert Transform as a Fourier Multiplier])[
  On the circle the Hilbert transform @eq-st-hilbert acts as the Fourier multiplier @eq-st-hilbert-fourier, sending $e^(i k x)$ to $-i thin op("sgn")(k) thin e^(i k x)$ and annihilating the mean $k = 0$. (a) Deduce its action on the real harmonics, $cos(k x) arrow.r sin(k x)$ and $sin(k x) arrow.r -cos(k x)$ for $k > 0$, and verify it on $e^(cos x) = I_0 (1) + 2 sum_(k gt.eq.slant 1) I_k (1) cos(k x)$ to recover $H{f}(y) = 2 sum_(k gt.eq.slant 1) I_k (1) sin(k y)$. (b) Show that $H$ is skew-adjoint and satisfies $H^2 = -(I - P_0)$, where $P_0$ projects onto the mean. (c) Explain why analytic periodic data give the super-geometric convergence seen in @etude-st-hilbert, whereas a function analytic only in a strip of half-width $delta$ about the real axis yields geometric convergence at rate $approx e^(-delta N\/2)$.
] <ex-st-hilbert-multiplier>

=== Computational Exercises

#exercise(title: [Optimal Truncation Length])[
  Extend the rescue story of @etude-st-rescue with `rescue_naive_vs_tailored`. (a) For $f(y) = sech(y)$ at fixed $N = 48$, sweep the truncation half-length $L in [4, 40]$ and plot the maximum inner-window error against $L$ on a semilogarithmic axis; locate the optimum and relate it to the $L$ at which $sech(L)$ first sits at truncation-noise level. (b) Repeat for a narrow pulse $sech(3 y)$ and a broad pulse $sech(y\/3)$, and confirm that the optimal $L$ tracks the support width of the pulse. (c) At the optimal $L$, plot the Chebyshev coefficients and confirm geometric decay; at the defensive choice $L = 30$, show that the decay has stalled.
] <ex-st-optimal-length>

#exercise(title: [Sideband Convergence and Breakdown])[
  Using `mathieu_sideband` and @etude-st-mathieu. (a) Reproduce the coefficient bar chart of $"ce"_(15)$ at $q = 10$ and confirm the carrier-plus-sideband cluster of @fig-st-mathieu. (b) At fixed $n = 15$ and $q = 10$, grow the sideband half-width $m$ from $1$ to $10$ and plot $|delta_m (q) - delta_("full")(q)|$, verifying the geometric approach to the high-$N$ reference. (c) For carriers $n in {3, 7, 15}$, sweep $q$ and record the value at which the $5 times 5$ box first fails to deliver three correct digits of $delta$; tabulate this threshold against $q\/n^2$ and confirm that breakdown sets in near $q\/n^2 = O(1)$.
] <ex-st-sideband-sweep>

#exercise(title: [Singularity Subtraction across Exponents])[
  Using `singularity_subtraction` and @etude-st-subtract. (a) Reproduce the coefficient-decay and error-versus-$N$ panels for @eq-st-subtract-f. (b) Replace the singular term by $(1 + x)^gamma$ for $gamma in {1\/2, 2\/3, 4\/3, 7\/3}$ and measure the algebraic decay exponent of the naive Chebyshev coefficients, comparing each against the predicted $n^(-(2 gamma + 1))$. (c) Confirm that subtracting the matching singular function restores geometric convergence in every case, and report the smallest $N$ at which each tricked expansion reaches a maximum error of $10^(-12)$.
] <ex-st-subtract-exponents>

#exercise(title: [Radiation Basis Parameter Study])[
  Using `radiation_scattering` and @etude-st-radiation. (a) Reproduce Boyd's reflection-coefficient table for $v = 1$ at $k in {0.3, 0.6, dots, 3.0}$ and overlay the closed form @eq-st-reflection. (b) At fixed $k = 0.6$, study the absolute error as a function of the rational-Chebyshev map parameter $ell$ and of $N$, and identify the $(ell, N)$ region of best accuracy. (c) Set $v = 2$, one of the reflectionless values $v = m (m + 1)$ with integer $m$, and confirm that the computed $cal(R)(k)$ collapses toward machine zero across the whole $k$ range while you monitor the asymptotic-constant drift of @fig-st-radiation.
] <ex-st-radiation-study>

#exercise(title: [Lanczos Economization versus Ioakimidis])[
  Combine `lanczos_economization` (@etude-st-lanczos) and `ioakimidis_root` (@etude-st-ioakimidis). (a) For the expensive determinant $D(lambda)$ of @etude-st-lanczos, reproduce the colleague-matrix root locations and study the root accuracy as the number of Chebyshev sample nodes $K$ grows. (b) Tabulate the number of expensive LU factorisations spent by the surrogate against a scan-and-bisect baseline for matrix sizes $M in {20, 100, 500}$, and confirm that the cost saving scales with the per-evaluation $O(M^3)$ work. (c) For a single simple root on an interval, compare the Ioakimidis quotient @eq-st-ioakimidis with the surrogate by plotting error against the number of function evaluations for both, and comment on which is preferable when several roots are wanted at once.
] <ex-st-lanczos-vs-ioak>

#exercise(title: [Hilbert Transform and the Analyticity Strip])[
  Using `hilbert_fourier` and @etude-st-hilbert. (a) Reproduce the $f(x) = e^(cos x)$ test and confirm super-geometric convergence, with the error bounded by the Bessel-coefficient tail $2 thin I_(N\/2)(1)$ of @fig-st-hilbert. (b) Apply the transform to the Poisson kernel $f(x) = (1 - r^2)\/(1 - 2 r cos x + r^2)$, analytic in a strip whose width shrinks as $r arrow.r 1^-$, for several $r in (0, 1)$, and verify geometric convergence at a rate that degrades with $r$. (c) Verify the operator identity $H^2 = -(I - P_0)$ numerically by applying the FFT multiplier twice and recovering the negative of the mean-subtracted signal to machine precision.
] <ex-st-hilbert-strip>

=== Project-Style Exercises

#exercise(title: [Symbolic Spectral Pipelines and Integer Swell])[
  The symbolic études @etude-st-bl and @etude-st-quartic are minimal instances of the symbolic--numeric pipelines now used for large differential--algebraic systems. (a) Push the Carrier--Pearson Galerkin construction of @eq-st-bl-coeffs to six, eight, and ten even moments with `symbolic_boundary_layer`, and track the digit-length of the integer coefficients in the common denominator as the order grows, the integer-swell phenomenon flagged in @sec-st-symbolic. (b) Identify the order at which exact-rational arithmetic should give way to high-precision floating point, and quantify the accuracy sacrificed by the switch. (c) For the quartic oscillator of @etude-st-quartic, raise the number of Galerkin moments and test the $N\/2$ reliability rule against the Bender--Orszag eigenvalues, reporting how many eigenvalues stay trustworthy at each order. The industrial-scale version of this pipeline is developed by @Stocco2023, @StoccoLarcherBertolazzi2024, and @BertolazziBiralDaLio2006.
] <ex-st-symbolic-pipeline>

#hint-for(<ex-st-symbolic-pipeline>)[Watch the bit-length of the common denominator @eq-st-bl-coeffs: it roughly doubles with each added even moment, so a ten-moment solve already carries hundred-digit integers. Switch to fifty-digit floating point once the exact-rational `simplify` step dominates the runtime.]

#exercise(title: [The Tau-Method for Time-Delay Systems])[
  The tau philosophy of @sec-st-tau, to find the exact solution of a nearby problem, has become a foundational tool in modern control theory. (a) Extend `tau_first_order` (@etude-st-tau) from the scalar ODE @eq-st-tau-original to a second-order boundary-value problem and to a linear eigenvalue problem, tracking the perturbation amplitude $|tau|$ and confirming its geometric decay with $N$. (b) Implement the shifted-Legendre tau discretisation of a scalar delay differential equation and reproduce numerically the equivalence with a Padé approximation of the delay transfer function established by @Provoost2024. (c) Compute the $cal(H)_2$ norm of a retarded delay system and observe the supergeometric convergence reported by @ProvoostMichielsH2norm2026, contrasting the sparsity of the tau matrices with the dense pseudospectral-collocation alternative. The founding document is Lanczos's 1938 paper @Lanczos1938.
] <ex-st-tau-delay>

#hint-for(<ex-st-tau-delay>)[For part (a), border the differentiation operator with the boundary rows exactly as in `tau_first_order`, but append one $tau_j thin T_(N - j)(x)$ forcing column for each boundary constraint you must satisfy; solve for the vector of $tau_j$ alongside the spectral coefficients.]

#exercise(title: [From Fourier to Hilbert Neural Operators])[
  The one-line Hilbert transform of @etude-st-hilbert has lately migrated into operator-learning architectures. (a) Implement the Fourier-multiplier transform @eq-st-hilbert-fourier as a differentiable layer in an automatic-differentiation framework, and check its gradient against a finite-difference reference. (b) Form the analytic signal $f + i thin H{f}$ and use it to separate instantaneous amplitude and phase on a phase-shifted travelling-wave test, the representation that the Fourier Neural Operator @LiKovachki2021 cannot express. (c) On a small causal or phase-sensitive PDE, reproduce the qualitative advantage of the Hilbert Neural Operator @HilbertNeuralOperator2025 over a matched Fourier Neural Operator, and relate the analytic-signal idea to the Hilbert--Huang line of work @HuangShenWind2010. Comment on what the experiment says about the place of analytical restructuring at the frontier of scientific machine learning.
] <ex-st-hilbert-operator>

#hint-for(<ex-st-hilbert-operator>)[The multiplier $-i thin op("sgn")(k)$ is a fixed diagonal in Fourier space, so the layer is just `ifft(m .* fft(x))` and autodiff differentiates the transform for free; the only care needed is to zero the $k = 0$ and Nyquist modes so the layer stays a real linear map.]

#exercise(title: [Tailored Spectral Collocation for Quasinormal Modes])[
  The chapter's guiding maxim, match the discretisation to the structure, reaches the research frontier in the spectral computation of black-hole and wormhole quasinormal modes. (a) Set up a Chebyshev-collocation eigenproblem for the quasinormal modes of a model effective potential on a compactified radial domain, imposing outgoing-wave (radiation) boundary conditions in the spirit of @etude-st-radiation. (b) Where the dispersion relation appears as a transcendental determinant, locate its roots with the colleague-matrix surrogate of @etude-st-lanczos in place of Newton iteration, and use the complex-plane detour of @etude-complex-detour to navigate the branch structure. (c) Reproduce the qualitative spectral morphology, a nearly equally spaced overdamped sequence punctuated by anomalous gaps, reported for Planck stars by @BaticDutykhScardigli2026, and contrast it with the wormhole spectra of @DutykhBaticWormholes2025; explain why low-order WKB cannot resolve these features.
] <ex-st-qnm>

#hint-for(<ex-st-qnm>)[Compactify the radial coordinate so that spatial infinity maps to a finite endpoint, factor the known oscillatory asymptotics $e^(i omega r_*)$ out of the wavefunction so the remaining envelope is smooth enough to collocate, and read the quasinormal frequencies off as the generalised eigenvalues of the resulting matrix pencil.]
