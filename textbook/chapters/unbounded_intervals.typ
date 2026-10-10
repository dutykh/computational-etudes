// textbook/chapters/unbounded_intervals.typ
// Chapter 20: Spectral Methods on Unbounded Intervals
// Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
// Email: denys.dutykh@ku.ac.ae
// Homepage: https://www.denys-dutykh.com/
// Last modified: April 2026

#import "../styles/template.typ": dropcap, etude-conclusion, idx, chapter-abstract, exercise, hint-for


= Spectral Methods on Unbounded Intervals <ch-unbounded>

#chapter-abstract(keywords: [Unbounded domains · Subgeometric convergence · Sinc functions · Hermite functions · Rational Chebyshev functions · Domain truncation])[
On an unbounded interval one never chooses the resolution $N$ alone; one must also choose a length scale --- a truncation length, a map parameter, or a basis scale --- and this second parameter is not an implementation detail but an explicit asymptotic modelling choice. This chapter explains why that choice is unavoidable, why it usually yields subgeometric rather than geometric convergence, and how different basis families respond to different behaviours at infinity. It distinguishes three broad strategies --- domain truncation, intrinsic unbounded-domain bases, and mapped finite-interval bases --- and develops the principal toolkits: sinc functions for exponential decay, Hermite functions for Gaussian decay, Laguerre functions on the half-line, rational Chebyshev functions on the whole line and on the semi-infinite interval, and the Weideman--Cloot sinh mapping. Throughout, the organising question is what the solution does at infinity and which basis knows that behaviour best, with coefficient-decay plots used to tune the scale parameter diagnostically rather than by trial and error. The chapter also treats infinity as either a behavioural or a numerical boundary condition, and offers restricted bases and asymptotic augmentation as remedies for algebraic or oscillatory tails.
]

#dropcap[On an unbounded interval#idx("unbounded interval"), one never chooses only the resolution $N$; one also chooses a length scale. This is the deepest teaching point of the chapter. The second parameter --- a truncation length, a map parameter#idx("map parameter"), or a basis scale --- is not an implementation detail that can be chosen by habit; it is an explicit asymptotic modelling choice. We will see why it is unavoidable, why it often imposes _subgeometric_ rather than geometric convergence, and how a few concrete families of basis functions (sinc, Hermite, Laguerre, rational Chebyshev $T B_n$ and $T L_n$) respond differently to different asymptotic regimes. The chapter's pedagogical spine is not a catalogue of bases but a single organising question: _what does the solution do at infinity, and which basis knows that behaviour best?_ The answer determines the method, and the method determines the cost.]

By the end of this chapter, you should be able to:

1. Explain why unbounded intervals usually lead to subgeometric convergence#idx("subgeometric convergence"), and why they force a scale parameter in addition to $N$.
2. Distinguish the three broad strategies for infinity --- domain truncation#idx("domain truncation"), intrinsic unbounded-domain bases, mapped finite-interval bases --- and move between them confidently.
3. Choose among sinc, Hermite, Laguerre, rational Chebyshev $T B_n$, $T L_n$, and the Weideman--Cloot sinh mapping based on the asymptotic behaviour of the target solution.
4. Tune a scale parameter $L$ or $h$ empirically and diagnostically, using coefficient-decay plots rather than trial-and-error alone.
5. Treat infinity either as a behavioural boundary condition#idx("behavioural boundary condition") or as a numerical one, and recognise when each viewpoint helps.
6. Recognise when standard bases fail because the tail is algebraic or oscillatory, and deploy restricted bases (odd-SB series) or asymptotic augmentation#idx("asymptotic augmentation") as remedies.

== Why One Parameter is Never Enough <sec-un-opening>

On a finite interval, the single parameter $N$ controls resolution; the interval itself is fixed. On an unbounded interval, neither is true. Even the simplest-possible setup --- domain truncation of a well-behaved, exponentially-decaying function $u(y) = sech(y)$ onto a large finite interval $[-L, L]$ --- cannot achieve machine precision by raising $N$ alone. The reason is crisp: the truncated approximation cannot see the tail beyond $|y| = L$, and that part of the function contributes a _domain-truncation error_ $E_("DT")(L) tilde.op u(L)$, which is independent of $N$. Only when $L$ is also increased does the total error descend.

This is the opening shock of the chapter and the source of Boyd#idx("Boyd")'s Rule-of-Thumb 14 @Boyd2000:

#block(stroke: 0.8pt + rgb(20, 45, 110), radius: 3pt, inset: 10pt,
  fill: rgb(248, 250, 254))[
  *Rule-of-Thumb 14 (penalties of unbounded intervals).*  On an infinite or semi-infinite interval:

  1. The rate of convergence is usually _subgeometric_ rather than geometric.
  2. In addition to the truncation $N$, one must choose a second numerical parameter such as a truncation length $L$, a map parameter, or a basis scale.
]

=== Computational Étude 20.1: When Spectral Convergence Stalls <etude-un-stalls>

We approximate $u(y) = sech(y)$ by Chebyshev collocation on $[-L, L]$ (Dirichlet truncation) and run three independent sweeps:

+ Fix $L = 6$, increase $N$. The error plateaus at roughly $e^(-L) approx 2.5 times 10^(-3)$: no amount of resolution can reduce it.
+ Fix $N = 32$, vary $L$. A V-shaped error curve reveals a sweet spot; both too-small $L$ (large truncation error) and too-large $L$ (too many Chebyshev points wasted outside the support of the function) are suboptimal.
+ Grow $L$ and $N$ together. The error descends geometrically, but subgeometrically: each doubling of $N$ buys _less_ than a doubling of the logarithmic error.

#figure(
  image("../figures/ch20/python/truncation_stalls.pdf", width: 100%),
  caption: [Étude 20.1: domain truncation of $sech(y)$ as three slices through one error surface (four panels in a 2$times$2 grid). _(a)_: $L$ fixed at $6$, $N$ varying --- the error plateaus at the domain-truncation level $e^(-L) approx 2.5 times 10^(-3)$ (navy dashed line). _(b)_: $N$ fixed at $32$, $L$ varying --- V-shaped error, with a broad but finite sweet spot near $L approx 6$. _(c)_: $N$ and $L$ growing together --- subgeometric descent, reaching $3 times 10^(-5)$ at $(N, L) = (96, 14)$. _(d)_: the *full two-dimensional error surface* $log_(10) E(N, L)$ as a heatmap (yellow = bad, dark = good), with the three sweeps overlaid as polylines. The unifying picture: the three plots in (a)--(c) are not three different experiments but three orthogonal cuts through the same surface --- a horizontal slice, a vertical slice, and a diagonal trace along the optimum.],
) <fig-un-stalls>

#etude-conclusion[
  The three sweeps are not three different experiments; they are three slices of the same two-parameter error surface $E(N, L)$. Sweep 1 fixes $L$ and proves that the floor $e^(-L)$ exists no matter how many points we add. Sweep 2 fixes $N$ and proves that there is an optimal $L$ --- making the window too large is just as harmful as making it too small. Sweep 3 grows $N$ and $L$ jointly along the optimum and confirms the convergence is _subgeometric_: a straight line on a $log E$ versus $sqrt(N)$ plot, not the straight line on a $log E$ versus $N$ plot we would get on a finite interval. This is the empirical face of Boyd's Rule-of-Thumb 14: an unbounded interval costs us a second knob to tune, and tuning it badly costs decimals.
]

Source files:
- `codes/python/ch20/truncation_stalls.py`
- `codes/matlab/ch20/truncation_stalls.m`
- `codes/julia/ch20/truncation_stalls.jl`

== Three Routes to Infinity <sec-un-three-routes>

Having seen the trap, we organise the landscape. @Boyd2000 identifies three broad families of approach to unbounded-domain spectral problems:

1. _Domain truncation._ Replace $[-infinity, +infinity]$ by a large but finite $[-L, L]$ and use standard Chebyshev or Fourier machinery; choose $L$ so that the solution at $|y| = L$ is negligible. Simple; requires both $N$ and $L$ to grow.
2. _Intrinsic unbounded-domain bases._ Use basis functions whose natural domain is already $[-infinity, +infinity]$ (sinc, Hermite) or $[0, +infinity]$ (Laguerre). No truncation is needed; the scale is carried by the grid spacing $h$ or a scaling factor $alpha$.
3. _Mapped finite-interval bases._ Change coordinates via $y = f(x)$ with $x in [-1, 1]$, then expand in Chebyshev polynomials of $x$. The resulting basis functions are _rational_ in $y$: the $T B_n$ family on $(-infinity, +infinity)$ and the $T L_n$ family on $[0, +infinity]$.

These categories are not mutually exclusive; mapping combined with domain truncation gives the Weideman--Cloot hybrid of @sec-un-weideman-cloot. The operating principle, however, is the same across all three:

#block(stroke: 0.8pt + rgb(20, 45, 110), radius: 3pt, inset: 10pt,
  fill: rgb(248, 250, 254))[
  *Principle (choose the basis that mirrors the tail).*  The asymptotic personality of the solution determines the method. Exponential decay favours Hermite, Laguerre, or $sech$-like truncation. Algebraic decay or asymptotic constants demand rational Chebyshev. Slow algebraic decay with non-decaying oscillation demands asymptotic augmentation.
]

== Domain Truncation is the Baseline <sec-un-truncation>

Domain truncation, for all its flaws, is the honest baseline. It costs no new basis functions, no new software, and no new conceptual commitment. It is also Boyd's "no-brain" strategy @Boyd2000. Two facts make it surprisingly effective in practice:

+ For an exponentially decaying solution, $E_("DT")(L) = |u(L)|$ falls as $e^(-L)$, which for $L = 20$ is already $10^(-9)$.
+ On a truncated interval, _Fourier_ is usually more efficient than Chebyshev. Chebyshev clusters points near $y = plus.minus L$ where the decaying solution is tiny; Fourier has uniform density throughout the interior where the solution lives. The Gibbs oscillation at $y = plus.minus L$ is bounded (Boyd 1988e) by $0.09 |u(L) - u(-L)| \/ N$, which is smaller than the truncation error $|u(L)|$ itself, hence harmless.

=== Computational Étude 20.2: Fourier Beats Chebyshev at the Far Boundary <etude-un-fourier>

We compare Fourier and Chebyshev domain truncation of $sech(y)$ on $[-L, L]$ at matched degrees of freedom ($M = 2N$ Fourier modes versus $N$ Chebyshev intervals). The interior error (on $[-L + 1, L - 1]$) settles the comparison: Fourier wins decisively at low to moderate $N$, while Chebyshev overtakes asymptotically only because of its superior pointwise approximation away from the boundary.

#figure(
  image("../figures/ch20/python/fourier_vs_chebyshev_truncation.pdf", width: 85%),
  caption: [Étude 20.2: Fourier vs Chebyshev domain truncation of $sech(y)$ on $[-L, L]$ with $L = 10$. Left: grid densities at $N = 32$ (Chebyshev, $2 N = 64$ Fourier) --- Fourier has uniformly denser interior points; Chebyshev wastes resolution near $y = plus.minus L$ where the solution is already below $10^(-4)$. Right: max-norm error on the interior $[-L + 1, L - 1]$ versus $N$ --- Fourier wins by $2$-$3$ orders of magnitude until Chebyshev's geometric tail overtakes at $N approx 100$.],
) <fig-un-fourier>

#etude-conclusion[
  At low to moderate resolution, Fourier on $[-L, L]$ wins by two to three orders of magnitude, simply because Chebyshev wastes its boundary clustering on a region where $sech(y)$ is already below the truncation floor. Chebyshev only overtakes Fourier asymptotically (here near $N approx 100$), once its geometric decay finally outpaces Fourier's polynomial-in-$N$ convergence on a non-periodic target. The practical lesson is unambiguous: when an *exponentially decaying* function is being truncated to a large $[-L, L]$, switch to Fourier and pretend the function is periodic --- the periodic image that wraps the wall is exponentially small, so it costs nothing.
]

Source files:
- `codes/python/ch20/fourier_vs_chebyshev_truncation.py`
- `codes/matlab/ch20/fourier_vs_chebyshev_truncation.m`
- `codes/julia/ch20/fourier_vs_chebyshev_truncation.jl`

=== Sponge Layers: A Wall Is Not Enough

When the unbounded physics propagates _waves_ rather than decays, truncation requires more care. A naive Dirichlet boundary at $y = plus.minus L$ acts as a perfect reflector; the waves that should exit the computation return as spurious interference. Two remedies are classical:

- _Robin / outgoing-wave boundary conditions._ For waves with asymptotic structure $e^(plus.minus i k y)$, the Sommerfeld-like condition $partial_y u - i k u = 0$ at $y = L$ annihilates the incoming mode and passes the outgoing mode through. This is sharp when $k$ is unique; it fails for broadband wave packets.
- _Sponge layers._ Add an artificial viscosity $nu(y)$ that is negligible in the interior and rises smoothly near $y = plus.minus L$. The outgoing waves are damped before they reflect. Boyd @Boyd2000 warns that the sponge must vary _more slowly_ than the wave; a discontinuous viscosity reflects almost as strongly as a Dirichlet wall. A good rule of thumb is to use a $C^infinity$ ramp and size the sponge at roughly three wavelengths.

== Sinc Functions: The Simplest Infinite-Interval Method <sec-un-sinc>

The Whittaker cardinal or _sinc_ expansion is the simplest pseudospectral method on $(-infinity, +infinity)$:
$ f(y) approx sum_(j = -N \/ 2)^(N \/ 2) f(j h) op("sinc")((y - j h) \/ h), quad op("sinc")(z) = sin(pi z) \/ (pi z). $ <eq-un-sinc>
Two parameters, $N$ and $h$, must be chosen together. The total error splits into:

- a _bandwidth_ error $E_W (h)$ that shrinks exponentially with $1 \/ h$ (because $f_W (y)$, the band-limited projection of $f$, is exponentially close to $f$ for smooth decaying targets);
- a _grid-span_ error $E_("DT")(N h) = cal(O)(f(N h \/ 2))$ that shrinks with the grid span $N h$.

Optimising the balance for $f(y) = sech(y)$ yields $h^star approx sqrt(pi^2 \/ (2 N))$ and the convergence rate $exp(-pi sqrt(N \/ 2))$ --- subgeometric, but very fast in practice.

=== Computational Étude 20.3: Two Masters, Span and Spacing <etude-un-sinc>

We demonstrate Boyd's subgeometric picture concretely: sweep $h$ at fixed $N$ to find a V-shaped optimum, then follow $h = sqrt(pi^2 \/ (2 N))$ to produce a linear-in-$sqrt(N)$ descent on the log-error plot.

#figure(
  image("../figures/ch20/python/sinc_two_masters.pdf", width: 100%),
  caption: [Étude 20.3: sinc expansion of $sech(y)$. _(a)_: error versus $h$ at fixed $N = 48$; the V-shape exhibits the trade-off between bandwidth error (left branch) and grid-span error (right branch), with the empirical optimum (navy dotted) and the theoretical $sqrt(pi^2 \/ (2 N))$ (teal dashed) agreeing tightly. _(b)_: the literal _two masters_ decomposition --- bandwidth error $E_W (h) tilde.op (4 \/ pi) e^(-pi^2 \/ (2 h))$ (navy dashed) and grid-span error $E_("DT")(N h) tilde.op 4 e^(-N h \/ 2)$ (coral dashed) shown separately, with their sum (teal solid) recovering the empirical V (grey); the two analytic curves cross at the empirical $h^star$. _(c)_: max-norm error vs $sqrt(N)$ at the optimal $h$; the curve is straight on this scale, confirming the subgeometric rate $exp(-pi sqrt(N \/ 2))$ (navy dotted guide). _(d)_: grid cartoon at $N = 8$ and $N = 32$ --- as $N$ quadruples, the span doubles _and_ the spacing halves.],
) <fig-un-sinc>

#etude-conclusion[
  The V-shaped error against $h$ is the diagrammatic content of "two masters": shrink $h$ and the bandwidth error vanishes but the grid-span $N h$ collapses; grow $h$ and the span is large enough to see the tail but the band cannot resolve the bulk. The minimum $h^star = sqrt(pi^2 \/ (2 N))$ tracks both errors at the same level, and the resulting joint convergence $exp(-pi sqrt(N \/ 2))$ is *subgeometric* but practically very fast --- the middle panel is a straight line on a $log E$-versus-$sqrt(N)$ plot, exactly the diagnostic shape of subgeometric decay. Sinc is the cheapest infinite-interval method, but its scale parameter $h$ must be tuned with $N$; using a fixed $h$ across resolutions silently throws decimals away.
]

Source files:
- `codes/python/ch20/sinc_two_masters.py`
- `codes/matlab/ch20/sinc_two_masters.m`
- `codes/julia/ch20/sinc_two_masters.jl`

== Hermite Functions: When the Basis Understands the Physics <sec-un-hermite>

The normalised Hermite functions
$ psi_n (y) = (2^n n! sqrt(pi))^(-1 \/ 2) H_n (y) e^(-y^2 \/ 2) $ <eq-un-hermite-def>
are the eigenfunctions of the quantum harmonic oscillator and the asymptotic eigenfunctions of several other important physical problems (Mathieu, prolate spheroidal, associated Legendre, Laplace tidal equation). Their convergence theory, developed by Hille (1939-1961) and sharpened in Boyd (1984a), is deep but practical. The key results @Boyd2000:

- For $f$ with _Gaussian_ decay ($f tilde.op e^(-p y^2)$), Hermite gives _geometric_ convergence --- the fastest possible.
- For _super-Gaussian_ decay ($f tilde.op e^(-p |y|^k)$, $k > 2$), Hermite is subgeometric with index $k \/ (2 (k - 1))$.
- For _sub-Gaussian_ decay (the common case, $k < 2$), Hermite is subgeometric with index $k \/ 2$.
- For _algebraic_ decay, Hermite converges only algebraically.

The lesson is: Hermite is a precision instrument, not a general-purpose tool. When the solution has Gaussian character and a well-matched scale, Hermite is spectacular; when it does not, Hermite wastes resolution on a decay rate the target does not share.

=== Scaling Matters: $alpha = sqrt(2 A)$

For a target $e^(-A y^2)$ with $A eq.not 1 \/ 2$, the unscaled basis $psi_n (y)$ will be poor because its built-in length scale ($sqrt(2)$) does not match the target's length scale ($1 \/ sqrt(A)$). The cure is to use the scaled basis
$ phi_n (y) = sqrt(alpha) psi_n (alpha y), quad alpha = sqrt(2 A), $ <eq-un-hermite-scale>
which restores geometric convergence. For $A = 1 \/ 2$, $alpha = 1$ and the match is automatic.

=== Computational Étude 20.4: The Cost of Width Mismatch <etude-un-hermite>

We expand $f(y) = e^(-A y^2)$ for four values of $A$ with and without scaling, and separately expand the quantum-oscillator ground state $psi_0 (y) = pi^(-1 \/ 4) e^(-y^2 \/ 2)$, which the Hermite basis represents exactly with a single coefficient.

#figure(
  image("../figures/ch20/python/hermite_width_mismatch.pdf", width: 100%),
  caption: [Étude 20.4: Hermite width mismatch (six panels in a 3$times$2 grid). _(a)_: cartoon view --- targets $e^(-A y^2)$ for $A in {0.1, 0.5, 2, 8}$ overlaid with the unscaled basis envelope $psi_0(y) = pi^(-1\/4) e^(-y^2 \/ 2)$ (dashed grey); $A = 0.5$ tracks the envelope exactly. _(b)_: coefficient decay at $A = 8$ --- unscaled (coral, slow algebraic decay) versus matched (navy, only $a_0 != 0$). _(c)_: unscaled ($alpha = 1$) convergence --- $A = 0.5$ is machine-precision exact by coincidence, all other $A$ converge only slowly. _(d)_: matched scaling $alpha = sqrt(2 A)$ --- all four curves collapse to machine precision at every $N$. _(e)_: optimal-$alpha$ scan at fixed $N = 16$ for three $A$ values; each curve has a sharp dip exactly at $alpha^star = sqrt(2 A)$ (dotted vertical lines). _(f)_: quantum oscillator ground state $psi_0$, machine-precision at $N = 0$ --- the ideal case where the basis matches the physics by construction.],
) <fig-un-hermite>

#etude-conclusion[
  The middle and right panels make the same point in two different ways. With $alpha = sqrt(2 A)$, the scaled basis is *literally* the eigenbasis of the target's Gaussian, so a single coefficient suffices --- the convergence curve is a flat line at machine precision. The unscaled basis, by contrast, only achieves this for the special value $A = 1 \/ 2$ where the widths happen to coincide; for other $A$, *every additional coefficient is paying interest on a width mismatch* and the convergence is no better than algebraic. The takeaway is operational: when using Hermite, the *scale* $alpha$ must be matched to the solution width before $N$ can do any work. This is the same lesson as Étude 20.1 (the second knob) seen through a different basis.
]

Source files:
- `codes/python/ch20/hermite_width_mismatch.py`
- `codes/matlab/ch20/hermite_width_mismatch.m`
- `codes/julia/ch20/hermite_width_mismatch.jl`

== Laguerre Functions on the Semi-Infinite Interval <sec-un-laguerre>

The Laguerre function#idx("Laguerre function")s $phi_n (y) = e^(-y \/ 2) L_n (y)$, $y in [0, +infinity)$, are close cousins of Hermite: product of an exponentially decaying envelope and an orthogonal polynomial. They are natural for problems with pure exponential decay at infinity --- above all, the hydrogen-atom radial equations and quantum-chemistry orbitals. Outside these physical contexts, Laguerre is more limited than the rational-Chebyshev $T L_n$ basis; Iranzo & Falqués (1992) document this in detail.

The key weakness, again, is asymptotic. Laguerre's envelope $e^(-y \/ 2)$ decays exponentially, so Laguerre fails on any function that decays _algebraically_ --- even the innocuous $1 \/ (1 + y)$. The rational $T L_n$ basis, by contrast, handles exponential, algebraic, and asymptote-to-a-constant cases with equal ease. We defer the formal definition of $T L_n$ to @sec-un-tln but already introduce them for the comparison.

=== Computational Étude 20.5: A Semi-Infinite Basis Meets a Semi-Infinite Problem <etude-un-lagvstln>

Two target functions on $[0, +infinity)$:
$ f_1 (y) = e^(-y) quad "(pure exponential)", quad f_2 (y) = 1 \/ (1 + y) quad "(algebraic $\sim 1 \/ y$)". $ <eq-un-lagtln>
Laguerre and $T L_n$ are compared on both. Laguerre reaches machine precision on $f_1$ at $N approx 20$, but plateaus at $10^(-3)$ on $f_2$ --- an accuracy floor no amount of resolution can break. $T L_n$ handles both, winning spectacularly on $f_2$.

#figure(
  image("../figures/ch20/python/laguerre_vs_tln.pdf", width: 85%),
  caption: [Étude 20.5: Laguerre vs rational $T L_n$ on $[0, +infinity)$. Left: pure exponential $f = e^(-y)$; both bases converge geometrically, Laguerre slightly faster thanks to basis-physics match. Right: algebraic $f = 1 \/ (1 + y)$; Laguerre plateaus near $10^(-3)$ (its exponential envelope cannot track $1 \/ y$), while $T L_n$ descends to machine precision by $N = 32$.],
) <fig-un-lagtln>

#etude-conclusion[
  Both bases live on $[0, +infinity)$; both are spectrally accurate for problems they are *designed* for. The two panels show that "designed for" is a sharp clause. Laguerre's $e^(-y \/ 2) L_n (y)$ envelope can model only exponential decay --- it has no room left for $1 \/ y$, so the right panel saturates at $10^(-3)$ no matter how many modes we add. The rational $T L_n$ basis carries no fixed envelope; its members asymptote to $plus.minus 1$ at infinity, so any decay rate (exponential, algebraic, or none at all) is reachable by recombining them. The general rule that emerges: *prefer rational Chebyshev as the default semi-infinite basis*; reserve Laguerre for problems where the physics literally produces $e^(-y \/ 2) L_n^k$ structure (notably hydrogenic radial equations).
]

Source files:
- `codes/python/ch20/laguerre_vs_tln.py`
- `codes/matlab/ch20/laguerre_vs_tln.m`
- `codes/julia/ch20/laguerre_vs_tln.jl`

== Rational Chebyshev Functions $T B_n (y)$ on $(-infinity, +infinity)$ <sec-un-tbn>

The central section of the chapter. The rational Chebyshev basis is defined by the map
$ y = ell x \/ sqrt(1 - x^2), quad x in [-1, 1], quad T B_n (y) = T_n (x) = cos(n t) quad "with" quad y = ell cot(t). $ <eq-un-tbn-def>
Geometrically: $T B_n$ is the image of the $n$-th Chebyshev polynomial under an _algebraic_ map from $[-1, 1]$ to $(-infinity, +infinity)$. Explicit first few:
$ T B_0 = 1, quad T B_1 = y \/ sqrt(1 + y^2), \
  T B_2 = (y^2 - 1) \/ (y^2 + 1), quad T B_3 = y(y^2 - 3) \/ (y^2 + 1)^(3 \/ 2). $ <eq-un-tbn-explicit>
The even-degree members are rational functions symmetric about $y = 0$; the odd-degree members are antisymmetric and carry a residual factor of $sqrt(ell^2 + y^2)$.

The decisive property: all $T B_n$ _asymptote to a constant_ at $y = plus.minus infinity$. This makes them exactly the right basis for functions that decay algebraically or approach a constant. Boyd @Boyd2000 summarises: _rational Chebyshev functions are the basis of choice for $f(y)$ that decay algebraically rather than exponentially, or which asymptote to a constant._

=== Computational Étude 20.6: The Function that Humiliates Hermite <etude-un-tbn>

Boyd's emblematic example:
$ f(y) = 1 \/ (1 + y^2). $ <eq-un-humiliate>
Bounded, smooth, symmetric --- but with simple poles at $y = plus.minus i$, so it decays only as $1 \/ y^2$. Hermite converges algebraically (exactly what the theorem of @sec-un-hermite predicts for algebraic decay); sinc does likewise. The $T B_n$ basis, by contrast, is exact at $N = 2$: because $f = 1 \/ 2 (T B_0 - T B_2)$ precisely for $ell = 1$.

#figure(
  image("../figures/ch20/python/tbn_humiliates_hermite.pdf", width: 85%),
  caption: [Étude 20.6: the function that humiliates Hermite. Left: $f(y) = 1 \/ (1 + y^2)$ --- an innocuous target. Right: three bases compared on a log-log error plot. Hermite converges algebraically ($"err" tilde.op 1 \/ sqrt(N)$, cf. Theorem 34 of @Boyd2000); sinc likewise. The $T B_n$ basis reaches machine precision at $N = 8$ because $f$ is literally in the span of $(T B_0, T B_2)$ for $ell = 1$. The whole chapter hinges on this comparison.],
) <fig-un-humiliate>

#etude-conclusion[
  The function $1 \/ (1 + y^2)$ is the reductio ad absurdum of "exponential bases for the unbounded line": Hermite and sinc both decay polynomially because the target decays only as $1 \/ y^2$, while $T B_n$ recovers it *exactly* at $N = 2$ thanks to the algebraic identity $1 \/ (1 + y^2) = (T B_0 - T B_2) \/ 2$ at $ell = 1$. This is more than a curiosity: it operationalises the Principle of @sec-un-three-routes by showing what happens when one ignores it. The lesson is sharp --- *the right basis for an algebraically decaying target is rational, not exponential*; using Hermite or sinc here is not "slower spectral convergence" but a category mismatch between basis and tail.
]

Source files:
- `codes/python/ch20/tbn_humiliates_hermite.py`
- `codes/matlab/ch20/tbn_humiliates_hermite.m`
- `codes/julia/ch20/tbn_humiliates_hermite.jl`

== Behavioural versus Numerical Boundary Conditions at Infinity <sec-un-behavioural>

Boundary conditions at infinity are fundamentally different from those at finite endpoints. On a finite interval, one imposes conditions like $u(-1) = u(+1) = 0$ _essentially_: the basis is modified to satisfy them exactly. At infinity, the requirement is almost always _behavioural_: we ask only that the solution be bounded (or have a specified asymptotic type), not that it equal a specified value. The rational-Chebyshev basis $T B_n$ and $T L_n$ both respect boundedness automatically, so the "condition at infinity" need not be imposed as an explicit constraint.

#block(stroke: 0.8pt + rgb(20, 45, 110), radius: 3pt, inset: 10pt,
  fill: rgb(248, 250, 254))[
  *Rule-of-Thumb 15 (behavioural conditions at infinity, @Boyd2000).*  If the desired solution of a differential equation is the only solution that is bounded at infinity, an unmodified series of rational functions $T B_n$ or $T L_n$ usually gives exponential accuracy.  If this fails, impose the constraint $u(infinity) = 0$ explicitly.
]

The rule is pragmatic, not theoretical: _try behavioural first_, and fall back on explicit constraints only if the computation misbehaves.

=== Computational Étude 20.7: Boundedness Is Sometimes Enough <etude-un-behavioural>

Boyd's semi-infinite Laguerre eigenvalue problem:
$ y u''(y) + (y + 1) u'(y) + lambda u(y) = 0, quad y in [0, +infinity), $ <eq-un-laguerre-evp>
with exact spectrum $lambda_n = n$ for $n = 0, 1, 2, dots$ and eigenfunctions $u_n (y) = e^(-y) L_n^1(y)$ (associated Laguerre). _Strategy A_: discretise directly on the $T L_n$ grid with no boundary constraint. _Strategy B_: change the unknown to $w = e^(y \/ 2) u(y)$, transforming the problem so that the unphysical blow-up solutions grow _exponentially_ rather than algebraically. Both in principle compute the same spectrum; in practice, Strategy B dominates.

#figure(
  image("../figures/ch20/python/behavioral_bc_laguerre.pdf", width: 85%),
  caption: [Étude 20.7: two strategies for Boyd's Laguerre eigenproblem. Left: count of "good" eigenvalues (within $5 %$ of an integer) versus $N$. Strategy A (naive, coral) reaches $33$ good eigenvalues at $N = 120$; Strategy B (behavioural recast, teal) reaches $49$. Right: the computed spectrum at $N = 40$. Strategy A (coral circles) shows complex-conjugate pairs with imaginary parts of order $0.01$ --- contamination of the real spectrum by slowly-decaying unphysical solutions. Strategy B (teal squares) returns clean real eigenvalues at $lambda_n = n$.],
) <fig-un-behavioural>

#etude-conclusion[
  Strategy A is what the unmodified rational basis offers: nothing is imposed at infinity, but the mathematics already requires the wanted eigenfunction to be the only one bounded there. It works well, recovering 33 clean eigenvalues at $N = 120$, but the unphysical solutions that grow only *algebraically* contaminate the spectrum with small imaginary parts at higher modes. Strategy B replaces the unknown by $w = e^(y \/ 2) u$, so the parasitic solutions now grow *exponentially* and are pushed far from the real axis: 49 clean eigenvalues at the same $N$, with no imaginary contamination. The lesson generalises: *behavioural conditions at infinity are usually enough, but if a slowly-growing parasite contaminates the matrix, change the unknown so that the parasite grows fast enough to be filtered out by the basis*.
]

Source files:
- `codes/python/ch20/behavioral_bc_laguerre.py`
- `codes/matlab/ch20/behavioral_bc_laguerre.m`
- `codes/julia/ch20/behavioral_bc_laguerre.jl`

== Slow Decay and Restricted Bases: the $S B_n$ Family <sec-un-sbn>

When the solution decays algebraically with _odd_ inverse powers of $y$ (like $tilde.op 1 \/ y$), the antisymmetric subfamily of $T B_n$ has the _wrong_ asymptotic type: each $T B_(2n + 1)$ carries a factor of $sqrt(ell^2 + y^2)$ that does not decay. The remedy is a companion family: the _sine_-based rational functions $S B_n (y)$ defined by
$ S B_n (y) = sin{(n + 1) op("arccot")(y \/ ell)}. $ <eq-un-sbn-def>
The $S B_n$ are rational in $y$ with odd-inverse-power asymptotics, and only the odd-index subfamily $S B_(2 k + 1)$ is needed for antisymmetric targets.

Boyd's Table 17.6 classifies four symmetry-and-asymptotics cases (symmetric/antisymmetric, even/odd inverse powers) into four restricted bases: $T B_(2 n)$, $T B_(2 n + 1)$, $S B_(2 n)$, and $S B_(2 n + 1)$. The practical upshot is that looking at the asymptotic behaviour of $u(y)$ before choosing the basis can save an order of magnitude in $N$.

=== Computational Étude 20.8: The Yoshida Jet <etude-un-yoshida>

The steady-state Yoshida jet#idx("Yoshida jet") in equatorial oceanography satisfies
$ v_(y y) - y^2 v = y, quad y in (-infinity, +infinity), $ <eq-un-yoshida>
with $v$ antisymmetric and decaying as $-1 \/ y$ at infinity. Odd-$S B_n$ collocation converges rapidly; Boyd's Table 17.7 quotes coefficients that match our computation to six decimal places at $N = 21$.

#figure(
  image("../figures/ch20/python/yoshida_jet.pdf", width: 85%),
  caption: [Étude 20.8: Yoshida jet velocity profile $v(y)$ solving $v'' - y^2 v = y$ with $v tilde.op -1 \/ y$. Left: the reference (21-mode $S B$ series) compared with low-$N$ truncations --- $N = 3$ is already accurate to three decimals in the interior. Right: convergence with $N$; the odd-$S B$ basis matches the algebraic $1 \/ y$ tail by construction, so the series descends rapidly to $10^(-5)$ by $N = 10$.],
) <fig-un-yoshida>

#etude-conclusion[
  The Yoshida jet is antisymmetric and decays as $-1 \/ y$ --- an *odd* inverse power that the standard rational $T B_n$ basis cannot represent efficiently because its odd-degree members carry a residual $sqrt(ell^2 + y^2)$. The odd-$S B_n$ basis is built precisely to match this asymptotic class; matching the basis to the symmetry class collapses the numerical problem to a few unknowns ($N = 3$ already accurate to three decimals, $N approx 10$ to five). The pattern recurs throughout the chapter: *taking a few minutes to classify the asymptotic personality of the solution pays back an order of magnitude in $N$*. A textbook $T B_n$ run on this same problem would converge much more slowly, despite using the same rational-Chebyshev infrastructure.
]

Source files:
- `codes/python/ch20/yoshida_jet.py`
- `codes/matlab/ch20/yoshida_jet.m`
- `codes/julia/ch20/yoshida_jet.jl`

== Rational Chebyshev Functions $T L_n (y)$ on $[0, +infinity)$ <sec-un-tln>

The semi-infinite analogue of $T B_n$ is defined by the map
$ y = ell (1 + x) \/ (1 - x), quad x in [-1, 1], quad T L_n (y) = T_n (x) = cos(n t) quad "with" quad y = ell cot^2(t \/ 2). $ <eq-un-tln-def>
First few:
$ T L_0 = 1, quad T L_1 = (y - ell) \/ (y + ell), quad T L_2 = (y^2 - 6 y ell + ell^2) \/ (y + ell)^2. $ <eq-un-tln-explicit>
The basis handles exponential, algebraic, and asymptote-to-a-constant decay at $y = infinity$ with equal ease, just like $T B_n$. Étude 20.5 already compared $T L_n$ to Laguerre and revealed $T L_n$'s strength on algebraic decay.

#figure(
  image("../figures/ch20/python/tln_map_illustration.pdf", width: 100%),
  caption: [The rational Chebyshev $T L_n$ map and basis on the half-line (four panels in a 2$times$2 grid). _(a)_: the algebraic semi-infinite map $y = ell, (1 + x) \/ (1 - x)$ for four values of the map parameter $ell$; the singularity at $x = 1$ pulls $x in (-1, 1)$ onto $y in [0, +infinity)$. Larger $ell$ delays the rise, scattering grid points further out. _(b)_: the corresponding collocation grids @eq-un-tln-def at fixed $N = 24$, plotted as horizontal tick rows over a generic $sech$-decaying target (grey). $ell = 1$ packs the grid near $y = 0$; $ell = 8$ pushes most nodes far out. _(c)_: the first five basis functions $T L_0, dots, T L_4$ at $ell = 2$, evaluated over $y in [0, 20]$. Each $T L_n$ asymptotes to $plus.minus 1$ at $y = +infinity$, making the basis natural for problems whose tail is algebraic or asymptote-to-a-constant. _(d)_: coefficient diagnostic --- $|a_n|$ vs $n$ for the $T L_n$ expansion of $sech(y)$ at $ell in {0.5, 2, 8}$; too-small $ell$ flattens early (domain truncation), too-large $ell$ slows the small-$n$ slope, and $ell approx 2$ descends fastest, foreshadowing Boyd's Rule-of-Thumb 16 (Étude 20.9).],
) <fig-tln-map>

== Choosing the Map Parameter $ell$ <sec-un-ell>

Every rational-Chebyshev calculation on an unbounded domain requires a choice of the map parameter $ell$. The usual advice is "take $ell$ of the same order of magnitude as the solution width", but this leaves students to guess. Boyd's Rule-of-Thumb 16 @Boyd2000 replaces guessing with diagnosis:

#block(stroke: 0.8pt + rgb(20, 45, 110), radius: 3pt, inset: 10pt,
  fill: rgb(248, 250, 254))[
  *Rule-of-Thumb 16 (optimising the map parameter).*  Plot the coefficients $a_j$ versus degree $j$ on a log-linear plot. If the graph _abruptly flattens_ at some $j = j^star < N$, then $ell$ is too small: the remaining coefficients contribute the domain-truncation-like tail. Increase $ell$ until the flattening is postponed to $j = N$.
]

=== Computational Étude 20.9: Read the Coefficients Before the Error <etude-un-ell>

We expand $sech(y)$ in the $T B_n$ basis at $N = 64$ for six values of $ell$ ranging from $0.5$ to $16$. The coefficient decay plots tell the whole story: too-small $ell$ produces an early flattening of $|a_n|$; too-large $ell$ produces a gentler small-$n$ slope (a pole of the mapped function is drifting toward $x = plus.minus 1$ in Chebyshev space). A broad valley of good $ell$ is clearly visible in the second panel, and it widens as $N$ grows.

#figure(
  image("../figures/ch20/python/ell_diagnostic.pdf", width: 85%),
  caption: [Étude 20.9: reading $ell$ from coefficient decay (four panels in a 2$times$2 grid). _(a)_: rolling-max envelope of $|a_n|$ for the three regimes $ell in {0.5, 2, 16}$ at $N = 64$ (small / good / large), with parity zig-zag suppressed --- $ell = 0.5$ flattens early, $ell = 16$ has a visibly gentler small-$n$ slope, $ell = 2$ descends cleanly to machine precision. _(b)_: the full six-$ell$ sweep as envelopes only, letting the reader compare every regime without overplotting. _(c)_: the scalar tail-size diagnostic $sum_(n > N \/ 2) |a_n|$ versus $ell$ at $N in {24, 48, 96}$ --- the valley of good $ell$ is broad and widens with $N$. _(d)_: extracting the valley centre and width from (c) as functions of $N$ on a log-log axis; the centre drifts mildly upward with $N$ and the width grows so that larger problems forgive coarser $ell$-tuning.],
) <fig-un-ell>

#etude-conclusion[
  Boyd's Rule-of-Thumb 16 promises a *quantitative* diagnostic for $ell$, and the left panel delivers it: at $ell = 0.5$ the coefficient curve breaks horizontally near $n approx 20$, signalling that the remaining 44 modes contribute the same kind of plateau-tail seen in Étude 20.1 for the truncation-window error. At $ell = 16$ the curve never breaks, but its small-$n$ slope is gentler --- a pole of the mapped function is drifting toward $x = plus.minus 1$ in Chebyshev space, wasting resolution on a different end of the spectrum. The middle values $ell in {1, 2, 4}$ produce monotone, geometric decay all the way to machine precision. The right panel turns the *qualitative shape* of the coefficient plot into a *scalar diagnostic* (the high-mode tail sum) whose minimum reliably lands inside the good-$ell$ valley; the valley *widens* with $N$, so larger problems forgive coarser tuning.
]

Source files:
- `codes/python/ch20/ell_diagnostic.py`
- `codes/matlab/ch20/ell_diagnostic.m`
- `codes/julia/ch20/ell_diagnostic.jl`

== Singular Endpoints via Composed Maps <sec-un-composed>

The rational-Chebyshev framework composes beautifully with the coordinate-transformation toolkit of @ch-coord-transforms. For a semi-infinite problem with a logarithmic singularity at $y = 0$ _and_ exponential decay at $y = infinity$, a single map $y = op("arcsinh")(e^z)$ handles both endpoints at once: for large $z$, $y tilde.op e^z \/ 2$ (linear in log-space, good for exponential decay); for $z arrow -infinity$, $y tilde.op e^z$ (exponentially fine spacing near $y = 0$, resolving the logarithmic singularity).

=== Computational Étude 20.10: One Global Expansion for $r K_1 (r)$ <etude-un-rk1>

The modified Bessel function $r K_1 (r)$, $r in [0, +infinity)$, has a logarithmic singularity at $r = 0$ ($r K_1 (r) = 1 + (r^2 \/ 2) log(r \/ 2) + dots$) and exponential decay at infinity. Standard library software patches _two_ expansions together: a power series with log terms for small $r$, and an asymptotic series in $1 \/ r$ for large $r$. Boyd's composed map#idx("composed map") replaces both with a single $T B_n$ expansion in $z$.

#figure(
  image("../figures/ch20/python/rK1_composed_map.pdf", width: 85%),
  caption: [Étude 20.10: global expansion of $r K_1 (r)$. Left: the function (navy) and the $N = 16$ $T B_n$ approximation in the composed coordinate (coral dashed); graphical agreement is within the line thickness. Right: max-norm error versus $N$ for $ell = 4$; a single global expansion reaches $10^(-10)$ at $N = 64$ with no piecewise switching logic.],
) <fig-un-rk1>

#etude-conclusion[
  Library implementations of $r K_1 (r)$ use *two* expansions glued together: a power series with logarithmic terms near $r = 0$, and an asymptotic series in $1 \/ r$ for large $r$. The composed map $y = op("arcsinh")(e^z)$ replaces both with a *single* $T B_n$ expansion in $z$, reaching $10^(-10)$ at $N = 64$. Two ideas combine: near $r = 0$, $z arrow -infinity$ provides exponentially fine spacing that resolves the logarithmic singularity; near $r = +infinity$, $z$ is roughly $log(r)$, which is the natural variable for exponential decay. The takeaway generalises beyond Bessel functions: *whenever a problem has fundamentally different asymptotic behaviours at its two ends, look for a coordinate that aligns both* --- a single global expansion will then beat any patchwork of local series.
]

Source files:
- `codes/python/ch20/rK1_composed_map.py`
- `codes/matlab/ch20/rK1_composed_map.m`
- `codes/julia/ch20/rK1_composed_map.jl`

== Oscillatory Tails and Asymptotic Augmentation <sec-un-oscillatory>

The last serious obstacle. A function that decays slowly _and_ oscillates at a non-decaying frequency --- the canonical example is $J_0 (y) tilde.op sqrt(2 \/ (pi y)) cos(y - pi \/ 4)$ --- is beyond the reach of standard rational-Chebyshev bases: no truncation $sum c_n T L_n (y)$ of finite length can resolve infinitely many wavelengths. The cure is _asymptotic augmentation_: write
$ sqrt(1 + y) J_0 (y) = a(y) cos(y - pi \/ 4) + b(y) sin(y - pi \/ 4), $ <eq-un-aug>
where $a(y)$ and $b(y)$ asymptote to constants at infinity. Both $a$ and $b$ are then expanded in $T L_n$ series; the oscillation lives entirely in the cosine and sine _carriers_, and the basis need only represent the slowly-varying amplitude and phase.

=== Computational Étude 20.11: Teach the Basis the Oscillation First <etude-un-j0>

We fit $g(y) = sqrt(1 + y) J_0 (y)$ by two methods: a naive $sum c_n T L_n (y)$ and the augmented $sum a_n T L_n (y) cos(y - pi \/ 4) + sum b_n T L_n (y) sin(y - pi \/ 4)$, both on a dense sample.

#figure(
  image("../figures/ch20/python/J0_oscillatory.pdf", width: 85%),
  caption: [Étude 20.11: asymptotic augmentation for $sqrt(1 + y) J_0 (y)$. Left: the target and its amplitude-phase decomposition at $N = 15$; $|a(y)|$ asymptotes to a constant $approx 0.89$ and $|phi(y)|$ to zero, confirming the WKB structure. Right: max-norm error on $[0, 50]$ versus $N$. The naive basis stays pinned near $1$ (it cannot resolve the oscillation); the augmented basis reaches $10^(-13)$ at $N = 20$.],
) <fig-un-j0>

#etude-conclusion[
  The naive curve in the right panel never moves: no $T L_n$ expansion of finite length can resolve infinitely many wavelengths, because the basis itself does not oscillate at the Bessel frequency. *Asymptotic augmentation* repairs this by *factoring out the oscillation* into known carriers $cos(y - pi \/ 4)$ and $sin(y - pi \/ 4)$, leaving the basis to represent only the slowly-varying amplitude $a(y)$ and phase $b(y)$. Both of these *do* asymptote to constants at infinity --- exactly the regime $T L_n$ was built for --- and the augmented method reaches $10^(-13)$ at $N = 20$. The general principle: *teach the basis the oscillation first, the amplitude afterwards*; whenever the WKB structure of the solution is known a priori, it should be built into the ansatz, not learned by the basis at extra cost.
]

Source files:
- `codes/python/ch20/J0_oscillatory.py`
- `codes/matlab/ch20/J0_oscillatory.m`
- `codes/julia/ch20/J0_oscillatory.jl`

== The Weideman--Cloot Sinh Mapping <sec-un-weideman-cloot>

A useful hybrid of mapping and domain truncation. Weideman & Cloot (1990) propose
$ y = sinh(ell t), quad t in [-pi, pi], quad y in [-sinh(ell pi), sinh(ell pi)] equiv [-y_("max"), y_("max")]. $ <eq-un-wc>
A standard Fourier pseudospectral method is then applied in $t$. Since $y_("max") = sinh(ell pi)$ grows _exponentially_ with $ell$, a modest logarithmic increase in $ell$ is enough to drive the domain-truncation error geometrically to zero. The Weideman--Cloot method is often more efficient than $T B_n$ at fixed $N$, particularly for $sech$-like functions with small to moderate $N$, but is also more sensitive to the choice of $ell$.

Chapter-final comparison (Boyd's Fig. 17.12): both methods reach $10^(-11)$ for $sech(y)$ at $N = 22$, but Weideman--Cloot's error valley in $ell$ is sharper (a factor of $30 %$ mistuning can cost two decimal places). Rational Chebyshev wins for robustness; Weideman--Cloot wins for raw efficiency when $ell$ is well-tuned.

#figure(
  image("../figures/ch20/python/weideman_cloot_map.pdf", width: 100%),
  caption: [The Weideman--Cloot sinh map @eq-un-wc and its grid behaviour (four panels in a 2$times$2 grid). _(a)_: the map $y = sinh(ell, t)$ for four values of $ell$; the reach $y_("max") = sinh(ell pi)$ is annotated and grows _exponentially_ in $ell$ (about $2$ at $ell = 0.5$, near $270$ at $ell = 2$). _(b)_: the resulting grids on the $y$-line at fixed $N = 64$, plotted as horizontal tick rows over a $sech(y)$ target. The uniform-in-$t$ Fourier grid maps to a non-uniform-in-$y$ grid that compresses near $y = 0$ (where $y'(0) = ell$ is small) and stretches into the tails. _(c)_: the same $sech(y)$ target shown in $y$-space (navy) and in $t$-space (coral, dashed, with $t$ rescaled to share the $y$-axis). In $t$-space the function becomes a smooth, periodic-friendly profile that Fourier represents efficiently --- this is exactly the trade-off that makes Weideman--Cloot work. _(d)_: Fourier-in-$t$ coefficient diagnostic --- $|hat f_k|$ of $sech(y(t))$ at fixed $N = 96$ for $ell in {0.5, 1, 1.5, 2}$. The optimal $ell$ is the one whose coefficients descend furthest before flattening, parallel to Étude 20.9's $ell$ diagnostic on the rational-Chebyshev side.],
) <fig-wc-map>

== A Decision Tree for Unbounded Intervals <sec-un-decision>

We close with a compact operational guide.

#block(stroke: 0.8pt + rgb(20, 45, 110), radius: 4pt, inset: 12pt,
  fill: rgb(248, 250, 254))[
  *Decision guide for unbounded intervals.*  Ask, in order:

  1. _Is the domain infinite or semi-infinite?_  $(-infinity, +infinity)$ uses $T B_n$, Hermite, or sinc; $[0, +infinity)$ uses $T L_n$, Laguerre.
  2. _What is the asymptotic personality of the solution?_
     - Exponential decay: Hermite / Laguerre / $T L_n$ all competitive, in order of decreasing specificity.
     - Gaussian decay on $(-infinity, +infinity)$: Hermite (geometric convergence).
     - Algebraic decay or asymptote to a constant: rational Chebyshev $T B_n$ / $T L_n$ (the only spectrally-accurate option).
     - Oscillatory non-decaying tail: asymptotic augmentation on top of $T L_n$.
  3. _Do I need a behavioural or explicit BC at infinity?_  Try behavioural; fall back on a change-of-unknown (like $w = e^(y \/ 2) u$) if a slowly-growing parasitic solution contaminates the matrix.
  4. _How do I choose the scale parameter?_  Plot the coefficients. If they flatten early, increase $ell$.
  5. _Does the solution have a singular or difficult endpoint?_  Compose two maps (like $y = op("arcsinh")(e^z)$) so a single expansion resolves both ends.
  6. _Do I need FFT speed, favourable timestep, or easy implementation?_  Rational Chebyshev is FFT-compatible and has spectral radius $cal(O)(1 \/ N)$ for the first derivative, matching finite differences. Hermite and Laguerre do not admit standard FFTs.
]

The chapter's closing maxim, to borrow Boyd's own phrasing: _on unbounded intervals, basis choice is asymptotic modelling_. If the basis knows the tail, the expansion converges geometrically; if it does not, no amount of resolution will save the computation.

== A Non-Exhaustive Literature Review <sec-un-lit-review>

The classical reference for unbounded-domain spectral methods is @Boyd2000 Chapter 17, from which this chapter descends directly. Boyd's own survey papers (1987a, 1987b, 1987c) develop the rational-Chebyshev $T B_n$ and $T L_n$ families, and Grosch & Orszag (1977) document the earliest serious comparison of logarithmic and algebraic maps.

Sinc methods have a literature of their own: Stenger (1981, 1993) wrote the definitive review and monograph; Lund & Bowers (1992) give a concrete implementation guide. The rate-of-convergence theory for Hermite functions is due to Hille (1939-1961), with sharpening by Bain (1978) and especially Boyd (1984a) for entire functions via steepest descent. For Laguerre, the quantum-chemistry literature is vast; Iranzo & Falqués (1992) give the only serious detailed comparison with rational Chebyshev. Weideman (1992) proves that rational-Chebyshev differentiation matrices have the same favourable skew-symmetric / symmetric structure as sinc and Hermite.

The Weideman--Cloot sinh map (Weideman & Cloot 1990, Cloot & Weideman 1992) is a workhorse in plasma physics and optics; Boyd (1994b) gives the most balanced comparison with rational Chebyshev. For oscillatory tails, Boyd (1987b, 1990d, 1995a) developed the asymptotic-augmentation framework that we used in Étude 20.11; the resulting "radiation function" bases are central to weakly nonlocal solitary-wave theory.

Among modern developments, the ultraspherical spectral method @OlverTownsend2013 @OlverTownsend2019 applies naturally to mapped bases and resolves some of the ill-conditioning that plagued earlier rational-basis implementations. In relativistic astrophysics, rational-Chebyshev-like methods are the standard tool for quasinormal-mode computations on unbounded domains, as documented in the Batić--Dutykh series cited in @ch-linear-eigen. For mapped Legendre bases applied to vortex stability in unbounded fluid domains, @LeeMarcus2023 provides a state-of-the-art demonstration of roughly three-fold efficiency over rational Chebyshev in specific cases.

== Summary and Exercises <sec-un-summary>

The exercises below progress from pencil-and-paper analysis of convergence rates and rational-basis identities, through numerical experiments that reproduce and extend the études of this chapter, to open-ended projects that reach into the research literature. The computational problems may be carried out in any of the book's three languages; the named scripts under `codes/python/ch20/`, `codes/julia/ch20/`, and `codes/matlab/ch20/` give a starting point.

=== Conceptual Exercises

#exercise[
  Derive Boyd's estimate $0.09 |u(L) - u(-L)| \/ N$ for the Gibbs overshoot of the Fourier truncation of a smooth decaying function on $[-L, L]$.
] <ex-un-gibbs-overshoot>

#exercise[
  Show that for $f(y) = sech(y)$, the optimal sinc spacing is $h^star = sqrt(pi^2 \/ (2 N))$ and the resulting error is $cal(O)(exp(-pi sqrt(N \/ 2)))$.
] <ex-un-sinc-optimal>

#hint-for(<ex-un-sinc-optimal>)[Write the total error as the sum of the bandwidth term $tilde.op e^(-pi^2 \/ (2 h))$, set by the strip half-width $pi \/ 2$ of $sech$, and the grid-span term $tilde.op e^(-N h)$, the tail beyond the outermost of the $N$ nodes on each side of the origin; balance the two masters by setting the exponents equal and solving for $h$, then substitute $h^star$ back into either term to obtain the joint rate.]

#exercise[
  Use Theorem 34 of @Boyd2000 to derive the algebraic decay rate of Hermite coefficients of $1 \/ (1 + y^2)$.
] <ex-un-hermite-decay>

#exercise[
  Compute the first three $T B_n$ in explicit closed form for $ell = 1$, and verify that $T B_0 - T B_2 = 2 \/ (1 + y^2)$.
] <ex-un-tbn-explicit-verify>

#exercise[
  Show that the $T L_n$ map $y = ell (1 + x) \/ (1 - x)$ sends a pole of $f(y)$ at $y = i y_0$ to $x = (i y_0 - ell) \/ (i y_0 + ell)$. Compute the imaginary part of this Chebyshev-plane pole and deduce the asymptotic rate of convergence for $f(y) = 1 \/ (ell^2 + y^2)$.
] <ex-un-tln-pole>

#hint-for(<ex-un-tln-pole>)[Substitute $y = i y_0$ into the inverse map $x = (y - ell) \/ (y + ell)$ and verify that the image has modulus one; the geometric convergence rate is then governed by the Bernstein ellipse passing through that image point, and for $f(y) = 1 \/ (ell^2 + y^2)$ the pole sits at $y_0 = ell$.]

#exercise[
  Explain, in one paragraph, why Hermite functions are the "right" basis for the quantum harmonic oscillator but the "wrong" basis for $1 \/ (1 + y^2)$.
] <ex-un-hermite-right-wrong>

#exercise[
  Show that the asymptotic expansion of $J_0 (y)$ at infinity motivates the amplitude-phase decomposition of Étude 20.11.
] <ex-un-j0-asymptotic>

=== Computational Exercises

#exercise[
  Re-run Étude 20.1 with $f(y) = y e^(-y^2)$ in place of $sech(y)$. How do the three sweeps change?
] <ex-un-rerun-stalls>

#exercise[
  In Étude 20.4, add a fifth case $A = 0$ (the constant function $f = 1$). What happens? Why?
] <ex-un-hermite-constant>

#exercise[
  For Étude 20.5, add a third test function $f_3 (y) = (1 + y)^(-3 \/ 2)$. Which basis wins?
] <ex-un-laguerre-third>

#exercise[
  Implement the $T B_n$ expansion of Étude 20.6 with $ell in {0.5, 1, 2, 4}$. Confirm Boyd's Rule-of-Thumb 16 by inspecting $|a_n|$.
] <ex-un-tbn-rot16>

#exercise[
  Extend Étude 20.8 to the full coefficient comparison against Boyd's Table 17.7: with your $N = 21$ computation, reproduce the table to six decimal places.
] <ex-un-yoshida-table>

#exercise[
  Implement the Weideman--Cloot sinh mapping and reproduce Boyd's Fig 17.12 for $sech(y)$.
] <ex-un-weideman-cloot>

#exercise[
  Apply the composed-map strategy of Étude 20.10 to $r K_0 (r)$ instead of $r K_1 (r)$. The logarithmic singularity is now at the function itself, not multiplicatively removed. How does the method respond?
] <ex-un-rk0-composed>

=== Project-Style Exercises

#exercise(title: [The Radiation-Function Basis])[
  Study Boyd (1990b) on weakly-nonlocal solitary waves and implement a $T B_n$ expansion augmented with one or two radiation functions. Compare with a naive $T B_n$ expansion for a test nonlocal soliton.
] <ex-un-radiation-basis>

#hint-for(<ex-un-radiation-basis>)[A weakly nonlocal solitary wave radiates a small-amplitude, non-decaying oscillatory tail that no decaying rational series can capture; append one or two non-decaying radiation functions to the $T B_n$ set and let the tail amplitude become an extra unknown, in the spirit of the asymptotic augmentation of Étude 20.11.]

#exercise(title: [Composition with @ch-coord-transforms])[
  For a semi-infinite problem with an endpoint singularity AND an oscillatory tail, build a hybrid: a composed map (Étude 20.10) plus asymptotic augmentation (Étude 20.11). Discuss the conditioning of the resulting matrix.
] <ex-un-hybrid-compose>

#hint-for(<ex-un-hybrid-compose>)[Apply the two devices in sequence: first the singularity-resolving composed coordinate of Étude 20.10 to handle the endpoint, then factor the oscillation into $cos$ and $sin$ carriers as in Étude 20.11, expanding the slowly-varying amplitude and phase in the mapped variable; track the condition number as the carrier frequency and map parameter vary.]

#exercise(title: [The Decision Guide as a Program])[
  Write a small classifier that, given a sampled representation of a target function, diagnoses its asymptotic type and recommends a basis. Include tests on the four test cases of this chapter.
] <ex-un-decision-classifier>

#hint-for(<ex-un-decision-classifier>)[Estimate the tail from the sampled data: regress $log |f(y)|$ against $y$, $y^2$, and $log y$ over the outer decades to separate exponential, Gaussian, and algebraic decay, and inspect the envelope of $|f|$ for a non-decaying oscillation; map each diagnosis onto the decision guide of @sec-un-decision.]
