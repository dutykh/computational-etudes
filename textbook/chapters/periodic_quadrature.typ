// textbook/chapters/periodic_quadrature.typ
// Chapter 16: Integration of Periodic Functions: Why the Trapezoidal Rule Becomes Spectral
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Email: denys.dutykh@ku.ac.ae
// Homepage: https://www.denys-dutykh.com/
// Last modified: April 2026

#import "../styles/template.typ": dropcap, num, format-table

// Enable equation numbering for this chapter
#set math.equation(numbering: "(1)")

= Integration of Periodic Functions: Why the Trapezoidal Rule Becomes Spectral <ch-periodic-trap>

#dropcap[The previous chapter delivered an uncomfortable verdict: polynomial exactness is a deceptive measure of practical accuracy. We are now ready to inspect the other side of the same coin. If exactness in the wrong space can mislead us, can exactness in the _right_ space make a humble formula look spectacular? The answer, drawn from two beautiful papers of Trefethen and Weideman @TrefethenWeideman2014 @Weideman2002, is yes, and the example is one most readers think they already know: the equispaced trapezoidal rule. Judged by piecewise-linear interpolation, this is a clumsy first-order method: even Simpson's rule beats it. But on a periodic interval the right yardstick is _trigonometric_ exactness, not piecewise-linear exactness, and once we change yardsticks the trapezoidal rule reveals itself as a Fourier spectral method in disguise. For an analytic periodic integrand its error decays geometrically; for an entire periodic integrand it decays faster than any geometric rate; even for a $C^infinity$ integrand that is not analytic, it decays faster than any algebraic rate. Conversely, for periodic integrands of merely finite smoothness, it decays at a rate dictated entirely by how many odd derivatives match across the period, not by the function's pointwise smoothness. This chapter develops the resulting five-class taxonomy of convergence rates, traces it back to Poisson's 1820s computation of the perimeter of an ellipse @Davis1959, and lands us squarely inside the worldview of Fourier spectral methods built up in @ch-fourier-grids, @ch-spectral-pde and @ch-fourier-pseudo. The grand moral of @ch-quadrature returns: choose the approximation space first, and then let the quadrature follow.]

By the end of this chapter, you should be able to:

1. State the periodic trapezoidal rule on $[0, 2 pi)$ and recognise it as the exact integral of the trigonometric interpolant of the data.
2. Prove that the rule is exact for trigonometric polynomials of degree $lt.eq.slant N - 1$, using the aliasing identity for $e^(i k theta)$ on the equispaced grid.
3. Identify five convergence classes for periodic integrands -- band-limited (exact), finite smoothness (algebraic), strip analyticity (geometric), entire (supergeometric), and $C^infinity$-but-not-analytic (subgeometric) -- and predict the class of a given integrand from its analytic structure.
4. Derive the geometric error bound $|I_N - I| lt.eq.slant 4 pi M / (e^(a N) - 1)$ for functions analytic in a strip of half-width $a$, by both the Fourier-series-with-aliasing argument and the residue-calculus argument with the characteristic function $m(theta) = -frac(1, 2) cot(N theta \/ 2)$ @TrefethenWeideman2014.
5. Use Weideman's explicit error formulas to predict the convergence rate of the trapezoidal rule on $1\/(a - cos x)$, $e^(cos x)$, and $|sin(x\/2)|^k$, and verify these predictions numerically.
6. Translate the periodic theorems to the real line: explain why $integral_(-infinity)^infinity e^(-x^2) f(x) dif x$ is computed to machine precision by an equispaced trapezoidal sum with only a dozen nodes.
7. Compute the Fourier coefficients of a smooth periodic function using the Fast Fourier Transform and explain why this is just the periodic trapezoidal rule applied to the defining integral.
8. Articulate the doubled-rate observation: trapezoidal quadrature converges asymptotically twice as fast as the underlying trigonometric interpolant, because half of the aliased modes still integrate exactly to zero.

// ============================================================================
== From Polynomial to Trigonometric Exactness <sec-trig-exactness-bridge>
// ============================================================================

=== The textbook paradox

Open any introductory calculus or numerical analysis textbook to the section on numerical integration and you will find the following error bound for the composite trapezoidal rule on $[0, 2 pi]$ with $N$ equal subintervals:
$ |I(f) - T_N (f)| lt.eq.slant frac((2 pi)^3, 12 N^2) max_([0, 2 pi]) |f''(x)|. $ <eq-trap-textbook>
For the periodic function $f(x) = e^(cos x)$ this bound, with $K = max |f''| = e$, gives $|I - T_(10)(f)| lt.eq.slant pi^3 e \/ 300 approx 0.28$. The actual error at $N = 10$ is approximately $3.5 times 10^(-9)$. The textbook estimate is correct but useless: it overestimates the actual error by a factor of about $10^8$. Weideman called this "like saying the distance between New York and London is less than $10^(11)$ miles" @Weideman2002.

This is more than a quirk. It is the symptom of a misalignment between the theory the bound is built from and the structure the integrand actually has. The classical bound treats $f$ as a generic $C^2$ function on a closed interval and asks: how accurately does a piecewise-linear interpolant approximate it? On the periodic interval $[0, 2 pi]$ with a periodic integrand this is the wrong question. The right question, as we shall see, asks how accurately a _trigonometric_ interpolant approximates $f$, and the answer is breathtaking.

=== The right exactness space is trigonometric

Here is the punchline of the entire chapter, stated in advance so that the reader can hold it in mind as the theorems unfold:

#block(
  fill: rgb("#142D6E").lighten(92%),
  stroke: (left: 3pt + rgb("#142D6E")),
  inset: (left: 12pt, y: 10pt, right: 10pt),
  width: 100%,
)[
*Thesis.* The periodic trapezoidal rule
$ I_N (f) = frac(2 pi, N) sum_(j=0)^(N-1) f(2 pi j \/ N) $ <eq-periodic-trap-N>
is _exact_ for every trigonometric polynomial of degree $lt.eq.slant N - 1$. The right exactness space is therefore not piecewise-linear functions but trigonometric polynomials. Once one accepts that change of yardstick, every spectacular convergence result for smooth periodic integrands falls out as a corollary.
]

Note the contrast with the previous chapter. There we had to abandon polynomial exactness as a metric because monomials misled us about practical accuracy. Here we _embrace_ exactness in a different polynomial space (trigonometric, not algebraic) and find that it perfectly predicts the spectacular performance of the trapezoidal rule. The grand moral of @ch-quadrature is therefore not "exactness is bad" but rather _exactness in the wrong space_ is bad. In the right space, exactness is everything.

=== A historical aside: Poisson, 1823

The first person to notice the phenomenon empirically appears to have been Poisson, in the 1820s @Davis1959. He computed the perimeter of an ellipse with semi-axes $1\/(2 pi)$ and $0.6 / (2 pi)$:
$ I = frac(1, 2 pi) integral_0^(2 pi) sqrt(1 - 0.36 sin^2 theta) dif theta = frac(2, pi) E(0.36) approx 0.9027799277721857..., $ <eq-poisson-ellipse>
where $E$ is the complete elliptic integral of the second kind. Exploiting the four-fold symmetry of the integrand, Poisson computed the trapezoidal sum $I_N$ for $N = 4, 8, dots, 20$, with just three nontrivial function evaluations at $N = 16$, and obtained ten correct digits. He proved, by what was essentially the Euler--Maclaurin argument, that his estimate was in error by less than $4.84 times 10^(-6)$. He could not have known it, but the actual error decays as $|I_N - I| = cal(O)(3^(-N))$, since the integrand has branch points in the complex $theta$-plane at $theta = plus.minus i log(3)$. Each new sample point brings about $log_10 (3) approx 0.48$ correct digits. In our first étude we reproduce Poisson's calculation and visualise its geometric convergence.

// ============================================================================
== Computational Étude 16.1: Poisson's Ellipse, the Original Paradox <sec-etude-poisson>
// ============================================================================

This étude reproduces Figure 1.1 of @TrefethenWeideman2014. We integrate $sqrt(1 - 0.36 sin^2 theta) \/ (2 pi)$ on $[0, 2 pi]$, exploiting the four-fold symmetry, and plot the absolute error against $N \/ 4$ on a semilogarithmic scale. The reference value is $(2\/pi) E(0.36)$ from `scipy.special.ellipe`.

```python
import numpy as np
from scipy.special import ellipe

def trapezoidal_periodic(f, N):
    """N-point periodic trapezoidal rule on [0, 2*pi)."""
    theta = 2.0 * np.pi * np.arange(N) / N
    return (2.0 * np.pi / N) * np.sum(f(theta))

integrand = lambda t: np.sqrt(1.0 - 0.36 * np.sin(t)**2)
I_exact = (2.0 / np.pi) * ellipe(0.36)
for N in [4, 8, 12, 16, 20]:
    I_N = trapezoidal_periodic(integrand, N) / (2.0 * np.pi)
    print(N // 4, abs(I_N - I_exact))
```

The equivalent MATLAB implementation:

```matlab
integrand = @(t) sqrt(1 - 0.36 * sin(t).^2);
I_exact = (2 / pi) * ellipke(0.36);   % Note: MATLAB convention
for N = 4:4:20
    theta = 2*pi*(0:N-1)/N;
    I_N = (2*pi/N) * sum(integrand(theta)) / (2*pi);
    fprintf('%d  %.4e\n', N/4, abs(I_N - I_exact));
end
```

The Julia implementation:

```julia
using SpecialFunctions
integrand(t) = sqrt(1 - 0.36 * sin(t)^2)
I_exact = (2 / π) * ellipe(0.36)
for N in 4:4:20
    θ = 2π * (0:N-1) / N
    I_N = (2π/N) * sum(integrand.(θ)) / (2π)
    println(N ÷ 4, "  ", abs(I_N - I_exact))
end
```

Running the Python version reproduces (a corrected version of) Poisson's original table:

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    {
      show table: format-table(auto, none, auto)
      table(
        columns: 3,
        align: (center, center, center),
        inset: (x: 1em, y: 0.6em),
        stroke: none,
        table.hline(stroke: 0.75pt + rgb("#142D6E")),
        table.header(
          table.cell(fill: rgb("#142D6E").lighten(85%))[*$N\/4$*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*$I_N$*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*$|I_N - I|$*],
        ),
        table.hline(stroke: 0.5pt + luma(180)),
        [1], [$0.9000000000000$], num[2.78e-3],
        [2], [$0.9027692569069$], num[1.07e-5],
        [3], [$0.9027798586496$], num[6.91e-8],
        [4], [$0.9027799272276$], num[5.45e-10],
        [5], [$0.9027799277674$], num[4.76e-12],
        [6], [$0.9027799277721$], num[4.45e-14],
        [7], [$0.9027799277722$], num[5.55e-16],
      )
    },
  ),
  caption: [Trapezoidal-rule approximations to $(2\/pi) E(0.36)$, exploiting the four-fold symmetry of Poisson's integrand. Each new sample point reduces the error by a factor of approximately $3$, in agreement with the geometric rate $cal(O)(3^(-N))$ predicted by the strip-analyticity theorem of @sec-strip-analyticity.],
) <tbl-poisson-ellipse>

The convergence is shown graphically in @fig-poisson-ellipse, alongside the theoretical envelope $3^(-N)$.

#figure(
  image("../figures/ch16/python/poisson_ellipse.pdf", width: 80%),
  caption: [Geometric convergence of the periodic trapezoidal rule on Poisson's ellipse integrand. Each independent sample reduces the error by approximately $3$, until machine precision is reached at $N\/4 approx 8$. The dashed line is the theoretical envelope $3^(-N)$, derived in @sec-strip-analyticity from the location of the branch points of the integrand in the complex plane.],
) <fig-poisson-ellipse>

The code generating @fig-poisson-ellipse is available in:
- `codes/python/ch16/trap_poisson_ellipse.py`
- `codes/matlab/ch16/trap_poisson_ellipse.m`
- `codes/julia/ch16/trap_poisson_ellipse.jl`

// ============================================================================
== Trigonometric Exactness via Fourier Series and Aliasing <sec-trig-exactness>
// ============================================================================

We now develop the algebraic framework that explains Poisson's spectacular accuracy. Throughout this section we work on $[0, 2 pi]$ and assume $f: [0, 2 pi] arrow bb(C)$ is $2 pi$-periodic and Lipschitz continuous, so that its Fourier series
$ f(theta) = sum_(j = -infinity)^infinity c_j e^(i j theta), quad c_j = frac(1, 2 pi) integral_0^(2 pi) e^(-i j theta) f(theta) dif theta $ <eq-fourier-series>
converges uniformly and absolutely. The Fourier coefficient $c_0 = (1\/(2 pi)) integral_0^(2 pi) f(theta) dif theta$ is exactly the mean of $f$, and the integral we want is $I(f) = 2 pi c_0$.

=== The aliasing identity

The starting point is a one-line observation about the discrete sum of a single Fourier mode on the equispaced grid $theta_k = 2 pi k \/ N$:
$ sum_(k=0)^(N-1) e^(i j theta_k) = sum_(k=0)^(N-1) e^(2 pi i j k \/ N) = cases(N \, & j "is an integer multiple of" N\,, 0 \, & "otherwise.") $ <eq-trap-aliasing>
This is a finite geometric series; when $j$ is _not_ a multiple of $N$ the ratio $e^(2 pi i j \/ N) eq.not 1$ and the sum telescopes to zero, while when $j = ell N$ for some integer $ell$ each term equals $1$ and the sum is $N$. Geometrically, the $N$ phasors $e^(2 pi i j k \/ N)$ either coincide or distribute themselves uniformly around the unit circle and cancel.

The identity @eq-trap-aliasing is the entire content of "aliasing" in this chapter. On the $N$-point grid $theta_k = 2 pi k \/ N$, the modes
$ e^(i j theta), e^(i (j + N) theta), e^(i (j + 2 N) theta), dots, e^(i (j - N) theta), e^(i (j - 2 N) theta), dots $
are _indistinguishable_ as discrete sequences. We have already met this phenomenon in @ch-fourier-grids in the context of trigonometric interpolation; here we shall see that it works _in our favour_ for integration.

=== The Fourier-series proof

Multiplying @eq-trap-aliasing by $c_j$ and summing over all integers $j$ gives the exact discrete sum
$ frac(2 pi, N) sum_(k=0)^(N-1) f(theta_k) = 2 pi sum_(ell = -infinity)^infinity c_(ell N). $ <eq-trap-aliased-sum>
Subtracting $I(f) = 2 pi c_0$ produces the master error formula:
$ T_N (f) - I(f) = 2 pi sum_(ell eq.not 0) c_(ell N). $ <eq-trap-error>
Equation @eq-trap-error is the rosetta stone of this chapter. It says, in plain language, that the trapezoidal-rule error is exactly the sum of all _aliases_ of the constant Fourier mode, that is, the Fourier coefficients $c_j$ at frequencies $j = N, -N, 2 N, -2 N, dots$. The faster these "shaded" coefficients decay, the smaller the error.

=== The band-limited exactness theorem

The simplest possible consequence of @eq-trap-error is the exactness statement we need:

#block(
  fill: rgb("#142D6E").lighten(92%),
  stroke: (left: 3pt + rgb("#142D6E")),
  inset: (left: 12pt, y: 10pt, right: 10pt),
  width: 100%,
)[
*Theorem 1 (Band-limited exactness).* Let $f$ be a trigonometric polynomial of degree $lt.eq.slant n$, that is, a finite linear combination of $e^(i j theta)$ with $|j| lt.eq.slant n$. Then for every $N > n$, the periodic trapezoidal rule is exact:
$ T_N (f) = I(f). $ <eq-trig-exactness>
]

#text(weight: "bold")[Proof.] All Fourier coefficients $c_j$ with $|j| > n$ are zero. The sum on the right of @eq-trap-error therefore contains only the terms with $|ell N| lt.eq.slant n$, and the smallest such nonzero index is $|ell N| = N > n$, so the sum is empty. #h(1fr) $square$

This is the exactness statement that justifies our thesis. Trigonometric polynomials of degree $lt.eq.slant N - 1$ are integrated exactly by the $N$-point trapezoidal rule. Note the contrast with the Newton--Cotes case studied in @ch-quadrature: there the polynomial exactness degree was an algebraic invariant of the formula, here it is a _trigonometric_ exactness degree, and the formula is a trivially equispaced sum rather than a carefully designed quadrature rule.

There is also an unexpected dividend hidden in @eq-trap-error: even when $f$ is a single mode $e^(i k theta)$ with $k$ between $N\/2$ and $N$ (so $f$ is _not_ band-limited at the Nyquist rate of two points per wavelength), the trapezoidal rule still gives the correct answer of zero. The reason is that $e^(i k theta)$ aliases on the grid to $e^(i (k - N) theta)$, which has degree $|k - N| < N\/2$ and is therefore integrated exactly to zero. _Quadrature only needs one point per wavelength, even though reconstruction needs two._ This is the cleanest explanation we know for the apparent paradox of the "one point per wavelength condition for integration", and it is the seed of the doubled-rate observation that will close the chapter.

// ============================================================================
== Computational Étude 16.2: Band-Limited Exactness and One-Mode Aliasing <sec-etude-band>
// ============================================================================

This étude verifies @eq-trig-exactness experimentally and demonstrates the one-mode aliasing dividend. We construct a random trigonometric polynomial of degree $m = 10$ and check that the trapezoidal error becomes machine zero as soon as $N > m$. We then take the single mode $cos(k theta)$ with $N = 16$ fixed and sweep $k$ from $0$ to $32$, observing that the only failures occur at $k = 16, 32, dots$ (the integer multiples of $N$).

```python
import numpy as np
rng = np.random.default_rng(seed=42)
m = 10
a = rng.standard_normal(m + 1)  # cos coefficients
b = rng.standard_normal(m + 1)  # sin coefficients

def trig_poly(theta):
    s = np.full_like(theta, a[0])
    for k in range(1, m + 1):
        s += a[k] * np.cos(k * theta) + b[k] * np.sin(k * theta)
    return s

I_exact = 2.0 * np.pi * a[0]   # only the constant term contributes
for N in range(4, 31):
    err = abs(trapezoidal_periodic(trig_poly, N) - I_exact)
    print(N, err)
```

The equivalent MATLAB implementation:

```matlab
rng(42); m = 10;
a = randn(m + 1, 1); b = randn(m + 1, 1);
trig_poly = @(t) a(1) + sum(a(2:end) .* cos((1:m)' * t) + ...
                            b(2:end) .* sin((1:m)' * t), 1);
I_exact = 2 * pi * a(1);
for N = 4:30
    theta = 2*pi*(0:N-1)/N;
    err = abs((2*pi/N) * sum(trig_poly(theta)) - I_exact);
    fprintf('%d  %.2e\n', N, err);
end
```

The Julia implementation:

```julia
using Random; Random.seed!(42)
m = 10
a = randn(m + 1); b = randn(m + 1)
trig_poly(θ) = a[1] + sum(a[k+1] * cos(k*θ) + b[k+1] * sin(k*θ) for k in 1:m)
I_exact = 2π * a[1]
for N in 4:30
    θ = 2π * (0:N-1) / N
    err = abs((2π/N) * sum(trig_poly.(θ)) - I_exact)
    println(N, "  ", err)
end
```

@fig-band-limited shows both halves of the experiment. In panel (a) the error drops from $cal(O)(1)$ to machine precision at exactly $N = m + 1 = 11$ and stays there. In panel (b) the trapezoidal rule integrates $cos(k theta)$ to machine precision for every $k$ except $k = 16$ and $k = 32$, even though for $8 < k < 16$ the mode is aliased (it cannot be _reconstructed_ from the samples), and for $16 < k < 32$ it is doubly aliased.

#figure(
  image("../figures/ch16/python/band_limited.pdf", width: 95%),
  caption: [(a) Trapezoidal error for a random trigonometric polynomial of degree $m = 10$. Once $N > m$, the rule is exact. (b) For fixed $N = 16$, the trapezoidal error in $integral_0^(2 pi) cos(k theta) dif theta$ as a function of $k$. The only failures are at $k = N$ and $k = 2 N$ (integer multiples of $N$). Modes $k$ with $N\/2 < k < N$ are aliased on the grid but still integrate exactly to zero, because the one-point-per-wavelength condition is sufficient for quadrature even when it is not sufficient for reconstruction.],
) <fig-band-limited>

The code generating @fig-band-limited is available in:
- `codes/python/ch16/trap_band_limited.py`
- `codes/matlab/ch16/trap_band_limited.m`
- `codes/julia/ch16/trap_band_limited.jl`

// ============================================================================
== The Convergence Taxonomy <sec-taxonomy>
// ============================================================================

@eq-trap-error tells us that the trapezoidal-rule error is controlled entirely by the rate at which the Fourier coefficients of $f$ decay at high frequencies. We can therefore organise periodic integrands into a hierarchy of convergence classes, each defined by a different decay rate of the $c_j$:

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    table(
      columns: (1.6fr, 1.4fr, 1.6fr),
      align: (left, left, left),
      inset: (x: 0.9em, y: 0.55em),
      stroke: none,
      table.hline(stroke: 0.75pt + rgb("#142D6E")),
      table.header(
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Class of $f$*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Coefficient decay*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Trapezoidal error*],
      ),
      table.hline(stroke: 0.5pt + luma(180)),
      [Trigonometric polynomial of degree $n$], [finite support: $c_j = 0$ for $|j| > n$], [exact for all $N > n$],
      [$C^k$ periodic, but $C^(k+1)$ jumps], [$|c_j| = cal(O)(|j|^(-k-1))$], [$cal(O)(N^(-k-1))$ algebraic],
      [Analytic in a strip $|"Im" theta| < a$], [$|c_j| = cal(O)(e^(-a |j|))$], [$cal(O)(e^(-a N))$ geometric],
      [Entire periodic (e.g. $e^(cos x)$)], [$|c_j| = cal(O)((c\/|j|)^(|j|))$], [supergeometric: $cal(O)((e\/(2 N))^N)$],
      [$C^infinity$ periodic but not analytic], [decays faster than every algebraic rate but slower than every $r^j$], [subgeometric: e.g.\u{a0}$cal(O)(e^(-c N^(2\/3)))$],
    ),
  ),
  caption: [The five-class taxonomy of convergence rates for the periodic trapezoidal rule. Each class is defined by the decay rate of the Fourier coefficients of the integrand, and that decay rate is in turn determined by the analytic structure of $f$ in the complex plane.],
) <tbl-taxonomy>

The five classes are not just a bookkeeping device; they correspond to genuinely distinct asymptotic regimes that one can _see_ on a single log plot. Each of the next five sections develops one of them in detail and supplies a computational étude that pins down the corresponding rate experimentally on one of Weideman's named functions @Weideman2002.

// ============================================================================
== Algebraic Convergence: Regularity of the Periodic Extension <sec-algebraic>
// ============================================================================

We start at the bottom of the staircase. A periodic function of merely finite smoothness has Fourier coefficients that decay only algebraically; the trapezoidal-rule error therefore decays algebraically too. The exact rate is governed by the _smoothness of the periodic extension_ rather than by the smoothness of the function on its open period. This distinction is the source of the most common student misconception about the trapezoidal rule, and we shall hammer it.

=== The Euler--Maclaurin formula

The classical tool for understanding algebraic convergence is the Euler--Maclaurin formula, which relates the trapezoidal rule of an arbitrary (not necessarily periodic) function to its true integral via the values of its odd derivatives at the endpoints:
$ T_N (f) - I(f) = sum_(k=1)^m frac(B_(2 k), (2 k)!) h^(2 k) [ f^((2 k - 1))(2 pi) - f^((2 k - 1))(0) ] + R_(m, N), $ <eq-euler-maclaurin>
where $h = 2 pi \/ N$, $B_(2 k)$ are the Bernoulli numbers, and $R_(m, N)$ is a remainder of order $h^(2 m + 2)$. The crucial point is the bracketed term: the formula picks up a correction whenever an _odd_ derivative of $f$ fails to match across the period. If $f$ is genuinely $2 pi$-periodic and analytic, all odd derivatives match at $0$ and $2 pi$, every bracket vanishes, and the Euler--Maclaurin sum is empty: this is the first hint that smoothly periodic integrands are integrated faster than any power of $1 \/ N$.

But what if the periodic _extension_ has only limited regularity? Suppose, for instance, that $f$ is $C^k$-periodic and continuous, but its $(k+1)$-st derivative has a jump at $0$ (and hence at every multiple of $2 pi$). Then in @eq-euler-maclaurin all the bracketed terms with $2 ell - 1 lt.eq.slant k$ vanish, but the bracket with $2 ell - 1 = k + 1$ (if $k$ is odd) or higher does not, and the leading error is of order $h^(k + 2) tilde N^(-k-2)$ when $k$ is odd, or order $h^(k+1) tilde N^(-k-1)$ when $k$ is even. The convergence is therefore _algebraic_: the rate is set by the highest-order odd derivative that remains continuous across the period.

=== The lesson

The pedagogical lesson of this section is worth stating in a box, because it routinely surprises students:

#block(
  fill: rgb("#142D6E").lighten(92%),
  stroke: (left: 3pt + rgb("#142D6E")),
  inset: (left: 12pt, y: 10pt, right: 10pt),
  width: 100%,
)[
*The dividing line is not "smooth" versus "non-smooth".* The dividing line is whether the periodic _extension_ is smooth. A function that is infinitely differentiable on $(0, 2 pi)$ but whose periodic extension has a jump in its first derivative at $0$ is integrated only at the algebraic rate $cal(O)(N^(-2))$, whereas a function that is $C^infinity$-periodic but not analytic is integrated at a subgeometric rate that can be much faster.
]

We test both halves of this statement in the next étude, using two of Weideman's classic examples.

// ============================================================================
== Computational Étude 16.3: Algebraic Decay on $|sin(x\/2)|^k$ <sec-etude-algebraic>
// ============================================================================

The two test functions are
$ f_2 (x) = abs(sin(x \/ 2)), quad f_3 (x) = abs(sin(x \/ 2))^3, $
both regarded as $2 pi$-periodic functions of $x$. Their exact integrals are $I_2 = 4$ and $I_3 = 8 \/ 3$ respectively. The function $f_2$ has a jump in its first derivative at $x = 0$ (the periodic extension is the absolute value $|sin(x\/2)|$, which has a corner there); $f_3$ has continuous first and second derivatives but a jump in its third derivative. From the Euler--Maclaurin analysis we therefore expect $cal(O)(N^(-2))$ convergence for $f_2$ and $cal(O)(N^(-4))$ for $f_3$. Weideman @Weideman2002 derives the explicit asymptotic constants
$ I(f_2) - T_N (f_2) tilde -frac(pi^2, 3) frac(1, N^2), quad I(f_3) - T_N (f_3) tilde -frac(pi^4, 30) frac(1, N^4). $ <eq-algebraic-rates>

```python
f2 = lambda x: np.abs(np.sin(x / 2))
f3 = lambda x: np.abs(np.sin(x / 2))**3
N_values = np.array([8, 16, 32, 64, 128, 256, 512, 1024, 2048])
err2, err3 = [], []
for N in N_values:
    err2.append(abs(trapezoidal_periodic(f2, N) - 4.0))
    err3.append(abs(trapezoidal_periodic(f3, N) - 8.0 / 3.0))
slope2 = np.polyfit(np.log10(N_values[2:]), np.log10(err2[2:]), 1)[0]
slope3 = np.polyfit(np.log10(N_values[2:]), np.log10(err3[2:]), 1)[0]
print(f"f2 slope: {slope2:.3f} (theory -2)")
print(f"f3 slope: {slope3:.3f} (theory -4)")
```

The equivalent MATLAB implementation:

```matlab
f2 = @(x) abs(sin(x/2));
f3 = @(x) abs(sin(x/2)).^3;
N_values = 2.^(3:11);
for k = 1:length(N_values)
    N = N_values(k); theta = 2*pi*(0:N-1)/N;
    err2(k) = abs((2*pi/N)*sum(f2(theta)) - 4);
    err3(k) = abs((2*pi/N)*sum(f3(theta)) - 8/3);
end
p2 = polyfit(log10(N_values(3:end)), log10(err2(3:end)), 1);
p3 = polyfit(log10(N_values(3:end)), log10(err3(3:end)), 1);
fprintf('f2 slope: %.3f, f3 slope: %.3f\n', p2(1), p3(1));
```

The Julia implementation:

```julia
f2(x) = abs(sin(x/2))
f3(x) = abs(sin(x/2))^3
N_values = [8, 16, 32, 64, 128, 256, 512, 1024, 2048]
err2 = [abs((2π/N)*sum(f2.(2π*(0:N-1)/N)) - 4)        for N in N_values]
err3 = [abs((2π/N)*sum(f3.(2π*(0:N-1)/N)) - 8/3)      for N in N_values]
# Linear fit gives slopes ~ -2 and -4
```

The numerical slopes from a linear fit through the asymptotic tail are $-2.000$ and $-4.001$, confirming the predicted rates to three digits. @fig-algebraic-decay shows the two convergence curves on a log--log plot.

#figure(
  image("../figures/ch16/python/algebraic_decay.pdf", width: 80%),
  caption: [Algebraic convergence on the periodic functions $f_2(x) = |sin(x\/2)|$ (slope $-2$, coral) and $f_3(x) = |sin(x\/2)|^3$ (slope $-4$, navy). The dashed envelopes are Weideman's explicit asymptotic constants $pi^2 \/ (3 N^2)$ and $pi^4 \/ (30 N^4)$ from @eq-algebraic-rates. The two slopes differ by a factor of two because $f_3$ has two more matching odd derivatives across the period than $f_2$.],
) <fig-algebraic-decay>

The code generating @fig-algebraic-decay is available in:
- `codes/python/ch16/trap_algebraic_decay.py`
- `codes/matlab/ch16/trap_algebraic_decay.m`
- `codes/julia/ch16/trap_algebraic_decay.jl`

// ============================================================================
== Geometric Convergence: Analyticity in a Strip <sec-strip-analyticity>
// ============================================================================

We now climb up the staircase to the most consequential rung: geometric convergence. The result of this section is the central theorem of the chapter.

=== Statement of the theorem

#block(
  fill: rgb("#142D6E").lighten(92%),
  stroke: (left: 3pt + rgb("#142D6E")),
  inset: (left: 12pt, y: 10pt, right: 10pt),
  width: 100%,
)[
*Theorem 2 (Strip-analyticity, Trefethen--Weideman 2014).* Let $f$ be $2 pi$-periodic, analytic, and bounded by $M$ in the open strip $abs("Im" theta) < a$ for some $a > 0$. Then for every $N gt.eq.slant 1$,
$ abs(I_N - I) lt.eq.slant frac(4 pi M, e^(a N) - 1) tilde 4 pi M e^(-a N), quad N arrow infinity, $ <eq-strip-thm>
and the constant $4 pi$ cannot be improved. The error therefore decays geometrically at rate $e^(-a N)$, where $a$ is the half-width of the maximal strip of analyticity of $f$.
]

The constant $a$ is the only quantity the practitioner needs to estimate from $f$. It is the distance from the real axis to the nearest singularity of $f$ in the complex $theta$-plane. For Poisson's ellipse integrand, $sqrt(1 - 0.36 sin^2 theta)$ has branch points where $sin^2 theta = 1\/0.36$, i.e. at $theta = plus.minus i log(3)$, so $a = log(3)$ and the predicted rate is $e^(-N log(3)) = 3^(-N)$. This is exactly the rate observed in @fig-poisson-ellipse.

We give two proofs of @eq-strip-thm. They illustrate two different argument styles, both useful elsewhere in spectral methods.

=== Proof by Fourier series and aliasing

This is the natural completion of the master error formula @eq-trap-error. The hypothesis that $f$ is bounded by $M$ in the strip $abs("Im" theta) < a$ allows us to estimate the Fourier coefficients of $f$ via Cauchy integrals on horizontal contours in the strip. Shifting the integration contour from the real axis $[0, 2 pi]$ down to $[0, 2 pi] - i a'$ (with $a' < a$ arbitrarily close to $a$), the periodicity of $f$ ensures that the contributions from the vertical sides of the rectangle cancel, and we are left with
$ c_j = frac(1, 2 pi) integral_(-i a' )^(2 pi - i a') e^(-i j theta) f(theta) dif theta. $
On this contour, $abs(e^(-i j theta)) = e^(-j a')$ for $j > 0$, so
$ abs(c_j) lt.eq.slant frac(1, 2 pi) dot 2 pi M e^(-j a') = M e^(-j a'). $
A symmetric argument with the contour shifted up gives the same bound for $j < 0$. Letting $a' arrow a$ and substituting into @eq-trap-error,
$ abs(I_N - I) = abs(2 pi sum_(ell eq.not 0) c_(ell N)) lt.eq.slant 2 pi sum_(ell eq.not 0) M e^(-|ell| N a) = 4 pi M sum_(ell = 1)^infinity e^(-ell N a) = frac(4 pi M, e^(a N) - 1), $
which is precisely @eq-strip-thm. The constant $4 pi$ is the right one because the leading term of the geometric series is exactly $4 pi M e^(-a N)$, and a one-line check on $f(theta) = e^(i N theta)$ shows that this leading constant cannot be reduced. #h(1fr) $square$

=== Proof by residue calculus

The second proof, due in essence to Davis @Davis1959, uses the _characteristic function_
$ m(theta) = -frac(1, 2) cot(N theta \/ 2) = -frac(1, 2) frac(1 + e^(-i N theta), 1 - e^(-i N theta)), $ <eq-cot-char>
which has simple poles at the equispaced points $theta_k = 2 pi k \/ N$ with residues equal to $-i \/ N$. By residue calculus, for any positively oriented contour $Gamma$ enclosing exactly the poles in $(0, 2 pi)$, we have $T_N (f) = integral_Gamma m(theta) f(theta) dif theta$. Choosing $Gamma$ to be the rectangle with vertices $pi \/ N - i a'$, $2 pi + pi \/ N - i a'$, $2 pi + pi \/ N + i a'$, $pi \/ N + i a'$, and using the fact that the contributions from the vertical sides cancel by periodicity, one obtains an explicit integral representation for $T_N (f) - I(f)$ that immediately yields the bound $|I_N - I| lt.eq.slant 4 pi M \/ (e^(a N) - 1)$. We omit the algebraic details (they can be found in @TrefethenWeideman2014, eqs. (3.19)--(3.23)), but emphasise that this proof never invokes Fourier series at all: the trapezoidal rule is recast as a contour integral, and the analyticity of $f$ in the strip lets us deform the contour into the half-planes where the cotangent is exponentially small. The same trick will reappear in spectral methods whenever we want to bound discrete sums by complex contour integrals.

=== Why $a$ matters

The half-width $a$ is the only handle the user has on the convergence rate of @eq-strip-thm. For functions with poles, $a$ is simply the imaginary part of the nearest pole to the real axis. For functions with branch points, $a$ is the imaginary part of the nearest branch point. For functions like $1\/(a - cos x)$ with $a > 1$, the singularities lie at $cos theta = a$, i.e. at $theta = plus.minus i log(a + sqrt(a^2 - 1))$, so $a_("strip") = log(a + sqrt(a^2 - 1))$ and the predicted rate is
$ e^(-N log(a + sqrt(a^2 - 1))) = (a + sqrt(a^2 - 1))^(-N) = (a - sqrt(a^2 - 1))^N, $
the second equality coming from $(a + sqrt(a^2 - 1))(a - sqrt(a^2 - 1)) = 1$. The next étude tests this prediction directly.

// ============================================================================
== Computational Étude 16.4: Geometric Convergence on the Poisson Kernel <sec-etude-poisson-kernel>
// ============================================================================

We take $f_4 (x) = 1 \/ (a - cos x)$ with $a = 2$. Its exact integral is $I = 2 pi \/ sqrt(3) approx 3.6275987284684357$, and the predicted geometric rate is $r^N$ with $r = 2 - sqrt(3) approx 0.2679$. By a partial-fraction expansion (Weideman 2002, eqs. (16)--(21)), one can derive the _exact_ trapezoidal error in closed form:
$ I(f_4) - T_N (f_4) = -8 pi frac(r, 1 - r^2) frac(r^N, 1 - r^N), $ <eq-poisson-kernel-exact>
which is the strict analogue of the strip bound @eq-strip-thm and even pins down the asymptotic constant.

```python
a = 2.0
f4 = lambda x: 1.0 / (a - np.cos(x))
I_exact = 2.0 * np.pi / np.sqrt(a**2 - 1.0)
r = a - np.sqrt(a**2 - 1.0)
for N in range(2, 31, 2):
    err = abs(trapezoidal_periodic(f4, N) - I_exact)
    theory = 8 * np.pi * r/(1 - r**2) * r**N / (1 - r**N)
    print(N, err, theory)
```

The equivalent MATLAB implementation:

```matlab
a = 2; f4 = @(x) 1 ./ (a - cos(x));
I_exact = 2*pi / sqrt(a^2 - 1);
r = a - sqrt(a^2 - 1);
for N = 2:2:30
    theta = 2*pi*(0:N-1)/N;
    err = abs((2*pi/N) * sum(f4(theta)) - I_exact);
    theory = 8*pi * r/(1 - r^2) * r^N / (1 - r^N);
    fprintf('%d  %.4e  %.4e\n', N, err, theory);
end
```

The Julia implementation:

```julia
a = 2.0; f4(x) = 1 / (a - cos(x))
I_exact = 2π / sqrt(a^2 - 1)
r = a - sqrt(a^2 - 1)
for N in 2:2:30
    θ = 2π * (0:N-1) / N
    err = abs((2π/N) * sum(f4.(θ)) - I_exact)
    theory = 8π * r/(1 - r^2) * r^N / (1 - r^N)
    println(N, "  ", err, "  ", theory)
end
```

@fig-poisson-kernel shows the result alongside a complex-plane picture of the singularities. The numerical and analytic curves agree to all displayed digits down to machine precision around $N approx 30$. The right-hand panel makes the relationship between the singularity geometry and the convergence rate explicit: the strip of analyticity is bounded by the horizontal lines $"Im" theta = plus.minus log(2 + sqrt(3))$, which contain the trapezoidal grid (red dots).

#figure(
  image("../figures/ch16/python/poisson_kernel.pdf", width: 95%),
  caption: [(a) Geometric decay of the trapezoidal error for $f_4(x) = 1 \/ (2 - cos x)$. The numerical errors (coral) and Weideman's exact closed-form error @eq-poisson-kernel-exact (navy dashed) agree to machine precision. (b) The integrand has poles in the complex plane at $theta = plus.minus i log(2 + sqrt(3)) = plus.minus 1.317i$ (and their $2 pi$-translates). The shaded region is the strip of analyticity; the trapezoidal nodes (coral dots) lie on the real axis. Each new node reduces the error by a factor of $2 - sqrt(3) approx 0.27$.],
) <fig-poisson-kernel>

The code generating @fig-poisson-kernel is available in:
- `codes/python/ch16/trap_poisson_kernel.py`
- `codes/matlab/ch16/trap_poisson_kernel.m`
- `codes/julia/ch16/trap_poisson_kernel.jl`

// ============================================================================
== Supergeometric Convergence: Entire Periodic Functions <sec-supergeometric>
// ============================================================================

If the integrand is analytic _everywhere_ in the complex plane, the strip bound @eq-strip-thm applies for every value of $a > 0$. The natural question is then how to choose $a$ to minimise the error. For an entire periodic function $f$ that grows in the imaginary direction like $e^(M(a))$ for some increasing function $M$, the optimal $a$ grows with $N$ and the resulting error decays _faster than any geometric rate_. This regime is called _supergeometric_ convergence.

The cleanest example is $f_5 (x) = e^(cos x)$, whose maximum on the strip $abs("Im" theta) < a$ is $e^(cosh a)$. Substituting into @eq-strip-thm,
$ abs(I_N - I) lt.eq.slant frac(4 pi e^(cosh a), e^(a N) - 1). $
A short calculus exercise shows that this bound is minimised by choosing $a$ close to $log(2 N)$, giving $cosh(a) approx N$ and
$ abs(I_N - I) lt.eq.slant frac(4 pi e^N, (2 N)^N - 1) tilde 4 pi (e \/ (2 N))^N, quad N arrow infinity. $ <eq-supergeometric-rate>
Each new node not only divides the error by a constant; it divides it by an _ever-growing_ factor, since $(e \/ (2 N))^N$ shrinks faster than any geometric rate $r^N$. The accuracy gain per sample point grows without bound.

We can see the same conclusion through the lens of the Bessel-function expansion of $e^(cos x)$:
$ e^(cos x) = I_0 (1) + 2 sum_(k=1)^infinity I_k (1) cos(k x), $
where $I_k$ is the modified Bessel function of the first kind. The exact integral is $I = 2 pi I_0 (1) approx 7.954926521012846$. By the master error formula @eq-trap-error and the symmetry $c_(-k) = c_k$,
$ I(f_5) - T_N (f_5) = -4 pi sum_(ell=1)^infinity (-1)^ell I_(ell N) (1) tilde -4 pi I_N (1) tilde 2 sqrt(2 pi \/ N) (e \/ (2 N))^N, $
where the asymptotics for $I_N (1)$ as $N arrow infinity$ come from the standard expansion of the modified Bessel function for large index. This is exactly @eq-supergeometric-rate up to the slowly varying prefactor $sqrt(2 pi \/ N)$.

// ============================================================================
== Computational Étude 16.5: Supergeometric Decay on $e^(cos x)$ <sec-etude-supergeometric>
// ============================================================================

We reproduce the convergence table from Section 4 of @TrefethenWeideman2014. The exact value is $2 pi I_0 (1)$, available from `scipy.special.iv(0, 1)` in Python.

```python
from scipy.special import iv
f5 = lambda x: np.exp(np.cos(x))
I_exact = 2.0 * np.pi * iv(0, 1.0)
for N in range(1, 17):
    I_N = trapezoidal_periodic(f5, N)
    err = abs(I_N - I_exact)
    theory = 2.0 * np.sqrt(2.0 * np.pi / N) * (np.e / (2.0 * N))**N
    print(N, I_N, err, theory)
```

The equivalent MATLAB implementation:

```matlab
f5 = @(x) exp(cos(x));
I_exact = 2 * pi * besseli(0, 1);
for N = 1:16
    theta = 2*pi*(0:N-1)/N;
    I_N = (2*pi/N) * sum(f5(theta));
    err = abs(I_N - I_exact);
    theory = 2 * sqrt(2*pi/N) * (exp(1)/(2*N))^N;
    fprintf('%d  %.16f  %.4e  %.4e\n', N, I_N, err, theory);
end
```

The Julia implementation:

```julia
using SpecialFunctions
f5(x) = exp(cos(x))
I_exact = 2π * besseli(0, 1)
for N in 1:16
    θ = 2π * (0:N-1) / N
    I_N = (2π/N) * sum(f5.(θ))
    err = abs(I_N - I_exact)
    theory = 2 * sqrt(2π/N) * (exp(1)/(2N))^N
    println(N, "  ", I_N, "  ", err, "  ", theory)
end
```

The numerical results agree with the asymptotic envelope @eq-supergeometric-rate to within a factor of two for every $N$ from $4$ onwards. By $N = 14$ the trapezoidal sum agrees with $2 pi I_0 (1)$ to roughly $14$ decimal digits, that is, almost all of double precision. @fig-supergeometric shows the error against $N$ on a semilogarithmic scale; note how the slope of the curve _steepens_ as $N$ grows, the diagnostic of supergeometric (faster than any geometric) decay.

#figure(
  image("../figures/ch16/python/supergeometric.pdf", width: 80%),
  caption: [Supergeometric decay of the trapezoidal-rule error for $f_5 (x) = e^(cos x)$. The dashed envelope is the Bessel-function asymptotic $2 sqrt(2 pi \/ N) (e \/ (2 N))^N$ derived above. The slope of the convergence curve becomes steeper with each new node, the signature of a faster-than-geometric decay; by $N = 14$ the trapezoidal sum agrees with $2 pi I_0 (1)$ to almost double-precision accuracy.],
) <fig-supergeometric>

The code generating @fig-supergeometric is available in:
- `codes/python/ch16/trap_supergeometric.py`
- `codes/matlab/ch16/trap_supergeometric.m`
- `codes/julia/ch16/trap_supergeometric.jl`

// ============================================================================
== Subgeometric Convergence: $C^infinity$ but Not Analytic <sec-subgeometric>
// ============================================================================

We have seen that band-limited integrands give exact answers, that integrands of finite smoothness give algebraic convergence, and that analytic-in-a-strip integrands give geometric or supergeometric convergence. The remaining slot in the taxonomy is curious: a function that is _smooth_ in the sense of $C^infinity$, but not analytic, must converge faster than any algebraic rate (because the Euler--Maclaurin sum @eq-euler-maclaurin terminates at no finite order) but cannot reach geometric convergence (because the Cauchy estimate of the Fourier coefficients fails as soon as the strip of analyticity collapses to the real line). The trapezoidal rule must therefore enter a regime that is genuinely intermediate: faster than $1\/N^k$ for every $k$, but slower than $r^N$ for every $r < 1$. Such intermediate rates are called _subgeometric_, and they look like
$ abs(I_N - I) tilde C exp(-c N^alpha) $
for some $0 < alpha < 1$.

Weideman @Weideman2002 supplies an explicit example: the function
$ f_6 (x) = cases(exp((cos x - 1) \/ (cos x + 1)) \, & 0 lt.eq.slant x lt.eq.slant 2 pi\, x eq.not pi\,, 0 \, & x = pi\,) $ <eq-f6-def>
which is $2 pi$-periodic and infinitely differentiable, but has an essential singularity at $x = pi$ and therefore lies in no strip of analyticity. By a careful analysis based on the Meijer $G$-function representation of the Fourier coefficients, Weideman proves that the trapezoidal error has the asymptotic envelope
$ abs(I (f_6) - T_N (f_6)) tilde 8 sqrt(pi \/ 3) exp(-frac(3, 2) N^(2 \/ 3)). $ <eq-subgeometric-rate>
The exponent $alpha = 2\/3$ is the diagnostic for the analytic structure of $f_6$. The exact integral is $I(f_6) = 2 e pi (1 - "erf"(1)) approx 2.6865868432934716$, computable from `scipy.special.erf`.

The visual diagnostic for subgeometric decay is to plot $log |I_N - I|$ against $N^(2 \/ 3)$ instead of $N$: if the rate is right, the curve becomes a straight line. This is the visual pun the étude is built around.

// ============================================================================
== Computational Étude 16.6: Subgeometric Decay on Weideman's $f_6$ <sec-etude-subgeometric>
// ============================================================================

```python
from scipy.special import erf
def f6(x):
    out = np.zeros_like(x)
    mask = (np.cos(x) + 1.0) > 1e-15
    out[mask] = np.exp((np.cos(x[mask]) - 1.0) / (np.cos(x[mask]) + 1.0))
    return out

I_exact = 2.0 * np.exp(1.0) * np.pi * (1.0 - erf(1.0))
N_vals = np.array([4, 6, 8, 10, 12, 16, 20, 24, 30, 40, 50, 60, 80, 100, 120, 160, 200])
errors = np.array([abs(trapezoidal_periodic(f6, N) - I_exact) for N in N_vals])
# Diagnostic plot: error vs N^{2/3}
import matplotlib.pyplot as plt
plt.semilogy(N_vals**(2/3), errors, 'o-')   # should be a straight line
```

The equivalent MATLAB implementation:

```matlab
f6 = @(x) (cos(x)+1 > 1e-15) .* exp((cos(x)-1) ./ max(cos(x)+1, 1e-15));
I_exact = 2 * exp(1) * pi * (1 - erf(1));
N_vals = [4 6 8 10 12 16 20 24 30 40 50 60 80 100 120 160 200];
for k = 1:length(N_vals)
    N = N_vals(k); theta = 2*pi*(0:N-1)/N;
    err(k) = abs((2*pi/N) * sum(f6(theta)) - I_exact);
end
semilogy(N_vals.^(2/3), err, 'o-')
```

The Julia implementation:

```julia
using SpecialFunctions
function f6(x)
    c = cos(x)
    (c + 1) > 1e-15 ? exp((c - 1) / (c + 1)) : 0.0
end
I_exact = 2 * exp(1) * π * (1 - erf(1))
N_vals = [4, 6, 8, 10, 12, 16, 20, 24, 30, 40, 50, 60, 80, 100, 120, 160, 200]
errors = [abs((2π/N)*sum(f6.(2π*(0:N-1)/N)) - I_exact) for N in N_vals]
```

@fig-subgeometric shows both diagnostic plots. In panel (a) the convergence curve plotted against $N$ looks suspicious: it bends downward but never settles into a straight line. In panel (b) the same data plotted against $N^(2 \/ 3)$ becomes nearly linear (the small wiggle is the oscillatory cosine factor that Weideman derives explicitly), confirming the rate $exp(-frac(3, 2) N^(2 \/ 3))$. The pedagogical point is precisely the visual transformation: by replotting the data on the right axis we can _read off_ the analytic structure of the integrand.

#figure(
  image("../figures/ch16/python/subgeometric.pdf", width: 95%),
  caption: [Subgeometric decay of the trapezoidal-rule error for the $C^infinity$-but-not-analytic function $f_6 (x) = exp((cos x - 1)\/(cos x + 1))$. (a) Error vs $N$: the curve bends downward but never becomes a straight line, in contrast to the geometric and supergeometric examples. (b) Error vs $N^(2\/3)$: the curve becomes nearly straight, confirming the predicted rate $exp(-frac(3, 2) N^(2\/3))$. The visual transformation between the two panels diagnoses the analytic structure of the integrand from the convergence data alone.],
) <fig-subgeometric>

The code generating @fig-subgeometric is available in:
- `codes/python/ch16/trap_subgeometric.py`
- `codes/matlab/ch16/trap_subgeometric.m`
- `codes/julia/ch16/trap_subgeometric.jl`

// ============================================================================
== The Trapezoidal Rule on the Real Line <sec-real-line>
// ============================================================================

The five-class taxonomy of @sec-taxonomy was developed for periodic integrands. The same arguments apply, with one substitution, to integrals over the entire real line. The substitution is to replace Fourier _series_ by Fourier _transforms_ and to take the limit $N arrow infinity$ in such a way that the grid spacing $h = 2 pi \/ N$ approaches zero.

Concretely, for a function $w$ on $bb(R)$ that decays sufficiently fast at $plus.minus infinity$, the (infinite) trapezoidal rule with step $h$ is
$ I_h = h sum_(k = -infinity)^infinity w(k h). $ <eq-real-line-trap>
Theorem 5.1 of @TrefethenWeideman2014 then says that if $w$ is analytic in the strip $abs("Im" z) < a$, decays uniformly to zero as $abs(x) arrow infinity$ in that strip, and has $L^1$ bound $M$ on every horizontal line $z = x + i b$ with $b in (-a, a)$, then
$ abs(I_h - I) lt.eq.slant frac(2 M, e^(2 pi a \/ h) - 1). $ <eq-real-line-thm>
The error decays as $cal(O)(e^(-2 pi a \/ h))$ as $h arrow 0$, which is the natural transcription of the periodic bound @eq-strip-thm with $N$ replaced by $2 pi \/ h$. For the famous example of Goodwin @Goodwin1949, $w(x) = e^(-x^2) \/ sqrt(pi)$, the integrand is entire (so $a$ can be taken arbitrarily large) and the error reaches machine precision at $h$ as large as $pi \/ 6$, that is, with only about $25$ nontrivial samples.

This result is the bridge from the periodic theory to all of the unbounded-domain trapezoidal-rule applications surveyed in Sections 5--7 of @TrefethenWeideman2014: inverse Laplace transforms via Hankel contours, special-function evaluation, and the double-exponential transformations of Takahasi--Mori and Stenger. We do not pursue these applications here, but the next étude makes the basic spectral convergence visible on the canonical Goodwin example.

// ============================================================================
== Computational Étude 16.7: Real-Line Trapezoidal Rule on the Gaussian <sec-etude-real-line>
// ============================================================================

We compute $integral_(-infinity)^infinity e^(-x^2) \/ sqrt(pi) dif x = 1$ by truncating the infinite sum @eq-real-line-trap to $|k| lt.eq.slant n$ with $n$ chosen large enough that the missing tail is below $10^(-300)$. We then plot the absolute error as a function of $h = 2 pi \/ N$ for $N = 1, 2, dots, 12$.

```python
def trapezoidal_real_line(w, h, n_max):
    k = np.arange(-n_max, n_max + 1)
    return h * np.sum(w(k * h))

w = lambda x: np.exp(-x**2) / np.sqrt(np.pi)
for N in range(1, 13):
    h = 2.0 * np.pi / N
    n_max = max(int(np.ceil(28.0 / h)), 30)   # ensure tail < 1e-300
    I_h = trapezoidal_real_line(w, h, n_max)
    print(N, h, I_h, abs(I_h - 1.0))
```

The equivalent MATLAB implementation:

```matlab
w = @(x) exp(-x.^2) / sqrt(pi);
for N = 1:12
    h = 2 * pi / N;
    n_max = max(ceil(28 / h), 30);
    k = -n_max:n_max;
    I_h = h * sum(w(k * h));
    fprintf('%d  %.6f  %.16f  %.4e\n', N, h, I_h, abs(I_h - 1));
end
```

The Julia implementation:

```julia
w(x) = exp(-x^2) / sqrt(π)
for N in 1:12
    h = 2π / N
    n_max = max(Int(ceil(28 / h)), 30)
    k = -n_max:n_max
    I_h = h * sum(w.(k * h))
    println(N, "  ", h, "  ", I_h, "  ", abs(I_h - 1))
end
```

The convergence is shown in @fig-real-line-gaussian. The semilogarithmic slope steepens noticeably as $h$ shrinks, reflecting the fact that the entire Gaussian admits arbitrarily large strip widths $a$, so the rate constant $2 pi a \/ h$ in the bound @eq-real-line-thm grows along the way. By $N = 12$ (i.e. $h = pi \/ 6$), the error is below $10^(-15)$.

#figure(
  image("../figures/ch16/python/real_line_gaussian.pdf", width: 80%),
  caption: [The trapezoidal rule on the real line, applied to $integral_(-infinity)^infinity e^(-x^2) \/ sqrt(pi) dif x = 1$. The reference dashed line is the simplest envelope $e^(-pi^2 \/ h)$ obtained by taking $a = pi \/ 2$ in the strip bound @eq-real-line-thm. The actual error decays even faster because the Gaussian is entire, allowing the optimal strip width $a$ to grow with $N$. By $N = 12$ the error is at machine precision.],
) <fig-real-line-gaussian>

The code generating @fig-real-line-gaussian is available in:
- `codes/python/ch16/trap_real_line.py`
- `codes/matlab/ch16/trap_real_line.m`
- `codes/julia/ch16/trap_real_line.jl`

// ============================================================================
== Computing Fourier Coefficients with the FFT <sec-fft>
// ============================================================================

We close the technical content of the chapter with the most consequential application of all. The Fourier coefficients of a $2 pi$-periodic function are themselves periodic integrals,
$ c_k = frac(1, 2 pi) integral_0^(2 pi) e^(-i k theta) f(theta) dif theta, $ <eq-fourier-def>
and the natural numerical estimator is the periodic trapezoidal rule applied to this integral:
$ hat(c)_k = frac(1, N) sum_(j=0)^(N-1) f(theta_j) e^(-2 pi i k j \/ N). $ <eq-fft-coeffs>
But the right-hand side of @eq-fft-coeffs is exactly the discrete Fourier transform (DFT) of the sample vector $(f(theta_0), dots, f(theta_(N-1)))$, divided by $N$. _The trapezoidal rule applied to the Fourier-coefficient integral is the FFT, up to a normalisation constant._

Every spectral convergence theorem we have proved in this chapter therefore transfers, for free, to the FFT. If $f$ is analytic in a strip of half-width $a$, then the FFT-computed coefficient $hat(c)_k$ satisfies
$ abs(hat(c)_k - c_k) = cal(O)(e^(-a N)), quad N arrow infinity, "for each fixed " k, $
and the entire vector of $N$ coefficients can be computed in $cal(O)(N log N)$ operations. This is the algorithmic foundation of every Fourier pseudospectral method we developed in @ch-spectral-pde, @ch-fourier-pseudo and @ch-polar. Whenever we wrote `np.fft.fft(u)/N` in those chapters to extract Fourier coefficients from samples, what we were _really_ doing was applying the periodic trapezoidal rule to the integral @eq-fourier-def. The spectacular accuracy of Fourier spectral discretisation is, at bottom, the spectacular accuracy of the periodic trapezoidal rule.

// ============================================================================
== Computational Étude 16.8: FFT Computation of Fourier Coefficients <sec-etude-fft>
// ============================================================================

We test @eq-fft-coeffs on the function $f(x) = 1\/(2 - cos x)$ from @sec-etude-poisson-kernel, whose Fourier coefficients are known in closed form:
$ c_k = frac(1, sqrt(3)) (2 - sqrt(3))^(|k|), quad k in bb(Z). $ <eq-poisson-coefficients>
We sample $f$ on $N = 32$ equispaced points, take the FFT, divide by $N$, and compare with the closed-form $c_k$.

```python
a = 2.0
f = lambda x: 1.0 / (a - np.cos(x))
N = 32
theta = 2.0 * np.pi * np.arange(N) / N
c_fft = np.fft.fft(f(theta)) / N
c_fft_sym = np.fft.fftshift(c_fft)        # reorder to k = -N/2..N/2 - 1
k_sym = np.arange(-N // 2, N // 2)
r = a - np.sqrt(a**2 - 1.0)
c_exact = (1.0 / np.sqrt(a**2 - 1.0)) * r ** np.abs(k_sym)
print("Max error in resolved band:", np.max(np.abs(c_fft_sym[1:-1] - c_exact[1:-1])))
```

The equivalent MATLAB implementation:

```matlab
a = 2; N = 32;
f = @(x) 1 ./ (a - cos(x));
theta = 2*pi*(0:N-1)/N;
c_fft = fft(f(theta)) / N;
c_fft_sym = fftshift(c_fft);
k_sym = -N/2 : N/2 - 1;
r = a - sqrt(a^2 - 1);
c_exact = (1 / sqrt(a^2 - 1)) * r .^ abs(k_sym);
disp(max(abs(c_fft_sym(2:end-1) - c_exact(2:end-1))))
```

The Julia implementation:

```julia
using FFTW
a = 2.0; N = 32
f(x) = 1 / (a - cos(x))
θ = 2π * (0:N-1) / N
c_fft = fft(f.(θ)) / N
c_fft_sym = fftshift(c_fft)
k_sym = -N÷2 : N÷2 - 1
r = a - sqrt(a^2 - 1)
c_exact = [(1 / sqrt(a^2 - 1)) * r^abs(k) for k in k_sym]
println("Max error: ", maximum(abs.(c_fft_sym[2:end-1] .- c_exact[2:end-1])))
```

@fig-fft-coefficients shows the result. In the resolved band $|k| lt.eq.slant 7$ or so, the FFT-computed coefficients agree with the exact formula to roughly machine precision. As $|k|$ approaches $N \/ 2 = 16$, the error grows because of the truncation and aliasing tail of the geometric series, in exact agreement with the master error formula @eq-trap-error applied to each $hat(c)_k$ individually.

#figure(
  image("../figures/ch16/python/fft_coefficients.pdf", width: 95%),
  caption: [(a) Fourier coefficients $|hat(f)_k|$ of $1\/(2 - cos x)$. The exact closed form $r^(|k|) \/ sqrt(3)$ (navy) and the FFT-computed values from $N = 32$ samples (coral dots) coincide on the plot. (b) The absolute error in each FFT coefficient. In the resolved band the error is at machine precision; near the Nyquist frequency $|k| = N\/2$ the error grows because of aliasing tails, but it remains far below the magnitude of the coefficient itself. _Every spectral coefficient computed by the FFT in this textbook is, secretly, an instance of the periodic trapezoidal rule._],
) <fig-fft-coefficients>

The code generating @fig-fft-coefficients is available in:
- `codes/python/ch16/trap_fft_coefficients.py`
- `codes/matlab/ch16/trap_fft_coefficients.m`
- `codes/julia/ch16/trap_fft_coefficients.jl`

// ============================================================================
== The Doubled-Rate Observation <sec-doubled-rate>
// ============================================================================

We have one more observation to make, and it is the perfect bookend to the band-limited exactness theorem of @sec-trig-exactness. Consider an analytic periodic integrand $f$. By Theorem 2 the trapezoidal-rule error converges as $cal(O)(e^(-a N))$, where $a$ is the half-width of the strip of analyticity. Now consider the corresponding _interpolation_ problem: the trigonometric interpolant $f_N$ of $f$ at the same $N$ points. By the standard analysis of Fourier interpolation (which we developed in @ch-fourier-grids), $f_N - f$ converges only as $cal(O)(e^(-a N \/ 2))$ in the $L^infinity$ norm: the rate is _half_ as fast as for the trapezoidal-rule integration.

The discrepancy is paid for by aliasing. When we interpolate, the aliased modes $c_(k + N), c_(k - N), dots$ are added to $c_k$ and produce a real reconstruction error at every grid point. When we integrate, the same aliased modes are added to $c_0$, and the master error formula @eq-trap-error tells us that they contribute the trapezoidal error. But _half_ of those modes have negative index and integrate exactly to zero by themselves; they cancel against their positive-index partners in the symmetric expression for the error and contribute only at second order. The net effect is that the integration error is asymptotically half the order of the interpolation error in the exponent, and so converges twice as fast.

Section 8 of @TrefethenWeideman2014 gives a rigorous statement of this folk theorem and includes a beautiful figure (their Figure 8.1) showing the two convergence rates side by side. The pedagogical takeaway for our purposes is the slogan:

#block(
  fill: rgb("#142D6E").lighten(92%),
  stroke: (left: 3pt + rgb("#142D6E")),
  inset: (left: 12pt, y: 10pt, right: 10pt),
  width: 100%,
)[
*The trapezoidal rule converges twice as fast as Fourier interpolation*, asymptotically, for analytic periodic integrands. This is because the aliased high-frequency modes that destroy interpolation at every grid point integrate (half of them) to zero -- the very fact that explained the band-limited exactness theorem of @sec-trig-exactness.
]

// ============================================================================
== Bridge to Fourier Spectral Methods <sec-trap-bridge>
// ============================================================================

We can now state the final synthesis of the chapter. Throughout @ch-fourier-grids, @ch-spectral-pde, @ch-fourier-pseudo and @ch-polar we built Fourier spectral methods on the foundational claim that an analytic periodic function is approximated to spectral accuracy by its trigonometric interpolant on $N$ equispaced points. We have now retrofitted that claim with a precise quadrature theorem. The two viewpoints --- Fourier interpolation in the time domain and the trapezoidal rule for the Fourier coefficient integrals in the frequency domain --- are not analogous, they are _identical_, and they share the same spectral convergence rate (with the trapezoidal rule winning by a factor of two in the exponent).

#block(
  fill: rgb("#142D6E").lighten(92%),
  stroke: (left: 3pt + rgb("#142D6E")),
  inset: (left: 12pt, y: 10pt, right: 10pt),
  width: 100%,
)[
*The Grand Moral, revisited.* The grand moral of @ch-quadrature was: choose the approximation space first, and let the quadrature follow. The grand moral of this chapter is the same, applied to a different space. For periodic integrands, the right approximation space is _trigonometric polynomials_, and the right quadrature is the equispaced trapezoidal rule. The spectacular accuracy of this combination is _not_ a numerical curiosity. It is the algebraic foundation of every Fourier spectral method ever written, and the engine that drives the FFT-based pseudospectral algorithms of @ch-fourier-pseudo. The trapezoidal rule, taught in elementary calculus as a humble first-order method, is in fact a Fourier spectral method in disguise.
]

This principle explains the central paradox we set out from. The textbook bound @eq-trap-textbook is wildly pessimistic for analytic periodic integrands not because the rule is overperforming, but because the bound is built from the wrong norm: it measures distance in $L^infinity$ over a finite interval, whereas the trapezoidal rule lives in the world of Fourier coefficients, where smooth periodic functions are exponentially small at high frequencies. Each of our five convergence classes is just a different rate of decay for the high-frequency tail of the Fourier expansion, and each produces a correspondingly precise quadrature rate.

A practical algorithmic checklist for the reader, distilling the chapter into four steps:

1. Put the periodic data on an equispaced grid.
2. Think in Fourier coefficients (the master error formula @eq-trap-error is everything).
3. Estimate the half-width $a$ of the strip of analyticity by locating the nearest complex singularity.
4. Then let the periodic trapezoidal rule do the work.

This is exactly the mentality of Fourier spectral methods, and it carries over verbatim to the FFT-based discretisations of the previous chapters.

// ============================================================================
== A Non-Exhaustive Literature Overview <sec-trap-literature>
// ============================================================================

The story of the periodic trapezoidal rule is one of the most surprising in numerical analysis. A formula introduced in elementary calculus as a humble first-order approximation turns out to be, on the right class of integrands, one of the most accurate quadrature methods ever devised. This realisation crystallised over almost two centuries, from Poisson's 1820s computation of the perimeter of an ellipse to the contemporary research programmes that exploit the same complex-analytic mechanism in solid-state physics, boundary integral equations, and rational approximation theory. The following overview traces this arc, structured around the master error formula @eq-trap-error and the five-class taxonomy of @sec-taxonomy that we developed in this chapter.

=== Origins: Poisson, Davis, Goodwin, and the Sinc-Methods Tradition

The first known instance of the periodic trapezoidal rule being used as a high-accuracy device dates to Siméon-Denis Poisson in the 1820s. Confronted with the integral
$ E(k) = integral_0^(pi \/ 2) sqrt(1 - k^2 sin^2 theta) dif theta $
for the perimeter of an ellipse, Poisson observed that the integrand is real-analytic and periodic on $[0, 2 pi]$, and that the equispaced trapezoidal rule applied to it converges spectacularly fast. With sixteen function values he obtained ten correct digits, an accuracy that no Newton--Cotes formula of comparable cost could match. Poisson did not have the language of complex analysis to explain his success; the mechanism (analyticity in a strip whose half-width is set by the imaginary part of the nearest branch point) would only be uncovered more than a century later. The historical attribution chain is reconstructed in detail in @TrefethenWeideman2014.

The first formal proof that periodic analytic integrands enjoy geometric trapezoidal convergence was given by Philip Davis in his 1959 contribution to the Madison symposium on numerical approximation @Davis1959. Davis showed that if $f$ is $2 pi$-periodic and analytic in an open strip of half-width $a > 0$, the trapezoidal-rule error decays as $cal(O)(e^(-a N))$, with the constant $a$ identified as the imaginary part of the nearest complex singularity. Davis's argument was the prototype of the contour-deformation reasoning that drives the strip-analyticity theorem of @sec-strip-analyticity.

In parallel, Goodwin's 1949 paper @Goodwin1949 showed that the same mechanism survives the passage from the periodic interval to the real line. Goodwin demonstrated that an equispaced trapezoidal sum applied to the Gaussian-decay integrand $e^(-x^2) f(x)$, with $f$ analytic in a strip, achieves machine precision with only a dozen samples. This is the historical seed of the real-line theorem we proved in @sec-real-line, and it remains the cleanest illustration of why entire integrands enjoy supergeometric convergence on $bb(R)$.

These two strands (Poisson--Davis on the periodic side and Goodwin on the real-line side) were synthesised into a comprehensive numerical analysis discipline by Frank Stenger and his collaborators over the following four decades. Stenger's 1993 monograph @Stenger1993 codifies the resulting body of knowledge under the umbrella of _sinc methods_, in which the equispaced trapezoidal rule for integration appears as one face of a broader theory that also encompasses sinc interpolation, sinc-Galerkin methods, and indefinite integration. The sinc-methods framework provides the unifying language in which the periodic and real-line theorems are recognised as two manifestations of a single complex-analytic mechanism.

=== The Modern Synthesis: Trefethen and Weideman

The contemporary picture of the periodic trapezoidal rule was established by two papers of Lloyd N. Trefethen and J. André C. Weideman that frame the entire chapter you have just read. The first, Weideman's 2002 _American Mathematical Monthly_ article @Weideman2002, articulated the textbook paradox in pedagogically pointed terms and introduced the explicit closed-form error formulas (for $1\/(a - cos x)$, $|sin(x \/ 2)|^k$, $e^(cos x)$) that we reproduced in @sec-taxonomy. Weideman's deliberate framing of the issue (the classical $cal(O)(N^(-2))$ bound for $f(x) = e^(cos x)$ overestimates the actual error by eight orders of magnitude at $N = 10$, "like saying that the distance between New York and London is less than $10^(11)$ miles") has become the canonical statement of the disagreement between piecewise-polynomial and trigonometric exactness.

The decisive synthesis came twelve years later in the Trefethen--Weideman 2014 _SIAM Review_ paper @TrefethenWeideman2014, which is the principal reference for almost every theorem in this chapter. The paper establishes the master error formula, the strip-analyticity theorem, the five-class convergence taxonomy (band-limited, finite smoothness, geometric, supergeometric, subgeometric), the doubled-rate observation comparing quadrature to interpolation, and the explicit translation to the real line via Fourier transforms. Equally importantly, the paper traces the historical pedigree of these ideas back through Davis, Goodwin, Poisson, and the Russian school of approximation theory, situating the trapezoidal rule as the beating heart of Fourier spectral methods rather than as an isolated curiosity. Trefethen's 2016 quadrature talk @Trefethen2016Quadrature subsequently distilled the same message for a broader audience, framing it as one of "ten things you should know about quadrature".

The conceptual shift these papers codified (from polynomial to trigonometric exactness) is the same shift that drives the entire programme of @ch-quadrature, @ch-fourier-grids, @ch-spectral-pde and @ch-fourier-pseudo. The trapezoidal rule becomes spectacular not because the formula is clever but because the function space in which we measure its error has been correctly identified.

=== The Real Line and Double-Exponential Transformations

Goodwin's real-line bound of @sec-real-line invites a natural question: what if the integrand is _not_ analytic in a strip around the real axis, but only in a strip around some other curve? The answer is the celebrated _double-exponential_ (DE) transformation introduced by Hidetosi Takahasi and Masatake Mori at the Kyoto Research Institute for Mathematical Sciences in 1974 @TakahasiMori1974. Their idea is to apply a nonlinear coordinate change such as $x = tanh(frac(pi, 2) sinh t)$ to a finite-interval integral with endpoint singularities. Under this map, the original integrand is pushed onto $bb(R)$ and acquires double-exponential decay $exp(-exp(|t|))$ at infinity, while the singularities at the endpoints of the original interval are pushed to $plus.minus infinity$. The transformed integrand satisfies the analyticity-and-decay hypotheses of the real-line trapezoidal theorem of @sec-real-line, and a few dozen equispaced samples produce machine precision regardless of the strength of the original endpoint singularities.

The DE family has been extended in the half-century since 1974 to oscillatory integrands, indefinite integrals, multiple integrals, Sinc--collocation discretisations of differential equations, and Fourier integrals. Mori and Sugihara's 2001 review @MoriSugihara2001 documents the full breadth of the programme, and modern implementations remain among the most reliable black-box quadrature routines available for difficult integrals on finite intervals. From the perspective of this chapter, the DE transformation is the cleanest demonstration of the principle that the trapezoidal rule's spectral accuracy is _portable_: any nonlinear change of variables that produces an analytic, rapidly decaying integrand on $bb(R)$ inherits the geometric convergence of @sec-real-line for free. This portability is also at the heart of Stenger's broader sinc-methods framework @Stenger1993.

=== Singular Boundary Integral Equations and the Zeta Correction

A second contemporary domain in which the periodic trapezoidal rule is the workhorse, but in which the strip-analyticity hypothesis fails dramatically, is the discretisation of boundary integral equations (BIEs) for the Laplace, Helmholtz, and Stokes equations. On a smooth closed contour the layer-potential kernel is naturally periodic, which suggests applying the periodic trapezoidal rule to the parametrised integral. The catch is that the kernel is logarithmically (in 2D) or algebraically (in 3D) singular when the integration variable approaches the target point. Under the five-class taxonomy of @sec-taxonomy the integrand falls into the finite-smoothness regime, and the trapezoidal rule converges only at a slow algebraic rate dictated by the highest odd derivative that remains continuous across the singular point.

Two divergent strategies were developed to recover spectral accuracy. Rainer Kress's 1990 Nyström method @Kress1990 splits the kernel into a smooth part and an explicitly singular part, integrates the singular part analytically using Fourier expansions, and applies the periodic trapezoidal rule to the smooth remainder. The Kress approach restores high-order convergence but produces dense interaction matrices that are difficult to combine with the Fast Multipole Method. Sharad Kapur and Vladimir Rokhlin's 1997 alternative @KapurRokhlin1997 takes a complementary route: instead of correcting the kernel globally, they leave the standard trapezoidal rule untouched everywhere except in a small stencil immediately adjacent to the singular target node, where the weights are replaced by precomputed corrections. This locally corrected scheme preserves the sparse structure required by fast solvers but only achieves modest order of convergence.

The breakthrough came in a series of papers by Bowei Wu and Per-Gunnar Martinsson published between 2020 and 2023, in which the two paradigms were unified through what the authors call the _zeta correction_ method @WuMartinsson2021 @WuMartinsson2023. The key analytical insight is that the asymptotic error coefficients of the punctured trapezoidal rule (in which the singular target node is simply omitted from the sum) can be evaluated in closed form using specific values of the Riemann zeta function in 2D and the Epstein zeta function in 3D. By matching these explicit error coefficients with a small fixed-size local stencil (typically 9 to 16 points), one obtains a correction whose support is local but whose order can be driven arbitrarily high. Numerical implementations have demonstrated up to 32nd-order convergence for 2D boundary integral equations and 9th-order for 3D surfaces, all while preserving the $cal(O)(N)$ cost structure required for Fast Multipole acceleration. The zeta correction is one of the cleanest contemporary illustrations of the principle that "trapezoidal rule + tiny analytic correction = spectral accuracy", and it has reframed the place of the periodic trapezoidal rule in computational electromagnetics, fluid dynamics, and acoustics. The same principle --- that a small, analytically informed correction can elevate a low-order quadrature rule to high or spectral accuracy --- appears in other guises as well: in complex-plane contour integration, where the trapezoidal rule on a deformed contour inherits exponential convergence from the analyticity of the integrand, and in the evaluation of Caputo integrals for fractional derivatives, where end-corrected trapezoidal rules can be highly effective @Fornberg2025[Sections 7.5 and 8.3--8.4].

=== Lightning Plus Polynomial: Spectral Accuracy on Non-Periodic Corner Domains

A third research thread, superficially distant from periodic quadrature but driven by the same complex-analytic philosophy, is the rational-function approach to functions with corner singularities. When standard polynomial expansions (Chebyshev, Legendre) are applied to a function such as $f(x) = x^alpha$ on a sector domain, the convergence rate collapses to a slow algebraic rate dictated by the strength of the singularity. Abinand Gopal and Lloyd Trefethen's 2019 _lightning solver_ @GopalTrefethen2019 bypasses polynomials entirely: it approximates the solution of a Laplace problem with corner singularities as the real part of a rational function whose poles are _preassigned_ in an exponentially clustered sequence tapering toward the singular corner. By fixing the pole locations a priori the difficult nonlinear optimisation of pole placement is replaced by a well-conditioned linear least-squares problem evaluated on the boundary, and the resulting approximation achieves root-exponential convergence on the original domain.

The lightning approach has one weakness: the exponentially clustered poles, so effective at absorbing the local singularity, are inefficient at modelling the smooth global background of the function across the interior of the domain. Forcing the clustered poles to do double duty results in catastrophic ill-conditioning of the least-squares system. Astrid Herremans, Daan Huybrechs, and Trefethen resolved this in their 2023 paper @HerremansHuybrechsTrefethen2023 by introducing the _lightning plus polynomial_ architecture: the rational basis is augmented with a low-degree polynomial sequence that absorbs the smooth global behaviour, while the local lightning poles are left exclusively to manage the microscopic defect at the corner. This division of labour stabilises the linear system and elevates the lightning method from a practical algorithm into a theoretically optimal spectral framework. Sharp convergence theorems were subsequently proved by Xiang, Yang, and Wu @XiangYangWu2024, who showed that the optimal taper of the clustered poles produces a root-exponential error envelope $cal(O)(exp(-C sqrt(N)))$, attaining the classical Vyacheslavov--Stahl bounds for rational approximation of functions with branch points.

The reason this material belongs in a chapter on periodic quadrature is that the underlying mechanism is identical to the one we have studied for the trapezoidal rule. Whether one is integrating a smooth periodic function on $[0, 2 pi]$, applying a zeta correction to a singular layer potential, or approximating a non-periodic function with corner branch points by a rational function, the asymptotic error is universally controlled by the analytic continuation of the integrand into the complex plane and the distance to the nearest singularity. The trapezoidal rule and the lightning method are dialects of the same complex-analytic language; what differs is only the choice of approximation space.

=== The Periodic Trapezoidal Rule in Solid-State Physics: Brillouin Zone Integration

A fourth domain in which the periodic trapezoidal rule reigns supreme, this time under a different name, is the computation of densities of states and optical conductivities in solid-state physics. In a periodic crystal lattice the macroscopic observables are obtained by integrating energy bands over the _Brillouin zone_, the unit cell of the reciprocal momentum lattice. Because the Brillouin zone is intrinsically periodic, the natural quadrature is the equispaced trapezoidal rule on a regular grid, and the rule has been the de facto standard for density-functional theory calculations since Hendrik Monkhorst and James Pack introduced it in 1976 under the name _Monkhorst--Pack grid_ @MonkhorstPack1976. From the perspective of @sec-taxonomy a Monkhorst--Pack grid is precisely the periodic trapezoidal rule in three dimensions: nothing more, nothing less.

The pathology that physicists confront is the singular nature of the density of states $delta(epsilon_(n bold(k)) - E)$, which is regularised in numerical practice by replacing the Dirac distribution with a smoothed kernel (Gaussian, Lorentzian, or Fermi--Dirac derivative) governed by a small _smearing parameter_ $eta > 0$. While the smearing technically restores analyticity, it does so at the cost of pulling the complex singularities of the integrand close to the real axis. In the language of @sec-strip-analyticity, the half-width $a$ of the analyticity strip shrinks proportionally to $eta$, and the geometric convergence of the trapezoidal rule degrades catastrophically: a fixed error tolerance requires $cal(O)(eta^(-1))$ nodes per spatial dimension, leading to $cal(O)(eta^(-3))$ scaling in three-dimensional materials. For low-temperature properties of strongly correlated systems studied via Dynamical Mean-Field Theory, the required values of $eta$ are infinitesimally small and the standard Monkhorst--Pack approach becomes computationally intractable.

A recent breakthrough is the iterated adaptive integration (IAI) algorithm of Jason Kaye, Sophie Beck, Alex Barnett, Lorenzo Van Muñoz, and Olivier Parcollet @KayeBeckBarnett2023. IAI exploits the observation that a heavily smeared density of states is not uniformly difficult: its sharp features reside on a low-dimensional Fermi surface, and uniform sampling wastes computational effort throughout the bulk of the Brillouin zone. The IAI algorithm uses recursive, dimension-by-dimension adaptive Gaussian quadrature to focus refinement exclusively on the narrow transport distributions near the Fermi surface, achieving polylogarithmic cost $cal(O)(log^3(eta^(-1)))$ as $eta arrow 0^+$, an exponential speedup over the equispaced Monkhorst--Pack grid. Ongoing follow-up work explores the alternative strategy of analytically continuing the integrand into the complex plane and deforming the integration contour away from the physical singularities, restoring the wide analyticity strip and recovering the full geometric convergence rate of the periodic trapezoidal rule. In both approaches the conceptual recipe is the same as in the rest of this chapter: locate the nearest singularity, choose the right approximation space, and let the trapezoidal rule do the work.

=== Closing Perspective: Accelerating Subgeometric Convergence

The fifth and final research thread in our overview concerns the rescue of subgeometric methods, the pathological tier at the bottom of the convergence taxonomy of @sec-taxonomy. Historically, the slowdown observed when applying Hermite or Laguerre quadrature to functions on unbounded domains was treated as an unavoidable penalty: the spectral methods community broadly accepted that Gauss--Hermite was inferior to truncated trapezoidal rules for large classes of integrands, as we saw in @sec-beyond-gauss of @ch-quadrature. A 2024 paper by Hao Hu and Haijun Yu @HuYu2024 has overturned this view by providing a rigorous error analysis of _scaled_ Hermite approximations in which the spatial scaling factor is dynamically chosen to balance the spatial truncation error against the frequency truncation error. They show that for functions originally exhibiting subgeometric Hermite convergence, optimal scaling accelerates the rate to full geometric exponential convergence, and for functions of merely finite smoothness it doubles the order of algebraic decay. The mechanism is structurally identical to the contour-tuning technique that Weideman and Trefethen had developed two decades earlier for inverse Laplace transforms via Talbot's method @WeidemanTrefethen2007: by dynamically optimising the geometric parameters of the integration contour as a function of the node count $N$, the trapezoidal-based Bromwich integral can be upgraded from a stubbornly subgeometric $cal(O)(exp(-c sqrt(N)))$ to a true geometric $cal(O)(exp(-c N))$. A parallel thread of work on supergeometric convergence of spectral collocation methods for weakly singular Volterra and Fredholm equations @HuangTangZhang2011 confirms that this accelerated-convergence phenomenon is universal across the broader landscape of spectral methods.

Across all five threads (the historical foundations, the Trefethen--Weideman synthesis, the double-exponential transformations of Takahasi--Mori, the zeta-corrected boundary integral equations of Wu--Martinsson, the lightning plus polynomial methods of Gopal--Trefethen and Herremans--Huybrechs--Trefethen, the Brillouin-zone integration algorithms of Kaye and collaborators, and the optimally scaled Hermite quadratures of Hu--Yu) the underlying mechanics remain consistent. Managing the distance to the nearest complex singularity, exploiting aliasing cancellations on equispaced grids, and dynamically tuning integration contours are the universal strategies required to unlock spectral convergence. The continuing success of these methods across machine learning, quantum mechanics, fluid dynamics, and computational electromagnetics underscores the foundational principle that closes both this chapter and the previous one: selecting the correct exactness space fundamentally overrides the algebraic degree of the underlying computational formula.

// ============================================================================
== Summary <sec-trap-summary>
// ============================================================================

This chapter has developed the theory of the periodic trapezoidal rule and its real-line cousin, guided by two papers of Trefethen and Weideman @TrefethenWeideman2014 @Weideman2002. The key findings are:

- *Trigonometric exactness.* The $N$-point periodic trapezoidal rule is exact for every trigonometric polynomial of degree $lt.eq.slant N - 1$. Even for single modes $cos(k theta)$ with $N\/2 < k < N$ -- aliased modes that cannot be reconstructed -- the rule still gives the correct answer of zero.
- *The master error formula.* For any Lipschitz periodic $f$, $T_N (f) - I(f) = 2 pi sum_(ell eq.not 0) c_(ell N)$. The trapezoidal-rule error is therefore controlled entirely by the high-frequency tail of the Fourier expansion of $f$.
- *Five-class taxonomy.* Periodic integrands fall into five convergence classes: band-limited (exact), finite smoothness (algebraic, with rate set by the regularity of the periodic extension), strip-analytic (geometric, $cal(O)(e^(-a N))$ with $a$ the strip half-width), entire (supergeometric, $(e \/ (2 N))^N$), and $C^infinity$-but-not-analytic (subgeometric, e.g.\u{a0}$exp(-c N^(2 \/ 3))$).
- *The strip-analyticity theorem.* If $f$ is analytic and bounded by $M$ in the open strip $|"Im" theta| < a$, then $|I_N - I| lt.eq.slant 4 pi M \/ (e^(a N) - 1)$. The constant $a$ is the only quantity the practitioner needs to estimate; it is the imaginary part of the nearest complex singularity of $f$.
- *Real-line generalisation.* For analytic, suitably decaying integrands on $bb(R)$, the equispaced trapezoidal rule with step $h = 2 pi \/ N$ gives an error of $cal(O)(e^(-2 pi a \/ h))$. This is the bridge from periodic to non-periodic spectral integration and the foundation of inverse-Laplace and Hankel-contour techniques.
- *FFT computation of Fourier coefficients.* Every Fourier coefficient computed by the FFT is, secretly, the periodic trapezoidal rule applied to the defining integral. Every Fourier spectral method in this textbook therefore inherits the convergence theory developed in this chapter.
- *The doubled-rate observation.* For analytic periodic integrands, the trapezoidal rule converges asymptotically twice as fast as the underlying trigonometric interpolant, because half of the aliased modes integrate exactly to zero.

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    table(
      columns: (0.6fr, 1.3fr, 1.1fr, 1.2fr),
      align: (center, left, left, left),
      inset: (x: 0.8em, y: 0.5em),
      stroke: none,
      table.hline(stroke: 0.75pt + rgb("#142D6E")),
      table.header(
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Étude*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Function*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Class*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Highlight*],
      ),
      table.hline(stroke: 0.5pt + luma(180)),
      [16.1], [Poisson's ellipse $sqrt(1 - 0.36 sin^2 theta)$], [geometric], [10 digits in 8 samples],
      [16.2], [random trig.\u{a0}polynomial of degree $m$], [exact for $N > m$], [aliased modes still integrate to zero],
      [16.3], [$f_2 = |sin(x\/2)|$, $f_3 = |sin(x\/2)|^3$], [algebraic], [slopes $-2$ and $-4$ on log-log],
      [16.4], [Poisson kernel $1\/(2 - cos x)$], [geometric], [exact closed-form error formula],
      [16.5], [$f_5 = e^(cos x)$], [supergeometric], [each new node adds more than one digit],
      [16.6], [Weideman's $f_6 = exp((cos x - 1)\/(cos x + 1))$], [subgeometric], [linear on $N^(2\/3)$ axis],
      [16.7], [Real-line Gaussian $e^(-x^2) \/ sqrt(pi)$], [geometric on $bb(R)$], [machine precision at $h = pi\/6$],
      [16.8], [FFT coefficients of $1\/(2 - cos x)$], [exact / geometric], [FFT $equiv$ trapezoidal rule for $hat(f)_k$],
    ),
  ),
  caption: [Summary of computational études in this chapter.],
) <tbl-ch16-summary>

The overarching lesson, parallel to that of @ch-quadrature: *the right exactness space matters more than the polynomial degree of the formula*. For periodic integrands the right exactness space is trigonometric, and the equispaced trapezoidal rule -- humble in algebraic terms -- is the natural quadrature in that space. Combined with the FFT, it is the engine of every Fourier spectral method we have developed in this textbook.

// ============================================================================
== Exercises <sec-trap-exercises>
// ============================================================================

*Exercise 16.1* (_Half-circles and ellipses_). Compute the perimeter of an ellipse with semi-axes $1$ and $0.8$ using the periodic trapezoidal rule. (a) Locate the branch points of the integrand $sqrt(1 - 0.36 sin^2 theta)$ in the complex plane. (b) Predict the geometric rate from your answer in (a). (c) Verify the rate numerically by plotting $log|I_N - I|$ against $N$ for $N = 4, 8, dots, 60$.

*Exercise 16.2* (_The Poisson kernel for varying $a$_). Repeat Étude 16.4 for $a = 1.1$, $a = 1.5$, $a = 2$, and $a = 5$. (a) For each $a$, predict the rate $r = a - sqrt(a^2 - 1)$ and the strip half-width $log(a + sqrt(a^2 - 1))$. (b) Verify the rates by linear regression of $log_(10) |I_N - I|$ on $N$. (c) Discuss the limiting behaviour as $a arrow 1^+$: what happens to the strip and to the convergence rate?

*Exercise 16.3* (_Match the rate to the function_). For each of the following functions on $[0, 2 pi]$, predict the convergence class and rate of the trapezoidal rule, then verify numerically: (a) $f(x) = sin(7 x) + cos(13 x)$; (b) $f(x) = abs(x - pi)$; (c) $f(x) = log(2 + cos x)$; (d) $f(x) = e^(2 cos x)$; (e) $f(x) = (cos x)^(1\/3)$ (with the principal branch).

*Exercise 16.4* (_The doubled-rate observation_). Take $f(x) = 1\/(2 - cos x)$. (a) Compute the Fourier interpolant $f_N$ at $N = 16$ Chebyshev-extreme points (no, equispaced points: this is a periodic function, use equispaced) and measure $||f_N - f||_infinity$. (b) Compute the trapezoidal-rule error $|T_N (f) - I|$. (c) Show that the interpolation error decays as $r^(N\/2)$ while the quadrature error decays as $r^N$, with $r = 2 - sqrt(3)$, by repeating for $N = 4, 8, 16, 32, 64$.

*Exercise 16.5* (_Goodwin's example beyond the Gaussian_). The Goodwin paper @Goodwin1949 also discusses $w(x) = e^(-x^2) f(x)$ for various smooth $f$. Take $f(x) = cos(x^3)$ and compute $integral_(-infinity)^infinity e^(-x^2) cos(x^3) dif x$ by the real-line trapezoidal rule. (a) Find the smallest $N$ that gives 10 correct digits. (b) Plot the convergence curve. (c) Compare with the result of Étude 15.10 (Gauss--Hermite vs.\u{a0}truncated quadrature) on the same integrand: which method is faster?

*Exercise 16.6* (_FFT vs.\u{a0}direct sum_). Verify in Python or Julia that the FFT-computed Fourier coefficients of $f(x) = 1 \/ (2 - cos x)$ at $N = 64$ agree with the closed form $r^(|k|) \/ sqrt(3)$ to machine precision in the resolved band $|k| lt.eq.slant 30$. Time the computation against a direct $cal(O)(N^2)$ sum and confirm the predicted $cal(O)(N log N)$ scaling.

*Exercise 16.7* (_Symmetry exploitation_). For real, even integrands (like Poisson's ellipse), the trapezoidal sum can be reduced to a sum over $N \/ 4 + 1$ values by combining the four-fold symmetry $f(theta) = f(-theta) = f(pi - theta) = f(pi + theta)$. (a) Derive the symmetry-reduced formula. (b) Verify that it gives the same numerical results as the unreduced sum on Poisson's ellipse. (c) Estimate the practical speedup for $N = 1024$.

*Exercise 16.8* (_The endpoint of subgeometric_). Construct a $C^infinity$-periodic function whose periodic extension has an essential singularity at one point of the period, and whose trapezoidal-rule error decays as $exp(-c N^(1\/2))$ instead of $exp(-c N^(2\/3))$. (Hint: replace $exp((cos x - 1)\/(cos x + 1))$ with a similar construction whose Fourier coefficients have a different decay structure.)

*Exercise 16.9* (_Reflection essay_). Write a two-page essay on the question: "Why is it that the trapezoidal rule, taught in elementary calculus as a crude first-order method, becomes a spectral method on periodic integrands?" Frame your answer in terms of approximation spaces, aliasing, and the geometry of complex singularities. The essay should _not_ rely on the Euler--Maclaurin formula at any point.
