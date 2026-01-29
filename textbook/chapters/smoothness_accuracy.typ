// textbook/chapters/smoothness_accuracy.typ
// Chapter 6: Smoothness and Spectral Accuracy
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Last modified: January 2026

#import "../styles/template.typ": dropcap, ii

// Enable equation numbering for this chapter
#set math.equation(numbering: "(1)")

= Smoothness and Spectral Accuracy <ch-smoothness>

#dropcap[In the previous chapter, we witnessed spectral methods achieving machine precision with remarkably few grid points. The test function $u(x) = 1\/(2 + sin(x))$ was differentiated to fourteen-digit accuracy using only fifty nodes, while finite difference methods of any fixed order would require millions of points to achieve comparable precision. But _why_? What is the source of this extraordinary accuracy?]

The answer lies in a beautiful chain of reasoning that connects the _smoothness_ of a function to the _accuracy_ of its numerical differentiation. This chapter develops the theoretical foundations that explain the spectacular performance of spectral methods. The argument proceeds in two steps:

1. *Smoothness implies rapid decay*: A smooth function has little energy at high wavenumbers, so its Fourier coefficients decay rapidly.

2. *Rapid decay implies small aliasing error*: When we sample a function on a discrete grid, high frequencies "fold" onto low frequencies through a phenomenon called aliasing. If the high-frequency coefficients are negligible, this aliasing causes negligible error.

These two steps, made precise by the theorems in this chapter, constitute the fundamental explanation for spectral accuracy. Understanding them provides insight not only into _why_ spectral methods work, but also into _when_ they work: the rate of convergence is controlled by the smoothness of the function being approximated.

== The Two Steps of Accuracy <sec-two-steps>

=== The Heuristic Picture

Consider a smooth, periodic function $f(x)$ on $[0, 2 pi)$. Its Fourier series representation is
$ f(x) = sum_(k = -infinity)^infinity hat(f)_k e^(ii k x), $ <eq-fourier-series>
where the Fourier coefficients are
$ hat(f)_k = frac(1, 2 pi) integral_0^(2 pi) f(x) e^(-ii k x) dif x. $ <eq-fourier-coeff>

The wavenumber $k$ corresponds to oscillations with frequency $|k|$ per period. A smooth function, by definition, changes gradually; it has no abrupt jumps or corners. Such a function cannot have significant energy at high frequencies, for high-frequency waves oscillate rapidly. Therefore, the coefficients $hat(f)_k$ must decay as $|k| arrow infinity$.

Now suppose we sample $f$ at $N$ equispaced points, producing the grid function $v_j = f(x_j)$ where $x_j = 2 pi j \/ N$. The discrete Fourier transform of this grid function involves only $N$ coefficients, corresponding to wavenumbers $|k| lt.eq.slant N\/2$. What happened to the higher wavenumbers?

The answer is _aliasing_. The discrete grid cannot distinguish between a wave with wavenumber $k$ and a wave with wavenumber $k + N$, since $e^(ii k x_j) = e^(ii (k + N) x_j)$ for all grid points $x_j$. High frequencies masquerade as low frequencies; they "fold" into the resolved spectral range. The discrete Fourier coefficient $tilde(f)_k$ is not equal to the continuous coefficient $hat(f)_k$, but rather to the sum of all coefficients that alias to $k$:
$ tilde(f)_k = sum_(j = -infinity)^infinity hat(f)_(k + j N). $ <eq-aliasing-preview>

This is the _aliasing formula_, and it reveals the key insight: if $hat(f)_k$ decays rapidly, the aliased contributions $hat(f)_(k plus.minus N), hat(f)_(k plus.minus 2 N), dots.h.c$ are negligible. The discrete coefficients $tilde(f)_k$ are then excellent approximations to the continuous coefficients $hat(f)_k$, and spectral methods achieve their remarkable accuracy.

=== Roadmap of This Chapter

We shall make this heuristic argument precise through four theorems:

- *Theorem 1* classifies functions by their smoothness and establishes the corresponding decay rates of their Fourier coefficients.
- *Theorem 2* (the aliasing formula) quantifies how high frequencies contaminate low frequencies upon discretization.
- *Theorems 3 and 4* combine these results to bound the errors in spectral interpolation and differentiation.

The following sections make this heuristic argument precise through four theorems that quantify the connection between smoothness and spectral accuracy.

=== The Density of Smooth Functions

Before diving into the theorems, we address a natural concern: the analysis that follows focuses on smooth functions, but what if our function of interest is merely continuous, or has limited regularity?

The answer lies in a classical approximation theorem. Smooth functions are _dense_ in the space of continuous functions: for any continuous function $f$ on a compact interval and any $epsilon > 0$, there exists a smooth function $g$ such that $||f - g||_infinity < epsilon$. For periodic functions on $[0, 2 pi]$, trigonometric polynomials provide such approximations; this is the content of the Weierstrass approximation theorem (1885).

The construction is elegant: convolve $f$ with a smooth "bump function" (a mollifier) of small support. The resulting function inherits the smoothness of the mollifier while approximating $f$ arbitrarily well. This technique, known as _mollification_, is fundamental in analysis and partial differential equations.

The implication for spectral methods is profound. Suppose we wish to approximate a merely continuous function $f$. We can:
1. Approximate $f$ by a smooth function $g$ with $||f - g||_infinity < epsilon\/2$;
2. Apply spectral methods to $g$, which (by Theorem 1) has rapidly decaying Fourier coefficients;
3. Achieve total error less than $epsilon$ with sufficiently many grid points.

This justifies our focus on smooth functions: spectral methods provide a _universal approximation framework_ for continuous functions, approaching any target through smooth intermediaries. The smoothness assumption is not a limitation but rather a pathway to understanding.

== The Decay of Fourier Coefficients <sec-decay>

=== Smoothness Classes

The following theorem, due in various parts to mathematicians from Riemann to Paley and Wiener, establishes the fundamental connection between function smoothness and spectral decay. We state it for functions on the real line $RR$; similar results hold for periodic functions.

#block(
  fill: rgb("#142D6E").lighten(92%),
  stroke: (left: 3pt + rgb("#142D6E")),
  inset: (left: 12pt, y: 10pt, right: 10pt),
  width: 100%,
)[
*Theorem 1* (Smoothness and Fourier decay). Let $u in L^2(RR)$ have Fourier transform $hat(u)$.

*(a)* If $u$ has $p - 1$ continuous derivatives in $L^2(RR)$ for some $p gt.eq.slant 0$ and a $p$-th derivative of bounded variation, then
$ hat(u)_k = O(|k|^(-p-1)) quad "as" |k| arrow infinity. $

*(b)* If $u$ has infinitely many continuous derivatives in $L^2(RR)$, then
$ hat(u)_k = O(|k|^(-m)) quad "as" |k| arrow infinity $
for every $m gt.eq.slant 0$.

*(c)* If there exists $a > 0$ such that $u$ can be extended to an analytic function in the complex strip $|"Im"(z)| < a$ with bounded $L^2$ norm along horizontal lines, then
$ hat(u)_k = O(e^(-a |k|)) quad "as" |k| arrow infinity. $

*(d)* If $u$ is entire (analytic throughout $CC$) and satisfies $|u(z)| = o(e^(a |z|))$ as $|z| arrow infinity$ for some $a > 0$, then $hat(u)$ has compact support: $hat(u)_k = 0$ for all $|k| > a$.
]

Parts (c) and (d) are known as the Paley--Wiener theorems. The theorem establishes a hierarchy of smoothness classes, each with its characteristic decay rate:

#figure(
  table(
    columns: (1fr, 1fr, 1fr),
    align: (left, center, center),
    inset: (x: 12pt, y: 8pt),
    stroke: none,
    table.hline(stroke: 1.5pt),
    table.header(
      [*Smoothness Class*],
      [*Decay Rate*],
      [*Example*],
    ),
    table.hline(stroke: 0.75pt),
    [Bounded variation], [$O(|k|^(-1))$], [Step function],
    [$p$ derivatives, $p$-th in BV], [$O(|k|^(-p-1))$], [B-splines],
    [$C^infinity$ (infinitely smooth)], [$O(|k|^(-m))$ for all $m$], [Bump functions],
    [Analytic in strip], [$O(e^(-a |k|))$], [$1\/(1 + x^2)$],
    [Entire], [Faster than any exponential], [$e^(-x^2)$],
    [Band-limited], [Compact support], [$"sinc"(x)$],
    table.hline(stroke: 1.5pt),
  ),
  caption: [The smoothness hierarchy and corresponding Fourier decay rates. More smoothness implies faster decay.],
) <tbl-smoothness-hierarchy>

=== Intuition via Integration by Parts

The algebraic decay rates in parts (a) and (b) can be understood through integration by parts. For a smooth periodic function $f$ with period $2 pi$, we have
$ hat(f)_k = frac(1, 2 pi) integral_0^(2 pi) f(x) e^(-ii k x) dif x = frac(1, ii k) dot frac(1, 2 pi) integral_0^(2 pi) f'(x) e^(-ii k x) dif x = frac(1, ii k) hat(f')_k. $
Each integration by parts gains a factor of $k^(-1)$ in the decay rate. If $f$ has $p$ derivatives, we can integrate by parts $p$ times:
$ hat(f)_k = frac(1, (ii k)^p) hat(f^((p)))_k. $ <eq-integration-by-parts>
If $f^((p))$ has bounded variation (so its Fourier coefficients are $O(k^(-1))$), then $hat(f)_k = O(k^(-p-1))$.

=== Intuition via Complex Analysis

The exponential decay in part (c) arises from complex analysis. If $f(z)$ is analytic in a strip $|"Im"(z)| < a$, we can shift the contour of integration in the Fourier integral:
$ hat(f)_k = integral_(-infinity)^infinity f(x) e^(-ii k x) dif x = integral_(-infinity)^infinity f(x + ii sigma) e^(-ii k (x + ii sigma)) dif x = e^(k sigma) integral_(-infinity)^infinity f(x + ii sigma) e^(-ii k x) dif x $
for any $|sigma| < a$. Taking $sigma = a - epsilon$ for small $epsilon > 0$ (and $k > 0$) gives $|hat(f)_k| lt.eq.slant C e^(-a k)$. The distance from the real axis to the nearest singularity controls the decay rate.

=== Computational Demonstration

@fig-decay-hierarchy illustrates the decay hierarchy for three representative functions on $[0, 2 pi]$:

1. $f_1(x) = |sin(x)|^3$: This function is $C^2$ but not $C^3$ (the third derivative has jump discontinuities). By Theorem 1(a), its Fourier coefficients decay as $O(k^(-4))$.

2. $f_2(x) = 1\/(1 + sin^2(x\/2))$: This function is smooth on the real line but has poles in the complex plane. It is analytic in a strip of finite width, so by Theorem 1(c), its coefficients decay geometrically: $O(c^(-k))$ for some $c > 1$.

3. $f_3(x) = exp(sin(x))$: This function is entire (analytic throughout $CC$). Its coefficients decay faster than any exponential.

#figure(
  image("../figures/ch06/python/decay_hierarchy.pdf", width: 95%),
  caption: [The decay hierarchy of Fourier coefficients. Functions with greater smoothness exhibit faster decay: algebraic ($O(k^(-4))$) for finite regularity, geometric ($O(e^(-alpha k))$) for strip-analytic functions, and super-geometric for entire functions. The plot shows $|hat(f)_k|$ versus wavenumber $k$ on a semi-log scale.],
) <fig-decay-hierarchy>

The following Python code computes the Fourier coefficients:

```python
def fourier_coefficients(f, N):
    """Compute Fourier coefficients via FFT."""
    x = np.linspace(0, 2*np.pi, N, endpoint=False)
    f_hat = np.fft.fft(f(x)) / N
    return np.abs(f_hat[:N//2])
```

The equivalent MATLAB implementation:

```matlab
function f_hat = fourier_coefficients(f, N)
    x = linspace(0, 2*pi, N+1);
    x = x(1:end-1);  % Remove endpoint for periodicity
    f_hat = abs(fft(f(x))) / N;
    f_hat = f_hat(1:N/2+1);
end
```

The code generating @fig-decay-hierarchy is available in:
- `codes/python/ch06_accuracy/fourier_decay.py`
- `codes/matlab/ch06_accuracy/fourier_decay.m`

== The Aliasing Phenomenon <sec-aliasing>

=== The Poisson Summation Formula

When we sample a continuous function at $N$ equispaced points, we lose information about frequencies beyond the _Nyquist frequency_ $N\/2$. This lost information does not simply disappear; it contaminates the lower frequencies through aliasing.

#block(
  fill: rgb("#142D6E").lighten(92%),
  stroke: (left: 3pt + rgb("#142D6E")),
  inset: (left: 12pt, y: 10pt, right: 10pt),
  width: 100%,
)[
*Theorem 2* (Aliasing formula). Let $u in L^2(RR)$ have a first derivative of bounded variation, and let $v$ be the grid function on $h ZZ$ defined by $v_j = u(x_j)$. Then for all $k in [-pi\/h, pi\/h]$,
$ hat(v)_k = sum_(j = -infinity)^infinity hat(u)_(k + 2 pi j \/ h). $ <eq-aliasing-formula>
]

For a periodic function sampled at $N$ equispaced points with spacing $h = 2 pi \/ N$, this becomes:
$ tilde(f)_k = sum_(j = -infinity)^infinity hat(f)_(k + j N). $ <eq-aliasing-periodic>

=== Spectral Folding

The aliasing formula has a vivid geometric interpretation: the frequency axis "folds" onto itself at the Nyquist frequency. All frequencies that differ by multiples of $N$ become indistinguishable on the discrete grid.

Consider a wave $e^(ii k x)$ sampled at $x_j = 2 pi j \/ N$. At the grid points,
$ e^(ii k x_j) = e^(2 pi ii k j \/ N) = e^(2 pi ii (k + N) j \/ N) = e^(ii (k + N) x_j). $
The discrete samples cannot distinguish wavenumber $k$ from wavenumber $k + N$.

@fig-aliasing-visualization demonstrates this phenomenon in three panels:
- *Left*: A function with high-frequency content sampled coarsely. The interpolant through the samples misses the true oscillations.
- *Center*: The frequency folding diagram, showing how wavenumbers outside $[-N\/2, N\/2]$ map into this range.
- *Right*: Comparison of the true spectrum (many points) with the aliased spectrum (few points).

#figure(
  image("../figures/ch06/python/aliasing_visualization.pdf", width: 95%),
  caption: [The aliasing phenomenon. _Left_: A function $f(x) = sin(x) + 0.5 sin(10 x)$ sampled at $N = 12$ points. The high-frequency component $sin(10 x)$ aliases to $sin(-2 x)$, producing a misleading interpolant. _Center_: Frequency folding diagram showing how wavenumbers $k = 10$ and $k = -10$ fold into the resolved range $[-6, 6]$. _Right_: True versus aliased spectrum for $f(x) = 1\/(1.5 + cos(x))$.],
) <fig-aliasing-visualization>

=== Why Smooth Functions Are Immune

For smooth functions, aliasing is harmless. If $hat(f)_k$ decays rapidly, the aliased contributions $hat(f)_(k plus.minus N), hat(f)_(k plus.minus 2 N), dots.h.c$ are negligible compared to $hat(f)_k$. The discrete coefficients $tilde(f)_k$ are then excellent approximations to the continuous coefficients $hat(f)_k$.

Quantitatively, the aliasing error is bounded by the tail of the Fourier series:
$ |tilde(f)_k - hat(f)_k| lt.eq.slant sum_(j eq.not 0) |hat(f)_(k + j N)| lt.eq.slant 2 sum_(m = N\/2)^infinity |hat(f)_m|. $ <eq-aliasing-error>

If $hat(f)_k = O(k^(-p-1))$, this sum is $O(N^(-p))$. If $hat(f)_k = O(e^(-a k))$, the sum is $O(e^(-a N\/2))$, which is exponentially small.

The code generating @fig-aliasing-visualization is available in:
- `codes/python/ch06_accuracy/aliasing_demo.py`
- `codes/matlab/ch06_accuracy/aliasing_demo.m`

== Accuracy of Spectral Differentiation <sec-accuracy>

=== From Fourier Coefficients to Differentiation Error

We now connect the decay of Fourier coefficients to the accuracy of spectral differentiation. The key observation is that differentiation in Fourier space is multiplication by $ii k$:
$ hat(f')_k = ii k hat(f)_k. $

For the discrete approximation, we have $hat(f')_"discrete" = ii k tilde(f)_k$. The error in the discrete derivative is therefore controlled by the aliasing error in the Fourier coefficients, amplified by the factor $k$.

#block(
  fill: rgb("#142D6E").lighten(92%),
  stroke: (left: 3pt + rgb("#142D6E")),
  inset: (left: 12pt, y: 10pt, right: 10pt),
  width: 100%,
)[
*Theorem 3* (Effect of discretization on the Fourier transform). Let $u in L^2(RR)$ have a first derivative of bounded variation, and let $v$ be the grid function on $h ZZ$ defined by $v_j = u(x_j)$. Then:

*(a)* If $u$ has $p - 1$ continuous derivatives with a $p$-th derivative of bounded variation, then
$ |hat(v)_k - hat(u)_k| = O(h^(p+1)) quad "as" h arrow 0. $

*(b)* If $u$ has infinitely many continuous derivatives, then
$ |hat(v)_k - hat(u)_k| = O(h^m) quad "as" h arrow 0 $
for every $m gt.eq.slant 0$.

*(c)* If $u$ is analytic in a strip $|"Im"(z)| < a$, then
$ |hat(v)_k - hat(u)_k| = O(e^(-pi (a - epsilon) \/ h)) quad "as" h arrow 0 $
for every $epsilon > 0$.
]

The corresponding result for differentiation error is:

#block(
  fill: rgb("#142D6E").lighten(92%),
  stroke: (left: 3pt + rgb("#142D6E")),
  inset: (left: 12pt, y: 10pt, right: 10pt),
  width: 100%,
)[
*Theorem 4* (Accuracy of spectral differentiation). Let $u in L^2(RR)$ have a $nu$-th derivative ($nu gt.eq.slant 1$) of bounded variation, and let $w$ be the $nu$-th spectral derivative of $u$ on the grid $h ZZ$. Then:

*(a)* If $u$ has $p - 1$ continuous derivatives with a $p$-th derivative of bounded variation ($p gt.eq.slant nu + 1$), then
$ |w_j - u^((nu))(x_j)| = O(h^(p - nu)) quad "as" h arrow 0. $

*(b)* If $u$ has infinitely many continuous derivatives, then
$ |w_j - u^((nu))(x_j)| = O(h^m) quad "as" h arrow 0 $
for every $m gt.eq.slant 0$.

*(c)* If $u$ is analytic in a strip $|"Im"(z)| < a$, then
$ |w_j - u^((nu))(x_j)| = O(e^(-pi (a - epsilon) \/ h)) quad "as" h arrow 0 $
for every $epsilon > 0$.

*(d)* If $u$ is band-limited with $hat(u)_k = 0$ for $|k| > pi\/h$, then
$ w_j = u^((nu))(x_j) $
exactly.
]

*Remark on the case $nu = 0$.* Theorem 4 requires $nu gt.eq.slant 1$ because the case $nu = 0$ is trivial: the "zeroth spectral derivative" is simply the sampled function values $w_j = v_j = u(x_j)$, so the pointwise error $|w_j - u(x_j)| = 0$ exactly. The non-trivial question for $nu = 0$ concerns _interpolation error_: how well does the trigonometric interpolant $p_(N)(x)$ approximate $u(x)$ at points _between_ grid points? This is controlled by Theorem 3, which bounds the Fourier coefficient error. The distinction is that Theorem 3 measures error in Fourier space, while Theorem 4 measures differentiation error in physical space at grid points.

=== Connection to Chapter 5

These theorems explain the convergence behavior observed in @sec-etude-convergence. The test function $u(x) = 1\/(2 + sin(x))$ is analytic on the real line, with nearest singularities in the complex plane at $x = -pi\/2 plus.minus ii dot "arcsinh"(2)$. The distance to the real axis is $a = "arcsinh"(2) approx 1.44$. By Theorem 4(c), the spectral differentiation error decays as $O(e^(-a N))$ for $N$ grid points on $[0, 2 pi)$.

The finite difference methods in Chapter 5 achieve only algebraic convergence ($O(N^(-2))$, $O(N^(-4))$, $O(N^(-6))$) because they use only local information, while spectral methods exploit global smoothness to achieve exponential convergence.

=== Computational Verification

@fig-convergence-rates demonstrates the convergence rates predicted by Theorem 4 for our three test functions:

- $|sin(x)|^3$: Algebraic convergence $O(N^(-3))$ because the function has limited smoothness.
- $1\/(1 + sin^2(x\/2))$: Geometric convergence $O(c^(-N))$ because the function is analytic in a strip.
- $exp(sin(x))$: Super-geometric convergence because the function is entire.

#figure(
  image("../figures/ch06/python/convergence_rates.pdf", width: 95%),
  caption: [Spectral differentiation convergence rates for three functions of increasing smoothness. The error is measured in the maximum norm $norm(f' - D bold(f))_infinity$. Functions with greater smoothness achieve faster convergence, exactly as predicted by Theorem 4.],
) <fig-convergence-rates>

The following Python code computes the differentiation error:

```python
def spectral_diff_error(f, f_deriv, N):
    """Compute max error in spectral differentiation."""
    D, x = spectral_diff_periodic(N)
    error = np.max(np.abs(D @ f(x) - f_deriv(x)))
    return error
```

The equivalent MATLAB implementation:

```matlab
function error = spectral_diff_error(f, f_deriv, N)
    [D, x] = spectral_diff_periodic(N);
    error = max(abs(D * f(x) - f_deriv(x)));
end
```

The code generating @fig-convergence-rates is available in:
- `codes/python/ch06_accuracy/convergence_rates.py`
- `codes/matlab/ch06_accuracy/convergence_rates.m`

== Summary

This chapter has established the theoretical foundations for spectral accuracy:

+ *Smoothness implies spectral decay*: The Fourier coefficients of a function decay at a rate determined by its smoothness. Algebraic decay for functions with finite regularity; exponential decay for analytic functions; super-exponential decay for entire functions.

+ *Aliasing is controlled by decay*: When a function is sampled on a discrete grid, high frequencies fold onto low frequencies. If the high-frequency coefficients are small, this aliasing introduces negligible error.

+ *Spectral methods exploit smoothness*: By using global information (all grid points), spectral methods capture the smoothness of a function and achieve convergence rates that match the decay of its Fourier coefficients.

+ *Exponential convergence for analytic functions*: For functions that are analytic in a strip containing the real axis, spectral methods achieve exponential convergence. The rate is controlled by the distance to the nearest singularity in the complex plane.

The theorems in this chapter provide not only convergence rates but also insight into _when_ spectral methods will excel (smooth functions) and when they may struggle (functions with limited regularity). In the next chapter, we extend these ideas to non-periodic problems on bounded domains, where Chebyshev methods take center stage.
