// textbook/chapters/smoothness_accuracy.typ
// Chapter 6: Smoothness and Spectral Accuracy
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Last modified: January 2026

#import "../styles/template.typ": dropcap

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
$ f(x) = sum_(k = -infinity)^infinity hat(f)_k e^(i k x), $ <eq-fourier-series>
where the Fourier coefficients are
$ hat(f)_k = frac(1, 2 pi) integral_0^(2 pi) f(x) e^(-i k x) d x. $ <eq-fourier-coeff>

The wavenumber $k$ corresponds to oscillations with frequency $|k|$ per period. A smooth function, by definition, changes gradually; it has no abrupt jumps or corners. Such a function cannot have significant energy at high frequencies, for high-frequency waves oscillate rapidly. Therefore, the coefficients $hat(f)_k$ must decay as $|k| arrow infinity$.

Now suppose we sample $f$ at $N$ equispaced points, producing the grid function $v_j = f(x_j)$ where $x_j = 2 pi j \/ N$. The discrete Fourier transform of this grid function involves only $N$ coefficients, corresponding to wavenumbers $|k| lt.eq N\/2$. What happened to the higher wavenumbers?

The answer is _aliasing_. The discrete grid cannot distinguish between a wave with wavenumber $k$ and a wave with wavenumber $k + N$, since $e^(i k x_j) = e^(i(k + N) x_j)$ for all grid points $x_j$. High frequencies masquerade as low frequencies; they "fold" into the resolved spectral range. The discrete Fourier coefficient $tilde(f)_k$ is not equal to the continuous coefficient $hat(f)_k$, but rather to the sum of all coefficients that alias to $k$:
$ tilde(f)_k = sum_(j = -infinity)^infinity hat(f)_(k + j N). $ <eq-aliasing-preview>

This is the _aliasing formula_, and it reveals the key insight: if $hat(f)_k$ decays rapidly, the aliased contributions $hat(f)_(k plus.minus N), hat(f)_(k plus.minus 2 N), dots.h.c$ are negligible. The discrete coefficients $tilde(f)_k$ are then excellent approximations to the continuous coefficients $hat(f)_k$, and spectral methods achieve their remarkable accuracy.

=== Roadmap of This Chapter

We shall make this heuristic argument precise through four theorems:

- *Theorem 1* classifies functions by their smoothness and establishes the corresponding decay rates of their Fourier coefficients.
- *Theorem 2* (the aliasing formula) quantifies how high frequencies contaminate low frequencies upon discretization.
- *Theorems 3 and 4* combine these results to bound the errors in spectral interpolation and differentiation.

The chapter culminates in a computational étude: the quantum harmonic oscillator, an eigenvalue problem whose solutions are entire functions. This example demonstrates spectral accuracy in action, achieving machine precision with modest computational effort.

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

*(a)* If $u$ has $p - 1$ continuous derivatives in $L^2(RR)$ for some $p gt.eq 0$ and a $p$-th derivative of bounded variation, then
$ hat(u)(k) = O(|k|^(-p-1)) quad "as" |k| arrow infinity. $

*(b)* If $u$ has infinitely many continuous derivatives in $L^2(RR)$, then
$ hat(u)(k) = O(|k|^(-m)) quad "as" |k| arrow infinity $
for every $m gt.eq 0$.

*(c)* If there exist $a, c > 0$ such that $u$ can be extended to an analytic function in the complex strip $|"Im"(z)| < a$ with bounded $L^2$ norm along horizontal lines, then
$ hat(u)(k) = O(e^(-a |k|)) quad "as" |k| arrow infinity. $

*(d)* If $u$ is entire (analytic throughout $CC$) and satisfies $|u(z)| = o(e^(a |z|))$ as $|z| arrow infinity$ for some $a > 0$, then $hat(u)$ has compact support: $hat(u)(k) = 0$ for all $|k| > a$.
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
    [Bounded variation], [$O(k^(-1))$], [Step function],
    [$p$ derivatives, $p$-th in BV], [$O(k^(-p-1))$], [B-splines],
    [$C^infinity$ (infinitely smooth)], [$O(k^(-m))$ for all $m$], [Bump functions],
    [Analytic in strip], [$O(e^(-a k))$], [$1\/(1 + x^2)$],
    [Entire], [Faster than any exponential], [$e^(-x^2)$],
    [Band-limited], [Compact support], [$"sinc"(x)$],
    table.hline(stroke: 1.5pt),
  ),
  caption: [The smoothness hierarchy and corresponding Fourier decay rates. More smoothness implies faster decay.],
) <tbl-smoothness-hierarchy>

=== Intuition via Integration by Parts

The algebraic decay rates in parts (a) and (b) can be understood through integration by parts. For a smooth periodic function $f$ with period $2 pi$, we have
$ hat(f)_k = frac(1, 2 pi) integral_0^(2 pi) f(x) e^(-i k x) d x = frac(1, i k) dot frac(1, 2 pi) integral_0^(2 pi) f'(x) e^(-i k x) d x = frac(1, i k) hat(f')_k. $
Each integration by parts gains a factor of $k^(-1)$ in the decay rate. If $f$ has $p$ derivatives, we can integrate by parts $p$ times:
$ hat(f)_k = frac(1, (i k)^p) hat(f^((p)))_k. $ <eq-integration-by-parts>
If $f^((p))$ has bounded variation (so its Fourier coefficients are $O(k^(-1))$), then $hat(f)_k = O(k^(-p-1))$.

=== Intuition via Complex Analysis

The exponential decay in part (c) arises from complex analysis. If $f(z)$ is analytic in a strip $|"Im"(z)| < a$, we can shift the contour of integration in the Fourier integral:
$ hat(f)(k) = integral_(-infinity)^infinity f(x) e^(-i k x) d x = integral_(-infinity)^infinity f(x + i sigma) e^(-i k (x + i sigma)) d x = e^(k sigma) integral_(-infinity)^infinity f(x + i sigma) e^(-i k x) d x $
for any $|sigma| < a$. Taking $sigma = a - epsilon$ for small $epsilon > 0$ (and $k > 0$) gives $|hat(f)(k)| lt.eq C e^(-a k)$. The distance from the real axis to the nearest singularity controls the decay rate.

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
$ hat(v)(k) = sum_(j = -infinity)^infinity hat(u)(k + 2 pi j \/ h). $ <eq-aliasing-formula>
]

For a periodic function sampled at $N$ equispaced points with spacing $h = 2 pi \/ N$, this becomes:
$ tilde(f)_k = sum_(j = -infinity)^infinity hat(f)_(k + j N). $ <eq-aliasing-periodic>

=== Spectral Folding

The aliasing formula has a vivid geometric interpretation: the frequency axis "folds" onto itself at the Nyquist frequency. All frequencies that differ by multiples of $N$ become indistinguishable on the discrete grid.

Consider a wave $e^(i k x)$ sampled at $x_j = 2 pi j \/ N$. At the grid points,
$ e^(i k x_j) = e^(2 pi i k j \/ N) = e^(2 pi i (k + N) j \/ N) = e^(i (k + N) x_j). $
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
$ |tilde(f)_k - hat(f)_k| lt.eq sum_(j eq.not 0) |hat(f)_(k + j N)| lt.eq 2 sum_(m = N\/2)^infinity |hat(f)_m|. $ <eq-aliasing-error>

If $hat(f)_k = O(k^(-p-1))$, this sum is $O(N^(-p))$. If $hat(f)_k = O(e^(-a k))$, the sum is $O(e^(-a N\/2))$, which is exponentially small.

The code generating @fig-aliasing-visualization is available in:
- `codes/python/ch06_accuracy/aliasing_demo.py`
- `codes/matlab/ch06_accuracy/aliasing_demo.m`

== Accuracy of Spectral Differentiation <sec-accuracy>

=== From Fourier Coefficients to Differentiation Error

We now connect the decay of Fourier coefficients to the accuracy of spectral differentiation. The key observation is that differentiation in Fourier space is multiplication by $i k$:
$ hat(f')(k) = i k hat(f)(k). $

For the discrete approximation, we have $hat(f')_"discrete" = i k tilde(f)_k$. The error in the discrete derivative is therefore controlled by the aliasing error in the Fourier coefficients, amplified by the factor $k$.

#block(
  fill: rgb("#142D6E").lighten(92%),
  stroke: (left: 3pt + rgb("#142D6E")),
  inset: (left: 12pt, y: 10pt, right: 10pt),
  width: 100%,
)[
*Theorem 3* (Effect of discretization on the Fourier transform). Let $u in L^2(RR)$ have a first derivative of bounded variation, and let $v$ be the grid function on $h ZZ$ defined by $v_j = u(x_j)$. Then:

*(a)* If $u$ has $p - 1$ continuous derivatives with a $p$-th derivative of bounded variation, then
$ |hat(v)(k) - hat(u)(k)| = O(h^(p+1)) quad "as" h arrow 0. $

*(b)* If $u$ has infinitely many continuous derivatives, then
$ |hat(v)(k) - hat(u)(k)| = O(h^m) quad "as" h arrow 0 $
for every $m gt.eq 0$.

*(c)* If $u$ is analytic in a strip $|"Im"(z)| < a$, then
$ |hat(v)(k) - hat(u)(k)| = O(e^(-pi (a - epsilon) \/ h)) quad "as" h arrow 0 $
for every $epsilon > 0$.
]

The corresponding result for differentiation error is:

#block(
  fill: rgb("#142D6E").lighten(92%),
  stroke: (left: 3pt + rgb("#142D6E")),
  inset: (left: 12pt, y: 10pt, right: 10pt),
  width: 100%,
)[
*Theorem 4* (Accuracy of spectral differentiation). Let $u in L^2(RR)$ have a $nu$-th derivative ($nu gt.eq 1$) of bounded variation, and let $w$ be the $nu$-th spectral derivative of $u$ on the grid $h ZZ$. Then:

*(a)* If $u$ has $p - 1$ continuous derivatives with a $p$-th derivative of bounded variation ($p gt.eq nu + 1$), then
$ |w_j - u^((nu))(x_j)| = O(h^(p - nu)) quad "as" h arrow 0. $

*(b)* If $u$ has infinitely many continuous derivatives, then
$ |w_j - u^((nu))(x_j)| = O(h^m) quad "as" h arrow 0 $
for every $m gt.eq 0$.

*(c)* If $u$ is analytic in a strip $|"Im"(z)| < a$, then
$ |w_j - u^((nu))(x_j)| = O(e^(-pi (a - epsilon) \/ h)) quad "as" h arrow 0 $
for every $epsilon > 0$.

*(d)* If $u$ is band-limited with $hat(u)(k) = 0$ for $|k| > pi\/h$, then
$ w_j = u^((nu))(x_j) $
exactly.
]

=== Connection to Chapter 5

These theorems explain the convergence behavior observed in @sec-etude-convergence. The test function $u(x) = 1\/(2 + sin(x))$ is analytic on the real line, with nearest singularities in the complex plane at $x = -pi\/2 plus.minus i dot "arcsinh"(2)$. The distance to the real axis is $a = "arcsinh"(2) approx 1.44$. By Theorem 4(c), the spectral differentiation error decays as $O(e^(-a N))$ for $N$ grid points on $[0, 2 pi)$.

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

== Computational Étude: The Quantum Harmonic Oscillator <sec-harmonic-oscillator>

We close this chapter with an example that demonstrates spectral accuracy for a physically important eigenvalue problem. The _quantum harmonic oscillator_ is described by the time-independent Schrödinger equation:
$ -u'' + x^2 u = lambda u, quad x in RR. $ <eq-harmonic-oscillator>

This equation arises throughout physics: in quantum mechanics (the harmonic potential), in vibration analysis (normal modes), and in probability theory (Hermite functions).

=== Exact Solution

The eigenvalues of @eq-harmonic-oscillator are well known:
$ lambda_n = 2 n + 1, quad n = 0, 1, 2, dots.h.c $

The corresponding eigenfunctions are the _Hermite functions_:
$ u_n(x) = H_n(x) e^(-x^2 \/ 2), $
where $H_n$ is the $n$-th Hermite polynomial. These functions decay like $e^(-x^2\/2)$ as $|x| arrow infinity$, which is faster than any polynomial. In fact, the Hermite functions are _entire_ functions (analytic throughout $CC$), so by Theorem 4, spectral methods should achieve super-geometric convergence.

=== Numerical Approach

Although the problem is posed on the infinite line $RR$, the rapid decay of the eigenfunctions allows us to truncate to a finite interval $[-L, L]$ for sufficiently large $L$. We impose homogeneous Dirichlet conditions $u(plus.minus L) = 0$, which are automatically satisfied (to machine precision) by the true eigenfunctions.

The numerical scheme uses the periodic spectral method from Chapter 5, rescaled to $[-L, L]$:

1. Set up $N$ equispaced grid points $x_j$ on $[-L, L]$.
2. Construct the second-derivative matrix $D^((2))_N$ rescaled by $(pi\/L)^2$.
3. Form the operator matrix $A = -D^((2))_N + "diag"(x_1^2, dots, x_N^2)$.
4. Compute eigenvalues of $A$ using standard linear algebra.

The core algorithm is remarkably simple. In Python:

```python
def harmonic_oscillator_eigenvalues(N, L):
    """Solve -u'' + x²u = λu on [-L, L]."""
    h = 2 * np.pi / N
    x = h * np.arange(N)
    x = L * (x - np.pi) / np.pi  # Map to [-L, L]

    # Second derivative matrix (rescaled)
    D2 = build_periodic_D2(N) * (np.pi / L)**2

    # Potential matrix
    V = np.diag(x**2)

    # Solve eigenvalue problem
    eigenvalues = np.linalg.eigvalsh(-D2 + V)
    return np.sort(eigenvalues)
```

The equivalent MATLAB implementation:

```matlab
function eigenvalues = harmonic_oscillator(N, L)
    h = 2*pi/N; x = h*(1:N)'; x = L*(x-pi)/pi;

    % Second derivative matrix
    column = [-pi^2/(3*h^2)-1/6, ...
        -.5*(-1).^(1:N-1)./sin(h*(1:N-1)/2).^2];
    D2 = (pi/L)^2 * toeplitz(column);

    eigenvalues = sort(eig(-D2 + diag(x.^2)));
end
```

=== Results

@fig-harmonic-oscillator shows the remarkable performance of the spectral method. The left panel displays the first five eigenfunctions computed with $N = 64$ points on $[-8, 8]$, along with the exact Hermite functions. The agreement is visually perfect.

#figure(
  image("../figures/ch06/python/harmonic_oscillator.pdf", width: 95%),
  caption: [The quantum harmonic oscillator eigenvalue problem. _Left_: First five eigenfunctions computed with $N = 64$ Chebyshev points on $[-8, 8]$ (dots) compared to exact Hermite functions (lines). _Right_: Eigenvalue error versus $N$ for the first four eigenvalues. The error decreases faster than any exponential, reaching machine precision around $N = 36$.],
) <fig-harmonic-oscillator>

The right panel shows the eigenvalue convergence. For $N = 36$ and $L = 8$, the first four eigenvalues are computed to approximately 13-digit accuracy:

#figure(
  table(
    columns: (auto, 1fr, 1fr),
    align: (center, center, center),
    inset: (x: 12pt, y: 8pt),
    stroke: none,
    table.hline(stroke: 1.5pt),
    table.header(
      [*$n$*],
      [*Computed $lambda_n$*],
      [*Error*],
    ),
    table.hline(stroke: 0.75pt),
    [$0$], [$0.99999999999996$], [$4 times 10^(-14)$],
    [$1$], [$3.00000000000003$], [$3 times 10^(-14)$],
    [$2$], [$4.99999999999997$], [$3 times 10^(-14)$],
    [$3$], [$6.99999999999999$], [$1 times 10^(-14)$],
    table.hline(stroke: 1.5pt),
  ),
  caption: [Computed eigenvalues of the harmonic oscillator with $N = 36$ and $L = 8$. The exact values are $lambda_n = 2n + 1$.],
) <tbl-harmonic-eigenvalues>

This is spectral accuracy in action. With just 36 grid points, we have computed eigenvalues to essentially machine precision. A finite difference method would require thousands of points to achieve comparable accuracy.

The code generating @fig-harmonic-oscillator is available in:
- `codes/python/ch06_accuracy/harmonic_oscillator.py`
- `codes/matlab/ch06_accuracy/harmonic_oscillator.m`

== Summary

This chapter has established the theoretical foundations for spectral accuracy:

+ *Smoothness implies spectral decay*: The Fourier coefficients of a function decay at a rate determined by its smoothness. Algebraic decay for functions with finite regularity; exponential decay for analytic functions; super-exponential decay for entire functions.

+ *Aliasing is controlled by decay*: When a function is sampled on a discrete grid, high frequencies fold onto low frequencies. If the high-frequency coefficients are small, this aliasing introduces negligible error.

+ *Spectral methods exploit smoothness*: By using global information (all grid points), spectral methods capture the smoothness of a function and achieve convergence rates that match the decay of its Fourier coefficients.

+ *Exponential convergence for analytic functions*: For functions that are analytic in a strip containing the real axis, spectral methods achieve exponential convergence. The rate is controlled by the distance to the nearest singularity in the complex plane.

The theorems in this chapter provide not only convergence rates but also insight into _when_ spectral methods will excel (smooth functions) and when they may struggle (functions with limited regularity). In the next chapter, we extend these ideas to non-periodic problems on bounded domains, where Chebyshev methods take center stage.
