// textbook/chapters/geometry_of_nodes.typ
// Chapter 4: The Geometry of Nodes
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Email: denys.dutykh@ku.ac.ae
// Homepage: https://www.denys-dutykh.com/
// Last modified: February 2026
#import "../styles/template.typ": dropcap, etude-conclusion, idx, chapter-abstract, exercise, hint-for

= The Geometry of Nodes <ch-geometry>

#chapter-abstract(keywords: [Polynomial interpolation · Runge phenomenon · Chebyshev points · Lebesgue constant · Barycentric interpolation · Potential theory])[
Polynomial interpolation underlies the pseudospectral methods for non-periodic problems, yet its success hinges entirely on where the nodes are placed. This chapter examines the paradox discovered by Carl Runge in 1901: interpolating a perfectly smooth function on an equispaced grid can make the approximation diverge wildly near the endpoints as the degree increases. After developing Lagrange interpolation and the Weierstrass approximation theorem, the chapter explains the Runge phenomenon through potential theory, showing how the limiting distribution of nodes governs whether the interpolant converges. The remedy is the family of Chebyshev points, which cluster near the boundaries as projections of equispaced points on a semicircle and so tame the growth of the interpolation error. The analysis is made quantitative through the Lebesgue constant, which measures interpolation stability, and through the numerically stable barycentric interpolation formula. Computational études on randomly placed nodes and on random angles reinforce, by experiment, why node geometry --- and not polynomial degree alone --- decides accuracy. The chapter thereby motivates the central role of Chebyshev grids in every later non-periodic spectral computation.
]

#dropcap[Before we can differentiate functions numerically using spectral methods, we must first understand how to represent them. Polynomial interpolation is the process of constructing a polynomial that passes through a given set of data points. It is the foundation upon which the pseudospectral methods for non-periodic problems are built (their periodic counterparts rest instead on trigonometric interpolation, as later chapters show). In this chapter, we explore a fascinating paradox: while polynomial interpolation seems entirely straightforward, the choice of interpolation node#idx("interpolation node")s determines whether the method succeeds brilliantly or fails catastrophically.]

The story begins with a surprising discovery by the German mathematician Carl Runge in 1901 @Runge1901. Attempting to approximate a simple, smooth function by interpolating polynomials, Runge found that increasing the polynomial degree made the approximation _worse_, not better. This counterintuitive phenomenon, now bearing his name, reveals deep connections between numerical analysis, complex analysis, and potential theory#idx("potential theory").

Our journey through this "étude in grid geometry" will explain the Runge phenomenon#idx("Runge phenomenon"), develop the theoretical framework of potential theory that predicts where interpolation succeeds or fails, and introduce the Chebyshev points#idx("Chebyshev points") that form the cornerstone of practical spectral methods.

== The Problem: Polynomial Interpolation

=== Lagrange Interpolation

Given $N + 1$ distinct points ${x_0, x_1, dots, x_N}$ on an interval $[a, b]$ and corresponding function values ${f_0, f_1, dots, f_N}$ where $f_j = f(x_j)$, there exists a unique polynomial $p_N (x)$ of degree at most $N$ such that
$ p_N (x_j) = f_j, quad j = 0, 1, dots, N. $ <eq-interpolation-condition>

This interpolating polynomial can be written explicitly using the _Lagrange formula_:
$ p_N (x) = sum_(k=0)^N f_k L_k (x), $ <eq-lagrange-formula>
where the _Lagrange basis polynomial#idx("Lagrange basis polynomial")s_ are
$ L_k (x) = product_(j = 0, j eq.not k)^N frac(x - x_j, x_k - x_j). $ <eq-lagrange-basis>

Each basis polynomial $L_(k) (x)$ has the cardinal property: it equals $1$ at $x_k$ and $0$ at all other nodes. This property ensures that substituting any node $x_j$ into the Lagrange formula yields $f_j$.

=== Interpolation versus Best Approximation

It is important to distinguish interpolation from best approximation. The _Weierstrass Approximation Theorem_, proved by Karl Weierstrass in 1885 @Weierstrass1885 when he was 70 years old, guarantees that any continuous function on $[a, b]$ can be uniformly approximated to arbitrary accuracy by polynomials: if $f$ is continuous on $[-1, 1]$ and $epsilon > 0$ is arbitrary, then there exists a polynomial $p$ such that $norm(f - p)_infinity < epsilon$.

The theorem was independently discovered at about the same time by Carl Runge, the same mathematician whose name we attached to the phenomenon of divergent equispaced interpolation. This coincidence is not accidental: Runge's deep investigations into polynomial approximation led him to both the positive existence result and the negative convergence phenomenon.

Weierstrass's original proof is remarkably elegant, connecting approximation theory to the heat equation. The idea is to extend $f(x)$ to a function $tilde(f)$ with compact support on the real line, then convolve it with the Gaussian kernel
$ phi(x) = frac(e^(-x^2 \/ 4t), sqrt(4 pi t)). $
This convolution solves the diffusion equation $partial u \/ partial t = partial^2 u \/ partial x^2$ and converges uniformly to $f$ as $t arrow 0$. Since the convolution is an entire function (analytic throughout the complex plane), it has a uniformly convergent Taylor series that can be truncated to yield polynomial approximations of arbitrary accuracy.

The theorem stimulated an extraordinary amount of mathematics in the early twentieth century, with many alternative proofs discovered in rapid succession: Picard @Picard1891, Lebesgue @Lebesgue1898 (in his first paper at age 23), Fejér @Fejer1900 (at age 20), Landau @Landau1908, Jackson @Jackson1911, and Bernstein @Bernstein1912, among others.

==== Bernstein Polynomials

The most constructive proof was given by Sergei Bernstein in 1912 @Bernstein1912. For $f in C([0, 1])$, the _Bernstein polynomial#idx("Bernstein polynomial")_ of degree $n$ is defined by
$ B_n (x) = sum_(k=0)^n f(k\/n) binom(n, k) x^k (1 - x)^(n-k). $ <eq-bernstein-polynomial>

Bernstein proved that $B_n (x) arrow f(x)$ uniformly as $n arrow infinity$. The convergence can be understood through a beautiful probabilistic interpretation: $B_n (x)$ represents the expected value of $f$ evaluated at the proportion of heads in $n$ independent coin flips, where each coin has probability $x$ of landing heads. As $n arrow infinity$, the law of large numbers ensures this proportion concentrates near $x$, and the expected value converges to $f(x)$.

Bernstein polynomials provide explicit approximations that converge for _any_ continuous function, though the convergence is typically slow: $O(1\/sqrt(n))$ for Lipschitz continuous functions. This is far inferior to the geometric convergence achieved by Chebyshev interpolation for analytic functions.

==== The Limits of Interpolation

Despite the Weierstrass theorem's guarantee that approximating polynomials exist, a fundamental negative result constrains how we can find them. The _Faber--Bernstein theorem_ @Faber1914 asserts that there is no fixed array of interpolation grids with $1, 2, 3, dots$ points that achieves convergence as $n arrow infinity$ for _all_ continuous functions.

This result might seem discouraging, but it has led to an unfortunate overemphasis on pathological functions in the numerical analysis literature. The practical reality is far more favorable: for functions with even modest smoothness, Chebyshev interpolation converges beautifully. The Weierstrass theorem applies to highly irregular functions like $sin(1\/x) sin(1\/sin(1\/x))$, which oscillates infinitely often near infinitely many points. Such functions are mathematical curiosities, not the smooth solutions to differential equations that arise in scientific computing.

Interpolation constructs a specific polynomial by enforcing exact agreement at the nodes. The hope is that as $N arrow infinity$, the interpolating polynomials $p_N$ converge uniformly to $f$. As we shall see, this hope is fulfilled for some node distributions but dramatically fails for others. The key insight of this chapter is that for smooth functions and well-chosen nodes, interpolation not only succeeds but achieves the optimal approximation rate up to a small constant factor.

=== Equispaced Nodes

The most natural choice of interpolation nodes is equally spaced points:
$ x_j = -1 + frac(2 j, N), quad j = 0, 1, dots, N. $ <eq-equispaced-nodes>

These nodes divide the interval $[-1, 1]$ into $N$ equal subintervals. For low-degree polynomials and well-behaved functions, equispaced interpolation works reasonably well. However, a fundamental problem emerges as $N$ increases.

== The Runge Phenomenon <sec-runge>

=== A Smooth but Troublesome Function

In 1901, Carl Runge @Runge1901 investigated the convergence of equispaced polynomial interpolation in general, and illustrated its central pitfall with a deceptively simple example:
$ f(x) = frac(1, 1 + 25 x^2). $ <eq-runge-function>

This _Runge function#idx("Runge function")_ is infinitely differentiable on the entire real line. Its graph is a smooth bell curve centered at the origin with maximum value $f(0) = 1$ and asymptotic decay to zero as $|x| arrow infinity$.

Runge discovered that polynomial interpolation on equispaced nodes#idx("equispaced nodes") _diverges_ for this function. Rather than improving with increasing $N$, the interpolating polynomials develop wild oscillations near the endpoints $x = plus.minus 1$.

=== Numerical Demonstration

@fig-runge-phenomenon illustrates this phenomenon. For $N = 6$, the interpolant reasonably approximates the function. But for $N = 10$ and $N = 14$, the approximation deteriorates dramatically at the boundaries, with oscillations that grow unboundedly as $N$ increases.

#figure(
  image("../figures/ch04/python/runge_phenomenon.pdf", width: 85%),
  caption: [The Runge phenomenon: polynomial interpolation of $f(x) = 1\/(1 + 25x^2)$ on equispaced nodes. As the polynomial degree increases, oscillations near the boundaries $x = plus.minus 1$ grow unboundedly. The exact function (solid) and interpolants (dashed) are shown for $N = 6$, $10$, and $14$.],
) <fig-runge-phenomenon>

The code generating @fig-runge-phenomenon is available in:
- `codes/python/ch04/runge_phenomenon.py`
- `codes/matlab/ch04/runge_phenomenon.m`
- `codes/julia/ch04/runge_phenomenon.jl`

=== Why Does This Happen?

The Runge phenomenon seems paradoxical: how can adding more information (more interpolation points) make the approximation worse? The answer lies in the _Lebesgue constant#idx("Lebesgue constant")_, which we will explore in @sec-lebesgue. But first, let us develop a deeper understanding through potential theory.

== Theoretical Explanation: Potential Theory <sec-potential-theory>

=== Singularities in the Complex Plane

The key to understanding the Runge phenomenon lies in extending our view from the real line to the complex plane. The Runge function has singularities (poles) where its denominator vanishes:
$ 1 + 25 z^2 = 0 quad arrow.r.double quad z = plus.minus frac(i, 5) = plus.minus 0.2 i. $

These poles lie on the imaginary axis, at a distance of only $0.2$ from the real interval $[-1, 1]$. This proximity is responsible for the divergence of equispaced interpolation.

=== The Potential Function

The convergence of polynomial interpolation is governed by a _potential function_ associated with the node distribution @SaffTotik1997. For a distribution with density $mu(x)$ on $[-1, 1]$, the potential at a point $z$ in the complex plane is:
$ phi(z) = - integral_(-1)^1 mu(x) ln |z - x| dif x. $ <eq-potential>

The _equipotential curve#idx("equipotential curve")s_ $phi(z) = "const"$ form closed contours around the interval $[-1, 1]$. The largest equipotential curve that does not enclose any singularity of $f$ determines the region of convergence.

=== Uniform versus Chebyshev Density

For equispaced nodes, the density is asymptotically uniform: $mu(x) = 1\/2$. The resulting equipotential curves are roughly elliptical, but they extend only a short distance from the interval before reaching the Runge function's poles at $plus.minus 0.2 i$.

For Chebyshev nodes (which we introduce in the next section), the density is:
$ mu(x) = frac(1, pi sqrt(1 - x^2)). $ <eq-chebyshev-density>

This density diverges at the endpoints, concentrating nodes near $x = plus.minus 1$. The corresponding potential simplifies to:
$ phi(z) = ln |z + sqrt(z^2 - 1)| - ln 2 = ln rho - ln 2, $ <eq-chebyshev-potential>
where $rho = |z + sqrt(z^2 - 1)|$. The equipotential curves are _Bernstein ellipse#idx("Bernstein ellipse")s_ with foci at $plus.minus 1$.

@fig-equipotential-curves compares the equipotential curves for both distributions, showing why Chebyshev interpolation succeeds where equispaced interpolation fails.

#figure(
  image("../figures/ch04/python/equipotential_curves.pdf", width: 95%),
  caption: [Equipotential curves in the complex plane. Left: uniform node density (equispaced nodes). Right: Chebyshev density, where equipotentials are Bernstein ellipses. The poles of the Runge function at $z = plus.minus 0.2 i$ are marked with crosses. For Chebyshev interpolation, the critical ellipse passing through the poles determines the convergence rate.],
) <fig-equipotential-curves>

The code generating @fig-equipotential-curves is available in:
- `codes/python/ch04/equipotential_curves.py`
- `codes/matlab/ch04/equipotential_curves.m`
- `codes/julia/ch04/equipotential_curves.jl`

== The Solution: Chebyshev Points <sec-chebyshev-points>

=== Definition and Geometric Construction

The _Chebyshev-Gauss-Lobatto#idx("Chebyshev-Gauss-Lobatto") points_ (often simply called Chebyshev points#idx("Chebyshev points")) are defined by:
$ x_j = cos(frac(j pi, N)), quad j = 0, 1, dots, N. $ <eq-chebyshev-points>

These points have a beautiful geometric interpretation: they are the projections onto the $x$-axis of $N + 1$ equally spaced points on the upper semicircle of the unit circle. @fig-chebyshev-circle illustrates this construction.

#figure(
  image("../figures/ch04/python/chebyshev_points_circle.pdf", width: 80%),
  caption: [Geometric construction of Chebyshev points. Points equally spaced on the upper semicircle (by angle $theta_j = j pi \/ N$) project vertically onto the Chebyshev nodes on the $x$-axis. Equal spacing on the circle maps to clustering near the endpoints $x = plus.minus 1$.],
) <fig-chebyshev-circle>

The code generating @fig-chebyshev-circle is available in:
- `codes/python/ch04/chebyshev_points_circle.py`
- `codes/matlab/ch04/chebyshev_points_circle.m`
- `codes/julia/ch04/chebyshev_points_circle.jl`

=== Node Density and Clustering

The geometric construction reveals why Chebyshev points cluster near the endpoints. Equal angular spacing on the circle corresponds to denser horizontal spacing near $x = plus.minus 1$. The node density is approximately:
$ rho(x) approx frac(N, pi sqrt(1 - x^2)), $ <eq-node-density>
which diverges as $x arrow plus.minus 1$. This clustering counteracts the growth of interpolation error at the boundaries.

=== Chebyshev Interpolation of the Runge Function

@fig-chebyshev-success demonstrates the dramatic improvement when using Chebyshev points instead of equispaced nodes for the Runge function. Where equispaced interpolation diverges, Chebyshev interpolation converges rapidly and uniformly.

#figure(
  image("../figures/ch04/python/chebyshev_success.pdf", width: 95%),
  caption: [Chebyshev interpolation of the Runge function. Left: the function and interpolants for $N = 6$, $10$, and $14$. Right: interpolation error on a logarithmic scale. Unlike equispaced interpolation, the error decreases rapidly with increasing $N$.],
) <fig-chebyshev-success>

Note that this formula produces nodes ordered from right to left: $x_0 = cos(0) = +1$ down to $x_N = cos(pi) = -1$. This ordering is natural for the cosine function and is the standard convention in spectral methods.

The code generating @fig-chebyshev-success is available in:
- `codes/python/ch04/chebyshev_success.py`
- `codes/matlab/ch04/chebyshev_success.m`
- `codes/julia/ch04/chebyshev_success.jl`

== Lagrange Basis Functions and Lebesgue Constants <sec-lebesgue>

=== The Interpolation Operator

To understand why node choice matters so dramatically, we formalize the interpolation process. For any continuous function $u in C([-1,1])$, there exists a unique interpolating polynomial $cal(P)_N (x) in bb(P)_N [bb(R)]$ of degree $N$ satisfying:
$ cal(P)_N (x_j) = u_j equiv u(x_j), quad j = 0, 1, dots, N. $
We denote this interpolating polynomial as $cal(I)_N [u]$, emphasizing its role as a _linear operator_ acting on continuous functions.

The _best approximating polynomial_ $cal(P)^* in bb(P)_N [bb(R)]$ is defined by:
$ norm(u - cal(P)^*)_infinity = inf_(cal(P) in bb(P)_N [bb(R)]) norm(u - cal(P))_infinity. $
In general the interpolant and the best approximation differ, $cal(I)_N [u] eq.not cal(P)^*$; they coincide, for instance, when $u$ is itself a polynomial of degree at most $N$, since interpolation then reproduces $u$ exactly and $u$ is trivially its own best approximation.

=== The Lebesgue Function

Examining the Lagrange basis functions more closely, we define the _Lebesgue function_ as the sum of absolute values of all basis polynomials:
$ Lambda_N (x) = sum_(k=0)^N |L_k (x)|. $ <eq-lebesgue-function>

The _Lebesgue constant_ is its maximum over the interval:
$ Lambda_N = max_(x in [-1, 1]) Lambda_N (x). $ <eq-lebesgue-constant>

=== Error Bounds via Operator Norm

The Lebesgue constant bounds the interpolation error through a beautiful argument. For any continuous $u$ with best approximation $cal(P)^*$:
$ norm(u - cal(I)_N [u])_infinity &= norm(u - cal(P)^* + cal(P)^* - cal(I)_N [u])_infinity \
&equiv norm(u - cal(P)^* + cal(I)_N [cal(P)^*] - cal(I)_N [u])_infinity \
&lt.eq.slant norm(u - cal(P)^*)_infinity + norm(cal(I)_N [cal(P)^*] - cal(I)_N [u])_infinity \
&lt.eq.slant norm(u - cal(P)^*)_infinity + norm(cal(I)_N) dot norm(cal(P)^* - u)_infinity \
&lt.eq.slant (1 + norm(cal(I)_N)) dot norm(u - cal(P)^*)_infinity. $ <eq-error-derivation>

The norm of the interpolation operator $cal(I)_N [dot]$ is defined as:
$ norm(cal(I)_N) := sup_(norm(u)_infinity = 1) norm(cal(I)_N [u])_infinity. $ <eq-operator-norm>
This operator norm is precisely the Lebesgue constant $Lambda_N$ for the node set ${x_j}_(j=0)^N$.

Thus, if $p_N^*$ is the best polynomial approximation of degree $N$ to $f$, and $p_N$ is the Lagrange interpolant, then:
$ norm(f - p_N)_infinity lt.eq.slant (1 + Lambda_N) norm(f - p_N^*)_infinity = (1 + Lambda_N) E_N (f), $ <eq-error-bound>
where $E_N (f)$ denotes the best approximation error.

This bound reveals the problem with equispaced nodes: even if $E_N (f)$ decreases with $N$ (as guaranteed by the Weierstrass theorem), the factor $(1 + Lambda_N)$ may grow so fast that the product increases.

=== Asymptotic Growth Rates

The growth rate of $Lambda_N$ depends critically on the node distribution:

- *Chebyshev nodes*: $Lambda_N^"Ch" = frac(2, pi) ln N + O(1)$ (logarithmic growth) @Brutman1978
- *Legendre nodes*#footnote[The _Legendre points_ are the roots of the Legendre polynomial $P_(N+1)$; in the Gauss--Lobatto variant used for interpolation on a closed interval, the interior nodes are the roots of $P'_N$ together with the endpoints $plus.minus 1$. They are developed fully in @ch-quadrature; here we need only their Lebesgue-constant growth.]: $Lambda_N^"Leg" = O(sqrt(N))$ (slow algebraic growth)
- *Equispaced nodes*: $Lambda_N^"eq" approx frac(2^(N+1), e N ln N)$ (exponential growth!) @Turetskii1940

The exponential growth of the Lebesgue constant for equispaced nodes explains the Runge phenomenon: even though the best approximation error decreases geometrically for smooth functions, the exponentially growing factor $(1 + Lambda_N)$ eventually dominates.

*Remark.* As shown by Vértesi @Vertesi1990, the Lebesgue constant for Chebyshev nodes is remarkably close to the smallest possible Lebesgue constant for _any_ node distribution:
$ Lambda_N^"Ch" = frac(2, pi) (ln N + gamma + ln frac(8, pi)) + o(1), $
$ Lambda_N^"min" = frac(2, pi) (ln N + gamma + ln frac(4, pi)) + o(1), $
where $gamma approx 0.5772156649...$ is the Euler--Mascheroni constant. The difference between these two expressions is only $frac(2, pi) ln 2 approx 0.44$, demonstrating that Chebyshev nodes are essentially optimal for polynomial interpolation.

=== Derivation of the Lebesgue Constant for Equi-Spaced Interpolation

We now derive the asymptotic formula for the Lebesgue constant of equi-spaced nodes. Let $L_k (x)$ denote the fundamental Lagrange polynomials and define the equidistant Lebesgue constant by
$ 1 + Lambda_N^"eq" = max_(x in [-1, 1]) sum_(k = 0)^N abs(L_k (x)). $

The interpretation is the usual one: if the interpolation error in the max norm satisfies $norm(f - p_N)_infinity lt.eq.slant epsilon$, then the data values fed into the interpolation formula may differ from the exact polynomial values $p_N (x_j)$ by at most $epsilon$ at each node $x_j$, and the factor $1 + Lambda_N^"eq"$ is the largest possible amplification of these perturbations by Lagrange interpolation#idx("Lagrange interpolation").

The size of the interval does not affect the value of the Lebesgue constant, so for convenience we work with equi-spaced nodes on $[0, N]$ instead of $[-1, 1]$:
$ x_j = j, quad j = 0, 1, dots, N. $

The Lebesgue constant on this interval is $1 + Lambda_N^"eq" = max_(x in [0, N]) sum_(k = 0)^N abs(L_k (x))$. By the symmetry $x arrow.r.bar N - x$, which maps node $j$ to node $N - j$, the Lebesgue function satisfies $lambda(N - x) = lambda(x)$, so the maximum over $[0, N]$ equals the maximum over $[0, N \/ 2]$. Furthermore, for equi-spaced nodes the Lebesgue function is largest in the first subinterval $[0, 1]$, where $x$ lies to one side of all the remaining nodes and the nodal polynomial $product_(j = 0)^N abs(x - x_j)$ is maximised (this is the same endpoint concentration that drives the Runge phenomenon). We may therefore restrict the optimisation to $x in [0, 1]$.

For these nodes the $k$-th fundamental Lagrange polynomial can be written as
$ L_k (x) = product_(j = 0, j eq.not k)^N frac(x - j, k - j) = frac(x (x - 1) dots.c (x - k + 1) (x - k - 1) dots.c (x - N), k (k - 1) dots.c 1 dot (-1) dots.c (k - N)). $

Hence
$ 1 + Lambda_N^"eq" = max_(x in [0, 1]) sum_(k = 0)^N abs(frac(x (x - 1) dots.c (x - k + 1) (x - k - 1) dots.c (x - N), k (k - 1) dots.c 1 dot (-1) dots.c (k - N))). $

We now multiply numerator and denominator of each term by $x - k$. This yields
$ x (x - 1) dots.c (x - k + 1) (x - k) (x - k - 1) dots.c (x - N) = x (1 - x) (2 - x) dots.c (N - x), $
so that
$ abs(L_k (x)) = frac(x (1 - x) (2 - x) dots.c (N - x), abs(k - x) k! (N - k)!). $

Consequently,
$ 1 + Lambda_N^"eq" = max_(x in [0, 1]) sum_(k = 0)^N frac(x (1 - x) (2 - x) dots.c (N - x), abs(k - x) k! (N - k)!). $

Using the Gamma function and repeated application of $Gamma(z + 1) = z Gamma(z)$, we obtain
$ (N - x) (N - 1 - x) dots.c (1 - x) = frac(Gamma(N + 1 - x), Gamma(1 - x)), $
and therefore
$ x (1 - x) (2 - x) dots.c (N - x) = x frac(Gamma(N + 1 - x), Gamma(1 - x)). $

Thus the exact expression becomes
$ 1 + Lambda_N^"eq" = max_(x in [0, 1]) x frac(Gamma(N + 1 - x), Gamma(1 - x)) sum_(k = 0)^N frac(1, abs(k - x) k! (N - k)!). $

Up to this point everything is exact. To estimate $Lambda_N^"eq"$ for large $N$, we introduce several standard approximations.

We first use the ratio asymptotics for the Gamma function:
$ frac(Gamma(N + 1 - x), Gamma(N + 1)) tilde N^(-x) quad (N arrow infinity), $
which gives
$ Gamma(N + 1 - x) tilde N! N^(-x). $

Next observe that, as $N$ grows, the binomial weights $1 \/ (k! (N - k)!)$ are sharply peaked near $k = N \/ 2$. In this region $abs(k - x) approx N \/ 2$ (the optimal $x$ will turn out to be small), so we approximate
$ frac(1, abs(k - x)) approx frac(2, N) $
and factor this out of the sum:
$ sum_(k = 0)^N frac(1, abs(k - x) k! (N - k)!) approx frac(2, N) sum_(k = 0)^N frac(1, k! (N - k)!). $

Using the binomial identity
$ sum_(k = 0)^N frac(N!, k! (N - k)!) = 2^N, $
we obtain
$ sum_(k = 0)^N frac(1, k! (N - k)!) = frac(1, N!) sum_(k = 0)^N binom(N, k) = frac(2^N, N!), $
and hence
$ 1 + Lambda_N^"eq" approx max_(x in [0, 1]) x frac(N! N^(-x), Gamma(1 - x)) dot frac(2, N) dot frac(2^N, N!). $

After cancelling $N!$ we arrive at
$ 1 + Lambda_N^"eq" approx max_(x in [0, 1]) frac(2^(N + 1), N) frac(x N^(-x), Gamma(1 - x)). $

To understand the dominant dependence on $N$ we analyse the factor $x N^(-x)$ for large $N$. Consider
$ phi(x) = x N^(-x) quad "for" quad x > 0. $

Then
$ ln phi(x) = ln x - x ln N, quad frac(dif, dif x) ln phi(x) = frac(1, x) - ln N. $

The maximum occurs where this derivative vanishes:
$ frac(1, x_*) - ln N = 0 quad arrow.r.double quad x_* = frac(1, ln N). $

Evaluating $phi$ at $x_*$ gives
$ phi(x_*) = x_* N^(-x_*) = frac(1, ln N) exp(-x_* ln N) = frac(1, ln N) e^(-1) = frac(1, e ln N). $

For large $N$ this maximiser lies in $(0, 1)$, so it is admissible. Moreover $Gamma(1 - x) arrow 1$ as $x arrow 0$, hence
$ Gamma(1 - x_*) approx 1, $
and the prefactor $1$ in $1 + Lambda_N^"eq"$ is negligible compared with the rapidly growing right-hand side. Thus
$ Lambda_N^"eq" approx frac(2^(N + 1), N) dot frac(1, e ln N) = frac(2^(N + 1), e N ln N) = O(frac(2^N, N ln N)). $ <eq-lebesgue-equispaced-asymptotic>

A slightly sharper estimate follows from the expansion
$ Gamma(1 - x) = 1 + gamma x + O(x^2), $
where $gamma$ is the Euler--Mascheroni constant. Substituting $x_* = 1 \/ ln N$ yields
$ Gamma(1 - x_*) = 1 + frac(gamma, ln N) + O(frac(1, (ln N)^2)), $
so at leading order in $1 \/ ln N$ we obtain
$ Lambda_N^"eq" approx frac(2^(N + 1), e N (gamma + ln N)). $

This recovers the classical asymptotic growth of the Lebesgue constant for equi-spaced interpolation nodes: it increases essentially like $2^N \/ (N ln N)$ as $N arrow infinity$.

=== Lebesgue Points and Multi-Dimensional Considerations

The nodes that _minimize_ the Lebesgue constant $Lambda_N arrow min$ for a given polynomial space $bb(P)_N [bb(R)]$ are called _Lebesgue points_. Finding these optimal nodes is a challenging optimization problem that has been solved only for moderate values of $N$ in one dimension.

In multiple dimensions, the situation becomes considerably more complex. The Lebesgue constant depends not only on the node distribution ${bold(x)_j}_(j=0)^N$ but also on the _shape of the domain_ $cal(U)$. The estimation of the Lebesgue constant in multi-dimensional non-Cartesian domains is a problem that remains essentially open today. The most precious information would be to find the node distribution over a triangle that minimizes the Lebesgue constant---knowledge that would be crucial for the design of new spectral elements @Hesthaven2007.

On _quadrilateral_ domains, the best known point distributions are tensor products of Gauss--Lobatto#footnote[Rehuel Lobatto (1797--1866) was a Dutch mathematician born into a Portuguese family.] or Chebyshev nodes. These tensor product constructions inherit the favorable properties of their one-dimensional counterparts.

On _triangular_ domains, the situation is more subtle. The _Fekete points_#footnote[Michael Fekete (1886--1957) was a Hungarian mathematician who completed his PhD under the supervision of Lipót Fejér. Fekete also provided private tutoring to János Neumann, known today as John von Neumann.] currently represent the best known choice. Fekete points are defined as the nodes that maximize the determinant of the Vandermonde matrix---equivalently, they maximize the volume of the interpolation simplex in function space. While not provably optimal in the Lebesgue sense, they exhibit excellent numerical properties and remain the standard choice for triangular spectral elements.

=== Visualization of Basis Functions

@fig-lagrange-basis compares the Lagrange basis polynomials for equispaced and Chebyshev nodes. For equispaced nodes, the basis functions near the endpoints develop large oscillations. For Chebyshev nodes, all basis functions remain well-behaved.

#figure(
  image("../figures/ch04/python/lagrange_basis.pdf", width: 95%),
  caption: [Lagrange basis polynomials $L_k (x)$ for $N = 10$. Left: equispaced nodes show large oscillations near the boundaries. Right: Chebyshev nodes yield bounded basis functions throughout the interval.],
) <fig-lagrange-basis>

@fig-lebesgue-functions shows the Lebesgue functions and the growth of Lebesgue constants for the three node distributions.

#figure(
  image("../figures/ch04/python/lebesgue_functions.pdf", width: 95%),
  caption: [Left: Lebesgue functions $Lambda_(N)(x)$ for $N = 10$ with equispaced, Legendre, and Chebyshev nodes. The equispaced case has large peaks near the boundaries. Right: Growth of Lebesgue constants with $N$. The equispaced constant grows exponentially, while Chebyshev grows only logarithmically.],
) <fig-lebesgue-functions>

The exponential growth of the equispaced Lebesgue constant dominates @fig-lebesgue-functions, making it difficult to distinguish the Legendre and Chebyshev curves. @fig-lebesgue-zoom provides a zoomed comparison of these two well-behaved distributions, clearly showing the $O(sqrt(N))$ growth of Legendre versus the superior $O(ln N)$ growth of Chebyshev.

#figure(
  image("../figures/ch04/python/lebesgue_constants_zoom.pdf", width: 75%),
  caption: [Lebesgue constants for Legendre and Chebyshev nodes (equispaced omitted for clarity). Chebyshev nodes achieve the slowest possible growth rate $O(ln N)$, while Legendre nodes grow as $O(sqrt(N))$. Both are far superior to equispaced nodes, but Chebyshev remains the optimal choice.],
) <fig-lebesgue-zoom>

The code generating these figures is available in:
- `codes/python/ch04/lagrange_basis.py`
- `codes/python/ch04/lebesgue_functions.py`
- `codes/python/ch04/lebesgue_constants_zoom.py`
- `codes/matlab/ch04/lagrange_basis.m`
- `codes/matlab/ch04/lebesgue_functions.m`
- `codes/julia/ch04/lagrange_basis.jl`
- `codes/julia/ch04/lebesgue_functions.jl`
- `codes/julia/ch04/lebesgue_constants_zoom.jl`

== Barycentric Interpolation

The Lagrange form is the theoretical bedrock of everything above, but it is not how one should evaluate an interpolant in practice. The _barycentric formula#idx("barycentric interpolation")_ is: it is as accurate as any known method, costs only $O(N)$ operations per evaluation once the nodes are fixed, and, for the Chebyshev points that concern us most, uses weights simple enough to write down by hand. Because the interpolant is the workhorse of every pseudospectral computation in the chapters that follow, the formula repays a careful development.

=== Why Not the Lagrange Formula?

The Lagrange formula, while mathematically elegant, is numerically awkward. Evaluating each basis polynomial $L_k (x)$ costs $O(N)$ operations, so a single evaluation of the interpolant costs $O(N^2)$; and adding one new node forces every basis polynomial to be recomputed from scratch. The individual products can also overflow or underflow, so that the sum is assembled from very large or very small terms and loses digits to cancellation.

=== The Barycentric Formula

The remedy is to factor the common _nodal polynomial_ $ell(x) = product_(j=0)^N (x - x_j)$ out of every Lagrange basis function. Writing $L_k (x) = ell(x) w_k \/ (x - x_k)$ identifies the _barycentric weights_
$ w_j = frac(1, product_(k eq.not j) (x_j - x_k)). $ <eq-barycentric-weights>
Substituting this into $p_N (x) = sum_k L_k (x) f_k$, and using the fact that the constant function $1$ is reproduced exactly, so that $ell(x) sum_k w_k \/ (x - x_k) = 1$, the factor $ell(x)$ cancels between the interpolant and this representation of unity and leaves the _barycentric interpolation formula_ @BerrutTrefethen2004:
$ p_N (x) = frac(sum_(j=0)^N frac(w_j, x - x_j) f_j, sum_(j=0)^N frac(w_j, x - x_j)). $ <eq-barycentric>
Once the weights are precomputed, each evaluation costs only $O(N)$ operations. Two properties make @eq-barycentric the method of choice.

- _Scale invariance._ Multiplying every weight $w_j$ by a common nonzero constant leaves $p_N (x)$ unchanged, since the constant cancels between numerator and denominator. Any factor shared by all the weights may therefore be discarded, which is exactly what renders the Chebyshev weights below so simple.
- _Numerical stability._ @Higham2004 proved that the second-form barycentric formula @eq-barycentric is forward stable for any set of well-separated nodes, the Chebyshev points among them, so the catastrophic cancellation that afflicts the naive Lagrange sum does not arise.

The apparent division by zero when $x$ meets a node $x_j$ is removable: there $p_N (x_j) = f_j$ by construction, and in floating point one simply detects the coincidence and returns the stored nodal value.

=== Weights for Chebyshev Points

Scale invariance pays off at once for the Chebyshev-Gauss-Lobatto points. Carrying out the product @eq-barycentric-weights and discarding the factor common to all $j$ leaves
$ w_j = (-1)^j delta_j, quad "where" delta_j = cases(1\/2 & "if" j = 0 "or" j = N, 1 & "otherwise".) $ <eq-chebyshev-weights>
The weights are, up to sign, all equal, with only the endpoints halved: no products to form, nothing to precompute, no cancellation to fear. This closed form, together with the approximation properties established above, is a large part of why the Chebyshev points are the standard choice for spectral methods; the barycentric formula @eq-barycentric is the tool that evaluates the resulting interpolants and, in the next chapters, differentiates them.

== Convergence Analysis <sec-convergence>

=== Geometric Convergence

For functions analytic in a region containing $[-1, 1]$, Chebyshev interpolation converges geometrically. The relevant region is the _Bernstein ellipse#idx("Bernstein ellipse")_ $cal(E)_rho$: the image of the circle $abs(z) = rho$ (with $rho > 1$) under the Joukowski map $z arrow.r.bar (z + z^(-1)) \/ 2$. It is the ellipse with foci at $plus.minus 1$ and semi-axes $(rho + rho^(-1)) \/ 2$ and $(rho - rho^(-1)) \/ 2$, so that $rho$ measures its size and $cal(E)_rho$ shrinks onto $[-1, 1]$ as $rho arrow.r 1^+$. If $f$ is analytic inside $cal(E)_rho$, then:
$ norm(f - p_N)_infinity lt.eq.slant frac(C, rho^N), $ <eq-geometric-convergence>
where $C$ depends on $f$ and $rho$ but not on $N$.

The parameter $rho$ is set by the location of the nearest singularity of $f$: the larger the ellipse of analyticity, the faster the decay. One caveat sharpens the statement. The clean bound @eq-geometric-convergence holds whenever $f$ is analytic _strictly inside_ $cal(E)_rho$; when the nearest singularity lies exactly _on_ $cal(E)_rho$, the type of that singularity can contribute a factor growing algebraically in $N$, so that the sharp rate is $rho^(-N)$ multiplied by a power of $N$ rather than $rho^(-N)$ alone @Boyd2000. The bare $O(rho^(-N))$ is thus best read as the exponential heart of the rate, not the whole of it.

=== The Runge Function Revisited

For the Runge function with poles at $z = plus.minus 0.2 i$, we can compute:
$ rho = |0.2 i + sqrt((0.2 i)^2 - 1)| = |0.2 i + i sqrt(1.04)| approx 1.22. $

Thus Chebyshev interpolation converges at rate $O(1.22^(-N))$; convergence is geometric but modest due to the poles being close to the real axis. In contrast, equispaced interpolation diverges because its effective $rho < 1$.

The divergence of equispaced interpolation can be understood through the error bound $norm(f - p_N)_infinity lt.eq.slant (1 + Lambda_N) E_N (f)$. For the Runge function, the best approximation error $E_N (f)$ still decreases geometrically---the function is analytic, after all. However, recall from @sec-lebesgue that the equispaced Lebesgue constant grows as $Lambda_N^"eq" approx 2^(N+1) \/ (e N ln N)$. This exponential growth eventually overwhelms the geometric decay of $E_N (f)$, causing the interpolation error $(1 + Lambda_N) E_N (f)$ to grow without bound.

=== Numerical Verification

@fig-convergence compares the convergence behavior for both node distributions. The Chebyshev error decreases geometrically, while the equispaced error grows without bound.

#figure(
  image("../figures/ch04/python/convergence_comparison.pdf", width: 80%),
  caption: [Convergence comparison for the Runge function. The maximum interpolation error is plotted against polynomial degree $N$ on a semilogarithmic scale. Chebyshev interpolation (green) converges at rate $O(rho^(-N))$ with $rho approx 1.22$. Equispaced interpolation (red) diverges exponentially.],
) <fig-convergence>

@fig-convergence-zoom shows the Chebyshev convergence alone, without the divergent equispaced curve that dominates the vertical scale. This reveals the beautiful geometric convergence: the error decreases by a factor of approximately $rho approx 1.22$ with each increase in polynomial degree. The theoretical rate matches the computed errors almost perfectly.

#figure(
  image("../figures/ch04/python/convergence_zoom.pdf", width: 75%),
  caption: [Chebyshev interpolation convergence for the Runge function (equispaced omitted for clarity). The error decreases geometrically at rate $O(rho^(-N))$ with $rho approx 1.22$, matching the theoretical prediction based on the pole locations at $z = plus.minus 0.2 i$.],
) <fig-convergence-zoom>

The code generating these figures is available in:
- `codes/python/ch04/convergence_comparison.py`
- `codes/python/ch04/convergence_zoom.py`
- `codes/matlab/ch04/convergence_comparison.m`
- `codes/julia/ch04/convergence_comparison.jl`
- `codes/julia/ch04/convergence_zoom.jl`

== Computational Étude 4.1: Random Nodes <sec-random-nodes>

Having studied the optimal Chebyshev distribution and the problematic equispaced distribution, a natural question arises: what happens if we choose interpolation nodes _randomly_? This question leads us into the realm of _experimental mathematics_, where we use computation to discover and conjecture mathematical relationships.

=== Motivation

The dramatically different behaviors of equispaced and Chebyshev nodes might suggest that the key factor is _regularity_ or _structure_ in node placement. Perhaps any "reasonable" arrangement would work? To test this hypothesis, we investigate the simplest unstructured choice: nodes drawn uniformly at random from $[-1, 1]$.

This étude is the first half of a controlled two-part experiment. Here we strip away all structure, drawing nodes uniformly at random, to ask whether the endpoint clustering of Chebyshev points is genuinely essential or merely one of many acceptable choices. Uniform sampling changes two things at once, the macroscopic _density_ of the nodes and the _separation_ between neighbours, so on its own it cannot tell us which of the two is decisive. The companion étude of @sec-random-chebyshev supplies the missing control: it restores the correct arcsine density while still sampling independently, and thereby isolates density from separation. Read together, the two études pin down which property of the Chebyshev grid actually matters. Along the way, the experiment illustrates the methodology of computational investigation, using numerical evidence to formulate conjectures about asymptotic behavior.

=== Experimental Setup

For each polynomial degree $N$ from $2$ to $30$, we generate $N + 1$ random points uniformly distributed on $[-1, 1]$, sort them, and compute the Lebesgue constant. We repeat this process $M = 200$ times to obtain statistical estimates of the mean, standard deviation, and range of $Lambda_N$ for random nodes.

=== Results

@fig-lebesgue-random shows the results of our Monte Carlo experiment. The left panel displays the growth of the Lebesgue constant with polynomial degree, including shaded regions showing the statistical variability. The right panel shows the distribution of $Lambda_N$ values for a fixed degree $N = 15$.

#figure(
  image("../figures/ch04/python/lebesgue_random_nodes.pdf", width: 95%),
  caption: [Lebesgue constant for random nodes. Left: growth with polynomial degree $N$, showing mean (solid line) and min-max range (light shading). Right: distribution of $log_(10)(Lambda_(15))$ over $500$ random trials on a logarithmic scale, revealing the enormous spread spanning many orders of magnitude. Vertical lines mark the mean, median, equispaced, and Chebyshev values.],
) <fig-lebesgue-random>

The key observations are striking:

+ *Explosive exponential growth*: Random nodes exhibit exponential growth of the Lebesgue constant, but surprisingly _faster_ than equispaced nodes. The _sample_ mean of $Lambda_N$ grows roughly as $e^(b N)$ with $b approx 1.0$, against $b approx 0.6$ for equispaced nodes, and the median (also marked in @fig-lebesgue-random) grows exponentially as well.

+ *Extreme variability*: There is enormous variability between different random samples. For $N = 15$, $Lambda_N$ can range from around $10^2$ to $10^{10}$ or more, spanning many orders of magnitude. This variability increases dramatically with $N$.

+ *Worse than equispaced*: Counter to naive intuition, random nodes perform _worse_ on average than the structured equispaced distribution. The reason is that random sampling occasionally produces clusters of nearby nodes, creating near-singular interpolation conditions. These worst-case realizations dominate the mean.

=== Empirical Asymptotic Formula

Fitting an exponential model to the mean Lebesgue constants yields an empirical formula:
$ Lambda_N^"rand" approx a dot e^(b N), $ <eq-random-asymptotic>
where the fitted constants are approximately $a approx 0.7$ and $b approx 1.0$, giving:
$ Lambda_N^"rand" approx 0.7 dot 2.7^N. $

Remarkably, the growth rate for random nodes is _faster_ than for equispaced nodes (which grow as $approx 1.8^N$). This counterintuitive result arises because random sampling occasionally produces clusters of nodes that are extremely close together, causing the Lebesgue constant to explode. The mean is dominated by these worst-case realizations.

A caveat is in order about the statistic itself. Because independent sampling can place two nodes arbitrarily close, the Lebesgue constant of a random node set has a very heavy right tail, and it is by no means clear that its population mean is finite; the sample mean over the $M$ trials is best read as a diagnostic of that tail rather than as an estimate of a well-defined expectation. It is for this reason that we lean on the _median_, which is stable and robust, when we speak of the "typical" growth rate, and the exponential fit above should be understood in that light.

The fundamental problem is twofold: (1) the lack of _clustering near the endpoints_, which only carefully designed distributions like Chebyshev points possess, and (2) the risk of _random clustering_ in the interior, which creates severely ill-conditioned interpolation problems.

#etude-conclusion[
  This étude provides strong numerical evidence for a fundamental principle: *the clustering of Chebyshev points near the endpoints is not merely convenient but essential when one seeks to push the polynomial degree toward the spectral limit*. Random nodes, despite their apparent "fairness", fail in two ways --- no boundary resolution to control Lagrange-basis oscillations, and risk of interior clustering that creates catastrophically ill-conditioned systems. The extreme variability in $Lambda_N$ is the most striking finding: equispaced nodes are suboptimal but *predictable*, while random nodes add a catastrophic-or-acceptable lottery on top. (Note: the severe boundary concentration of Chebyshev grids is the price of pushing the polynomial degree to the spectral limit; high-order finite-difference methods on mildly clustered or equispaced grids reach moderate accuracy without it @Fornberg2025.) The étude is also a small lesson in computational mathematics --- a natural hypothesis is tested, unexpected behaviour discovered ("random is *worse* than equispaced"), and the empirical asymptotics quantified by data fitting.
]

The code generating @fig-lebesgue-random is available in:
- `codes/python/ch04/lebesgue_random_nodes.py`
- `codes/matlab/ch04/lebesgue_random_nodes.m`
- `codes/julia/ch04/lebesgue_random_nodes.jl`

== Computational Étude 4.2: Random Angles on the Circle <sec-random-chebyshev>

The previous étude revealed that uniform random nodes on $[-1, 1]$ perform poorly because they lack the endpoint clustering of Chebyshev points. This raises a natural follow-up question: what if we generate random points that _do_ cluster near the endpoints? Specifically, what happens if we generate random angles uniformly on the semicircle and project them via cosine, mimicking the construction of Chebyshev nodes?

=== The Arcsine Distribution

Recall that Chebyshev--Gauss--Lobatto nodes are defined as $x_j = cos(j pi / N)$ for $j = 0, 1, dots, N$. The angles $theta_j = j pi / N$ are _equispaced_ on $[0, pi]$, and the cosine projection creates the characteristic endpoint clustering.

Our experiment generates _random_ angles uniformly on $[0, pi]$ and applies the same cosine projection:
$ theta_k tilde "Uniform"(0, pi), quad x_k = cos(theta_k). $

If $theta$ is uniformly distributed on $[0, pi]$, then $x = cos(theta)$ has the _arcsine distribution_ with probability density function:
$ f(x) = frac(1, pi sqrt(1 - x^2)), quad x in [-1, 1]. $ <eq-arcsine-pdf>
This is precisely the Chebyshev weight function! The arcsine distribution naturally concentrates probability mass near $plus.minus 1$, providing the endpoint clustering we seek.

=== Implementation

We conduct the same Monte Carlo experiment as before: for each $N$ from $2$ to $30$, we generate $200$ random realizations and compute the Lebesgue constant for each.

=== Surprising Results

@fig-lebesgue-random-chebyshev displays the results, which are profoundly counterintuitive.

#figure(
  image("../figures/ch04/python/lebesgue_random_chebyshev.pdf", width: 95%),
  caption: [Lebesgue constant for "random Chebyshev" nodes (uniform angles projected via cosine). Left: growth with polynomial degree, showing explosive exponential behavior with enormous variability. Right: distribution of $log_(10)(Lambda_(15))$ over $500$ trials. Despite having the "correct" arcsine marginal distribution, these nodes perform _catastrophically worse_ than any other distribution studied.],
) <fig-lebesgue-random-chebyshev>

The key observations are striking and instructive:

+ *Catastrophic exponential growth*: Random Chebyshev nodes exhibit the _worst_ performance of any distribution we have studied. The mean Lebesgue constant grows explosively, reaching astronomical values ($10^(20)$ or higher) by $N = 30$.

+ *Extreme variability*: The spread across realizations is staggering, spanning $10$ to $20$ orders of magnitude. Some realizations produce reasonable Lebesgue constants (the minimum values grow polynomially), while others are catastrophically large.

+ *Worse than uniform random*: Counter to our hypothesis, having the "correct" marginal distribution (arcsine) actually makes things _worse_, not better.

=== Understanding the Failure

Why does this happen? The answer lies in the distinction between _marginal distributions_ and _joint distributions_.

When we generate $N + 1$ random angles uniformly on $[0, pi]$, there is no guarantee of minimum separation. Two angles $theta_i$ and $theta_j$ can be arbitrarily close. When this happens, the corresponding nodes $x_i = cos(theta_i)$ and $x_j = cos(theta_j)$ become nearly identical, creating a _near-singular_ interpolation problem.

In contrast, deterministic Chebyshev nodes have angles that are exactly $pi / N$ apart. This _guaranteed minimum separation_ ensures that no two nodes can cluster together.

The mathematical insight is profound: *the success of Chebyshev nodes comes not from the arcsine marginal distribution, but from the deterministic structure that guarantees minimum node separation*.

=== The Lesson

This étude teaches a crucial lesson about interpolation theory:

+ *Marginal distribution is not sufficient*: Having the "right" probability distribution for individual nodes does not guarantee good collective behavior.

+ *Structure matters*: The deterministic, equispaced placement of angles in true Chebyshev nodes provides essential guarantees that random sampling cannot replicate.

+ *Minimum separation is key*: The $O(1/N)$ minimum spacing between Chebyshev angles translates to guaranteed node separation, preventing the catastrophic clustering that plagues random approaches.

This finding reinforces a fundamental principle: in numerical analysis, _structured_ methods often outperform _random_ methods precisely because structure provides guarantees that randomness cannot.

*Remark (repulsion restores the marginal's promise).* The failure documented here is one of _independent_ sampling, not of randomness as such. A point process that preserves the arcsine marginal but forbids near-coincidences supplies exactly the separation that independent sampling lacks. The natural such model is a _determinantal point process#idx("determinantal point process")_, whose points repel one another; on $[-1, 1]$ one can construct one whose marginal is precisely the arcsine (equilibrium) density @BardenetHardy2020. With the near-coincidences suppressed, the pathological clustering behind @fig-lebesgue-random-chebyshev is removed and the growth of the Lebesgue constant is tamed accordingly. This is the exception that proves the rule: what the deterministic Chebyshev construction provides for free, a guaranteed minimum spacing, is exactly what a competitive random model must engineer. Density sets the stage; separation decides the outcome.

#etude-conclusion[
  The result is sobering: *reproducing the arcsine marginal distribution of Chebyshev nodes through random sampling is not enough*. The Lebesgue constants for random cosine-projected points grow *even faster* than for uniform random nodes, with variability spanning many orders of magnitude. The root cause is the absence of *guaranteed minimum separation*: true Chebyshev points space their angles exactly $pi \/ N$ apart, a rigid structure that independent random sampling cannot replicate. The étude reinforces the lesson of Étude 4.1: what makes Chebyshev nodes special is their *deterministic architecture*, not merely their density profile. Two random configurations may share the same average spacing yet differ catastrophically in their worst-case gaps. When an interpolation grid can be chosen freely, structured node sets with provable separation guarantees are the safe default; and when a random model is genuinely wanted, the separation must be engineered in, not left to chance, as the remark above makes precise.
]

The code generating @fig-lebesgue-random-chebyshev is available in:
- `codes/python/ch04/lebesgue_random_chebyshev.py`
- `codes/matlab/ch04/lebesgue_random_chebyshev.m`
- `codes/julia/ch04/lebesgue_random_chebyshev.jl`

== Practical Guidelines and Outlook

=== When to Use Which Grid

Based on the theory developed in this chapter, we can state clear guidelines:

+ *Always prefer Chebyshev points* for polynomial interpolation on a finite interval. The $O(ln N)$ growth of the Lebesgue constant ensures stability.

+ *Avoid equispaced nodes* for high-degree polynomial interpolation. The exponential growth of $Lambda_N$ makes this a losing proposition.

+ *Consider the singularity structure* of the function being approximated. The location of complex singularities determines the convergence rate.

=== Mapping to General Intervals

Chebyshev points are defined on $[-1, 1]$, but can be mapped to any interval $[a, b]$ via the linear transformation:
$ x_"phys" = frac(a + b, 2) + frac(b - a, 2) x_"ref", $ <eq-interval-mapping>
where $x_"ref" in [-1, 1]$ is the reference coordinate.

=== Preview of Differentiation Matrices

The Chebyshev points introduced in this chapter will play a central role in the differentiation matrices we develop in subsequent chapters. The clustering of nodes near the boundaries, far from being a peculiarity, is precisely what enables accurate spectral differentiation. The connection between node distribution, interpolation accuracy, and differentiation stability is one of the beautiful unifying themes of spectral methods.

In the next chapter, we will see how to convert function values at Chebyshev points into accurate approximations of derivatives, building on the geometric insights developed here.

== A non-exhaustive literature overview

The theory of polynomial interpolation and node geometry presented in this chapter draws from a rich intellectual tradition spanning two centuries, from foundational existence theorems to modern algorithmic refinements. This overview traces the principal threads.

*The theoretical bedrock.* The genesis of polynomial approximation theory lies in the work of Karl Weierstrass @Weierstrass1885, who established that continuous functions on a closed interval can be uniformly approximated by polynomials to arbitrary precision. This existence result was subsequently given constructive proofs by several mathematicians in rapid succession: Picard @Picard1891 using integral operators, Lebesgue @Lebesgue1898 through geometric methods, Fejér @Fejer1900 via Cesàro summation of Fourier series, Landau @Landau1908 through convolution with rational kernels, Jackson @Jackson1911 using trigonometric approximation, and most elegantly by Bernstein @Bernstein1912 through the probabilistic construction that now bears his name. However, the negative result of Faber @Faber1914 demonstrated that no single triangular array of interpolation nodes can guarantee convergence for _all_ continuous functions, shifting the focus from universal convergence to convergence for smooth functions on well-chosen grids.

*The Runge phenomenon and its quantification.* The pathology of equispaced interpolation was discovered by Runge @Runge1901, who showed that the seemingly innocuous function $(1 + 25 x^2)^(-1)$ defies polynomial interpolation on uniform grids. The mechanism of failure was quantified through the Lebesgue constant: Turetskii @Turetskii1940 derived the asymptotic formula $Lambda_N^"eq" tilde 2^(N+1) \/ (e N ln N)$ for equispaced nodes, establishing that error amplification grows exponentially. In contrast, Brutman @Brutman1978 proved that the Lebesgue constant for Chebyshev nodes grows only logarithmically, and Vértesi @Vertesi1990 showed that this logarithmic rate is asymptotically optimal among all possible node distributions. Detailed numerical experiments by Luttmann and Rivlin @LuttmannRivlin1965 corroborated these theoretical predictions and revealed fine structural properties of the Lebesgue function. A comprehensive modern survey of Lebesgue constants across different node distributions is provided by Smith @Smith2006.

*The potential-theoretic framework.* The deep explanation for why certain node distributions succeed lies in logarithmic potential theory, developed comprehensively by Saff and Totik @SaffTotik1997. The key insight is that the convergence of polynomial interpolation is governed by equipotential curves in the complex plane: the equilibrium measure $mu_"eq" (x) = 1 \/ (pi sqrt(1 - x^2))$ generates Bernstein ellipses as its level curves, and any node distribution that asymptotically follows this density (such as Chebyshev or Legendre nodes) will guarantee convergence for functions analytic inside the corresponding ellipse. This framework unifies the disparate observations about node geometry into a single elegant theory.

*Algorithmic realization.* The translation of these theoretical insights into stable numerical algorithms required further innovation. Berrut and Trefethen @BerrutTrefethen2004 revived and popularized the barycentric form of Lagrange interpolation, which reduces evaluation cost to $O(N)$ operations while avoiding the catastrophic cancellation inherent in the classical Lagrange formula. Higham @Higham2004 provided a rigorous forward error analysis proving that the barycentric formula is numerically stable whenever the Lebesgue constant is moderate, thereby confirming that the theoretical stability of Chebyshev nodes translates faithfully into finite-precision arithmetic. These algorithmic developments underpin the practical implementations used throughout this book; see @Trefethen2000 and @Boyd2000 for comprehensive treatments.

*Multi-dimensional frontiers.* While one-dimensional interpolation theory is well established, the extension to higher dimensions remains an active area of research. On triangular domains, the Fekete points (which maximize the Vandermonde determinant) provide the natural generalization of Chebyshev nodes and are now standard in high-order spectral element methods @Hesthaven2007. On quadrilateral domains, a breakthrough came with the _Padua points_, discovered by De Marchi, Caliari, and Vianello @DeMarchi2005: the first explicitly known set of unisolvent nodes on the square $[-1, 1]^2$ with a Lebesgue constant growing only as $O(ln^2 N)$, a bound subsequently established by Bos, Caliari, De Marchi, Vianello, and Xu @Bos2006 @Bos2007. These points are generated by the self-intersections of a Lissajous curve and offer a practical alternative to tensor-product grids for bivariate interpolation.

== Summary

This chapter has established that, within the polynomial interpolation framework, the choice of nodes is decisive for spectral accuracy. (We note that alternative frameworks, such as radial basis function interpolation, can achieve spectral accuracy on more general node layouts, including scattered points in multiple dimensions; see the discussion in Chapter 3.) The key findings are:

+ *The Runge phenomenon*: High-degree polynomial interpolation on equispaced nodes diverges near the boundaries, with errors growing exponentially as $N$ increases.

+ *Potential theory explanation*: The convergence rate of polynomial interpolation is governed by the location of the function's singularities in the complex plane and by the potential-theoretic properties of the node distribution.

+ *Chebyshev points*: Nodes distributed as $x_j = cos(j pi \/ N)$ cluster near the boundaries and yield a Lebesgue constant that grows only as $O(ln N)$, ensuring stable, spectrally accurate interpolation.

+ *Lebesgue constant*: This quantity measures the worst-case amplification of interpolation error. Its slow logarithmic growth for Chebyshev nodes, versus exponential growth for equispaced nodes, is the geometric foundation of spectral methods.

+ *Barycentric interpolation*: The barycentric formula provides a numerically stable $O(N)$ algorithm for evaluating the interpolant at any point, avoiding the explicit construction of the Lagrange polynomials.

+ *Structure over randomness*: Random nodes with the correct marginal distribution perform catastrophically worse than deterministic Chebyshev nodes. The guaranteed minimum separation of deterministic grids, not the density profile alone, is what prevents ill-conditioning.

These geometric insights motivate the differentiation matrices developed in the next chapter, where function values at Chebyshev points are converted into accurate spectral approximations of derivatives.

== Exercises <sec-geometry-of-nodes-exercises>

The exercises below progress from pencil-and-paper properties of polynomial interpolation and node geometry, through numerical experiments that reproduce and extend the études of this chapter, to open-ended projects that reach into the research literature on optimal node distributions. The computational problems may be carried out in any of the book's three languages; the named scripts under `codes/` give a starting point.

=== Conceptual Exercises

#exercise(title: [Cardinality and Uniqueness of the Interpolant])[
  Let ${x_0, dots, x_N}$ be $N + 1$ distinct nodes and let $L_k (x)$ be the Lagrange basis polynomials of @eq-lagrange-basis. (a) Prove the cardinal property $L_k (x_j) = delta_(j k)$ directly from the product formula, and deduce that the Lagrange formula @eq-lagrange-formula satisfies the interpolation conditions @eq-interpolation-condition. (b) Prove that the interpolating polynomial of degree at most $N$ is unique: if two such polynomials agree at all $N + 1$ nodes, their difference is a polynomial of degree at most $N$ with $N + 1$ roots, hence identically zero. (c) Conclude that ${L_0, dots, L_N}$ is a basis of $bb(P)_N [bb(R)]$.
] <ex-geo-cardinal-uniqueness>

#exercise(title: [Partition of Unity and a Lower Bound on the Lebesgue Function])[
  (a) By interpolating the constant function $f equiv 1$, which lies in $bb(P)_N [bb(R)]$ and is therefore reproduced exactly, show that the Lagrange basis forms a partition of unity, $sum_(k=0)^N L_k (x) = 1$ for all $x$. (b) Deduce that the Lebesgue function @eq-lebesgue-function satisfies $Lambda_N (x) gt.eq.slant 1$ pointwise, with equality wherever all the $L_k (x)$ share the same sign. (c) Show that $Lambda_N (x_j) = 1$ at every node, so the Lebesgue constant @eq-lebesgue-constant obeys $Lambda_N gt.eq.slant 1$ for every node set, its maximum being attained between nodes.
] <ex-geo-partition-unity>

#exercise(title: [The Near-Best Approximation Inequality])[
  Let $cal(I)_N [u]$ be the interpolation operator and $cal(P)^*$ the best uniform approximation to $u$ in $bb(P)_N [bb(R)]$. (a) Using the projection property $cal(I)_N [cal(P)^*] = cal(P)^*$, justify each step of the derivation @eq-error-derivation and conclude the near-best bound @eq-error-bound, $norm(u - cal(I)_N [u])_infinity lt.eq.slant (1 + Lambda_N) norm(u - cal(P)^*)_infinity$. (b) Explain why the operator norm @eq-operator-norm equals the Lebesgue constant. (c) Using the equispaced asymptotics @eq-lebesgue-equispaced-asymptotic, show that this bound permits divergence even when the best-approximation error $E_N (u)$ decays geometrically.
] <ex-geo-near-best>

#exercise(title: [Bernstein Ellipse and the Convergence Rate])[
  The convergence rate of Chebyshev interpolation is set by the largest Bernstein ellipse $cal(E)_rho$ free of singularities of $f$, as in @eq-geometric-convergence. The Joukowski map $z = (w + w^(-1)) \/ 2$ sends the circle $|w| = rho$ to $cal(E)_rho$, and the parameter is recovered as $rho = |z + sqrt(z^2 - 1)|$, consistent with the Chebyshev potential @eq-chebyshev-potential. (a) For a pole at $z = i alpha$ with $alpha > 0$, show that $rho = alpha + sqrt(1 + alpha^2)$. (b) Apply this to the Runge function @eq-runge-function, whose poles sit at $z = plus.minus 0.2 i$, to recover $rho approx 1.22$. (c) Show that $rho > 1$ for every $alpha > 0$, so Chebyshev interpolation of any such function converges geometrically, however close the pole lies to the interval.
] <ex-geo-bernstein-ellipse>

#hint-for(<ex-geo-bernstein-ellipse>)[Substitute $z = i alpha$ into $rho = |z + sqrt(z^2 - 1)|$ and simplify $sqrt((i alpha)^2 - 1) = i sqrt(alpha^2 + 1)$; the modulus of $i (alpha + sqrt(alpha^2 + 1))$ is the positive real number $alpha + sqrt(alpha^2 + 1)$.]

#exercise(title: [The Arcsine Density from Equispaced Angles])[
  The Chebyshev--Gauss--Lobatto nodes @eq-chebyshev-points arise by projecting equispaced angles through the cosine map. (a) Treating $theta$ as uniform on $[0, pi]$ and setting $x = cos theta$, use the change-of-variables formula for densities to derive the arcsine density @eq-arcsine-pdf, $f(x) = 1 \/ (pi sqrt(1 - x^2))$, and identify it with the Chebyshev weight @eq-chebyshev-density. (b) Show that the implied count of nodes in a small interval $[x, x + dif x]$ reproduces the clustering density @eq-node-density. (c) Explain in one or two sentences why this same density appears as the equilibrium measure of $[-1, 1]$ in the potential @eq-potential.
] <ex-geo-arcsine-density>

#exercise(title: [Barycentric Weights for Chebyshev Points])[
  Starting from the general barycentric weights @eq-barycentric-weights, $w_j = 1 \/ product_(k eq.not j) (x_j - x_k)$, derive the closed form @eq-chebyshev-weights for the Chebyshev--Gauss--Lobatto nodes $x_j = cos(j pi \/ N)$. (a) Observe that the weights enter the barycentric formula @eq-barycentric only through ratios, so any common multiplicative constant may be dropped. (b) Establish the alternating-sign structure $w_j = (-1)^j delta_j$ with the endpoint halving $delta_0 = delta_N = 1 \/ 2$. (c) Explain why these weights#idx("barycentric weights") are immune to the overflow that afflicts the raw products $product_(k eq.not j) (x_j - x_k)$ when $N$ is large.
] <ex-geo-bary-weights>

#exercise(title: [Chebyshev Points as Extrema of the Chebyshev Polynomial])[
  Let $T_N (x) = cos(N arccos x)$ be the degree-$N$ Chebyshev polynomial. (a) Show that the Chebyshev--Gauss--Lobatto points @eq-chebyshev-points are exactly the points of $[-1, 1]$ where $T_N (x) = plus.minus 1$, that is, the $N - 1$ interior extrema of $T_N$ together with the two endpoints. (b) Conclude that these nodes are the zeros of $(1 - x^2) T'_N (x)$. (c) Contrast them with the Chebyshev points of the first kind, the zeros of $T_N$, given by $x_j = cos((2 j + 1) pi \/ (2 N))$ for $j = 0, dots, N - 1$, and note that neither set places two nodes at the same location.
] <ex-geo-cheb-extrema>

#exercise(title: [Bernstein Polynomials Reproduce Linear Functions])[
  Consider the Bernstein polynomial @eq-bernstein-polynomial $B_n (x) = sum_(k=0)^n f(k\/n) binom(n, k) x^k (1 - x)^(n - k)$ on $[0, 1]$. (a) Using the binomial theorem, prove the two identities $B_n [1] = 1$ and $B_n [x] = x$, so that constants and linear functions are reproduced exactly. (b) Interpret $B_n (x)$ as the expectation of $f(S_n \/ n)$, where $S_n$ counts heads in $n$ independent tosses of a coin with heads-probability $x$, and use the law of large numbers to argue uniform convergence $B_n arrow.r f$. (c) Explain why the resulting rate $O(1 \/ sqrt(n))$ for Lipschitz $f$ is far slower than the geometric rate @eq-geometric-convergence of Chebyshev interpolation, so that Bernstein polynomials are prized for theory rather than for computation.
] <ex-geo-bernstein-linear>

=== Computational Exercises

#exercise(title: [Runge Phenomenon: Equispaced vs Chebyshev])[
  Interpolate $f(x) = 1\/(1 + 25x^2)$ on $[-1, 1]$ using (a) equispaced nodes and (b) Chebyshev nodes, for $N = 8, 16, 24, 32$. (a) Plot the interpolants and the pointwise error $|f(x) - p_N (x)|$ on a fine evaluation grid. (b) Plot the maximum interpolation error versus $N$ on a semilogarithmic scale for both node distributions. (c) Explain the divergence for equispaced nodes in terms of the potential-theoretic results of @sec-potential-theory. The script `runge_phenomenon` is a useful starting point.
] <ex-geo-runge>

#exercise(title: [Lagrange Interpolation Error and Analyticity])[
  For the function $f(x) = e^x$ on $[-1, 1]$: (a) Compute the Chebyshev interpolant $p_N (x)$ for $N = 4, 8, 12, 16, 20$ and measure the maximum error $norm(f - p_N)_infinity$. (b) Plot the error versus $N$ on a semilogarithmic scale and verify exponential convergence. (c) The convergence rate is governed by the parameter $rho$ of the Bernstein ellipse inside which $f$ is analytic. For $f(x) = e^x$, this function is entire, so $rho arrow infinity$ and the convergence is super-exponential. Verify that the error decays faster than $C rho^(-N)$ for any fixed $rho > 1$.
] <ex-geo-analytic-convergence>

#exercise(title: [Convergence Rate and the Distance to the Pole])[
  Consider the family $f_alpha (x) = 1 \/ (alpha^2 + x^2)$ on $[-1, 1]$, whose poles lie at $z = plus.minus i alpha$. For $alpha in {1, 0.5, 0.2, 0.1, 0.05}$: (a) Compute the Chebyshev interpolant for $N$ up to a few hundred and measure $norm(f_alpha - p_N)_infinity$. (b) Fit the geometric decay $norm(f_alpha - p_N)_infinity tilde.op C rho^(-N)$ and compare the measured rate against the prediction $rho = alpha + sqrt(1 + alpha^2)$ implied by @eq-geometric-convergence. (c) Plot the measured and predicted $rho$ against $alpha$ and confirm that $rho arrow.r 1^+$ as the pole approaches the interval. The script `convergence_zoom` is a useful starting point.
] <ex-geo-pole-sweep>

#exercise(title: [Lebesgue Constant Computation])[
  For $N = 2, 4, 8, 16, 32, 64$, compute the Lebesgue constant $Lambda_N = max_(x in [-1,1]) sum_(j=0)^N |L_j (x)|$ for (a) equispaced nodes, (b) Chebyshev nodes of the first kind, and (c) Chebyshev nodes of the second kind (Chebyshev--Lobatto). Evaluate the maximum over a fine grid of 10000 points. Plot all three curves on a semilogarithmic scale and verify the asymptotic growth rates $Lambda_N^("equi") tilde.op 2^(N+1) \/ (e N ln N)$ and $Lambda_N^("Cheb") tilde.op (2\/pi) ln N$. The script `lebesgue_functions` is a useful starting point.
] <ex-geo-lebesgue-compute>

#exercise(title: [Barycentric versus Naive Lagrange Evaluation])[
  Implement two evaluators for the Chebyshev interpolant: the naive Lagrange sum @eq-lagrange-formula built from the basis @eq-lagrange-basis, and the barycentric formula @eq-barycentric with the Chebyshev weights @eq-chebyshev-weights. (a) For $f(x) = 1 \/ (1 + 16 x^2)$ and growing $N$, compare the maximum evaluation error of the two formulas and show that the naive form loses accuracy through cancellation while the barycentric form does not. (b) Count operations and time both evaluators on a fixed set of $M$ output points, confirming the $O(N)$ versus $O(N^2)$ cost per evaluation point. (c) Verify experimentally that rescaling the barycentric weights by a common constant leaves the result unchanged. The script `runge_phenomenon` provides a naive Lagrange routine to start from.
] <ex-geo-barycentric-stability>

#exercise(title: [Potential-Theoretic Equilibrium Measure])[
  The equilibrium measure for a compact set $K subset.eq RR$ minimises the logarithmic energy. For $K = [-1, 1]$, the equilibrium measure has the density $rho(x) = 1\/(pi sqrt(1 - x^2))$, the arcsine distribution. (a) Generate $N = 500$ random points from this distribution using the inverse CDF method $x_j = cos(pi U_j)$ where $U_j tilde.op "Uniform"(0, 1)$. Histogram the result and overlay the theoretical density. (b) Compute the Lebesgue constant of these random arcsine-distributed nodes and compare it to the Lebesgue constant of the deterministic Chebyshev nodes. (c) Explain why deterministic Chebyshev nodes vastly outperform random samples from the same density.
] <ex-geo-equilibrium-measure>

=== Project-Style Exercises

#exercise(title: [Optimal Interpolation Points by Optimization])[
  The nodes that minimise the Lebesgue constant $Lambda_N$ for a given polynomial degree are the Lebesgue points#idx("Lebesgue points"). (a) For $N$ up to about ten, set up the minimisation of $Lambda_N$ over the interior node positions, with the endpoints fixed at $plus.minus 1$, and solve it numerically using the equioscillation of the optimal Lebesgue function. (b) Compare the optimal constants against the Chebyshev values and against the lower bound of @Vertesi1990, confirming that the gap to Chebyshev is the small additive constant $(2 \/ pi) ln 2 approx 0.44$ reported in this chapter. (c) Discuss why, given this tiny gap, Chebyshev points are preferred in practice over the harder-won optimal points.
] <ex-geo-lebesgue-points>

#hint-for(<ex-geo-lebesgue-points>)[At the optimum the local maxima of the Lebesgue function between consecutive nodes are all equal (equioscillation); exploit these equal-ripple conditions, or a direct minimax optimiser, rather than a smooth gradient method, since $Lambda_N$ is only piecewise smooth in the node positions.]

#exercise(title: [Fekete and Padua Points in Two Dimensions])[
  Multivariate interpolation forces a choice of node geometry on the domain. (a) On the square $[-1, 1]^2$, implement the Padua points#idx("Padua points") of @DeMarchi2005 and verify numerically the $O(ln^2 N)$ growth of their Lebesgue constant established by @Bos2006 @Bos2007, contrasting them with tensor-product Chebyshev grids. (b) On the triangle, compute Fekete points#idx("Fekete points") by maximising the Vandermonde determinant for low to moderate degree and tabulate the resulting Lebesgue constants, as used in spectral-element practice @Hesthaven2007. (c) Discuss what remains open about minimising the Lebesgue constant on the triangle.
] <ex-geo-multidim-points>

#hint-for(<ex-geo-multidim-points>)[The Padua points are the boundary contacts and self-intersections of a single Lissajous curve, so they can be written in closed form before any optimisation; the Fekete computation, by contrast, is a genuine maximisation of $|det V|$ over node positions and is best warm-started from a structured grid.]

#exercise(title: [Potential Theory for General Node Families])[
  Use the logarithmic potential framework of @SaffTotik1997 to predict convergence for arbitrary node densities#idx("equilibrium measure"). (a) For a one-parameter family of densities interpolating between the uniform density and the Chebyshev density @eq-chebyshev-density, evaluate the potential @eq-potential numerically and plot its equipotential curves in the complex plane. (b) For a test function with known complex singularities, predict the region of convergence as the largest singularity-free equipotential and compare it against measured interpolation errors. (c) Relate the limiting case to the Chebyshev potential @eq-chebyshev-potential and the Bernstein-ellipse rate @eq-geometric-convergence.
] <ex-geo-potential-general>

#hint-for(<ex-geo-potential-general>)[Convergence at a point $z$ is governed by whether the potential there exceeds its value on the interval $[-1, 1]$; the interpolant converges inside the largest equipotential level curve that encloses no singularity of $f$, so the decisive quantity is the potential difference between the nearest singularity and the interval.]

#exercise(title: [Minimum Separation and Random Node Families])[
  The two computational études of this chapter, in Sections #ref(<sec-random-nodes>, supplement: none) and #ref(<sec-random-chebyshev>, supplement: none), traced the failure of random nodes to the absence of a guaranteed minimum separation. (a) Construct node families with tunable separation, for example jittered Chebyshev points, van der Corput or other low-discrepancy sequences, and Chebyshev angles perturbed by noise of controlled amplitude. (b) Measure how the Lebesgue constant depends on the smallest node gap and seek an empirical law linking the two. (c) Relate the findings to the observation of @Fornberg2025 that mild clustering already suffices for high-order accuracy, and propose a practical criterion for deciding when a near-Chebyshev grid is good enough.
] <ex-geo-min-separation>

#hint-for(<ex-geo-min-separation>)[Hold the marginal density fixed at the arcsine law and vary only the minimum gap, so that any change in the Lebesgue constant is attributable to separation rather than to density; plotting $ln Lambda_N$ against the reciprocal of the smallest gap is a good first diagnostic.]
