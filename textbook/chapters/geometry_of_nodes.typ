// textbook/chapters/geometry_of_nodes.typ
#import "../styles/template.typ": dropcap

= The Geometry of Nodes

#dropcap[Before we can differentiate functions numerically using spectral methods, we must first understand how to represent them. Polynomial interpolation---the process of constructing a polynomial that passes through a given set of data points---is the foundation upon which pseudospectral methods are built. In this chapter, we explore a fascinating paradox: while polynomial interpolation seems entirely straightforward, the choice of interpolation nodes determines whether the method succeeds brilliantly or fails catastrophically.]

The story begins with a surprising discovery by the German mathematician Carl Runge in 1901. Attempting to approximate a simple, smooth function by interpolating polynomials, Runge found that increasing the polynomial degree made the approximation _worse_, not better. This counterintuitive phenomenon, now bearing his name, reveals deep connections between numerical analysis, complex analysis, and potential theory.

Our journey through this "étude in grid geometry" will explain the Runge phenomenon, develop the theoretical framework of potential theory that predicts where interpolation succeeds or fails, and introduce the Chebyshev points that form the cornerstone of practical spectral methods.

== The Problem: Polynomial Interpolation

=== Lagrange Interpolation

Given $N + 1$ distinct points ${x_0, x_1, dots, x_N}$ on an interval $[a, b]$ and corresponding function values ${f_0, f_1, dots, f_N}$ where $f_j = f(x_j)$, there exists a unique polynomial $p_N (x)$ of degree at most $N$ such that
$ p_N (x_j) = f_j, quad j = 0, 1, dots, N. $ <eq-interpolation-condition>

This interpolating polynomial can be written explicitly using the _Lagrange formula_:
$ p_N (x) = sum_(k=0)^N f_k L_k (x), $ <eq-lagrange-formula>
where the _Lagrange basis polynomials_ are
$ L_k (x) = product_(j = 0, j eq.not k)^N frac(x - x_j, x_k - x_j). $ <eq-lagrange-basis>

Each basis polynomial $L_(k) (x)$ has the cardinal property: it equals $1$ at $x_k$ and $0$ at all other nodes. This property ensures that substituting any node $x_j$ into the Lagrange formula yields $f_j$.

=== Interpolation versus Best Approximation

It is important to distinguish interpolation from best approximation. The Weierstrass Approximation Theorem guarantees that any continuous function on $[a, b]$ can be uniformly approximated to arbitrary accuracy by polynomials. A constructive proof was given by Sergei Bernstein in 1912 using what are now called _Bernstein polynomials_---explicit polynomial approximations that converge uniformly to any continuous function, though typically more slowly than interpolation-based methods. However, this theorem says nothing about _which_ polynomials achieve this approximation or how to construct them efficiently.

Interpolation constructs a specific polynomial by enforcing exact agreement at the nodes. The hope is that as $N arrow infinity$, the interpolating polynomials $p_N$ converge uniformly to $f$. As we shall see, this hope is fulfilled for some node distributions but dramatically fails for others.

=== Equispaced Nodes

The most natural choice of interpolation nodes is equally spaced points:
$ x_j = -1 + frac(2 j, N), quad j = 0, 1, dots, N. $ <eq-equispaced-nodes>

These nodes divide the interval $[-1, 1]$ into $N$ equal subintervals. For low-degree polynomials and well-behaved functions, equispaced interpolation works reasonably well. However, a fundamental problem emerges as $N$ increases.

== The Runge Phenomenon

=== A Smooth but Troublesome Function

In 1901, Carl Runge studied the interpolation of a deceptively simple function:
$ f(x) = frac(1, 1 + 25 x^2). $ <eq-runge-function>

This _Runge function_ is infinitely differentiable on the entire real line. Its graph is a smooth bell curve centered at the origin with maximum value $f(0) = 1$ and asymptotic decay to zero as $|x| arrow infinity$.

Runge discovered that polynomial interpolation on equispaced nodes _diverges_ for this function. Rather than improving with increasing $N$, the interpolating polynomials develop wild oscillations near the endpoints $x = plus.minus 1$.

=== Numerical Demonstration

@fig-runge-phenomenon illustrates this phenomenon. For $N = 6$, the interpolant reasonably approximates the function. But for $N = 10$ and $N = 14$, the approximation deteriorates dramatically at the boundaries, with oscillations that grow unboundedly as $N$ increases.

#figure(
  image("../figures/ch04/python/runge_phenomenon.pdf", width: 85%),
  caption: [The Runge phenomenon: polynomial interpolation of $f(x) = 1\/(1 + 25x^2)$ on equispaced nodes. As the polynomial degree increases, oscillations near the boundaries $x = plus.minus 1$ grow unboundedly. The exact function (solid) and interpolants (dashed) are shown for $N = 6$, $10$, and $14$.],
) <fig-runge-phenomenon>

The following Python code computes the Lagrange interpolant for a given set of nodes:

```python
def lagrange_interpolate(x_nodes, f_nodes, x_eval):
    """Evaluate the Lagrange interpolating polynomial."""
    n = len(x_nodes)
    p_eval = np.zeros_like(x_eval)

    for k in range(n):
        L_k = np.ones_like(x_eval)
        for j in range(n):
            if j != k:
                L_k *= (x_eval - x_nodes[j]) / (x_nodes[k] - x_nodes[j])
        p_eval += f_nodes[k] * L_k

    return p_eval
```

The equivalent MATLAB implementation:

```matlab
function p = lagrange_interp(x_nodes, f_nodes, x_eval)
    n = length(x_nodes);
    p = zeros(size(x_eval));

    for k = 1:n
        L_k = ones(size(x_eval));
        for j = 1:n
            if j ~= k
                L_k = L_k .* (x_eval - x_nodes(j)) / (x_nodes(k) - x_nodes(j));
            end
        end
        p = p + f_nodes(k) * L_k;
    end
end
```

The code generating @fig-runge-phenomenon is available in:
- `codes/python/ch04_geometry_of_nodes/runge_phenomenon.py`
- `codes/matlab/ch04_geometry_of_nodes/runge_phenomenon.m`

=== Why Does This Happen?

The Runge phenomenon seems paradoxical: how can adding more information (more interpolation points) make the approximation worse? The answer lies in the _Lebesgue constant_, which we will explore in @sec-lebesgue. But first, let us develop a deeper understanding through potential theory.

== Theoretical Explanation: Potential Theory <sec-potential-theory>

=== Singularities in the Complex Plane

The key to understanding the Runge phenomenon lies in extending our view from the real line to the complex plane. The Runge function has singularities (poles) where its denominator vanishes:
$ 1 + 25 z^2 = 0 quad arrow.r.double quad z = plus.minus frac(i, 5) = plus.minus 0.2 i. $

These poles lie on the imaginary axis, at a distance of only $0.2$ from the real interval $[-1, 1]$. This proximity is responsible for the divergence of equispaced interpolation.

=== The Potential Function

The convergence of polynomial interpolation is governed by a _potential function_ associated with the node distribution. For a distribution with density $mu(x)$ on $[-1, 1]$, the potential at a point $z$ in the complex plane is:
$ phi(z) = - integral_(-1)^1 mu(x) ln |z - x| dif x. $ <eq-potential>

The _equipotential curves_ $phi(z) = "const"$ form closed contours around the interval $[-1, 1]$. The largest equipotential curve that does not enclose any singularity of $f$ determines the region of convergence.

=== Uniform versus Chebyshev Density

For equispaced nodes, the density is asymptotically uniform: $mu(x) = 1\/2$. The resulting equipotential curves are roughly elliptical, but they extend only a short distance from the interval before reaching the Runge function's poles at $plus.minus 0.2 i$.

For Chebyshev nodes (which we introduce in the next section), the density is:
$ mu(x) = frac(1, pi sqrt(1 - x^2)). $ <eq-chebyshev-density>

This density diverges at the endpoints, concentrating nodes near $x = plus.minus 1$. The corresponding potential simplifies to:
$ phi(z) = ln |z + sqrt(z^2 - 1)| - ln 2 = ln rho - ln 2, $ <eq-chebyshev-potential>
where $rho = |z + sqrt(z^2 - 1)|$. The equipotential curves are _Bernstein ellipses_ with foci at $plus.minus 1$.

@fig-equipotential-curves compares the equipotential curves for both distributions, showing why Chebyshev interpolation succeeds where equispaced interpolation fails.

#figure(
  image("../figures/ch04/python/equipotential_curves.pdf", width: 95%),
  caption: [Equipotential curves in the complex plane. Left: uniform node density (equispaced nodes). Right: Chebyshev density, where equipotentials are Bernstein ellipses. The poles of the Runge function at $z = plus.minus 0.2 i$ are marked with crosses. For Chebyshev interpolation, the critical ellipse passing through the poles determines the convergence rate.],
) <fig-equipotential-curves>

The code generating @fig-equipotential-curves is available in:
- `codes/python/ch04_geometry_of_nodes/equipotential_curves.py`
- `codes/matlab/ch04_geometry_of_nodes/equipotential_curves.m`

== The Solution: Chebyshev Points <sec-chebyshev-points>

=== Definition and Geometric Construction

The _Chebyshev-Gauss-Lobatto points_ (often simply called Chebyshev points) are defined by:
$ x_j = cos(frac(j pi, N)), quad j = 0, 1, dots, N. $ <eq-chebyshev-points>

These points have a beautiful geometric interpretation: they are the projections onto the $x$-axis of $N + 1$ equally spaced points on the upper semicircle of the unit circle. @fig-chebyshev-circle illustrates this construction.

#figure(
  image("../figures/ch04/python/chebyshev_points_circle.pdf", width: 80%),
  caption: [Geometric construction of Chebyshev points. Points equally spaced on the upper semicircle (by angle $theta_j = j pi \/ N$) project vertically onto the Chebyshev nodes on the $x$-axis. Equal spacing on the circle maps to clustering near the endpoints $x = plus.minus 1$.],
) <fig-chebyshev-circle>

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

The Chebyshev nodes in Python and MATLAB:

```python
def chebyshev_nodes(N):
    """Chebyshev-Gauss-Lobatto points on [-1, 1]."""
    j = np.arange(N + 1)
    return np.cos(j * np.pi / N)
```

```matlab
% Chebyshev-Gauss-Lobatto points
j = 0:N;
x_cheb = cos(j * pi / N);
```

Note that this formula produces nodes ordered from right to left: $x_0 = cos(0) = +1$ down to $x_N = cos(pi) = -1$. This ordering is natural for the cosine function and is the standard convention in spectral methods.

The code generating @fig-chebyshev-success is available in:
- `codes/python/ch04_geometry_of_nodes/chebyshev_success.py`
- `codes/matlab/ch04_geometry_of_nodes/chebyshev_success.m`

== Lagrange Basis Functions and Lebesgue Constants <sec-lebesgue>

=== The Interpolation Operator

To understand why node choice matters so dramatically, we formalize the interpolation process. For any continuous function $u in C([-1,1])$, there exists a unique interpolating polynomial $cal(P)_N (x) in bb(P)_N [bb(R)]$ of degree $N$ satisfying:
$ cal(P)_N (x_j) = u_j equiv u(x_j), quad j = 0, 1, dots, N. $
We denote this interpolating polynomial as $cal(I)_N [u]$, emphasizing its role as a _linear operator_ acting on continuous functions.

The _best approximating polynomial_ $cal(P)^* in bb(P)_N [bb(R)]$ is defined by:
$ norm(u - cal(P)^*)_infinity = inf_(cal(P) in bb(P)_N [bb(R)]) norm(u - cal(P))_infinity. $
While it is not obliged that $cal(P)^* equiv cal(I)_N [u]$, since $cal(P)^* in bb(P)_N [bb(R)]$, we necessarily have $cal(P)^* equiv cal(I)_N [u]$ when interpolating a polynomial of degree at most $N$.

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

- *Chebyshev nodes*: $Lambda_N^"Ch" = frac(2, pi) ln N + O(1)$ (logarithmic growth)
- *Legendre nodes*: $Lambda_N^"Leg" = O(sqrt(N))$ (slow algebraic growth)
- *Equispaced nodes*: $Lambda_N^"eq" approx frac(2^(N+1), e N ln N)$ (exponential growth!)

The exponential growth of the Lebesgue constant for equispaced nodes explains the Runge phenomenon: even though the best approximation error decreases geometrically for smooth functions, the exponentially growing factor $(1 + Lambda_N)$ eventually dominates.

*Remark.* As shown by Vértesi @Vertesi1990, the Lebesgue constant for Chebyshev nodes is remarkably close to the smallest possible Lebesgue constant for _any_ node distribution:
$ Lambda_N^"Ch" = frac(2, pi) (ln N + gamma + ln frac(8, pi)) + o(1), $
$ Lambda_N^"min" = frac(2, pi) (ln N + gamma + ln frac(4, pi)) + o(1), $
where $gamma approx 0.5772156649...$ is the Euler--Mascheroni constant. The difference between these two expressions is only $frac(2, pi) ln 2 approx 0.44$, demonstrating that Chebyshev nodes are essentially optimal for polynomial interpolation.

=== Lebesgue Points and Multi-Dimensional Considerations

The nodes that _minimize_ the Lebesgue constant $Lambda_N arrow min$ for a given polynomial space $bb(P)_N [bb(R)]$ are called _Lebesgue points_. Finding these optimal nodes is a challenging optimization problem that has been solved only for moderate values of $N$ in one dimension.

In multiple dimensions, the situation becomes considerably more complex. The Lebesgue constant depends not only on the node distribution ${bold(x)_j}_(j=0)^N$ but also on the _shape of the domain_ $cal(U)$. The estimation of the Lebesgue constant in multi-dimensional non-Cartesian domains is a problem that remains essentially open today. The most precious information would be to find the node distribution over a triangle that minimizes the Lebesgue constant---knowledge that would be crucial for the design of new spectral elements @Hesthaven2007.

On _quadrilateral_ domains, the best known point distributions are tensor products of Gauss--Lobatto#footnote[Rehuel Lobatto (1797--1866) was a Dutch mathematician born into a Portuguese family.] or Chebyshev nodes. These tensor product constructions inherit the favorable properties of their one-dimensional counterparts.

On _triangular_ domains, the situation is more subtle. The _Fekete points_#footnote[Michael Fekete (1886--1957) was a Hungarian mathematician who completed his PhD under the supervision of Lipót Fejér. Fekete also provided private tutoring to János Neumann, known today as John von Neumann.] currently represent the best known choice. Fekete points are defined as the nodes that maximize the determinant of the Vandermonde matrix---equivalently, they maximize the volume of the interpolation simplex in function space. While not provably optimal in the Lebesgue sense, they exhibit excellent numerical properties and remain the standard choice for triangular spectral elements.

=== Visualization of Basis Functions

@fig-lagrange-basis compares the Lagrange basis polynomials for equispaced and Chebyshev nodes. For equispaced nodes, the basis functions near the endpoints develop large oscillations. For Chebyshev nodes, all basis functions remain well-behaved.

#figure(
  image("../figures/ch04/python/lagrange_basis.pdf", width: 95%),
  caption: [Lagrange basis polynomials $L_k(x)$ for $N = 10$. Left: equispaced nodes show large oscillations near the boundaries. Right: Chebyshev nodes yield bounded basis functions throughout the interval.],
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
- `codes/python/ch04_geometry_of_nodes/lagrange_basis.py`
- `codes/python/ch04_geometry_of_nodes/lebesgue_functions.py`
- `codes/python/ch04_geometry_of_nodes/lebesgue_constants_zoom.py`
- `codes/matlab/ch04_geometry_of_nodes/lagrange_basis.m`
- `codes/matlab/ch04_geometry_of_nodes/lebesgue_functions.m`

== Barycentric Interpolation

=== Numerical Stability Issues

The Lagrange formula, while mathematically elegant, is numerically problematic. Computing each basis polynomial $L_k (x)$ requires $N$ multiplications and divisions, and the resulting values can be very large or very small, leading to catastrophic cancellation when summed.

=== The Barycentric Formula

A more stable and efficient approach is the _barycentric interpolation formula_:
$ p_N (x) = frac(sum_(j=0)^N frac(w_j, x - x_j) f_j, sum_(j=0)^N frac(w_j, x - x_j)), $ <eq-barycentric>
where the _barycentric weights_ are:
$ w_j = frac(1, product_(k eq.not j) (x_j - x_k)). $ <eq-barycentric-weights>

This formula requires only $O(N)$ operations to evaluate $p_N(x)$ once the weights are precomputed.

=== Weights for Chebyshev Points

For Chebyshev-Gauss-Lobatto points, the weights have a particularly simple form:
$ w_j = (-1)^j delta_j, quad "where" delta_j = cases(1\/2 & "if" j = 0 "or" j = N, 1 & "otherwise".) $ <eq-chebyshev-weights>

This simplicity, combined with the excellent approximation properties, makes Chebyshev points the standard choice for spectral methods.

== Convergence Analysis <sec-convergence>

=== Geometric Convergence

For functions analytic in a region containing $[-1, 1]$, Chebyshev interpolation converges geometrically. If $f$ is analytic inside the Bernstein ellipse $cal(E)_rho$ with parameter $rho > 1$, then:
$ norm(f - p_N)_infinity lt.eq.slant frac(C, rho^N), $ <eq-geometric-convergence>
where $C$ depends on $f$ and $rho$ but not on $N$.

The parameter $rho$ is determined by the location of the nearest singularity: if the closest singularity to $[-1, 1]$ lies on the ellipse $cal(E)_rho$, then convergence is at rate $O(rho^(-N))$.

=== The Runge Function Revisited

For the Runge function with poles at $z = plus.minus 0.2i$, we can compute:
$ rho = |0.2i + sqrt((0.2i)^2 - 1)| = |0.2i + i sqrt(1.04)| approx 1.22. $

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
- `codes/python/ch04_geometry_of_nodes/convergence_comparison.py`
- `codes/python/ch04_geometry_of_nodes/convergence_zoom.py`
- `codes/matlab/ch04_geometry_of_nodes/convergence_comparison.m`

== Computational Experiment: Random Nodes <sec-random-nodes>

Having studied the optimal Chebyshev distribution and the problematic equispaced distribution, a natural question arises: what happens if we choose interpolation nodes _randomly_? This question leads us into the realm of _experimental mathematics_, where we use computation to discover and conjecture mathematical relationships.

=== Motivation

The dramatically different behaviors of equispaced and Chebyshev nodes might suggest that the key factor is _regularity_ or _structure_ in node placement. Perhaps any "reasonable" arrangement would work? To test this hypothesis, we investigate the simplest unstructured choice: nodes drawn uniformly at random from $[-1, 1]$.

This experiment serves two purposes. First, it tests whether the special clustering of Chebyshev points near the endpoints is truly essential, or whether it is merely one of many acceptable distributions. Second, it demonstrates the methodology of computational experimentation---using numerical evidence to formulate conjectures about asymptotic behavior.

=== Experimental Setup

For each polynomial degree $N$ from $2$ to $30$, we generate $N + 1$ random points uniformly distributed on $[-1, 1]$, sort them, and compute the Lebesgue constant. We repeat this process $M = 200$ times to obtain statistical estimates of the mean, standard deviation, and range of $Lambda_N$ for random nodes.

The Python implementation uses NumPy's random number generation:

```python
def random_nodes(N, rng=None):
    """Generate N+1 random nodes, sorted, on [-1, 1]."""
    if rng is None:
        rng = np.random.default_rng()
    return np.sort(rng.uniform(-1, 1, N + 1))

def monte_carlo_lebesgue(N, M=200):
    """Compute M samples of Lebesgue constant for random nodes."""
    samples = np.zeros(M)
    for m in range(M):
        x_rand = random_nodes(N)
        samples[m] = lebesgue_constant(x_rand)
    return samples
```

The equivalent MATLAB code:

```matlab
% Generate M samples for polynomial degree N
samples = zeros(M, 1);
for m = 1:M
    x_rand = sort(2 * rand(N+1, 1) - 1);  % Uniform on [-1, 1]
    samples(m) = max(lebesgue_function(x_rand, x_fine));
end
```

=== Results

@fig-lebesgue-random shows the results of our Monte Carlo experiment. The left panel displays the growth of the Lebesgue constant with polynomial degree, including shaded regions showing the statistical variability. The right panel shows the distribution of $Lambda_N$ values for a fixed degree $N = 15$.

#figure(
  image("../figures/ch04/python/lebesgue_random_nodes.pdf", width: 95%),
  caption: [Lebesgue constant for random nodes. Left: growth with polynomial degree $N$, showing mean (solid line) and min-max range (light shading). Right: distribution of $log_(10)(Lambda_(15))$ over $500$ random trials on a logarithmic scale, revealing the enormous spread spanning many orders of magnitude. Vertical lines mark the mean, median, equispaced, and Chebyshev values.],
) <fig-lebesgue-random>

The key observations are striking:

+ *Explosive exponential growth*: Random nodes exhibit exponential growth of the Lebesgue constant, but surprisingly _faster_ than equispaced nodes. The mean $Lambda_N$ grows roughly as $e^(b N)$ with $b approx 1.0$, compared to $b approx 0.6$ for equispaced.

+ *Extreme variability*: There is enormous variability between different random samples. For $N = 15$, $Lambda_N$ can range from around $10^2$ to $10^{10}$ or more, spanning many orders of magnitude. This variability increases dramatically with $N$.

+ *Worse than equispaced*: Counter to naive intuition, random nodes perform _worse_ on average than the structured equispaced distribution. The reason is that random sampling occasionally produces clusters of nearby nodes, creating near-singular interpolation conditions. These worst-case realizations dominate the mean.

=== Empirical Asymptotic Formula

Fitting an exponential model to the mean Lebesgue constants yields an empirical formula:
$ Lambda_N^"rand" approx a dot e^(b N), $ <eq-random-asymptotic>
where the fitted constants are approximately $a approx 0.7$ and $b approx 1.0$, giving:
$ Lambda_N^"rand" approx 0.7 dot 2.7^N. $

Remarkably, the growth rate for random nodes is _faster_ than for equispaced nodes (which grow as $approx 1.8^N$). This counterintuitive result arises because random sampling occasionally produces clusters of nodes that are extremely close together, causing the Lebesgue constant to explode. The mean is dominated by these worst-case realizations.

The fundamental problem is twofold: (1) the lack of _clustering near the endpoints_, which only carefully designed distributions like Chebyshev points possess, and (2) the risk of _random clustering_ in the interior, which creates severely ill-conditioned interpolation problems.

=== Discussion

This computational experiment provides strong numerical evidence for a fundamental principle: *the clustering of Chebyshev points near the endpoints is not merely convenient but essential*. Random nodes, despite their apparent "fairness" in covering the interval, fail in two ways: they do not provide the boundary resolution needed to control Lagrange basis oscillations, and they risk interior clustering that creates catastrophically ill-conditioned systems.

The enormous variability in $Lambda_N$ for random nodes is perhaps the most striking finding. While equispaced nodes are suboptimal, they at least provide _predictable_ (if exponentially growing) behavior. Random nodes introduce an additional layer of uncertainty---any given random realization might be acceptable or catastrophic.

The experiment illustrates the power of computational mathematics. By systematic numerical investigation, we have:
- Tested a natural hypothesis (random nodes might work)
- Discovered unexpected behavior (random is _worse_ than equispaced)
- Quantified the asymptotic growth through data fitting
- Reinforced our understanding of why structured distributions matter

The code generating @fig-lebesgue-random is available in:
- `codes/python/ch04_geometry_of_nodes/lebesgue_random_nodes.py`
- `codes/matlab/ch04_geometry_of_nodes/lebesgue_random_nodes.m`

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
