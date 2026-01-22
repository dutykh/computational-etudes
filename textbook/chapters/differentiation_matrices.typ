// textbook/chapters/differentiation_matrices.typ
// Chapter 5: Differentiation Matrices
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Last modified: January 2026

#import "../styles/template.typ": dropcap

// Enable equation numbering for this chapter
#set math.equation(numbering: "(1)")

= Differentiation Matrices

#dropcap[In the previous chapter, we mastered the art of polynomial interpolation---constructing polynomials that pass exactly through a set of data points. We discovered that the choice of nodes determines whether interpolation succeeds or fails, with Chebyshev points emerging as the optimal choice for non-periodic problems. Now we take the next logical step: having represented a function as an interpolating polynomial, how do we _differentiate_ it? The answer leads us to one of the most elegant structures in numerical analysis: the differentiation matrix.]

The remarkable insight of pseudospectral methods is that differentiation can be accomplished by a single matrix-vector multiplication. Given function values $bold(u) = (u_0, u_1, dots, u_N)^top$ at the grid points, we can approximate the derivative values $bold(u)' = (u'_0, u'_1, dots, u'_N)^top$ as
$ bold(u)' approx D bold(u), $
where $D$ is the _differentiation matrix_. This matrix encapsulates the entire differentiation process: interpolate, differentiate, evaluate.

This chapter develops the theory and practice of constructing differentiation matrices. We begin with familiar finite difference approximations, viewing them through the lens of matrix algebra. We then take the crucial limiting step: what happens when the stencil extends to include _all_ grid points? The answer reveals that spectral methods are not a separate species from finite differences, but rather their natural culmination---the limiting case as the stencil width grows without bound.

== From Interpolation to Differentiation

=== The Matrix Perspective

Let us begin with a fundamental observation. Given $N + 1$ distinct nodes ${x_0, x_1, dots, x_N}$ and function values ${u_0, u_1, dots, u_N}$, the unique interpolating polynomial of degree at most $N$ can be written using the Lagrange formula:
$ p_N (x) = sum_(j=0)^N u_j L_j (x), $
where $L_j (x)$ are the Lagrange basis polynomials from @sec-lebesgue.

To approximate the derivative of $u$ at the nodes, we simply differentiate the interpolating polynomial:
$ p'_N (x) = sum_(j=0)^N u_j L'_j (x). $

Evaluating this at node $x_i$ yields:
$ p'_N (x_i) = sum_(j=0)^N u_j L'_j (x_i). $ <eq-deriv-at-node>

This is a linear combination of the function values! The coefficients $L'_j (x_i)$ depend only on the node locations, not on the function values. We can therefore write @eq-deriv-at-node in matrix form:
$ bold(u)' approx D bold(u), quad "where" quad D_(i j) = L'_j (x_i). $ <eq-diff-matrix-def>

The entry $D_(i j)$ represents the weight given to the function value at $x_j$ when approximating the derivative at $x_i$.

=== The Barycentric Formula for Differentiation Matrix Entries

For general (non-equispaced) nodes, the differentiation matrix entries can be computed using the barycentric weights $w_j$ introduced in @sec-lebesgue:
$ D_(i j) = frac(w_j \/ w_i, x_i - x_j) quad "for" quad i eq.not j. $ <eq-diff-entry-bary>

The diagonal entries are determined by the important _consistency condition_: the derivative of a constant function is zero. Since $D$ applied to the constant vector $(1, 1, dots, 1)^top$ must yield zero, we have
$ sum_(j=0)^N D_(i j) = 0 quad arrow.r.double quad D_(i i) = - sum_(j eq.not i) D_(i j). $ <eq-diff-diagonal>

This row-sum property provides a convenient way to compute the diagonal entries and serves as a useful sanity check.

== Finite Difference Matrices <sec-fd-matrices>

=== The Local Approach

Before tackling the full spectral differentiation matrix, let us review the familiar territory of finite differences. The classical approach approximates the derivative using only _nearby_ function values---a local stencil.

The simplest example is the _second-order central difference_:
$ u'(x_i) approx frac(u_(i+1) - u_(i-1), 2h), $ <eq-fd2>
where $h$ is the grid spacing. This formula is "second-order accurate" because the error is $O(h^2)$ for smooth functions, as can be verified by Taylor expansion.

For higher accuracy, we can include more neighbors. The _fourth-order central difference_ uses five points:
$ u'(x_i) approx frac(-u_(i+2) + 8 u_(i+1) - 8 u_(i-1) + u_(i-2), 12h). $ <eq-fd4>

The _sixth-order central difference_ extends to seven points:
$ u'(x_i) approx frac(u_(i+3) - 9 u_(i+2) + 45 u_(i+1) - 45 u_(i-1) + 9 u_(i-2) - u_(i-3), 60h). $ <eq-fd6>

=== Matrix View of Finite Differences

These formulas can be assembled into differentiation matrices. For the second-order scheme @eq-fd2, assuming periodic boundary conditions on $N$ equispaced points with spacing $h = 2 pi \/ N$, the matrix is _tridiagonal_ (plus corner entries for periodicity):
$ D^((2)) = frac(1, 2h) mat(
  0, 1, 0, dots.c, 0, -1;
  -1, 0, 1, dots.c, 0, 0;
  0, -1, 0, dots.c, 0, 0;
  dots.v, dots.down, dots.down, dots.down, dots.down, dots.v;
  1, 0, 0, dots.c, -1, 0
). $ <eq-fd2-matrix>

The fourth-order scheme @eq-fd4 produces a _pentadiagonal_ matrix with bandwidth 5, and the sixth-order scheme yields a _heptadiagonal_ matrix with bandwidth 7.

@fig-fd-bandwidth illustrates this progression. As the order of accuracy increases, the stencil widens and the matrix bandwidth grows.

#figure(
  image("../figures/ch05/python/fd_matrix_bandwidth.pdf", width: 95%),
  caption: [Sparsity patterns of differentiation matrices for $N = 20$ grid points. Left: second-order finite differences (tridiagonal, bandwidth 3). Center: fourth-order finite differences (pentadiagonal, bandwidth 5). Right: spectral method (dense, bandwidth $N$). The progression illustrates the key insight: spectral methods are the limiting case of finite differences as the stencil width extends to the full domain.],
) <fig-fd-bandwidth>

=== The Limit Question

This observation leads to a fundamental question: _if widening the stencil improves accuracy, what is the ultimate limit?_

The answer is a stencil that uses _all_ $N$ grid points. The resulting matrix is _dense_---every point influences the derivative at every other point. This is precisely the spectral differentiation matrix.

From this perspective, spectral methods are not a separate species from finite differences. They are the natural endpoint of the progression: as we increase the order of accuracy by including more neighbors, we eventually include all points. The sparse banded matrix becomes dense, and algebraic convergence transforms into spectral (exponential) convergence.

This unifying viewpoint, emphasized by Fornberg @Fornberg1996, provides both conceptual clarity and practical guidance. It suggests that finite differences and spectral methods lie on a continuum, with the choice of stencil width being a tunable parameter balancing accuracy against computational cost.

== The Periodic Spectral Differentiation Matrix <sec-spectral-periodic>

=== Periodic Problems and Equispaced Nodes

A remarkable simplification occurs for _periodic_ problems. On a domain like $[0, 2 pi)$ with periodic boundary conditions, the natural interpolation basis consists of trigonometric polynomials (sines and cosines) rather than algebraic polynomials.

For periodic problems, the optimal node distribution is _equispaced_:
$ x_j = frac(2 pi j, N), quad j = 0, 1, dots, N - 1. $ <eq-periodic-nodes>

This is in stark contrast to the non-periodic case studied in @sec-chebyshev-points, where equispaced nodes lead to the Runge phenomenon. The difference lies in the boundary: periodic functions have no endpoints, so there is no need for the clustering that Chebyshev nodes provide. The periodicity itself regularizes the problem.

=== Derivation via the Fourier Basis

The spectral differentiation matrix for periodic problems can be derived by considering the trigonometric interpolant and differentiating it analytically. For $N$ equispaced points (with $N$ even), the interpolating trigonometric polynomial is
$ p(x) = sum_(j=0)^(N-1) u_j phi_j (x), $
where $phi_j (x)$ is the _periodic cardinal function_ (or discrete Dirichlet kernel) satisfying $phi_j (x_k) = delta_(j k)$.

The periodic cardinal function can be written as a sum of complex exponentials:
$ phi_j (x) = frac(1, N) sum_(k=-N\/2+1)^(N\/2) e^(i k (x - x_j)). $
This sum can be evaluated in closed form. Writing $theta = (x - x_j)\/2$ and using the geometric series, we obtain:
$ phi_j (x) = frac(sin(N theta \/ 2), N sin(theta \/ 2)) = frac(sin(N(x - x_j)\/2), N sin((x - x_j)\/2)). $ <eq-cardinal-periodic>

To find the differentiation matrix, we need to compute $phi'_j (x_m)$ for all $m$. Let $xi = (x - x_j)\/2$ for brevity. Then @eq-cardinal-periodic becomes $phi_j = sin(N xi) \/ (N sin xi)$. Applying the quotient rule:
$ frac(d phi_j, d xi) = frac(N cos(N xi) sin xi - sin(N xi) cos xi, N sin^2 xi). $
Since $d xi \/ d x = 1\/2$, we obtain:
$ phi'_j (x) = frac(1, 2) dot frac(cos(N xi), sin xi) - frac(1, 2) dot frac(sin(N xi) cos xi, N sin^2 xi). $ <eq-cardinal-deriv>
We evaluate this at the grid points $x = x_m$.

#block(inset: (left: 1em))[
  *Case 1: Diagonal entries* ($m = j$). \
  When $x = x_j$, we have $xi = 0$, so both numerator and denominator of @eq-cardinal-deriv vanish. Applying L'Hôpital's rule (twice) yields $phi'_j (x_j) = 0$. Geometrically, the cardinal function is symmetric about its peak at $x_j$, so its derivative must vanish there.
]

#block(inset: (left: 1em))[
  *Case 2: Off-diagonal entries* ($m eq.not j$). \
  When $x = x_m$ with $m eq.not j$, we have $xi = pi(m-j)\/N$. At these points:
  - $sin(N xi) = sin(pi(m-j)) = 0$ (the second term in the numerator vanishes),
  - $cos(N xi) = cos(pi(m-j)) = (-1)^(m-j)$.

  Substituting into @eq-cardinal-deriv:
  $ phi'_j (x_m) = frac(1, 2) frac((-1)^(m-j), tan(pi(m-j)\/N)) = frac((-1)^(m-j), 2) cot frac((m-j) pi, N). $
]

Since $D_(m k) = phi'_k (x_m)$, the spectral differentiation matrix entries are:
$ D_(j k) = cases(
  display(frac(1, 2) (-1)^(j-k) cot(frac((j - k) pi, N))) & "if" j eq.not k,
  0 & "if" j = k.
) $ <eq-spectral-periodic>

=== Properties of the Spectral Matrix

The matrix defined by @eq-spectral-periodic has several remarkable properties:

+ *Skew-symmetry*: $D^top = -D$. This reflects the fact that differentiation is an antisymmetric operator in appropriate function spaces.

+ *Zero diagonal*: $D_(j j) = 0$. The derivative approximation at any point depends only on neighboring values, not on the value at the point itself.

+ *Toeplitz structure*: $D_(j k)$ depends only on the difference $j - k$. This means the matrix has constant diagonals---a consequence of the translation-invariance of differentiation on periodic domains.

+ *Circulant*: Due to the periodic boundary conditions, the matrix is actually _circulant_: each row is a cyclic shift of the previous row. Circulant matrices can be diagonalized by the discrete Fourier transform, enabling $O(N log N)$ matrix-vector products via the FFT.

+ *Dense*: Unlike finite difference matrices, every off-diagonal entry is nonzero. This is the price we pay for spectral accuracy.

+ *Exactness for trigonometric polynomials*: If $u(x)$ is a trigonometric polynomial of degree at most $N\/2$, then $D bold(u)$ gives the _exact_ derivative values.

@fig-spectral-structure visualizes these properties for $N = 16$.

#figure(
  image("../figures/ch05/python/spectral_matrix_structure.pdf", width: 95%),
  caption: [Structure of the periodic spectral differentiation matrix for $N = 16$. _Left_: heatmap showing the matrix entries. The Toeplitz (constant diagonal) structure is evident, as is the skew-symmetry ($D^top = -D$) with positive values (red) appearing where negative values (blue) appear in the transpose. _Right_: the first row entries $D_(0,k)$ (blue dots) plotted against the column offset $k = j - i$. The dashed red curve shows the continuous cotangent function $1/2 dot (-1)^k cot(k pi \/ N)$ from @eq-spectral-periodic, confirming that the discrete matrix entries lie exactly on this curve. The antisymmetric pattern ($D_(0,k) = -D_(0,-k)$) reflects the skew-symmetry of~$D$.],
) <fig-spectral-structure>

=== Python Implementation

The following Python code constructs the periodic spectral differentiation matrix:

```python
import numpy as np

def spectral_diff_periodic(N):
    """
    Construct the periodic spectral differentiation matrix.

    Parameters:
        N : int - Number of grid points (should be even)

    Returns:
        D : ndarray (N, N) - Differentiation matrix
        x : ndarray (N,) - Grid points on [0, 2π)
    """
    h = 2 * np.pi / N
    x = h * np.arange(N)
    D = np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            if i != j:
                D[i, j] = 0.5 * ((-1) ** (i - j)) / np.tan((i - j) * np.pi / N)

    return D, x
```

The equivalent MATLAB implementation:

```matlab
function [D, x] = spectral_diff_periodic(N)
    h = 2 * pi / N;
    x = h * (0:N-1)';
    D = zeros(N, N);

    for i = 1:N
        for j = 1:N
            if i ~= j
                D(i, j) = 0.5 * ((-1)^(i-j)) / tan((i-j) * pi / N);
            end
        end
    end
end
```

The code implementing these algorithms is available in:
- `codes/python/ch05_differentiation_matrices/spectral_matrix_periodic.py`
- `codes/matlab/ch05_differentiation_matrices/spectral_matrix_periodic.m`

== Fornberg's Recursive Algorithm <sec-fornberg>

=== Motivation: Robustness for Arbitrary Grids

The explicit formula @eq-spectral-periodic is elegant but "brittle"---it applies only to equispaced periodic grids. What if we need differentiation weights for non-equispaced nodes? Or for non-periodic problems? Or for higher-order derivatives?

The direct approach would be to solve the Vandermonde system that arises from polynomial interpolation. However, Vandermonde matrices are notoriously ill-conditioned, especially for large $N$. We need a more robust algorithm.

In 1988, Bengt Fornberg published a remarkable recursive algorithm @Fornberg1988 that computes finite difference weights for _any_ node distribution and _any_ derivative order. The algorithm is:
- Numerically stable even for large stencils
- Efficient: $O(M N^2)$ for all derivatives up to order $M$ on $N$ nodes
- Elegant: the recursion has a beautiful mathematical structure

=== The Recursive Algorithm

The key insight is that we can build up the weights incrementally. Let $delta_m^((k, n))$ denote the weight for node $x_k$ when approximating the $m$-th derivative using nodes $x_0, x_1, dots, x_n$ (a stencil of $n + 1$ points), evaluated at some point $xi$.

The algorithm proceeds level by level:
- *Level 0* ($n = 0$): A single node. For interpolation ($m = 0$), the weight is simply $1$.
- *Level 1* ($n = 1$): Two nodes. Weights are computed from level 0.
- *Level $n$*: Weights for $n + 1$ nodes are computed from level $n - 1$.

The recursion formulas are:

For $k < n$ (weights for nodes already in the previous stencil):
$ delta_m^((k, n)) = frac((xi - x_n) delta_m^((k, n-1)) - m delta_(m-1)^((k, n-1)), x_n - x_k). $ <eq-fornberg-old>

For $k = n$ (weight for the newly added node):
$ delta_m^((n, n)) = frac(c_n, c_(n-1)) (m delta_(m-1)^((n-1, n-1)) - (xi - x_(n-1)) delta_m^((n-1, n-1))), $ <eq-fornberg-new>
where $c_n = product_(ell=0)^(n-1) (x_n - x_ell)$.

The base case is $delta_0^((0, 0)) = 1$, and we define $delta_m^((k, n)) = 0$ when $m < 0$ or $m > n$.

=== The Stencil Pyramid

@fig-stencil-pyramid visualizes the recursive structure. Each level of the pyramid contains weights for one stencil size. Arrows show the data dependencies: weights at level $n$ depend only on weights from level $n - 1$.

#figure(
  image("../figures/ch05/python/stencil_pyramid.pdf", width: 85%),
  caption: [Fornberg's recursive algorithm visualized as a "stencil pyramid." Each box represents a weight $delta_m^((k, n))$ for node index $k$ in a stencil of $n + 1$ nodes. The algorithm builds from top (single node) to bottom (full stencil). Green arrows indicate updates to existing node weights via @eq-fornberg-old; red arrows show computation of new node weights via @eq-fornberg-new. The recursive structure ensures numerical stability.],
) <fig-stencil-pyramid>

The pyramid structure explains why the algorithm is numerically stable. Each weight is computed as a simple combination of weights from the previous level, avoiding the accumulation of rounding errors that plagues direct methods.

=== Python Implementation

The following code implements Fornberg's algorithm, translated from the MATLAB version by Toby Driscoll @Driscoll2007:

```python
import numpy as np

def fdweights(xi, x, m):
    """
    Compute finite difference weights using Fornberg's algorithm.

    Parameters:
        xi : float - Point where derivative is approximated
        x  : array - Node positions
        m  : int - Derivative order (0=interpolation, 1=first, etc.)

    Returns:
        w  : array - Weights for the m-th derivative approximation
    """
    x = np.asarray(x, dtype=float)
    n = len(x) - 1
    w = np.zeros(n + 1)
    x_shifted = x - xi  # Translate evaluation point to origin

    for k in range(n + 1):
        w[k] = _weight(x_shifted, m, n, k)

    return w

def _weight(x, m, j, k):
    """Recursive weight computation (evaluation point at 0)."""
    if m < 0 or m > j:
        return 0.0
    elif m == 0 and j == 0:
        return 1.0
    else:
        if k < j:
            c = (x[j] * _weight(x, m, j-1, k)
                 - m * _weight(x, m-1, j-1, k)) / (x[j] - x[k])
        else:
            beta = np.prod(x[j-1] - x[:j-1]) / np.prod(x[j] - x[:j])
            c = beta * (m * _weight(x, m-1, j-1, j-1)
                        - x[j-1] * _weight(x, m, j-1, j-1))
        return c
```

This algorithm is the universal tool for computing differentiation matrix entries on any grid.

The code implementing this algorithm is available in:
- `codes/python/ch05_differentiation_matrices/fdweights.py`
- `codes/matlab/ch05_differentiation_matrices/fdweights.m`

== Computational Étude: The Rational Trigonometric Test <sec-etude-convergence>

We now conduct a computational experiment that reveals the dramatic difference between finite difference and spectral accuracy. The goal is not merely to verify theoretical convergence rates, but to develop intuition for _why_ spectral methods achieve such remarkable precision.

=== Choosing a Test Function

The choice of test function is crucial. We need a function that is smooth and periodic (so that both finite difference and spectral methods apply), yet not so simple that it masks the differences between methods. A pure trigonometric function like $sin(k x)$ would be inappropriate: since $sin(k x)$ is an eigenfunction of both the continuous and discrete differentiation operators, spectral methods would give exact results trivially, revealing nothing about their convergence behavior.

Instead, we choose the _rational trigonometric function_:
$ u(x) = frac(1, 2 + sin(x)) quad "on" quad [0, 2 pi). $ <eq-test-function>
This function satisfies our requirements admirably. It is infinitely differentiable (analytic) on the entire real line, and clearly periodic with period $2 pi$. Most importantly, it has an _infinite_ Fourier series---unlike $sin(x)$ or finite trigonometric polynomials, its spectral content extends to all frequencies, ensuring there are no "lucky cancellations" in our numerical differentiation.

The exact derivative, computed by the quotient rule, is:
$ u'(x) = - frac(cos(x), (2 + sin(x))^2). $ <eq-test-derivative>

There is a deeper reason for choosing this particular function. Although $u(x)$ is smooth on the real line, it has singularities in the _complex plane_. The denominator $2 + sin(x)$ vanishes when $sin(x) = -2$, which has no real solutions but does have complex solutions at $x = -pi\/2 plus.minus i dot "arcsinh"(2)$. The distance from the real axis to these nearest singularities is $d = "arcsinh"(2) approx 1.44$. As we shall see, this distance controls the convergence rate of spectral methods---a beautiful connection to potential theory that we explored in @sec-potential-theory.

=== The Experiment

Our experimental design is straightforward. For a sequence of grid sizes $N = 4, 6, 8, 10, dots, 64$, we construct four differentiation matrices: three finite difference matrices of orders 2, 4, and 6 (using Fornberg's algorithm from the previous section), and the periodic spectral differentiation matrix given by @eq-spectral-periodic.

For each matrix $D$, we sample the test function at the $N$ equispaced grid points to form the vector $bold(u) = (u(x_0), u(x_1), dots, u(x_(N-1)))^top$, compute the numerical derivative $D bold(u)$, and compare it to the exact derivative values $bold(u)'_"exact" = (u'(x_0), u'(x_1), dots, u'(x_(N-1)))^top$. The error is measured in the maximum norm:
$ epsilon_N = norm(bold(u)'_"exact" - D bold(u))_infinity = max_(0 lt.eq.slant j lt.eq.slant N-1) |u'(x_j) - (D bold(u))_j|. $

=== Results and Interpretation

@fig-convergence-diff displays the results on a semi-logarithmic plot. The visual contrast between finite difference and spectral methods is striking and reveals a fundamental distinction in their convergence behavior.

#figure(
  image("../figures/ch05/python/convergence_comparison.pdf", width: 85%),
  caption: [Differentiation error comparison: finite differences versus spectral method. The test function is $u(x) = 1\/(2 + sin(x))$ on the periodic domain $[0, 2pi)$. Finite difference methods (FD2, FD4, FD6) exhibit algebraic convergence with the expected rates $O(N^(-2))$, $O(N^(-4))$, $O(N^(-6))$. The spectral method achieves geometric (exponential) convergence, reaching machine precision around $N = 20$. Dashed lines show theoretical convergence rates.],
) <fig-convergence-diff>

The finite difference methods exhibit _algebraic_ (polynomial) convergence. On a log-log plot, their error curves would appear as straight lines; on our semi-log plot, they curve downward with decreasing slope. The second-order method (FD2) shows error decreasing as $O(N^(-2))$, which is equivalent to $O(h^2)$ since $h = 2 pi \/ N$. The fourth-order method (FD4) improves to $O(N^(-4))$, and the sixth-order method (FD6) achieves $O(N^(-6))$. These rates match the theoretical predictions from Taylor series analysis: a $p$-th order finite difference scheme has truncation error $O(h^p)$.

The spectral method, in contrast, exhibits _geometric_ (exponential) convergence. Its error curve appears as a straight line on the semi-log plot, indicating that $log epsilon_N$ decreases linearly with $N$. Mathematically, $epsilon_N = O(c^(-N))$ for some constant $c > 1$. The method reaches machine precision ($approx 10^(-15)$) at around $N = 20$---with just 20 grid points, we have computed the derivative to the full precision available in double-precision arithmetic!

To appreciate what this means in practice, consider trying to achieve 15-digit accuracy with finite differences. For the second-order method, we would need to solve $C N^(-2) approx 10^(-15)$. Even with a modest constant $C approx 1$, this requires $N approx 10^(7.5) approx 30$ million grid points. For a problem in three dimensions, this would mean $30^3 times 10^(21)$ unknowns---utterly impractical. The spectral method achieves the same accuracy with 20 points.

=== Discussion: Why Does Spectral Win?

The difference between algebraic and spectral convergence is fundamental:

*Algebraic convergence* ($O(N^(-p))$): Doubling $N$ reduces the error by a factor of $2^p$. This is good, but the improvement is polynomial.

*Spectral convergence* ($O(c^(-N))$): Doubling $N$ _squares_ the error (roughly). This exponential improvement is why spectral methods can achieve machine precision with modest grid sizes.

The source of spectral convergence is the _analyticity_ of the function being approximated. For the test function @eq-test-function, the nearest singularities in the complex plane are at distance approximately $"arcsinh"(2) approx 1.44$ from the real axis. Potential theory (cf.~@sec-potential-theory) tells us that the convergence rate is controlled by this distance: the interpolation error decreases like $rho^(-N)$ where $rho = e^d$ and $d$ is the distance to the nearest singularity.

The code generating @fig-convergence-diff is available in:
- `codes/python/ch05_differentiation_matrices/convergence_comparison.py`
- `codes/matlab/ch05_differentiation_matrices/convergence_comparison.m`

== Higher-Order Derivatives

=== Second Derivatives: Squaring vs. Direct Construction

For second derivatives, we have two options:
1. *Matrix squaring*: $D^((2)) = D dot D$, where $D$ is the first-derivative matrix
2. *Direct construction*: Build $D^((2))$ directly using second-derivative weights

Matrix squaring is simpler and often sufficient. The resulting matrix $(D dot D)_(i j)$ represents the combined effect of differentiating twice. For smooth functions, this gives accurate second derivatives.

For the periodic spectral case, the second-derivative matrix has a closed form:
$ D^((2))_(j k) = cases(
  display(- frac(1, 2) (-1)^(j-k) / sin^2((j - k) pi \/ N)) & "if" j eq.not k,
  display(- frac(pi^2 (N^2 + 2), 3 (2 pi)^2)) & "if" j = k.
) $

However, matrix squaring $D^2$ is often accurate enough and more convenient.

=== Fornberg's Algorithm for Higher Derivatives

A key advantage of Fornberg's algorithm is that it handles any derivative order $m$ with no additional complexity. The same recursive structure computes interpolation weights ($m = 0$), first-derivative weights ($m = 1$), second-derivative weights ($m = 2$), and so on.

This generality is valuable when solving PDEs that involve mixed derivatives or high-order terms.

== Looking Ahead: The Non-Periodic Case

The periodic spectral matrix @eq-spectral-periodic relies crucially on the periodicity of the problem. For _non-periodic_ problems on finite intervals, such as boundary value problems with Dirichlet or Neumann conditions, we need a different approach.

The solution, which we develop in the next chapter, is to use _Chebyshev differentiation matrices_. These combine the Chebyshev nodes from @sec-chebyshev-points with differentiation matrix techniques similar to those presented here. The resulting matrices are dense but not Toeplitz, reflecting the non-uniform node distribution.

The key formulas are similar in spirit to @eq-spectral-periodic but more complex due to the clustering of Chebyshev nodes near the boundaries. The Chebyshev differentiation matrix will become our primary tool for spectral solutions of differential equations on bounded domains.

== Summary

This chapter has established differentiation matrices as the computational heart of pseudospectral methods:

+ *Matrix representation*: Differentiation of interpolating polynomials can be expressed as matrix-vector multiplication: $bold(u)' approx D bold(u)$.

+ *Finite differences as sparse approximations*: Classical FD methods correspond to sparse (banded) differentiation matrices. Higher-order FD methods use wider stencils and denser matrices.

+ *Spectral methods as the limit*: When the stencil extends to all grid points, we obtain the dense spectral differentiation matrix. This limiting viewpoint unifies finite difference and spectral approaches.

+ *Periodic spectral matrix*: For periodic problems on equispaced grids, the matrix entries are given by the cotangent formula @eq-spectral-periodic.

+ *Fornberg's algorithm*: A stable, efficient recursive algorithm computes differentiation weights for arbitrary node distributions and derivative orders.

+ *Spectral accuracy*: For smooth functions, spectral methods achieve exponentially fast convergence, vastly outperforming any fixed-order finite difference scheme.

The differentiation matrix is our passport to spectral solutions of differential equations. In the chapters ahead, we will use these matrices to solve boundary value problems, initial value problems, and eventually the partial differential equations that motivated our journey.
