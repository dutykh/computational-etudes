// textbook/chapters/chebyshev_differentiation.typ
// Chapter 6: Chebyshev Differentiation Matrices
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Last modified: January 2026

#import "../styles/template.typ": dropcap

// Enable equation numbering for this chapter
#set math.equation(numbering: "(1)")

= Chebyshev Differentiation Matrices <ch-chebyshev>

#dropcap[Having developed the theory of differentiation matrices for periodic problems in the previous chapter, we now face a new challenge: what happens when the domain is _bounded_? Many problems in science and engineering---heat conduction, wave propagation, fluid dynamics---are posed on finite intervals with boundary conditions at the endpoints. For such problems, the elegant trigonometric framework of Fourier methods must give way to something new.]

The key insight of this chapter is that polynomial interpolation, when done correctly, provides the foundation for spectral methods on bounded domains. The "correct" approach, as we discovered in @ch-geometry, requires carefully chosen interpolation nodes. Equispaced points lead to the Runge phenomenon; Chebyshev points do not. This chapter builds on that foundation to construct the _Chebyshev differentiation matrix_---the non-periodic analog of the Fourier differentiation matrix.

== The Non-Periodic Challenge <sec-nonperiodic>

=== From Periodic to Bounded Domains

In the previous chapter, we exploited the periodicity of the domain $[0, 2 pi)$ to construct differentiation matrices using equispaced nodes and trigonometric interpolation. The resulting matrices had beautiful structure: circulant, with entries determined by a simple cotangent formula.

For problems on a bounded interval like $[-1, 1]$, we face several new difficulties:

1. *No periodicity*: The function values at the endpoints are independent, not related by periodicity.

2. *The Runge phenomenon*: Equispaced nodes lead to wild oscillations near the boundaries, as we demonstrated dramatically in @sec-runge.

3. *Boundary conditions*: Physical problems typically impose conditions at $x = plus.minus 1$, which must be incorporated into the differentiation process.

The solution to these challenges comes from our study of interpolation theory: the Chebyshev-Gauss-Lobatto points
$ x_j = cos(j pi \/ N), quad j = 0, 1, dots, N, $ <eq-cheb-nodes>
cluster near the boundaries, counteracting the Runge phenomenon and enabling spectral accuracy.

=== Grid Comparison: Equispaced vs. Chebyshev

@fig-grid-comparison illustrates the fundamental difference between equispaced and Chebyshev grids. The equispaced grid distributes points uniformly across $[-1, 1]$, while the Chebyshev grid clusters points near the boundaries according to the projection from a circle.

#figure(
  image("../figures/ch06/python/grid_comparison.pdf", width: 95%),
  caption: [Comparison of equispaced and Chebyshev-Gauss-Lobatto grids for $N = 16$ intervals. Top: equispaced points are distributed uniformly. Middle: Chebyshev points cluster near the boundaries. Bottom left: the circle projection interpretation---Chebyshev points are the projections of equally-spaced points on a semicircle. Bottom right: comparison of grid spacing near the left boundary, showing the $O(N^(-2))$ clustering of Chebyshev points.],
) <fig-grid-comparison>

The boundary clustering is not a minor detail---it is essential for stability. Near the boundaries, where polynomial interpolation tends to oscillate wildly, the Chebyshev grid provides many closely-spaced points to control the behavior. Near the center, where interpolation is naturally well-behaved, fewer points suffice.

The spacing of Chebyshev points near the boundary is $O(N^(-2))$, compared to $O(N^(-1))$ for equispaced points. This denser boundary clustering has important implications for time-stepping in differential equations, as we shall see in later chapters.

== The Chebyshev Differentiation Matrix <sec-cheb-matrix>

=== Differentiation via Interpolation

The construction of the Chebyshev differentiation matrix follows the same principle as in the periodic case: interpolate, then differentiate.

Given function values $bold(v) = (v_0, v_1, dots, v_N)^top$ at the Chebyshev points @eq-cheb-nodes, we first construct the unique interpolating polynomial $p(x)$ of degree at most $N$ satisfying $p(x_j) = v_j$. The derivative approximation at the grid points is then
$ bold(w) = D_N bold(v), quad "where" quad w_i = p'(x_i). $

The matrix $D_N$ is the _Chebyshev differentiation matrix_. Its entries can be computed from the Lagrange basis polynomials:
$ (D_N)_(i j) = L'_j (x_i), $
where $L_j$ is the Lagrange interpolating polynomial centered at $x_j$.

=== Explicit Formulas

The entries of the Chebyshev differentiation matrix can be written explicitly. Define the weights
$ c_j = cases(2 & "if" j = 0 "or" j = N, 1 & "otherwise.") $

Then the matrix entries are:

*Off-diagonal entries* ($i eq.not j$):
$ (D_N)_(i j) = frac(c_i, c_j) frac((-1)^(i+j), x_i - x_j). $ <eq-cheb-offdiag>

*Diagonal entries* (interior nodes, $j = 1, dots, N-1$):
$ (D_N)_(j j) = -frac(x_j, 2(1 - x_j^2)). $ <eq-cheb-diag>

*Corner entries*:
$ (D_N)_(0 0) = frac(2 N^2 + 1, 6), quad (D_N)_(N N) = -frac(2 N^2 + 1, 6). $ <eq-cheb-corner>

=== The Negative Sum Trick

While the formulas @eq-cheb-diag and @eq-cheb-corner give exact expressions for the diagonal entries, direct evaluation can be numerically unstable. A more robust approach uses the _negative sum trick_: since the derivative of a constant function must be zero, each row of $D_N$ must sum to zero:
$ (D_N)_(j j) = - sum_(k eq.not j) (D_N)_(j k). $ <eq-negative-sum>

This ensures that $D_N bold(1) = bold(0)$ to machine precision, where $bold(1)$ is the vector of all ones.

The following code implements the Chebyshev differentiation matrix:

```python
def cheb_matrix(N):
    """Chebyshev differentiation matrix."""
    if N == 0:
        return np.array([[0.0]]), np.array([1.0])

    # Chebyshev-Gauss-Lobatto points
    x = np.cos(np.pi * np.arange(N + 1) / N)

    # Weights: c_0 = c_N = 2, others = 1
    c = np.ones(N + 1)
    c[0], c[N] = 2.0, 2.0

    # Off-diagonal entries
    X = np.tile(x, (N + 1, 1))
    dX = X - X.T
    C = np.outer(c, 1.0 / c)
    sign = np.outer((-1.0)**np.arange(N + 1), (-1.0)**np.arange(N + 1))

    with np.errstate(divide='ignore'):
        D = C * sign / (-dX)

    # Negative sum trick for diagonal
    np.fill_diagonal(D, 0.0)
    D[np.diag_indices(N + 1)] = -np.sum(D, axis=1)

    return D, x
```

```matlab
function [D, x] = cheb_matrix(N)
% Chebyshev differentiation matrix.
    if N == 0
        D = 0; x = 1; return
    end

    % Chebyshev-Gauss-Lobatto points
    x = cos(pi * (0:N)' / N);

    % Weights: c_0 = c_N = 2, others = 1
    c = ones(N+1, 1);
    c(1) = 2; c(N+1) = 2;

    % Off-diagonal entries
    X = repmat(x, 1, N+1);
    dX = X - X';
    C = c * (1 ./ c');
    sign = (-1).^((0:N)' + (0:N));

    D = C .* sign ./ (-dX);

    % Negative sum trick for diagonal
    D(1:N+2:end) = 0;
    D(1:N+2:end) = -sum(D, 2);
end
```

== Small-$N$ Examples <sec-small-n>

=== Hand Calculations

To develop intuition, let us compute the smallest Chebyshev matrices by hand.

For $N = 1$, we have two nodes: $x_0 = 1$ and $x_1 = -1$. The only polynomial passing through two points is a line, and its derivative is constant:
$ D_1 = mat(1/2, -1/2; 1/2, -1/2). $

For $N = 2$, we have three nodes: $x_0 = 1$, $x_1 = 0$, and $x_2 = -1$. The middle row of $D_2$ is particularly illuminating:
$ D_2 = mat(3/2, -2, 1/2; 1/2, 0, -1/2; -1/2, 2, -3/2). $

The middle row $(1/2, 0, -1/2)$ is exactly the _centered finite difference_ formula! This reveals a beautiful connection: at interior points where the grid happens to be locally symmetric, the spectral method reduces to the familiar finite difference formula.

== Matrix Structure and Properties <sec-matrix-structure>

=== The Dense Matrix

@fig-cheb-matrix-structure visualizes the structure of the Chebyshev differentiation matrix for $N = 16$.

#figure(
  image("../figures/ch06/python/cheb_matrix_structure.pdf", width: 95%),
  caption: [Structure of the Chebyshev differentiation matrix for $N = 16$. Left: heatmap showing the matrix entries, with red indicating positive values and blue indicating negative. The large corner entries $(D_N)_(0 0)$ and $(D_N)_(N N)$ are visible. Right: row profiles showing boundary row (red) and interior row (green). The boundary row has large entries reflecting the $O(N^2)$ corner values.],
) <fig-cheb-matrix-structure>

Unlike the sparse banded matrices of finite difference methods, the Chebyshev differentiation matrix is _dense_: every entry is generally nonzero. This is the price we pay for spectral accuracy---information from every grid point contributes to the derivative at every other point.

=== Cardinal Functions

The columns of $D_N$ have a natural interpretation: column $j$ contains the derivatives of the $j$th Lagrange cardinal function evaluated at all the grid points. @fig-cheb-cardinal illustrates this connection.

#figure(
  image("../figures/ch06/python/cheb_cardinal.pdf", width: 95%),
  caption: [Chebyshev cardinal functions (Lagrange basis polynomials). Left: several cardinal functions for $N = 10$, each peaking at value $1$ at its corresponding node and vanishing at all others. Right: a single cardinal function with tangent lines at the grid points---the slopes of these tangent lines are precisely the entries in the corresponding column of the differentiation matrix.],
) <fig-cheb-cardinal>

== Demonstration: The Witch of Agnesi <sec-witch>

=== A Smooth Test Function

To demonstrate spectral differentiation in action, we use the _Witch of Agnesi_:
$ u(x) = frac(1, 1 + 4 x^2), $
with exact derivative
$ u'(x) = frac(-8 x, (1 + 4 x^2)^2). $

This function is smooth and analytic on $[-1, 1]$, with poles at $x = plus.minus i\/2$ in the complex plane. The distance from $[-1, 1]$ to the nearest singularity determines the rate of exponential convergence.

@fig-cheb-diff-demo shows the function and its spectral derivative approximation for $N = 10$ and $N = 20$ grid points.

#figure(
  image("../figures/ch06/python/cheb_diff_demo.pdf", width: 95%),
  caption: [Chebyshev spectral differentiation of the Witch of Agnesi $u(x) = 1\/(1 + 4x^2)$. Top row: function values at Chebyshev points. Bottom row: comparison of exact and spectral derivatives, with maximum errors indicated. The error decreases exponentially with $N$.],
) <fig-cheb-diff-demo>

The exponential convergence is evident from @tab-witch-errors, which shows the maximum differentiation error for various values of $N$.

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    table(
      columns: 4,
      align: (center, center, center, center),
      inset: (x: 1em, y: 0.6em),
      stroke: none,
      table.hline(stroke: 0.75pt + rgb("#142D6E")),
      table.header(
        table.cell(fill: rgb("#142D6E").lighten(85%))[*$N$*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Max error*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*$N$*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Max error*],
      ),
      table.hline(stroke: 0.5pt + luma(180)),
      [4], [$8.5 times 10^(-1)$], [20], [$1.9 times 10^(-3)$],
      [6], [$4.8 times 10^(-1)$], [24], [$3.3 times 10^(-4)$],
      [8], [$2.4 times 10^(-1)$], [28], [$5.6 times 10^(-5)$],
      [10], [$1.2 times 10^(-1)$], [32], [$9.4 times 10^(-6)$],
      [12], [$5.3 times 10^(-2)$], [36], [$1.5 times 10^(-6)$],
      [14], [$2.4 times 10^(-2)$], [40], [$2.5 times 10^(-7)$],
      [16], [$1.0 times 10^(-2)$], [50], [$5.8 times 10^(-9)$],
    ),
  ),
  caption: [Maximum differentiation error for the Witch of Agnesi $u(x) = 1\/(1 + 4x^2)$. The error decreases exponentially with $N$.],
) <tab-witch-errors>

The following code demonstrates spectral differentiation:

```python
def differentiate_witch(N):
    """Differentiate u = 1/(1+4x²) using Chebyshev spectral method."""
    D, x = cheb_matrix(N)

    # Function and exact derivative
    u = 1.0 / (1.0 + 4.0 * x**2)
    du_exact = -8.0 * x / (1.0 + 4.0 * x**2)**2

    # Spectral derivative
    du_spectral = D @ u

    # Maximum error
    max_error = np.max(np.abs(du_spectral - du_exact))
    return x, du_spectral, du_exact, max_error
```

```matlab
function [x, du_spectral, du_exact, max_error] = differentiate_witch(N)
% Differentiate u = 1/(1+4x²) using Chebyshev spectral method.
    [D, x] = cheb_matrix(N);

    % Function and exact derivative
    u = 1 ./ (1 + 4*x.^2);
    du_exact = -8*x ./ (1 + 4*x.^2).^2;

    % Spectral derivative
    du_spectral = D * u;

    % Maximum error
    max_error = max(abs(du_spectral - du_exact));
end
```

== Spectral Convergence <sec-convergence>

=== Four Functions of Increasing Smoothness

The rate of spectral convergence depends critically on the smoothness of the function being differentiated. To illustrate this, we examine four test functions with different regularity:

1. $|x|^(5\/2)$: Second derivative continuous, third derivative in bounded variation. Expected convergence: $O(N^(-2.5))$.

2. $e^(-1\/(1-x^2))$: The "bump function," infinitely differentiable ($C^oo$) but not analytic at $x = plus.minus 1$. Expected convergence: faster than any algebraic rate, but not exponential.

3. $tanh(5 x)$: Analytic on $[-1, 1]$ with poles at $x = plus.minus i pi \/ 10$. Expected convergence: exponential, $O(rho^(-N))$ with $rho = 1 + pi\/10 approx 1.31$.

4. $x^8$: Polynomial of degree 8. Expected convergence: exact for $N gt.eq.slant 8$.

@fig-convergence-waterfall displays the maximum differentiation error versus $N$ for these four functions.

#figure(
  image("../figures/ch06/python/convergence_waterfall.pdf", width: 95%),
  caption: [Spectral convergence for four functions of increasing smoothness. Top left: $|x|^(5\/2)$ shows algebraic convergence $O(N^(-2.5))$ due to limited smoothness. Top right: the bump function $e^(-1\/(1-x^2))$ achieves superalgebraic (faster than any power) but not exponential convergence. Bottom left: $tanh(5x)$ demonstrates exponential convergence until machine precision. Bottom right: the polynomial $x^8$ is differentiated exactly for $N gt.eq.slant 8$.],
) <fig-convergence-waterfall>

The message is clear: spectral methods achieve their promised exponential convergence only for analytic functions. For less smooth functions, convergence is still rapid but algebraic, with the rate determined by the degree of smoothness.

== Summary

This chapter has developed the Chebyshev differentiation matrix for spectral methods on bounded domains. Key points include:

- *Chebyshev points* @eq-cheb-nodes provide the optimal node distribution, clustering near boundaries to avoid the Runge phenomenon.

- *The negative sum trick* @eq-negative-sum ensures numerical stability by computing diagonal entries from row sums.

- *Matrix structure*: The Chebyshev differentiation matrix is dense, with $O(N^2)$ corner entries.

- *Spectral convergence* depends on smoothness: exponential for analytic functions, algebraic for functions with limited regularity.

In the next chapter, we will use these differentiation matrices to solve boundary value problems---the natural next step in applying spectral methods to differential equations.
