// textbook/chapters/chebyshev_differentiation.typ
// Chapter 7: Chebyshev Differentiation Matrices
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Email: denys.dutykh@ku.ac.ae
// Homepage: https://www.denys-dutykh.com/
// Last modified: February 2026

#import "../styles/template.typ": dropcap, num, format-table, etude-conclusion, idx, chapter-abstract

// Enable equation numbering for this chapter

= Chebyshev Differentiation Matrices <ch-chebyshev>

#chapter-abstract(keywords: [Chebyshev differentiation matrix · Chebyshev--Gauss--Lobatto nodes · Bounded domains · Spectral convergence · Negative-sum trick · Boundary clustering])[
Fourier differentiation relies on periodicity, yet many problems in science and engineering live on bounded intervals with conditions imposed at the endpoints. This chapter constructs the Chebyshev differentiation matrix, the non-periodic counterpart of the Fourier matrix, built upon the Chebyshev--Gauss--Lobatto nodes --- the cosines of equally spaced angles --- which cluster toward the boundaries and thereby evade the Runge phenomenon. Explicit formulas are derived for every entry: the off-diagonal terms involving alternating signs and node differences, the diagonal entries obtained from an interior formula or the numerically safer negative-sum trick, and the special corner values at the endpoints. Small-$N$ examples make the construction concrete, and the matrix's structure and properties --- its nilpotent action on polynomials and its lack of symmetry --- are analysed. Two computational études close the chapter: differentiating the Witch of Agnesi, whose nearby complex poles temper the rate of convergence, and a direct demonstration of spectral convergence on an analytic function. The Chebyshev differentiation matrix established here becomes the workhorse for the boundary-value and eigenvalue problems of the chapters that follow.
]

#dropcap[Having developed the theory of differentiation matrices for periodic problems in the previous chapter, we now face a new challenge: what happens when the domain is _bounded_? Many problems in science and engineering (_e.g._ heat conduction, wave propagation, fluid dynamics) are posed on finite intervals with boundary conditions at the endpoints. For such problems, the elegant trigonometric framework of Fourier methods must give way to something new.]

The key insight of this chapter is that polynomial interpolation, when done correctly, provides the foundation for spectral methods on bounded domains. The "correct" approach, as we discovered in @ch-geometry, requires carefully chosen interpolation nodes. Equispaced points lead to the Runge phenomenon; Chebyshev points do not. This chapter builds on that foundation to construct the _Chebyshev differentiation matrix#idx("Chebyshev differentiation matrix")_ --- the non-periodic analog of the Fourier differentiation matrix. The construction follows the framework developed by Gottlieb and Orszag @GottliebOrszag1977 and presented in detail by Trefethen @Trefethen2000.

== The Non-Periodic Challenge <sec-nonperiodic>

=== From Periodic to Bounded Domains

In the previous chapter, we exploited the periodicity of the domain $[0, 2 pi)$ to construct differentiation matrices using equispaced nodes and trigonometric interpolation. The resulting matrices had beautiful structure: circulant, with entries determined by a simple cotangent formula.

For problems on a bounded interval like $[-1, 1]$, we face several new difficulties:

1. *No periodicity*: The function values at the endpoints are independent, not related by periodicity.

2. *The Runge phenomenon*: Equispaced nodes lead to wild oscillations near the boundaries, as we demonstrated dramatically in @sec-runge.

3. *Boundary conditions*: Physical problems typically impose conditions at $x = plus.minus 1$, which must be incorporated into the differentiation process.

The solution to these challenges comes from our study of interpolation theory: the Chebyshev-Gauss-Lobatto points
$ x_j = cos((j pi) / N), quad j = 0, 1, dots, N, $ <eq-cheb-nodes>
cluster near the boundaries, counteracting the Runge phenomenon and enabling spectral accuracy.

=== Grid Comparison: Equispaced vs. Chebyshev

@fig-grid-comparison illustrates the fundamental difference between equispaced and Chebyshev grids. The equispaced grid distributes points uniformly across $[-1, 1]$, while the Chebyshev grid clusters points near the boundaries according to the projection from a circle.

#figure(
  image("../figures/ch07/python/grid_comparison.pdf", width: 95%),
  caption: [Comparison of equispaced and Chebyshev-Gauss-Lobatto grids for $N = 16$ intervals. Top: equispaced points are distributed uniformly. Middle: Chebyshev points cluster near the boundaries. Bottom left: the circle projection interpretation: Chebyshev points are the projections of equally-spaced points on a semicircle. Bottom right: comparison of grid spacing near the left boundary, showing the $O(N^(-2))$ clustering of Chebyshev points.],
) <fig-grid-comparison>

The code generating @fig-grid-comparison is available in:
- `codes/python/ch07/cheb_grid_comparison.py`
- `codes/matlab/ch07/cheb_grid_comparison.m`
- `codes/julia/ch07/cheb_grid_comparison.jl`

The boundary clustering is not a minor detail, it is essential for stability. Near the boundaries, where polynomial interpolation tends to oscillate wildly, the Chebyshev grid provides many closely-spaced points to control the behavior. Near the center, where interpolation is naturally well-behaved, fewer points suffice.

The spacing of Chebyshev points near the boundary is $O(N^(-2))$, compared to $O(N^(-1))$ for equispaced points. This denser boundary clustering has important implications for time-stepping in differential equations, as we shall see in later chapters. In particular, this $O(N^(-2))$ clustering imposes a severe time-stepping restriction: explicit schemes for advection require $Delta t tilde O(N^(-2))$, and for diffusion $Delta t tilde O(N^(-4))$ --- the so-called CFL bottleneck of Chebyshev methods @Trefethen2000.

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

These formulas are derived in full by Gottlieb and Orszag @GottliebOrszag1977, Trefethen @Trefethen2000 [Chapter 6], and Boyd @Boyd2000. The matrix entries are:

*Off-diagonal entries* ($i eq.not j$):
$ (D_N)_(i j) = frac(c_i, c_j) frac((-1)^(i+j), x_i - x_j). $ <eq-cheb-offdiag>

*Diagonal entries* (interior nodes, $j = 1, dots, N-1$):
$ (D_N)_(j j) = -frac(x_j, 2(1 - x_j^2)). $ <eq-cheb-diag>

*Corner entries*:
$ (D_N)_(0 0) = frac(2 N^2 + 1, 6), quad (D_N)_(N N) = -frac(2 N^2 + 1, 6). $ <eq-cheb-corner>

=== The Negative Sum Trick

While the formulas @eq-cheb-diag and @eq-cheb-corner give exact expressions for the diagonal entries, direct evaluation can be numerically unstable. A more robust approach uses the _negative sum trick_: since the derivative of a constant function must be zero, each row of $D_N$ must sum to zero:
$ (D_N)_(j j) = - sum_(k eq.not j) (D_N)_(j k). $ <eq-negative-sum>

This ensures that $D_N bold(1) = bold(0)$ to machine precision, where $bold(1)$ is the vector of all ones. The numerical importance of this device was analyzed by Baltensperger and Berrut @BaltenspergerBerrut1999, who showed that naive evaluation of the diagonal formulas can lose several digits of accuracy; the comprehensive study by Weideman and Reddy @WeidemanReddy2000 confirmed the negative sum trick as a mandatory ingredient of robust spectral codes.

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

The Julia implementation:

```julia
function cheb_matrix(N)
    N == 0 && return (zeros(1, 1), [1.0])

    # Chebyshev-Gauss-Lobatto points
    x = [cos(π * j / N) for j in 0:N]

    # Weights: c_0 = c_N = 2, others = 1
    c = ones(N + 1)
    c[1] = 2.0; c[N+1] = 2.0

    # Off-diagonal entries
    X = repeat(x', N + 1, 1)
    dX = X .- X'
    C = c * (1.0 ./ c')
    sign_mat = ((-1.0) .^ (0:N)) * ((-1.0) .^ (0:N))'

    D = C .* sign_mat ./ (-dX .+ I(N + 1))
    D = D .- Diagonal(diag(D))

    # Negative sum trick for diagonal
    for i in 1:N+1
        D[i, i] = -sum(D[i, :])
    end
    return D, x
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
  image("../figures/ch07/python/cheb_matrix_structure.pdf", width: 95%),
  caption: [Structure of the Chebyshev differentiation matrix for $N = 16$. Left: heatmap showing the matrix entries, with red indicating positive values and blue indicating negative. The large corner entries $(D_N)_(0 0)$ and $(D_N)_(N N)$ are visible. Right: row profiles showing boundary row (red) and interior row (green). The boundary row has large entries reflecting the $O(N^2)$ corner values.],
) <fig-cheb-matrix-structure>

The code generating @fig-cheb-matrix-structure and @fig-cheb-interior-row is available in:
- `codes/python/ch07/cheb_matrix_structure.py`
- `codes/matlab/ch07/cheb_matrix_structure.m`
- `codes/julia/ch07/cheb_matrix_structure.jl`

The row profile plot in @fig-cheb-matrix-structure shows both boundary and interior rows together, but the vastly different scales make the interior row structure difficult to discern. @fig-cheb-interior-row isolates the interior row to reveal its characteristic structure.

#figure(
  image("../figures/ch07/python/cheb_matrix_interior_row.pdf", width: 70%),
  caption: [Interior row profile of the Chebyshev differentiation matrix ($i = 9$, $N = 16$). Unlike the boundary row with entries of order $O(N^2)$, interior rows have much smaller magnitudes. The largest entry in absolute value is the diagonal element $D_(9,9)$, with oscillating signs typical of differentiation stencils.],
) <fig-cheb-interior-row>

Several features are visible in @fig-cheb-interior-row. The largest entry in magnitude is the diagonal element $(D_N)_(i i)$, which for interior nodes is relatively small (recall that diagonal entries vanish for truly centered points, while off-center interior points have $O(1)$ diagonal entries). The off-diagonal entries exhibit alternating signs, reminiscent of finite difference stencils, with the magnitude decaying as we move away from the diagonal. Notably, the entries near the row ends (columns $j = 0$ and $j = N$) remain significant, reflecting the global nature of spectral differentiation: even distant boundary points contribute to the derivative at interior locations.

Unlike the sparse banded matrices of finite difference methods, the Chebyshev differentiation matrix is _dense_: every entry is generally nonzero. This is the price we pay for spectral accuracy---information from every grid point contributes to the derivative at every other point.

=== Cardinal Functions

The columns of $D_N$ have a natural interpretation: column $j$ contains the derivatives of the $j$th Lagrange cardinal function evaluated at all the grid points. @fig-cheb-cardinal illustrates this connection.

#figure(
  image("../figures/ch07/python/cheb_cardinal.pdf", width: 95%),
  caption: [Chebyshev cardinal functions (Lagrange basis polynomials). Left: several cardinal functions for $N = 10$, each peaking at value $1$ at its corresponding node and vanishing at all others. Right: a single cardinal function with tangent lines at the grid points---the slopes of these tangent lines are precisely the entries in the corresponding column of the differentiation matrix.],
) <fig-cheb-cardinal>

The code generating @fig-cheb-cardinal is available in:
- `codes/python/ch07/cheb_cardinal.py`
- `codes/matlab/ch07/cheb_cardinal.m`
- `codes/julia/ch07/cheb_cardinal.jl`

The barycentric form of Lagrange interpolation, highlighted by Berrut and Trefethen @BerrutTrefethen2004, provides an alternative route to computing these derivatives that decouples the node distribution from the algorithmic structure.

== Computational Étude 7.1: The Witch of Agnesi <sec-witch>

=== A Smooth Test Function

To demonstrate spectral differentiation in action, we use the _Witch of Agnesi#idx("Witch of Agnesi")_:
$ u(x) = frac(1, 1 + 4 x^2), $
with exact derivative
$ u'(x) = frac(-8 x, (1 + 4 x^2)^2). $

This function is smooth and analytic on $[-1, 1]$, with poles at $x = plus.minus i\/2$ in the complex plane. The distance from $[-1, 1]$ to the nearest singularity determines the rate of exponential convergence. This geometric picture is a direct application of the Bernstein ellipse theory developed in @ch-smoothness; see also Trefethen @Trefethen2000 [Chapter 8] for a thorough treatment.

@fig-cheb-diff-demo shows the function and its spectral derivative approximation for $N = 10$ and $N = 20$ grid points.

#figure(
  image("../figures/ch07/python/cheb_diff_demo.pdf", width: 95%),
  caption: [Chebyshev spectral differentiation of the Witch of Agnesi $u(x) = 1\/(1 + 4x^2)$. Top row: function values at Chebyshev points. Bottom row: comparison of exact and spectral derivatives, with maximum errors indicated. The error decreases exponentially with $N$.],
) <fig-cheb-diff-demo>

The code generating @fig-cheb-diff-demo is available in:
- `codes/python/ch07/cheb_diff_demo.py`
- `codes/matlab/ch07/cheb_diff_demo.m`
- `codes/julia/ch07/cheb_diff_demo.jl`

The exponential convergence is evident from @tab-witch-errors, which shows the maximum differentiation error for various values of $N$.

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    {
      show table: format-table(auto, auto, auto, auto)
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
        [4], num[8.5e-1], [20], num[1.9e-3],
        [6], num[4.8e-1], [24], num[3.3e-4],
        [8], num[2.4e-1], [28], num[5.6e-5],
        [10], num[1.2e-1], [32], num[9.4e-6],
        [12], num[5.3e-2], [36], num[1.5e-6],
        [14], num[2.4e-2], [40], num[2.5e-7],
        [16], num[1.0e-2], [50], num[5.8e-9],
      )
    },
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

The Julia implementation:

```julia
function differentiate_witch(N)
    # Differentiate u = 1/(1+4x²) using Chebyshev spectral method.
    D, x = cheb_matrix(N)

    # Function and exact derivative
    u = @. 1.0 / (1.0 + 4.0 * x^2)
    du_exact = @. -8.0 * x / (1.0 + 4.0 * x^2)^2

    # Spectral derivative
    du_spectral = D * u

    # Maximum error
    max_error = maximum(abs.(du_spectral .- du_exact))
    return x, du_spectral, du_exact, max_error
end
```

#etude-conclusion[
  The Witch of Agnesi $u(x) = 1 \/ (1 + 4 x^2)$ has poles at $x = plus.minus i \/ 2$. Applying the exact formula $rho = y + sqrt(y^2 + 1)$ with $y = 1\/2$ gives the Bernstein-ellipse parameter $rho = (1 + sqrt(5)) \/ 2 approx 1.618$ --- the golden ratio. (The naïve approximation $rho approx 1 + y = 1.5$ is off by roughly 7% for singularities this close to the real axis.) The error in the table accordingly decreases by roughly a factor of $1.618$ per unit of $N$, and by $N = 50$ has dropped to $approx 6 times 10^(-9)$. The figure confirms that *even at $N = 10$* the spectral derivative is visually indistinguishable from the exact curve, despite the relatively close singularities. The étude highlights a key difference from the periodic case: the Chebyshev grid *clusters near the boundaries* of $[-1, 1]$, which is exactly where the Lagrange basis polynomials would oscillate most if equispaced nodes were used. *Boundary clustering is a built-in stabiliser*, allowing spectral accuracy without the periodicity assumption.
]

== Computational Étude 7.2: Spectral Convergence <sec-convergence>

=== Four Functions of Increasing Smoothness

The rate of spectral convergence#idx("spectral convergence") depends critically on the smoothness of the function being differentiated. The theoretical framework developed in @ch-smoothness applies equally to Chebyshev methods on bounded domains: Theorem 1 establishes that smoother functions have more rapidly decaying spectral coefficients, and Theorem 4 translates this decay into bounds on the differentiation error. The convergence rates we report below measure how fast the maximum error $norm(D_N bold(v) - u'(bold(x)))_infinity$ decreases as $N$ increases.

To illustrate these principles, we examine four test functions with different regularity:

1. $|x|^(5\/2)$: This function has fractional Hölder regularity $s = 5\/2$. The first two derivatives,
$ (|x|^(5\/2))' = frac(5, 2)|x|^(3\/2) op("sgn")(x), quad (|x|^(5\/2))'' = frac(15, 4)|x|^(1\/2), $
are continuous, but the third derivative $(15\/8)|x|^(-1\/2) op("sgn")(x)$ is unbounded at $x = 0$. By the generalization of Theorem 1 to non-integer regularity, the Chebyshev coefficients decay as $O(n^(-s-1)) = O(n^(-7\/2))$, and by Theorem 4 the differentiation error is $O(N^(-s)) = O(N^(-5\/2)) = O(N^(-2.5))$. The fractional exponent in the function directly determines the convergence rate.

2. $e^(-1\/(1-x^2))$: The "bump function" is infinitely differentiable ($C^oo$) on $[-1, 1]$, vanishing together with all its derivatives at $x = plus.minus 1$. However, it is not analytic at the endpoints: no Taylor series converges there. By Theorem 1(b), $C^oo$ functions have coefficients that decay faster than any algebraic rate, $O(n^(-m))$ for all $m$. The differentiation error thus converges superalgebraically (faster than any power of $N$) but not exponentially.

3. $tanh(5 x)$: This function is real-analytic on $[-1, 1]$, but $tanh(z)$ has poles at $z = plus.minus i pi \/ 2$ in the complex plane. Scaling by $5$ moves the nearest poles to $z = plus.minus i pi \/ 10$. By Theorem 1(c), analyticity in a strip of half-width $a = pi\/10$ yields Chebyshev coefficients decaying as $O(e^(-a n))$. For the Chebyshev ellipse with parameter $rho$, the convergence rate is $O(rho^(-N))$ where $rho approx 1 + a = 1 + pi\/10 approx 1.31$.

4. $x^8$: A polynomial of degree $8$ is its own interpolant for $N gt.eq.slant 8$. This is the bounded-domain analog of band-limited functions: since $x^8$ has only finitely many nonzero Chebyshev coefficients (up to $T_8$), the spectral derivative is _exact_ once $N$ is large enough to resolve all terms. This corresponds to Theorem 4(d) for the periodic case.

@fig-convergence-waterfall displays the maximum differentiation error versus $N$ for these four functions.

#figure(
  image("../figures/ch07/python/convergence_waterfall.pdf", width: 95%),
  caption: [Convergence of the spectral differentiation error $norm(D_N bold(v) - u'(bold(x)))_infinity$ as a function of $N$ for four test functions. Top left: $|x|^(5\/2)$ shows algebraic convergence $O(N^(-2.5))$ due to limited smoothness (Hölder regularity $5\/2$). Top right: the bump function $e^(-1\/(1-x^2))$ achieves superalgebraic (faster than any power) but not exponential convergence, consistent with $C^oo$ but non-analytic behavior. Bottom left: $tanh(5x)$ demonstrates exponential convergence until machine precision, as expected for analytic functions. Bottom right: the polynomial $x^8$ is differentiated exactly for $N gt.eq.slant 8$.],
) <fig-convergence-waterfall>

#etude-conclusion[
  The figure is a *vivid catalogue of the four convergence regimes* predicted by the theory of @ch-smoothness. Algebraic ($|x|^(5 \/ 2)$, slope $-2.5$ matching the Hölder exponent), super-algebraic ($C^infinity$ bump, faster than any power but slower than exponential), geometric ($tanh(5 x)$, textbook straight-line descent to machine precision), and *exact* (the polynomial $x^8$ once $N gt.eq.slant 8$). The practical message is that the *promised* exponential convergence of spectral methods materialises only for analytic functions; for less smooth targets, convergence is still rapid but algebraic, with the rate set by the regularity. The hierarchy is not just theoretical --- monitoring the convergence rate during a computation lets one *diagnose* the regularity of the solution, a technique invaluable for validating numerical results and detecting hidden singularities.
]

The code generating @fig-convergence-waterfall is available in:
- `codes/python/ch07/cheb_convergence.py`
- `codes/matlab/ch07/cheb_convergence.m`
- `codes/julia/ch07/cheb_convergence.jl`

== A non-exhaustive literature overview

The intellectual roots of Chebyshev spectral differentiation reach back to Fourier's 1822 treatise on heat conduction @Fourier1822, which demonstrated that functions --- even those with discontinuities --- could be represented as infinite sums of trigonometric modes. While Fourier series are naturally suited to periodic domains, their direct application to bounded intervals destroys uniform convergence at the boundaries. The resolution came with the work of Cornelius Lanczos#idx("Lanczos") in the 1930s @Lanczos1938, who recognized that Chebyshev polynomials --- defined as $T_n (x) = cos(n arccos x)$ --- inherit the favourable convergence properties of Fourier series while providing a robust basis for non-periodic functions on $[-1, 1]$. This insight remained largely theoretical until the rediscovery of the Fast Fourier Transform by Cooley and Tukey @Cooley1965 made the transformation between physical and spectral space computationally feasible in $O(N log N)$ operations. The subsequent formalization by Gottlieb and Orszag @GottliebOrszag1977 established the mathematical rigour that underpins modern practice.

The transition from theory to reliable computation required several algorithmic refinements. Baltensperger and Berrut @BaltenspergerBerrut1999 provided a careful analysis of rounding errors in the pseudospectral differentiation matrix, demonstrating that the naive evaluation of diagonal entries can lose several digits of accuracy and establishing the negative sum trick as essential practice. Weideman and Reddy @WeidemanReddy2000 consolidated these insights into a widely-used MATLAB differentiation matrix suite that made spectral methods accessible to a broad community of practitioners. Berrut and Trefethen @BerrutTrefethen2004 highlighted the extraordinary stability of the barycentric form of Lagrange interpolation, providing an alternative algorithmic pathway that decouples the choice of nodes from the structure of the differentiation algorithm. Trefethen's _Spectral Methods in MATLAB_ @Trefethen2000 remains the most accessible introduction to the pseudospectral approach, presenting the key ideas through compact, readable programs.

Several comprehensive textbooks have codified the theory for practitioners. The foundational monograph by Gottlieb and Orszag @GottliebOrszag1977 provides the rigorous mathematical framework, including convergence theorems and stability analysis. Boyd @Boyd2000 offers an encyclopaedic treatment of Chebyshev and Fourier methods with an emphasis on practical algorithms and an extensive catalogue of convergence rates. Fornberg @Fornberg1996 gives a unified perspective on pseudospectral and finite-difference methods, while Canuto _et al._ @Canuto2006 develop the functional-analytic foundations in full generality.

In recent years, Chebyshev spectral methods have found striking applications beyond their traditional home in fluid dynamics. Batic, Dutykh, and Beek @Batic2025a employed Chebyshev differentiation matrices to compute the quasinormal mode frequencies of noncommutative geometry-inspired wormholes, a problem in general relativity where the high accuracy of spectral methods is essential for resolving the small imaginary parts of the complex eigenfrequencies. Melia _et al._ @Melia2026 developed hardware-accelerated hierarchical Poincaré--Steklov solvers that exploit the structure of spectral discretizations on modern GPU architectures, achieving significant speedups for large-scale elliptic problems.

Contemporary research is extending Chebyshev differentiation in several new directions. Zayernouri _et al._ @Zayernouri2024 have developed spectral and spectral element methods for fractional differential equations, exploiting the global nature of Chebyshev polynomials to handle the non-local operators that arise in anomalous diffusion and memory effects. Tasnin @Tasnin2025 has explored the integration of Chebyshev spectral methods with physics-informed neural networks, embedding spectral layers into deep learning architectures to overcome the spectral bias that limits standard networks in resolving high-frequency features. These developments, together with emerging work on quantum algorithms for spectral linear systems, suggest that Chebyshev differentiation matrices will remain central tools in computational science for the foreseeable future.

== Summary

This chapter has developed the Chebyshev differentiation matrix for spectral methods on bounded domains. Key points include:

- *Chebyshev points* @eq-cheb-nodes provide the optimal node distribution, clustering near boundaries to avoid the Runge phenomenon.

- *The negative sum trick* @eq-negative-sum ensures numerical stability by computing diagonal entries from row sums.

- *Matrix structure*: The Chebyshev differentiation matrix is dense, with $O(N^2)$ corner entries.

- *Spectral convergence* depends on smoothness: exponential for analytic functions, algebraic for functions with limited regularity.

In the next chapter, we will use these differentiation matrices to solve boundary value problems---the natural next step in applying spectral methods to differential equations.

== Exercises <sec-chebyshev-differentiation-exercises>

*Exercise 7.1* (_Explicit Construction and Verification of $D_N$_). Construct the Chebyshev differentiation matrix $D_N$ explicitly for $N = 4$. (a) Write out the $5 times 5$ matrix and verify the entries using the formulas from this chapter. (b) Verify the negative sum trick: confirm that $sum_(j=0)^N D_(i j) = 0$ for each row $i$. (c) Test the matrix by differentiating $f(x) = T_3 (x)$ (the third Chebyshev polynomial) and verifying that the result is $3 U_2 (x)$ evaluated at the Chebyshev points, where $U_2$ is the Chebyshev polynomial of the second kind.

*Exercise 7.2* (_Condition Number of $D_N$_). Compute the 2-norm condition number $kappa(D_N) = ||D_N||_2 ||D_N^(-1)||_2$ of the Chebyshev differentiation matrix for $N = 4, 8, 16, 32, 64$. (a) Plot $kappa(D_N)$ versus $N$ on a log-log scale and verify the scaling $kappa(D_N) = cal(O)(N^(2))$. (b) Compute the condition number of $D_N^2$ and verify $kappa(D_N^2) = cal(O)(N^4)$. (c) Explain why this rapid growth of condition numbers motivates the preconditioning techniques and alternative formulations discussed in later chapters.

*Exercise 7.3* (_Chebyshev Interpolation of $|x|$ and Convergence Analysis_). Interpolate $f(x) = |x|$ at $N + 1$ Chebyshev points for $N = 8, 16, 32, 64, 128, 256$. (a) Plot the interpolant and the error. Note the slow convergence caused by the cusp at $x = 0$. (b) Measure the maximum error and verify algebraic convergence $cal(O)(N^(-1))$. (c) Now interpolate $g(x) = x |x|$ ($C^1$ but not $C^2$) and verify that the convergence improves to $cal(O)(N^(-2))$. This illustrates the direct connection between regularity and spectral convergence rate.

*Exercise 7.4* (_Clenshaw--Curtis vs Gauss--Chebyshev Quadrature_). Implement both Clenshaw--Curtis quadrature (based on Chebyshev--Lobatto nodes) and Gauss--Chebyshev quadrature (based on Chebyshev roots). (a) Compute $integral_(-1)^1 e^x dif x$ and $integral_(-1)^1 1\/(1 + 25x^2) dif x$ with both rules for $N = 4, 8, 16, 32, 64$ and plot the error versus $N$. (b) Verify that both methods achieve spectral convergence for analytic integrands. (c) For the polynomial $f(x) = x^(2N)$, show that Gauss quadrature with $N + 1$ nodes is exact while Clenshaw--Curtis requires $N$ nodes to be exact. Explain the factor-of-two difference in terms of polynomial exactness.
