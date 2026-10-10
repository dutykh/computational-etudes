// textbook/chapters/mise_en_bouche.typ
// Chapter 3: Mise en Bouche
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Email: denys.dutykh@ku.ac.ae
// Homepage: https://www.denys-dutykh.com/
// Last modified: February 2026
#import "../styles/template.typ": dropcap, num, format-table, etude-conclusion, idx, chapter-abstract, exercise, hint-for

= Mise en Bouche

#chapter-abstract(keywords: [Method of weighted residuals · Collocation · Galerkin method · Spectral approximation · Residual minimisation · Trial functions])[
Conceived as an intellectual appetizer, this chapter distils the essential mechanics of spectral methods into transparent, hand-computable examples with only two or three unknowns. Rather than confronting high-degree expansions at once, the reader watches every step of the approximation unfold by direct calculation. The unifying framework is the Method of Weighted Residuals, which represents the unknown solution as a finite sum of basis functions and forces the residual to vanish in some averaged sense; different choices of weighting recover the collocation (pseudospectral) and Galerkin families that recur throughout the book. A first worked étude constructs a trial function that satisfies the boundary conditions automatically and enforces the residual at interior collocation points, converting a differential equation into a small, explicitly solvable linear system. A second étude contrasts collocation with Galerkin on the same problem, exposing their differing accuracy and bookkeeping. The chapter thereby builds the intuition --- how spectral methods turn calculus into algebra --- that the systematic Fourier and Chebyshev machinery of later chapters will formalise and scale.
]

#dropcap[In French cuisine, a _mise en bouche_ is a small appetizer offered by the chef to stimulate the palate before the main courses arrive. In this chapter we offer a similar intellectual appetizer: a compact, self-contained taste of spectral methods that illuminates the essential mechanics before we develop the full machinery of Fourier and Chebyshev approximation.]

Rather than jumping immediately to high-degree polynomials with $N = 100$, we perform hand calculations with just $N = 2$ or $N = 3$ unknowns. This low-dimensional setting makes every step transparent. We can verify each formula by direct computation and gain intuition that will guide us through the more sophisticated developments to come.

The techniques presented here follow the classical exposition in @Boyd2000, adapted to our pedagogical goals. The unifying framework of the Method of Weighted Residuals was first systematized in the seminal review by @Finlayson1966 and later expanded in @Finlayson1972. This framework connects the collocation#idx("collocation") (pseudospectral) approach we favor in this book with the Galerkin method#idx("Galerkin method")s that dominate finite element analysis.

== The Method of Weighted Residuals

=== Series Expansions and the Residual Function

The central idea of spectral methods is to approximate the unknown function $u(x)$ by a finite sum of basis functions:
$ u(x) approx u_N (x) = sum_(n=0)^N a_n phi_n (x), $ <eq-series-expansion>
where ${phi_n (x)}$ are known basis functions and ${a_n}$ are unknown coefficients to be determined.

When we substitute this approximation into a differential equation
$ cal(L) u = f(x), $
where $cal(L)$ is a linear differential operator, the result is generally not zero. The _residual#idx("residual") function_ measures this discrepancy:
$ R(x; a_0, a_1, dots, a_N) = cal(L) u_N - f. $ <eq-residual>

For the exact solution, $R(x) equiv 0$. The challenge is to choose the coefficients ${a_n}$ so that the residual is as small as possible. Different spectral methods correspond to different strategies for minimizing this residual.

=== Two Minimization Strategies

The two most important strategies are:

+ *Collocation (Pseudospectral) Method*: Force the residual to be exactly zero at a set of carefully chosen points ${x_j}$, called collocation points:
  $ R(x_j; a_0, dots, a_N) = 0, quad j = 1, 2, dots, N+1. $
  This gives $N+1$ equations for $N+1$ unknowns.

+ *Galerkin Method*: Require the residual to be orthogonal to each basis function (a concept originating with @Galerkin1915) in the sense of a weighted inner product#idx("weighted inner product"):
  $ integral_(-1)^1 R(x) phi_k (x) w(x) dif x = 0, quad k = 0, 1, dots, N, $
  where $w(x)$ is a weight function (often $w(x) = 1$ for polynomial bases).

Both methods convert the differential equation into a system of algebraic equations for the unknown coefficients. The collocation approach is simpler to implement and handles nonlinear terms easily, while the Galerkin approach often provides better global accuracy in weighted norms.

== Computational Étude 3.1: A First Collocation Example

We illustrate the collocation method with a complete worked example that can be verified by hand calculation.

=== Problem Statement

Consider the boundary value problem on $[-1, 1]$:
$ u''(x) - (4 x^2 + 2) u(x) = 0, quad -1 lt.eq.slant x lt.eq.slant 1, $ <eq-bvp1>
with boundary conditions
$ u(-1) = 1, quad u(1) = 1. $

=== The Exact Solution

The exact solution is
$ u_"exact" (x) = e^(x^2 - 1). $ <eq-exact1>

Let us verify this claim. The first derivative is
$ u'_"exact" (x) = 2 x e^(x^2 - 1). $

The second derivative is
$ u''_"exact" (x) = (2 + 4 x^2) e^(x^2 - 1) = (4 x^2 + 2) u_"exact" (x). $

Substituting into the ODE:
$ u''_"exact" - (4 x^2 + 2) u_"exact" = (4 x^2 + 2) u_"exact" - (4 x^2 + 2) u_"exact" = 0. checkmark $

The boundary conditions are satisfied:
$ u_"exact" (plus.minus 1) = e^(1 - 1) = e^0 = 1. checkmark $

=== Trial Function

To satisfy the boundary conditions automatically, we write the approximation in a form that equals $1$ at $x = plus.minus 1$ regardless of the coefficient values. A convenient choice is:
$ u_2 (x) = 1 + (1 - x^2)(a_0 + a_1 x + a_2 x^2). $ <eq-trial1>

The factor $(1 - x^2)$ vanishes at the endpoints, so
$ u_2 (plus.minus 1) = 1 + 0 dot (dots.c) = 1 $
for any values of $a_0$, $a_1$, $a_2$. We have three undetermined coefficients.

Expanding the trial function#idx("trial function"):
$ u_2 (x) = 1 + a_0 + a_1 x + a_2 x^2 - a_0 x^2 - a_1 x^3 - a_2 x^4 $
$ = (1 + a_0) + a_1 x + (a_2 - a_0) x^2 - a_1 x^3 - a_2 x^4. $

=== Computing the Residual

The residual is
$ R(x; a_0, a_1, a_2) = u_(2)^('') (x) - (4 x^2 + 2) u_2 (x). $

Computing the second derivative of $u_2$:
$ u_(2)^(') (x) = a_1 + 2(a_2 - a_0) x - 3 a_1 x^2 - 4 a_2 x^3, $
$ u_(2)^('') (x) = 2(a_2 - a_0) - 6 a_1 x - 12 a_2 x^2. $

Substituting into the residual and simplifying (a calculation best verified with computer algebra), the residual is a polynomial of degree six in $x$ with coefficients that depend linearly on $a_0$, $a_1$, $a_2$:
$ R(x) = &(-2 - 4 a_0 + 2 a_2) - 8 a_1 x + (-4 - 2 a_0 - 14 a_2) x^2 - 2 a_1 x^3 \
&+ (4 a_0 - 6 a_2) x^4 + 4 a_1 x^5 + 4 a_2 x^6. $

=== Collocation Conditions

We have three unknowns, so we choose three collocation points in the interior of the interval. A simple choice is
$ x_1 = -1/2, quad x_2 = 0, quad x_3 = 1/2. $

Setting the residual to zero at these points gives three linear equations:

*At $x = -1\/2$*:
$ R(-1/2) = -17/4 a_0 + 33/8 a_1 - 25/16 a_2 - 3 = 0. $

*At $x = 0$*:
$ R(0) = -4 a_0 + 2 a_2 - 2 = 0. $

*At $x = 1\/2$*:
$ R(1/2) = -17/4 a_0 - 33/8 a_1 - 25/16 a_2 - 3 = 0. $

=== Solving the System

From the second equation:
$ -4 a_0 + 2 a_2 = 2 quad arrow.r.double quad a_2 = 1 + 2 a_0. $

Adding the first and third equations (the $a_1$ terms cancel):
$ -17/2 a_0 - 25/8 a_2 = 6. $

Substituting $a_2 = 1 + 2 a_0$:
$ -17/2 a_0 - 25/8 (1 + 2 a_0) = 6 $
$ -17/2 a_0 - 25/8 - 25/4 a_0 = 6 $
$ -(34/4 + 25/4) a_0 = 6 + 25/8 $
$ -59/4 a_0 = 73/8 $
$ a_0 = -73/118. $

Then
$ a_2 = 1 + 2 dot (-73/118) = 1 - 146/118 = -28/118 = -14/59. $

Subtracting the first equation from the third:
$ -33/4 a_1 = 0 quad arrow.r.double quad a_1 = 0. $

The vanishing of $a_1$ reflects the symmetry of the problem: both the differential equation and the boundary conditions are symmetric about $x = 0$, so the solution must be an even function. An odd coefficient like $a_1$ would break this symmetry.

=== The Approximate Solution

Substituting the coefficients back:
$ u_2 (x) = 1 + (1 - x^2) (-73/118 - 14/59 x^2). $

After simplification, this becomes the even polynomial:
$ u_2 (x) = 14/59 x^4 + 45/118 x^2 + 45/118. $ <eq-approx1>

The boundary conditions are satisfied:
$ u_2 (plus.minus 1) = 14/59 + 45/118 + 45/118 = 28/118 + 90/118 = 118/118 = 1. checkmark $

=== Error Analysis

The following table compares the exact and approximate solutions at several points:

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
          table.cell(fill: rgb("#142D6E").lighten(85%))[*$x$*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*$u_"exact" (x)$*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*$u_"approx" (x)$*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*Error*],
        ),
        table.hline(stroke: 0.5pt + luma(180)),
        num[-1], num[1.00000], num[1.00000], num[0.00000],
        num[-0.5], num[0.47237], num[0.49153], num[-0.01916],
        num[0], num[0.36788], num[0.38136], num[-0.01348],
        num[0.5], num[0.47237], num[0.49153], num[-0.01916],
        num[1], num[1.00000], num[1.00000], num[0.00000],
      )
    },
  ),
  caption: [Comparison of exact and three-coefficient collocation approximation.],
) <tab-error1>

The maximum error is approximately $2 times 10^(-2)$, which is remarkably good for such a low-order approximation. @fig-collocation-example1 shows the solutions graphically.

#figure(
  image("../figures/ch03/python/collocation_example1.pdf", width: 95%),
  caption: [Left: exact solution $u(x) = e^(x^2 - 1)$ compared with the three-coefficient collocation approximation. The collocation points $x = -1\/2, 0, 1\/2$ are marked with squares. Right: the error $u_"exact" - u_"approx"$.],
) <fig-collocation-example1>

The code generating @fig-collocation-example1 is available in:
- `codes/python/ch03/collocation_example1.py`
- `codes/matlab/ch03/collocation_example1.m`
- `codes/julia/ch03/collocation_example1.jl`

#etude-conclusion[
  With just three free coefficients, the collocation approximation reproduces $u(x) = e^(x^2 - 1)$ to within $2 times 10^(-2)$ across the entire interval. The error plot reveals a smooth, oscillatory pattern with zeros precisely at the three collocation points --- exactly as the method guarantees --- and symmetric peaks between them. The vanishing of $a_1$ is a direct consequence of the even symmetry of the differential equation and boundary conditions; symmetry is doing free work for us. The key lesson is that the *method of weighted residuals#idx("method of weighted residuals")*, even in its simplest collocation form, converts a differential equation into a small algebraic system that can be solved by hand. The trial-function strategy of building boundary conditions into the basis --- the factor $(1 - x^2)$ that vanishes at the endpoints --- is a recurring theme in spectral methods. In later chapters we will replace this ad-hoc polynomial basis with Chebyshev expansions and the hand-chosen points with optimal Gauss--Lobatto distributions, but the core principle on display here will remain unchanged.
]

== Computational Étude 3.2: Collocation versus Galerkin

To compare the two main approaches to the Method of Weighted Residuals, a comparison famously analyzed by @Villadsen1967, we consider a second example where both methods can be applied with explicit hand calculations.

=== Problem Statement

Consider the reaction-diffusion equation on $[-1, 1]$:
$ (dif^2 u) / (dif x^2) - 4 u = -1, quad -1 lt.eq.slant x lt.eq.slant 1, $ <eq-bvp2>
with homogeneous Dirichlet boundary conditions:
$ u(-1) = 0, quad u(1) = 0. $

=== The Exact Solution

The homogeneous equation $u'' - 4u = 0$ has the general solution $u_h = A cosh(2x) + B sinh(2x)$. A particular solution for the constant forcing $-1$ is $u_p = 1\/4$. The general solution is therefore
$ u(x) = A cosh(2x) + B sinh(2x) + 1/4. $

The boundary condition $u(-1) = 0$ gives $A cosh(2) - B sinh(2) + 1\/4 = 0$.
The boundary condition $u(1) = 0$ gives $A cosh(2) + B sinh(2) + 1\/4 = 0$.

Adding these equations: $2 A cosh(2) + 1\/2 = 0$, so $A = -1\/(4 cosh(2))$.
Subtracting: $2 B sinh(2) = 0$, so $B = 0$.

The exact solution is:
$ u_"exact" (x) = 1/4 (1 - (cosh(2x)) / (cosh(2))). $ <eq-exact2>

The maximum value occurs at $x = 0$:
$ u_"exact" (0) = 1/4 (1 - 1 / cosh(2)) approx 0.1835. $

=== Trial Function and Basis

Since the boundary conditions are homogeneous, we choose basis functions that automatically vanish at $x = plus.minus 1$. For a symmetric problem like this one, we use even functions:
$ phi_0 (x) = 1 - x^2, quad phi_1 (x) = (1 - x^2) x^2 = x^2 - x^4. $

Both functions vanish at $x = plus.minus 1$ and are even in $x$. Our trial function is:
$ u_1 (x) = a_0 phi_0 (x) + a_1 phi_1 (x) = a_0 (1 - x^2) + a_1 (x^2 - x^4). $ <eq-trial2>

=== The Residual

The operator is $cal(L) u = u'' - 4u$. We compute:
$ phi_(0)^('') = -2, quad phi_(1)^('') = 2 - 12 x^2. $

Applying the operator to each basis function:
$ cal(L) phi_0 = -2 - 4(1 - x^2) = -6 + 4 x^2, $
$ cal(L) phi_1 = (2 - 12 x^2) - 4(x^2 - x^4) = 2 - 16 x^2 + 4 x^4. $

The residual is:
$ R(x) = a_0 cal(L) phi_0 + a_1 cal(L) phi_1 - (-1) $
$ = a_0 (-6 + 4 x^2) + a_1 (2 - 16 x^2 + 4 x^4) + 1. $

=== Collocation Method

With two unknowns, we need two collocation points. Due to symmetry, we can use points in $[0, 1)$. We choose $x_1 = 0$ and $x_2 = 0.5$.

*At $x = 0$*:
$ cal(L) phi_0 (0) = -6, quad cal(L) phi_1 (0) = 2. $
$ R(0) = -6 a_0 + 2 a_1 + 1 = 0 quad arrow.r.double quad 6 a_0 - 2 a_1 = 1. $

*At $x = 0.5$*:
$ cal(L) phi_0 (0.5) = -6 + 4 dot 0.25 = -5, $
$ cal(L) phi_1 (0.5) = 2 - 16 dot 0.25 + 4 dot 0.0625 = 2 - 4 + 0.25 = -1.75. $
$ R(0.5) = -5 a_0 - 1.75 a_1 + 1 = 0 quad arrow.r.double quad 5 a_0 + 1.75 a_1 = 1. $

Solving the system:
$ cases(6 a_0 - 2 a_1 = 1, 5 a_0 + 1.75 a_1 = 1) $

Multiply the first equation by $5$ and the second by $6$:
$ 30 a_0 - 10 a_1 = 5, quad 30 a_0 + 10.5 a_1 = 6. $

Subtracting: $20.5 a_1 = 1$, so $a_1 approx 0.04878$.

Substituting back: $6 a_0 = 1 + 2 dot 0.04878 = 1.09756$, so $a_0 approx 0.1829$.

The collocation solution is:
$ u_"coll" (x) approx 0.1829 (1 - x^2) + 0.0488 (x^2 - x^4). $

=== Galerkin Method

The Galerkin conditions require the residual to be orthogonal to each basis function:
$ integral_(-1)^1 R(x) phi_k (x) dif x = 0, quad k = 0, 1. $

This gives a symmetric matrix system $bold(A) bold(a) = bold(b)$ where:
$ A_(i j) = integral_(-1)^1 cal(L) phi_j (x) dot phi_i (x) dif x, $
$ b_i = integral_(-1)^1 (-1) dot phi_i (x) dif x = -integral_(-1)^1 phi_i (x) dif x. $

Computing the integrals (using standard formulas for powers of $x$):

For $phi_0 = 1 - x^2$:
$ integral_(-1)^1 phi_0 (x) dif x = integral_(-1)^1 (1 - x^2) dif x = [x - x^3 / 3]_(-1)^1 = 2 - 2/3 = 4/3. $

For $phi_1 = x^2 - x^4$:
$ integral_(-1)^1 phi_1 (x) dif x = [x^3 / 3 - x^5 / 5]_(-1)^1 = 2/3 - 2/5 = 4/15. $

The right-hand side is:
$ b_0 = -4/3, quad b_1 = -4/15. $

The matrix entries require more computation. Using computer algebra or careful integration:
$ A_(0 0) = -104/15, quad A_(0 1) = A_(1 0) = -8/7, quad A_(1 1) = -328/315. $

Solving the system yields:
$ a_0 approx 0.1832, quad a_1 approx 0.0550. $

The Galerkin method requires numerical integration to assemble the system.

=== Comparison

The following table compares the two methods at the central point $x = 0$:

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    {
      show table: format-table(none, auto, auto)
      table(
        columns: 3,
        align: (left, center, center),
        inset: (x: 1em, y: 0.6em),
        stroke: none,
        table.hline(stroke: 0.75pt + rgb("#142D6E")),
        table.header(
          table.cell(fill: rgb("#142D6E").lighten(85%))[*Method*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*$u(0)$*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*Absolute Error*],
        ),
        table.hline(stroke: 0.5pt + luma(180)),
        [Exact], num[0.1835], num[0],
        [Collocation ($N = 2$)], num[0.1829], num[0.0006],
        [Galerkin ($N = 2$)], num[0.1832], num[0.0003],
      )
    },
  ),
  caption: [Comparison of spectral approximation#idx("spectral approximation")s at the central maximum.],
) <tab-comparison>

For this problem, the Galerkin method is more accurate both at the central point and in a global sense. This is consistent with the error plot in @fig-collocation-vs-galerkin, which shows the Galerkin error (green) remaining closer to zero across the entire interval. The Galerkin method minimizes the error in a root-mean-square sense, which typically leads to better overall accuracy for smooth problems.

@fig-collocation-vs-galerkin shows both approximate solutions compared to the exact solution, along with the error profiles.

#figure(
  image("../figures/ch03/python/collocation_vs_galerkin.pdf", width: 95%),
  caption: [Left: exact solution compared with collocation and Galerkin approximations. The collocation points $x = 0$ and $x = 0.5$ are marked. Right: error profiles for both methods.],
) <fig-collocation-vs-galerkin>

The code generating @fig-collocation-vs-galerkin is available in:
- `codes/python/ch03/collocation_vs_galerkin.py`
- `codes/matlab/ch03/collocation_vs_galerkin.m`
- `codes/julia/ch03/collocation_vs_galerkin.jl`

#etude-conclusion[
  Both methods achieve remarkable accuracy with only two free coefficients: collocation misses the exact maximum at $x = 0$ by $6 times 10^(-4)$, Galerkin by $3 times 10^(-4)$. The error profiles confirm that Galerkin provides *uniformly* smaller errors across the entire interval --- a consequence of its global minimisation principle, which enforces orthogonality of the residual against the basis functions rather than demanding zero residual at a few discrete points. The practical trade-off, however, is instructive: collocation needs only pointwise evaluations of the differential operator, while Galerkin demands inner-product integrals (here evaluated analytically, but in general requiring numerical quadrature). This asymmetry --- *collocation simpler, Galerkin more accurate* --- is a persistent theme throughout spectral methods. For the nonlinear and high-dimensional problems of later chapters, *pseudospectral collocation will dominate* precisely because it avoids the integral-assembly step, and the accuracy gap narrows rapidly as $N$ increases.
]

== Conclusions and Questions

These simple examples illuminate several important features of spectral methods:

+ *Automatic satisfaction of boundary conditions*. By choosing trial functions that vanish at the boundaries (or equal the prescribed boundary values), we eliminate the boundary conditions from the algebraic system.

+ *Symmetry exploitation*. When the problem has symmetry, the solution inherits that symmetry. In our examples, the vanishing of odd coefficients ($a_1 = 0$) reflects the even symmetry of the exact solution.

+ *Simplicity of implementation*. Even with hand calculations, we can obtain remarkably accurate approximations with just a few coefficients.

+ *Trade-offs between methods*. Collocation is simpler to implement (no integrals to compute), while Galerkin typically provides better accuracy by minimizing a global error measure.

These examples raise important questions that we will address in subsequent chapters:

- *What is the optimal choice of basis functions?* Using simple powers of $x$ works for small $N$, but becomes numerically unstable for large $N$. Chebyshev and Fourier bases are far superior.

- *What are the optimal collocation points?* Our ad hoc choices $x = -1\/2, 0, 1\/2$ worked well, but there exist optimal point distributions (Gauss and Gauss--Lobatto points) derived from orthogonal polynomial theory.

- *How fast does the error decrease as $N$ increases?* For smooth solutions, spectral methods achieve exponential convergence (the error decreases like $c^(-N)$ for some $c > 1$), which is dramatically faster than the algebraic convergence $O(N^(-p))$ of finite difference and finite element methods.

The following chapters will develop the theory and algorithms needed to answer these questions and to apply spectral methods to a wide range of problems.

== A Broader Perspective

For demanding students who wish to understand how spectral methods fit into the wider landscape of numerical analysis, this section provides a high-level comparison with other approaches, particularly finite element methods. The discussion follows @Boyd2000[Section 1.3].

=== Local versus Global Basis Functions

The fundamental distinction between spectral methods and finite element methods lies in the _support_ of the basis functions. In finite element methods, the computational domain is divided into many small sub-intervals (or triangles, tetrahedra in higher dimensions), and the basis functions $phi_n (x)$ are _local_: they are polynomials of fixed, low degree (typically linear or quadratic) that are non-zero only over one or two adjacent elements.

In contrast, spectral methods use _global_ basis functions. Each $phi_n (x)$ is a polynomial (or trigonometric polynomial) of potentially high degree that is non-zero (except at isolated points) over the entire computational domain. This global character is both the source of spectral methods' power and the reason for some of their limitations.

@fig-local-vs-global illustrates this distinction schematically. The finite element basis function (left) has compact support and contributes to the solution only locally. The spectral basis function (right) influences the solution everywhere.

#figure(
  grid(
    columns: 2,
    gutter: 2em,
    // Left: Local (FEM) basis function
    align(center)[
      #box(
        width: 3.8cm,
        height: 2.5cm,
        stroke: 0.5pt + luma(180),
        radius: 3pt,
        inset: 0.3em,
      )[
        #place(bottom + left, dy: -0.3em, dx: 0.3em)[
          #line(start: (0pt, 0pt), end: (3.2cm, 0pt), stroke: 0.8pt + luma(100))
        ]
        #place(bottom + left, dy: -0.3em, dx: 0.3em)[
          // Element divisions
          #for i in range(5) {
            place(dx: i * 0.8cm, line(start: (0pt, 0pt), end: (0pt, 3pt), stroke: 0.5pt + luma(100)))
          }
        ]
        #place(bottom + left, dy: -0.5em, dx: 0.3em)[
          // Hat function (piecewise linear)
          #curve(
            stroke: 1.5pt + rgb("#142D6E"),
            curve.move((0cm, 0pt)),
            curve.line((0.8cm, 0pt)),
            curve.line((1.6cm, 1.5cm)),
            curve.line((2.4cm, 0pt)),
            curve.line((3.2cm, 0pt)),
          )
        ]
        #place(top + center, dy: 0.2em)[
          #text(size: 8pt, fill: luma(80))[Local support]
        ]
      ]
      #v(0.3em)
      #text(size: 9pt)[Finite element basis]
    ],
    // Right: Global (spectral) basis function
    align(center)[
      #box(
        width: 3.8cm,
        height: 2.5cm,
        stroke: 0.5pt + luma(180),
        radius: 3pt,
        inset: 0.3em,
      )[
        #place(bottom + left, dy: -0.3em, dx: 0.3em)[
          #line(start: (0pt, 0pt), end: (3.2cm, 0pt), stroke: 0.8pt + luma(100))
        ]
        #place(bottom + left, dy: -0.5em, dx: 0.3em)[
          // Chebyshev polynomial T_5, sampled densely so it renders as a smooth
          // global oscillation rather than a piecewise-linear sawtooth.
          #let npts = 60
          #curve(
            stroke: 1.5pt + rgb("#142D6E"),
            curve.move((0cm, (0.9 + 0.55 * calc.cos(5 * calc.acos(-1.0))) * 1cm)),
            ..range(1, npts + 1).map(i => {
              let t = -1.0 + 2.0 * i / npts
              let y = 0.9 + 0.55 * calc.cos(5 * calc.acos(t))
              curve.line((i / npts * 3.2 * 1cm, y * 1cm))
            }),
          )
        ]
        #place(top + center, dy: 0.2em)[
          #text(size: 8pt, fill: luma(80))[Global support]
        ]
      ]
      #v(0.3em)
      #text(size: 9pt)[Spectral basis]
    ],
  ),
  caption: [Schematic comparison of local (finite element) and global (spectral) basis functions. The finite element "hat function" is non-zero only over two adjacent elements, while the spectral basis function oscillates smoothly across the entire domain.],
) <fig-local-vs-global>

=== Refinement Strategies

When a numerical approximation is insufficiently accurate, there are three fundamentally different strategies to improve it, illustrated schematically in @fig-refinement-strategies:

+ *$h$-refinement*: Subdivide each element into smaller pieces, reducing the mesh spacing $h$ uniformly throughout the domain. This increases the number of elements while keeping the polynomial degree fixed.

+ *$r$-refinement* (adaptive): Redistribute the mesh points, clustering them in regions where the solution has steep gradients or other features requiring high resolution. The total number of degrees of freedom remains roughly constant.

+ *$p$-refinement*: Keep the mesh fixed while increasing $p$, the polynomial degree within each element. For a single-element domain, this is precisely what spectral methods do. This approach is mathematically equivalent to the $p$-version of the finite element method established by @Babuska1981.

#figure(
  grid(
    columns: 3,
    gutter: 1.5em,
    // h-refinement
    align(center)[
      #stack(
        dir: ttb,
        spacing: 0.5em,
        // Original mesh
        box(
          width: 2.8cm,
          height: 1.2cm,
          stroke: 0.8pt + luma(150),
        )[
          #place(center + horizon)[
            #grid(columns: 2, gutter: 0pt,
              box(width: 1.4cm, height: 1.2cm, stroke: (right: 0.5pt + luma(150))),
              box(width: 1.4cm, height: 1.2cm),
            )
          ]
        ],
        // Arrow
        text(size: 14pt)[↓],
        // Refined mesh
        box(
          width: 2.8cm,
          height: 1.2cm,
          stroke: 0.8pt + luma(150),
        )[
          #place(center + horizon)[
            #grid(columns: 4, gutter: 0pt,
              box(width: 0.7cm, height: 1.2cm, stroke: (right: 0.5pt + luma(150))),
              box(width: 0.7cm, height: 1.2cm, stroke: (right: 0.5pt + luma(150))),
              box(width: 0.7cm, height: 1.2cm, stroke: (right: 0.5pt + luma(150))),
              box(width: 0.7cm, height: 1.2cm),
            )
          ]
        ],
      )
      #v(0.3em)
      #text(size: 9pt, style: "italic")[$h$-refinement] \
      #text(size: 8pt, fill: luma(100))[Smaller elements]
    ],
    // r-refinement
    align(center)[
      #stack(
        dir: ttb,
        spacing: 0.5em,
        // Original mesh
        box(
          width: 2.8cm,
          height: 1.2cm,
          stroke: 0.8pt + luma(150),
        )[
          #place(center + horizon)[
            #grid(columns: 4, gutter: 0pt,
              box(width: 0.7cm, height: 1.2cm, stroke: (right: 0.5pt + luma(150))),
              box(width: 0.7cm, height: 1.2cm, stroke: (right: 0.5pt + luma(150))),
              box(width: 0.7cm, height: 1.2cm, stroke: (right: 0.5pt + luma(150))),
              box(width: 0.7cm, height: 1.2cm),
            )
          ]
        ],
        // Arrow
        text(size: 14pt)[↓],
        // Adapted mesh
        box(
          width: 2.8cm,
          height: 1.2cm,
          stroke: 0.8pt + luma(150),
        )[
          #place(center + horizon)[
            #grid(columns: 4, gutter: 0pt,
              box(width: 1.1cm, height: 1.2cm, stroke: (right: 0.5pt + luma(150))),
              box(width: 0.35cm, height: 1.2cm, stroke: (right: 0.5pt + luma(150))),
              box(width: 0.35cm, height: 1.2cm, stroke: (right: 0.5pt + luma(150))),
              box(width: 1.0cm, height: 1.2cm),
            )
          ]
        ],
      )
      #v(0.3em)
      #text(size: 9pt, style: "italic")[$r$-refinement] \
      #text(size: 8pt, fill: luma(100))[Adaptive redistribution]
    ],
    // p-refinement
    align(center)[
      #stack(
        dir: ttb,
        spacing: 0.5em,
        // Low degree
        box(
          width: 2.8cm,
          height: 1.2cm,
          stroke: 0.8pt + luma(150),
        )[
          #place(bottom + left, dx: 0.1cm, dy: -0.2cm)[
            #curve(
              stroke: 1.2pt + rgb("#142D6E"),
              // Degree-2 polynomial: a single smooth arch.
              curve.move((0cm, 0.3cm)),
              ..range(1, 25).map(i => {
                let t = -1.0 + 2.0 * i / 24
                curve.line((i / 24 * 2.6 * 1cm, (0.7 - 0.4 * t * t) * 1cm))
              }),
            )
          ]
          #place(top + right, dx: -0.15cm, dy: 0.15cm)[
            #text(size: 7pt, fill: luma(100))[$p = 2$]
          ]
        ],
        // Arrow
        text(size: 14pt)[↓],
        // High degree
        box(
          width: 2.8cm,
          height: 1.2cm,
          stroke: 0.8pt + luma(150),
        )[
          #place(bottom + left, dx: 0.1cm, dy: -0.2cm)[
            #curve(
              stroke: 1.2pt + rgb("#142D6E"),
              // Degree-8 polynomial (T_8), sampled densely so it reads as a
              // smooth high-degree oscillation rather than a sawtooth.
              curve.move((0cm, 0.8cm)),
              ..range(1, 49).map(i => {
                let t = -1.0 + 2.0 * i / 48
                curve.line((i / 48 * 2.6 * 1cm, (0.5 + 0.3 * calc.cos(8 * calc.acos(t))) * 1cm))
              }),
            )
          ]
          #place(top + right, dx: -0.15cm, dy: 0.15cm)[
            #text(size: 7pt, fill: luma(100))[$p = 8$]
          ]
        ],
      )
      #v(0.3em)
      #text(size: 9pt, style: "italic")[$p$-refinement] \
      #text(size: 8pt, fill: luma(100))[Higher polynomial degree]
    ],
  ),
  caption: [Three strategies for improving accuracy in numerical methods. Spectral methods employ $p$-refinement: increasing the polynomial degree while using a single element (or few elements) spanning the entire domain.],
) <fig-refinement-strategies>

Spectral methods can be viewed as the extreme form of $p$-refinement: a single element spans the entire domain, and accuracy is improved solely by increasing the polynomial degree. This strategy is devastatingly effective when the solution is smooth, but struggles when the solution has discontinuities or sharp gradients.

=== Trade-offs: Sparse versus Full Matrices

The choice between local and global basis functions entails fundamental trade-offs:

*Finite element advantages:*
- _Sparse matrices_: Since each basis function is non-zero over only a few elements, the stiffness matrix has mostly zero entries. Sparse matrix solvers can exploit this structure, dramatically reducing computational cost for large systems.
- _Geometric flexibility_: The small elements (triangles, tetrahedra) can be fitted to irregularly shaped domains like automobile bodies or aircraft wings.

*Finite element disadvantages:*
- _Low accuracy per degree of freedom_: Each basis function is a polynomial of low degree, so many elements are needed for high accuracy.

*Spectral method advantages:*
- _High accuracy for smooth problems_: The high-degree global polynomials capture smooth solutions with far fewer degrees of freedom.
- _Efficiency with iterative solvers_: When fast iterative methods are used, spectral methods can be much more efficient than low-order methods for many problem classes.

*Spectral method disadvantages:*
- _Full matrices_: The global basis functions create dense matrices where most entries are non-zero.
- _Geometric limitations_: Spectral methods are most natural on simple domains (intervals, rectangles, disks) and require more sophisticated techniques for irregular geometries.

For problems with smooth solutions on regular domains (many important problems in fluid dynamics, quantum mechanics, and wave propagation fall into this category), the accuracy advantage of spectral methods often outweighs the matrix structure disadvantage.

=== Spectral Element Methods

A natural question arises: can we combine the geometric flexibility of finite elements with the high accuracy of spectral methods? The answer is yes, through _spectral element methods_, as introduced by @Patera1984.

In spectral element methods, the domain is subdivided into elements (as in finite elements), but within each element, the polynomial degree $p$ is chosen to be moderately high, typically $p = 6$ to $8$. This hybrid approach inherits several advantages:

- The element subdivision provides geometric flexibility and matrix sparsity.
- The high polynomial degree within each element provides spectral-like accuracy.
- The theoretical framework is essentially the same as for global spectral methods.

Spectral element codes are typically written so that $p$ is a user-adjustable parameter, allowing practitioners to balance accuracy and cost for their specific application. We will not develop spectral element methods in detail in this book, but the reader should be aware that much of the theory developed for global spectral methods transfers directly to the spectral element context.

=== Radial Basis Functions and RBF-FD

A parallel line of development, largely absent from the classical spectral methods literature but increasingly prominent in modern computational science, is the use of Radial Basis Functions (RBFs). In this approach, the unknown function is approximated as a linear combination of translates $phi(||bold(x) - bold(x)_j||)$ of a single radially symmetric kernel (such as the Gaussian or multiquadric), centered at the data points. When the nodes are placed on a grid and the kernel is taken to the "flat" limit, the RBF interpolant reduces exactly to the polynomial pseudospectral interpolant studied in this book; in this sense, pseudospectral methods are a special case of the RBF framework. The decisive advantage of RBFs is their ability to achieve spectral accuracy on _scattered_ nodes in _any number of dimensions_, without requiring tensor-product grids or simple geometries. The RBF-FD (radial basis function--generated finite differences) method carries this further by evaluating local RBF interpolants on small node stencils, thereby recovering the sparsity of finite difference methods while retaining much higher accuracy than traditional FD stencils. The evolution FD $arrow.r$ PS $arrow.r$ RBF $arrow.r$ RBF-FD represents a unifying intellectual arc, developed at length in @Fornberg2025. We do not pursue RBF methods in this book, but the reader should be aware that they offer a natural path from the global polynomial methods developed here toward flexible, high-order computation on complex domains.

=== The Convergence of Methods at High Order

Perhaps the most profound insight from the comparison between finite element and spectral methods is this: _for sufficiently high polynomial degree, the two approaches become essentially equivalent_.

Low-order finite elements (linear, quadratic) can be derived, justified, and implemented without knowledge of Fourier or Chebyshev convergence theory. However, as the polynomial degree increases, ad hoc approaches become increasingly ill-conditioned and numerically unstable, a phenomenon rigorously analyzed in @GottliebOrszag1977. The only practical way to implement well-behaved high-order finite elements (say, sixth order or higher) is to use the technology of spectral methods: Chebyshev or Legendre basis functions, Gaussian quadrature, and the convergence theory we will develop in subsequent chapters.

Thus, the question "Are finite elements or spectral methods better?" becomes somewhat artificial for high-order approximations. The real question is: _Does the problem at hand require high-order accuracy, or is second or fourth order sufficient?_

When the solution is smooth and high accuracy is needed, the spectral/high-order approach is clearly superior. When the solution has discontinuities, shocks, or boundary layers, or when the geometry is highly irregular, low-order methods with adaptive mesh refinement may be more practical. The wise practitioner chooses the tool appropriate to the problem.

== A non-exhaustive literature overview

The Method of Weighted Residuals (MWR), which forms the backbone of the techniques explored in this chapter, is not a single invention but a synthesis of ideas that spanned the early 20th century. The unifying framework was formalized in the landmark review by @Finlayson1966, who demonstrated that diverse approximation schemes---including the method of moments, the subdomain method, and the Galerkin method---could all be rigorously classified by their choice of weight function. This taxonomy transformed numerical analysis from a collection of ad hoc recipes into a coherent mathematical discipline.

The specific tension between *collocation* and *Galerkin* methods discussed in this chapter was a central debate in the computational community for decades. While Galerkin methods offered theoretical optimality in energy norms, their implementation was hampered by the high cost of numerical integration. This bottleneck was resolved by the seminal work of @Villadsen1967 in chemical engineering. They introduced *orthogonal collocation*, proving that by selecting collocation points as the roots of orthogonal polynomials, one could achieve the accuracy of Galerkin methods with the computational efficiency of point-wise evaluation. This insight is the direct mathematical ancestor of the pseudospectral methods we utilize in this text.

The rigorous analysis of these methods was cemented in the 1970s. @Orszag1971 and @GottliebOrszag1977 established the error estimates for "pseudospectral" methods, demonstrating that the aliasing errors introduced by collocation are generally of the same order as the truncation errors, thus validating the use of the Fast Fourier Transform for nonlinear problems in fluid dynamics.

The "broader perspective" of local versus global bases mirrors the historical development of the *Finite Element Method* (FEM). While classical FEM relies on mesh refinement ($h$-version), @Babuska1981 pioneered the $p$-version of the finite element method, proving that increasing the polynomial degree $p$ on a fixed mesh yields exponential convergence for smooth solutions. This philosophy evolved into the *Spectral Element Method* (SEM) introduced by @Patera1984, which combines the geometric flexibility of finite elements with the high-order accuracy of spectral expansions---a synthesis that remains the gold standard for high-fidelity simulation in complex geometries.

In the contemporary era (2020--2026), spectral methods are undergoing a renaissance at the intersection of machine learning and nonlocal physics. @Meuris2023 and @Ngueabou2025 review how Deep Operator Networks (DeepONet) and Physics-Informed Neural Networks (PINNs) are being used to "learn" optimal spectral bases, overcoming the curse of dimensionality in high-dimensional PDEs. Furthermore, the extension of spectral methods to *fractional differential equations*, as detailed in the recent monograph by @Zayernouri2024, and the development of spectral solvers for *nonlocal diffusion on bounded domains* @Mustapha2025, illustrate the enduring adaptability of the spectral approach to the frontiers of mathematical physics.

== Summary

This chapter has provided a first taste of spectral methods through low-dimensional, hand-computable examples:

+ *Method of Weighted Residuals*: Spectral methods approximate the solution as a finite sum of basis functions and determine the unknown coefficients by minimizing the residual --- either pointwise (collocation) or in a weighted integral sense (Galerkin).

+ *Collocation simplicity*: Forcing the residual to vanish at selected points converts a differential equation into a small algebraic system, with no integrals to evaluate.

+ *Galerkin accuracy*: Requiring orthogonality of the residual against basis functions typically yields better global accuracy, at the cost of computing integrals.

+ *Boundary conditions*: Choosing trial functions that automatically satisfy the boundary data eliminates boundary conditions from the algebraic system.

+ *Symmetry exploitation*: When the problem possesses symmetry, the solution inherits that symmetry, and odd coefficients vanish automatically --- reducing the effective size of the algebraic system.

+ *Local versus global bases*: Finite element methods use low-degree, compactly supported basis functions on many small elements, producing sparse matrices. Spectral methods use high-degree global polynomials on the entire domain, producing dense but small matrices. For smooth solutions, the spectral approach requires far fewer degrees of freedom.

The examples in this chapter were intentionally simple. The chapters that follow will replace ad hoc polynomial bases with Chebyshev and Fourier expansions, and hand-chosen collocation points with optimal distributions derived from orthogonal polynomial theory.

== Exercises <sec-mise-en-bouche-exercises>

The exercises below progress from pencil-and-paper properties of the Method of Weighted Residuals, through numerical experiments that reproduce and extend the two études of this chapter, to open-ended projects that reach into the current research literature. The computational problems may be carried out in any of the book's three languages; the named scripts under `codes/` give a starting point.

=== Conceptual Exercises

#exercise(title: [The Residual is Affine in the Coefficients])[
  Let $cal(L)$ be a linear differential operator and let $u_N$ be the truncated expansion @eq-series-expansion. (a) Show that the residual @eq-residual is an affine function of the coefficient vector $bold(a) = (a_0, dots, a_N)^top$, namely $R(x; bold(a)) = sum_(n=0)^N a_n cal(L) phi_n (x) - f(x)$. (b) Deduce that both the collocation conditions and the Galerkin conditions yield a linear system $bold(A) bold(a) = bold(b)$, and give explicit formulas for $bold(A)$ and $bold(b)$ in each case. (c) Identify precisely where the linearity of $cal(L)$ is used, and describe what changes when $cal(L)$ is nonlinear.
] <ex-meb-residual-affine>

#exercise(title: [The Weight-Function Taxonomy])[
  @Finlayson1966 unified the Method of Weighted Residuals by classifying methods according to the weight functions ${w_k}$ in the orthogonality conditions $integral_(-1)^1 R(x) thin w_k (x) dif x = 0$. Starting from this single template, derive the weights that reproduce: (a) the collocation method, as point evaluation at nodes $x_k$; (b) the subdomain method#idx("subdomain method"), in which the residual integrates to zero over each of $N + 1$ subintervals; (c) the least-squares method, which minimises $integral_(-1)^1 R^2 dif x$ over the coefficients; and (d) the Galerkin method. State $w_k$ explicitly in each case, and for least squares show that $w_k = partial R \/ partial a_k$.
] <ex-meb-weight-taxonomy>

#exercise(title: [Collocation as Interpolation of the Residual])[
  For the boundary-value problem @eq-bvp1 with the trial function @eq-trial1, the residual @eq-residual is a polynomial of degree six in $x$, yet only three collocation conditions are imposed. (a) Show that, for a linear problem discretised with $N + 1$ basis functions and $N + 1$ collocation points, the collocation conditions are equivalent to requiring that the polynomial interpolant of the residual through the collocation points vanish identically. (b) Explain why the residual itself need not vanish between the collocation points, and relate its size there to the interpolation error of the residual. (c) Argue that this is why the placement of the collocation points matters even when their number is fixed.
] <ex-meb-collocation-interpolation>

#hint-for(<ex-meb-collocation-interpolation>)[The residual minus its interpolant through the $N + 1$ points is a polynomial vanishing at those points, so it is controlled by the standard interpolation-error estimate; driving the residual down between nodes is therefore a question of where the nodes are placed.]

#exercise(title: [Galerkin as Energy Minimisation])[
  Write @eq-bvp2 in its self-adjoint, positive-definite form $-u'' + 4 u = 1$ on $[-1, 1]$ with $u(plus.minus 1) = 0$, and set $cal(M) u = -u'' + 4 u$. (a) Using integration by parts, show that $cal(M)$ is symmetric and positive definite for the inner product $chevron.l f, g chevron.r = integral_(-1)^1 f(x) g(x) dif x$ on functions vanishing at $x = plus.minus 1$. (b) Define the energy $J[v] = 1/2 chevron.l cal(M) v, v chevron.r - chevron.l 1, v chevron.r$ and show that, restricted to a trial space such as @eq-trial2, the stationarity conditions $partial J \/ partial a_k = 0$ are exactly the Galerkin equations; conclude that the Galerkin coefficients minimise $J$ over the trial space and that the Galerkin solution is the best approximation to the exact solution in the energy norm $||v||_cal(M) = chevron.l cal(M) v, v chevron.r^(1\/2)$. (c) Use this best-approximation property to explain why the Galerkin error in @fig-collocation-vs-galerkin lies uniformly below the collocation error.
] <ex-meb-galerkin-energy>

#hint-for(<ex-meb-galerkin-energy>)[On the trial space $J$ is a quadratic in $bold(a)$ whose Hessian is the Galerkin matrix and whose gradient set to zero is the Galerkin system. For the best-approximation claim, write the error as exact minus Galerkin and use the $cal(M)$-orthogonality of the Galerkin residual to the trial space.]

#exercise(title: [Symmetry Forces Even Solutions])[
  In @eq-bvp1 both the operator and the boundary conditions are invariant under the reflection $x arrow.r.bar -x$, and the collocation solution @eq-approx1 turned out to be even, with the odd coefficient $a_1 = 0$. Make this precise. (a) Writing $(cal(R) u)(x) = u(-x)$, show that if $cal(L)$ commutes with $cal(R)$, the forcing $f$ is even, and the boundary conditions are symmetric, then the exact solution is even. (b) Prove that if, in addition, the trial basis in @eq-series-expansion and the set of collocation points are symmetric about the origin, then the collocation system forces every odd coefficient to vanish. (c) Explain how this observation halves the effective number of unknowns.
] <ex-meb-symmetry-even>

#exercise(title: [Lifting Inhomogeneous Boundary Conditions])[
  The trial function @eq-trial1 builds the boundary data into the basis through the factor $(1 - x^2)$, which vanishes at $x = plus.minus 1$. Generalise this device. (a) For Dirichlet data $u(-1) = alpha$ and $u(1) = beta$, construct an affine boundary-lifting function $B(x)$ with $B(-1) = alpha$ and $B(1) = beta$, and write $u_N (x) = B(x) + sum_(n=0)^N a_n (1 - x^2) psi_n (x)$ so that the boundary conditions hold for every choice of coefficients. (b) Show that the residual @eq-residual then depends on the data only through $cal(L) B$, so that inhomogeneous boundary conditions act as an additional forcing term. (c) Adapt the construction to a Neumann condition $u' (1) = gamma$.
] <ex-meb-bc-lifting>

#exercise(title: [Method of Weighted Residuals with Legendre Polynomials])[
  Consider $-u'' = 1$ on $[-1, 1]$ with $u(-1) = u(1) = 0$. The exact solution is $u(x) = (1 - x^2) \/ 2$. (a) Expand the solution as $u_N (x) = sum_(k=0)^N a_k P_k (x)$ where $P_k$ are Legendre polynomials and enforce the boundary conditions on the coefficients. (b) Apply Galerkin projection using the Legendre inner product $chevron.l f, g chevron.r = integral_(-1)^1 f(x) g(x) dif x$. (c) Show that the Galerkin solution with $N = 2$ recovers the exact solution. Why?
] <ex-meb-legendre-galerkin>

#exercise(title: [The Monomial Basis and the Hilbert Matrix])[
  The chapter cautions that the monomial basis $phi_n (x) = x^n$ becomes numerically unstable as $N$ grows. Quantify this on $[0, 1]$. (a) Show that the Galerkin mass matrix with entries $M_(i j) = integral_0^1 x^i x^j dif x$ is the Hilbert matrix#idx("Hilbert matrix") $M_(i j) = 1 \/ (i + j + 1)$. (b) Look up or estimate how the condition number of the $N times N$ Hilbert matrix grows with $N$, and explain what this implies for solving the coefficient system in floating-point arithmetic. (c) Explain qualitatively why an orthogonal basis (Legendre or Chebyshev) replaces this matrix by one that is diagonal or well-conditioned, motivating the developments of later chapters.
] <ex-meb-hilbert-matrix>

=== Computational Exercises

#exercise(title: [Refining the First Étude])[
  Reproduce and extend the first collocation example, the problem @eq-bvp1 with exact solution @eq-exact1, using the script `collocation_example1`. (a) Re-derive the three-coefficient approximation @eq-approx1 and confirm the error table @tab-error1. (b) Replace the ad hoc points $x = -1\/2, 0, 1\/2$ by $N$ interior Chebyshev points (defined in @sec-chebyshev-points) and the trial function @eq-trial1 by $u_N (x) = 1 + (1 - x^2) sum_(k=0)^(N-1) a_k T_k (x)$, and assemble the collocation system for general $N$. (c) Tabulate the maximum error for $N = 2, 4, 6, 8, 10$ and confirm that the convergence is exponential, consistent with the entire, smooth character of @eq-exact1.
] <ex-meb-refine-etude1>

#exercise(title: [Collocation versus Galerkin Convergence])[
  Extend the second étude, the reaction-diffusion problem @eq-bvp2 with exact solution @eq-exact2, using the script `collocation_vs_galerkin`. (a) Implement both collocation and Galerkin for the symmetric trial space $phi_k (x) = (1 - x^2) x^(2 k)$, $k = 0, dots, N - 1$, generalising @eq-trial2. (b) For $N = 1, 2, 3, 4, 5$ compute the maximum error of each method and plot both curves against $N$ on a semilogarithmic scale, reproducing and extending @fig-collocation-vs-galerkin. (c) Quantify the roughly constant factor by which the Galerkin error lies below the collocation error, and discuss whether that gap persists as $N$ grows.
] <ex-meb-coll-gal-sweep>

#exercise(title: [The Four Weighting Strategies in Practice])[
  On the reaction-diffusion problem @eq-bvp2 with the two-function basis @eq-trial2, implement all four weighting strategies classified in @Finlayson1966: collocation, the subdomain method, least squares, and Galerkin. (a) Assemble and solve the $2 times 2$ system for each. (b) Tabulate $u(0)$ and the maximum error of each method against the exact solution @eq-exact2, extending @tab-comparison. (c) Rank the four methods by accuracy and by implementation cost, noting which require quadrature and over what integrand, and comment on whether the ranking matches the discussion in the chapter.
] <ex-meb-four-weightings>

#exercise(title: [Collocation for a Variable-Coefficient BVP])[
  Consider the boundary value problem $-(a(x) u')' = f(x)$ on $[-1, 1]$ with $u(-1) = u(1) = 0$, where $a(x) = 1 + x^2\/2$. (a) Choose $f(x)$ so that $u(x) = (1 - x^2) cos(pi x)$ is the exact solution. (b) Solve the problem using collocation with polynomial trial functions $phi_k (x) = (1 - x^2) x^(k-1)$ for $k = 1, dots, N$, using Chebyshev points as collocation nodes. (c) Tabulate the maximum error for $N = 4, 6, 8, 10, 12$ and determine whether the convergence is algebraic or exponential.
] <ex-meb-variable-coefficient>

#exercise(title: [Effect of Collocation Point Placement])[
  Solve $-u'' + u = e^x$ on $[-1, 1]$ with $u(-1) = u(1) = 0$ using collocation with trial functions $phi_k (x) = (1 - x^2) T_(k-1) (x)$. Compare the accuracy for three choices of $N = 8$ interior collocation points: (a) equispaced, (b) Chebyshev, (c) Legendre--Gauss--Lobatto. Explain the observed differences in terms of the Lebesgue constant.
] <ex-meb-point-placement>

#exercise(title: [Convergence Rate: Polynomial versus Trigonometric Bases])[
  Solve $-u'' = e^x$ on $[0, 1]$ with $u(0) = u(1) = 0$, whose exact solution is $u(x) = 1 + (e - 1) x - e^x$, using (a) polynomial trial functions $phi_k (x) = x(1 - x) x^(k-1)$ and (b) trigonometric trial functions $phi_k (x) = sin(k pi x)$, in both cases by Galerkin projection. Plot the maximum error against $N$ on a semilogarithmic scale and explain why the polynomial basis converges exponentially while the sine basis converges only algebraically. Relate the sine rate to the jump in the second derivative of the odd periodic extension of $u$ at the endpoints, where $u'' (0) = -1$ and $u'' (1) = -e$ do not vanish, and contrast it with the entire, analytic character of $u$ that the polynomial basis exploits.
] <ex-meb-poly-vs-trig>

=== Project-Style Exercises

#exercise(title: [Spectral Element Method])[
  The spectral element method#idx("spectral element method") of @Patera1984 marries the geometric flexibility of finite elements with spectral accuracy; it is the $p$-refinement strategy of @fig-refinement-strategies carried onto a mesh of a few elements, in the spirit of @Babuska1981. (a) Split $[-1, 1]$ into the two elements $[-1, 0]$ and $[0, 1]$ and represent the solution of a boundary-value problem by a moderate-degree polynomial on each, enforcing the equation by collocation inside each element and $C^0$ continuity at the shared node. (b) Verify spectral convergence as the per-element degree $p$ increases. (c) Introduce a sharp interior layer (for instance through a small diffusion coefficient) and study how clustering the element interface near the layer, a form of the $r$-refinement of @fig-refinement-strategies, improves accuracy over a single global expansion.
] <ex-meb-spectral-element>

#hint-for(<ex-meb-spectral-element>)[Assemble each element's collocation system independently, then couple the two by interface conditions: equality of the solution value and of the one-sided first derivatives at the shared node, exactly as the global trial functions of this chapter were forced to meet the boundary data.]

#exercise(title: [Radial Basis Function Collocation])[
  As the chapter notes, pseudospectral collocation is the flat-kernel limit of radial basis function#idx("radial basis function") (RBF) interpolation, and RBF methods carry spectral accuracy to scattered nodes in any dimension; the arc FD $arrow.r$ PS $arrow.r$ RBF $arrow.r$ RBF-FD is developed in @Fornberg2025. (a) Solve a one-dimensional boundary-value problem by RBF collocation, approximating $u(x) approx sum_j lambda_j thin phi(|x - x_j| \/ epsilon.alt)$ with a Gaussian or multiquadric kernel of shape parameter $epsilon.alt$, imposing the equation and the boundary conditions at the centres. (b) Study the accuracy and the matrix condition number as functions of $epsilon.alt$, exhibiting the flat-limit ($epsilon.alt arrow.r 0$) ill-conditioning and its trade-off against accuracy. (c) Verify numerically that, on a uniform grid in the flat limit, the RBF derivative approaches the pseudospectral derivative of this book.
] <ex-meb-rbf-collocation>

#hint-for(<ex-meb-rbf-collocation>)[Form the interpolation matrix $bold(Phi)_(i j) = phi(|x_i - x_j| \/ epsilon.alt)$ and a second matrix from $cal(L)$ applied to each kernel; collocating $cal(L) u = f$ at the interior centres and $u = $ data at the boundary centres gives a square system for the weights $lambda_j$. Sweep $epsilon.alt$ on a logarithmic scale to trace the accuracy-conditioning trade-off.]

#exercise(title: [Residual Minimisation by Neural Networks])[
  A physics-informed neural network#idx("physics-informed neural network") (PINN) is a Method of Weighted Residuals in which the trial function in @eq-series-expansion is replaced by a neural network and the residual @eq-residual is driven down by stochastic optimisation at sampled collocation points; operator-learning variants such as DeepONet instead learn a solution map across families of problems. The reviews @Meuris2023 and @Ngueabou2025 connect these ideas to classical spectral bases. (a) Train a small network to solve the reaction-diffusion problem @eq-bvp2, using a loss that sums the squared residual at interior points and a penalty enforcing the boundary conditions, and compare against the exact solution @eq-exact2 and the two-coefficient collocation and Galerkin solutions. (b) Replace the generic network by an architecture whose final layer is a spectral (Chebyshev or Fourier) expansion, and study whether this spectral bias accelerates training and improves accuracy. (c) Discuss, with reference to @Meuris2023, the regimes in which a learned basis can outperform the fixed polynomial bases of this chapter.
] <ex-meb-pinn>

#hint-for(<ex-meb-pinn>)[The PINN loss is a discrete weighted-residual functional: residual collocation at the interior samples plus a boundary penalty is precisely least-squares collocation with the network as trial function. A well-conditioned spectral final layer keeps a shallow network close to the classical linear solve.]

#exercise(title: [Fractional and Nonlocal Operators])[
  Spectral methods extend naturally to operators that are nonlocal in space. Replace the local second derivative in @eq-bvp2 by a fractional derivative#idx("fractional derivative") of Caputo or Riemann--Liouville type of order $alpha in (1, 2)$, following the monograph @Zayernouri2024, or by a nonlocal diffusion operator on a bounded domain in the sense of @Mustapha2025. (a) Discretise the resulting boundary-value problem by a weighted-residual method (collocation or Galerkin) in a polynomial basis, taking care with the weak singularity of the fractional kernel. (b) Validate the scheme on a manufactured solution with a known fractional derivative, such as a suitable power of $(1 - x^2)$. (c) Contrast the sparsity and structure of the resulting system matrix with the dense but small matrices of the local problems in this chapter, and discuss how non-locality reshapes the cost of the method.
] <ex-meb-fractional-nonlocal>
