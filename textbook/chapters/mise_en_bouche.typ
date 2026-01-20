// textbook/chapters/mise_en_bouche.typ
#import "../styles/template.typ": dropcap

= Mise en Bouche

#dropcap[In French cuisine, a _mise en bouche_ is a small appetizer offered by the chef to stimulate the palate before the main courses arrive. In this chapter we offer a similar intellectual appetizer: a compact, self-contained taste of spectral methods that illuminates the essential mechanics before we develop the full machinery of Fourier and Chebyshev approximation.]

Rather than jumping immediately to high-degree polynomials with $N = 100$, we perform hand calculations with just $N = 2$ or $N = 3$ unknowns. This low-dimensional setting makes every step transparent. We can verify each formula by direct computation and gain intuition that will guide us through the more sophisticated developments to come.

The techniques presented here follow the classical exposition in @Boyd2000, adapted to our pedagogical goals. The Method of Weighted Residuals provides the unifying framework that connects the collocation (pseudospectral) approach we favor in this book with the Galerkin methods that dominate finite element analysis.

== The Method of Weighted Residuals

=== Series Expansions and the Residual Function

The central idea of spectral methods is to approximate the unknown function $u(x)$ by a finite sum of basis functions:
$ u(x) approx u_N (x) = sum_(n=0)^N a_n phi_n (x), $ <eq-series-expansion>
where ${phi_n (x)}$ are known basis functions and ${a_n}$ are unknown coefficients to be determined.

When we substitute this approximation into a differential equation
$ cal(L) u = f(x), $
where $cal(L)$ is a linear differential operator, the result is generally not zero. The _residual function_ measures this discrepancy:
$ R(x; a_0, a_1, dots, a_N) = cal(L) u_N - f. $ <eq-residual>

For the exact solution, $R(x) equiv 0$. The challenge is to choose the coefficients ${a_n}$ so that the residual is as small as possible. Different spectral methods correspond to different strategies for minimizing this residual.

=== Two Minimization Strategies

The two most important strategies are:

+ *Collocation (Pseudospectral) Method*: Force the residual to be exactly zero at a set of carefully chosen points ${x_j}$, called collocation points:
  $ R(x_j; a_0, dots, a_N) = 0, quad j = 1, 2, dots, N+1. $
  This gives $N+1$ equations for $N+1$ unknowns.

+ *Galerkin Method*: Require the residual to be orthogonal to each basis function in the sense of a weighted inner product:
  $ integral_(-1)^1 R(x) phi_k (x) w(x) dif x = 0, quad k = 0, 1, dots, N, $
  where $w(x)$ is a weight function (often $w(x) = 1$ for polynomial bases).

Both methods convert the differential equation into a system of algebraic equations for the unknown coefficients. The collocation approach is simpler to implement and handles nonlinear terms easily, while the Galerkin approach often provides better global accuracy in weighted norms.

== A First Collocation Example

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

Expanding the trial function:
$ u_2 (x) = 1 + a_0 + a_1 x + a_2 x^2 - a_0 x^2 - a_1 x^3 - a_2 x^4 $
$ = (1 + a_0) + a_1 x + (a_2 - a_0) x^2 - a_1 x^3 - a_2 x^4. $

=== Computing the Residual

The residual is
$ R(x; a_0, a_1, a_2) = u_2''(x) - (4 x^2 + 2) u_2 (x). $

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

The implementation of this approximation is straightforward. In Python:

```python
def u_approx(x):
    """Evaluate the collocation approximation."""
    a0, a1, a2 = -73/118, 0, -14/59
    return 1 + (1 - x**2) * (a0 + a1*x + a2*x**2)
```

The equivalent MATLAB implementation:

```matlab
% Collocation coefficients
a0 = -73/118;  a1 = 0;  a2 = -14/59;

% Approximate solution (anonymous function)
u_approx = @(x) 1 + (1 - x.^2) .* (a0 + a1*x + a2*x.^2);
```

=== Error Analysis

The following table compares the exact and approximate solutions at several points:

#figure(
  table(
    columns: 4,
    align: (center, center, center, center),
    stroke: 0.5pt,
    inset: 6pt,
    [$x$], [$u_"exact" (x)$], [$u_"approx" (x)$], [Error],
    [$-1$], [$1.00000$], [$1.00000$], [$0.00000$],
    [$-0.5$], [$0.47237$], [$0.49153$], [$-0.01916$],
    [$0$], [$0.36788$], [$0.38136$], [$-0.01348$],
    [$0.5$], [$0.47237$], [$0.49153$], [$-0.01916$],
    [$1$], [$1.00000$], [$1.00000$], [$0.00000$],
  ),
  caption: [Comparison of exact and three-coefficient collocation approximation.],
) <tab-error1>

The maximum error is approximately $2 times 10^(-2)$, which is remarkably good for such a low-order approximation. @fig-collocation-example1 shows the solutions graphically.

#figure(
  image("../figures/ch03/python/collocation_example1.pdf", width: 95%),
  caption: [Left: exact solution $u(x) = e^(x^2 - 1)$ compared with the three-coefficient collocation approximation. The collocation points $x = -1\/2, 0, 1\/2$ are marked with squares. Right: the error $u_"exact" - u_"approx"$.],
) <fig-collocation-example1>

The code that generated this figure is available in both Python and MATLAB:
- `codes/python/ch03_mise_en_bouche/collocation_example1.py`
- `codes/matlab/ch03_mise_en_bouche/collocation_example1.m`

== Collocation versus Galerkin

To compare the two main approaches to the Method of Weighted Residuals, we consider a second example where both methods can be applied with explicit hand calculations.

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
$ u_"exact" (0) = 1/4 (1 - 1 / cosh(2)) approx 0.1834. $

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

The following Python code assembles and solves the collocation system:

```python
def solve_collocation():
    """Solve the collocation system at x = 0 and x = 0.5."""
    # Operator L[φ] = φ'' - 4φ applied to basis functions
    L_phi0 = lambda x: -2 - 4*(1 - x**2)
    L_phi1 = lambda x: (2 - 12*x**2) - 4*(x**2 - x**4)

    # Build system matrix at collocation points
    A = np.array([[L_phi0(0.0), L_phi1(0.0)],
                  [L_phi0(0.5), L_phi1(0.5)]])
    b = np.array([-1.0, -1.0])  # RHS from f = -1
    return np.linalg.solve(A, b)
```

The equivalent MATLAB implementation:

```matlab
% Operator L = d²/dx² - 4 applied to basis functions
L_phi0 = @(x) -2 - 4*(1 - x.^2);
L_phi1 = @(x) (2 - 12*x.^2) - 4*(x.^2 - x.^4);

% Build and solve collocation system
A_coll = [L_phi0(0.0), L_phi1(0.0);
          L_phi0(0.5), L_phi1(0.5)];
coeffs = A_coll \ [-1; -1];
```

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

The Galerkin method requires numerical integration to assemble the system. In Python:

```python
def solve_galerkin():
    """Solve using Galerkin: ⟨R, φₖ⟩ = 0 for k = 0, 1."""
    from scipy import integrate

    # Matrix entries: A_{ij} = ∫ L[φⱼ] φᵢ dx
    A00, _ = integrate.quad(lambda x: L_phi0(x) * phi0(x), -1, 1)
    A01, _ = integrate.quad(lambda x: L_phi1(x) * phi0(x), -1, 1)
    A10, _ = integrate.quad(lambda x: L_phi0(x) * phi1(x), -1, 1)
    A11, _ = integrate.quad(lambda x: L_phi1(x) * phi1(x), -1, 1)

    A = np.array([[A00, A01], [A10, A11]])
    b0, _ = integrate.quad(lambda x: -phi0(x), -1, 1)
    b1, _ = integrate.quad(lambda x: -phi1(x), -1, 1)
    return np.linalg.solve(A, [b0, b1])
```

The equivalent MATLAB implementation uses the built-in `integral` function:

```matlab
% Matrix entries: A_{ij} = ∫ L[φⱼ] φᵢ dx
A00 = integral(@(x) L_phi0(x) .* phi0(x), -1, 1);
A01 = integral(@(x) L_phi1(x) .* phi0(x), -1, 1);
A10 = integral(@(x) L_phi0(x) .* phi1(x), -1, 1);
A11 = integral(@(x) L_phi1(x) .* phi1(x), -1, 1);

A_gal = [A00, A01; A10, A11];

% RHS: b_i = ∫ f φᵢ dx where f = -1
b0 = integral(@(x) -phi0(x), -1, 1);
b1 = integral(@(x) -phi1(x), -1, 1);
coeffs = A_gal \ [b0; b1];
```

=== Comparison

The following table compares the two methods at the central point $x = 0$:

#figure(
  table(
    columns: 3,
    align: (left, center, center),
    stroke: 0.5pt,
    inset: 6pt,
    [*Method*], [$u(0)$], [*Absolute Error*],
    [Exact], [$0.1835$], [$0$],
    [Collocation ($N = 2$)], [$0.1829$], [$0.0006$],
    [Galerkin ($N = 2$)], [$0.1832$], [$0.0003$],
  ),
  caption: [Comparison of spectral approximations at the central maximum.],
) <tab-comparison>

For this problem, the Galerkin method is more accurate both at the central point and in a global sense. This is consistent with the error plot in @fig-collocation-vs-galerkin, which shows the Galerkin error (green) remaining closer to zero across the entire interval. The Galerkin method minimizes the error in a root-mean-square sense, which typically leads to better overall accuracy for smooth problems.

@fig-collocation-vs-galerkin shows both approximate solutions compared to the exact solution, along with the error profiles.

#figure(
  image("../figures/ch03/python/collocation_vs_galerkin.pdf", width: 95%),
  caption: [Left: exact solution compared with collocation and Galerkin approximations. The collocation points $x = 0$ and $x = 0.5$ are marked. Right: error profiles for both methods.],
) <fig-collocation-vs-galerkin>

The code that generated this figure is available in both Python and MATLAB:
- `codes/python/ch03_mise_en_bouche/collocation_vs_galerkin.py`
- `codes/matlab/ch03_mise_en_bouche/collocation_vs_galerkin.m`

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
          // Chebyshev-like polynomial (oscillates everywhere)
          #curve(
            stroke: 1.5pt + rgb("#142D6E"),
            curve.move((0cm, 0.7cm)),
            curve.line((0.4cm, 1.4cm)),
            curve.line((0.8cm, 0.5cm)),
            curve.line((1.2cm, 1.2cm)),
            curve.line((1.6cm, 0.3cm)),
            curve.line((2.0cm, 1.5cm)),
            curve.line((2.4cm, 0.4cm)),
            curve.line((2.8cm, 1.3cm)),
            curve.line((3.2cm, 0.6cm)),
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
  caption: [Schematic comparison of local (finite element) and global (spectral) basis functions. The finite element "hat function" is non-zero only over two adjacent elements, while the spectral basis function oscillates across the entire domain.],
) <fig-local-vs-global>

=== Refinement Strategies

When a numerical approximation is insufficiently accurate, there are three fundamentally different strategies to improve it, illustrated schematically in @fig-refinement-strategies:

+ *$h$-refinement*: Subdivide each element into smaller pieces, reducing the mesh spacing $h$ uniformly throughout the domain. This increases the number of elements while keeping the polynomial degree fixed.

+ *$r$-refinement* (adaptive): Redistribute the mesh points, clustering them in regions where the solution has steep gradients or other features requiring high resolution. The total number of degrees of freedom remains roughly constant.

+ *$p$-refinement*: Keep the mesh fixed while increasing $p$, the polynomial degree within each element. For a single-element domain, this is precisely what spectral methods do.

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
              curve.move((0cm, 0.3cm)),
              curve.line((1.3cm, 0.7cm)),
              curve.line((2.6cm, 0.2cm)),
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
              curve.move((0cm, 0.5cm)),
              curve.line((0.35cm, 0.8cm)),
              curve.line((0.7cm, 0.3cm)),
              curve.line((1.05cm, 0.7cm)),
              curve.line((1.4cm, 0.25cm)),
              curve.line((1.75cm, 0.75cm)),
              curve.line((2.1cm, 0.35cm)),
              curve.line((2.45cm, 0.6cm)),
              curve.line((2.6cm, 0.4cm)),
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

A natural question arises: can we combine the geometric flexibility of finite elements with the high accuracy of spectral methods? The answer is yes, through _spectral element methods_.

In spectral element methods, the domain is subdivided into elements (as in finite elements), but within each element, the polynomial degree $p$ is chosen to be moderately high, typically $p = 6$ to $8$. This hybrid approach inherits several advantages:

- The element subdivision provides geometric flexibility and matrix sparsity.
- The high polynomial degree within each element provides spectral-like accuracy.
- The theoretical framework is essentially the same as for global spectral methods.

Spectral element codes are typically written so that $p$ is a user-adjustable parameter, allowing practitioners to balance accuracy and cost for their specific application. We will not develop spectral element methods in detail in this book, but the reader should be aware that much of the theory developed for global spectral methods transfers directly to the spectral element context.

=== The Convergence of Methods at High Order

Perhaps the most profound insight from the comparison between finite element and spectral methods is this: _for sufficiently high polynomial degree, the two approaches become essentially equivalent_.

Low-order finite elements (linear, quadratic) can be derived, justified, and implemented without knowledge of Fourier or Chebyshev convergence theory. However, as the polynomial degree increases, ad hoc approaches become increasingly ill-conditioned and numerically unstable. The only practical way to implement well-behaved high-order finite elements (say, sixth order or higher) is to use the technology of spectral methods: Chebyshev or Legendre basis functions, Gaussian quadrature, and the convergence theory we will develop in subsequent chapters.

Thus, the question "Are finite elements or spectral methods better?" becomes somewhat artificial for high-order approximations. The real question is: _Does the problem at hand require high-order accuracy, or is second or fourth order sufficient?_

When the solution is smooth and high accuracy is needed, the spectral/high-order approach is clearly superior. When the solution has discontinuities, shocks, or boundary layers, or when the geometry is highly irregular, low-order methods with adaptive mesh refinement may be more practical. The wise practitioner chooses the tool appropriate to the problem.
