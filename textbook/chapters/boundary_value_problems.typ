// textbook/chapters/boundary_value_problems.typ
// Chapter 8: Boundary Value Problems
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Email: denys.dutykh@ku.ac.ae
// Homepage: https://www.denys-dutykh.com/
// Last modified: February 2026

#import "../styles/template.typ": dropcap, num, format-table, etude-conclusion, idx

// Enable equation numbering for this chapter

= Boundary Value Problems <ch-bvp>

#dropcap[With the Chebyshev differentiation matrix in hand, we are ready to tackle one of the most important classes of problems in applied mathematics: boundary value problem#idx("boundary value problem")s (BVPs). Unlike initial value problems, where we march forward in time from given initial conditions, BVPs impose constraints at multiple locations, typically at the boundaries of the domain. This spatial coupling makes BVPs inherently global, and spectral methods are ideally suited to exploit this structure.]

This chapter demonstrates how to transform differential equations into linear algebra problems. The differentiation matrix $D_N$ from @ch-chebyshev becomes the workhorse, and its square $D_N^2$ handles second-order equations. Imposing boundary conditions requires a simple "matrix surgery", removing rows and columns corresponding to boundary points. The result is a systematic approach that handles linear, variable-coefficient, and even nonlinear problems with remarkable ease.

== Second Derivatives and Matrix Squaring <sec-second-deriv>

=== The Need for $D^2$

Most BVPs in physics involve second-order derivatives: the heat equation, the wave equation, the Poisson equation#idx("Poisson equation"), and many others all feature $u_(x x)$. These classical equations are treated extensively in the foundational spectral methods texts by Gottlieb and Orszag @GottliebOrszag1977 and Canuto _et al._ @Canuto2006. To apply spectral collocation, we need a second derivative matrix.

Two approaches present themselves:

1. *Direct formulas*: Derive explicit expressions for $(D^2)_(i j)$ analogous to the first derivative formulas in @ch-chebyshev. These exist but are complex. The definitive treatment is the MATLAB Differentiation Matrix Suite of Weideman and Reddy @WeidemanReddy2000, which provides algorithmically stable recurrence relations for arbitrary-order derivatives on Chebyshev and Hermite grids.

2. *Matrix squaring*: Simply compute $D^2 = D times D$. This is $O(N^3)$ but entirely adequate for spectral $N$-values (typically $N < 200$).

We adopt the second approach for its simplicity. Trefethen @Trefethen2000 advocates this pedagogical choice, noting that optimised BLAS routines make the $O(N^3)$ matrix multiplication extremely efficient for the moderate $N$ typical of spectral methods. Baltensperger and Berrut @BaltenspergerBerrut1999 showed that for $N < 200$, the conditioning of the spectral operator ($O(N^4)$ for second derivatives) dominates the error budget long before the method of matrix construction becomes the limiting factor. The product $D_N^2$ gives us the second derivative matrix: if $bold(v)$ contains function values at the Chebyshev points, then $D_N^2 bold(v)$ approximates the second derivative values.

== Imposing Boundary Conditions <sec-boundary-conditions>

=== Dirichlet Conditions: Matrix Stripping

The most common boundary conditions specify the function values at the endpoints:
$ u(-1) = alpha, quad u(1) = beta. $
These are _Dirichlet conditions_. Boyd @Boyd2000 classifies boundary conditions into "numerical" (imposed through the discrete operator) and "behavioural" (satisfied by the basis functions automatically); matrix stripping belongs to the numerical category.

With Chebyshev points ordered as $x_0 = 1$, $x_1$, ..., $x_(N-1)$, $x_N = -1$, the boundary conditions fix $v_0 = beta$ and $v_N = alpha$. The differential equation need only be enforced at the _interior_ points $x_1, dots, x_(N-1)$.

The implementation is straightforward:
- *Remove* the first and last rows of the equation (we don't need the ODE at boundary points).
- *Remove* the first and last columns (the boundary values are known, not unknowns).

The result is an $(N-1) times (N-1)$ interior system. Denoting the entries of $D_N^2$ by $(D_N^2)_(i j)$ with $i, j = 0, 1, dots, N$, we define
$ tilde(D)^2_(i j) = (D_N^2)_(i j), quad i, j = 1, dots, N-1. $

For homogeneous conditions ($alpha = beta = 0$), the boundary terms vanish and we simply solve $tilde(D)^2 bold(u)_("int") = bold(f)_("int")$.

=== Neumann and Robin Conditions: Boundary Bordering <sec-boundary-bordering>

Matrix stripping works well for Dirichlet conditions, but _Neumann_ conditions ($u'("boundary") = g$) or _Robin_ conditions ($alpha u + beta u' = gamma$) require a different approach. Since the boundary condition involves the derivative, we cannot simply remove the boundary unknowns.

The _boundary bordering#idx("boundary bordering")_ technique keeps the full $(N + 1) times (N + 1)$ system and _replaces_ boundary rows with the appropriate conditions. This technique, codified by Boyd @Boyd2000, provides a unified framework for imposing Dirichlet, Neumann, Robin, and even nonlinear boundary conditions within the same linear algebra structure. Specifically:
- *Neumann* $u'(x_0) = g$: replace row $0$ of $D^2$ with row $0$ of $D_N$ (the first derivative matrix), and set the corresponding right-hand side entry to $g$.
- *Robin* $alpha u(x_0) + beta u'(x_0) = gamma$: replace row $0$ with $alpha bold(e)_0^top + beta D_N [0, :]$, where $bold(e)_0$ is the first unit vector.
- *Dirichlet* $u(x_0) = alpha$: replace row $0$ with $bold(e)_0^top$ and set the right-hand side to $alpha$.

This approach handles any combination of boundary conditions uniformly. We demonstrate it in @sec-mixed-bc below.

== Computational Étude 8.1: The Poisson Equation <sec-poisson-1d>

=== A Model Problem

Consider the one-dimensional Poisson equation:
$ u_(x x) = sin(pi x) + 2 cos(2 pi x), quad x in (-1, 1), quad u(plus.minus 1) = 0. $ <eq-poisson-1d>

This has the exact solution
$ u(x) = -frac(sin(pi x), pi^2) + frac(1 - cos(2 pi x), 2 pi^2). $

=== Spectral Solution

The solution procedure is direct:
1. Construct $D_N$ and compute $D_N^2$.
2. Extract the interior submatrix $tilde(D)^2$ (rows and columns $1, dots, N-1$ of $D_N^2$).
3. Evaluate the right-hand side at interior points.
4. Solve the linear system $tilde(D)^2 bold(u)_("int") = bold(f)_("int")$.
5. Embed the result in the full vector with boundary values.

@fig-poisson-1d shows the solution and demonstrates spectral convergence.

#figure(
  image("../figures/ch08/python/poisson_1d.pdf", width: 95%),
  caption: [Solution of the 1D Poisson equation @eq-poisson-1d. Left: numerical solution (circles) compared with exact solution (line) for $N = 16$. Right: exponential convergence of the maximum error, reaching machine precision by $N = 24$.],
) <fig-poisson-1d>

@tbl-poisson-convergence quantifies the convergence rate. The error decreases by roughly two orders of magnitude for each increment of two in $N$, the hallmark of spectral (exponential) convergence. The theoretical basis for this exponential convergence rate was established by Gottlieb and Orszag @GottliebOrszag1977, who proved that for analytic functions the spectral approximation error decays faster than any polynomial in $1\/N$; the monograph by Canuto _et al._ @Canuto2006 provides a comprehensive modern treatment. Machine precision is reached by $N = 24$, after which further refinement cannot improve the result due to floating-point arithmetic limitations.

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    {
      show table: format-table(auto, auto)
      table(
        columns: (auto, 1fr),
        align: (center, center),
        inset: (x: 1em, y: 0.6em),
        stroke: none,
        table.hline(stroke: 0.75pt + rgb("#142D6E")),
        table.header(
          table.cell(fill: rgb("#142D6E").lighten(85%))[*$N$*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*Max Error*],
        ),
        table.hline(stroke: 0.5pt + luma(180)),
        [4], num[5.78e-1],
        [6], num[4.44e-2],
        [8], num[3.01e-3],
        [10], num[1.33e-4],
        [12], num[7.55e-6],
        [16], num[7.84e-9],
        [20], num[3.62e-12],
        [24], num[6.07e-16],
        [32], num[3.75e-16],
        [48], num[7.91e-16],
      )
    },
  ),
  caption: [Convergence of the Chebyshev spectral method for the Poisson equation @eq-poisson-1d. The error decreases exponentially with $N$ until machine precision ($approx 2.2 times 10^(-16)$) is reached at $N = 24$.],
) <tbl-poisson-convergence>

The implementation is remarkably concise:

```python
def solve_poisson(N, f):
    """Solve u_xx = f(x) with u(±1) = 0."""
    D, x = cheb_matrix(N)
    D2 = D @ D

    # Extract interior system
    D2_int = D2[1:N, 1:N]
    f_int = f(x[1:N])

    # Solve
    u_int = np.linalg.solve(D2_int, f_int)

    # Assemble full solution
    u = np.zeros(N + 1)
    u[1:N] = u_int
    return x, u
```

```matlab
function [x, u] = solve_poisson(N, f)
% Solve u_xx = f(x) with u(±1) = 0.
    [D, x] = cheb_matrix(N);
    D2 = D * D;

    % Extract interior system
    D2_int = D2(2:N, 2:N);
    f_int = f(x(2:N));

    % Solve
    u_int = D2_int \ f_int;

    % Assemble full solution
    u = zeros(N+1, 1);
    u(2:N) = u_int;
end
```

The Julia implementation:

```julia
function solve_poisson(N, f)
    # Solve u_xx = f(x) with u(±1) = 0.
    D, x = cheb_matrix(N)
    D2 = D * D

    # Extract interior system
    D2_int = D2[2:N, 2:N]
    f_int = f.(x[2:N])

    # Solve and assemble
    u_int = D2_int \ f_int
    u = zeros(N + 1)
    u[2:N] = u_int
    return x, u
end
```

#etude-conclusion[
  The convergence data reveal the defining advantage of spectral methods: *exponential* convergence. The error drops by roughly two orders of magnitude for each increment of two in $N$, and machine precision is attained by $N = 24$ --- 23 interior unknowns. A second-order finite-difference scheme would need thousands of grid points to match this, and could never reach $10^(-16)$ because of its algebraic $cal(O)(h^2)$ rate. Equally striking is the simplicity: the entire solution procedure reduces to *form $D_N^2$, strip the boundary rows and columns, solve the linear system*. This "matrix surgery" converts a continuous BVP into a small dense linear-algebra problem that takes microseconds. The dense coupling --- every collocation point influences every other --- is precisely what encodes the global polynomial approximation that makes spectral convergence possible.
]

The code generating @fig-poisson-1d is available in:
- `codes/python/ch08/bvp_linear.py`
- `codes/matlab/ch08/bvp_linear.m`
- `codes/julia/ch08/bvp_linear.jl`

== Computational Étude 8.2: Variable Coefficient Problems <sec-variable-coeff>

=== The Airy-Type Equation

Variable coefficients pose no additional difficulty for spectral methods. Consider:
$ u_(x x) - (1 + x^2) u = 1, quad x in (-1, 1), quad u(plus.minus 1) = 0. $ <eq-variable-coeff>

The variable coefficient#idx("variable coefficient") $(1 + x^2)$ becomes a diagonal matrix. The discretized operator is:
$ L = D_N^2 - "diag"(1 + x^2). $

After extracting the interior system, we solve $tilde(L) bold(u)_("int") = bold(1)_("int")$.

@fig-variable-coeff compares this solution with the constant-coefficient case.

#figure(
  image("../figures/ch08/python/variable_coeff.pdf", width: 95%),
  caption: [Variable coefficient BVP @eq-variable-coeff. Left: comparison of solutions for variable coefficient $(1 + x^2)$ and constant coefficient $(1)$. The variable coefficient reduces the solution amplitude, especially near the boundaries where $1 + x^2$ is largest. Right: verification of the constant-coefficient case against its exact solution.],
) <fig-variable-coeff>

The implementation requires only a minor modification to the Poisson solver:

```python
def solve_variable_coeff(N, coeff_func):
    """Solve u_xx - c(x)*u = 1 with u(±1) = 0."""
    D, x = cheb_matrix(N)
    D2 = D @ D

    # Build operator L = D² - diag(c(x))
    c = coeff_func(x)
    L = D2 - np.diag(c)

    # Extract interior system
    L_int = L[1:N, 1:N]
    rhs = np.ones(N - 1)

    # Solve and assemble
    u_int = np.linalg.solve(L_int, rhs)
    u = np.zeros(N + 1)
    u[1:N] = u_int
    return x, u
```

```matlab
function [x, u] = solve_variable_coeff(N, coeff_func)
% Solve u_xx - c(x)*u = 1 with u(±1) = 0.
    [D, x] = cheb_matrix(N);
    D2 = D * D;

    % Build operator L = D² - diag(c(x))
    c = coeff_func(x);
    L = D2 - diag(c);

    % Extract interior system
    L_int = L(2:N, 2:N);
    rhs = ones(N-1, 1);

    % Solve and assemble
    u_int = L_int \ rhs;
    u = zeros(N+1, 1);
    u(2:N) = u_int;
end
```

The Julia implementation:

```julia
function solve_variable_coeff(N, coeff_func)
    # Solve u_xx - c(x)*u = 1 with u(±1) = 0.
    D, x = cheb_matrix(N)
    D2 = D * D

    # Build operator L = D² - diag(c(x))
    c = coeff_func.(x)
    L = D2 - Diagonal(c)

    # Extract interior system
    L_int = L[2:N, 2:N]
    rhs = ones(N - 1)

    # Solve and assemble
    u_int = L_int \ rhs
    u = zeros(N + 1)
    u[2:N] = u_int
    return x, u
end
```

#etude-conclusion[
  The variable coefficient $(1 + x^2)$ is absorbed by a *single diagonal multiplication*: $L = D_N^2 - "diag"(1 + x^2)$ --- a direct consequence of collocation. There is no quadrature, no basis-function integral, no algorithmic complexity beyond the constant-coefficient case. The spectral rate is fully preserved because $(1 + x^2)$ is itself a low-degree polynomial that the Chebyshev grid represents exactly. Beyond the algorithmic point, the figure also yields a *physical* insight that pure theory obscures: the spatially varying stiffness suppresses the solution amplitude near the boundaries where the coefficient is largest --- spatially inhomogeneous damping that would be hard to predict from the equation alone, and a small demonstration of computation as a tool for building physical intuition.
]

The code generating @fig-variable-coeff is available in:
- `codes/python/ch08/bvp_variable_coeff.py`
- `codes/matlab/ch08/bvp_variable_coeff.m`
- `codes/julia/ch08/bvp_variable_coeff.jl`

== Computational Étude 8.3: Mixed Boundary Conditions <sec-mixed-bc>

=== Boundary Bordering in Practice

Not all problems come with Dirichlet conditions at both endpoints. Many physical situations involve _Neumann_ conditions (prescribed flux or insulation) or a mix of condition types. Consider the steady-state heat equation with a fixed temperature at one end and an insulated boundary at the other:
$ u_(x x) = -cos(pi x \/ 2), quad x in (-1, 1), quad u(-1) = 0, quad u'(1) = 0. $ <eq-mixed-bc>
The Dirichlet condition $u(-1) = 0$ fixes the temperature at the left boundary, while the Neumann condition $u'(1) = 0$ models perfect insulation at the right boundary (zero heat flux).

The exact solution is
$ u(x) = frac(4, pi^2) cos(pi x \/ 2) + frac(2, pi) (x + 1), $
which one can verify satisfies the ODE and both boundary conditions.

Following the boundary bordering approach from @sec-boundary-bordering, we keep the full $(N + 1) times (N + 1)$ system and replace boundary rows:
+ Start with $L = D_N^2$ and $bold(f) = f(bold(x))$.
+ *Row $0$* (Neumann at $x_0 = 1$): replace $L[0, :] = D_N [0, :]$ and set $f_0 = 0$.
+ *Row $N$* (Dirichlet at $x_N = -1$): replace $L[N, :] = bold(e)_N^top$ and set $f_N = 0$.
+ Solve $L bold(u) = bold(f)$.

In Python:

```python
def solve_mixed_bc(N):
    """Solve u_xx = f(x) with u(-1) = 0, u'(1) = 0."""
    D2, D, x = cheb_second_derivative_matrix(N)

    L = D2.copy()
    rhs = rhs_function(x)

    # Neumann at x[0] = 1: replace with derivative row
    L[0, :] = D[0, :]
    rhs[0] = 0.0

    # Dirichlet at x[N] = -1: replace with identity row
    L[N, :] = 0.0
    L[N, N] = 1.0
    rhs[N] = 0.0

    return x, np.linalg.solve(L, rhs)
```

The equivalent MATLAB implementation:

```matlab
function [x, u] = solve_mixed_bc(N, rhs_func)
    [D, x] = cheb_matrix(N);
    D2 = D * D;

    L = D2;
    rhs = rhs_func(x);

    % Neumann at x(1) = 1: replace with derivative row
    L(1, :) = D(1, :);
    rhs(1) = 0;

    % Dirichlet at x(N+1) = -1: replace with identity row
    L(N+1, :) = 0;
    L(N+1, N+1) = 1;
    rhs(N+1) = 0;

    u = L \ rhs;
end
```

And in Julia:

```julia
function solve_mixed_bc(N)
    D2, D, x = cheb_second_derivative_matrix(N)

    L = copy(D2)
    rhs = rhs_function.(x)

    # Neumann at x[1] = 1: replace with derivative row
    L[1, :] = D[1, :]
    rhs[1] = 0.0

    # Dirichlet at x[N+1] = -1: replace with identity row
    L[N+1, :] .= 0.0
    L[N+1, N+1] = 1.0
    rhs[N+1] = 0.0

    return x, L \ rhs
end
```

@fig-mixed-bc shows the solution and convergence behavior.

#figure(
  image("../figures/ch08/python/mixed_bc.pdf", width: 95%),
  caption: [Mixed boundary conditions @eq-mixed-bc. Left: the solution rises from $u(-1) = 0$ (Dirichlet) to a maximum at $x = 1$ where the slope is zero (Neumann, insulated). Right: spectral convergence reaches approximately $10^(-14)$ by $N = 16$, after which round-off effects from the large entries in $D_N$ near the boundary cause a slight accuracy plateau.],
) <fig-mixed-bc>

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    {
      show table: format-table(auto, auto)
      table(
        columns: (auto, 1fr),
        align: (center, center),
        inset: (x: 1em, y: 0.6em),
        stroke: none,
        table.hline(stroke: 0.75pt + rgb("#142D6E")),
        table.header(
          table.cell(fill: rgb("#142D6E").lighten(85%))[*$N$*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*Max Error*],
        ),
        table.hline(stroke: 0.5pt + luma(180)),
        [4],  num[1.46e-2],
        [6],  num[2.04e-4],
        [8],  num[1.70e-6],
        [10], num[9.40e-9],
        [12], num[3.69e-11],
        [16], num[8.06e-14],
        [20], num[1.25e-13],
        [24], num[3.11e-14],
        [32], num[3.79e-13],
        [48], num[1.80e-12],
      )
    },
  ),
  caption: [Convergence for the mixed BC problem @eq-mixed-bc. The error decreases spectrally until $N approx 16$, reaching $approx 10^(-14)$. Unlike the pure Dirichlet case (@tbl-poisson-convergence), which reaches machine precision at $10^(-16)$, the Neumann condition enforced through the first derivative matrix $D_N$ introduces a slightly higher noise floor. This is because the entries of $D_N$ at the boundary grow as $cal(O)(N^2)$, amplifying round-off errors.],
) <tbl-mixed-bc-convergence>

#etude-conclusion[
  The *boundary bordering* technique provides a unified framework for Dirichlet, Neumann, and Robin conditions within the same linear-algebra structure. Rather than redesigning the discretisation for each boundary type, one simply replaces the appropriate rows of the full $(N + 1) times (N + 1)$ system: an identity row for Dirichlet, a first-derivative row for Neumann, a linear combination for Robin. This flexibility is a hallmark of collocation; Galerkin methods often need different trial-function spaces for different boundary types. The convergence data also expose a subtle limitation: the Neumann case has a higher noise floor ($approx 10^(-14)$) than the pure Dirichlet case ($approx 10^(-16)$), because $D_N$ entries at the boundary grow as $cal(O)(N^2)$ and amplify round-off in the boundary row. The convergence *rate* is exponential for every BC type; the achievable *floor* depends on the type.
]

The code generating @fig-mixed-bc is available in:
- `codes/python/ch08/bvp_mixed_bc.py`
- `codes/matlab/ch08/bvp_mixed_bc.m`
- `codes/julia/ch08/bvp_mixed_bc.jl`

== Computational Étude 8.4: The Bratu Equation <sec-bratu>

=== A Classic Nonlinear Problem

The Bratu equation#idx("Bratu equation") models combustion and thermal explosion. The equation was originally derived by Frank-Kamenetskii @FrankKamenetskii1955 in the context of thermal ignition theory, where the exponential nonlinearity models the Arrhenius temperature dependence of reaction rates:
$ u_(x x) + lambda e^u = 0, quad x in (-1, 1), quad u(plus.minus 1) = 0. $ <eq-bratu>

This equation exhibits a _turning point_ phenomenon: solutions exist only for $lambda lt.eq.slant lambda_c$, where $lambda_c approx 0.878$ for the domain $[-1, 1]$. Boyd @Boyd1986Bratu produced the definitive spectral analysis of this bifurcation, using arc-length continuation to resolve the turning point singularity and compute $lambda_c$ to high precision. Above this critical value, no solution exists. For $lambda = 0.5$ (well below the critical value), a unique solution exists.

The nonlinearity $e^u$ prevents direct linear algebra. Instead, we use _Newton iteration#idx("Newton iteration")_: linearize, solve, update, repeat.

=== Newton Iteration

To solve a nonlinear equation, we cannot simply invert a matrix. Instead, we employ an _iterative_ strategy: start from an initial guess, and progressively refine it until the solution is accurate enough.

Newton's method is the most widely used approach for nonlinear equations --- a comprehensive treatment, including convergence theory and affine invariance, can be found in Deuflhard @Deuflhard2011. The key idea is _linearization_: at each iteration, we approximate the nonlinear problem by a linear one, solve that linear problem exactly, and use the result to improve our current approximation. Concretely, suppose we have a current approximation $bold(u)^((k))$. We seek a correction $delta bold(u)$ such that $bold(u)^((k)) + delta bold(u)$ satisfies the equation more closely. Expanding the residual to first order:
$ bold(F)(bold(u)^((k)) + delta bold(u)) approx bold(F)(bold(u)^((k))) + J(bold(u)^((k))) delta bold(u) = bold(0), $
where $J = partial bold(F) \/ partial bold(u)$ is the _Jacobian matrix_. Setting this to zero gives the Newton step:
$ J(bold(u)^((k))) delta bold(u) = -bold(F)(bold(u)^((k))). $ <eq-newton-step>

For our Bratu problem, the discrete residual is $bold(F)(bold(u)) = tilde(D)^2 bold(u) + lambda e^(bold(u))$, and the Jacobian is:
$ J(bold(u)) = tilde(D)^2 + lambda "diag"(e^(bold(u))). $
Note that the Jacobian changes at every iteration because $"diag"(e^(bold(u)))$ depends on the current solution.

The complete algorithm is:

+ Choose an initial guess $bold(u)^((0))$ (we use $bold(u)^((0)) = bold(0)$).
+ For $k = 0, 1, 2, dots$:
  + Compute the residual $bold(F)^((k)) = tilde(D)^2 bold(u)^((k)) + lambda e^(bold(u)^((k)))$.
  + Compute the Jacobian $J^((k)) = tilde(D)^2 + lambda "diag"(e^(bold(u)^((k))))$.
  + Solve the linear system $J^((k)) delta bold(u) = -bold(F)^((k))$ for the Newton correction $delta bold(u)$.
  + Update: $bold(u)^((k+1)) = bold(u)^((k)) + delta bold(u)$.
  + If the stopping criterion is satisfied, terminate.

*Stopping criterion.* A critical practical detail is deciding when to stop iterating. Common choices include:
- _Correction norm_: $norm(delta bold(u))_infinity < "tol"$, i.e., the Newton step becomes negligibly small.
- _Residual norm_: $norm(bold(F)^((k)))_infinity < "tol"$, i.e., the equation is satisfied to the desired accuracy.
- _Relative correction_: $norm(delta bold(u))_infinity \/ norm(bold(u)^((k)))_infinity < "tol"$, which accounts for the magnitude of the solution.

In our implementation, we use the correction norm $norm(delta bold(u))_infinity < 10^(-10)$ as the stopping criterion, combined with a maximum iteration count of 50 as a safety guard against non-convergence. This is a natural choice: when the correction becomes smaller than $10^(-10)$, the solution has stabilized to at least 10 significant digits.

*Convergence rate.* Newton's method enjoys _quadratic convergence_ when sufficiently close to a solution: the error at each step is roughly the square of the error at the previous step. This means that once the iteration begins to converge, the number of correct digits approximately doubles with each step. In practice, for $lambda = 0.5$, Newton's method converges in about 5 to 8 iterations starting from $bold(u)^((0)) = bold(0)$, as shown in @fig-bratu.

In the contemporary era, interest has shifted to fractional-order generalisations of the Bratu equation. Salama _et al._ @Salama2025 have extended spectral collocation using shifted Lucas polynomials to solve fractional Bratu equations, modelling non-local memory effects in combustion. Gheorghiu @Gheorghiu2020 revisited the classical problem to benchmark the accuracy of modern Newton--Kantorovich solvers, achieving spectral accuracy with as few as 8--16 grid points.

@fig-bratu shows the solution and convergence behavior.

#figure(
  image("../figures/ch08/python/bratu.pdf", width: 95%),
  caption: [Nonlinear BVP: the Bratu equation @eq-bratu with $lambda = 0.5$. Left: solutions for different values of $N$. Right: convergence history comparing Newton iteration (quadratic convergence) with fixed-point iteration (linear convergence). Newton typically converges in $5$--$8$ iterations.],
) <fig-bratu>

The Newton iteration is straightforward to implement:

```python
def solve_bratu_newton(N, lam=0.5, tol=1e-10, max_iter=50):
    """Solve u_xx + λ*exp(u) = 0 with u(±1) = 0 using Newton iteration."""
    D, x = cheb_matrix(N)
    D2 = D @ D
    D2_int = D2[1:N, 1:N]

    u = np.zeros(N - 1)  # Initial guess

    for k in range(max_iter):
        exp_u = np.exp(u)
        # Residual: F = D²u + λ*exp(u)
        F = D2_int @ u + lam * exp_u
        # Jacobian: J = D² + λ*diag(exp(u))
        J = D2_int + lam * np.diag(exp_u)
        # Newton step
        delta_u = np.linalg.solve(J, -F)
        u = u + delta_u

        if np.max(np.abs(delta_u)) < tol:
            break

    # Assemble full solution
    u_full = np.zeros(N + 1)
    u_full[1:N] = u
    return x, u_full, k + 1
```

```matlab
function [x, u_full, iterations] = solve_bratu_newton(N, lam, tol, max_iter)
% Solve u_xx + λ*exp(u) = 0 with u(±1) = 0 using Newton iteration.
    if nargin < 4, max_iter = 50; end
    if nargin < 3, tol = 1e-10; end
    if nargin < 2, lam = 0.5; end

    [D, x] = cheb_matrix(N);
    D2 = D * D;
    D2_int = D2(2:N, 2:N);

    u = zeros(N-1, 1);  % Initial guess

    for k = 1:max_iter
        exp_u = exp(u);
        % Residual: F = D²u + λ*exp(u)
        F = D2_int * u + lam * exp_u;
        % Jacobian: J = D² + λ*diag(exp(u))
        J = D2_int + lam * diag(exp_u);
        % Newton step
        delta_u = J \ (-F);
        u = u + delta_u;

        if max(abs(delta_u)) < tol
            break
        end
    end

    % Assemble full solution
    u_full = zeros(N+1, 1);
    u_full(2:N) = u;
    iterations = k;
end
```

The Julia implementation:

```julia
function solve_bratu_newton(N; lam=0.5, tol=1e-10, max_iter=50)
    # Solve u_xx + λ*exp(u) = 0 with u(±1) = 0 using Newton iteration.
    D, x = cheb_matrix(N)
    D2 = D * D
    D2_int = D2[2:N, 2:N]

    u = zeros(N - 1)  # Initial guess

    for k in 1:max_iter
        exp_u = exp.(u)
        # Residual and Jacobian
        F = D2_int * u + lam * exp_u
        J = D2_int + lam * Diagonal(exp_u)
        # Newton step
        delta_u = J \ (-F)
        u .+= delta_u

        maximum(abs.(delta_u)) < tol && break
    end

    # Assemble full solution
    u_full = zeros(N + 1)
    u_full[2:N] = u
    return x, u_full
end
```

#etude-conclusion[
  The Bratu equation shows how spectral collocation converts a *nonlinear* PDE into a nonlinear algebraic system directly amenable to Newton iteration. The Jacobian $J = tilde(D)^2 + lambda "diag"(e^(bold(u)))$ is simply the linear spectral operator augmented by a diagonal matrix that updates at each step --- the same structural simplicity as the variable-coefficient case extended into the nonlinear regime. Newton's *quadratic* convergence is on full display: once the iterates enter the basin of attraction, the number of correct digits roughly doubles at each step, with $5$--$8$ iterations sufficing from a zero initial guess. The right panel contrasts this with a Picard fixed-point iteration whose linear rate demands many more iterations --- the cost of assembling and factoring the Jacobian is more than compensated by the rapid residual reduction. Note: for $lambda < lambda_c$ Bratu has *two* branches; the zero initial guess naturally selects the lower one. Reaching the upper branch requires continuation or a carefully chosen warm start.
]

The code generating @fig-bratu is available in:
- `codes/python/ch08/bvp_nonlinear.py`
- `codes/matlab/ch08/bvp_nonlinear.m`
- `codes/julia/ch08/bvp_nonlinear.jl`

== Computational Étude 8.5: Eigenvalue Problems <sec-eigenvalue>

=== Resolution Limits

The eigenvalue problem for the Laplacian,
$ u_(x x) = lambda u, quad x in (-1, 1), quad u(plus.minus 1) = 0, $ <eq-eigenvalue>
has exact eigenvalues $lambda_n = -(n pi \/ 2)^2$ and eigenfunctions $u_n (x) = sin(n pi (x + 1) \/ 2)$.

Spectral methods compute these eigenvalues with spectral accuracy, _but only for the low modes_. High-frequency modes require many points per wavelength (ppw) for accurate representation. The rule of thumb is:
$ "Need" gt.eq.slant pi " points per wavelength for spectral accuracy." $
This resolution criterion was established in the foundational work of Kreiss and Oliger @KreissOliger1972, who demonstrated the superior efficiency of spectral methods over finite differences for wave propagation. Recent dispersion analysis by Gkanos _et al._ @Gkanos2025 confirms that while spectral methods minimise phase lag significantly better than finite differences, the "$pi$ points" rule is an asymptotic stability limit rather than an accuracy guarantee for large domains.

@fig-eigenvalue demonstrates this resolution limit.

#figure(
  image("../figures/ch08/python/eigenvalue_problem.pdf", width: 95%),
  caption: [Eigenvalue problem @eq-eigenvalue for $N = 36$. Six eigenmodes are shown with their points per wavelength (ppw) and eigenvalue errors. Low modes (high ppw) are resolved to near machine precision. As ppw drops below $pi approx 3.14$, accuracy degrades rapidly. This is not a bug but a fundamental resolution limit.],
) <fig-eigenvalue>

@tbl-eigenvalue-accuracy quantifies this resolution effect. For the first few modes, the points per wavelength is large and the spectral method delivers near machine precision. As the mode number increases and the ppw drops toward $pi approx 3.14$, accuracy degrades sharply. Beyond this threshold, the computed eigenvalues bear little resemblance to the exact values.

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    {
      show table: format-table(auto, auto, auto, auto, auto, auto)
      table(
        columns: (auto, auto, auto, auto, auto, auto),
        align: (center, center, center, center, center, center),
        inset: (x: 0.7em, y: 0.6em),
        stroke: none,
        table.hline(stroke: 0.75pt + rgb("#142D6E")),
        table.header(
          table.cell(fill: rgb("#142D6E").lighten(85%))[*$n$*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*ppw*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*$lambda_n$ (exact)*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*$lambda_n$ (numerical)*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*Absolute Error*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*Relative Error*],
        ),
        table.hline(stroke: 0.5pt + luma(180)),
        [1],  num[72.00], num[-2.4674],    num[-2.4674],    num[5.68e-14], num[2.30e-14],
        [5],  num[14.40], num[-61.6850],   num[-61.6850],   num[2.77e-13], num[4.49e-15],
        [10], num[7.20],  num[-246.7401],  num[-246.7401],  num[5.56e-11], num[2.25e-13],
        [15], num[4.80],  num[-555.1652],  num[-555.1653],  num[1.97e-5],  num[3.56e-8],
        [20], num[3.60],  num[-986.9604],  num[-988.0301],  num[1.07e0],   num[1.08e-3],
        [23], num[3.13],  num[-1305.2552], num[-1270.1437], num[3.51e1],   num[2.69e-2],
        [25], num[2.88],  num[-1542.1257], num[-1567.3682], num[2.52e1],   num[1.64e-2],
        [30], num[2.40],  num[-2220.6610], num[-5860.9103], num[3.64e3],   num[1.64e0],
        [35], num[2.06],  num[-3022.5663], num[-79900.13],  num[7.69e4],   num[2.54e1],
      )
    },
  ),
  caption: [Eigenvalue accuracy for $N = 36$. The points per wavelength (ppw) equals $2 N \/ n$. For modes with ppw $gt.eq.slant pi$, both absolute and relative errors are near machine precision. The transition at ppw $approx pi$ is remarkably sharp: mode 23 (ppw $= 3.13 approx pi$) already shows $approx 3%$ relative error, and beyond this threshold accuracy collapses completely.],
) <tbl-eigenvalue-accuracy>

The eigenvalue computation uses standard linear algebra:

```python
def compute_laplacian_eigenvalues(N):
    """Compute eigenvalues of u_xx = λu with u(±1) = 0."""
    D, x = cheb_matrix(N)
    D2 = D @ D
    D2_int = D2[1:N, 1:N]

    # Compute eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eig(D2_int)

    # Sort by magnitude (most negative first)
    idx = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Exact eigenvalues: λ_n = -(nπ/2)²
    n = np.arange(1, N)
    exact_eigenvalues = -(n * np.pi / 2)**2

    return eigenvalues, eigenvectors, exact_eigenvalues
```

```matlab
function [eigenvalues, eigenvectors, exact_eigenvalues] = compute_laplacian_eigenvalues(N)
% Compute eigenvalues of u_xx = λu with u(±1) = 0.
    [D, x] = cheb_matrix(N);
    D2 = D * D;
    D2_int = D2(2:N, 2:N);

    % Compute eigenvalues and eigenvectors
    [eigenvectors, Lambda] = eig(D2_int);
    eigenvalues = diag(Lambda);

    % Sort by magnitude (most negative first)
    [eigenvalues, idx] = sort(eigenvalues);
    eigenvectors = eigenvectors(:, idx);

    % Exact eigenvalues: λ_n = -(nπ/2)²
    n = (1:N-1)';
    exact_eigenvalues = -(n * pi / 2).^2;
end
```

The Julia implementation:

```julia
function compute_laplacian_eigenvalues(N)
    # Compute eigenvalues of u_xx = λu with u(±1) = 0.
    D, x = cheb_matrix(N)
    D2 = D * D
    D2_int = D2[2:N, 2:N]

    # Compute eigenvalues and eigenvectors
    eig = eigen(D2_int)
    idx = sortperm(real.(eig.values))
    eigenvalues = real.(eig.values[idx])
    eigenvectors = real.(eig.vectors[:, idx])

    # Exact eigenvalues: λ_n = -(nπ/2)²
    exact_eigenvalues = [-(n * π / 2)^2 for n in 1:N-1]
    return eigenvalues, eigenvectors, exact_eigenvalues
end
```

#etude-conclusion[
  The eigenvalue computation exposes a *fundamental resolution limit* that is intrinsic to all discrete methods but especially transparent in the spectral setting. Low modes with many points-per-wavelength (ppw) reach machine precision: mode $n = 1$ with ppw $= 72$ achieves $2.3 times 10^(-14)$ relative error, and mode $n = 10$ with ppw $= 7.2$ retains $10^(-13)$. The transition at ppw $approx pi$ is remarkably abrupt: mode $n = 23$ (ppw $= 3.13$) is already $approx 3 %$ off, and beyond this threshold the computed eigenvalues become entirely spurious. The $(N - 1) times (N - 1)$ matrix necessarily produces $N - 1$ eigenvalues but only the first $approx 2 N \/ pi$ are physically meaningful; the rest are *polluted* artefacts. The practical rule of thumb: *discard any mode with fewer than $pi$ points per wavelength*.
]

The code generating @fig-eigenvalue is available in:
- `codes/python/ch08/bvp_eigenvalue.py`
- `codes/matlab/ch08/bvp_eigenvalue.m`
- `codes/julia/ch08/bvp_eigenvalue.jl`

== Computational Étude 8.6: Two-Dimensional Poisson Problem <sec-2d>

=== Tensor Products

For problems on a rectangle $[-1, 1]^2$, we use _tensor product_ grids: Chebyshev points in both $x$ and $y$ directions. The total number of grid points is $(N + 1)^2$.

#block(
  width: 100%,
  inset: (x: 1.2em, y: 0.9em),
  stroke: (left: 2pt + rgb("#142D6E").lighten(50%)),
  fill: rgb("#142D6E").lighten(95%),
)[
  *Remark (Padua points).* Tensor product grids are not the only option for bivariate polynomial interpolation on the square. The _Padua points_, discovered by De Marchi, Caliari, and Vianello @DeMarchi2005, are the first known example (and to date the only one) of a unisolvent point set on $[-1, 1]^2$ with _minimal growth_ of the Lebesgue constant, proven to be $cal(O)(log^2 n)$ by Bos, Caliari, De Marchi, Vianello, and Xu @Bos2006. For total polynomial degree $n$, Padua points require only $(n + 1)(n + 2) \/ 2$ points, roughly half the $(n + 1)^2$ needed by the tensor product grid#idx("tensor product grid"). The points admit an elegant geometric construction: they lie exactly on the self-intersections and boundary contacts of a generating Lissajous curve on $[-1, 1]^2$, which gives rise to four distinct families obtained by successive 90-degree rotations @Bos2007. An efficient implementation is available as Algorithm 886 @Caliari2008. While we use tensor product grids throughout this chapter for their simplicity and natural compatibility with the Kronecker product#idx("Kronecker product") structure of differential operators, Padua points offer a more economical alternative for pure interpolation and approximation problems.
]

The 2D Laplacian operator is built using _Kronecker products_:
$ L = I times.o D^2 + D^2 times.o I, $ <eq-laplacian-2d>
where $I$ is the identity matrix and $times.o$ denotes the Kronecker product.

@fig-tensor-grid shows the tensor product grid structure.

#figure(
  image("../figures/ch08/python/tensor_grid.pdf", width: 55%),
  caption: [Chebyshev tensor product grid for $N = 16$. Points cluster near all four boundaries, providing the boundary resolution needed for spectral accuracy in two dimensions.],
) <fig-tensor-grid>

=== 2D Poisson Problem

Consider the 2D Poisson equation:
$ u_(x x) + u_(y y) = -2 pi^2 sin(pi x) sin(pi y), quad u = 0 "on boundary." $ <eq-poisson-2d>

The exact solution is $u(x, y) = sin(pi x) sin(pi y)$.

@fig-poisson-2d shows the solution.

#figure(
  image("../figures/ch08/python/poisson_2d.pdf", width: 95%),
  caption: [Solution of the 2D Poisson equation @eq-poisson-2d. Left: 3D surface plot of the solution. Right: contour plot with Chebyshev grid points overlaid. For $N = 16$, the maximum error is approximately $10^(-12)$.],
) <fig-poisson-2d>

The Kronecker product formulation leads to concise code:

```python
def solve_poisson_2d(N, f):
    """Solve ∇²u = f on [-1,1]² with u = 0 on boundary."""
    D, x = cheb_matrix(N)
    D2 = D @ D
    D2_int = D2[1:N, 1:N]
    x_int = x[1:N]

    # Build 2D Laplacian: L = I ⊗ D² + D² ⊗ I
    I = np.eye(N - 1)
    L = np.kron(I, D2_int) + np.kron(D2_int, I)

    # Right-hand side on interior grid
    X, Y = np.meshgrid(x_int, x_int)
    F = f(X, Y).flatten(order='F')

    # Solve and reshape
    u_vec = np.linalg.solve(L, F)
    U_int = u_vec.reshape((N-1, N-1), order='F')

    # Embed in full grid with zero boundary
    U = np.zeros((N+1, N+1))
    U[1:N, 1:N] = U_int
    return np.meshgrid(x, x), U
```

```matlab
function [grids, U] = solve_poisson_2d(N, f)
% Solve ∇²u = f on [-1,1]² with u = 0 on boundary.
    [D, x] = cheb_matrix(N);
    D2 = D * D;
    D2_int = D2(2:N, 2:N);
    x_int = x(2:N);

    % Build 2D Laplacian: L = I ⊗ D² + D² ⊗ I
    I = eye(N - 1);
    L = kron(I, D2_int) + kron(D2_int, I);

    % Right-hand side on interior grid
    [X, Y] = meshgrid(x_int, x_int);
    F = f(X, Y);

    % Solve and reshape
    u_vec = L \ F(:);
    U_int = reshape(u_vec, N-1, N-1);

    % Embed in full grid with zero boundary
    U = zeros(N+1, N+1);
    U(2:N, 2:N) = U_int;
    [grids{1}, grids{2}] = meshgrid(x, x);
end
```

The Julia implementation:

```julia
function solve_poisson_2d(N, f)
    # Solve ∇²u = f on [-1,1]² with u = 0 on boundary.
    D, x = cheb_matrix(N)
    D2 = D * D
    D2_int = D2[2:N, 2:N]
    x_int = x[2:N]
    n_int = N - 1

    # Build 2D Laplacian: L = I ⊗ D² + D² ⊗ I
    Imat = I(n_int)
    L = kron(Imat, D2_int) + kron(D2_int, Imat)

    # Right-hand side on interior grid
    F = [f(x_int[j], x_int[i]) for i in 1:n_int, j in 1:n_int]

    # Solve and reshape
    u_vec = L \ vec(F)
    U_int = reshape(u_vec, n_int, n_int)

    # Embed in full grid with zero boundary
    U = zeros(N + 1, N + 1)
    U[2:N, 2:N] = U_int
    return x, U
end
```

The sparsity pattern of the 2D Laplacian operator, shown in @fig-laplacian-sparsity, reveals the Kronecker product structure.

#figure(
  image("../figures/ch08/python/laplacian_sparsity.pdf", width: 55%),
  caption: [Sparsity pattern of the 2D Laplacian matrix for $N = 16$. The $(N-1)^2 times (N-1)^2$ matrix shows the characteristic block structure arising from the Kronecker product formulation @eq-laplacian-2d.],
) <fig-laplacian-sparsity>

Despite the visual impression of sparsity in @fig-laplacian-sparsity, a closer look reveals an important structural property: the Chebyshev second derivative matrix $D^2$ is _fully dense_ (every collocation point couples to every other). The apparent sparsity arises entirely from the Kronecker product structure. The term $I times.o D^2$ contributes $(N - 1)$ dense diagonal blocks of size $(N - 1) times (N - 1)$, while $D^2 times.o I$ scatters $(N - 1) times (N - 1)$ identity-sized blocks across the off-diagonal positions. By inclusion-exclusion (subtracting the $(N - 1)^2$ diagonal entries counted twice), the exact number of nonzero entries is
$ "nnz"(L) = (N - 1)^2 (2(N - 1) - 1) = (N - 1)^2 (2 N - 3), $
so that each row contains exactly $2(N - 1) - 1$ nonzeros.

@tbl-laplacian-sparsity quantifies how the sparsity evolves with $N$.

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    {
      show table: format-table(auto, none, auto, auto, auto)
      table(
        columns: (auto, auto, auto, auto, auto),
        align: (center, center, center, center, center),
        inset: (x: 1em, y: 0.6em),
        stroke: none,
        table.hline(stroke: 0.75pt + rgb("#142D6E")),
        table.header(
          table.cell(fill: rgb("#142D6E").lighten(85%))[*$N$*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*Matrix Size*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*Nonzeros*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*Density (%)*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*nnz/row*],
        ),
        table.hline(stroke: 0.5pt + luma(180)),
        [8],  [$49 times 49$],     num[637],    num[26.5],  [13],
        [12], [$121 times 121$],    num[2541],   num[17.4],  [21],
        [16], [$225 times 225$],    num[6525],   num[12.9],  [29],
        [20], [$361 times 361$],    num[13357],  num[10.2],  [37],
        [24], [$529 times 529$],    num[23805],  num[8.5],   [45],
        [32], [$961 times 961$],    num[58621],  num[6.3],   [61],
      )
    },
  ),
  caption: [Sparsity statistics for the 2D Laplacian matrix $L = I times.o D^2 + D^2 times.o I$. While the density decreases as $cal(O)(1\/N)$, each row has $2(N - 1) - 1$ nonzeros, growing linearly with $N$. For comparison, a standard second-order finite difference stencil on the same grid has exactly 5 nonzeros per row regardless of $N$. This denser coupling is the price paid for spectral (exponential) convergence.],
) <tbl-laplacian-sparsity>

#etude-conclusion[
  The 1D-to-2D step shows the power of the Kronecker construction $L = I times.o D^2 + D^2 times.o I$, which assembles the 2D Laplacian from two copies of the 1D second-derivative matrix *without any new algorithmic ideas*. The resulting $(N - 1)^2 times (N - 1)^2$ system inherits the spectral accuracy of its 1D building blocks: at $N = 16$ the maximum error is $approx 10^(-12)$, matching the 1D rate. The price is *cost scaling*: a direct dense solve is $cal(O)(N^6)$ in 2D. The matrix has $2 (N - 1) - 1$ nonzeros per row (linear in $N$), in contrast to the fixed-five 5-point stencil of second-order FD. This is the cost of dense global coupling. For moderate $N$ (up to $approx 32$) the direct approach is fine; for larger problems exploit the tensor-product structure (e.g.\ matrix diagonalisation) to drop to $cal(O)(N^3)$.
]

The code generating @fig-poisson-2d is available in:
- `codes/python/ch08/bvp_2d_poisson.py`
- `codes/matlab/ch08/bvp_2d_poisson.m`
- `codes/julia/ch08/bvp_2d_poisson.jl`

== Computational Étude 8.7: The Helmholtz Equation <sec-helmholtz>

=== Near-Resonance Behavior

The Helmholtz equation#idx("Helmholtz equation") models wave phenomena:
$ u_(x x) + u_(y y) + k^2 u = f(x, y), quad u = 0 "on boundary." $ <eq-helmholtz>

When $k^2$ approaches an eigenvalue of the Laplacian, the system becomes _nearly resonant_ and the solution amplitude grows dramatically. The near-singularity of the Helmholtz operator near resonance is a manifestation of the _pollution effect_ analysed by Babuska and Sauter @BabuskaSauter1997, who showed that maintaining accuracy at high wavenumbers requires $N tilde k^(3\/2)$ or even $N tilde k^2$ rather than the naïve $N approx k\/pi$. Moiola and Spence @MoiolaSpence2014 further clarified the sign-indefinite nature of the Helmholtz operator, explaining why standard iterative solvers stall near resonance.

For a localized Gaussian forcing $f(x, y) = e^(-20[(x - 0.3)^2 + (y + 0.4)^2])$ and $k = 7$, we are near resonance with the $(2, 4)$ mode (theoretical $k approx 7.02$).

@fig-helmholtz shows the characteristic modal structure of the near-resonant solution.

#figure(
  image("../figures/ch08/python/helmholtz.pdf", width: 95%),
  caption: [Helmholtz equation @eq-helmholtz with $k = 7$, near resonance with the $(2, 4)$ eigenmode ($k_(2,4) approx 7.02$). Left: 3D surface showing the wave-like structure. Right: contour plot with forcing location marked. The solution exhibits the characteristic pattern of the $(2, 4)$ eigenmode.],
) <fig-helmholtz>

The Helmholtz solver modifies the 2D Laplacian by adding the $k^2 I$ term:

```python
def solve_helmholtz(N, k, f):
    """Solve ∇²u + k²u = f on [-1,1]² with u = 0 on boundary."""
    D, x = cheb_matrix(N)
    D2 = D @ D
    D2_int = D2[1:N, 1:N]
    x_int = x[1:N]

    # Build Helmholtz operator: L = I ⊗ D² + D² ⊗ I + k²I
    I = np.eye(N - 1)
    n_int = (N - 1)**2
    L = np.kron(I, D2_int) + np.kron(D2_int, I) + k**2 * np.eye(n_int)

    # Right-hand side
    X, Y = np.meshgrid(x_int, x_int)
    F = f(X, Y).flatten(order='F')

    # Solve
    u_vec = np.linalg.solve(L, F)
    U_int = u_vec.reshape((N-1, N-1), order='F')

    U = np.zeros((N+1, N+1))
    U[1:N, 1:N] = U_int
    return np.meshgrid(x, x), U
```

```matlab
function [grids, U] = solve_helmholtz(N, k, f)
% Solve ∇²u + k²u = f on [-1,1]² with u = 0 on boundary.
    [D, x] = cheb_matrix(N);
    D2 = D * D;
    D2_int = D2(2:N, 2:N);
    x_int = x(2:N);

    % Build Helmholtz operator: L = I ⊗ D² + D² ⊗ I + k²I
    I = eye(N - 1);
    n_int = (N - 1)^2;
    L = kron(I, D2_int) + kron(D2_int, I) + k^2 * eye(n_int);

    % Right-hand side
    [X, Y] = meshgrid(x_int, x_int);
    F = f(X, Y);

    % Solve
    u_vec = L \ F(:);
    U_int = reshape(u_vec, N-1, N-1);

    U = zeros(N+1, N+1);
    U(2:N, 2:N) = U_int;
    [grids{1}, grids{2}] = meshgrid(x, x);
end
```

The Julia implementation:

```julia
function solve_helmholtz(N, k, f)
    # Solve ∇²u + k²u = f on [-1,1]² with u = 0 on boundary.
    D, x = cheb_matrix(N)
    D2 = D * D
    D2_int = D2[2:N, 2:N]
    x_int = x[2:N]
    n_int = N - 1

    # Build Helmholtz operator: L = I ⊗ D² + D² ⊗ I + k²I
    Imat = I(n_int)
    L = kron(Imat, D2_int) + kron(D2_int, Imat) + k^2 * I(n_int^2)

    # Right-hand side
    F = [f(x_int[j], x_int[i]) for i in 1:n_int, j in 1:n_int]

    # Solve
    u_vec = L \ vec(F)
    U_int = reshape(u_vec, n_int, n_int)

    U = zeros(N + 1, N + 1)
    U[2:N, 2:N] = U_int
    return x, U
end
```

#etude-conclusion[
  When the wavenumber $k = 7$ sits close to the eigenvalue $k_(2, 4) approx 7.02$ of the Laplacian, the operator $nabla^2 + k^2 I$ becomes *nearly singular*, and even a modest localised Gaussian source produces a response with amplitude far exceeding the source itself. The solution locks onto the spatial pattern of the $(2, 4)$ eigenmode --- the one whose eigenvalue is nearest to $-k^2$. The near-resonance#idx("near-resonance") regime is a *stringent test* of numerical accuracy: low-order methods need very fine meshes to resolve the large-amplitude oscillatory structure, and the ill-conditioning further degrades accuracy. The spectral method handles this gracefully ($N = 24$ suffices) because the Chebyshev basis represents the smooth eigenmode structure with few points and the direct solver sidesteps the iterative-solver difficulties for sign-indefinite systems. *Robustness near resonance* is a quiet but real advantage of spectral methods.
]

The code generating @fig-helmholtz is available in:
- `codes/python/ch08/bvp_helmholtz.py`
- `codes/matlab/ch08/bvp_helmholtz.m`
- `codes/julia/ch08/bvp_helmholtz.jl`

== Computational Étude 8.8: The Quantum Harmonic Oscillator Revisited <sec-harmonic-oscillator>

In @sec-harmonic-oscillator-fourier, we solved the quantum harmonic oscillator eigenvalue problem
$ -u'' + x^2 u = lambda u, quad x in RR $ <eq-harmonic-oscillator>
using the periodic (Fourier) spectral method. We now revisit this same problem using the Chebyshev spectral methods developed in this chapter, and compare the two approaches. The comparison of Fourier (equispaced) and Chebyshev (clustered) grids for problems on unbounded domains is a recurring theme in spectral methods; Boyd @Boyd2000 discusses the interplay between domain truncation and basis choice, noting that when eigenfunctions decay rapidly, Fourier grids can outperform Chebyshev grids despite the latter's theoretical optimality on bounded domains.

=== Chebyshev Approach

Recall that the exact eigenvalues are $lambda_n = 2n + 1$ and the eigenfunctions are Hermite functions that decay like $e^(-x^2\/2)$. We truncate to $[-L, L]$ with homogeneous Dirichlet conditions $u(plus.minus L) = 0$, which are automatically satisfied by the true eigenfunctions for large $L$.

The Chebyshev method proceeds as follows:

1. Compute the Chebyshev--Gauss--Lobatto points $t_j = cos(j pi \/ N)$ on $[-1, 1]$ and the differentiation matrix $D_N$.
2. Map to the physical domain: $x_j = L t_j$, giving points on $[-L, L]$.
3. Form the second-derivative matrix $D^((2)) = D_N^2 \/ L^2$ (the factor $1\/L^2$ accounts for the chain rule).
4. Strip boundary rows and columns (matrix stripping, as in @sec-boundary-conditions) to obtain the interior operator.
5. Assemble $A = -D^((2))_"int" + "diag"(x_1^2, dots, x_(N-1)^2)$.
6. Compute eigenvalues of $A$.

In Python, using the Chebyshev matrix from @sec-cheb-matrix:

```python
def solve_chebyshev(N, L):
    """Solve -u'' + x²u = λu on [-L, L] using Chebyshev."""
    D, t = cheb_matrix(N)
    x = L * t                      # Map to [-L, L]
    D2 = (D @ D) / L**2            # Second derivative

    # Strip boundaries (matrix stripping)
    D2_int = D2[1:N, 1:N]
    x_int = x[1:N]

    A = -D2_int + np.diag(x_int**2)
    eigenvalues, eigenvectors = np.linalg.eig(A)
    idx = np.argsort(np.real(eigenvalues))
    return np.real(eigenvalues[idx]), np.real(eigenvectors[:, idx]), x
```

The equivalent MATLAB implementation:

```matlab
function [eigenvalues, eigvecs, x] = solve_chebyshev(N, L)
    [D, t] = cheb_matrix(N);
    x = L * t;                     % Map to [-L, L]
    D2 = (D * D) / L^2;           % Second derivative

    % Strip boundaries
    D2_int = D2(2:N, 2:N);
    x_int = x(2:N);

    A = -D2_int + diag(x_int.^2);
    [eigvecs, E] = eig(A);
    eigenvalues = sort(diag(E));
end
```

And in Julia:

```julia
function solve_chebyshev(N, L)
    D, t = cheb_matrix(N)
    x = L * t                      # Map to [-L, L]
    D2 = (D * D) / L^2            # Second derivative

    # Strip boundaries
    D2_int = D2[2:N, 2:N]
    x_int = x[2:N]

    A = -D2_int + Diagonal(x_int .^ 2)
    eig_result = eigen(Matrix(A))
    return sort(real.(eig_result.values)), eig_result.vectors, x
end
```

=== Results

@fig-harmonic-oscillator shows the Chebyshev method in action. The left panel displays the first five eigenfunctions computed with $N = 64$ Chebyshev points on $[-8, 8]$, overlaid with the exact Hermite functions. The agreement is visually perfect.

#figure(
  image("../figures/ch08/python/harmonic_oscillator.pdf", width: 95%),
  caption: [The quantum harmonic oscillator solved with Chebyshev spectral methods. _Left_: First five eigenfunctions computed with $N = 64$ Chebyshev--Gauss--Lobatto points on $[-8, 8]$ (dots) compared to exact Hermite functions (lines). _Right_: Eigenvalue convergence comparison between the Fourier and Chebyshev methods for $lambda_0$ and $lambda_1$. Both methods achieve spectral convergence, but the Fourier approach reaches machine precision with fewer points.],
) <fig-harmonic-oscillator>

=== Fourier vs Chebyshev: A Comparison

The right panel of @fig-harmonic-oscillator compares the convergence of the Fourier and Chebyshev methods on the same problem with $L = 8$. @tbl-harmonic-comparison presents the eigenvalue errors for both approaches.

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    {
      show table: format-table(auto, auto, auto, auto, auto)
      table(
        columns: (auto, 1fr, 1fr, 1fr, 1fr),
        align: (center, center, center, center, center),
        inset: (x: 1em, y: 0.6em),
        stroke: none,
        table.hline(stroke: 0.75pt + rgb("#142D6E")),
        table.header(
          table.cell(fill: rgb("#142D6E").lighten(85%))[*$N$*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*Fourier $lambda_0$*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*Chebyshev $lambda_0$*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*Fourier $lambda_1$*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*Chebyshev $lambda_1$*],
        ),
        table.hline(stroke: 0.5pt + luma(180)),
        [12], num[2.19e-2], num[3.39e-1], num[1.72e-1], num[1.96e0],
        [18], num[3.00e-5], num[3.87e-2], num[6.44e-4], num[2.54e-1],
        [24], num[2.37e-9], num[1.20e-3], num[9.84e-8], num[1.49e-2],
        [30], num[1.77e-14], num[1.44e-5], num[7.54e-13], num[2.95e-4],
        [36], num[6.00e-15], num[7.38e-8], num[6.66e-15], num[2.17e-6],
        [42], num[2.11e-14], num[1.75e-10], num[2.80e-14], num[6.82e-9],
        [48], num[1.72e-14], num[2.01e-13], num[1.02e-14], num[1.01e-11],
      )
    },
  ),
  caption: [Eigenvalue errors for the quantum harmonic oscillator: Fourier vs Chebyshev methods with $L = 8$. Both methods achieve machine precision, but the Fourier method converges faster for this particular problem.],
) <tbl-harmonic-comparison>

The comparison reveals an interesting phenomenon: the Fourier method converges _faster_ than the Chebyshev method for this problem, despite being designed for periodic functions. This can be understood by considering the nature of the eigenfunctions and the grid point distributions:

- *Eigenfunctions decay rapidly.* The Hermite functions are effectively zero beyond $|x| approx 5$ for the lowest modes. On $[-8, 8]$, the solution is essentially zero near the boundaries, making it _effectively periodic_.

- *Equispaced points sample uniformly.* The Fourier method places points uniformly across $[-L, L]$, providing resolution where the eigenfunctions are non-negligible.

- *Chebyshev points cluster near boundaries.* The Chebyshev--Gauss--Lobatto points accumulate near $x = plus.minus L$, precisely where the eigenfunctions are vanishingly small. This "wastes" resolution on regions where there is nothing to resolve.

The Chebyshev method remains more general: it handles any non-periodic boundary value problem, including those where the solution has significant structure near the boundaries. The Fourier method exploits the specific structure of this problem. This comparison illustrates that the best numerical method depends on the problem at hand; matching the method to the solution structure can yield significant efficiency gains.

#etude-conclusion[
  The side-by-side comparison evaluates two spectral strategies on the *same* problem. Both reach machine precision, but Fourier needs $N = 30$ while Chebyshev needs $N approx 48$ --- a $60 %$ overhead traceable to Chebyshev's boundary clustering, which deploys resolution near $x = plus.minus L$ where the Hermite eigenstates have already decayed to zero. This is the same issue identified in Étude 5.3: the eigenfunctions are effectively periodic on a wide window because they vanish exponentially at both endpoints, so Fourier *matches the structure* and Chebyshev wastes some resolution. The Chebyshev method nevertheless remains the more general tool: for BVPs whose solutions have significant boundary structure (the generic case), the clustered points are exactly where they are needed. *Matching the method to the solution structure pays off; generality comes at a cost.*
]

The code generating @fig-harmonic-oscillator is available in:
- `codes/python/ch08/harmonic_oscillator.py`
- `codes/matlab/ch08/harmonic_oscillator.m`
- `codes/julia/ch08/harmonic_oscillator.jl`

== A non-exhaustive literature overview

The application of spectral methods to boundary value problems is a narrative of overcoming the rigidity of the Fast Fourier Transform to handle non-periodic constraints. The foundational shift from the integral-based Galerkin formulations of the early 20th century --- pioneered by Galerkin @Galerkin1915 and formalised by Lanczos @Lanczos1938 --- to the differential point-wise approach used in this chapter occurred in the 1970s. The seminal works of Kreiss and Oliger @KreissOliger1972 and Orszag @Orszag1971 established the viability of the collocation (pseudospectral) method, demonstrating that satisfying the differential equation exactly at grid points could yield exponential convergence without the computational burden of convolution sums. The construction of the differentiation matrix, the algebraic heart of this chapter, was systematised by Weideman and Reddy @WeidemanReddy2000 in their MATLAB Differentiation Matrix Suite, which provided not only the explicit algorithms for Chebyshev and Hermite matrices but also the rigorous error analysis justifying the "negative sum trick" to mitigate round-off errors. The debate between computing second derivatives via matrix squaring ($D^2$) versus explicit recurrence formulas was settled for practical purposes by Baltensperger and Berrut @BaltenspergerBerrut1999, who showed that for the moderate resolutions ($N < 200$) typical of spectral methods, the operational simplicity of matrix squaring outweighs theoretical stability concerns.

The technique of boundary bordering --- replacing rows of the differentiation matrix with boundary constraints --- was codified by Boyd @Boyd2000 in his influential text on Chebyshev and Fourier spectral methods. While alternative methods like basis recombination (constructing basis functions that inherently satisfy the boundary conditions) offer better conditioning for high-order problems, as demonstrated by Shen @Shen1995, boundary bordering remains the standard for general-purpose solvers due to its flexibility in handling Robin and mixed conditions. The chapter's "matrix surgery" approach --- stripping boundary rows and columns for Dirichlet problems --- follows the pedagogical tradition established by Gottlieb and Orszag @GottliebOrszag1977 and refined by Trefethen @Trefethen2000, who emphasised the elegance of reducing differential equations to small, dense linear systems solvable by direct methods.

For nonlinear problems, the Bratu equation serves as the historical benchmark. Originally derived by Frank-Kamenetskii @FrankKamenetskii1955 in the context of thermal ignition theory, this equation was given its definitive spectral treatment by Boyd @Boyd1986Bratu, who used Chebyshev spectral methods with arc-length continuation to resolve the turning point singularity and compute the critical parameter $lambda_c approx 3.513830719$ to high precision. In the contemporary era (2020--2026), the focus has shifted to fractional-order generalisations. Salama _et al._ @Salama2025 have extended spectral collocation using shifted Lucas polynomials to solve fractional Bratu equations, modelling non-local memory effects in combustion. Gheorghiu @Gheorghiu2020 revisited the classical problem to benchmark the accuracy of modern Newton--Kantorovich solvers, reconfirming the method's unrivalled efficiency for smooth nonlinearities. A comprehensive treatment of Newton methods for nonlinear PDEs, including convergence theory and affine invariance, can be found in the monograph by Deuflhard @Deuflhard2011.

In the realm of wave propagation, the Helmholtz equation reveals the subtle limitations of spectral accuracy. The "$pi$ points per wavelength" rule is necessary but not sufficient for high-frequency problems due to the _pollution effect_, a phenomenon rigorously analysed by Babuska and Sauter @BabuskaSauter1997. Recent work by Gkanos _et al._ @Gkanos2025 utilises dispersion analysis to quantify these errors in complex acoustical environments, while Galkowski and Spence @GalkowskiSpence2025 have applied semiclassical analysis to the high-frequency Helmholtz equation, providing theoretical upper bounds on the grid densities required to resolve scattering from non-convex obstacles. Moiola and Spence @MoiolaSpence2014 clarified the sign-indefinite nature of the Helmholtz operator, explaining why standard iterative solvers stall near resonance and motivating the development of complex-shifted Laplacian preconditioners.

For two-dimensional problems, the chapter's tensor product approach using Kronecker products is efficient ($O(N^3)$) but restricted to rectangular geometries. The mention of Padua points connects to a major breakthrough in approximation theory: De Marchi, Caliari, and Vianello @DeMarchi2005 discovered these as the first explicit unisolvent point set on the square with optimal Lebesgue constant growth ($O(log^2 N)$), later proven rigorously by Bos _et al._ @Bos2006. The elegant geometric construction via Lissajous curves @Bos2007 and the efficient Algorithm 886 @Caliari2008 make Padua points a practical alternative to tensor product grids when interpolation accuracy, rather than differential operator structure, is the primary concern.

Finally, the field is currently witnessing a paradigm shift towards quantum algorithms for spectral methods. Liu _et al._ @Liu2025quantum have recently proposed a Quantum Spectral Method for non-periodic BVPs, encoding Chebyshev basis coefficients into quantum states to achieve polylogarithmic computational complexity ($O("poly"(log N))$) for solving the linear systems. This suggests that the "matrix surgery" described in this chapter --- currently an $O(N^3)$ operation on classical silicon --- could one day become an $O("poly"(log N))$ operation on a quantum processor. The treatment of non-local boundary conditions, as explored by Mustapha _et al._ @Mustapha2025 for fractional differential equations, further extends the boundaries of what spectral methods can address, requiring operational matrices derived from shifted Jacobi or Lucas polynomials to handle integral constraints rather than local derivative conditions.

== Summary

This chapter has demonstrated the power and simplicity of spectral collocation for boundary value problems:

- *Matrix stripping* handles Dirichlet boundary conditions by removing boundary rows and columns.

- *Variable coefficients* become diagonal matrices---no additional complexity.

- *Nonlinear problems* yield to Newton iteration, with the Jacobian easily computed from the differentiation matrix.

- *Eigenvalue problems* reveal the resolution limits of spectral methods: $gt.eq.slant pi$ points per wavelength are needed.

- *Tensor products* extend the method to higher dimensions via Kronecker products.

The key insight throughout is that spectral methods transform differential equations into _dense linear algebra_. The matrices are small (because spectral accuracy requires few points), and the resulting systems can be solved directly. This directness, no iteration and no mesh refinement, is the hallmark of spectral methods for smooth problems.

== Exercises <sec-boundary-value-problems-exercises>

*Exercise 8.1* (_Spectral Collocation for a Variable-Coefficient ODE_). Solve the boundary value problem $-u'' + (1 + x^2) u = f(x)$ on $[-1, 1]$ with $u(-1) = u(1) = 0$, where $f(x)$ is chosen so that $u(x) = (1 - x^2) e^(-x^2)$ is the exact solution. (a) Implement Chebyshev spectral collocation with matrix stripping. (b) Plot the maximum error versus $N$ for $N = 8, 12, 16, 20, 24, 32$ on a semilogarithmic scale. (c) Verify exponential convergence and determine how many points are needed for 12-digit accuracy.

*Exercise 8.2* (_Newton Iteration for a Nonlinear BVP_). Solve the Bratu problem $-u'' = lambda e^u$ on $[0, 1]$ with $u(0) = u(1) = 0$, for $lambda = 1$. (a) Formulate the Newton iteration for the spectral collocation system. The Jacobian is $J = D_N^2 + lambda "diag"(e^(bold(u)))$. (b) Starting from $bold(u)^((0)) = bold(0)$, iterate until $||bold(u)^((k+1)) - bold(u)^((k))||_infinity < 10^(-12)$. Report the number of iterations required. (c) Plot the converged solution and verify its accuracy by substituting back into the ODE. (d) Investigate what happens as $lambda arrow lambda_c approx 3.51$: how many iterations does Newton's method require, and does it still converge?

*Exercise 8.3* (_Boundary Bordering vs Matrix Stripping_). For the BVP $-u'' = e^x$ on $[-1, 1]$ with $u(-1) = 0$, $u(1) = 1$: (a) Solve using matrix stripping (remove boundary rows/columns, modify RHS). (b) Solve using boundary bordering (replace first and last rows of $D_N^2$ with the boundary conditions). (c) Compare the condition numbers of the two resulting linear systems for $N = 8, 16, 32, 64$ and discuss which approach is better conditioned.

*Exercise 8.4* (_2D Poisson Equation with Tensor Products_). Solve $-Delta u = f(x, y)$ on $[-1, 1]^2$ with $u = 0$ on the boundary, where $f$ is chosen so that $u(x, y) = sin(pi x) sin(pi y)$ is the exact solution. (a) Construct the 2D Laplacian using Kronecker products: $L = D_2 times.o I + I times.o D_2$. (b) Apply boundary stripping to both dimensions. (c) Measure the maximum error for $N = 8, 12, 16, 20$ and verify exponential convergence. (d) Time the computation for $N = 16, 24, 32, 40$ and verify the $cal(O)(N^3)$ scaling of the direct solve.

*Exercise 8.5* (_Convergence Study: Semilogarithmic Plot_). Consider the eigenvalue problem $-u'' = lambda u$ on $[-1, 1]$ with $u(-1) = u(1) = 0$. The exact eigenvalues are $lambda_k = (k pi \/ 2)^2$. (a) Compute the first 10 numerical eigenvalues using spectral collocation for $N = 16, 32, 48, 64$. (b) For each $N$, plot the relative error $|lambda_k^("num") - lambda_k^("exact")| \/ lambda_k^("exact")$ versus $k$ on a semilogarithmic scale. (c) Determine the "resolution limit": for each $N$, what is the largest $k$ for which the relative error is less than $10^(-6)$? Verify the rule of thumb $k_(max) approx N \/ pi$.
