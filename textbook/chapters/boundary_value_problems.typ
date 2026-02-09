// textbook/chapters/boundary_value_problems.typ
// Chapter 7: Boundary Value Problems
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Last modified: January 2026

#import "../styles/template.typ": dropcap, num, format-table

// Enable equation numbering for this chapter
#set math.equation(numbering: "(1)")

= Boundary Value Problems <ch-bvp>

#dropcap[With the Chebyshev differentiation matrix in hand, we are ready to tackle one of the most important classes of problems in applied mathematics: boundary value problems (BVPs). Unlike initial value problems, where we march forward in time from given initial conditions, BVPs impose constraints at multiple locations, typically at the boundaries of the domain. This spatial coupling makes BVPs inherently global, and spectral methods are ideally suited to exploit this structure.]

This chapter demonstrates how to transform differential equations into linear algebra problems. The differentiation matrix $D_N$ from @ch-chebyshev becomes the workhorse, and its square $D_N^2$ handles second-order equations. Imposing boundary conditions requires a simple "matrix surgery", removing rows and columns corresponding to boundary points. The result is a systematic approach that handles linear, variable-coefficient, and even nonlinear problems with remarkable ease.

== Second Derivatives and Matrix Squaring <sec-second-deriv>

=== The Need for $D^2$

Most BVPs in physics involve second-order derivatives: the heat equation, the wave equation, the Poisson equation, and many others all feature $u_(x x)$. To apply spectral collocation, we need a second derivative matrix.

Two approaches present themselves:

1. *Direct formulas*: Derive explicit expressions for $(D^2)_(i j)$ analogous to the first derivative formulas in @ch-chebyshev. These exist but are complex.

2. *Matrix squaring*: Simply compute $D^2 = D times D$. This is $O(N^3)$ but entirely adequate for spectral $N$-values (typically $N < 200$).

We adopt the second approach for its simplicity. The product $D_N^2$ gives us the second derivative matrix: if $bold(v)$ contains function values at the Chebyshev points, then $D_N^2 bold(v)$ approximates the second derivative values.

== Imposing Boundary Conditions <sec-boundary-conditions>

=== Dirichlet Conditions: Matrix Stripping

The most common boundary conditions specify the function values at the endpoints:
$ u(-1) = alpha, quad u(1) = beta. $
These are _Dirichlet conditions_.

With Chebyshev points ordered as $x_0 = 1$, $x_1$, ..., $x_(N-1)$, $x_N = -1$, the boundary conditions fix $v_0 = beta$ and $v_N = alpha$. The differential equation need only be enforced at the _interior_ points $x_1, dots, x_(N-1)$.

The implementation is straightforward:
- *Remove* the first and last rows of the equation (we don't need the ODE at boundary points).
- *Remove* the first and last columns (the boundary values are known, not unknowns).

The result is an $(N-1) times (N-1)$ interior system. Denoting the entries of $D_N^2$ by $(D_N^2)_(i j)$ with $i, j = 0, 1, dots, N$, we define
$ tilde(D)^2_(i j) = (D_N^2)_(i j), quad i, j = 1, dots, N-1. $

For homogeneous conditions ($alpha = beta = 0$), the boundary terms vanish and we simply solve $tilde(D)^2 bold(u)_("int") = bold(f)_("int")$.

=== Neumann and Robin Conditions: Boundary Bordering <sec-boundary-bordering>

Matrix stripping works well for Dirichlet conditions, but _Neumann_ conditions ($u'("boundary") = g$) or _Robin_ conditions ($alpha u + beta u' = gamma$) require a different approach. Since the boundary condition involves the derivative, we cannot simply remove the boundary unknowns.

The _boundary bordering_ technique keeps the full $(N + 1) times (N + 1)$ system and _replaces_ boundary rows with the appropriate conditions:
- *Neumann* $u'(x_0) = g$: replace row $0$ of $D^2$ with row $0$ of $D_N$ (the first derivative matrix), and set the corresponding right-hand side entry to $g$.
- *Robin* $alpha u(x_0) + beta u'(x_0) = gamma$: replace row $0$ with $alpha bold(e)_0^top + beta D_N [0, :]$, where $bold(e)_0$ is the first unit vector.
- *Dirichlet* $u(x_0) = alpha$: replace row $0$ with $bold(e)_0^top$ and set the right-hand side to $alpha$.

This approach handles any combination of boundary conditions uniformly. We demonstrate it in @sec-mixed-bc below.

== Linear BVP: The Poisson Equation <sec-poisson-1d>

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

@tbl-poisson-convergence quantifies the convergence rate. The error decreases by roughly two orders of magnitude for each increment of two in $N$, the hallmark of spectral (exponential) convergence. Machine precision is reached by $N = 24$, after which further refinement cannot improve the result due to floating-point arithmetic limitations.

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

== Variable Coefficient Problems <sec-variable-coeff>

=== The Airy-Type Equation

Variable coefficients pose no additional difficulty for spectral methods. Consider:
$ u_(x x) - (1 + x^2) u = 1, quad x in (-1, 1), quad u(plus.minus 1) = 0. $ <eq-variable-coeff>

The variable coefficient $(1 + x^2)$ becomes a diagonal matrix. The discretized operator is:
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

== Mixed Boundary Conditions <sec-mixed-bc>

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

The complete source code is available in `codes/python/ch08/bvp_mixed_bc.py`, `codes/matlab/ch08/bvp_mixed_bc.m`, and `codes/julia/ch08/bvp_mixed_bc.jl`.

== Nonlinear BVP: The Bratu Equation <sec-bratu>

=== A Classic Nonlinear Problem

The Bratu equation models combustion and thermal explosion:
$ u_(x x) + lambda e^u = 0, quad x in (-1, 1), quad u(plus.minus 1) = 0. $ <eq-bratu>

This equation exhibits a _turning point_ phenomenon: solutions exist only for $lambda lt.eq.slant lambda_c$, where $lambda_c approx 0.878$ for the domain $[-1, 1]$. Above this critical value, no solution exists. For $lambda = 0.5$ (well below the critical value), a unique solution exists.

The nonlinearity $e^u$ prevents direct linear algebra. Instead, we use _Newton iteration_: linearize, solve, update, repeat.

=== Newton Iteration

To solve a nonlinear equation, we cannot simply invert a matrix. Instead, we employ an _iterative_ strategy: start from an initial guess, and progressively refine it until the solution is accurate enough.

Newton's method is the most widely used approach for nonlinear equations. The key idea is _linearization_: at each iteration, we approximate the nonlinear problem by a linear one, solve that linear problem exactly, and use the result to improve our current approximation. Concretely, suppose we have a current approximation $bold(u)^((k))$. We seek a correction $delta bold(u)$ such that $bold(u)^((k)) + delta bold(u)$ satisfies the equation more closely. Expanding the residual to first order:
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

== Eigenvalue Problems <sec-eigenvalue>

=== Resolution Limits

The eigenvalue problem for the Laplacian,
$ u_(x x) = lambda u, quad x in (-1, 1), quad u(plus.minus 1) = 0, $ <eq-eigenvalue>
has exact eigenvalues $lambda_n = -(n pi \/ 2)^2$ and eigenfunctions $u_n (x) = sin(n pi (x + 1) \/ 2)$.

Spectral methods compute these eigenvalues with spectral accuracy, _but only for the low modes_. High-frequency modes require many points per wavelength (ppw) for accurate representation. The rule of thumb is:
$ "Need" gt.eq.slant pi " points per wavelength for spectral accuracy." $

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

== Two-Dimensional Problems <sec-2d>

=== Tensor Products

For problems on a rectangle $[-1, 1]^2$, we use _tensor product_ grids: Chebyshev points in both $x$ and $y$ directions. The total number of grid points is $(N + 1)^2$.

#block(
  width: 100%,
  inset: (x: 1.2em, y: 0.9em),
  stroke: (left: 2pt + rgb("#142D6E").lighten(50%)),
  fill: rgb("#142D6E").lighten(95%),
)[
  *Remark (Padua points).* Tensor product grids are not the only option for bivariate polynomial interpolation on the square. The _Padua points_, discovered at the University of Padua by De Marchi, Caliari, and Vianello (2005), are the first known example (and to date the only one) of a unisolvent point set on $[-1, 1]^2$ with _minimal growth_ of the Lebesgue constant, proven to be $cal(O)(log^2 n)$ by Bos, Caliari, De Marchi, Vianello, and Xu (2006). For total polynomial degree $n$, Padua points require only $(n + 1)(n + 2) \/ 2$ points, roughly half the $(n + 1)^2$ needed by the tensor product grid. The points admit an elegant geometric construction: they lie exactly on the self-intersections and boundary contacts of a generating Lissajous curve on $[-1, 1]^2$, which gives rise to four distinct families obtained by successive 90-degree rotations. While we use tensor product grids throughout this chapter for their simplicity and natural compatibility with the Kronecker product structure of differential operators, Padua points offer a more economical alternative for pure interpolation and approximation problems.
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

== The Helmholtz Equation <sec-helmholtz>

=== Near-Resonance Behavior

The Helmholtz equation models wave phenomena:
$ u_(x x) + u_(y y) + k^2 u = f(x, y), quad u = 0 "on boundary." $ <eq-helmholtz>

When $k^2$ approaches an eigenvalue of the Laplacian, the system becomes _nearly resonant_ and the solution amplitude grows dramatically.

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

== Computational Étude: The Quantum Harmonic Oscillator Revisited <sec-harmonic-oscillator>

In @sec-harmonic-oscillator-fourier, we solved the quantum harmonic oscillator eigenvalue problem
$ -u'' + x^2 u = lambda u, quad x in RR $ <eq-harmonic-oscillator>
using the periodic (Fourier) spectral method. We now revisit this same problem using the Chebyshev spectral methods developed in this chapter, and compare the two approaches.

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

The code generating @fig-harmonic-oscillator is available in:
- `codes/python/ch08/harmonic_oscillator.py`
- `codes/matlab/ch08/harmonic_oscillator.m`
- `codes/julia/ch08/harmonic_oscillator.jl`

== Summary

This chapter has demonstrated the power and simplicity of spectral collocation for boundary value problems:

- *Matrix stripping* handles Dirichlet boundary conditions by removing boundary rows and columns.

- *Variable coefficients* become diagonal matrices---no additional complexity.

- *Nonlinear problems* yield to Newton iteration, with the Jacobian easily computed from the differentiation matrix.

- *Eigenvalue problems* reveal the resolution limits of spectral methods: $gt.eq.slant pi$ points per wavelength are needed.

- *Tensor products* extend the method to higher dimensions via Kronecker products.

The key insight throughout is that spectral methods transform differential equations into _dense linear algebra_. The matrices are small (because spectral accuracy requires few points), and the resulting systems can be solved directly. This directness, no iteration and no mesh refinement, is the hallmark of spectral methods for smooth problems.
