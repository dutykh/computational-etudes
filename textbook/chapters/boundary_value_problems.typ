// textbook/chapters/boundary_value_problems.typ
// Chapter 7: Boundary Value Problems
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Last modified: January 2026

#import "../styles/template.typ": dropcap

// Enable equation numbering for this chapter
#set math.equation(numbering: "(1)")

= Boundary Value Problems <ch-bvp>

#dropcap[With the Chebyshev differentiation matrix in hand, we are ready to tackle one of the most important classes of problems in applied mathematics: boundary value problems (BVPs). Unlike initial value problems, where we march forward in time from given initial conditions, BVPs impose constraints at multiple locations---typically at the boundaries of the domain. This spatial coupling makes BVPs inherently global, and spectral methods are ideally suited to exploit this structure.]

This chapter demonstrates how to transform differential equations into linear algebra problems. The differentiation matrix $D_N$ from @ch-chebyshev becomes the workhorse, and its square $D_N^2$ handles second-order equations. Imposing boundary conditions requires a simple "matrix surgery"---removing rows and columns corresponding to boundary points. The result is a systematic approach that handles linear, variable-coefficient, and even nonlinear problems with remarkable ease.

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

The result is an $(N-1) times (N-1)$ interior system:
$ tilde(D)^2 = D_N^2 [1 : N-1, 1 : N-1]. $

For homogeneous conditions ($alpha = beta = 0$), the boundary terms vanish and we simply solve $tilde(D)^2 bold(u)_("int") = bold(f)_("int")$.

== Linear BVP: The Poisson Equation <sec-poisson-1d>

=== A Model Problem

Consider the one-dimensional Poisson equation:
$ u_(x x) = sin(pi x) + 2 cos(2 pi x), quad x in (-1, 1), quad u(plus.minus 1) = 0. $ <eq-poisson-1d>

This has the exact solution
$ u(x) = -frac(sin(pi x), pi^2) + frac(1 - cos(2 pi x), 2 pi^2). $

=== Spectral Solution

The solution procedure is direct:
1. Construct $D_N$ and compute $D_N^2$.
2. Extract the interior submatrix $tilde(D)^2 = D_N^2 [1 : N-1, 1 : N-1]$.
3. Evaluate the right-hand side at interior points.
4. Solve the linear system $tilde(D)^2 bold(u)_("int") = bold(f)_("int")$.
5. Embed the result in the full vector with boundary values.

@fig-poisson-1d shows the solution and demonstrates spectral convergence.

#figure(
  image("../figures/ch08/python/poisson_1d.pdf", width: 95%),
  caption: [Solution of the 1D Poisson equation @eq-poisson-1d. Left: numerical solution (circles) compared with exact solution (line) for $N = 16$. Right: exponential convergence of the maximum error, reaching machine precision by $N = 24$.],
) <fig-poisson-1d>

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

== Nonlinear BVP: The Bratu Equation <sec-bratu>

=== A Classic Nonlinear Problem

The Bratu equation models combustion and thermal explosion:
$ u_(x x) + lambda e^u = 0, quad x in (-1, 1), quad u(plus.minus 1) = 0. $ <eq-bratu>

This equation exhibits a _turning point_ phenomenon: solutions exist only for $lambda lt.eq.slant lambda_c$, where $lambda_c approx 0.878$ for the domain $[-1, 1]$. Above this critical value, no solution exists. For $lambda = 0.5$ (well below the critical value), a unique solution exists.

The nonlinearity $e^u$ prevents direct linear algebra. Instead, we use _Newton iteration_: linearize, solve, update, repeat.

=== Newton Iteration

Newton's method for the discrete system $F(bold(u)) = D^2 bold(u) + lambda e^(bold(u)) = bold(0)$ proceeds as follows:

1. Compute the residual: $bold(F) = D^2 bold(u) + lambda e^(bold(u))$.
2. Compute the Jacobian: $J = D^2 + lambda "diag"(e^(bold(u)))$.
3. Solve for the Newton step: $J delta bold(u) = -bold(F)$.
4. Update: $bold(u) arrow.l bold(u) + delta bold(u)$.

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

== Eigenvalue Problems <sec-eigenvalue>

=== Resolution Limits

The eigenvalue problem for the Laplacian,
$ u_(x x) = lambda u, quad x in (-1, 1), quad u(plus.minus 1) = 0, $ <eq-eigenvalue>
has exact eigenvalues $lambda_n = -(n pi \/ 2)^2$ and eigenfunctions $u_n (x) = sin(n pi (x + 1) \/ 2)$.

Spectral methods compute these eigenvalues with spectral accuracy---_for the low modes_. High-frequency modes require many points per wavelength (ppw) for accurate representation. The rule of thumb is:
$ "Need" gt.eq.slant pi " points per wavelength for spectral accuracy." $

@fig-eigenvalue demonstrates this resolution limit.

#figure(
  image("../figures/ch08/python/eigenvalue_problem.pdf", width: 95%),
  caption: [Eigenvalue problem @eq-eigenvalue for $N = 36$. Six eigenmodes are shown with their points per wavelength (ppw) and eigenvalue errors. Low modes (high ppw) are resolved to near machine precision. As ppw drops below $pi approx 3.14$, accuracy degrades rapidly. This is not a bug but a fundamental resolution limit.],
) <fig-eigenvalue>

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

== Two-Dimensional Problems <sec-2d>

=== Tensor Products

For problems on a rectangle $[-1, 1]^2$, we use _tensor product_ grids: Chebyshev points in both $x$ and $y$ directions. The total number of grid points is $(N + 1)^2$.

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

The sparsity pattern of the 2D Laplacian operator, shown in @fig-laplacian-sparsity, reveals the Kronecker product structure.

#figure(
  image("../figures/ch08/python/laplacian_sparsity.pdf", width: 55%),
  caption: [Sparsity pattern of the 2D Laplacian matrix for $N = 16$. The $(N-1)^2 times (N-1)^2$ matrix shows the characteristic block structure arising from the Kronecker product formulation @eq-laplacian-2d.],
) <fig-laplacian-sparsity>

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

== Summary

This chapter has demonstrated the power and simplicity of spectral collocation for boundary value problems:

- *Matrix stripping* handles Dirichlet boundary conditions by removing boundary rows and columns.

- *Variable coefficients* become diagonal matrices---no additional complexity.

- *Nonlinear problems* yield to Newton iteration, with the Jacobian easily computed from the differentiation matrix.

- *Eigenvalue problems* reveal the resolution limits of spectral methods: $gt.eq.slant pi$ points per wavelength are needed.

- *Tensor products* extend the method to higher dimensions via Kronecker products.

The key insight throughout is that spectral methods transform differential equations into _dense linear algebra_. The matrices are small (because spectral accuracy requires few points), and the resulting systems can be solved directly. This directness---no iteration, no mesh refinement---is the hallmark of spectral methods for smooth problems.
