// textbook/chapters/advanced_boundary_conditions.typ
// Chapter 13: Advanced Boundary Conditions for Spectral Methods
// Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
// Email: denys.dutykh@ku.ac.ae
// Homepage: https://www.denys-dutykh.com/
// Last modified: February 2026

#import "../styles/template.typ": dropcap, num, format-table

// Enable equation numbering for this chapter
#set math.equation(numbering: "(1)")

= Advanced Boundary Conditions for Spectral Methods <ch-advanced-bc>

#dropcap[The spectral machinery assembled in the preceding chapters has been deployed on problems with a single, simple type of boundary condition: homogeneous Dirichlet conditions $u = 0$ at the endpoints, imposed by stripping the first and last rows and columns of the differentiation matrix. This "matrix surgery" is elegant and efficient, but it covers only a narrow slice of the boundary conditions encountered in practice. Inhomogeneous Dirichlet data, Neumann and Robin conditions, time-dependent boundaries, nonlinear boundary constraints, and --- most strikingly --- boundary conditions that depend on the unknown eigenvalue all arise in applications ranging from heat transfer to black hole physics. In this chapter, we develop two systematic strategies for handling these diverse boundary conditions within the Chebyshev spectral framework, and we demonstrate them through a sequence of increasingly sophisticated computational études.]

By the end of this chapter, you should be able to:

1. Distinguish between _Method I_ (restricting interpolants to satisfy boundary conditions) and _Method II_ (replacing equations to enforce boundary conditions), and know when each is appropriate.
2. Apply the _lifting technique_ to reduce inhomogeneous Dirichlet conditions to homogeneous ones, including the case of time-dependent boundary data.
3. Impose Neumann, Robin, and mixed boundary conditions via the _tau approach_ (boundary row replacement) for both static and time-dependent problems.
4. Handle _nonlinear boundary conditions_ by coupling Newton iteration with the spectral discretisation.
5. Extend boundary condition imposition to _two-dimensional_ problems using Kronecker products.
6. Formulate and solve _quadratic eigenvalue problems_ arising from frequency-dependent boundary conditions, with application to the quasinormal modes of black holes.

== Two Strategies for Boundary Conditions <sec-two-strategies>

The literature on spectral methods identifies two fundamentally different philosophies for imposing boundary conditions @Trefethen2000 @Boyd2000 @Canuto2006 @GottliebOrszag1977. Understanding the mathematical distinction between these philosophies is essential for ensuring stability and achieving optimal complexity in advanced physical simulations.

=== Method I: Restricting the Interpolant

The first strategy modifies the solution space so that the interpolating polynomial _automatically_ satisfies the boundary conditions. In the simplest case --- homogeneous Dirichlet conditions $u(plus.minus 1) = 0$ --- this amounts to stripping the boundary rows and columns from the differentiation matrix, as we did throughout @ch-bvp. The resulting $(N - 1) times (N - 1)$ interior system operates on a subspace of polynomials that vanish at the endpoints.

For inhomogeneous Dirichlet conditions $u(a) = alpha$, $u(b) = beta$, the restriction is achieved through _lifting_: one defines a simple, continuous function $w(x)$ satisfying the boundary data and solves for the remainder $v = u - w$ subject to homogeneous conditions. The original solution is then reconstructed as $u = v + w$. This technique has proven highly adaptable, finding extended, state-of-the-art use in fractional diffusion models where classical trace theorems fail, as well as in global lifting formulations for reduced-order parameter-dependent modelling @Mustapha2025 @Ervin2020.

*Advantages*: Method I produces a smaller system (the boundary unknowns are eliminated), and the solution automatically satisfies the boundary conditions to machine precision.

*Limitations*: Method I requires an explicit lifting function, which may be difficult or impossible to construct when the boundary condition depends on the unknown (as in eigenvalue problems with frequency-dependent boundary data) or when the condition involves derivatives.

=== Method II: Replacing Equations (the Tau Approach)

The second strategy keeps the full $(N + 1) times (N + 1)$ system and _replaces_ the boundary rows of the differential operator with equations encoding the boundary conditions. This is a collocation variant of the classical _tau method_ introduced by Lanczos @Lanczos1938 in the context of spectral Galerkin approximations. This boundary bordering technique has since been rigorously codified for arbitrary grid distributions @Fornberg1998, and remains the premier method for discretising time-delay systems and models featuring eigenvalue-dependent boundary operators @Provoost2024.

For a Dirichlet condition $u(x_0) = alpha$ at the first Chebyshev point $x_0 = 1$, row $0$ of the operator matrix $L$ is replaced by the unit vector $bold(e)_0^top$, and the corresponding right-hand side entry is set to $alpha$. For a Neumann condition $u'(x_0) = g$, row $0$ is replaced by the first row of the differentiation matrix $D_N [0, :]$, with right-hand side $g$. For a general linear boundary condition
$ alpha_0 u(x_b) + alpha_1 u'(x_b) = gamma, $ <eq-general-linear-bc>
the replacement row is $alpha_0 bold(e)_b^top + alpha_1 D_N [b, :]$, where $bold(e)_b$ selects the boundary node. This was already previewed in @sec-boundary-bordering of @ch-bvp.

*Advantages*: Method II handles any linear (and, with Newton iteration, nonlinear) boundary condition uniformly. It extends naturally to boundary conditions that depend on the eigenvalue parameter, making it indispensable for certain spectral eigenvalue problems.

*Limitations*: The full $(N + 1) times (N + 1)$ system is slightly larger, and the boundary rows may affect the condition number of the operator matrix.

=== When to Use Which

A useful rule of thumb:
- *Method I* is ideal for problems with standard Dirichlet data where a lifting function is easy to write down. It produces a cleaner, smaller system.
- *Method II* is necessary whenever the boundary condition involves derivatives (Neumann, Robin), depends on unknown parameters (eigenvalue problems), or changes in time in a way that couples to the solution.

In practice, many problems can be solved by either method. The computational études in this chapter illustrate both approaches and highlight the situations where Method II becomes essential.

== Inhomogeneous Dirichlet Conditions and Lifting <sec-lifting>

=== The Lifting Technique

Consider a second-order boundary value problem
$ cal(L) u = f(x), quad x in (a, b), quad u(a) = alpha, quad u(b) = beta, $ <eq-inhom-bvp>
where $cal(L)$ is a linear differential operator. The boundary data $alpha$ and $beta$ may be nonzero, so the homogeneous matrix-stripping approach from @ch-bvp does not directly apply.

The lifting technique proceeds in three steps:
1. *Choose a lifting function* $w(x)$ satisfying $w(a) = alpha$ and $w(b) = beta$. The simplest choice is the linear interpolant:
$ w(x) = alpha frac(b - x, b - a) + beta frac(x - a, b - a). $ <eq-lifting-function>
2. *Define the remainder* $v(x) = u(x) - w(x)$. Since $v(a) = 0$ and $v(b) = 0$, the function $v$ satisfies homogeneous Dirichlet conditions.
3. *Solve the modified problem* $cal(L) v = f(x) - cal(L) w(x)$ with homogeneous conditions, then reconstruct $u = v + w$.

The modified right-hand side $f - cal(L) w$ accounts for the effect of the lifting function on the differential equation. For the linear lifting function @eq-lifting-function, the second derivative $w''(x) = 0$, so only lower-order terms in $cal(L)$ contribute to the correction. For operators of the form $cal(L) u = u''$, the right-hand side is simply $f(x)$, unchanged.

When the boundary data is time-dependent, $alpha = alpha(t)$ and $beta = beta(t)$, the lifting function becomes $w(x, t)$ and an additional source term $-w_t$ appears in the modified equation. This is explored in Étude 13.3.

=== Computational Étude 13.1: Inhomogeneous Poisson with Lifting <etude-inhom-lifting>

We solve the inhomogeneous Poisson equation
$ u_(x x) = e^(4 x), quad -1 < x < 1, quad u(-1) = 0, quad u(1) = 1. $ <eq-inhom-poisson>
The lifting function is $w(x) = (x + 1) \/ 2$, which satisfies $w(-1) = 0$ and $w(1) = 1$. Since $w''(x) = 0$, the remainder $v = u - w$ satisfies
$ v_(x x) = e^(4 x), quad v(plus.minus 1) = 0, $
which is a standard homogeneous Dirichlet problem.

The implementation constructs the Chebyshev differentiation matrix, extracts the interior submatrix of $D^2$, solves for $v$, and reconstructs $u = v + w$. In Python:

```python
D, x = cheb_matrix(N)
D2 = D @ D
# Interior system (rows/cols 1..N-1)
D2_int = D2[1:N, 1:N]
f_int = np.exp(4 * x[1:N])
v_int = np.linalg.solve(D2_int, f_int)
# Reconstruct full solution with lifting
w = (x + 1) / 2
u = np.concatenate([[0], v_int, [0]]) + w
```

In MATLAB:

```matlab
[D, x] = cheb_matrix(N);
D2 = D^2;
D2_int = D2(2:N, 2:N);
f_int = exp(4 * x(2:N));
v_int = D2_int \ f_int;
w = (x + 1) / 2;
u = [0; v_int; 0] + w;
```

In Julia:

```julia
D, x = cheb_matrix(N)
D2 = D^2
D2_int = D2[2:N, 2:N]
f_int = exp.(4 .* x[2:N])
v_int = D2_int \ f_int
w = (x .+ 1) ./ 2
u = [0; v_int; 0] .+ w
```

#figure(
  image("../figures/ch13/python/inhom_lifting.pdf", width: 95%),
  caption: [Inhomogeneous Poisson equation @eq-inhom-poisson solved by lifting. _Left_: the numerical solution (circles) agrees with the reference solution for $N = 16$. _Right_: the maximum error decreases exponentially with $N$, reaching machine precision by $N approx 20$.],
) <fig-inhom-lifting>

The code generating @fig-inhom-lifting is available in:
- `codes/python/ch13/bc_inhom_lifting.py`
- `codes/matlab/ch13/bc_inhom_lifting.m`
- `codes/julia/ch13/bc_inhom_lifting.jl`

=== Discussion

The convergence plot in @fig-inhom-lifting illustrates a hallmark of Chebyshev spectral methods: the error decreases geometrically with $N$, reaching machine precision near $N = 20$ for this smooth problem. This exponential convergence rate stands in sharp contrast to finite-difference or finite-element methods, where the error decreases only algebraically. The lifting technique is both conceptually transparent and computationally trivial --- the linear interpolant $w(x) = (x + 1) \/ 2$ has vanishing second derivative, so the modified right-hand side is identical to the original. This simplicity, however, depends critically on the Dirichlet character of the boundary data: the lifting function must satisfy the boundary values exactly, and such a function is easy to construct only when the values themselves are prescribed.

The lifting technique is simple and effective, but it is limited to Dirichlet conditions where the boundary values are explicitly known. For derivative conditions, we need Method II.

== Neumann and Robin Conditions via the Tau Method <sec-neumann-robin>

=== Neumann Conditions: Derivative Rows

A Neumann boundary condition specifies the _derivative_ of the solution at a boundary point:
$ u'(x_b) = g. $
Since the boundary value $u(x_b)$ is unknown, we cannot strip it from the system. Instead, we keep the full $(N + 1) times (N + 1)$ operator matrix and replace the boundary row with the corresponding row of the first derivative matrix $D_N$.

Concretely, for the problem $u'' = f(x)$ with $u'(-1) = g_L$ and $u(1) = g_R$, the discrete system is
$ A bold(u) = bold(b), $
where $A$ is constructed from $D_N^2$ with two modifications:
- Row $0$ (corresponding to $x_0 = 1$) is replaced by $bold(e)_0^top$ (Dirichlet), with $b_0 = g_R$.
- Row $N$ (corresponding to $x_N = -1$) is replaced by $D_N [N, :]$ (Neumann), with $b_N = g_L$.

=== Robin Conditions: Combining Evaluation and Derivative Rows

A Robin (or impedance) boundary condition is a linear combination of the function value and its derivative:
$ alpha u(x_b) + beta u'(x_b) = gamma. $ <eq-robin-bc>
The corresponding replacement row is $alpha bold(e)_b^top + beta D_N [b, :]$, with right-hand side $gamma$.

Robin conditions arise naturally in heat transfer (convective boundary), acoustics (impedance boundary), and electromagnetics (absorbing boundary). Recent rigorous mathematical analyses demonstrate that Robin formulations are absolutely essential for capturing bilateral heat exchange in coupled fluid-solid environments, as well as for modelling frequency-dependent local reacting surfaces in high-fidelity acoustic simulations @Clarte2020 @Vorlander2022. The parameter ratio $alpha \/ beta$ controls the character of the condition: when $beta = 0$, it reduces to Dirichlet; when $alpha = 0$, to Neumann. Thus Robin conditions interpolate smoothly between the two classical types.

=== Computational Étude 13.2: Helmholtz Equation with Robin Boundaries <etude-helmholtz-robin>

We solve the Helmholtz equation
$ -u_(x x) + k^2 u = f(x), quad -1 < x < 1, $ <eq-helmholtz-robin>
with a Robin condition at the left boundary and a Dirichlet condition at the right:
$ u'(-1) - alpha u(-1) = 0, quad u(1) = 0. $ <eq-helmholtz-robin-bc>
The parameter $alpha gt.eq.slant 0$ controls the boundary character: $alpha = 0$ gives a pure Neumann condition (zero flux), while $alpha arrow infinity$ approaches Dirichlet (zero value). We choose $k = 2$ and $f(x) = e^(- x^2)$ as a smooth Gaussian source.

The operator matrix is $L = -D_N^2 + k^2 I$. The tau modifications are:
- Row $0$ ($x = 1$): replace with $bold(e)_0^top$, set $b_0 = 0$ (Dirichlet).
- Row $N$ ($x = -1$): replace with $D_N [N, :] - alpha bold(e)_N^top$, set $b_N = 0$ (Robin).

In Python:

```python
D, x = cheb_matrix(N)
D2 = D @ D
k = 2.0
L = -D2 + k**2 * np.eye(N + 1)
f = np.exp(-x**2)
# Dirichlet at x = 1 (row 0)
L[0, :] = 0; L[0, 0] = 1; f[0] = 0
# Robin at x = -1 (row N)
L[N, :] = D[N, :] - alpha * np.eye(N + 1)[N, :]
f[N] = 0
u = np.linalg.solve(L, f)
```

In MATLAB:

```matlab
[D, x] = cheb_matrix(N);
D2 = D^2;
k = 2;
L = -D2 + k^2 * eye(N + 1);
f = exp(-x.^2);
% Dirichlet at x = 1 (row 1)
L(1, :) = 0; L(1, 1) = 1; f(1) = 0;
% Robin at x = -1 (row N+1)
L(N+1, :) = D(N+1, :) - alpha * eye(N+1, N+1);
L(N+1, :) = L(N+1, :);  % Robin row
f(N+1) = 0;
u = L \ f;
```

In Julia:

```julia
D, x = cheb_matrix(N)
D2 = D^2
k = 2.0
L = -D2 + k^2 * I(N + 1)
f = exp.(-x.^2)
# Dirichlet at x = 1 (row 1)
L[1, :] .= 0; L[1, 1] = 1; f[1] = 0
# Robin at x = -1 (row N+1)
L[N+1, :] = D[N+1, :] - α * I(N+1)[N+1, :]
f[N+1] = 0
u = L \ f
```

@fig-helmholtz-robin shows the solution for several values of $alpha$. As $alpha$ increases from $0$ (Neumann) to large values, the solution transitions from having a nonzero value at $x = -1$ with zero derivative, to being clamped near zero at both boundaries. @fig-robin-conditioning shows that the condition number of the discrete operator depends on $alpha$ but remains moderate for practical values.

#figure(
  image("../figures/ch13/python/helmholtz_robin.pdf", width: 95%),
  caption: [Helmholtz equation @eq-helmholtz-robin with Robin boundary condition @eq-helmholtz-robin-bc for $k = 2$ and $N = 32$. _Left_: solutions for $alpha = 0$ (Neumann), $alpha = 1$, $alpha = 5$, and $alpha = 20$. As $alpha$ increases, the solution is progressively suppressed at $x = -1$. _Right_: spectral convergence of the maximum error for $alpha = 5$.],
) <fig-helmholtz-robin>

#figure(
  image("../figures/ch13/python/robin_conditioning.pdf", width: 95%),
  caption: [Condition number of the discrete Helmholtz operator with Robin boundary condition as a function of $alpha$ for $N = 16$, $32$, and $64$. The condition number is moderate for all practical values of $alpha$ and grows polynomially with $N$, as expected for spectral differentiation matrices.],
) <fig-robin-conditioning>

The code generating @fig-helmholtz-robin and @fig-robin-conditioning is available in:
- `codes/python/ch13/bc_helmholtz_robin.py`
- `codes/matlab/ch13/bc_helmholtz_robin.m`
- `codes/julia/ch13/bc_helmholtz_robin.jl`

=== Discussion

The family of solutions displayed in @fig-helmholtz-robin reveals how Robin conditions interpolate continuously between the Neumann and Dirichlet limits. For $alpha = 0$, the solution attains a nonzero value at $x = -1$ with zero slope, whereas for large $alpha$ the solution is increasingly pinned near zero at the left boundary. This smooth parametric transition is difficult to appreciate from the analytical theory alone, where the two limiting cases are typically treated as distinct problems with different solution structures. The computation makes the interpolation explicit and quantitative.

@fig-robin-conditioning addresses a practical concern that often arises when boundary rows of the operator matrix are replaced: the worry that the resulting system may be ill-conditioned. The data show that the condition number grows polynomially with $N$ (as expected for Chebyshev differentiation matrices) and depends only mildly on $alpha$, confirming that the tau approach is numerically robust for Robin conditions across a wide range of parameters.

== Allen--Cahn Equation: Metastability and Boundary Effects <sec-allen-cahn-bc>

The Allen--Cahn equation @AllenCahn1979
$ u_t = epsilon u_(x x) + u - u^3, quad -1 < x < 1, $ <eq-allen-cahn-bc>
models phase separation in binary alloys and provides a rich testbed for boundary condition techniques. The nonlinear reaction term $u - u^3$ has three equilibria: $u = -1$ and $u = 1$ (stable) and $u = 0$ (unstable). Solutions develop sharp transition layers (interfaces) between the two stable phases, and these interfaces exhibit _metastability_ --- they persist for exponentially long times before suddenly annihilating @Trefethen2000. The asymptotic projection methods detailing this exponentially slow drift, scaling as $cal(O)(e^(-C\/epsilon))$, were rigorously established by Carr and Pego and Fusco and Hale @CarrPego1989 @Alikakos1991. Contemporary spectral methods are now actively exploring how non-local fractional Laplacian operators fundamentally alter this time-to-collision, enhancing the metastability of the phase boundaries.

=== Metastability under Fixed Dirichlet Conditions

With fixed boundary conditions $u(-1, t) = -1$ and $u(1, t) = 1$, we apply Method I. The lifting function $w(x) = x$ satisfies the boundary data, and the remainder $v(x, t) = u(x, t) - x$ satisfies
$ v_t = epsilon (v_(x x) + 0) + (v + x) - (v + x)^3, quad v(plus.minus 1, t) = 0. $
The spatial operator is discretised with the interior Chebyshev system (matrix stripping), and time integration is performed with a simple explicit scheme. The initial condition $u(x, 0) = 0.53 x + 0.47 sin(-1.5 pi x)$ creates three internal transition layers that persist until approximately $t approx 45$, when two interfaces annihilate simultaneously.

=== Time-Dependent Boundaries and Interface Dynamics

When the boundary data varies in time, Method I requires updating the lifting function at each time step, introducing additional source terms. Method II offers a simpler alternative: at each time step, simply _overwrite_ the boundary values of the solution vector.

With time-dependent conditions $u(-1, t) = -1$ and $u(1, t) = 1 + sin^2(t \/ 5)$, the oscillating right boundary pumps energy into the system, shifting the equilibrium interface location from $x = 0$ to approximately $x approx -0.4$ and causing the interfaces to annihilate earlier, near $t approx 30$.

=== Computational Étude 13.3: Allen--Cahn with Fixed and Driven Boundaries <etude-allen-cahn-bc>

We implement both scenarios using Chebyshev spatial discretisation with $N = 80$ and $epsilon = 0.01$.

*Part A (fixed boundaries, Method I)*: The lifting function $w(x) = x$ reduces the problem to homogeneous Dirichlet conditions on $v$. We march in time with a fourth-order exponential time differencing scheme (ETDRK4) or, more simply, with the backward Euler / forward Euler (IMEX) splitting: treat $epsilon u_(x x)$ implicitly and $u - u^3$ explicitly.

*Part B (driven boundaries, Method II)*: The full system including boundaries is advanced in time, and at the end of each time step the boundary values are overwritten:
$ u_0 = 1 + sin^2(t \/ 5), quad u_N = -1. $

In Python (Part A sketch):

```python
D, x = cheb_matrix(N)
D2 = D @ D
D2_int = D2[1:N, 1:N]
w = x.copy()  # lifting function
v = u0 - w    # initial remainder
dt = 5 * N**(-2)
A = np.eye(N - 1) - eps * dt * D2_int
for step in range(n_steps):
    vx = v[1:N] + w[1:N]
    rhs = v[1:N] + dt * (vx - vx**3)
    v[1:N] = np.linalg.solve(A, rhs)
    u = v + w
```

In MATLAB:

```matlab
[D, x] = cheb_matrix(N);
D2 = D^2;
D2_int = D2(2:N, 2:N);
w = x;
v = u0 - w;
dt = 5 * N^(-2);
A = eye(N - 1) - eps * dt * D2_int;
[LA, UA] = lu(A);
for step = 1:n_steps
    vx = v(2:N) + w(2:N);
    rhs = v(2:N) + dt * (vx - vx.^3);
    v(2:N) = UA \ (LA \ rhs);
    u = v + w;
end
```

In Julia:

```julia
D, x = cheb_matrix(N)
D2 = D^2
D2_int = D2[2:N, 2:N]
w = copy(x)
v = u0 .- w
dt = 5 * N^(-2)
A = I(N - 1) - eps * dt * D2_int
F = lu(A)
for step in 1:n_steps
    vx = v[2:N] .+ w[2:N]
    rhs = v[2:N] .+ dt .* (vx .- vx.^3)
    v[2:N] = F \ rhs
    u = v .+ w
end
```

#figure(
  image("../figures/ch13/python/allen_cahn_fixed.pdf", width: 95%),
  caption: [Allen--Cahn equation @eq-allen-cahn-bc with fixed Dirichlet conditions $u(plus.minus 1) = plus.minus 1$ and $epsilon = 0.01$. The space-time plot shows three initial transition layers that persist for a long time (metastability) before two interfaces suddenly annihilate near $t approx 45$, leaving only the single central interface required by the boundary data.],
) <fig-allen-cahn-fixed>

#figure(
  image("../figures/ch13/python/allen_cahn_driven.pdf", width: 95%),
  caption: [Allen--Cahn equation with time-dependent right boundary $u(1, t) = 1 + sin^2(t\/5)$. The oscillating boundary pumps energy into the domain, shifting the equilibrium interface to $x approx -0.4$ and causing earlier annihilation near $t approx 30$.],
) <fig-allen-cahn-driven>

The code generating @fig-allen-cahn-fixed and @fig-allen-cahn-driven is available in:
- `codes/python/ch13/bc_allen_cahn.py`
- `codes/matlab/ch13/bc_allen_cahn.m`
- `codes/julia/ch13/bc_allen_cahn.jl`

=== Discussion

The space-time plots in @fig-allen-cahn-fixed and @fig-allen-cahn-driven provide a vivid computational demonstration of _metastability_, a phenomenon whose theoretical analysis involves delicate exponential asymptotics that are far from intuitive. The fixed-boundary simulation (@fig-allen-cahn-fixed) shows three transition layers that persist essentially unchanged for dozens of time units before abruptly annihilating --- a behaviour that no linear stability analysis would predict. The timescale of this persistence grows exponentially as $epsilon arrow 0$, making direct numerical simulation the most accessible route to studying the phenomenon quantitatively.

Comparing the two boundary scenarios highlights how boundary conditions can control interface dynamics. The oscillating right boundary in @fig-allen-cahn-driven breaks the symmetry of the equilibrium configuration, shifts the surviving interface to the left, and accelerates the annihilation event. From a methodological standpoint, the two simulations also contrast the two boundary imposition strategies: Method I (lifting) for the fixed case and Method II (direct overwrite) for the time-dependent case. The latter is simpler to implement when boundaries vary in time, since it avoids recomputing source terms associated with the time derivative of the lifting function.

== Nonlinear Boundary Conditions <sec-nonlinear-bc>

=== Radiation Boundary Conditions

In heat transfer, the Stefan--Boltzmann law dictates that a body radiates energy proportional to the _fourth power_ of its temperature. At a boundary, this leads to the nonlinear Robin condition
$ u_x(1, t) + sigma u(1, t)^4 = 0, $ <eq-radiation-bc>
where $sigma > 0$ is a radiation parameter. This condition represents the full Stefan--Boltzmann radiation law, a formulation critical for accurately modelling high-temperature transient heat conduction, magnetohydrodynamic (MHD) flows, and highly coupled thermodynamic boundary layers. Because of the $u^4$ nonlinearity, neither Method I (lifting) nor a simple tau row replacement suffices; Newton--Raphson iteration must be actively coupled to the spectral discretisation.

=== Newton Iteration at the Boundary

For an implicit time step, the discrete system couples the interior PDE with the nonlinear boundary condition. We solve this coupled system using Newton iteration.

At each time step, the unknowns are the nodal values $bold(u) = (u_0, u_1, dots, u_N)^top$. The residual vector $bold(F)(bold(u))$ has three components:
- *Row $0$ (Dirichlet at $x = 0$)*: $F_0 = u_0 - 1$.
- *Rows $1, dots, N - 1$ (heat equation)*: $F_i = u_i - u_i^("old") - Delta t [D^2 bold(u)]_i$.
- *Row $N$ (radiation at $x = 1$)*: $F_N = [D bold(u)]_N + sigma u_N^4$.

The Jacobian $J = partial bold(F) \/ partial bold(u)$ is assembled by differentiating each component. The Newton update is $bold(u)^((k+1)) = bold(u)^((k)) - J^(-1) bold(F)(bold(u)^((k)))$, converging quadratically from a good initial guess (the previous time step).

=== Computational Étude 13.4: Radiative Cooling at a Boundary <etude-radiative-bc>

We solve the heat equation on $[0, 1]$ with a hot wall at $x = 0$ and a radiating surface at $x = 1$:
$ u_t = u_(x x), quad u(0, t) = 1, quad u_x(1, t) + sigma u(1, t)^4 = 0. $ <eq-radiative-problem>
The initial condition is $u(x, 0) = 1$ (uniform hot state). As heat radiates from the right boundary, the temperature drops and the system approaches a steady state where the conductive flux $u_x(1)$ balances the radiative loss $sigma u(1)^4$.

We map $[0, 1]$ to $[-1, 1]$ using $x = (xi + 1) \/ 2$ (with scaled derivatives $D_x = 2 D_xi$) and discretise with $N = 32$ Chebyshev points. The backward Euler time step gives a nonlinear system at each step, solved by Newton iteration (typically 3--4 iterations to machine precision).

In Python:

```python
D, xi = cheb_matrix(N)
D = 2 * D  # scale for [0, 1]
D2 = D @ D
dt = 0.01
for step in range(n_steps):
    u_old = u.copy()
    for newton in range(10):
        F = np.zeros(N + 1)
        F[0] = u[0] - 1                            # Dirichlet
        F[1:N] = u[1:N] - u_old[1:N] - dt * D2[1:N] @ u
        F[N] = D[N] @ u + sigma * u[N]**4           # Radiation
        J = np.zeros((N + 1, N + 1))
        J[0, 0] = 1
        J[1:N, :] = np.eye(N + 1)[1:N] - dt * D2[1:N]
        J[N, :] = D[N]; J[N, N] += 4 * sigma * u[N]**3
        u -= np.linalg.solve(J, F)
        if np.linalg.norm(F) < 1e-12:
            break
```

In MATLAB:

```matlab
[D, xi] = cheb_matrix(N);
D = 2 * D;
D2 = D^2;
dt = 0.01;
for step = 1:n_steps
    u_old = u;
    for newton = 1:10
        F = zeros(N+1, 1);
        F(1) = u(1) - 1;
        F(2:N) = u(2:N) - u_old(2:N) - dt * D2(2:N,:) * u;
        F(N+1) = D(N+1,:) * u + sigma * u(N+1)^4;
        J = zeros(N+1);
        J(1, 1) = 1;
        J(2:N, :) = eye(N+1, N+1) - dt * D2(2:N,:);
        J(2:N, :) = J(2:N, :);
        J(N+1, :) = D(N+1,:);
        J(N+1, N+1) = J(N+1, N+1) + 4*sigma*u(N+1)^3;
        u = u - J \ F;
        if norm(F) < 1e-12, break; end
    end
end
```

In Julia:

```julia
D, ξ = cheb_matrix(N)
D = 2D
D2 = D^2
dt = 0.01
for step in 1:n_steps
    u_old = copy(u)
    for newton in 1:10
        F = zeros(N + 1)
        F[1] = u[1] - 1
        F[2:N] = u[2:N] - u_old[2:N] - dt * D2[2:N, :] * u
        F[N+1] = D[N+1, :]' * u + σ * u[N+1]^4
        J = zeros(N + 1, N + 1)
        J[1, 1] = 1.0
        J[2:N, :] = I(N + 1)[2:N, :] - dt * D2[2:N, :]
        J[N+1, :] = D[N+1, :]
        J[N+1, N+1] += 4σ * u[N+1]^3
        u .-= J \ F
        norm(F) < 1e-12 && break
    end
end
```

#figure(
  image("../figures/ch13/python/radiative_cooling.pdf", width: 95%),
  caption: [Radiative cooling problem @eq-radiative-problem with $sigma = 1$ and $N = 32$. _Left_: temperature profiles at $t = 0, 0.1, 0.5, 1, 5$ showing the approach to steady state. The right boundary cools as heat radiates away. _Right_: boundary temperature $u(1, t)$ versus time for $sigma = 0.1$, $1$, and $10$, demonstrating that stronger radiation leads to faster and deeper cooling.],
) <fig-radiative-cooling>

The code generating @fig-radiative-cooling is available in:
- `codes/python/ch13/bc_radiative.py`
- `codes/matlab/ch13/bc_radiative.m`
- `codes/julia/ch13/bc_radiative.jl`

=== Discussion

The temperature profiles in @fig-radiative-cooling illustrate a competition between conduction and radiation that pure analysis renders as an implicit algebraic equation for the steady-state boundary temperature but cannot easily display as a dynamical process. The computation reveals the full transient relaxation: the boundary temperature drops rapidly at first, then approaches the steady state on a timescale that depends strongly on the radiation parameter $sigma$. For large $sigma$, the Stefan--Boltzmann $u^4$ law creates an effective "stiffness" at the boundary --- small changes in temperature produce large changes in radiative flux --- which the Newton iteration handles naturally with quadratic convergence.

This etude also demonstrates a key advantage of Method II for nonlinear boundary conditions. The Jacobian of the coupled interior-boundary system includes the derivative of the $u^4$ radiation term, which enters only in a single matrix entry ($J_(N, N)$). The Newton iteration typically converges in three to four steps per time step, confirming that the previous time step provides an excellent initial guess. This coupling between a nonlinear algebraic boundary constraint and a linear interior PDE is a recurring pattern in engineering applications, and the tau approach combined with Newton iteration provides a systematic framework for such problems.

== Two-Dimensional Problems with Mixed Boundary Data <sec-2d-mixed-bc>

=== Kronecker Products with Boundary Replacement

In two dimensions, the spectral Laplacian on $[-1, 1]^2$ is assembled via Kronecker products as described in @ch-spectral-pde:
$ L = I_y times.o D_x^2 + D_y^2 times.o I_x, $ <eq-2d-laplacian>
where $times.o$ denotes the Kronecker product. The solution vector is formed by stacking the columns (or rows) of the two-dimensional grid values into a single vector of length $(N + 1)^2$. While purely Dirichlet conditions can be solved with highly efficient $cal(O)(N^3)$ complexity via the matrix diagonalisation methods pioneered by Haidvogel and Zang @HaidvogelZang1979, the full Kronecker product assembly paired with the Tau method remains necessary to resolve mixed or piecewise boundary conditions without resorting to complex iterative solvers.

Boundary conditions are imposed by Method II: for each grid point on the boundary of the square, the corresponding row of $L$ is replaced by the appropriate condition. For Dirichlet data, the row becomes a unit vector selecting that grid point, with the right-hand side set to the prescribed boundary value. This requires careful bookkeeping to identify which entries of the vectorised solution correspond to boundary nodes, but the principle is identical to the one-dimensional case.

=== Computational Étude 13.5: Laplace Equation with Piecewise Boundary Data <etude-laplace-2d>

We solve the Laplace equation on the unit square following the setup from Trefethen's Program 36 @Trefethen2000:
$ u_(x x) + u_(y y) = 0, quad -1 < x, y < 1, $ <eq-laplace-2d>
with piecewise boundary data:
$ u(x, y) = cases(
  sin^4(pi x) & "if" y = 1 "and" -1 < x < 0,
  frac(1, 5) sin(3 pi y) & "if" x = 1,
  0 & "otherwise."
) $ <eq-laplace-2d-bc>

The boundary data is discontinuous at two corner points, $(0, 1)$ and $(1, 1)$, which limits the convergence rate near these corners. Away from the discontinuities, the spectral method resolves the solution with high accuracy.

The implementation assembles the Kronecker product Laplacian @eq-2d-laplacian, identifies boundary nodes by their position in the vectorised grid, replaces the corresponding rows with identity rows, and solves the resulting linear system.

In Python:

```python
D, x = cheb_matrix(N)
D2 = D @ D
I = np.eye(N + 1)
L = np.kron(I, D2) + np.kron(D2, I)
y = x.copy()
XX, YY = np.meshgrid(x, y)
xx, yy = XX.flatten(), YY.flatten()
# Right-hand side (Laplace: all zeros for interior)
rhs = np.zeros((N + 1)**2)
# Identify and impose boundary conditions
for k in range((N + 1)**2):
    if xx[k] == -1 or xx[k] == 1 or yy[k] == -1 or yy[k] == 1:
        L[k, :] = 0; L[k, k] = 1
        if yy[k] == 1 and xx[k] < 0:
            rhs[k] = np.sin(np.pi * xx[k])**4
        elif xx[k] == 1:
            rhs[k] = 0.2 * np.sin(3 * np.pi * yy[k])
u = np.linalg.solve(L, rhs).reshape(N + 1, N + 1)
```

In MATLAB:

```matlab
[D, x] = cheb_matrix(N);
D2 = D^2;
I = eye(N + 1);
L = kron(I, D2) + kron(D2, I);
[XX, YY] = meshgrid(x, x);
xx = XX(:); yy = YY(:);
rhs = zeros((N+1)^2, 1);
for k = 1:(N+1)^2
    if abs(xx(k)) == 1 || abs(yy(k)) == 1
        L(k, :) = 0; L(k, k) = 1;
        if yy(k) == 1 && xx(k) < 0
            rhs(k) = sin(pi * xx(k))^4;
        elseif xx(k) == 1
            rhs(k) = 0.2 * sin(3*pi*yy(k));
        end
    end
end
u = reshape(L \ rhs, N+1, N+1);
```

In Julia:

```julia
D, x = cheb_matrix(N)
D2 = D^2
Id = I(N + 1)
L = kron(Id, D2) + kron(D2, Id)
XX = [xi for yi in x, xi in x]
YY = [yi for yi in x, xi in x]
xx, yy = vec(XX), vec(YY)
rhs = zeros((N + 1)^2)
L = Matrix(L)
for k in 1:(N+1)^2
    if abs(xx[k]) ≈ 1 || abs(yy[k]) ≈ 1
        L[k, :] .= 0; L[k, k] = 1
        if yy[k] ≈ 1 && xx[k] < 0
            rhs[k] = sin(π * xx[k])^4
        elseif xx[k] ≈ 1
            rhs[k] = 0.2 * sin(3π * yy[k])
        end
    end
end
u = reshape(L \ rhs, N + 1, N + 1)
```

#figure(
  image("../figures/ch13/python/laplace_2d_mixed.pdf", width: 95%),
  caption: [Laplace equation @eq-laplace-2d on $[-1, 1]^2$ with piecewise Dirichlet boundary data @eq-laplace-2d-bc, solved with $N = 24$. The surface plot shows the smooth interior solution interpolating between the $sin^4(pi x)$ hump on the top-left boundary and the $sin(3 pi y) \/ 5$ oscillation on the right boundary. The solution satisfies $Delta u = 0$ to spectral accuracy in the interior.],
) <fig-laplace-2d>

The code generating @fig-laplace-2d is available in:
- `codes/python/ch13/bc_laplace_2d.py`
- `codes/matlab/ch13/bc_laplace_2d.m`
- `codes/julia/ch13/bc_laplace_2d.jl`

=== Discussion

The surface plot in @fig-laplace-2d demonstrates that the Chebyshev spectral method resolves the two-dimensional Laplacian with high accuracy in the interior, even when the boundary data is only piecewise smooth. The discontinuities at the corners $(0, 1)$ and $(1, 1)$ limit the global convergence rate to algebraic, as the Gibbs phenomenon propagates along the boundary edges. Away from these corners, however, the solution is analytic and the spectral method achieves its full exponential convergence rate.

From an implementation perspective, the main challenge in two dimensions is the bookkeeping required to identify boundary nodes in the vectorised (Kronecker product) system. Each row of the $(N + 1)^2 times (N + 1)^2$ operator corresponding to a boundary grid point must be replaced by an identity row, with the right-hand side set to the prescribed boundary value. The logic is straightforward but error-prone when the grid ordering (column-major versus row-major) is not handled consistently. This etude provides a template that extends directly to more complex two-dimensional problems with mixed Dirichlet, Neumann, and Robin conditions on different segments of the boundary.

== Frequency-Dependent Boundary Conditions and Quadratic Eigenvalue Problems <sec-qep>

The études so far have treated boundary conditions where the prescribed data --- values, derivatives, or nonlinear expressions --- are independent of the solution's eigenvalue or frequency. We now turn to a fundamentally different situation: boundary conditions that _depend on the unknown eigenvalue_. This makes Method I (lifting) impossible, since the lifting function would need to know the eigenvalue in advance. Method II (tau row replacement) handles this naturally, at the cost of producing a _polynomial eigenvalue problem_ rather than a standard one. Specifically, this results in a Quadratic Eigenvalue Problem (QEP), which requires specialised linearisation techniques --- such as mapping the system to companion matrices --- to be solved reliably, as extensively surveyed by Tisseur and Meerbergen @Tisseur2001.

=== Physical Motivation: Quasinormal Modes of Black Holes

When a black hole is perturbed --- for example, by the inspiral and merger of a binary system --- it "rings" at characteristic complex frequencies called _quasinormal modes_ (QNMs). The real part of each frequency determines the oscillation rate, while the imaginary part determines the damping rate. These QNMs are the gravitational-wave signature of the newly formed black hole settling into its final stationary state. For classical reviews, see Nollert @Nollert1999, Kokkotas and Schmidt @KokkotasSchmidt1999, Berti, Cardoso, and Starinets @Berti2009, and Konoplya and Zhidenko @Konoplya2011. More recently, the direct analytical mapping of quantum mechanical spectral techniques to QNM boundary value problems has been comprehensively detailed by Hatsuda and Kimura @HatsudaKimura2021.

The computation of QNMs reduces to a one-dimensional eigenvalue problem in the radial direction. In the so-called _tortoise coordinate_ $x in (-infinity, +infinity)$, the linearised perturbation equation for a massless scalar field on a black hole background takes the Schrödinger-like form
$ frac(d^2 psi, d x^2) + (omega^2 - V(x)) psi = 0, $ <eq-qnm-ode>
where $omega$ is the complex quasinormal frequency and $V(x)$ is an effective potential that depends on the black hole geometry. The key difficulty lies in the boundary conditions.

Recent advancements in computational relativity rely heavily on Chebyshev spectral collocation to resolve the quasinormal modes of exotic and alternative spacetime geometries. By casting the metric perturbation equations into a matrix-based QEP, spectral methods have recently been utilised to uncover the non-perturbative spectra of rapidly rotating black holes in shift-symmetric gravity @CanoRuiperez2024. Furthermore, recent spectral analyses have demonstrated the physical stability of noncommutative geometry-inspired dirty black holes and Morris--Thorne wormholes, successfully identifying previously undetected overdamped modes indicative of rapid decay without oscillation in near-extremal configurations @Batic2024b @Batic2025a. These studies confirm that, for complex, frequency-dependent boundaries where WKB approximations prove inadequate, the spectral Tau method remains the definitive numerical tool.

=== The Pöschl--Teller Potential

For pedagogical clarity, we work with the Pöschl--Teller potential @PoschlTeller1933
$ V(x) = V_0 op("sech")^2(x), $ <eq-poschl-teller>
which approximates the Regge--Wheeler potential of a Schwarzschild black hole near the peak of the potential barrier. The parameter $V_0$ controls the height of the barrier; for the Schwarzschild case with angular momentum number $ell$, the correspondence is $V_0 = ell(ell + 1) \/ 4$. This model potential has the crucial advantage that exact QNM frequencies are known analytically @FerrariMashhoon1984, providing a benchmark for our numerical method. The exact fundamental quasinormal frequency for $V_0 = 2$ (corresponding to $ell = 1$ in the Schwarzschild case) is
$ omega_0 = sqrt(V_0 - 1\/4) - i\/2 = sqrt(7) \/ 2 - i \/ 2. $ <eq-qnm-exact>
More generally, the overtone frequencies are $omega_n = sqrt(V_0 - (n + 1\/2)^2) - i(n + 1\/2)$ for $n = 0, 1, 2, dots$

The Pöschl--Teller potential decays exponentially as $x arrow plus.minus infinity$, mimicking the essential physics: the potential barrier near the black hole's photon sphere, with flat asymptotics representing the horizon (left) and spatial infinity (right). In the context of the Schwarzschild black hole, the QNM spectrum and its physical interpretation have been extensively studied using spectral methods by Batic and Dutykh @Batic2024 @Batic2024b @Batic2024c @Batic2025 @Batic2025a @Batic2025b.

=== Why Method I Fails: The Eigenvalue Appears in the BCs

The physical boundary conditions for QNMs encode a precise asymptotic behaviour:
- *At the horizon* ($x arrow -infinity$): the wave must be purely _ingoing_, $psi tilde e^(-i omega x)$.
- *At spatial infinity* ($x arrow +infinity$): the wave must be purely _outgoing_, $psi tilde e^(i omega x)$.

Differentiating these asymptotic forms gives the boundary conditions
$ psi'(-L) + i omega psi(-L) = 0, quad psi'(L) - i omega psi(L) = 0, $ <eq-qnm-bc>
where $L$ is a large truncation length chosen so that $V(plus.minus L) approx 0$.

The crucial observation is that _the unknown eigenvalue $omega$ appears in the boundary conditions_. This makes Method I impossible: we cannot construct a lifting function or restrict the basis without knowing $omega$ in advance. Method II, on the other hand, naturally accommodates the $omega$-dependent boundary rows, producing a _polynomial_ eigenvalue problem.

=== Quadratic Eigenvalue Problem Formulation

To cast the problem in standard form, we introduce the substitution $lambda = -i omega$, so that $omega = i lambda$ and $omega^2 = -lambda^2$. The ODE @eq-qnm-ode becomes
$ psi'' - V(x) psi - lambda^2 psi = 0, $ <eq-qnm-lambda>
and the boundary conditions @eq-qnm-bc become
$ psi'(-L) - lambda psi(-L) = 0, quad psi'(L) + lambda psi(L) = 0. $ <eq-qnm-bc-lambda>

Discretising on $N + 1$ Chebyshev points mapped to $[-L, L]$, the domain scaling gives $D = D_N \/ L$ (the differentiation matrix scaled for the interval $[-L, L]$). The discrete system takes the form of a _quadratic eigenvalue problem_ (QEP) @Tisseur2001:
$ (M_0 + lambda M_1 + lambda^2 M_2) bold(psi) = bold(0), $ <eq-qep>
where the coefficient matrices are assembled row by row:

*Interior rows* ($i = 1, dots, N - 1$):
$ [M_0]_i = [D^2 - op("diag")(V(x_j))]_i, quad [M_1]_i = 0, quad [M_2]_i = [-I]_i. $

*Boundary row $0$* ($x_0 = L$, the outgoing condition):
$ [M_0]_(0, j) = D_(0, j), quad [M_1]_(0, j) = delta_(0 j), quad [M_2]_(0, j) = 0. $

*Boundary row $N$* ($x_N = -L$, the ingoing condition):
$ [M_0]_(N, j) = D_(N, j), quad [M_1]_(N, j) = -delta_(N j), quad [M_2]_(N, j) = 0. $

=== Linearisation and Solution

The QEP @eq-qep can be solved by _companion linearisation_: introducing the auxiliary variable $bold(phi) = lambda bold(psi)$, we obtain the $2(N + 1) times 2(N + 1)$ generalised eigenvalue problem
$ underbrace(mat(M_0, M_1; bold(0), I), A) vec(bold(psi), bold(phi)) = lambda underbrace(mat(bold(0), -M_2; I, bold(0)), B) vec(bold(psi), bold(phi)). $ <eq-linearised-qep>

This is solved by a standard generalised eigenvalue solver. The computed eigenvalues $lambda$ are converted back to quasinormal frequencies via $omega = i lambda$. Physical QNMs are identified by the criterion $op("Im")(omega) < 0$, corresponding to modes that decay in time (as required by the outgoing/ingoing radiation conditions).

The spectrum also contains _spurious_ eigenvalues arising from the domain truncation and the linearisation. These are readily identified: they do not converge as $N$ or $L$ is increased, while the physical QNMs converge exponentially.

=== Computational Étude 13.6: Quasinormal Modes of a Black Hole <etude-qnm>

We compute the QNM spectrum for the Pöschl--Teller potential @eq-poschl-teller with $V_0 = 2$, using $N = 80$ Chebyshev points on $[-L, L]$ with $L = 8$.

In Python:

```python
D, x_cheb = cheb_matrix(N)
D = D / L                    # scale for [-L, L]
x = L * x_cheb               # physical coordinates
D2 = D @ D
V = V0 / np.cosh(x)**2
I = np.eye(N + 1)
Z = np.zeros((N + 1, N + 1))
# Assemble QEP matrices
M0 = D2 - np.diag(V)
M1 = np.zeros_like(M0)
M2 = -I.copy()
# Boundary row 0: psi' + lambda*psi = 0
M0[0, :] = D[0, :]; M1[0, :] = 0; M1[0, 0] = 1; M2[0, :] = 0
# Boundary row N: psi' - lambda*psi = 0
M0[N, :] = D[N, :]; M1[N, :] = 0; M1[N, N] = -1; M2[N, :] = 0
# Companion linearisation
A = np.block([[M0, M1], [Z, I]])
B = np.block([[Z, -M2], [I, Z]])
eigenvalues = scipy.linalg.eig(A, B, right=False)
omega = 1j * eigenvalues
```

In MATLAB:

```matlab
[D, x_cheb] = cheb_matrix(N);
D = D / L;
x = L * x_cheb;
D2 = D^2;
V = V0 ./ cosh(x).^2;
I = eye(N + 1);
% Assemble QEP coefficient matrices
M0 = D2 - diag(V);
M1 = zeros(N+1);
M2 = -I;
M0(1,:) = D(1,:); M1(1,1) = 1; M2(1,:) = 0;
M0(N+1,:) = D(N+1,:); M1(N+1,N+1) = -1; M2(N+1,:) = 0;
% MATLAB has polyeig for polynomial eigenvalue problems
lambda = polyeig(M0, M1, M2);
omega = 1i * lambda;
```

In Julia:

```julia
D, x_cheb = cheb_matrix(N)
D = D / L
x = L * x_cheb
D2 = D^2
V = V0 ./ cosh.(x).^2
Id = Matrix(1.0I, N + 1, N + 1)
Z = zeros(N + 1, N + 1)
M0 = D2 - diagm(V)
M1 = zeros(N + 1, N + 1)
M2 = -copy(Id)
M0[1, :] = D[1, :]; M1[1, 1] = 1.0; M2[1, :] .= 0
M0[N+1, :] = D[N+1, :]; M1[N+1, N+1] = -1.0; M2[N+1, :] .= 0
A = [M0 M1; Z Id]
B = [Z -M2; Id Z]
λs = eigvals(A, B)
ω = 1im .* λs
```

The fundamental QNM frequency computed with $N = 80$ and $L = 8$ is
$ omega_0^("num") approx 1.3227 - 0.5002 i, $
with an error $|omega_0^("num") - omega_0^("exact")| approx 2.3 times 10^(-4)$. The accuracy is limited by the domain truncation at $L = 8$; increasing $L$ to $9$ reduces the error to $approx 7.5 times 10^(-5)$. To reach machine precision, one must increase both $N$ and $L$ simultaneously, balancing spectral resolution against the domain extent needed to accommodate the decaying eigenfunction tails.

#figure(
  image("../figures/ch13/python/qnm_spectrum.pdf", width: 95%),
  caption: [Quasinormal mode spectrum for the Pöschl--Teller potential @eq-poschl-teller with $V_0 = 2$, $N = 80$, $L = 8$. Physical QNMs (filled circles) cluster in the lower half of the complex $omega$-plane ($op("Im")(omega) < 0$), while spurious modes (grey dots) appear scattered. The fundamental mode $omega_0 approx 1.323 - 0.500 i$ agrees with the exact Ferrari--Mashhoon result @FerrariMashhoon1984 to approximately four significant digits (limited by domain truncation).],
) <fig-qnm-spectrum>

#figure(
  image("../figures/ch13/python/qnm_eigenfunctions.pdf", width: 95%),
  caption: [Eigenfunctions of the first two quasinormal modes. _Left_: the fundamental mode ($n = 0$), showing a single antinode centred on the potential peak. _Right_: the first overtone ($n = 1$), with one nodal crossing. The oscillatory tails at large $|x|$ reflect the outgoing/ingoing wave character imposed by the boundary conditions @eq-qnm-bc.],
) <fig-qnm-eigenfunctions>

#figure(
  image("../figures/ch13/python/qnm_convergence.pdf", width: 95%),
  caption: [Convergence of the fundamental QNM frequency $omega_0$. _Left_: error $|omega_0^("num") - omega_0^("exact")|$ versus $N$ for fixed $L = 8$, showing initial spectral convergence that saturates at $approx 2 times 10^(-4)$ due to domain truncation. _Right_: error versus truncation length $L$ for fixed $N = 80$, showing exponential decay until $L approx 9$, beyond which the grid becomes too coarse to resolve the solution and the error increases again.],
) <fig-qnm-convergence>

The code generating @fig-qnm-spectrum, @fig-qnm-eigenfunctions, and @fig-qnm-convergence is available in:
- `codes/python/ch13/bc_qnm_poschl_teller.py`
- `codes/matlab/ch13/bc_qnm_poschl_teller.m`
- `codes/julia/ch13/bc_qnm_poschl_teller.jl`

=== Discussion

The QNM spectrum in @fig-qnm-spectrum reveals a feature that has no counterpart in classical Sturm--Liouville theory: the eigenvalues are _complex_, with negative imaginary parts encoding the damping rate of each mode. The physical QNMs form an ordered sequence in the lower half of the complex $omega$-plane, while spurious eigenvalues --- artefacts of the finite domain truncation and the companion linearisation --- appear scattered and do not converge with increasing $N$ or $L$. Distinguishing physical from spurious modes by checking convergence is a practical necessity that the computation makes straightforward.

The eigenfunctions shown in @fig-qnm-eigenfunctions exhibit oscillatory tails at large $|x|$ that reflect the outgoing and ingoing wave character imposed by the boundary conditions @eq-qnm-bc. These tails grow without bound as $|x| arrow infinity$ (since the modes are not $L^2$-normalizable), which is precisely why the truncation length $L$ must be chosen with care. @fig-qnm-convergence quantifies this trade-off: for fixed $L$, the error in $omega_0$ saturates at a level determined by the domain truncation, while for fixed $N$, there is an optimal $L$ beyond which the Chebyshev grid becomes too coarse to resolve the increasingly oscillatory tails. This interplay between spectral resolution and domain extent is a recurring theme in unbounded-domain computations, and the QNM problem provides an instructive example where both parameters must be balanced simultaneously.

== Eigenvalue Problems with Mixed Boundary Conditions <sec-eigenvalue-mixed-bc>

=== Sturm--Liouville with Neumann/Robin Endpoints

The QNM problem of the preceding section is exotic in that the boundary conditions depend on the eigenvalue itself. But even the classical Sturm--Liouville eigenproblem
$ -frac(d, d x) (p(x) frac(d u, d x)) + q(x) u = lambda w(x) u, quad a < x < b, $ <eq-sturm-liouville>
becomes nontrivial when the boundary conditions involve derivatives. With Dirichlet conditions at both endpoints, the eigenvalues are computed by stripping the differentiation matrix (Method I), as demonstrated in @ch-bvp. But a Neumann condition $u'(b) = 0$ or a Robin condition $alpha u(b) + beta u'(b) = 0$ requires Method II: the boundary row of the operator matrix must be replaced, and the eigenvalue problem solved on the full $(N + 1) times (N + 1)$ system.

=== Computational Étude 13.7: Vibrating String with a Free End <etude-vibrating-string>

We compute the eigenmodes of a vibrating string that is fixed at the left end and free at the right:
$ -u_(x x) = lambda u, quad 0 < x < 1, quad u(0) = 0, quad u'(1) = 0. $ <eq-vibrating-string>
The exact eigenvalues and eigenfunctions are
$ lambda_n = ((2 n - 1) pi \/ 2)^2, quad u_n(x) = sin((2 n - 1) pi x \/ 2), quad n = 1, 2, 3, dots $ <eq-vibrating-exact>
The fundamental mode is a quarter-wavelength standing wave, and successive overtones add half-wavelengths.

We map $[0, 1]$ to $[-1, 1]$ and construct the scaled second derivative matrix. The tau modifications are:
- Row $0$ ($x = 1$): replace with $D[0, :]$ (Neumann condition $u'(1) = 0$).
- Row $N$ ($x = 0$): replace with $bold(e)_N^top$ (Dirichlet condition $u(0) = 0$).

The eigenvalue problem becomes $A bold(u) = lambda bold(u)$, where $A$ is the modified $-D^2$ matrix.

In Python:

```python
D, xi = cheb_matrix(N)
D = 2 * D  # scale for [0, 1]
D2 = D @ D
A = -D2.copy()
# Neumann at x = 1 (row 0): u'(1) = 0
A[0, :] = D[0, :]
# Dirichlet at x = 0 (row N): u(0) = 0
A[N, :] = 0; A[N, N] = 1
eigenvalues, eigenvectors = np.linalg.eig(A)
idx = np.argsort(np.real(eigenvalues))
lam = np.real(eigenvalues[idx])
```

In MATLAB:

```matlab
[D, xi] = cheb_matrix(N);
D = 2 * D;
D2 = D^2;
A = -D2;
A(1, :) = D(1, :);           % Neumann at x = 1
A(N+1, :) = 0; A(N+1, N+1) = 1;  % Dirichlet at x = 0
[V, Lam] = eig(A);
[lam, idx] = sort(diag(Lam));
V = V(:, idx);
```

In Julia:

```julia
D, ξ = cheb_matrix(N)
D = 2D
D2 = D^2
A = -D2
A[1, :] = D[1, :]
A[N+1, :] .= 0; A[N+1, N+1] = 1.0
F = eigen(A)
idx = sortperm(real.(F.values))
λ = real.(F.values[idx])
```

#figure(
  image("../figures/ch13/python/vibrating_string.pdf", width: 95%),
  caption: [Eigenmodes of the vibrating string @eq-vibrating-string with $N = 36$. _Left_: the first six eigenmodes (solid lines) agree with the exact solutions $sin((2n-1) pi x \/ 2)$ (dashed). _Right_: eigenvalue errors $|lambda_n^("num") - lambda_n^("exact")|$ versus mode number, showing that the first $approx 2 N \/ 3$ eigenvalues are resolved to high accuracy before the resolution limit is reached.],
) <fig-vibrating-string>

The code generating @fig-vibrating-string is available in:
- `codes/python/ch13/bc_vibrating_string.py`
- `codes/matlab/ch13/bc_vibrating_string.m`
- `codes/julia/ch13/bc_vibrating_string.jl`

=== Discussion

The eigenvalue comparison in @fig-vibrating-string confirms that the tau approach resolves the first $approx 2 N \/ 3$ eigenvalues to high accuracy, a rule of thumb consistent with the resolution limits observed for spectral eigenvalue problems in @ch-bvp. The remaining eigenvalues suffer from under-resolution: the corresponding eigenfunctions oscillate too rapidly for the $N$-point Chebyshev grid to represent faithfully. This resolution boundary is sharp and predictable, making it straightforward to determine how large $N$ must be to capture a desired number of modes.

The mixed Dirichlet--Neumann boundary conditions in this problem produce the quarter-wavelength mode family @eq-vibrating-exact, with eigenvalues $(2 n - 1)^2 pi^2 \/ 4$ that are distinct from the full-wavelength Dirichlet--Dirichlet spectrum $n^2 pi^2$ studied in earlier chapters. The ability to switch between boundary condition types by replacing a single row of the operator matrix --- without modifying the interior discretisation --- illustrates the modularity and flexibility of the tau approach. This same row-replacement strategy extends to more general Sturm--Liouville problems with Robin endpoints, periodic conditions, or even boundary conditions coupling the values at the two endpoints.

== A Non-Exhaustive Literature Overview

The treatment of boundary conditions within spectral collocation methods represents a critical mathematical juncture where theoretical elegance must inevitably reconcile with the complex demands of physical reality. While the foundational principles of spectral approximation thrive in domains governed by simple, homogeneous Dirichlet constraints, the transition to advanced computational modelling introduces a profound layer of complexity. Modern physical systems demand the resolution of inhomogeneous data, mixed derivative constraints, highly nonlinear energy fluxes, dynamic moving interfaces, and, most intricately, boundary operators that depend inherently upon the unknown eigenvalue of the system.

=== The Evolution of Spectral Boundary Imposition: Tau Methods and Lifting

The methodological imposition of boundary conditions in spectral methods has evolved through decades of rigorous numerical analysis, branching into two dominant philosophies: the restriction of the interpolant (Method I) and the replacement of discrete equations (Method II).

The genesis of the latter approach lies in the seminal work of Cornelius Lanczos @Lanczos1938, who introduced the classical "Tau method" within the context of Chebyshev polynomial expansions. Lanczos recognised a fundamental limitation in purely periodic or homogeneous orthogonal bases: they could not easily accommodate arbitrary boundary constraints without destroying the uniform convergence properties of the expansion. By modifying the highest-order expansion coefficients to explicitly satisfy the boundary constraints, the Tau method circumvented these restrictive assumptions, effectively introducing a spectral alternative to the local Taylor-series approximations inherent in classical finite difference schemes.

This modal formulation of the Tau method remained the standard until the 1970s when the rediscovery of the Fast Fourier Transform (FFT) catalysed a shift toward physical-space techniques. The definitive monograph by Gottlieb and Orszag @GottliebOrszag1977 formally established the "pseudospectral" or collocation method, demonstrating that boundary conditions could be imposed directly on the spatial grid by replacing the boundary rows of the discrete differential operator with equations encoding the exact boundary constraints. This technique, commonly referred to as "boundary bordering" or "tau row replacement," provides a uniform mathematical framework for handling Neumann, Robin, and mixed conditions. The encyclopaedic treatment by Boyd @Boyd2000 further codified these strategies, categorising the implementation of boundary conditions into basis recombination (a pure Galerkin approach), penalty methods, and collocation tau row replacement. The comprehensive monographs by Canuto _et al._ @Canuto2006 and Shen, Tang, and Wang @ShenTangWang2011 give rigorous error analysis for various boundary condition formulations in both Galerkin and collocation settings. Trefethen's _Spectral Methods in MATLAB_ @Trefethen2000, particularly Chapter 13, is the direct inspiration for the material in this chapter. His Programs 32--37 demonstrate lifting, tau row replacement, and mixed-dimension discretisations (Fourier in one direction, Chebyshev in the other) with remarkable economy of code.

Conversely, for systems governed by inhomogeneous Dirichlet data, the "lifting technique" (Method I) emerged as an elegant and computationally efficient alternative to tau row replacement. By constructing an explicit, smooth lifting function that perfectly satisfies the prescribed boundary data, the inhomogeneous boundary value problem is mathematically reduced to a homogeneous one, allowing the interior system to be solved on a restricted subspace of polynomials that vanish at the endpoints. Recent literature has expanded the utility of the lifting technique far beyond simple linear interpolants. In the context of fractional diffusion equations and non-local operators, lifting techniques have been utilised to bypass the breakdown of classical trace theorems, enabling the formulation of spectral fractional Laplacians with non-zero boundary constraints @Mustapha2025 @Ervin2020.

=== Nonlinear Boundary Constraints: The Stefan--Boltzmann Law and Iterative Coupling

When boundary conditions transcend linear combinations of values and derivatives, the spectral framework must be intimately coupled with iterative root-finding algorithms. A canonical and highly challenging example arises in thermodynamics and heat transfer, where the Stefan--Boltzmann law dictates that an absolutely black body radiates energy proportional to the fourth power of its absolute surface temperature. The resulting nonlinear Robin-type condition, typically expressed as $alpha u + nu nabla u + sigma u^4 = g$, introduces extreme stiffness into the boundary layer and cannot be accommodated by simple matrix surgery. Instead, the algebraic nonlinear boundary constraint must be solved simultaneously with the interior partial differential equation (PDE). In a discrete spectral framework, this is achieved by incorporating Newton--Raphson iterations directly into the implicit time-stepping loop, treating the boundary nodes as fully coupled nonlinear algebraic variables.

The mathematical justification for utilising Robin conditions in thermal flows has been subjected to rigorous analysis. Recent studies comparing ideal transmission conditions (enforcing absolute temperature and flux continuity across an interface) with simplified Robin approximations demonstrate that the latter can be justified mathematically for high convection velocities, provided the thermal conductivity of the surrounding fluid is significantly lower than that of the solid body @Clarte2020. However, in highly coupled thermodynamic systems --- such as modelling transient heat conduction in composite materials, or assessing the thermal shielding of hypersonic aerospace vehicles --- the full nonlinear Stefan--Boltzmann condition must be retained to accurately model bilateral heat exchange, as the surrounding medium is actively heated by the radiation.

=== Time-Dependent Boundaries, Phase Separation, and Allen--Cahn Metastability

In the realm of time-dependent boundary conditions and interior phase dynamics, the Allen--Cahn equation @AllenCahn1979 serves as a premier mathematical testbed. Originally formulated to model the microscopic dynamics of antiphase boundary motion and domain coarsening in binary alloys, the nonlinear reaction-diffusion equation is renowned for exhibiting "metastability." This phenomenon describes a state where sharp transition layers (interfaces) separating stable phases persist in an unstable equilibrium for exponentially long times before suddenly accelerating and annihilating.

The rigorous asymptotic analysis of this slow evolution was pioneered in the late 1980s and 1990s by Carr and Pego, as well as Fusco and Hale @CarrPego1989 @Alikakos1991. By utilising spectral estimates and dynamic manifold projections, they proved that the interface velocity scales at an exponentially small rate, approximately $cal(O)(e^(-C\/epsilon))$, where $epsilon$ dictates the interface width. Because the interface motion is so infinitesimally slow, standard low-order numerical methods often fail; artificial numerical dissipation can easily overwhelm the true physical dynamics, causing the interfaces to collapse prematurely. Spectral methods, virtually free of dispersive and dissipative errors, are therefore the optimal choice for resolving these metastable states @Trefethen2000. Spectral discretisations for Allen--Cahn and the related Cahn--Hilliard equation are analysed in detail by Shen and Yang @ShenYang2010, who develop energy-stable schemes that preserve the variational structure of the continuous problem.

=== Multidimensional Geometries and Kronecker Product Algebra

The extension of advanced boundary conditions to two- and three-dimensional domains fundamentally relies on the algebraic properties of Kronecker products. For elliptic equilibrium problems --- such as the Poisson, Laplace, and Helmholtz equations on rectangular domains --- the tensor-product grid yields a discrete high-order Laplacian formulated as the Kronecker sum of one-dimensional second-derivative operators:
$ L = I times.o D_x^2 + D_y^2 times.o I. $

The historical challenge of applying global spectral methods to multidimensional domains was the computational complexity of inverting the resulting dense matrices. For a 2D grid of size $N times N$, the full Kronecker matrix is of size $N^2 times N^2$, making direct Gaussian elimination an $cal(O)(N^6)$ operation. To circumvent this, Haidvogel and Zang @HaidvogelZang1979 introduced the matrix diagonalisation method. By exploiting the eigenvalue and eigenvector decomposition of the 1D differentiation matrices, the multidimensional system is decoupled into a sequence of independent one-dimensional problems, thereby reducing the computational complexity to a highly efficient $cal(O)(N^3)$. This approach was further refined by Shen @Shen1995, who constructed specific recombined Legendre and Chebyshev basis functions that satisfied homogeneous Dirichlet boundary conditions automatically, leading to sparse, banded matrices and optimal $cal(O)(N^2)$ complexity for rectangular and polar geometries.

However, a significant limitation of the matrix diagonalisation technique is its rigidity regarding boundary conditions. When non-Dirichlet boundary conditions --- such as Neumann, Robin, or piecewise mixed constraints --- are introduced, the matrix diagonalisation method breaks down, requiring messy row operations or costly iterative procedures that destroy the optimal complexity. In these scenarios, the full Kronecker product assembly, combined with the Tau method (row replacement), remains the most robust and generalised approach.

=== Frequency-Dependent Boundaries and the Quadratic Eigenvalue Problem

A profound shift in the mathematical structure of spectral methods occurs when the boundary condition itself depends directly upon the unknown eigenvalue (or frequency) of the system. This paradigm frequently arises in dissipative acoustics, structural vibrations, fluid-structure interactions, and relativistic astrophysics. Because the eigenvalue $lambda$ (or $omega$) is unknown _a priori_, the lifting technique is rendered mathematically impossible; the lifting function would require prior knowledge of the eigenvalue to satisfy the boundary constraint. Consequently, the Tau method (boundary row replacement) becomes the exclusive pathway forward.

Replacing the boundary rows of the discrete differential operator with frequency-dependent algebraic equations fundamentally transforms the standard linear eigenvalue problem ($A bold(v) = lambda bold(v)$) into a Polynomial Eigenvalue Problem (PEP), most commonly the Quadratic Eigenvalue Problem (QEP). The definitive mathematical survey on the properties, conditioning, and numerical solution techniques for QEPs was provided by Tisseur and Meerbergen @Tisseur2001. The standard methodology for solving the QEP in a spectral collocation framework involves "linearisation" --- casting the quadratic system of size $N times N$ into a generalised linear eigenvalue problem of size $2N times 2N$ using companion matrices.

While companion linearisation is the most widely adopted technique due to its direct implementation in standard software libraries (_e.g._, MATLAB's `polyeig`), it introduces significant numerical challenges. The conditioning of the linearised matrix polynomials can be highly sensitive to the scaling of the constituent matrices. Poor scaling can lead to severe backward errors and large variations in the pseudospectrum @Tisseur2001. For large-scale multidimensional problems, transforming the system into a dense $2N times 2N$ matrix is computationally prohibitive. Consequently, advanced iterative methods such as the Jacobi--Davidson method, inexact residual iteration methods, and Krylov subspace techniques have been developed to target specific regions of the quadratic spectrum without computing the full dense eigensystem.

=== Modern Frontiers: Quasinormal Modes in Black Hole Perturbation Theory

The most striking and mathematically rigorous contemporary application of spectral methods with frequency-dependent boundary conditions lies in gravitational wave astronomy and black hole perturbation theory. When a black hole is subjected to an external perturbation --- such as the asymmetric collapse of a star, or the inspiral and merger of a binary system --- it undergoes a period of exponentially damped oscillation known as the "ringdown" phase. The frequencies of this emitted gravitational radiation are strictly complex; the real part dictates the oscillation rate, while the imaginary part governs the exponential decay. These characteristic, intrinsic frequencies are known as Quasinormal Modes (QNMs).

The computation of QNMs is essentially an eigenvalue problem governed by master wave equations, such as the Regge--Wheeler or Zerilli equations for static Schwarzschild black holes, and the Teukolsky equation for rotating Kerr black holes. The effective potential barrier in these Schrödinger-like equations decays to zero at both spatial infinity and at the event horizon. Consequently, the physical boundary conditions demand purely outgoing waves at spatial infinity and purely ingoing waves at the horizon. Because the wave speed and phase are inextricably tied to the unknown complex frequency $omega$, the boundary conditions inherently depend on the eigenvalue, culminating directly in a Quadratic Eigenvalue Problem (QEP).

The historical approach to calculating QNM spectra relied heavily on the WKB (Wentzel--Kramers--Brillouin) approximation and Leaver's continued fraction method. However, the modern era of computational relativity has seen a massive paradigm shift toward Chebyshev spectral collocation methods. This is driven by their geometric convergence rates and their unique ability to compute entire sectors of the spectrum simultaneously --- capturing not only the fundamental modes but also highly damped overtones that govern the early ringdown phase. An unconventional but highly influential review by Hatsuda and Kimura @HatsudaKimura2021 explicitly demonstrated how varied techniques from quantum mechanics, including Borel summations and Padé approximants, could be mapped directly onto spectral problems in black hole perturbation theory.

Recent breakthroughs (2024--2026) have established the pseudo-spectral method as the premier numerical tool for investigating alternative theories of gravity. The spectral method has been successfully deployed to compute QNMs for rapidly rotating black holes in shift-symmetric Einstein-scalar-Gauss--Bonnet theory, capturing the full non-perturbative spectrum where previous analytical methods were strictly limited to quadratic order perturbations @CanoRuiperez2024. Furthermore, a significant body of state-of-the-art research by Batic, Dutykh, and collaborators has aggressively pushed the boundaries of spectral QNM computation into exotic and alternative spacetime geometries. In 2024, they deployed a unified spectral approach to compute the QNMs of Lee--Wick black holes, validating the spectral method's robustness against complex, higher-derivative gravity modifications @Batic2024c. Subsequent landmark research in 2025 applied Chebyshev spectral collocation to noncommutative geometry-inspired Schwarzschild black holes @Batic2024b. By transforming the perturbation equations into an eigenvalue problem over a compact domain using the Chebyshev roots grid, this work fundamentally challenged previous claims in the literature regarding the instability of these quantum-corrected black holes. The spectral analysis revealed remarkable stability under scalar, electromagnetic, and gravitational perturbations. Crucially, the spectral QEP formulation identified entirely new "overdamped" modes --- indicative of rapid decay without oscillation --- that become particularly prominent in near-extremal regimes @Batic2024b. Further extending this rigorous framework, Batic, Dutykh, and Beek (2025) applied the spectral method to noncommutative geometry-inspired Morris--Thorne (Bronnikov--Ellis) wormholes @Batic2025a. The spectral approach corrected critical inaccuracies previously reported by WKB approximations, resolving scalar, electromagnetic, and vector-type gravitational perturbations with high precision. The spectral analysis validated that, for large rescaled mass parameters, the QNMs of the noncommutative wormhole transition smoothly to those of the classical Schwarzschild wormhole, solidifying the spectral method as the gold standard for mapping the gravitational-wave signatures of the universe @Batic2025a. These advanced computational études underscore the indispensable nature of the Tau method, global polynomial approximation, and QEP linearisations in deciphering the fundamental physics of the cosmos.

Beyond the methods presented here, several advanced boundary condition techniques merit mention. _Penalty_ and _simultaneous approximation term_ (SAT) methods impose boundary conditions weakly, adding penalisation terms to the differential operator rather than replacing equations. These methods are particularly important for hyperbolic systems where strong imposition can lead to instabilities. _Integral_ (or _nonlocal_) boundary conditions, where the constraint involves an integral of the solution over the domain, arise in certain inverse problems and can be handled by appending a quadrature-weighted row to the tau system. _Fractional_ boundary conditions, involving fractional-order derivatives at the boundary, are an emerging topic in the modelling of anomalous diffusion and non-local transport.

== Summary <sec-advanced-bc-summary>

This chapter has developed two systematic strategies for imposing boundary conditions in Chebyshev spectral methods, progressing from elementary Dirichlet data to frequency-dependent eigenvalue conditions:

+ *Method I (lifting)* reduces inhomogeneous Dirichlet conditions to homogeneous ones by subtracting a simple function satisfying the boundary data. The technique is efficient and exact, but limited to problems where the boundary values are explicitly known.

+ *Method II (tau / row replacement)* keeps the full system and replaces boundary rows with equations encoding the boundary conditions. This handles Neumann, Robin, nonlinear, time-dependent, and eigenvalue-dependent conditions uniformly.

+ *Nonlinear boundary conditions* (such as the Stefan--Boltzmann radiation law) couple the boundary algebraic equation to the interior PDE, requiring Newton iteration within the time-stepping loop.

+ *Two-dimensional problems* with mixed boundary data are handled by Method II applied to the Kronecker product Laplacian, with careful bookkeeping to identify boundary nodes in the vectorised system.

+ *Frequency-dependent boundary conditions* make Method I impossible and Method II essential. The unknown eigenvalue in the boundary condition produces a _quadratic eigenvalue problem_ (QEP), solved by companion linearisation and a generalised eigenvalue solver.

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    table(
      columns: (0.6fr, 1.2fr, 1fr, 1fr),
      align: (center, left, left, left),
      inset: (x: 0.8em, y: 0.5em),
      stroke: none,
      table.hline(stroke: 0.75pt + rgb("#142D6E")),
      table.header(
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Étude*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Problem*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Key Technique*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Method*],
      ),
      table.hline(stroke: 0.5pt + luma(180)),
      [13.1], [$u'' = e^(4 x)$, inhomogeneous Dirichlet], [Lifting function], [Method I],
      [13.2], [$-u'' + k^2 u = f$, Robin BC], [Tau rows: $alpha bold(e)^top + beta D[N, :]$], [Method II],
      [13.3], [Allen--Cahn with driven BCs], [Lifting vs direct overwrite], [Methods I & II],
      [13.4], [Heat eq. with $u_x + sigma u^4 = 0$], [Newton at boundary], [Method II],
      [13.5], [$nabla^2 u = 0$, piecewise Dirichlet (2D)], [Kronecker + row replacement], [Method II],
      [13.6], [QNM: $psi'' + (omega^2 - V) psi = 0$], [QEP linearisation], [Method II],
      [13.7], [$-u'' = lambda u$, Dirichlet + Neumann], [Tau eigenvalue problem], [Method II],
    ),
  ),
  caption: [Summary of computational études in this chapter.],
) <tbl-ch13-summary>

== Exercises <sec-advanced-bc-exercises>

*Exercise 13.1* (_Pure Neumann and compatibility_). Consider the Poisson equation $u_(x x) = cos(pi x)$ on $[-1, 1]$ with Neumann conditions $u'(plus.minus 1) = 0$ at both endpoints. (a) Verify that the compatibility condition $integral_(-1)^1 f(x) d x = 0$ is satisfied. (b) The solution is unique only up to an additive constant. Construct the discrete system using tau rows for both Neumann conditions, and observe that the matrix is singular. (c) Add the zero-mean constraint $sum_j w_j u_j = 0$ (using Clenshaw--Curtis quadrature weights $w_j$) as an additional row, replacing one of the interior equations. Solve the resulting bordered system and verify the solution against the exact answer $u(x) = -cos(pi x) \/ pi^2 + C$.

*Exercise 13.2* (_Fourier--Chebyshev wave tank_). Implement the two-dimensional wave equation $u_(t t) = u_(x x) + u_(y y)$ on the rectangular domain $[-3, 3] times [-1, 1]$ with periodic boundary conditions in $x$ (Fourier discretisation) and Neumann conditions $u_y(x, plus.minus 1, t) = 0$ in $y$ (Chebyshev with tau rows). Use leapfrog time stepping. Start with a Gaussian pulse $u(x, y, 0) = e^(-8((x + 3\/2)^2 + y^2))$ and set $u_t(x, y, 0)$ for a rightward-travelling wave. Investigate the CFL stability limit $Delta t tilde 5 \/ (N_x + N_y^2)$ experimentally.

*Exercise 13.3* (_QNM domain truncation_). In Étude 13.6, the truncation length $L$ controls how well the computational domain approximates the infinite line. For $V_0 = 2$, the exact fundamental frequency is given by @eq-qnm-exact. Fix $N = 80$ and compute $|omega_0^("num") - omega_0^("exact")|$ for $L = 2, 3, 4, dots, 12$. Plot the error versus $L$ on a semilogarithmic scale and explain the observed exponential convergence in terms of the decay rate of the Pöschl--Teller potential and the QNM eigenfunction.

*Exercise 13.4* (_Allen--Cahn metastability lifetime_). Modify Étude 13.3 to measure the _metastability lifetime_ $T(epsilon)$, defined as the time at which the solution $u(x, t)$ first becomes monotone (single interface). Compute $T(epsilon)$ for $epsilon = 0.02, 0.015, 0.01, 0.007, 0.005$ and plot $log T$ versus $1 \/ epsilon$. Verify that $T tilde e^(C \/ epsilon)$ for some constant $C$ that depends on the initial interface configuration.
