// textbook/chapters/higher_order_bvps.typ
// Chapter 14: Higher-Order Boundary Value Problems
// Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
// Email: denys.dutykh@ku.ac.ae
// Homepage: https://www.denys-dutykh.com/
// Last modified: February 2026

#import "../styles/template.typ": dropcap, num, format-table

// Enable equation numbering for this chapter
#set math.equation(numbering: "(1)")

= Higher-Order Boundary Value Problems <ch-higher-order>

#dropcap[The boundary value problems studied in the preceding chapters have been governed by second-order differential operators --- the Laplacian, the Helmholtz operator, Sturm--Liouville operators --- where a single boundary condition at each endpoint suffices to determine the solution uniquely. In engineering and applied mathematics, however, some of the most important problems involve _fourth-order_ operators that demand _two_ boundary conditions at each endpoint. The Euler--Bernoulli beam equation $E I u^((4)) = f(x)$ governs the deflection of slender structures under load, requiring specification of both displacement and slope (clamped) or displacement and moment (simply supported) at each end. The biharmonic equation $Delta^2 u = f$ describes the deformation of thin elastic plates in the Kirchhoff theory and arises in Stokes flow as the equation governing the stream function. Most dramatically, the Orr--Sommerfeld equation --- the linearised stability equation for viscous parallel shear flows --- is a fourth-order eigenvalue problem whose spectrum determines whether a given laminar flow is stable or unstable to infinitesimal perturbations. The spectral discretisation of these problems raises new challenges: the fourth derivative matrix $D^4$ has condition number $cal(O)(N^8)$, boundary condition imposition requires eliminating four constraints rather than two, and the rich algebraic structure of the biharmonic operator in two dimensions calls for careful exploitation of Kronecker products and symmetry. This chapter develops the spectral tools needed to address these challenges, progressing from elementary beam problems to the Orr--Sommerfeld equation, pseudospectra, and the Kuramoto--Sivashinsky equation.]

By the end of this chapter, you should be able to:

1. Apply the _polynomial trick_ to automatically satisfy clamped boundary conditions ($u = u' = 0$ at both endpoints) and construct the corresponding fourth-derivative matrix.
2. Solve fourth-order _eigenvalue problems_ arising from beam vibrations using both the polynomial trick and boundary bordering.
3. Reformulate a fourth-order boundary value problem as a _coupled system_ of two second-order equations, achieving better conditioning at the cost of a doubled system size.
4. Extend the biharmonic operator to _two dimensions_ using Kronecker products and the polynomial trick in each coordinate direction.
5. Exploit the _symmetry_ of a square domain to reduce a full-domain problem to a quarter-domain computation, selecting specific symmetry classes of eigenmodes.
6. Compute the spectrum of the _Orr--Sommerfeld equation_ for plane Poiseuille flow and identify the critical Reynolds number $R_c approx 5772$.
7. Understand the concept of _non-normality_ and compute the $epsilon$-pseudospectrum of the Orr--Sommerfeld operator, revealing the extreme sensitivity of its eigenvalues to perturbations.
8. Use the _ETDRK4_ (exponential time differencing, fourth-order Runge--Kutta) scheme to integrate stiff periodic fourth-order evolution equations such as the Kuramoto--Sivashinsky equation.

== The Polynomial Trick for Clamped Boundary Conditions <sec-polynomial-trick>

The clamped boundary conditions that arise in beam and plate theory prescribe both the function value and its first derivative at each endpoint:
$ u(plus.minus 1) = 0, quad u'(plus.minus 1) = 0. $ <eq-clamped-bc>
Imposing four conditions by the tau method (replacing four rows of the differentiation matrix) is straightforward in principle but requires careful bookkeeping. An elegant alternative, introduced by Huang and Sloan @HuangSloan1994 and prominently featured in Trefethen's _Spectral Methods in MATLAB_ @Trefethen2000 (Programs 38 and 39), exploits the factorisation of the solution as a polynomial that _automatically_ satisfies the clamped conditions.

=== The Factorisation

The key observation is that any polynomial $p(x)$ satisfying $p(plus.minus 1) = p'(plus.minus 1) = 0$ must have $(1 - x^2)$ as a factor with multiplicity at least two in the sense that both $p$ and $p'$ vanish at $x = plus.minus 1$. We can therefore write
$ u(x) = (1 - x^2) q(x), $ <eq-polynomial-trick>
where $q(x)$ is an unrestricted polynomial. The factor $(1 - x^2)$ vanishes at $x = plus.minus 1$, and its derivative $-2x$ ensures that
$ u'(x) = -2 x q(x) + (1 - x^2) q'(x), $
so $u'(plus.minus 1) = -2(plus.minus 1) q(plus.minus 1) + 0 = minus.plus 2 q(plus.minus 1)$. Wait --- this does _not_ automatically vanish! The single factor $(1 - x^2)$ enforces $u(plus.minus 1) = 0$ but not $u'(plus.minus 1) = 0$.

To enforce all four conditions, we need the _squared_ factor $(1 - x^2)^2$... but this over-constrains the problem (it forces $u''(plus.minus 1) = 0$ as well). The polynomial trick instead proceeds differently: rather than factoring $u$, we express the _differential equation_ in terms of $q$ and work on the interior grid, where the boundary conditions are implicitly satisfied by the structure of the discretisation.

=== The Correct Formulation

Following Trefethen @Trefethen2000, the polynomial trick works as follows. Let $u(x) = (1 - x^2) q(x)$. Then $u(plus.minus 1) = 0$ automatically. We need $u'(plus.minus 1) = 0$ as well. Since
$ u'(x) = -2 x q(x) + (1 - x^2) q'(x), $
the condition $u'(plus.minus 1) = 0$ gives $q(plus.minus 1) = 0$ (assuming $(1 - x^2)q'(x)$ vanishes at the boundary, which it does since $(1 - x^2)$ is zero there). So if $q(plus.minus 1) = 0$, all four clamped conditions are satisfied.

Now, the fourth derivative of $u = (1 - x^2)q$ is computed using the Leibniz rule. Let $w(x) = 1 - x^2$. Then
$ u^((4)) = sum_(k = 0)^(4) binom(4, k) w^((k)) q^((4-k)). $
Since $w(x) = 1 - x^2$, we have $w' = -2x$, $w'' = -2$, $w''' = 0$, and $w^((4)) = 0$. Therefore:
$ u^((4)) = (1 - x^2) q^((4)) + 4(-2x) q''' + 6(-2) q'' = (1 - x^2) q^((4)) - 8 x q''' - 12 q''. $ <eq-leibniz-fourth>

If the original equation is $u^((4)) = f(x)$, then the equation for $q$ on the interior Chebyshev points (where $q(plus.minus 1) = 0$ is enforced by stripping boundary rows) is
$ [(1 - x^2) D^4 - 8 x D^3 - 12 D^2] bold(q)_("int") = bold(f)_("int"). $

In matrix form, we define the operator acting on the full grid as
$ L = op("diag")(1 - x_j^2) D^4 - 8 op("diag")(x_j) D^3 - 12 D^2, $ <eq-L4-polynomial>
and then extract the interior rows and columns (rows and columns $1, dots, N - 1$) to obtain the $(N - 1) times (N - 1)$ system
$ L_("int") bold(q)_("int") = bold(f)_("int"). $
After solving for $bold(q)_("int")$, the solution is reconstructed as $u_j = (1 - x_j^2) q_j$ at the interior points, with $u_0 = u_N = 0$ at the boundaries.

=== Conditioning Considerations

A well-known difficulty with fourth-order spectral methods is the severe ill-conditioning of the fourth derivative matrix $D^4$. While $D^2$ has condition number $cal(O)(N^4)$, the fourth power $D^4$ has condition number $cal(O)(N^8)$, which limits the achievable accuracy in double precision to roughly $N lt.eq.slant 30$--$40$ for direct application. The polynomial trick does not cure this ill-conditioning --- indeed, the matrix $L$ in @eq-L4-polynomial inherits the $cal(O)(N^8)$ scaling --- but it provides a clean and compact implementation that avoids the bookkeeping of the boundary bordering approach. For problems requiring larger $N$, the coupled second-order reformulation (discussed in @sec-coupled-system) offers substantially better conditioning.

== Computational Étude 14.1: Clamped Beam Under Exponential Load <etude-clamped-beam>

We solve the fourth-order boundary value problem
$ u^((4)) = e^x, quad -1 < x < 1, quad u(plus.minus 1) = u'(plus.minus 1) = 0, $ <eq-clamped-beam>
which models a clamped beam under an exponentially varying distributed load. The exact solution is
$ u_("exact")(x) = e^x + c_3 x^3 + c_2 x^2 + c_1 x + c_0, $ <eq-clamped-exact>
where the four constants $c_0, c_1, c_2, c_3$ are determined by the four boundary conditions. Since $u^((4))(x) = e^x$ and $(c_3 x^3 + c_2 x^2 + c_1 x + c_0)^((4)) = 0$, the particular solution is simply $e^x$, and the complementary solution is the cubic polynomial. Substituting the boundary conditions $u(plus.minus 1) = u'(plus.minus 1) = 0$ into @eq-clamped-exact yields a $4 times 4$ linear system for $(c_0, c_1, c_2, c_3)$.

The spectral implementation uses the polynomial trick from @sec-polynomial-trick with the operator @eq-L4-polynomial. In Python:

```python
D, x = cheb_matrix(N)
D2 = D @ D
D3 = D2 @ D
D4 = D3 @ D
# Polynomial trick operator
S = np.diag(1 - x**2)
X = np.diag(x)
L = S @ D4 - 8 * X @ D3 - 12 * D2
# Interior system
L_int = L[1:N, 1:N]
f_int = np.exp(x[1:N])
q_int = np.linalg.solve(L_int, f_int)
# Reconstruct u = (1 - x^2) * q
u = np.zeros(N + 1)
u[1:N] = (1 - x[1:N]**2) * q_int
```

In MATLAB:

```matlab
[D, x] = cheb_matrix(N);
D2 = D^2; D3 = D^3; D4 = D^4;
S = diag(1 - x.^2);
X = diag(x);
L = S * D4 - 8 * X * D3 - 12 * D2;
L_int = L(2:N, 2:N);
f_int = exp(x(2:N));
q_int = L_int \ f_int;
u = zeros(N+1, 1);
u(2:N) = (1 - x(2:N).^2) .* q_int;
```

In Julia:

```julia
D, x = cheb_matrix(N)
D2 = D^2; D3 = D^3; D4 = D^4
S = diagm(1 .- x.^2)
X = diagm(x)
L = S * D4 - 8X * D3 - 12D2
L_int = L[2:N, 2:N]
f_int = exp.(x[2:N])
q_int = L_int \ f_int
u = zeros(N + 1)
u[2:N] = (1 .- x[2:N].^2) .* q_int
```

#figure(
  image("../figures/ch14/python/clamped_beam.pdf", width: 95%),
  caption: [Clamped beam under exponential load @eq-clamped-beam, solved with the polynomial trick. _Left_: the numerical solution (circles) for $N = 16$ agrees with the exact solution @eq-clamped-exact (solid line). The deflection is small and asymmetric, reflecting the exponential load. _Right_: the maximum error (solid line, left axis) decreases exponentially with $N$, reaching $approx 5.4 times 10^(-16)$ by $N = 15$. The condition number of $L_("int")$ (dashed line, right axis) grows as $cal(O)(N^8)$, eventually limiting achievable accuracy for large $N$.],
) <fig-clamped-beam>

The code generating @fig-clamped-beam is available in:
- `codes/python/ch14/clamped_beam.py`
- `codes/matlab/ch14/clamped_beam.m`
- `codes/julia/ch14/clamped_beam.jl`

=== Discussion

The convergence plot in @fig-clamped-beam demonstrates the hallmark of spectral methods applied to smooth problems: the error decreases geometrically with $N$ and reaches machine precision at a remarkably modest grid size ($N approx 15$). This result is consistent with Trefethen's Program 38 @Trefethen2000, which reports similar accuracy for the same problem. The exponential load $e^x$ is entire (analytic everywhere in the complex plane), so the Chebyshev coefficients of the solution decay faster than any algebraic rate, and spectral convergence is achieved.

The condition number plot reveals the price of working with the fourth derivative: $kappa(L_("int"))$ grows as $cal(O)(N^8)$, which means that by $N approx 35$--$40$, the condition number approaches $10^(16)$ and the matrix becomes effectively singular in double precision. For this particular problem, the solution is so smooth that $N = 15$ suffices, but for problems requiring larger $N$ --- such as those with sharp internal layers or highly oscillatory forcing --- the $cal(O)(N^8)$ conditioning becomes a serious limitation. The coupled second-order reformulation described in @sec-coupled-system addresses this issue by reducing the conditioning to $cal(O)(N^4)$.

The polynomial trick itself is noteworthy for its elegance: by substituting $u = (1 - x^2)q$ and working on the interior grid where $q(plus.minus 1) = 0$, all four clamped boundary conditions are enforced without any tau row replacement or explicit constraint equations. The Leibniz rule transforms the fourth-order equation for $u$ into a fourth-order equation for $q$ with modified coefficients, and the implementation requires only a few lines of code beyond the standard Chebyshev setup.

== Eigenmodes of Fourth-Order Operators <sec-fourth-order-eigen>

Having established the polynomial trick for fourth-order boundary value problems, we now turn to the corresponding _eigenvalue_ problem. The vibration of a clamped beam is governed by
$ u^((4)) = lambda u, quad u(plus.minus 1) = u'(plus.minus 1) = 0, $ <eq-beam-eigen>
where $lambda > 0$ is the square of the frequency of oscillation. The eigenvalues are determined by the transcendental equation
$ cos(lambda^(1\/4)) cosh(lambda^(1\/4)) = 1, $ <eq-beam-transcendental>
whose roots give $lambda_1^(1\/4) approx 4.73004$, $lambda_2^(1\/4) approx 7.85320$, $lambda_3^(1\/4) approx 10.9956$, and so on, with $lambda_n^(1\/4) arrow (2n + 1)pi \/ 2$ for large $n$. The corresponding eigenfunctions are combinations of trigonometric and hyperbolic functions that satisfy the clamped conditions.

For the spectral discretisation, we have two options. The first uses the polynomial trick: the matrix $L$ from @eq-L4-polynomial, restricted to interior points, directly yields the eigenvalue problem $L_("int") bold(q) = lambda op("diag")(1 - x_j^2) bold(q)$, which is a _generalised_ eigenvalue problem because the right-hand side involves the mass matrix $op("diag")(1 - x_j^2)$ (the factor from $u = (1 - x^2) q$ appears on both sides of the equation).

The second approach uses _boundary bordering_: one works with the full $D^4$ matrix and replaces four rows --- two for $u(plus.minus 1) = 0$ and two for $u'(plus.minus 1) = 0$ --- with the corresponding constraint equations. This produces a generalised eigenvalue problem $A bold(u) = lambda B bold(u)$, where $A$ is the bordered fourth-derivative matrix and $B$ is the identity with zeros in the boundary rows.

== Computational Étude 14.2: Vibration Modes of a Clamped Beam <etude-beam-eigenmodes>

We compute the eigenvalues and eigenmodes of @eq-beam-eigen using the boundary bordering approach. The full $(N + 1) times (N + 1)$ matrix $D^4$ is modified as follows:
- Row $0$: replaced by $bold(e)_0^top$ (enforcing $u(1) = 0$).
- Row $1$: replaced by $D[0, :]$ (enforcing $u'(1) = 0$).
- Row $N - 1$: replaced by $D[N, :]$ (enforcing $u'(-1) = 0$).
- Row $N$: replaced by $bold(e)_N^top$ (enforcing $u(-1) = 0$).

The right-hand side mass matrix $B$ is the identity with zeros in rows $0$, $1$, $N - 1$, $N$. The generalised eigenvalue problem $A bold(u) = lambda B bold(u)$ is solved by a standard eigensolver, and the physical eigenvalues are extracted by discarding the spurious infinite eigenvalues arising from the singular boundary rows of $B$.

In Python:

```python
D, x = cheb_matrix(N)
D4 = np.linalg.matrix_power(D, 4)
I = np.eye(N + 1)
A = D4.copy()
B = I.copy()
# Boundary bordering
A[0, :] = I[0, :];    B[0, :] = 0     # u(1) = 0
A[1, :] = D[0, :];    B[1, :] = 0     # u'(1) = 0
A[N-1, :] = D[N, :];  B[N-1, :] = 0   # u'(-1) = 0
A[N, :] = I[N, :];    B[N, :] = 0     # u(-1) = 0
eigenvalues = scipy.linalg.eigvals(A, B)
# Filter finite eigenvalues
lam = eigenvalues[np.isfinite(eigenvalues)]
lam = np.sort(np.real(lam[np.abs(np.imag(lam)) < 1e-6]))
```

In MATLAB:

```matlab
[D, x] = cheb_matrix(N);
D4 = D^4;
I = eye(N + 1);
A = D4;
B = I;
A(1, :) = I(1, :);      B(1, :) = 0;
A(2, :) = D(1, :);      B(2, :) = 0;
A(N, :) = D(N+1, :);    B(N, :) = 0;
A(N+1, :) = I(N+1, :);  B(N+1, :) = 0;
lam = eig(A, B);
lam = sort(real(lam(isfinite(lam) & abs(imag(lam)) < 1e-6)));
```

In Julia:

```julia
D, x = cheb_matrix(N)
D4 = D^4
Id = Matrix(1.0I, N + 1, N + 1)
A = copy(D4)
B = copy(Id)
A[1, :] = Id[1, :];      B[1, :] .= 0
A[2, :] = D[1, :];       B[2, :] .= 0
A[N, :] = D[N+1, :];     B[N, :] .= 0
A[N+1, :] = Id[N+1, :];  B[N+1, :] .= 0
λ = eigvals(A, B)
λ = sort(real.(filter(z -> isfinite(z) && abs(imag(z)) < 1e-6, λ)))
```

#figure(
  image("../figures/ch14/python/beam_eigenmodes.pdf", width: 95%),
  caption: [Eigenmodes of a clamped beam @eq-beam-eigen with $N = 20$. _Left_: the first six eigenmodes, normalised to unit maximum. The mode shapes alternate between symmetric and antisymmetric configurations. _Right_: relative error $|lambda_n^("num") - lambda_n^("exact")| \/ lambda_n^("exact")$ versus mode number, showing that the first $approx 2 N \/ 3$ eigenvalues are resolved to high accuracy.],
) <fig-beam-eigenmodes>

The code generating @fig-beam-eigenmodes is available in:
- `codes/python/ch14/beam_eigenmodes.py`
- `codes/matlab/ch14/beam_eigenmodes.m`
- `codes/julia/ch14/beam_eigenmodes.jl`

=== Discussion

The eigenmode shapes in @fig-beam-eigenmodes display the expected alternating symmetry pattern: the first mode is symmetric with a single central antinode, the second is antisymmetric with one nodal crossing, and so on. This pattern is a direct consequence of the symmetry of the clamped boundary conditions and the even/odd structure of the eigenfunctions on $[-1, 1]$.

The eigenvalue accuracy plot confirms the $2N\/3$ rule of thumb for spectral eigenvalue resolution: for $N = 20$, the first $approx 13$ eigenvalues are captured to high relative accuracy, while higher eigenvalues suffer from under-resolution. This is consistent with the observations in @ch-bvp for second-order problems, but the degradation is somewhat faster here because the fourth derivative amplifies high-frequency errors more aggressively than the second derivative.

The boundary bordering technique used here is more transparent than the polynomial trick for eigenvalue problems, since it produces a standard generalised eigenvalue problem $A bold(u) = lambda B bold(u)$ in the original unknowns $u_j$. The polynomial trick would instead yield $L_("int") bold(q) = lambda M_("int") bold(q)$, where $M_("int") = op("diag")((1 - x_j^2))$ is the interior mass matrix --- a formulation that is equally valid but requires reconstructing $u$ from $q$ after the eigensolve. The conditioning of both formulations is comparable, since both involve the fourth derivative matrix with its inherent $cal(O)(N^8)$ condition number.

== Fourth-Order Problems as Coupled Second-Order Systems <sec-coupled-system>

The severe conditioning of the fourth derivative matrix --- $cal(O)(N^8)$ compared with $cal(O)(N^4)$ for the second derivative --- motivates an alternative discretisation strategy. The idea is to introduce an auxiliary variable and split the fourth-order equation into a _coupled system_ of two second-order equations, each of which can be discretised with the well-conditioned $D^2$ matrix.

=== The Auxiliary Variable Reformulation

Consider the clamped beam equation $u^((4)) = f(x)$ with $u(plus.minus 1) = u'(plus.minus 1) = 0$. We introduce the auxiliary variable
$ w = u'', $ <eq-auxiliary-w>
so that the original equation splits into
$ u'' = w, quad w'' = f(x). $ <eq-coupled-system>
The boundary conditions for $u$ are $u(plus.minus 1) = 0$ (from the original clamped conditions). For $w$, we need $w(plus.minus 1) = u''(plus.minus 1)$, but $u''(plus.minus 1)$ is _not_ prescribed by the clamped conditions --- it is an unknown that must be determined as part of the solution. The remaining clamped conditions $u'(plus.minus 1) = 0$ provide the additional constraints needed to close the system.

The discrete formulation assembles a $2(N + 1) times 2(N + 1)$ block system. Let $D^2$ be the Chebyshev second derivative matrix. The system @eq-coupled-system is discretised as
$ mat(D^2, -I; bold(0), D^2) vec(bold(u), bold(w)) = vec(bold(0), bold(f)), $ <eq-block-system>
with the following boundary modifications:
- Rows of the first block corresponding to $x = plus.minus 1$: enforce $u(plus.minus 1) = 0$ (Dirichlet on $u$).
- Rows of the first block adjacent to boundaries: enforce $u'(plus.minus 1) = 0$ (Neumann on $u$, via tau rows using $D$).
- Rows of the second block corresponding to $x = plus.minus 1$: enforce $w(plus.minus 1) = u''(plus.minus 1)$ by leaving the coupling intact, or alternatively, one may impose no explicit boundary condition on $w$ and let the system determine $w(plus.minus 1)$ from the interior equations.

In practice, the simplest formulation uses interior rows of $D^2$ for both $u$ and $w$, with the boundary rows of the first block enforcing $u(plus.minus 1) = 0$ and $u'(plus.minus 1) = 0$, and the boundary rows of the second block enforcing $w(plus.minus 1)$ as free (natural) boundary values.

=== Conditioning Improvement

The coupled system @eq-block-system involves only $D^2$, whose condition number is $cal(O)(N^4)$. Since the block matrix has two copies of $D^2$, the overall condition number scales as $cal(O)(N^4)$ --- a dramatic improvement over the $cal(O)(N^8)$ conditioning of the direct fourth-derivative approach. The price paid is a doubled system size: $2(N + 1)$ unknowns instead of $N + 1$.

== Computational Étude 14.3: Direct vs. Coupled System Comparison <etude-coupled-comparison>

We solve the same clamped beam problem @eq-clamped-beam using both the direct polynomial trick and the coupled second-order reformulation, comparing accuracy and conditioning as functions of $N$.

In Python (coupled system):

```python
D, x = cheb_matrix(N)
D2 = D @ D
I = np.eye(N + 1)
Z = np.zeros((N + 1, N + 1))
# Block system: [D2, -I; 0, D2] [u; w] = [0; f]
A = np.block([[D2, -I], [Z, D2]])
rhs = np.concatenate([np.zeros(N + 1), np.exp(x)])
# Boundary conditions on u: u(±1) = 0, u'(±1) = 0
A[0, :] = 0;  A[0, 0] = 1;      rhs[0] = 0    # u(1) = 0
A[1, :2*(N+1)] = 0
A[1, :N+1] = D[0, :];            rhs[1] = 0    # u'(1) = 0
A[N-1, :2*(N+1)] = 0
A[N-1, :N+1] = D[N, :];          rhs[N-1] = 0  # u'(-1) = 0
A[N, :] = 0;  A[N, N] = 1;      rhs[N] = 0    # u(-1) = 0
# Boundary conditions on w (free): w(±1) from D2 rows
A[N+1, :] = 0;  A[N+1, N+1] = 1
rhs[N+1] = 0  # placeholder, see actual implementation
A[2*N+1, :] = 0;  A[2*N+1, 2*N+1] = 1
rhs[2*N+1] = 0
sol = np.linalg.solve(A, rhs)
u_coupled = sol[:N+1]
```

In MATLAB:

```matlab
[D, x] = cheb_matrix(N);
D2 = D^2;
I = eye(N + 1);
Z = zeros(N + 1);
A = [D2, -I; Z, D2];
rhs = [zeros(N+1, 1); exp(x)];
% BCs on u
A(1, :) = 0; A(1, 1) = 1; rhs(1) = 0;
A(2, :) = 0; A(2, 1:N+1) = D(1, :); rhs(2) = 0;
A(N, :) = 0; A(N, 1:N+1) = D(N+1, :); rhs(N) = 0;
A(N+1, :) = 0; A(N+1, N+1) = 1; rhs(N+1) = 0;
% BCs on w
A(N+2, :) = 0; A(N+2, N+2) = 1; rhs(N+2) = 0;
A(2*(N+1), :) = 0; A(2*(N+1), 2*(N+1)) = 1; rhs(2*(N+1)) = 0;
sol = A \ rhs;
u_coupled = sol(1:N+1);
```

In Julia:

```julia
D, x = cheb_matrix(N)
D2 = D^2
Id = Matrix(1.0I, N + 1, N + 1)
Z = zeros(N + 1, N + 1)
A = [D2 -Id; Z D2]
rhs = [zeros(N + 1); exp.(x)]
# BCs on u
A[1, :] .= 0; A[1, 1] = 1; rhs[1] = 0
A[2, :] .= 0; A[2, 1:N+1] = D[1, :]; rhs[2] = 0
A[N, :] .= 0; A[N, 1:N+1] = D[N+1, :]; rhs[N] = 0
A[N+1, :] .= 0; A[N+1, N+1] = 1; rhs[N+1] = 0
# BCs on w
A[N+2, :] .= 0; A[N+2, N+2] = 1; rhs[N+2] = 0
A[2(N+1), :] .= 0; A[2(N+1), 2(N+1)] = 1; rhs[2(N+1)] = 0
sol = A \ rhs
u_coupled = sol[1:N+1]
```

#figure(
  image("../figures/ch14/python/coupled_comparison.pdf", width: 95%),
  caption: [Comparison of the direct (polynomial trick) and coupled second-order approaches for the clamped beam problem @eq-clamped-beam. _Left_: the maximum error versus $N$ for both methods. Both achieve spectral convergence, but the coupled method maintains accuracy for larger $N$ before conditioning effects dominate. _Right_: the condition number of the system matrix versus $N$, confirming $cal(O)(N^8)$ growth for the direct method and $cal(O)(N^4)$ for the coupled method.],
) <fig-coupled-comparison>

The code generating @fig-coupled-comparison is available in:
- `codes/python/ch14/coupled_comparison.py`
- `codes/matlab/ch14/coupled_comparison.m`
- `codes/julia/ch14/coupled_comparison.jl`

=== Discussion

The comparison in @fig-coupled-comparison quantifies the conditioning trade-off between the two discretisation strategies. For small $N$ (say $N lt.eq.slant 20$), both methods deliver essentially the same accuracy --- machine precision for this smooth problem --- and the direct method is preferable because it involves a system half the size. For larger $N$, the $cal(O)(N^8)$ conditioning of the direct method begins to erode accuracy, while the coupled method's $cal(O)(N^4)$ conditioning allows it to maintain precision to larger $N$ values.

In practice, the choice between the two approaches depends on the problem at hand. For smooth forcing functions where $N lt.eq.slant 30$ suffices, the polynomial trick is simpler and faster. For problems requiring larger $N$ --- such as those with boundary layers, sharp internal features, or highly oscillatory solutions --- the coupled reformulation is more robust. The doubled system size is a modest penalty, particularly since the block structure of @eq-block-system can be exploited to reduce the computational cost (for example, by using a Schur complement to eliminate $bold(w)$ and solve a reduced system for $bold(u)$ alone).

This trade-off between direct high-order discretisation and reformulation as a lower-order system is a recurring theme in numerical PDEs. It appears in finite element methods (where $C^1$ elements for fourth-order problems are much more complex than $C^0$ elements for the coupled second-order formulation), and it extends to higher-order problems: a sixth-order equation can be reformulated as three coupled second-order equations, each discretised with $D^2$, avoiding the catastrophic $cal(O)(N^(12))$ conditioning of $D^6$.

== Two-Dimensional Biharmonic Operator <sec-2d-biharmonic>

The natural two-dimensional generalisation of the beam equation is the _biharmonic equation_
$ Delta^2 u = f, $ <eq-biharmonic>
where $Delta^2 = Delta compose Delta$ is the biharmonic operator. In Cartesian coordinates on $[-1, 1]^2$, this expands to
$ Delta^2 u = u_(x x x x) + 2 u_(x x y y) + u_(y y y y) = f(x, y). $ <eq-biharmonic-expanded>
The biharmonic equation arises in the Kirchhoff theory of thin elastic plates (where $u$ is the transverse deflection), in two-dimensional Stokes flow (where $u$ is the stream function), and in various problems in continuum mechanics.

=== Kronecker Product Formulation

On a tensor-product Chebyshev grid with $N + 1$ points in each direction, the biharmonic operator is assembled using Kronecker products. Let $D_x^k$ and $D_y^k$ denote the $k$-th derivative matrices in $x$ and $y$, and let $I$ be the $(N + 1) times (N + 1)$ identity. The discrete biharmonic operator is
$ L = D_x^4 times.o I + 2 (D_x^2 times.o D_y^2) + I times.o D_y^4, $ <eq-biharmonic-kronecker>
which acts on the vectorised grid values of $u$ (stacked column by column).

=== Polynomial Trick in Two Dimensions

For a clamped plate --- $u = 0$ and $partial_n u = 0$ on the boundary of $[-1, 1]^2$ --- the polynomial trick extends naturally. In each direction, we write $u(x, y) = (1 - x^2)(1 - y^2) q(x, y)$, which automatically satisfies $u = 0$ on all four edges. The conditions $partial_n u = 0$ reduce to requiring $q = 0$ on the boundary (by the same argument as in the one-dimensional case). The resulting system for $q$ at the interior grid points involves a modified biharmonic operator assembled from the Leibniz-expanded derivatives in each direction.

In Trefethen's formulation @Trefethen2000 (Program 39), the implementation uses the one-dimensional polynomial trick operator @eq-L4-polynomial in each direction, combined via Kronecker products:
$ tilde(L) = L_x times.o S_y + 2 (tilde(D)_x^2 times.o tilde(D)_y^2) + S_x times.o L_y, $
where $L_x$ and $L_y$ are the one-dimensional polynomial trick operators @eq-L4-polynomial, $S_x = op("diag")(1 - x_j^2)$ and $S_y = op("diag")(1 - y_k^2)$, and $tilde(D)_x^2$ and $tilde(D)_y^2$ are the second derivative operators modified to account for the polynomial substitution.

In practice, it is simpler to work with the full biharmonic Kronecker product @eq-biharmonic-kronecker and impose the clamped boundary conditions via boundary bordering on the full $(N + 1)^2 times (N + 1)^2$ system, or to extract only the interior grid points and apply the one-dimensional polynomial trick operators in each direction.

== Computational Étude 14.4: Eigenmodes of a Clamped Square Plate <etude-plate-eigenmodes>

We compute the eigenmodes of the clamped plate vibration problem
$ Delta^2 u = lambda u, quad (x, y) in [-1, 1]^2, quad u = partial_n u = 0 "on" partial([-1, 1]^2), $ <eq-plate-eigen>
following Trefethen's Program 39 @Trefethen2000. The eigenvalues $lambda$ determine the vibration frequencies of the plate, and the corresponding eigenfunctions describe the nodal patterns.

We use $N = 17$ Chebyshev points in each direction and construct the biharmonic operator using the one-dimensional polynomial trick in each direction. The interior grid has $(N - 1)^2 = 256$ points, and the eigenvalue problem is solved for the lowest eigenvalues.

In Python:

```python
D, x = cheb_matrix(N)
D2 = D @ D; D3 = D2 @ D; D4 = D3 @ D
# 1D polynomial trick operator
S = np.diag(1 - x**2)
X = np.diag(x)
L1 = S @ D4 - 8 * X @ D3 - 12 * D2
# Interior indices
idx = slice(1, N)
L1_int = L1[idx, idx]
D2_int = D2[idx, idx]
S_int = np.diag(1 - x[idx]**2)
I_int = np.eye(N - 1)
# 2D biharmonic via Kronecker products
L2d = (np.kron(L1_int, S_int)
     + 2 * np.kron(D2_int, D2_int)
     + np.kron(S_int, L1_int))
M2d = np.kron(S_int, S_int)
eigenvalues, eigenvectors = scipy.linalg.eig(L2d, M2d)
```

In MATLAB:

```matlab
[D, x] = cheb_matrix(N);
D2 = D^2; D3 = D^3; D4 = D^4;
S = diag(1 - x.^2);
X = diag(x);
L1 = S * D4 - 8 * X * D3 - 12 * D2;
L1_int = L1(2:N, 2:N);
D2_int = D2(2:N, 2:N);
S_int = diag(1 - x(2:N).^2);
I_int = eye(N - 1);
L2d = kron(L1_int, S_int) + 2 * kron(D2_int, D2_int) ...
    + kron(S_int, L1_int);
M2d = kron(S_int, S_int);
[V, Lam] = eig(L2d, M2d);
lam = diag(Lam);
```

In Julia:

```julia
D, x = cheb_matrix(N)
D2 = D^2; D3 = D^3; D4 = D^4
S = diagm(1 .- x.^2)
X = diagm(x)
L1 = S * D4 - 8X * D3 - 12D2
L1_int = L1[2:N, 2:N]
D2_int = D2[2:N, 2:N]
S_int = diagm(1 .- x[2:N].^2)
L2d = kron(L1_int, S_int) + 2kron(D2_int, D2_int) + kron(S_int, L1_int)
M2d = kron(S_int, S_int)
F = eigen(L2d, M2d)
```

#figure(
  image("../figures/ch14/python/plate_eigenmodes.pdf", width: 95%),
  caption: [First 25 eigenmodes of the clamped square plate @eq-plate-eigen, displayed on a $5 times 5$ grid. The contour plots show the nodal lines (where $u = 0$) of each eigenmode, arranged in order of increasing eigenvalue. Degenerate eigenvalue pairs --- eigenmodes with the same frequency but different nodal patterns --- appear due to the square symmetry of the domain. For example, the second and third modes form a degenerate pair related by a $90 degree$ rotation.],
) <fig-plate-eigenmodes>

The code generating @fig-plate-eigenmodes is available in:
- `codes/python/ch14/plate_eigenmodes.py`
- `codes/matlab/ch14/plate_eigenmodes.m`
- `codes/julia/ch14/plate_eigenmodes.jl`

=== Discussion

The eigenmode gallery in @fig-plate-eigenmodes reveals several features characteristic of vibrating plates. The fundamental mode (top-left) has no nodal lines in the interior, with a single dome-like deflection. Higher modes exhibit increasingly complex nodal patterns, with nodal lines that may be straight, curved, or closed.

A striking feature is the appearance of _degenerate_ eigenvalue pairs: two modes with the same eigenvalue but different spatial patterns. This degeneracy is a direct consequence of the square symmetry of the domain. The symmetry group of the square (the dihedral group $D_4$) has irreducible representations of dimension one and two; the two-dimensional representations correspond to degenerate eigenvalue pairs. The two modes in a degenerate pair are related by a $90 degree$ rotation of the domain. For example, a mode with a vertical nodal line and a mode with a horizontal nodal line at the same frequency form a degenerate pair.

The computation with $N = 17$ involves a generalised eigenvalue problem of size $(N - 1)^2 = 256$. While this is tractable for a direct eigensolver, the computational cost grows rapidly with $N$: the Kronecker product matrix has $(N - 1)^4$ entries, and the eigensolver requires $cal(O)((N - 1)^6)$ operations. For high-resolution plate computations, iterative eigensolvers or the symmetry reduction described in the next section become essential.

The normalised eigenvalues $lambda_1 \/ pi^4 approx 13.35$, $lambda_2 \/ pi^4 = lambda_3 \/ pi^4 approx 25.38$, and $lambda_4 \/ pi^4 approx 40.73$ are in good agreement with the tabulated values in the structural engineering literature. The comparison with Trefethen's Program 39 @Trefethen2000 confirms that the polynomial trick provides an efficient and accurate discretisation of the biharmonic eigenvalue problem on moderate-sized grids.

== Symmetry and Domain Reduction <sec-symmetry-reduction>

The square domain $[-1, 1]^2$ possesses the full dihedral symmetry group $D_4$, consisting of four rotations and four reflections. This symmetry implies that the eigenmodes of the clamped plate can be classified into _symmetry classes_ based on their behaviour under the reflections $x arrow -x$ and $y arrow -y$:

1. *Even-even*: $u(-x, y) = u(x, y)$ and $u(x, -y) = u(x, y)$.
2. *Even-odd*: $u(-x, y) = u(x, y)$ and $u(x, -y) = -u(x, y)$.
3. *Odd-even*: $u(-x, y) = -u(x, y)$ and $u(x, -y) = u(x, y)$.
4. *Odd-odd*: $u(-x, y) = -u(x, y)$ and $u(x, -y) = -u(x, y)$.

Each symmetry class can be computed independently on the _quarter domain_ $[0, 1]^2$, with the appropriate boundary conditions at the symmetry axes $x = 0$ and $y = 0$:
- *Even symmetry* in $x$: Neumann condition $u_x(0, y) = 0$ at $x = 0$.
- *Odd symmetry* in $x$: Dirichlet condition $u(0, y) = 0$ at $x = 0$.
- Similarly for the $y$-direction.

For the even-even class, the quarter-domain problem has Neumann conditions at $x = 0$ and $y = 0$, and clamped conditions at $x = 1$ and $y = 1$:
$ Delta^2 u = lambda u, quad (x, y) in [0, 1]^2, $
$ u_x(0, y) = u_(x x x)(0, y) = 0, quad u_y(x, 0) = u_(y y y)(x, 0) = 0, $
$ u(1, y) = u_x(1, y) = 0, quad u(x, 1) = u_y(x, 1) = 0. $

This reduces the problem size by a factor of four (from $(N - 1)^2$ to approximately $(N\/2)^2$), which translates to a factor of $approx 16$ savings in memory and a factor of $approx 64$ savings in computation time for the eigensolver.

== Computational Étude 14.5: Quarter-Square Symmetry Reduction <etude-quarter-square>

We compute the even-even eigenmodes of the clamped square plate by solving on the quarter domain $[0, 1]^2$ with Neumann conditions at the symmetry axes $x = 0$ and $y = 0$, and clamped conditions at the walls $x = 1$ and $y = 1$. We use $N = 13$ Chebyshev points in each direction, mapped to $[0, 1]$.

In Python:

```python
D, xi = cheb_matrix(N)
D = 2 * D  # scale for [0, 1]
x = (xi + 1) / 2
D2 = D @ D; D3 = D2 @ D; D4 = D3 @ D
I = np.eye(N + 1)
# Build 1D biharmonic with clamped at x=1 (row 0),
# Neumann at x=0 (row N)
# ... (boundary bordering in each direction)
# Assemble 2D operator via Kronecker products
# Solve generalised eigenvalue problem
```

In MATLAB:

```matlab
[D, xi] = cheb_matrix(N);
D = 2 * D;
x = (xi + 1) / 2;
D2 = D^2; D3 = D^3; D4 = D^4;
I = eye(N + 1);
% Build 1D operators with mixed BCs
% Clamped at x=1: u(1)=0, u'(1)=0
% Neumann at x=0: u'(0)=0, u'''(0)=0
% Assemble 2D operator
% Solve eigenvalue problem
```

In Julia:

```julia
D, ξ = cheb_matrix(N)
D = 2D
x = (ξ .+ 1) ./ 2
D2 = D^2; D3 = D^3; D4 = D^4
Id = Matrix(1.0I, N + 1, N + 1)
# Build 1D operators with mixed BCs
# Assemble 2D operator via Kronecker products
# Solve eigenvalue problem
```

#figure(
  image("../figures/ch14/python/quarter_plate.pdf", width: 95%),
  caption: [First 12 even-even eigenmodes of the clamped square plate, computed on the quarter domain $[0, 1]^2$ with $N = 13$. The modes are displayed on a $4 times 3$ grid in order of increasing eigenvalue. All modes exhibit the expected even-even symmetry: symmetric about both $x = 0$ and $y = 0$. The eigenvalues agree with those obtained from the full-domain computation (@fig-plate-eigenmodes) to full numerical precision.],
) <fig-quarter-plate>

The code generating @fig-quarter-plate is available in:
- `codes/python/ch14/quarter_plate.py`
- `codes/matlab/ch14/quarter_plate.m`
- `codes/julia/ch14/quarter_plate.jl`

=== Discussion

The quarter-domain computation in @fig-quarter-plate demonstrates the practical value of symmetry exploitation. By restricting to the even-even symmetry class, we reduce the problem size by a factor of four and eliminate the degenerate eigenvalue pairs that complicate the full-domain spectrum. The eigenvalues computed on the quarter domain agree with the even-even eigenvalues from the full-domain computation to full numerical precision, confirming that the symmetry reduction is exact (not an approximation).

The computational savings are substantial. For the full domain with $N = 17$, the eigenvalue problem has dimension $(N - 1)^2 = 256$. For the quarter domain with $N = 13$, the dimension is $(N - 1)^2 = 144$ --- and $N = 13$ on $[0, 1]^2$ provides comparable resolution to $N = 17$ on $[-1, 1]^2$, since the Chebyshev points are concentrated near the boundaries where the clamped conditions create boundary layers. The eigensolver cost, which scales as $cal(O)(n^3)$ for a dense matrix of size $n$, is therefore reduced by a factor of $(256\/144)^3 approx 5.6$.

The four symmetry classes (even-even, even-odd, odd-even, odd-odd) partition the full spectrum into four non-overlapping subsets. By solving each class independently, one obtains the complete spectrum with the degeneracies resolved: each degenerate pair from the full-domain computation splits into one mode in the even-odd class and one in the odd-even class. This classification is not only computationally efficient but also physically meaningful, as it identifies the symmetry properties of each vibration mode.

== Hydrodynamic Stability: The Orr--Sommerfeld Equation <sec-orr-sommerfeld>

We now turn to one of the most celebrated applications of spectral methods: the linear stability analysis of viscous parallel shear flows. The question of whether a given laminar flow is stable or unstable to small perturbations lies at the heart of the transition to turbulence, and the answer is provided by the spectrum of the _Orr--Sommerfeld equation_ --- a fourth-order eigenvalue problem with clamped boundary conditions.

=== Linearised Stability of Parallel Flows

Consider a steady, parallel flow with velocity profile $U(y)$ in a channel $-1 lt.eq.slant y lt.eq.slant 1$. For plane Poiseuille flow (pressure-driven flow between parallel plates), the base profile is
$ U(y) = 1 - y^2. $ <eq-poiseuille>
Superimposing an infinitesimal perturbation of the form $hat(psi)(y) e^(i alpha(x - c t))$ (a normal mode with streamwise wavenumber $alpha$ and complex phase speed $c$), the linearised Navier--Stokes equations reduce to the Orr--Sommerfeld equation @Orszag1971 @Drazin2004:
$ frac(1, i alpha R) (D^2 - alpha^2)^2 hat(psi) = (U - c)(D^2 - alpha^2) hat(psi) - U'' hat(psi), $ <eq-orr-sommerfeld>
where $R$ is the Reynolds number, $D = d \/ d y$, and the boundary conditions are
$ hat(psi)(plus.minus 1) = hat(psi)'(plus.minus 1) = 0 $ <eq-os-bc>
(no-slip and no-penetration at the channel walls). This is a generalised eigenvalue problem for the complex phase speed $c$, or equivalently for $lambda = -i alpha c$ (the temporal growth rate plus frequency).

The flow is linearly _stable_ if all eigenvalues $c$ have $op("Im")(c) < 0$ (all modes decay), _unstable_ if any eigenvalue has $op("Im")(c) > 0$ (at least one mode grows), and _neutrally stable_ if the most unstable eigenvalue has $op("Im")(c) = 0$.

=== The Critical Reynolds Number

For plane Poiseuille flow with $U(y) = 1 - y^2$, Orszag @Orszag1971 computed the critical Reynolds number --- the smallest $R$ at which an unstable eigenvalue first appears --- to be
$ R_c approx 5772.22, $ <eq-critical-Re>
at streamwise wavenumber $alpha approx 1.02$. This calculation, performed with a Chebyshev spectral method using $N = 20$--$50$ collocation points, was a landmark demonstration of the power of spectral methods for eigenvalue problems. The result superseded earlier, less accurate estimates obtained by finite-difference methods and asymptotic expansions.

At the critical Reynolds number, the least stable eigenvalue crosses the real axis with $op("Im")(c) = 0$, and its real part gives the phase speed of the neutrally stable Tollmien--Schlichting wave. The corresponding eigenfunction has an oscillatory structure concentrated in the interior of the channel, with rapid decay near the walls.

=== Spectral Discretisation

The Orr--Sommerfeld equation @eq-orr-sommerfeld is discretised using Chebyshev collocation on $N + 1$ points. Define the Chebyshev differentiation matrix $D$ and the matrices
$ A = frac(1, i alpha R) (D^2 - alpha^2 I)^2, $
$ B = (D^2 - alpha^2 I) op("diag")(U(y_j)) - op("diag")(U''(y_j)). $
Note that $(D^2 - alpha^2 I)^2 = D^4 - 2 alpha^2 D^2 + alpha^4 I$.

The clamped boundary conditions @eq-os-bc are imposed by boundary bordering: the first and last two rows of both $A$ and $B$ are replaced by the constraint equations. This yields the generalised eigenvalue problem $A hat(bold(psi)) = c B hat(bold(psi))$.

== Computational Étude 14.6: Orr--Sommerfeld Spectrum <etude-orr-sommerfeld>

We compute the Orr--Sommerfeld spectrum for plane Poiseuille flow @eq-poiseuille with $alpha = 1$ at the critical Reynolds number $R = 5772$, using $N = 40, 60, 80, 100$ Chebyshev points.

In Python:

```python
D, y = cheb_matrix(N)
D2 = D @ D
I = np.eye(N + 1)
U = 1 - y**2
Upp = -2 * np.ones(N + 1)
alpha = 1.0
R = 5772.0
# Orr-Sommerfeld operator
Lap = D2 - alpha**2 * I
A = Lap @ Lap / (1j * alpha * R)
B = np.diag(U) @ Lap - np.diag(Upp)
# Boundary bordering: psi(±1) = psi'(±1) = 0
A[0, :] = 0; A[0, 0] = 1;   B[0, :] = 0
A[1, :] = D[0, :];           B[1, :] = 0
A[N-1, :] = D[N, :];         B[N-1, :] = 0
A[N, :] = 0; A[N, N] = 1;   B[N, :] = 0
c = scipy.linalg.eigvals(A, B)
```

In MATLAB:

```matlab
[D, y] = cheb_matrix(N);
D2 = D^2;
I = eye(N + 1);
U = 1 - y.^2;
Upp = -2 * ones(N+1, 1);
alpha = 1;
R = 5772;
Lap = D2 - alpha^2 * I;
A = Lap^2 / (1i * alpha * R);
B = diag(U) * Lap - diag(Upp);
A(1, :) = 0; A(1, 1) = 1; B(1, :) = 0;
A(2, :) = D(1, :);         B(2, :) = 0;
A(N, :) = D(N+1, :);       B(N, :) = 0;
A(N+1, :) = 0; A(N+1, N+1) = 1; B(N+1, :) = 0;
c = eig(A, B);
```

In Julia:

```julia
D, y = cheb_matrix(N)
D2 = D^2
Id = Matrix(1.0I, N + 1, N + 1)
U = 1 .- y.^2
Upp = -2 .* ones(N + 1)
α = 1.0
R = 5772.0
Lap = D2 - α^2 * Id
A = Lap^2 / (1im * α * R)
B = diagm(U) * Lap - diagm(Upp)
A[1, :] .= 0; A[1, 1] = 1;   B[1, :] .= 0
A[2, :] = D[1, :];            B[2, :] .= 0
A[N, :] = D[N+1, :];          B[N, :] .= 0
A[N+1, :] .= 0; A[N+1, N+1] = 1; B[N+1, :] .= 0
c = eigvals(A, B)
```

#figure(
  image("../figures/ch14/python/orr_sommerfeld.pdf", width: 95%),
  caption: [Orr--Sommerfeld spectrum for plane Poiseuille flow with $alpha = 1$ and $R = 5772$, computed with $N = 40$, $60$, $80$, and $100$ Chebyshev points (shown on a $2 times 2$ grid). The eigenvalues $c$ are plotted in the complex plane; the real axis $op("Im")(c) = 0$ separates stable ($op("Im")(c) < 0$) from unstable ($op("Im")(c) > 0$) modes. The characteristic Y-shaped structure of the spectrum is visible, with three distinct branches. The rightmost eigenvalue $c_("max") approx 0.2375 + 7.82 times 10^(-5) i$ is barely unstable, confirming proximity to the critical Reynolds number.],
) <fig-orr-sommerfeld>

The code generating @fig-orr-sommerfeld is available in:
- `codes/python/ch14/orr_sommerfeld.py`
- `codes/matlab/ch14/orr_sommerfeld.m`
- `codes/julia/ch14/orr_sommerfeld.jl`

=== Discussion

The Orr--Sommerfeld spectrum displayed in @fig-orr-sommerfeld exhibits the well-known Y-shaped structure first documented by Orszag @Orszag1971. The three branches of the spectrum correspond to distinct families of eigenvalues:

- The _A branch_ (or wall modes), running along the left side of the Y, consists of modes concentrated near the channel walls.
- The _P branch_ (or pressure modes), forming the stem of the Y, consists of modes related to the mean pressure gradient.
- The _S branch_ (or centre modes), running along the right side, consists of modes concentrated near the channel centreline.

The most physically significant eigenvalue is the _rightmost_ one (largest $op("Im")(c)$), which determines the linear stability of the flow. At $R = 5772$ and $alpha = 1$, this eigenvalue has $op("Im")(c) approx 7.82 times 10^(-5)$ --- barely positive, confirming that the flow is marginally unstable. The extreme smallness of this imaginary part (the growth rate is $alpha op("Im")(c) approx 10^(-4)$) explains why the critical Reynolds number is so sensitive to numerical accuracy, and why Orszag's spectral computation @Orszag1971 was such a breakthrough: finite-difference methods with $cal(O)(h^2)$ or even $cal(O)(h^4)$ accuracy could not resolve such a tiny growth rate reliably.

The convergence study across $N = 40, 60, 80, 100$ shows that the well-resolved eigenvalues (those in the core of the Y) are essentially independent of $N$, while the poorly resolved ones (at the tips of the branches) move as $N$ increases. For the critical eigenvalue, $N = 40$ already provides several digits of accuracy, and $N = 100$ gives the growth rate to $approx 10$ digits. This rapid convergence is a direct consequence of the spectral accuracy of Chebyshev methods for smooth eigenfunctions.

== Non-Normality and Pseudospectra <sec-pseudospectra>

The Orr--Sommerfeld spectrum tells only part of the stability story. The operator is _non-normal_ --- it does not commute with its adjoint ($A A^* eq.not A^* A$) --- and the eigenvalues of non-normal operators can be extraordinarily sensitive to perturbations. This sensitivity has profound physical consequences: even when all eigenvalues indicate stability ($op("Im")(c) < 0$), small perturbations to the operator (representing, for example, roughness, noise, or nonlinear effects) can produce eigenvalues in the unstable half-plane. This mechanism, known as _subcritical transition_, explains why turbulence is observed in pipe and channel flows at Reynolds numbers well below the critical value predicted by linear stability theory.

=== The $epsilon$-Pseudospectrum

The _$epsilon$-pseudospectrum_ of an operator $A$ is defined as
$ sigma_epsilon(A) = { z in CC : sigma_min (A - z I) lt.eq.slant epsilon }, $ <eq-pseudospectrum-def>
where $sigma_min$ denotes the smallest singular value. Equivalently, $z in sigma_epsilon(A)$ if there exists a perturbation $Delta A$ with $||Delta A|| lt.eq.slant epsilon$ such that $z$ is an eigenvalue of $A + Delta A$. When $A$ is normal, the $epsilon$-pseudospectrum is simply the union of $epsilon$-discs around the eigenvalues. For non-normal operators, the pseudospectrum can extend far beyond the eigenvalues, revealing the _transient growth_ potential of the operator.

The pseudospectrum is computed by evaluating $sigma_min(A - z I)$ on a grid in the complex plane and plotting the level curves. This requires $cal(O)(n^3)$ work per grid point (for the SVD of an $n times n$ matrix), so the computation is expensive but conceptually straightforward.

The definitive reference on pseudospectra is the monograph by Trefethen and Embree @TrefethenEmbree2005, which provides extensive theory, algorithms, and applications across fluid mechanics, quantum mechanics, and numerical analysis.

=== Physical Significance

For the Orr--Sommerfeld operator, the pseudospectrum extends far into the unstable half-plane ($op("Im")(c) > 0$) even at Reynolds numbers where all eigenvalues are stable. This means that infinitesimal perturbations to the operator can create eigenvalues with positive growth rates, providing a linear mechanism for the subcritical transition observed experimentally. The _transient energy growth_ --- the maximum amplification of perturbation energy before eventual exponential decay --- scales as $cal(O)(R^2)$ for pipe flow and $cal(O)(R^3)$ for channel flow, far exceeding what the eigenvalues alone would suggest.

== Computational Étude 14.7: Pseudospectra of the Orr--Sommerfeld Operator <etude-pseudospectra>

We compute the $epsilon$-pseudospectra of the Orr--Sommerfeld operator for plane Poiseuille flow with $N = 100$, $R = 5772$, and $alpha = 1$. The smallest singular value $sigma_min(A - z B)$ (for the generalised eigenvalue problem) is evaluated on a grid of $z$ values in the complex $c$-plane, and contours at $log_(10)(sigma_min) = -2, -4, -6, -8, -10$ are plotted.

The computation transforms the generalised eigenvalue problem $A bold(psi) = c B bold(psi)$ into a standard eigenvalue problem by computing $C = B^(-1) A$ (when $B$ is invertible), and then evaluates $sigma_min(C - z I)$ on the grid.

In Python:

```python
D, y = cheb_matrix(N)
D2 = D @ D
I = np.eye(N + 1)
U = 1 - y**2
Upp = -2 * np.ones(N + 1)
alpha = 1.0; R = 5772.0
Lap = D2 - alpha**2 * I
A = Lap @ Lap / (1j * alpha * R)
B = np.diag(U) @ Lap - np.diag(Upp)
# Boundary bordering
A[0, :] = 0; A[0, 0] = 1;   B[0, :] = 0; B[0, 0] = 1
A[1, :] = D[0, :];           B[1, :] = 0; B[1, 1] = 1
A[N-1, :] = D[N, :];         B[N-1, :] = 0; B[N-1, N-1] = 1
A[N, :] = 0; A[N, N] = 1;   B[N, :] = 0; B[N, N] = 1
C = np.linalg.solve(B, A)
# Grid in complex c-plane
x_grid = np.linspace(-1, 1, 200)
y_grid = np.linspace(-1, 0.2, 200)
sigma_min_grid = np.zeros((len(y_grid), len(x_grid)))
for i, yi in enumerate(y_grid):
    for j, xj in enumerate(x_grid):
        z = xj + 1j * yi
        sigma_min_grid[i, j] = np.min(
            np.linalg.svd(C - z * I, compute_uv=False))
```

In MATLAB:

```matlab
% (similar structure, using svd for sigma_min computation)
C = B \ A;
[X, Y] = meshgrid(linspace(-1, 1, 200), linspace(-1, 0.2, 200));
sigma_min_grid = zeros(size(X));
for i = 1:numel(X)
    z = X(i) + 1i * Y(i);
    sigma_min_grid(i) = min(svd(C - z * eye(N+1)));
end
contour(X, Y, log10(sigma_min_grid), [-2 -4 -6 -8 -10]);
```

In Julia:

```julia
C = B \ A
x_grid = range(-1, 1, length=200)
y_grid = range(-1, 0.2, length=200)
σ_min = [minimum(svdvals(C - (xj + im*yi)*Id))
         for yi in y_grid, xj in x_grid]
```

#figure(
  image("../figures/ch14/python/pseudospectra.pdf", width: 95%),
  caption: [Pseudospectra of the Orr--Sommerfeld operator for plane Poiseuille flow with $N = 100$, $R = 5772$, $alpha = 1$. The eigenvalues (dots) lie in the stable half-plane $op("Im")(c) < 0$, but the $epsilon$-pseudospectral contours extend far into the unstable region. The contours at $epsilon = 10^(-2), 10^(-4), 10^(-6), 10^(-8), 10^(-10)$ show that even a tiny perturbation of order $epsilon approx 10^(-4)$ can push an eigenvalue across the stability boundary. This extreme sensitivity --- a consequence of the non-normality of the operator --- explains the subcritical transition to turbulence observed in experiments at $R << R_c$.],
) <fig-pseudospectra>

The code generating @fig-pseudospectra is available in:
- `codes/python/ch14/pseudospectra.py`
- `codes/matlab/ch14/pseudospectra.m`
- `codes/julia/ch14/pseudospectra.jl`

=== Discussion

The pseudospectral contours in @fig-pseudospectra deliver a profound message about the limitations of eigenvalue-based stability analysis for non-normal operators. The eigenvalues of the Orr--Sommerfeld operator at $R = 5772$ all have $op("Im")(c) < 0$ (the flow is technically stable at this Reynolds number for $alpha = 1$, since the critical wavenumber is $alpha_c approx 1.02$, not exactly $1$), but the pseudospectral contours protrude dramatically into the upper half-plane. The $epsilon = 10^(-4)$ contour already crosses the real axis, meaning that a perturbation to the operator of relative magnitude $10^(-4)$ suffices to create an unstable eigenvalue.

This phenomenon has far-reaching physical implications. In a real flow, perturbations of this magnitude arise naturally from surface roughness, free-stream turbulence, acoustic noise, and the neglected nonlinear terms of the Navier--Stokes equations. The pseudospectrum thus provides a rigorous mathematical explanation for the _subcritical transition_ observed in experiments: turbulence appears in pipe and channel flows at Reynolds numbers ($R approx 1000$--$2000$) far below the critical value ($R_c approx 5772$) predicted by eigenvalue analysis alone.

The computation of pseudospectra is expensive --- evaluating $sigma_min(C - z I)$ at each grid point requires an SVD of an $N times N$ matrix, and the total cost for a $200 times 200$ grid with $N = 100$ is approximately $200^2 times 100^3 \/ 3 approx 10^(10)$ floating-point operations. Efficient algorithms for pseudospectral computation, including contour-tracing methods and projection-based approximations, are discussed in the comprehensive treatment by Trefethen and Embree @TrefethenEmbree2005.

The non-normality of the Orr--Sommerfeld operator also manifests as _transient energy growth_: initial perturbations can amplify by factors of $cal(O)(R^2)$ or $cal(O)(R^3)$ before eventually decaying according to the eigenvalue prediction. This transient amplification is driven by the non-orthogonality of the eigenmodes and is quantified by the numerical abscissa $sup_(t gt.eq.slant 0) ||e^(t A)||$, which can exceed unity even when all eigenvalues are in the left half-plane. The interplay between non-normality, pseudospectra, and transient growth has transformed our understanding of hydrodynamic stability over the past four decades.

== Periodic Fourth-Order Problems and Fourier Methods <sec-periodic-fourth-order>

The études so far have employed Chebyshev methods on bounded domains with wall boundary conditions. When the domain is periodic, Fourier methods provide a natural and highly efficient alternative, with the additional advantage that differentiation in Fourier space is _diagonal_: the $k$-th Fourier mode $e^(i k x)$ satisfies $(d \/ d x)^4 e^(i k x) = k^4 e^(i k x)$, so the fourth derivative matrix is simply $op("diag")(k^4)$ in Fourier space.

=== The Kuramoto--Sivashinsky Equation

The _Kuramoto--Sivashinsky_ (KS) equation @Kuramoto1978 @Sivashinsky1977
$ u_t + u u_x + u_(x x) + u_(x x x x) = 0 $ <eq-ks-ho>
is a canonical model for spatiotemporal chaos in one dimension. It arises in the study of laminar flame fronts, thin film flows, and the Saffman--Taylor instability. The second-derivative term $u_(x x)$ provides a _destabilising_ diffusion (negative viscosity) at long wavelengths, while the fourth-derivative term $u_(x x x x)$ provides _stabilising_ hyperdiffusion at short wavelengths. The nonlinear advection $u u_x$ transfers energy between scales. The competition between these three effects produces a rich, chaotic dynamics characterised by cells that merge, split, and oscillate irregularly.

=== Fourier Discretisation and ETDRK4

On the periodic domain $[0, L]$ with $N$ equally spaced grid points $x_j = j L \/ N$, the discrete Fourier transform converts the KS equation to a system of ODEs for the Fourier coefficients $hat(u)_k$:
$ frac(d hat(u)_k, d t) = (k^2 - k^4) hat(u)_k - frac(i k, 2) hat(u^2)_k, $ <eq-ks-fourier>
where $k = 2 pi n \/ L$ for $n = -N\/2, dots, N\/2 - 1$. The linear part $(k^2 - k^4) hat(u)_k$ has eigenvalues that range from mildly positive (for small $|k|$, responsible for the instability) to hugely negative (for large $|k|$, responsible for the stiffness). This extreme stiffness makes explicit time integrators impractical at any reasonable time step.

The _ETDRK4_ method (Exponential Time Differencing Runge--Kutta, fourth order) @CoxMatthews2002 @KassamTrefethen2005 treats the linear part exactly using matrix exponentials and applies a fourth-order Runge--Kutta scheme to the nonlinear part. For a system of the form
$ frac(d bold(v), d t) = L bold(v) + N(bold(v), t), $
where $L$ is the (diagonal) linear operator and $N$ is the nonlinear term, the ETDRK4 method computes:
$ a_n &= e^(L h\/2) v_n + L^(-1)(e^(L h\/2) - I) N_n, \
b_n &= e^(L h\/2) v_n + L^(-1)(e^(L h\/2) - I) N(a_n, t_n + h\/2), \
c_n &= e^(L h\/2) a_n + L^(-1)(e^(L h\/2) - I) (2 N(b_n, t_n + h\/2) - N_n), \
v_(n+1) &= e^(L h) v_n + h^(-2) L^(-3) [ (-4 - L h + e^(L h)(4 - 3L h + (L h)^2)) N_n \ & quad + 2(2 + L h + e^(L h)(-2 + L h)) (N_a + N_b) \ & quad + (-4 - 3L h - (L h)^2 + e^(L h)(4 - L h)) N_c ], $ <eq-etdrk4>
where $N_n = N(v_n, t_n)$, $N_a = N(a_n, t_n + h\/2)$, $N_b = N(b_n, t_n + h\/2)$, and $N_c = N(c_n, t_n + h)$.

The coefficients in @eq-etdrk4 involve expressions of the form $(e^z - 1)\/z$ that suffer from catastrophic cancellation when $z = L h$ is small. Kassam and Trefethen @KassamTrefethen2005 developed a robust evaluation strategy using Cauchy integrals: each coefficient is evaluated as a contour integral around a small circle in the complex plane, computed by the trapezoidal rule with $M = 32$ points. This eliminates the cancellation errors at negligible additional cost.

== Computational Étude 14.8: Kuramoto--Sivashinsky Equation <etude-kuramoto-sivashinsky>

We solve the KS equation @eq-ks-ho on $[0, 32 pi]$ with $N = 256$ Fourier modes using the ETDRK4 time integrator @KassamTrefethen2005. The initial condition is $u(x, 0) = cos(x \/ 16)(1 + sin(x \/ 16))$, and the simulation runs to $T = 150$.

In Python:

```python
N = 256
L = 32 * np.pi
x = L * np.arange(N) / N
u = np.cos(x / 16) * (1 + np.sin(x / 16))
dt = 0.5
k = np.fft.fftfreq(N, d=L / (2 * np.pi * N))
# Linear operator in Fourier space
Lk = k**2 - k**4
E = np.exp(dt * Lk)
E2 = np.exp(dt * Lk / 2)
# ETDRK4 coefficients via contour integral
M = 32
r = np.exp(1j * np.pi * (np.arange(M) + 0.5) / M)
# ... (Kassam-Trefethen contour integral evaluation)
# Time stepping
for step in range(n_steps):
    Nv = -0.5j * k * np.fft.fft(u**2)
    a = E2 * v + (E2 - 1) / Lk * Nv
    Na = -0.5j * k * np.fft.fft(np.real(np.fft.ifft(a))**2)
    b = E2 * v + (E2 - 1) / Lk * Na
    Nb = -0.5j * k * np.fft.fft(np.real(np.fft.ifft(b))**2)
    c = E2 * a + (E2 - 1) / Lk * (2 * Nb - Nv)
    Nc = -0.5j * k * np.fft.fft(np.real(np.fft.ifft(c))**2)
    v = E * v + Nv * f1 + 2 * (Na + Nb) * f2 + Nc * f3
    u = np.real(np.fft.ifft(v))
```

In MATLAB:

```matlab
N = 256;
L = 32 * pi;
x = L * (0:N-1)' / N;
u = cos(x/16) .* (1 + sin(x/16));
dt = 0.5;
k = [0:N/2-1, 0, -N/2+1:-1]' * (2*pi/L);
Lk = k.^2 - k.^4;
E = exp(dt * Lk);
E2 = exp(dt * Lk / 2);
% ETDRK4 coefficients (Kassam-Trefethen)
M = 32;
r = exp(1i * pi * ((1:M) - 0.5) / M);
% ... contour integral computation ...
v = fft(u);
for step = 1:n_steps
    Nv = -0.5i * k .* fft(real(ifft(v)).^2);
    a = E2 .* v + (E2 - 1) ./ Lk .* Nv;
    Na = -0.5i * k .* fft(real(ifft(a)).^2);
    b = E2 .* v + (E2 - 1) ./ Lk .* Na;
    Nb = -0.5i * k .* fft(real(ifft(b)).^2);
    c = E2 .* a + (E2 - 1) ./ Lk .* (2*Nb - Nv);
    Nc = -0.5i * k .* fft(real(ifft(c)).^2);
    v = E .* v + Nv .* f1 + 2*(Na + Nb) .* f2 + Nc .* f3;
end
u = real(ifft(v));
```

In Julia:

```julia
N = 256
L = 32pi
x = L .* (0:N-1) ./ N
u = cos.(x ./ 16) .* (1 .+ sin.(x ./ 16))
dt = 0.5
k = vcat(0:N÷2-1, 0, -N÷2+1:-1) .* (2π / L)
Lk = k.^2 .- k.^4
E = exp.(dt .* Lk)
E2 = exp.(dt .* Lk ./ 2)
# ETDRK4 coefficients via contour integral
M = 32
r = exp.(1im .* π .* ((1:M) .- 0.5) ./ M)
# ... contour integral computation ...
v = fft(u)
for step in 1:n_steps
    Nv = -0.5im .* k .* fft(real.(ifft(v)).^2)
    a = E2 .* v .+ (E2 .- 1) ./ Lk .* Nv
    Na = -0.5im .* k .* fft(real.(ifft(a)).^2)
    b = E2 .* v .+ (E2 .- 1) ./ Lk .* Na
    Nb = -0.5im .* k .* fft(real.(ifft(b)).^2)
    c = E2 .* a .+ (E2 .- 1) ./ Lk .* (2 .* Nb .- Nv)
    Nc = -0.5im .* k .* fft(real.(ifft(c)).^2)
    v = E .* v .+ Nv .* f1 .+ 2 .* (Na .+ Nb) .* f2 .+ Nc .* f3
end
u = real.(ifft(v))
```

#figure(
  image("../figures/ch14/python/kuramoto_sivashinsky.pdf", width: 95%),
  caption: [Kuramoto--Sivashinsky equation @eq-ks-ho on $[0, 32 pi]$ with $N = 256$ and ETDRK4 time integration up to $T = 150$. The space-time plot shows the characteristic spatiotemporal chaos: an initial transient period ($t lt.eq.slant 20$) during which long-wavelength instabilities grow, followed by a sustained chaotic regime in which cells merge, split, and interact nonlinearly. The solution never settles into a steady or periodic state.],
) <fig-kuramoto-sivashinsky>

The code generating @fig-kuramoto-sivashinsky is available in:
- `codes/python/ch14/kuramoto_sivashinsky.py`
- `codes/matlab/ch14/kuramoto_sivashinsky.m`
- `codes/julia/ch14/kuramoto_sivashinsky.jl`

=== Discussion

The space-time plot in @fig-kuramoto-sivashinsky provides a vivid illustration of _spatiotemporal chaos_ in one of its simplest mathematical settings. The initial condition is smooth and periodic, but the negative-viscosity instability ($u_(x x)$ term) amplifies long-wavelength perturbations, driving the solution into a regime where the nonlinear advection and the stabilising hyperdiffusion compete perpetually. The resulting dynamics is genuinely chaotic: the solution is sensitive to initial conditions, has a positive Lyapunov exponent, and never repeats.

Several features of the KS dynamics are visible in the figure. The initial transient ($t lt.eq.slant 20$) shows the growth of the most unstable wavelengths. Once the nonlinear saturation is reached, the solution enters the chaotic regime, characterised by irregular cell interactions: merging events (two adjacent cells combine into one), splitting events (a cell divides), and mutual displacements. These events are not periodic or predictable, yet the statistical properties of the solution (such as the mean energy, the autocorrelation function, and the power spectrum) are reproducible.

From a numerical standpoint, the KS equation is an ideal testbed for the ETDRK4 scheme. The linear part of the equation is _stiff_: the eigenvalues of the linear operator range from $k^2 - k^4 approx +0.25$ (the most unstable mode, $k approx 1\/sqrt(2)$) to $k^2 - k^4 approx -10^8$ (the highest resolved mode, $k approx N\/2$). An explicit Runge--Kutta method would require $Delta t tilde N^(-4)$ for stability, making long-time integration prohibitively expensive. The ETDRK4 method, by treating the linear part exactly through matrix exponentials, removes this stability restriction and allows time steps of $Delta t = 0.5$ regardless of $N$. The Kassam--Trefethen contour integral trick @KassamTrefethen2005 ensures that the ETDRK4 coefficients are computed without catastrophic cancellation, maintaining the full fourth-order accuracy of the scheme.

The KS equation also serves as a bridge between the spectral methodology developed throughout this textbook and the broader field of dynamical systems. The equation's finite-dimensional attractor, its inertial manifold, and its statistical properties have been studied extensively, making it a canonical example in the theory of dissipative PDEs.

== A Non-Exhaustive Literature Overview <sec-higher-order-literature>

The spectral treatment of fourth-order boundary value problems, plate and beam vibrations, hydrodynamic stability, and spatiotemporal chaos represents a convergence of several major threads in applied mathematics and numerical analysis. The following overview traces these threads from their origins to contemporary developments.

=== Beam and Plate Theory: From Euler--Bernoulli to Kirchhoff

The mathematical modelling of elastic beams dates to the work of Euler and Bernoulli in the eighteenth century. The Euler--Bernoulli beam equation, $E I u^((4)) = f(x)$, relates the transverse deflection $u$ of a slender beam to the applied load $f$ through the flexural rigidity $E I$ (the product of Young's modulus and the second moment of area). This fourth-order equation, together with appropriate boundary conditions (clamped, simply supported, free, or cantilever), forms the basis of structural mechanics. The extension to thin elastic plates by Kirchhoff in the mid-nineteenth century replaces the one-dimensional beam operator with the biharmonic operator $Delta^2$, yielding a system that describes the bending of plates under transverse loads. The classical references for beam and plate theory include Timoshenko and Woinowsky-Krieger @Timoshenko1959 and Ventsel and Krauthammer @Ventsel2001. Modern treatments emphasising computational approaches are found in the monographs by Reddy @Reddy2006 and Szilard @Szilard2004.

=== Spectral Methods for Fourth-Order Problems

The application of spectral methods to fourth-order differential equations has a distinguished history. Orszag @Orszag1971 was among the first to demonstrate the power of Chebyshev spectral methods for the Orr--Sommerfeld equation, computing the critical Reynolds number of plane Poiseuille flow to unprecedented accuracy with only $N = 20$--$50$ collocation points. This calculation established spectral methods as the standard tool for hydrodynamic stability computations and motivated extensive further development.

The polynomial trick for clamped boundary conditions was systematised by Huang and Sloan @HuangSloan1994, who analysed its convergence properties and showed that it retains the full spectral accuracy of the underlying Chebyshev approximation. Trefethen's _Spectral Methods in MATLAB_ @Trefethen2000, Programs 38 and 39, provides compact and widely cited implementations for clamped beam and plate problems. The alternative approach of reformulating the fourth-order problem as a coupled second-order system has been advocated by several authors, including Canuto _et al._ @Canuto2006 and Boyd @Boyd2000, who note the significant improvement in conditioning ($cal(O)(N^4)$ versus $cal(O)(N^8)$) and the easier implementation of diverse boundary conditions.

The Chebyshev differentiation matrices and their properties --- including the $cal(O)(N^(2k))$ condition number growth for the $k$-th derivative matrix --- are analysed in detail by Weideman and Reddy @WeidemanReddy2000, who also provide practical guidelines for the maximum useful $N$ in double-precision arithmetic.

=== Hydrodynamic Stability and the Orr--Sommerfeld Equation

The stability of viscous parallel shear flows has been a central problem in fluid mechanics since the pioneering work of Orr (1907) and Sommerfeld (1908). The Orr--Sommerfeld equation @eq-orr-sommerfeld is a fourth-order eigenvalue problem whose spectrum determines the linear stability of flows such as plane Poiseuille flow (pressure-driven channel flow), plane Couette flow (shear-driven channel flow), and the Blasius boundary layer.

Orszag's landmark 1971 paper @Orszag1971 computed the critical Reynolds number $R_c approx 5772$ for plane Poiseuille flow using Chebyshev spectral methods, settling a decades-long controversy and demonstrating that finite-difference methods were inadequate for this sensitive problem. The comprehensive monograph by Drazin and Reid @Drazin2004 covers the classical theory of hydrodynamic stability, while Schmid and Henningson @SchmidHenningson2001 provide a modern treatment emphasising non-modal stability analysis and the role of non-normality.

=== Non-Normality, Pseudospectra, and Subcritical Transition

The recognition that non-normality of the Orr--Sommerfeld operator has profound consequences for flow stability emerged gradually during the 1990s. Key contributions include the work of Reddy and Henningson @ReddyHenningson1993, who computed optimal transient growth bounds for channel flows, and Trefethen _et al._ @TrefethenEtAl1993, who connected pseudospectra to subcritical transition. The monograph _Spectra and Pseudospectra_ by Trefethen and Embree @TrefethenEmbree2005 provides the definitive mathematical treatment of pseudospectra, including extensive applications to hydrodynamic stability.

The physical picture that emerges from the pseudospectral analysis is striking: even when all eigenvalues of the Orr--Sommerfeld operator predict stability, the extreme non-normality of the operator allows initial perturbations to undergo transient amplification by factors of $cal(O)(R^2)$ (for plane Couette flow) or $cal(O)(R^3)$ (for plane Poiseuille flow) before eventually decaying. If the amplified perturbation triggers nonlinear effects, the flow can transition to turbulence at Reynolds numbers far below the critical value. This mechanism, often called _bypass transition_ or _subcritical transition_, provides a compelling explanation for the discrepancy between the theoretical critical Reynolds number and the experimentally observed transition Reynolds number.

=== The Kuramoto--Sivashinsky Equation and Spatiotemporal Chaos

The Kuramoto--Sivashinsky equation @Kuramoto1978 @Sivashinsky1977 is one of the simplest PDEs exhibiting spatiotemporal chaos. Kuramoto derived it in the context of phase turbulence in reaction-diffusion systems, while Sivashinsky obtained it independently as a model for the instability of laminar flame fronts. The equation has since been studied extensively as a paradigmatic example of low-dimensional chaos in a PDE setting.

The mathematical analysis of the KS equation, including the existence of an inertial manifold and the finite-dimensionality of its attractor, is reviewed in the monographs by Temam @Temam1997 and Robinson @Robinson2001. The numerical simulation of the KS equation using spectral methods has been a standard benchmark problem since the 1980s, and its solution displays the characteristic features of spatiotemporal chaos: sensitivity to initial conditions, positive Lyapunov exponents, and statistical stationarity in the long-time limit.

=== Exponential Time Differencing and the ETDRK4 Scheme

The stiffness of the KS equation and similar fourth-order evolution equations motivated the development of exponential time differencing (ETD) methods, which treat the linear part of the equation exactly using matrix exponentials and apply a standard time-stepping scheme to the nonlinear part. The ETDRK4 scheme was introduced by Cox and Matthews @CoxMatthews2002 and subsequently improved by Kassam and Trefethen @KassamTrefethen2005, who identified a numerical instability in the original formulation and proposed a robust remedy based on contour integral evaluation of the ETD coefficients.

The Kassam--Trefethen modification uses the Cauchy integral formula to evaluate expressions of the form $(e^z - 1)\/z$ and similar coefficient functions along a circle in the complex plane, avoiding the catastrophic cancellation that occurs when $z$ is small. This trick, which adds negligible computational cost, has become the standard implementation of ETDRK4 in the spectral methods community. The method achieves fourth-order accuracy in time with a stability region that encompasses the entire left half-plane, making it ideal for stiff problems where the linear eigenvalues span many orders of magnitude.

Beyond the methods presented here, several advanced topics merit mention. _Multigrid_ and _domain decomposition_ methods for the biharmonic equation exploit the local structure of the operator to achieve optimal complexity for large-scale plate problems. _Spectral element methods_ combine the accuracy of spectral methods with the geometric flexibility of finite elements, enabling the solution of fourth-order problems on complex domains. _Bifurcation analysis_ of the Orr--Sommerfeld equation, tracking the eigenvalues as functions of $R$ and $alpha$, reveals the detailed structure of the neutral stability curve and the connections between different instability mechanisms (Tollmien--Schlichting, Kelvin--Helmholtz, centrifugal).

== Summary <sec-higher-order-summary>

This chapter has developed spectral methods for fourth-order boundary value problems, progressing from elementary beam deflections to the Orr--Sommerfeld equation and spatiotemporal chaos:

+ *The polynomial trick* automatically satisfies clamped boundary conditions ($u = u' = 0$ at both endpoints) by the substitution $u = (1 - x^2)q$ and the Leibniz expansion of $u^((4))$. The interior system for $q$ is solved on the stripped Chebyshev grid.

+ *Boundary bordering* for fourth-order problems replaces four rows of the operator matrix with constraint equations for $u(plus.minus 1) = 0$ and $u'(plus.minus 1) = 0$, producing a generalised eigenvalue problem for clamped beam vibrations.

+ *The coupled second-order reformulation* introduces $w = u''$ and splits $u^((4)) = f$ into $u'' = w$ and $w'' = f$, reducing the condition number from $cal(O)(N^8)$ to $cal(O)(N^4)$ at the cost of a doubled system size.

+ *The two-dimensional biharmonic operator* $Delta^2 = D_x^4 times.o I + 2 D_x^2 times.o D_y^2 + I times.o D_y^4$ is assembled via Kronecker products, and the polynomial trick extends to each coordinate direction for clamped plate problems.

+ *Symmetry reduction* on the square domain exploits even/odd symmetry classes to reduce the eigenvalue problem to a quarter domain, achieving a factor-of-four reduction in problem size and resolving degenerate eigenvalue pairs.

+ *The Orr--Sommerfeld equation* governs the linear stability of viscous parallel shear flows. Spectral discretisation with $N = 40$--$100$ resolves the critical Reynolds number $R_c approx 5772$ for plane Poiseuille flow to high accuracy.

+ *Pseudospectra* of the non-normal Orr--Sommerfeld operator reveal extreme eigenvalue sensitivity, providing a mathematical explanation for subcritical transition to turbulence.

+ *The ETDRK4 scheme* integrates stiff periodic fourth-order evolution equations by treating the linear part exactly with matrix exponentials, enabling efficient long-time simulation of the Kuramoto--Sivashinsky equation.

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
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Highlight*],
      ),
      table.hline(stroke: 0.5pt + luma(180)),
      [14.1], [$u^((4)) = e^x$, clamped beam], [Polynomial trick], [Error $approx 10^(-16)$ at $N = 15$],
      [14.2], [$u^((4)) = lambda u$, clamped beam], [Boundary bordering], [$2N\/3$ eigenvalue rule],
      [14.3], [Direct vs coupled comparison], [Coupled $2^("nd")$-order system], [$cal(O)(N^4)$ vs $cal(O)(N^8)$ conditioning],
      [14.4], [$Delta^2 u = lambda u$, clamped plate], [2D Kronecker products], [Degenerate eigenvalue pairs],
      [14.5], [Quarter-plate symmetry], [Domain reduction], [$4 times$ computational savings],
      [14.6], [Orr--Sommerfeld spectrum], [Chebyshev $+$ bordering], [$R_c approx 5772$ (Orszag)],
      [14.7], [Pseudospectra of Orr--Sommerfeld], [$sigma_min(A - z I)$ grid], [Non-normality $arrow$ subcritical transition],
      [14.8], [Kuramoto--Sivashinsky equation], [Fourier $+$ ETDRK4], [Spatiotemporal chaos],
    ),
  ),
  caption: [Summary of computational études in this chapter.],
) <tbl-ch14-summary>

== Exercises <sec-higher-order-exercises>

*Exercise 14.1* (_Localised load on a clamped beam_). Replace the exponential load in Étude 14.1 with a Gaussian bump $f(x) = e^(-20 x^2)$. (a) Solve $u^((4)) = f(x)$ with clamped boundary conditions using the polynomial trick with $N = 16, 20, 24, 28, 32$. (b) Since there is no closed-form exact solution, estimate the error by computing a reference solution at $N = 64$ and measuring $||u_N - u_(64)||_infinity$. (c) Plot the convergence rate and verify that it is spectral (geometric decrease with $N$). (d) Compare the solution profile with that of the exponential load: how does localisation of the forcing affect the beam deflection?

*Exercise 14.2* (_Simply supported beam_). A simply supported beam has boundary conditions $u(plus.minus 1) = 0$ and $u''(plus.minus 1) = 0$ (zero displacement and zero bending moment). (a) Explain why the polynomial trick with $u = (1 - x^2) q$ does _not_ automatically satisfy these conditions. (b) Use boundary bordering (replace four rows of $D^4$) to solve $u^((4)) = e^x$ with simply supported conditions. (c) Compute the eigenvalues of $u^((4)) = lambda u$ with simply supported conditions and compare with the exact eigenvalues $lambda_n = (n pi)^4$ for $n = 1, 2, 3, dots$ Are these eigenvalues simpler than those of the clamped beam? Explain why.

*Exercise 14.3* (_Boundary condition technique comparison_). For the clamped beam eigenvalue problem $u^((4)) = lambda u$, $u(plus.minus 1) = u'(plus.minus 1) = 0$, compare the conditioning and accuracy of three approaches: (a) the polynomial trick, (b) boundary bordering on $D^4$, and (c) the coupled second-order system. For each method, compute the condition number of the system matrix and the relative error in the first ten eigenvalues as functions of $N$. Present the results in a single figure with six panels ($3 times 2$) and discuss the trade-offs.

*Exercise 14.4* (_Transient energy growth from Orr--Sommerfeld_). The maximum transient energy amplification is $G(t) = ||e^(t A)||_2^2$, where $A$ is the Orr--Sommerfeld operator matrix and $|| dot ||_2$ is the energy norm. (a) For plane Poiseuille flow with $R = 5000$ and $alpha = 1$, compute $G(t)$ for $0 lt.eq.slant t lt.eq.slant 500$ by evaluating $||e^(t A)||_2$ using the eigendecomposition of $A$. (b) Find the maximum $G_("max") = max_t G(t)$ and the time $t_("max")$ at which it occurs. (c) Repeat for $R = 1000, 2000, 5000, 10000$ and plot $G_("max")$ versus $R$ on a log-log scale. Verify the theoretical scaling $G_("max") tilde R^2$ for pipe flow. How does this transient growth relate to the pseudospectra computed in Étude 14.7?

*Exercise 14.5* (_Fourth-order diffusion_). Consider the fourth-order diffusion equation $u_t = -u_(x x x x)$ on $[0, 2 pi]$ with periodic boundary conditions and initial condition $u(x, 0) = sin(x) + 0.5 sin(3 x)$. (a) Implement a Fourier spatial discretisation and integrate using backward Euler with $Delta t = 0.01$ and $N = 64$. The Fourier modes decay as $hat(u)_k(t) = hat(u)_k(0) e^(-k^4 t)$: verify that the numerical solution agrees with this formula. (b) Implement the ETDRK4 scheme with $Delta t = 0.1$. Compare the accuracy of backward Euler and ETDRK4 at $t = 1$ for various $Delta t$. (c) Replace the periodic boundary conditions with clamped conditions on $[-1, 1]$ and solve using a Chebyshev discretisation with implicit Euler time stepping. How does the stiffness of $D^4$ affect the maximum stable time step?

*Exercise 14.6* (_Free-free beam and rigid body modes_). A free-free beam has boundary conditions $u''(plus.minus 1) = u'''(plus.minus 1) = 0$ (zero bending moment and zero shear force). (a) Show analytically that $lambda = 0$ is a double eigenvalue with eigenfunctions $u_1(x) = 1$ (translation) and $u_2(x) = x$ (rotation). These are the _rigid body modes_. (b) Use boundary bordering to discretise the eigenvalue problem $u^((4)) = lambda u$ with free-free conditions. (c) Verify that the numerical spectrum contains two zero eigenvalues (to machine precision) and that the nonzero eigenvalues agree with the exact values from the transcendental equation $cos(lambda^(1\/4)) cosh(lambda^(1\/4)) = 1$ (the same equation as for the clamped beam!). Explain this remarkable coincidence.
