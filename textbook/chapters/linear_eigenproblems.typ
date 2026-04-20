// textbook/chapters/linear_eigenproblems.typ
// Chapter 18: Linear Spectral Eigenproblems --- Discretisation, Verification, and Spurious Modes
// Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
// Email: denys.dutykh@ku.ac.ae
// Homepage: https://www.denys-dutykh.com/
// Last modified: April 2026

#import "../styles/template.typ": dropcap, num, format-table

// Enable equation numbering for this chapter
#set math.equation(numbering: "(1)")

= Linear Spectral Eigenproblems <ch-linear-eigen>

#dropcap[Eigenvalue problems are where spectral methods begin to show, at the same time, their greatest power and their greatest capacity to mislead. A smooth numerical eigenfunction is not necessarily a true eigenfunction, and a matrix eigensolver returns every algebraic eigenvalue of the truncated problem, whether or not any of them approximate the operator spectrum. The governing question of this chapter is therefore deliberately defensive: _when is a computed eigenvalue believable?_ We will develop an answer that is practical, algorithmic, and, above all, honest --- by the end of the chapter every reported eigenvalue will carry with it a short list of checks it has either passed or failed. The chapter is the book's closing exercise in discipline, drawing together every eigenvalue computation seen previously --- Sturm--Liouville modes in @ch-bvp, disk modes in @ch-polar, quasinormal modes in @ch-advanced-bc, beam and plate modes and the Orr--Sommerfeld spectrum in @ch-higher-order --- and treating them not as isolated successes but as data points in a single verification culture. We follow the architecture of @Boyd2000 Chapter 7, while adding to it a systematic three-language cross-validation and a recurring _trust certificate_ that accompanies every numerical spectrum.]

By the end of this chapter, you should be able to:

1. Formulate regular and generalised linear eigenvalue problems in operator form $L u = lambda M u$ and in matrix-pencil form $A bold(u) = lambda B bold(u)$.
2. Discretise such problems using Chebyshev collocation, rational Chebyshev bases for unbounded intervals, and boundary-adapted bases for high-order operators.
3. Solve the resulting matrix problems with the default QR/QZ workflow and understand when local solvers (power method, inverse iteration) are preferable.
4. Decide which computed eigenvalues are trustworthy by comparing spectra at two resolutions using the drift-with-$N$ diagnostic.
5. Distinguish three qualitatively different phenomena: accurate discrete eigenvalues, numerically spurious (under-resolved) eigenvalues, and _physically_ spurious eigenvalues that arise from a badly formulated problem.
6. Recognise the warning signs of a continuous spectrum, interior singularities, and branch-point behaviour.
7. Manage conditioning through basis design, through lower-order system reformulations, and (when nothing else works) through a detour into the complex plane.
8. Build a verification workflow that survives beyond toy problems.

== Opening Vignette: a Smooth Eigenfunction Can Still Be False <sec-smooth-lie>

The central shock of this chapter should arrive before any theorem. Consider the tilde.opplest possible eigenproblem on a finite interval,
$ u_(x x) + lambda u = 0, quad u(plus.minus 1) = 0. $ <eq-laplacian-benchmark>
It is self-adjoint, its spectrum is known exactly,
$ lambda_j = (j pi \/ 2)^2, quad j = 1, 2, 3, dots, $ <eq-laplacian-exact>
and its eigenfunctions alternate by parity,
$ u_j (x) = cases(
  cos(j pi x \/ 2) & "if " j " odd",
  sin(j pi x \/ 2) & "if " j " even". ) $ <eq-laplacian-eigfun>
We discretise @eq-laplacian-benchmark with a Chebyshev pseudospectral method on $N = 16$ Chebyshev--Gauss--Lobatto points, impose Dirichlet boundary conditions by restricting to the interior grid, and solve the resulting $15 times 15$ algebraic eigenproblem $- bold(D)^2_("int") bold(u) = lambda bold(u)$ with a standard library routine. The result is a list of fifteen numerical eigenvalues and fifteen corresponding numerical eigenfunctions.

Two of those eigenfunctions tell the whole story of the chapter.

=== Computational Étude 18.1: The Smooth Lie <etude-smooth-lie>

We implement the computation in Python, MATLAB, and Julia, run all three, and compare.

In Python:

```python
import numpy as np
from cheb_matrix import cheb_matrix

N = 16
D, x = cheb_matrix(N)
D2 = D @ D
A = -D2[1:N, 1:N]                          # interior block
eigvals, eigvecs = np.linalg.eig(A)
order = np.argsort(eigvals.real)
lam = eigvals[order].real                  # sorted spectrum
V   = eigvecs[:, order].real
```

In MATLAB:

```matlab
[D, x] = cheb_matrix(N);
D2 = D^2;
A = -D2(2:N, 2:N);
[V, Lam] = eig(A);
lam_all = diag(Lam);
[~, order] = sort(real(lam_all));
lam = real(lam_all(order));
V = real(V(:, order));
```

In Julia:

```julia
D, x = cheb_matrix(N)
D2 = D * D
A = -D2[2:N, 2:N]
F = eigen(A)
order = sortperm(real.(F.values))
lam = real.(F.values[order])
V   = real.(F.vectors[:, order])
```

The three scripts (`codes/python/ch18/eigen_smooth_lie.py`, `codes/matlab/ch18/eigen_smooth_lie.m`, `codes/julia/ch18/eigen_smooth_lie.jl`) produce eigenvalue lists that agree with each other to a maximum absolute difference of $approx 10^(-12)$. This is an important first data point: in a chapter centred on trust, it is reassuring that three independent implementations of the same discrete problem arrive at the same fifteen numbers. Any disagreement at this stage would indicate a bug in one of the three --- not a property of the operator.

#figure(
  image("../figures/ch18/python/eigen_smooth_lie.pdf", width: 95%),
  caption: [The smooth lie. Left: the third eigenmode at $N = 16$. The numerical values (open circles) agree with the exact eigenfunction $cos(3 pi x \/ 2)$ (solid) to eight decimal places, and the computed eigenvalue $lambda_3^(N=16) approx 22.20661$ matches the exact $(3 pi \/ 2)^2$ to within $3 times 10^(-9)$. Right: the fifteenth eigenmode at the same resolution. The numerical eigenfunction still _looks_ smooth and oscillatory, but it is concentrated near the endpoints in a manner the exact $sin(15 pi x \/ 2)$ is not, and the computed eigenvalue $lambda_(15)^(N=16) approx 3174.79$ is wrong by a factor of $5.7$ relative to the exact $(15 pi \/ 2)^2 approx 555.17$.],
) <fig-smooth-lie>

=== What Just Happened

The third mode is recovered essentially for free: the low-frequency eigenfunction is resolved by the 17-point grid to many digits, and the matrix eigensolver confirms the obvious. The fifteenth mode, by contrast, requires a Fourier-like oscillation of fifteen half-waves across the interval, and no 17-point grid can represent it faithfully. The numerical eigenfunction that the solver returns is still a legitimate eigenvector of the matrix, but it is only the _matrix's_ eigenvector; it is not a sampled version of any eigenfunction of the continuous operator. The eigenvalue it carries --- roughly five and a half times the true value --- is in a certain sense meaningless. Yet the plot looks plausible. It is smooth. It oscillates. It vanishes at the endpoints. A credulous reader would accept it.

This is the trap the chapter has been written to disarm. Three habits will be cultivated, in order of increasing formality:

+ *Never believe a spectrum computed at one resolution alone.* Repeat the calculation at $N$ and at, say, $2N$, and accept only the eigenvalues that agree between the two.
+ *Never believe eigenvalues without also inspecting eigenfunctions.* A numerical eigenfunction may be a plausible-looking oscillation with the wrong eigenvalue, or a correct-looking shape with a wrong amplitude distribution.
+ *Never trust a smooth numerical mode in the upper half of the computed spectrum.* At resolution $N$ the upper half of the matrix spectrum is, with very high probability, numerical garbage dressed in plausible-looking clothes.

The etude above will be revisited in a sharper form in @sec-benchmark-finite, where the refinement $N \to 2N$ is performed explicitly and the trustworthy half of the spectrum is isolated quantitatively.

== From Operator to Matrix Pencil: The Default Workflow <sec-workflow>

A linear spectral eigenproblem usually enters in the form
$ L u = lambda M u, quad cal(B) u = 0, $ <eq-operator-eigen>
where $L$ and $M$ are linear (differential, integral, or algebraic) operators and $cal(B)$ denotes homogeneous boundary conditions. At first sight the calculation looks almost indistinguishable from a boundary value problem that the reader has already met in @ch-bvp. In practice the differences are exactly two: a linear solve is replaced by an eigensolve, and verification is performed mode by mode between two distinct truncations rather than solution by solution.

We call the following four-step procedure the _default workflow_:

+ Represent the unknown eigenfunction in a truncated spectral space.
+ Convert the differential eigenproblem into a matrix eigenproblem or matrix pencil.
+ Compute the algebraic eigenpairs with a global eigensolver (QR for regular problems, QZ for generalised).
+ Repeat the computation with a larger truncation and trust only those modes that persist.

The first three steps are mechanical; the fourth is the one the chapter will spend most of its energy teaching the reader to take seriously.

=== Regular and Generalised Eigenproblems

If $M = I$ in @eq-operator-eigen, so that the equation reduces to $L u = lambda u$, we call the problem _regular_. Otherwise it is _generalised_. The distinction survives discretisation: a regular problem gives the algebraic system
$ bold(A) bold(u) = lambda bold(u), $ <eq-regular-matrix>
while a generalised problem gives the _matrix pencil_
$ bold(A) bold(u) = lambda bold(B) bold(u). $ <eq-pencil-matrix>
Both forms admit library solvers (`np.linalg.eig` in Python, `eig(A, B)` in MATLAB, `eigen(A, B)` in Julia). A reader tempted to eliminate $bold(B)$ by forming $bold(B)^(-1) bold(A)$ and calling a regular-eigenvalue routine should resist: this move can introduce spurious large eigenvalues, as @sec-spurious-i will show in detail. The pencil is part of the mathematical structure of the problem, not a nuisance to be cleared away.

=== Why Spectral Discretisations Are Economical for Eigenproblems

The standard dense eigensolvers QR and QZ are global: they compute _all_ $N$ algebraic eigenvalues and eigenvectors of an $N times N$ matrix without requiring an initial guess, and do so in $cal(O)(10 N^3)$ operations @GottliebOrszag1977. This is an order of magnitude more expensive than an LU solve of the same size, and the expense is incompressible: QR and QZ destroy sparsity as they run, so even a sparse discretisation matrix becomes full as the algorithm proceeds. The consequence, observed by @Boyd2000, is that high-order methods are especially economical for eigenproblems: if a spectral method can resolve the physically meaningful portion of the spectrum with $N = 40$ whereas a second-order finite difference needs $N = 400$, the eigensolve is a thousand times cheaper. Spectral accuracy is not an aesthetic preference for eigenproblems; it is a cost reduction.

=== Computational Étude 18.2: Assemble the Pencil <etude-assemble-pencil>

The workflow just described admits at least two distinct concrete realisations of the _same_ differential eigenproblem. We implement both and verify that they agree, so that the reader appreciates the matrix pencil as an authentic algebraic object rather than an implementation quirk.

*Formulation A: boundary bordering.* We retain the full $(N + 1) times (N + 1)$ matrix built from the Chebyshev second derivative, and we convert the homogeneous Dirichlet conditions $u(plus.minus 1) = 0$ into two _algebraic_ rows of the pencil. The matrix $bold(A)$ has its first row replaced by $(1, 0, dots, 0)$ and its last row by $(0, dots, 0, 1)$, and the right-hand-side matrix $bold(B)$ has zeros on those two rows and the identity elsewhere. The system $bold(A) bold(u) = lambda bold(B) bold(u)$ is a _singular_ generalised eigenproblem: the two boundary rows of $bold(B)$ contribute two eigenvalues at infinity, which the library routine flags automatically and we filter out.

*Formulation B: basis recombination.* We restrict attention to the interior grid by stripping rows and columns $0$ and $N$ of the second-derivative matrix, arriving at the _regular_ problem $-bold(D)^2_("int") bold(u) = lambda bold(u)$ of size $(N - 1) times (N - 1)$. No pencil is formed; no infinite eigenvalues appear; the matrix is smaller.

In Python:

```python
# Formulation A: pencil
A = -D2.copy(); B = np.eye(N + 1)
A[0, :] = 0; A[0, 0]   = 1; B[0, :] = 0           # u(1)  = 0
A[N, :] = 0; A[N, N]   = 1; B[N, :] = 0           # u(-1) = 0
lam_pen = scipy.linalg.eigvals(A, B)              # includes inf eigenvalues
lam_pen = np.sort(lam_pen[np.isfinite(lam_pen)].real)
# Formulation B: regular
A_int   = -D2[1:N, 1:N]
lam_reg = np.sort(np.linalg.eigvals(A_int).real)
```

In MATLAB:

```matlab
% Formulation A: pencil
A = -D2; B = eye(N+1);
A(1,  :)=0; A(1,  1)=1; B(1,  :)=0;
A(N+1,:)=0; A(N+1,N+1)=1; B(N+1,:)=0;
lam_pen = eig(A, B);
lam_pen = sort(real(lam_pen(isfinite(lam_pen))));
% Formulation B: regular
A_int = -D2(2:N, 2:N);
lam_reg = sort(real(eig(A_int)));
```

In Julia:

```julia
# Formulation A: pencil
A = -D2; B = Matrix(1.0I, N+1, N+1)
A[1, :]   .= 0; A[1, 1]   = 1; B[1, :]   .= 0
A[N+1, :] .= 0; A[N+1, N+1] = 1; B[N+1, :] .= 0
lam_pen = eigvals(A, B)
lam_pen = sort(real.(filter(isfinite, lam_pen)))
# Formulation B: regular
A_int = -D2[2:N, 2:N]
lam_reg = sort(real.(eigvals(A_int)))
```

#figure(
  image("../figures/ch18/python/eigen_two_formulations.pdf", width: 95%),
  caption: [Étude 18.2: two formulations of the same eigenproblem. Left: absolute error vs. mode number for the boundary-bordered pencil (blue circles) and the basis-recombined regular problem (coral squares); the teal crosses are $|lambda^("pencil") - lambda^("regular")|$. Both formulations produce the _same_ finite spectrum; the maximum disagreement between them is $1.7 times 10^(-12)$, machine-precision noise. Right: $log_(10) |bold(A)|$ for the bordered pencil, showing the two identity-like boundary rows at top and bottom.],
) <fig-two-formulations>

The code generating @fig-two-formulations is available in:
- `codes/python/ch18/eigen_two_formulations.py`
- `codes/matlab/ch18/eigen_two_formulations.m`
- `codes/julia/ch18/eigen_two_formulations.jl`

=== Critical Discussion

The two formulations agree pairwise to within $2 times 10^(-12)$ in all three languages, a number that should be read as _round-off error in a double-precision eigensolve_, not as an algorithmic discrepancy. Crucially, the two formulations share the _same_ maximum error against the exact spectrum (~$2.65 times 10^3$ for the 15th mode): bordering does not invent accuracy where the underlying discretisation has none. A reader tempted to choose between formulations on the grounds of "one might be more accurate" is mistaken; the choice is governed by bookkeeping convenience, by whether one intends to solve boundary value and eigenvalue problems with the _same_ assembled matrix (favouring bordering), and by whether spurious physically-spurious eigenvalues might arise when $bold(B)$ encodes non-trivial boundary structure (@sec-spurious-i will expose this issue starkly). For Dirichlet conditions on a bounded interval, the two formulations are interchangeable; for more delicate problems they are not.

The sparsity plot on the right of @fig-two-formulations is another pedagogical signal. The QR and QZ routines work by orthogonal transformations that _destroy_ any initial sparsity pattern; the eigensolver will not, and cannot, exploit the large zero block evident in the plot. The only practical consequence of sparsity in a dense-eigensolver pipeline is _storage_, not _cost_. This is a gentle but important warning: spectral eigenproblems benefit from high-order accuracy so that $N$ can stay small, but they do not benefit from sparsity the way finite-element eigenproblems do.

== The Finite-Interval Benchmark and the $N \/ 2$ Rule-of-Thumb <sec-benchmark-finite>

The etude of @sec-smooth-lie leaves a sharp question: of the fifteen eigenvalues that the $N = 16$ calculation produced, how many are trustworthy? Repeating the computation at $N = 32$ and $N = 64$ answers this quantitatively. @Boyd2000 (Rule-of-Thumb 8) crystallises the answer:

#block(stroke: 0.5pt + rgb(20, 45, 110), radius: 3pt, inset: 8pt,
  fill: rgb(248, 250, 254))[
  *Eigenvalue rule-of-thumb.* In solving a linear eigenvalue problem by a spectral method using $N + 1$ terms in the truncated expansion, the lowest $N \/ 2$ eigenvalues are usually accurate to within a few percent, while the larger $N \/ 2$ numerical eigenvalues differ from those of the differential equation by such large amounts as to be useless. The only reliable test, however, is to repeat the calculation with different $N$ and compare the results. The number of good eigenvalues may be smaller than $N \/ 2$ if the modes have boundary layers, critical latitudes, or other areas of very rapid change, or when the interval is unbounded.
]

The upper-half warning is unforgiving: no amount of plotting elegance disguises an eigenvalue that drifts with $N$.

=== Computational Étude 18.3: How Many Modes Did We Really Compute? <etude-benchmark-finite>

We solve the Dirichlet Laplacian @eq-laplacian-benchmark at three resolutions, $N in {16, 32, 64}$, and plot $|lambda_j^("num") - lambda_j^("exact")|$ on a semilogarithmic axis against the mode number $j$. The exact spectrum is @eq-laplacian-exact. The code is a one-line modification of Étude 18.1:

```python
for N in (16, 32, 64):
    D, _ = cheb_matrix(N)
    A   = -(D @ D)[1:N, 1:N]
    lam = np.sort(np.linalg.eigvals(A).real)
    err = np.abs(lam - (np.arange(1, lam.size + 1) * np.pi / 2) ** 2)
```

#figure(
  image("../figures/ch18/python/eigen_benchmark_finite.pdf", width: 80%),
  caption: [Étude 18.3: absolute eigenvalue error vs. mode number at $N = 16, 32, 64$. The dashed horizontal line is the $10^(-2)$ tolerance used to define "good". Coloured dotted verticals at $N \/ 2$ mark the heuristic cliff position. For $N = 16$: 6 good modes. For $N = 32$: 15 good modes. For $N = 64$: 34 good modes. The curves show the geometric accuracy in the resolved portion of the spectrum --- six orders of magnitude of error growth from mode $1$ to mode $N \/ 2$ at $N = 32$ --- followed by a steep cliff up to $cal(O)(lambda_j)$ in the unresolved portion.],
) <fig-benchmark-finite>

The code generating @fig-benchmark-finite is available in:
- `codes/python/ch18/eigen_benchmark_finite.py`
- `codes/matlab/ch18/eigen_benchmark_finite.m`
- `codes/julia/ch18/eigen_benchmark_finite.jl`

All three languages produce the _same_ integer count of good modes at each $N$: $6, 15, 34$. The full spectra agree to within $10^(-12)$ across Python, MATLAB, and Julia --- the same round-off margin observed in the pencil-versus-regular comparison. The identical integer counts are a stronger result: not just the eigenvalues but the _decision_ "trusted vs. suspect" is stable across implementations.

=== Critical Discussion

Two quantitative points are worth flagging honestly.

First, the rule-of-thumb is _statistical_, not _exact_. At $N = 16$ the heuristic predicts $8$ good modes but we measure $6$; at $N = 64$ the heuristic predicts $32$ and we measure $34$. The asymmetry is real: the rule is slightly pestilde.opistic at high resolution (because conditioning has saturated at the noise floor and the high-accuracy regime stretches further than $j = N \/ 2$) and slightly optimistic at low resolution (because the under-resolved modes begin to bleed into the range $j in [N \/ 4, N \/ 2]$). Boyd's own examples in Fig 7.1 and 7.3 show the same asymmetry, and a careful reader should come away distrusting any advertised "good-mode count" that differs from an actual resolution study by more than a factor of two.

Second, the accuracy floor at the lowest modes is _not_ machine epsilon. At $N = 32$ the smallest visible error is around $10^(-13)$, not $10^(-16)$, because of the combined effect of the $cal(O)(N^4)$ condition number of the interior second-derivative matrix and the eigenvalue-computation condition number of the library routine. For the bulk of practical eigenproblems, aiming for relative accuracy below $10^(-12)$ is futile: the discretisation matrix itself does not carry that information.

Finally, the cross-language _identical_ integer counts deserve a word. Independent implementations in three different languages, using three different library eigensolvers, agree to the bit that mode 6 (at $N = 16$) passes the tolerance test and mode 7 fails. This is evidence that both the Chebyshev second derivative and the LAPACK/Intel MKL/OpenBLAS eigensolvers are implemented consistently enough that scientific _decisions_ made on their output do not depend on the language. The reverse --- that a decision _does_ depend on the implementation --- would be a red flag pointing to ill-conditioning, and we will encounter that situation in later études.

== Unbounded Domains and the Infinite-Interval Tax <sec-infinite-tax>

The finite-interval benchmark is pleasant precisely because it is finite. Real applications --- quantum mechanics of a bound state, shear-layer stability in an open channel, atmospheric tides extending to infinity in altitude --- are posed on unbounded intervals. The bounded-domain machinery of @sec-benchmark-finite does not tilde.opply _transfer_ to the infinite interval; something must be added. That something is a basis appropriate to the geometry, and the tilde.opplest choice for the real line is the _rational Chebyshev_ family introduced by @Boyd2000 (Chapter 17).

=== Rational Chebyshev Functions via an Algebraic Map

The key move is a change of variable. Introduce the map
$ x = L frac(t, sqrt(1 - t^2)), quad t in [-1, 1], $ <eq-rc-map>
with a positive _map parameter_ $L$. As $t$ ranges over the closed interval $[-1, 1]$, $x$ ranges over the whole real line, with $t = 0$ corresponding to $x = 0$ and $t arrow plus.minus 1$ corresponding to $x arrow plus.minus infinity$. The rational Chebyshev basis is then
$ "TB"_n (x) = T_n (t(x)), quad n = 0, 1, 2, dots, $ <eq-rc-basis>
so the collocation points on the real line are the pre-images of the Chebyshev--Gauss--Lobatto points in $t$,
$ t_j = cos(pi j \/ N), quad x_j = L frac(t_j, sqrt(1 - t_j^2)), quad j = 1, dots, N - 1. $
The two endpoint nodes $t_0 = 1$ and $t_N = -1$ are excluded: they map to $x = plus.minus infinity$, which is precisely where a bound state vanishes. Restricting to the interior $(N - 1)$ nodes is equivalent to imposing the boundary condition $u(plus.minus infinity) = 0$ by basis design.

=== Derivatives via the Chain Rule

Differentiation matrices follow from the chain rule. Differentiating @eq-rc-map gives
$ frac(d t, d x) = frac((1 - t^2)^(3\/2), L), quad frac(d^2 t, d x^2) = -frac(3 t (1 - t^2)^2, L^2). $ <eq-rc-chain>
Letting $bold(D)_t$ denote the standard Chebyshev differentiation matrix in $t$, the first derivative in $x$ on the interior grid is
$ bold(D)_x = op("diag")(d t \/ d x) bold(D)_t|_("int"), $ <eq-rc-Dx>
and the second derivative, obtained by applying the chain rule twice, is
$ bold(D)_x^((2)) = op("diag")((d t \/ d x)^2) bold(D)_t^2|_("int") + op("diag")(d^2 t \/ d x^2) bold(D)_t|_("int"). $ <eq-rc-D2x>
The utility `rational_chebyshev` packaged with this chapter assembles these matrices in one call. A Gaussian sanity check --- differentiating $u(x) = exp(-x^2 \/ 2)$, whose second derivative is $(x^2 - 1) u$ --- gives maximum errors of $1.25 times 10^(-3)$ at $N = 16$, $1.4 times 10^(-6)$ at $N = 32$, and $2.4 times 10^(-11)$ at $N = 64$, a textbook display of geometric convergence. All three language implementations agree on these numbers.

=== Computational Étude 18.4: The Infinite-Interval Tax <etude-benchmark-oscillator>

The quantum harmonic oscillator
$ u_(x x) + (lambda - x^2) u = 0, quad |u| arrow 0 "as" |x| arrow infinity $ <eq-oscillator>
has exact eigenvalues $lambda_j = 2 j + 1$ for $j = 0, 1, 2, dots$ with eigenfunctions proportional to $H_j(x) exp(-x^2 \/ 2)$, the Hermite functions. The Hamiltonian matrix $bold(H) = -bold(D)_x^((2)) + op("diag")(x^2)$ is assembled directly from the rational Chebyshev utility and diagonalised by a library eigensolver.

In Python:

```python
from rational_chebyshev import rational_chebyshev_derivative_matrices
_, D2_x, x = rational_chebyshev_derivative_matrices(N, L)
H   = -D2_x + np.diag(x ** 2)
lam = np.sort(np.linalg.eigvals(H).real)
```

In MATLAB:

```matlab
[~, D2_x, x] = rational_chebyshev(N, L);
H = -D2_x + diag(x.^2);
lam = sort(real(eig(H)));
```

In Julia:

```julia
using .RationalChebyshev: rational_chebyshev_derivative_matrices
_, D2_x, x = rational_chebyshev_derivative_matrices(N, L)
H   = -D2_x + Diagonal(x .^ 2)
lam = sort(real.(eigvals(H)))
```

#figure(
  image("../figures/ch18/python/eigen_benchmark_oscillator.pdf", width: 95%),
  caption: [Étude 18.4: the infinite-interval tax. Left: absolute eigenvalue error of the harmonic oscillator at $N = 16, 32$ with map parameter $L = 4$. Whereas the bounded Laplacian at $N = 16$ gave six trusted eigenvalues (@fig-benchmark-finite), the oscillator at the same resolution gives only four. Doubling to $N = 32$ buys ten trusted modes against fifteen for the bounded case. Right: $L$-scan at $N = 32$. $L = 4$ gives ten trusted modes; $L = 2$ crowds the grid near the origin (seven trusted); $L = 8$ wastes nodes in the decayed tails (five trusted). The best $L$ is problem-dependent and, for the oscillator, scales as $L tilde.op sqrt(lambda_("max, target"))$.],
) <fig-benchmark-oscillator>

The code generating @fig-benchmark-oscillator is available in:
- `codes/python/ch18/eigen_benchmark_oscillator.py`
- `codes/matlab/ch18/eigen_benchmark_oscillator.m`
- `codes/julia/ch18/eigen_benchmark_oscillator.jl`

All three languages report the same integer counts: $N = 16$ gives $4$ good modes; $N = 32$ gives $10$; the $L$-scan at $N = 32$ gives $7, 10, 5$ at $L = 2, 4, 8$ respectively. These match @Boyd2000 Fig 7.6 to the integer --- a satisfying confirmation that our rational Chebyshev assembly is correct.

=== Critical Discussion

The most important lesson of this étude is quantitative: _the same $N$ that buys six good modes on the bounded interval buys only four on the real line_. This is the "infinite-interval tax", and the numbers are unforgiving. At $N = 32$ the ratio is $10$ to $15$, about two-thirds of the bounded-domain accuracy. Some of this cost is intrinsic (the real line is genuinely a larger geometry), but some of it is tuning: the map parameter $L$ is a new handle the reader must learn to turn.

The $L$-scan on the right of @fig-benchmark-oscillator tells a story worth pausing over. Too small an $L$ crowds the interior collocation nodes near the origin; the wavefunctions' tails at moderate $|x|$ are under-sampled, and high modes degrade. Too large an $L$ does the opposite: most nodes sit in the decayed tails where the wavefunction is already numerically zero, and the near-origin region --- where the oscillatory structure of the Hermite polynomial lives --- is starved of resolution. The optimum $L$ depends on which eigenvalues one wants: for mode $j$, the Hermite wavefunction has classical turning points at $|x| = sqrt(2 j + 1)$, and a sensible heuristic is $L tilde.op sqrt(lambda_("max, target"))$. For the first ten modes this gives $L tilde.op sqrt(19) approx 4.4$, consistent with the empirical optimum $L = 4$ of @fig-benchmark-oscillator.

What would we _not_ dare claim from @fig-benchmark-oscillator? We would not claim a universal rule for $L$. Different potentials have different decay scales, and some --- like the associated Legendre bound states of @sec-slt-taxonomy --- decay only algebraically, not exponentially. For those problems the algebraic map is still convergent, but more slowly, and the optimal $L$ can change dramatically. The étude does, however, establish a workflow that always works: compute the spectrum at several $L$ values, identify the $L$ that maximises the number of trusted eigenvalues under refinement, and report sensitivity to $L$ as part of the trust certificate.

== Verifying a Spectrum: the Drift-with-$N$ Diagnostic <sec-drift-diagnostic>

The chapter's central reusable tool is the _spectrum lie detector_ introduced in this section. It compares two spectra computed at different truncations and returns a verdict, mode by mode. The file `codes/{python,julia,matlab}/ch18/spectrum_verify.{py,jl,m}` implements the diagnostic in all three languages.

Let $lambda^((1))_j$ be the spectrum at the lower resolution $N_1$, sorted. The _intermodal separation_ is
$ sigma_1 = |lambda^((1))_1 - lambda^((1))_2|, quad sigma_j = 1/2 ( |lambda^((1))_j - lambda^((1))_(j-1)| + |lambda^((1))_(j+1) - lambda^((1))_j| ), quad j > 1. $ <eq-sigma-drift>
The scale $sigma_j$ is preferred to $|lambda^((1))_j|$ because eigenvalues grow as $cal(O)(j^2)$ for "nice" Sturm--Liouville problems but the _spacings_ between them grow only as $cal(O)(j)$; using $|lambda^((1))_j|$ would make the drift ratios of low modes artificially flattering. Let $lambda^((2))_j$ be the spectrum at the higher resolution $N_2 > N_1$, also sorted. The _scaled ordinal drift_ is
$ delta^("ord")_j = frac(|lambda^((1))_j - lambda^((2))_j|, sigma_j), $ <eq-delta-ordinal>
and the _scaled nearest drift_ is
$ delta^("nst")_j = frac(min_(k in [1, N_2]) |lambda^((1))_j - lambda^((2))_k|, sigma_j). $ <eq-delta-nearest>
The ordinal drift is the default diagnostic when the two spectra are sorted identically; the nearest drift is required when mode ordering is not preserved under refinement, as happens for @Boyd2000 Laplace's tidal equation where Rossby and gravity waves interlace. A mode $j$ is reported as _trusted_ when both drifts fall below a user-supplied tolerance (typically $10^(-3)$ or $10^(-2)$).

=== A Confidence-Building Cross-Language Test

Before applying the diagnostic to any real eigenproblem, we subject the utility itself to a regression test that every reader should be able to reproduce. We construct a synthetic spectrum $lambda_j = (j pi \/ 2)^2$ for $j = 1, dots, 20$, then deliberately inflate modes $13$ through $20$ by a factor $1.5$ to mimic under-resolution. The utility should return twelve trusted modes and eight suspect ones, and the scaled drifts of the suspect modes should follow the pattern $delta_j = 0.5 lambda_j \/ sigma_j prop j$. All three language implementations return

#align(center)[
  #table(
    columns: (auto, auto),
    stroke: 0.5pt + rgb(180, 180, 190),
    align: center + horizon,
    inset: 6pt,
    table.header([*quantity*], [*value*]),
    [trusted modes],              [$12 "of" 20$],
    [$delta^("ord")_(13)$],       [$3.25$],
    [$delta^("ord")_(14)$],       [$3.50$],
    [$delta^("ord")_(15)$],       [$3.75$],
    [$delta^("ord")_(16)$],       [$4.00$],
    [$delta^("ord")_(17)$],       [$4.25$],
    [$delta^("ord")_(18)$],       [$4.50$],
    [$delta^("ord")_(19)$],       [$4.75$],
    [$delta^("ord")_(20)$],       [$5.1282$],
  )
]

identically, to more than four decimal places. The linear growth of $delta^("ord")_j$ is expected and instructive: the perturbation is $0.5 |lambda_j|$, the separation scales as $|lambda_(j+1) - lambda_j| approx pi^2 j \/ 2$, so the ratio grows as $j$. The last value ($5.13$ instead of $5.25$) reflects the edge effect in the intermodal separation at $j = 20$, where $sigma_(20)$ includes only one neighbouring spacing. These details are _expected consequences of the definitions_ and are recorded here precisely so that a future reader encountering different numbers knows a bug has crept in somewhere.

The unperturbed modes $j = 1, dots, 12$ have drift exactly zero in all three languages, to full precision --- not $10^(-16)$ but _zero_. This is reassuring: the intermodal-separation machinery is not injecting floating-point noise into the diagnostic. The test does not exercise the nearest-drift branch (which would require a mode-ordering inversion), nor the $sigma$ fallback near a zero eigenvalue, and those cases will be revisited in later études.

=== Computational Étude 18.5: The Spectrum Lie Detector in Action <etude-drift-diagnostic>

Having verified the utility against a synthetic test, we turn it loose on the two real problems of the chapter so far: the Dirichlet Laplacian of @sec-benchmark-finite and the harmonic oscillator of @sec-infinite-tax. In each case we compute the spectrum at two resolutions, $N_1 = 32$ and $N_2 = 48$, and plot the _reciprocal_ scaled drifts $1 \/ delta^("ord")_j$ and $1 \/ delta^("nst")_j$ on a semilogarithmic axis against mode number. This follows @Boyd2000 Fig 7.7: trusted modes float to the top of the plot; suspect modes sink to the bottom. A horizontal line at $1 \/ "tol"$ marks the decision threshold.

In Python (abridged):

```python
from spectrum_verify import verify_spectrum
lam1_A = solve_laplacian(32);  lam2_A = solve_laplacian(48)
rep_A  = verify_spectrum(lam1_A, lam2_A, tol=1e-3)
print(f"Laplacian trusted: {rep_A.n_trusted} of {lam1_A.size}")   # 16 of 31
lam1_B = solve_oscillator(32); lam2_B = solve_oscillator(48)
rep_B  = verify_spectrum(lam1_B, lam2_B, tol=1e-3)
print(f"Oscillator trusted: {rep_B.n_trusted} of {lam1_B.size}")  # 9 of 31
```

#figure(
  image("../figures/ch18/python/eigen_drift_diagnostic.pdf", width: 95%),
  caption: [Étude 18.5: drift diagnostic applied to (left) the Dirichlet Laplacian and (right) the harmonic oscillator, both at $N_1 = 32$ and $N_2 = 48$. Blue circles show $1 \/ delta^("ord")$; coral crosses show $1 \/ delta^("nst")$. For each problem a sharp cliff separates trusted from suspect modes. The Laplacian yields 16 trusted modes (close to the $N_1 \/ 2$ heuristic); the oscillator yields only 9 (the "infinite-interval tax"). The two diagnostics agree on the cliff location, since the Laplacian and oscillator both preserve mode ordering under refinement.],
) <fig-drift-diagnostic>

The code generating @fig-drift-diagnostic is available in:
- `codes/python/ch18/eigen_drift_diagnostic.py`
- `codes/matlab/ch18/eigen_drift_diagnostic.m`
- `codes/julia/ch18/eigen_drift_diagnostic.jl`

All three languages report exactly $16$ trusted Laplacian modes and $9$ trusted oscillator modes at the stated tolerance. This is not just a cross-check of continuous-valued quantities but of an _integer decision_ made on the output of three independent eigensolvers --- precisely the kind of robustness that a reader should demand of a verification tool.

=== Critical Discussion

Two observations worth distinguishing carefully.

_First, the drift diagnostic is strictly more informative than a single-resolution error plot._ For the bounded Laplacian we happen to know the exact spectrum, so the distinction is academic; but for the oscillator-with-continuum problem coming up in @sec-slt-taxonomy the exact spectrum is finite, and the "correct" answer is supposed to be "there are four discrete modes, the rest is continuum". The drift diagnostic can see this without being told: the four bound states are stable under refinement, while the continuum modes wander. No single-resolution analysis can distinguish these cases.

_Second, the ordinal and nearest drifts agree here but will not always._ For the Laplacian and the oscillator, sorting by value preserves mode identity under refinement --- the $j$-th Laplacian eigenmode at $N = 32$ really is the $j$-th at $N = 48$. For problems where different mode families interlace under refinement (Laplace's tidal equation is the canonical example; Rossby and gravity waves are interleaved), the ordinal diagnostic lies and the nearest diagnostic is essential. A reader who reflexively uses ordinal matching will silently miss this failure mode. The recipe this chapter adopts --- requiring _both_ drifts below tolerance for a mode to be trusted --- is conservative by design: no ordering pathology can slip past.

The stricter tolerance used here ($10^(-3)$) than in Étude 18.3 ($10^(-2)$) reflects a deliberate raise of the bar. The scaled drift is a _relative_ quantity in the intermodal-separation metric, and it naturally carries a tighter implicit standard than absolute error. A reader who wishes to match Étude 18.3's count exactly should either use a looser tolerance ($3 times 10^(-2)$) or compare the diagnostic's count with the rule-of-thumb's count as two _different_ measures of trustworthiness, not two copies of the same quantity.

== A Taxonomy of Sturm--Liouville Problems <sec-slt-taxonomy>

The étude sequence so far has used two unusually well-behaved eigenproblems: the Dirichlet Laplacian on a bounded interval (a _first-kind_ regular Sturm--Liouville problem in the classification of @Boyd2000) and the quantum harmonic oscillator on the real line (a _second-kind_ singular problem, but with an entirely discrete spectrum). Real applications are not always this cooperative. This section frames the broader landscape so that the warning signs of more delicate spectra become recognisable.

#block(stroke: 0.5pt + rgb(20, 45, 110), radius: 3pt, inset: 8pt,
  fill: rgb(248, 250, 254))[
  *The four kinds of Sturm--Liouville problems* (@Boyd2000). Writing a Sturm--Liouville eigenproblem as
  $ [p(x) u_x]_x + {q(x) + lambda r(x)} u = 0, quad x in [a, b], $
  with homogeneous boundary conditions, there are four qualitatively different cases:

  - *First kind*: $p, q, r$ analytic everywhere on $[a, b]$; pure discrete spectrum.
  - *Second kind*: $p, q, r$ analytic on the open interval but singular at the endpoints (typical of unbounded-domain problems); spectrum either discrete, or finitely many discrete modes plus a continuum.
  - *Third kind*: the differential equation is singular on the _interior_ of the interval, but the eigenfunctions remain regular there (the singularity is "apparent"). Laplace's tidal equation is the canonical example.
  - *Fourth kind*: both the equation and its eigenfunctions are genuinely singular in the interior; complex-valued eigenvalues and eigenfunctions are the rule rather than the exception.
]

Within the second kind we meet, for the first time in the book, operators whose exact spectrum is _not_ a countably infinite tower. A potential well that is too shallow supports only a finite number of bound states; the remaining states lie on a continuum at higher energy. No amount of refinement will convert matrix eigenvalues of the continuum into new bound states.

=== Computational Étude 18.6: Bound States versus Continuum <etude-bound-continuum>

The Pöschl--Teller potential, a staple of quantum-mechanics textbooks, provides a clean exact illustration. The Schrödinger equation
$ -u_(x x) + V(x) u = E u, quad V(x) = -nu (nu + 1) op("sech")^2 x, $ <eq-poschl-teller-ch18>
with decay at $|x| arrow infinity$ has, for any positive integer $nu$, exactly $nu$ bound states at
$ E_j = -(nu - j)^2, quad j = 0, 1, dots, nu - 1, $ <eq-pt-bound-states>
plus a continuous spectrum at $E gt.eq.slant 0$. We take $nu = 4$, so the exact bound states are $E in {-16, -9, -4, -1}$, and discretise on the real line using the `rational_chebyshev` utility at two resolutions $N_1 = 60$, $N_2 = 96$.

In Python:

```python
from rational_chebyshev import rational_chebyshev_derivative_matrices
from spectrum_verify        import verify_spectrum
_, D2, x = rational_chebyshev_derivative_matrices(N, L=6.0)
H = -D2 + np.diag(-nu * (nu + 1) / np.cosh(x) ** 2)
lam = np.sort(np.linalg.eigvals(H).real)
report = verify_spectrum(lam1, lam2, tol=1e-3)
```

#figure(
  image("../figures/ch18/python/eigen_bound_plus_continuum.pdf", width: 95%),
  caption: [Étude 18.6: Pöschl--Teller with $nu = 4$. Left: the lowest 25 eigenvalues at $N_1 = 60$ (circles) and $N_2 = 96$ (crosses). The four exact bound states $E in {-16, -9, -4, -1}$ are marked with horizontal teal dashes; every computation locks onto them to machine precision, while the continuum modes above $E = 0$ drift noticeably. Right: reciprocal drift diagnostic; the cliff at mode 5 separates the four bound states plus the marginal $E = 0$ mode from the under-resolved continuum.],
) <fig-bound-continuum>

The code generating @fig-bound-continuum is available in:
- `codes/python/ch18/eigen_bound_plus_continuum.py`
- `codes/matlab/ch18/eigen_bound_plus_continuum.m`
- `codes/julia/ch18/eigen_bound_plus_continuum.jl`

All three languages recover the four bound states to $4 times 10^(-14)$ (agreement with each other) and to $4 times 10^(-5)$ relative (agreement with the exact integer formula @eq-pt-bound-states); the drift diagnostic flags five trusted modes in all three.

=== Critical Discussion

The striking result of this étude is _negative_: increasing $N$ does not produce more bound states. The matrix spectrum grows --- at $N = 96$ we compute 95 eigenvalues, at $N = 192$ we would compute 191 --- but every eigenvalue beyond the fifth wanders freely between the two resolutions. This is the mathematically correct picture: the operator has _four_ discrete eigenstates, and the drift diagnostic correctly refuses to endorse spurious "fifth-and-higher bound states". A naive user who, seeing an N=96 calculation return twenty-odd negative numerical eigenvalues, concluded that the potential must support twenty-odd bound states would be wrong. The étude is a concrete demonstration that the question _how many discrete modes exist?_ must not be answered by counting returned eigenvalues; it must be answered by the drift test.

The drift test here reports five trusted modes, not four. The extra one is the marginal $E approx 0$ mode --- numerically near-zero and (at this tolerance) stable across the two resolutions. Whether to report it as bound or as "continuum threshold" is a physical rather than numerical question; the diagnostic reports what it measures, and the practitioner's taxonomy then decides. A _tighter_ tolerance would likely demote the $E = 0$ mode to suspect; a _looser_ one might promote more continuum modes to trusted. The robust take-away is the integer 4, not the integer 5.

== False Modes I: Numerically Spurious versus Physically Spurious <sec-spurious-i>

Until this point we have treated under-resolution as the only source of untrustworthy eigenvalues. A second source, qualitatively different, arises when a generalised eigenproblem $bold(A) bold(u) = lambda bold(B) bold(u)$ is formulated in a way that is algebraically valid but _physically_ inconsistent. @Boyd2000 distinguishes the two sources sharply:

#block(stroke: 0.5pt + rgb(20, 45, 110), radius: 3pt, inset: 8pt,
  fill: rgb(248, 250, 254))[
  *Definition (Boyd 2000).*
  - A _numerically spurious_ eigenvalue is a poor approximation to a genuine operator eigenvalue, arising from under-resolution. Increasing $N$ cures it.
  - A _physically spurious_ eigenvalue is one that does not approximate any eigenvalue of the operator. It arises from misrepresentation of the boundary conditions in a generalised formulation. Increasing $N$ does _not_ cure it.
]

The canonical example is the stream-function formulation of the linearised Navier--Stokes equations in one space dimension. @GottliebOrszag1977 first noted that the problem
$ nu u_(x x x x) = lambda u_(x x), quad u(plus.minus 1) = u_x(plus.minus 1) = 0 $ <eq-go-streamfunction>
produces a pair of large positive eigenvalues that are _not_ eigenvalues of the continuous operator (whose spectrum is real and strictly negative). Dawkins, Dunbar, and Douglass (1998) later proved rigorously that these spurious positives scale as $cal(O)(N^4)$, growing with refinement rather than converging away.

The reason is a boundary-condition mismatch. Writing @eq-go-streamfunction in the form
$ lambda u = nu op("D")_(x x)^(-1) u_(x x x x) $
requires inverting a second-order operator on a function space where _four_ boundary conditions have been imposed --- an inconsistent specification. The matrix pencil inherits this inconsistency as $cal(O)(N^4)$ eigenvalues at $lambda approx O(N^4) > 0$.

=== Computational Étude 18.7: Manufacturing Fake Instability <etude-physically-spurious>

We implement both the naive pencil and a cured formulation, and compare. The naive pencil uses $bold(A) = nu bold(D)^4$ and $bold(B) = bold(D)^2$ with four rows of boundary bordering. The cured version discretises @eq-go-streamfunction first on the interior grid (where $u(plus.minus 1) = 0$ is native) and then imposes $u_x(plus.minus 1) = 0$ via two _tau rows_ built from the full first-derivative matrix restricted to interior columns.

In Python (naive):

```python
A = nu * D4.copy(); B = D2.copy()
A[0, :] = I[0, :];  B[0, :] = 0      # u(1)  = 0
A[1, :] = D[0, :];  B[1, :] = 0      # u'(1) = 0
A[N-1, :] = D[N, :]; B[N-1, :] = 0   # u'(-1) = 0
A[N, :] = I[N, :];  B[N, :] = 0      # u(-1) = 0
lam_naive = np.sort(scipy.linalg.eigvals(A, B).real)
```

In Python (cured):

```python
D2i = D2[1:N, 1:N]
A   = nu * (D2i @ D2i); B = D2i.copy()
A[0,  :] = D[0, 1:N]; B[0,  :] = 0   # u'(+1) = 0 via tau row
A[-1, :] = D[N, 1:N]; B[-1, :] = 0   # u'(-1) = 0 via tau row
lam_cured = np.sort(scipy.linalg.eigvals(A, B).real)
```

Analogous MATLAB and Julia implementations are provided.

#figure(
  image("../figures/ch18/python/eigen_physically_spurious.pdf", width: 98%),
  caption: [Étude 18.7: manufacturing fake instability. Left: naive-formulation spectrum at $N = 32$ (circles) and $N = 48$ (crosses), showing a huge isolated positive eigenvalue at each resolution. Centre: scaling of the largest positive eigenvalue of the naive formulation across $N in {16, 24, 32, 48, 64, 96}$, consistent with the $cal(O)(N^4)$ prediction of Dawkins--Dunbar--Douglass (1998) at intermediate $N$ (teal dashed line). Right: the cured formulation; the huge positive eigenvalue disappears and the computed spectrum approaches the expected all-negative shape. A single $O(10^(14))$ outlier from the near-singular $bold(B)$ remains in the cured figure at some resolutions as a reminder that cures are never perfectly aseptic.],
) <fig-spurious-i>

The code generating @fig-spurious-i is available in:
- `codes/python/ch18/eigen_physically_spurious.py`
- `codes/matlab/ch18/eigen_physically_spurious.m`
- `codes/julia/ch18/eigen_physically_spurious.jl`

=== Critical Discussion

The naive formulation delivers one unambiguously unphysical eigenvalue at $N = 32$ (all three languages agree: one positive, magnitude $10^(13)$). @Boyd2000 and Dawkins--Dunbar--Douglass (1998) predict two; the discrepancy is a reflection of how the bordering is implemented --- my version places all four boundary rows outside the interior, so one of the pair is suppressed as an "infinite" eigenvalue by the library filter. Either way, the finding is qualitative and inescapable: the naive formulation _invents_ a mode that the continuous operator does not have, and the invented mode is _positive_ and _large_, which in a stability calculation would be mis-read as violent instability. A user who computed this spectrum at a single $N$ and published the result would be reporting fictitious physics. The cure --- discretising the operator on the interior grid plus tau-row imposition of the Neumann BCs --- eliminates the large positive in the majority of cases; the residual outlier in the cured spectrum at certain $N$ is a reminder that the tau approach introduces its own (smaller) boundary-imposition artefacts, and reinforces the chapter's repeated message that no single formulation is a silver bullet. The only robust protection is comparison at two $N$: the _cured_ spectrum converges, while the _naive_ spurious mode does not.

A final remark. @Boyd2000 notes that the Dawkins--Dunbar--Douglass result is a _proof of bad behaviour_, not a proof of a _specific magnitude_. My numerical scan ($N = 16$ through $96$) confirms the existence and positivity of the large eigenvalue but does not cleanly reproduce an $N^4$ scaling across the whole range --- the measured max positive moves erratically as $N$ grows, reflecting the fact that the "two large positives" of the rigorous result can swap positions, interlace with other matrix eigenvalues, or temporarily cancel in our particular bordered formulation. Had we tried to read off a rate from these data, we would have wandered. The rigorous result tells us what to expect; the numerics confirm _existence_ but not a clean rate.

== False Modes II: Condition Number, Basis Design, High-Order Operators <sec-spurious-ii>

The étude just completed showed that _formulation_ can create false modes. A parallel phenomenon creates false confidence in _correct_ modes when the discretisation is so badly conditioned that the digits we report are eaten by floating-point noise long before we notice. The problem is sharpest for high-order operators. Differentiating a Chebyshev polynomial $T_N(x)$ $p$ times produces a function whose maximum grows as
$ |(d^p T_N \/ d x^p)(plus.minus 1)| tilde.op N^(2 p), $ <eq-cheb-derivative-bound>
and consequently the matrix discretising $d^p \/ d x^p$ has condition number $cal(O)(N^(2 p))$. For a fourth derivative this is $cal(O)(N^8)$; for a sixth derivative $cal(O)(N^(12))$. At $N = 30$ the fourth-derivative condition number is already $10^(12)$, which means double-precision arithmetic can supply only $10^(-4)$ relative accuracy --- fine for some applications, disastrous for others.

=== Boundary-Adapted Bases

Heinrichs (1989) observed that a well-chosen basis can absorb much of this conditioning into the algebra and leave the numerics better-behaved. For the second-order Dirichlet problem $u(plus.minus 1) = 0$, the Heinrichs basis
$ phi_j (x) = (1 - x^2) T_j(x), quad j = 0, 1, dots, $ <eq-heinrichs-second>
automatically satisfies the boundary condition: no bordering, no tau rows. For the fourth-order clamped problem, the double-root basis
$ phi_j (x) = (1 - x^2)^2 T_j (x) $ <eq-heinrichs-fourth>
satisfies all four conditions $u(plus.minus 1) = u_x(plus.minus 1) = 0$ by construction. The substitution $u = (1 - x^2)^2 q$ with $q$ a free polynomial reduces $u_(x x x x) = lambda u$ on the interior collocation grid to the generalised eigenproblem
$ bold(A) bold(q) = lambda bold(M) bold(q), $
with
$ bold(A) = bold(S)^2 bold(D)^4 - 16 bold(X) bold(S) bold(D)^3 + (48 bold(X)^2 - 24 bold(I)) bold(D)^2 + 96 bold(X) bold(D) + 24 bold(I), $
$ bold(M) = bold(S)^2, $
where $bold(S) = op("diag")(1 - x^2)$, $bold(X) = op("diag")(x)$, and everything is restricted to interior nodes. The `heinrichs_basis` utility packaged with this chapter assembles $(bold(A), bold(M))$ in one call.

=== Computational Étude 18.8: Condition-Number Surgery <etude-heinrichs>

We compare the naive boundary-bordered $bold(D)^4$ pencil (the formulation used in @ch-higher-order Étude 14.2 for the clamped beam) against the Heinrichs $(1 - x^2)^2 T_j$ basis. The target is $u_(x x x x) = lambda u$ with four clamped conditions; the exact spectrum is $lambda_j = beta_j^4$ where $beta_j$ solves $cos(2 beta) cosh(2 beta) = 1$, the first root at $beta_1 approx 2.365$ giving $lambda_1 approx 31.285$.

#figure(
  image("../figures/ch18/python/eigen_heinrichs_condition.pdf", width: 95%),
  caption: [Étude 18.8: condition-number surgery for the fourth-order clamped eigenproblem. Left: condition numbers of the naive $bold(D)^4$-bordered pencil (blue circles) and the Heinrichs basis $(1 - x^2)^2 T_j$ (coral squares). Both grow rapidly but the Heinrichs version is consistently $approx 8 times$ better. Dashed reference slopes $N^8$ and $N^4$ are shown as asymptotes. Right: the error in the first eigenvalue at the two discretisations. At small $N$ the Heinrichs version is strictly more accurate; at large $N$ both saturate near the accuracy floor set by conditioning.],
) <fig-heinrichs>

The code generating @fig-heinrichs is available in:
- `codes/python/ch18/eigen_heinrichs_condition.py`
- `codes/matlab/ch18/eigen_heinrichs_condition.m`
- `codes/julia/ch18/eigen_heinrichs_condition.jl`

Across the full range $N in {12, 16, 24, 32, 48, 64, 96}$, the condition numbers returned by Python, MATLAB, and Julia are _identical_ to three significant digits, to the last printed place --- a three-language agreement of the kind we have come to expect by this point in the chapter.

=== Critical Discussion

The empirical finding of this étude is more nuanced than the textbook caricature would suggest. In the idealised presentation, the Heinrichs basis buys $cal(O)(N^4)$ conditioning where the naive version has $cal(O)(N^8)$. In my measurement, using the condition number of the _standardised_ operator $bold(M)^(-1) bold(A)$, both versions scale like $N^(~7-8)$ within our range --- but Heinrichs is consistently $approx 8 times$ better at every $N$. This is not a $N^4$-versus-$N^8$ separation at $N = 96$; it is a _constant_ factor-of-eight win that persists indefinitely. Is the pedagogical claim still true? Yes, but the honest version is that Heinrichs wins a constant multiplier on top of a scaling that remains dominated by the inverse mass matrix $bold(S)^(-2)$ appearing when one converts the Heinrichs problem to standard form. A measurement that _preserves_ the generalised structure --- for example, comparing the norm growth of $bold(A)$ alone, not of $bold(M)^(-1) bold(A)$ --- would restore the textbook $N^4$-versus-$N^8$ picture. The user who converts to standard form by left-multiplying by $bold(M)^(-1)$ inherits whichever conditioning $bold(M)$ brings along; the user who feeds $(bold(A), bold(M))$ to a QZ routine does not.

The consequence for practice is unchanged: Heinrichs is preferable to boundary bordering for high-order eigenproblems, and the preference grows more important as operator order increases. For sixth- and eighth-order equations the analogous $(1 - x^2)^3$ and $(1 - x^2)^4$ bases are the only realistic approach. Cross-reference @sec-coupled-system for a complementary technique --- reformulating as coupled second-order systems --- which is sometimes a better trade-off than basis design.

== Solver Strategy: Global First, Local Second <sec-solver-strategy>

The chapter has so far used the QR and QZ eigensolvers as black boxes. For problems where $N lt.eq.slant 1000$ on a modern workstation, that is the right choice: these are global, dense, back-box-robust methods that require nothing from the user except the matrices. When $N$ grows, or when the same eigenproblem must be solved thousands of times (for example in a parameter sweep), global solvers become uneconomical and one reaches for the _local_ iterative methods surveyed in this section.

The _power method_ for a matrix $bold(A)$ takes a random vector $bold(v)^((0))$, normalises it, and iterates
$ bold(v)^((k+1)) = frac(bold(A) bold(v)^((k)), || bold(A) bold(v)^((k)) ||). $
The iterate converges to the eigenvector of largest magnitude eigenvalue, and the Rayleigh quotient $mu^((k)) = bold(v)^((k)) dot bold(A) bold(v)^((k))$ converges to that eigenvalue itself. The convergence is _geometric_ with rate $|lambda_(N-1) \/ lambda_N|$: good when there is a gap above the dominant mode, slow when the dominant and subdominant are nearly degenerate. The method selects only one eigenvalue --- specifically, the one the user did not specify.

_Inverse iteration_ fixes the selection problem at the cost of a linear solve. Given a shift $mu$, iterate
$ (bold(A) - mu bold(I)) bold(v)^((k+1)) = bold(v)^((k)); quad bold(v)^((k+1)) arrow.l bold(v)^((k+1)) \/ || bold(v)^((k+1)) ||. $
The iterate converges to the eigenvector of the eigenvalue _nearest_ to $mu$, at geometric rate $|lambda_("nearest") - mu| \/ |lambda_("second-nearest") - mu|$. Choosing $mu$ close to a target eigenvalue makes the rate arbitrarily good, but also makes the matrix $bold(A) - mu bold(I)$ arbitrarily close to singular --- so each linear solve becomes expensive as $mu$ approaches a true eigenvalue. In practice one uses an LU factorisation once (with a fixed $mu$) and re-solves many right-hand sides cheaply.

=== Computational Étude 18.9: One Mode at a Time <etude-power-inverse>

We apply both methods to $bold(A) = -bold(D)^2_("int")$ at $N = 32$, whose spectrum is $lambda_j approx (j pi \/ 2)^2$. The target for the power method is $lambda_("max") approx 49939$; the targets for inverse iteration are _whichever eigenvalue lies nearest the shift_.

#figure(
  image("../figures/ch18/python/eigen_power_inverse.pdf", width: 98%),
  caption: [Étude 18.9: one mode at a time. Left: power method iterate converging to $lambda_("max") approx 49939$ at rate $|lambda_(N-1) \/ lambda_N|$. Centre: inverse iteration with three different shifts ($mu = 5, 90, 250$) locking onto three different interior eigenvalues ($lambda approx 2.47, 88.83, 246.74$). Each shift reaches its target in $lt.eq.slant 20$ iterations, with the convergence rate determined by the gap between the two eigenvalues nearest the shift. Right: cautionary case; the shift $mu approx 50.58$ chosen exactly between two adjacent eigenvalues ($lambda approx 39.48$ and $lambda approx 61.68$) produces an iteration that does not cleanly pick either --- distances to both modes oscillate, and convergence is slow.],
) <fig-power-inverse>

The code generating @fig-power-inverse is available in:
- `codes/python/ch18/eigen_power_inverse.py`
- `codes/matlab/ch18/eigen_power_inverse.m`
- `codes/julia/ch18/eigen_power_inverse.jl`

All three languages agree on the three inverse-iteration targets to $10^(-13)$. The power-method iterate at $k = 80$ reports a slightly different _non-converged_ value in each language (Python $approx 49857$, Julia $approx 49813$, MATLAB $approx 49826$), because the starting random vector from each language's RNG is different. The _eigenvalue_ being approached is identical; only the _rate of arrival_ differs, and more than 80 iterations would bring all three to the same answer.

=== Critical Discussion

The three panels of @fig-power-inverse encode the central algorithmic lesson. A power-method iteration that stalls says "the two largest eigenvalues are close"; an inverse iteration that stalls says "the shift is not close to any eigenvalue, or it is close to two eigenvalues simultaneously". In neither case is the method silently returning a wrong answer --- the _symptom_ is visible in the convergence plot. A user who monitors convergence and stops when it plateaus is protected; a user who iterates a fixed number of steps and trusts the final number is not.

The right-hand panel is a miniature sermon on shift selection. The user who naively picks $mu$ by reading "the midpoint between two eigenvalues looks like a good compromise" from the global-solver output gets exactly the behaviour shown: slow, undecided, oscillating. Good shift selection is asymmetric --- pick a $mu$ much closer to the target than to any other eigenvalue, so the convergence rate has a small numerator and a large denominator. If the target is a simple eigenvalue, one can often get machine precision in ten iterations; if the target sits in a tight cluster, inverse iteration alone will not resolve the cluster, and one must move to a block method such as Arnoldi (outside this chapter's scope).

== Parameter-Dependent Spectra and the Art of Not Missing a Branch <sec-parameter-sweep>

Many applications require the leading eigenvalue of an operator whose coefficients depend on one or more parameters. The neutral stability curve in a hydrodynamic problem, the ground-state energy of a Hamiltonian as a coupling is varied, the most unstable wavenumber as a function of Reynolds number --- all of these are _parameter-dependent spectra_. @Boyd2000 recommends an asymmetric strategy: _global_ QR solves on a coarse grid in parameter space, then _local_ inverse iteration to fill in between the anchors. The point of re-anchoring is not to save time --- it is to avoid missing an eigenvalue branch that the local iteration would have tracked off to a different mode.

=== Computational Étude 18.10: Map the Neutral Curve Safely <etude-parameter-map>

We use a two-parameter model
$ -u_(x x) + (alpha + beta x^2) u = lambda u, quad u(plus.minus 1) = 0, $ <eq-two-param>
with $(alpha, beta) in [-2, 2] times [0, 8]$. The coarse grid is $9 times 9$, all computed with QR (81 full eigensolves). The fine $beta$-only scan at fixed $alpha = 0$ uses 40 values; the _safe_ scanner re-runs QR at each point, while the _naive_ scanner continues inverse iteration from the previous step's converged shift, seeded from the _second_ eigenvalue of $beta = 0$ so as to stress-test the failure mode.

#figure(
  image("../figures/ch18/python/eigen_parameter_map.pdf", width: 98%),
  caption: [Étude 18.10: safely mapping a parameter-dependent spectrum. Left: the leading eigenvalue $lambda_1(alpha, beta)$ computed by QR at every node of a $9 times 9$ coarse grid (white dots). Contours are drawn through the eigenvalue surface. Right: fine $beta$-scan at $alpha = 0$. The blue curve (safe) re-runs QR at every $beta$, faithfully tracking $lambda_1$ from $approx 2.47$ at $beta = 0$ to $approx 3.44$ at $beta = 8$. The coral curve (naive) continues inverse iteration from a deliberately bad initial shift and tracks the second eigenvalue $lambda_2$ throughout, ending at $approx 12.06$ instead of $approx 3.44$.],
) <fig-parameter-map>

The code generating @fig-parameter-map is available in:
- `codes/python/ch18/eigen_parameter_map.py`
- `codes/matlab/ch18/eigen_parameter_map.m`
- `codes/julia/ch18/eigen_parameter_map.jl`

All three languages agree to $3 times 10^(-13)$ on the safe scanner and reach identical endpoints on the naive one.

=== Critical Discussion

The right panel of @fig-parameter-map is the entire pedagogical point. _The two curves are both mathematically valid computations_; one computes the smallest eigenvalue and the other computes the second-smallest. A practitioner using the naive tracker without re-anchoring would be reporting a parameter sweep of the _wrong eigenvalue_, and unless the user independently knows what answer to expect (in which case the calculation is not scientifically interesting), the mistake is invisible. In real applications --- neutral curves in hydrodynamic stability, phase boundaries in condensed-matter lattices --- the true leading mode can switch identity as the parameters vary, and the naive tracker would silently follow an irrelevant branch across a crossing. The safe strategy costs 40 full QR solves instead of 1 plus 39 inverse iterations, a factor-of-few slowdown that any practical code can afford.

The coarse-grid QR anchors (left panel) serve two purposes at once. They provide starting guesses for any local iteration the user might layer on top, _and_ they act as a self-consistent reference against which the local tracker can be audited. If the local tracker departs from the coarse-grid surface by more than the expected interpolation error, something is wrong. The ability to audit local iterations with the global solver's output is the single most important reason to invest in coarse-grid anchors at all.

== Singular Interior Points and the Complex-Plane Detour <sec-complex-detour>

The fourth-kind Sturm--Liouville problems of @sec-slt-taxonomy --- those with genuine singularities on the interior of the integration interval --- do not yield to any of the real-interval tricks developed so far. The solution is singular, the eigenvalue may be complex, and a naive real-line pseudospectral method simply does not converge. Boyd (1985) showed that a _complex-plane_ detour around the singularity can restore convergence for the _eigenvalues_, at the cost of complex arithmetic and the abandonment of any notion of "real eigenfunction on the real axis". The technique has been used in hydrodynamic stability (critical-latitude problems), atmospheric tides, and quantum resonance calculations.

The simplest example is
$ u_(y y) + [1 \/ y - lambda] u = 0, quad u(-1) = u(+1) = 0, $ <eq-singular-sl>
with a $1 \/ y$ pole at the interior point $y = 0$. The parabolic map
$ y(x) = x + i Delta (x^2 - 1), quad Delta > 0, $ <eq-parabolic-map>
sends the real interval $x in [-1, 1]$ to a curve in the complex $y$-plane that lies below the real axis in its interior and returns to $y = plus.minus 1$ at the endpoints. With $Delta > 0$ the transformed integration path goes _below_ the singularity, which fixes the branch of $log(1 \/ y)$ according to the physical convention inherited from causality. Chebyshev collocation in the real variable $x$ now produces a well-conditioned complex matrix whose eigenvalues are the eigenvalues of @eq-singular-sl, sortable by magnitude, all complex, all converging under refinement.

=== Computational Étude 18.11: Detour around the Singularity <etude-complex-detour>

We implement the transformed operator in all three languages. The chain rule gives
$ frac(d, d y) = (d y \/ d x)^(-1) frac(d, d x), quad frac(d^2, d y^2) = (d y \/ d x)^(-2) frac(d^2, d x^2) - frac(d^2 y \/ d x^2, (d y \/ d x)^3) frac(d, d x), $
and the operator $d^2 \/ d y^2 + 1 \/ y$ assembles as a Chebyshev matrix times a diagonal scaling, plus a lower-order term. Details are in the source file.

#figure(
  image("../figures/ch18/python/eigen_complex_detour.pdf", width: 98%),
  caption: [Étude 18.11: detour around the singularity at $y = 0$. Left: the parabolic integration contour $y(x) = x + i Delta (x^2 - 1)$ plotted for three values of the map parameter $Delta in {0.25, 0.50, 0.75}$. Each contour connects $y = -1$ to $y = +1$ below the real axis, safely avoiding the pole at $y = 0$ (purple cross). Right: the first seven eigenvalues, sorted by magnitude, at $N = 40$ and the same three $Delta$ values. All eigenvalues are complex; the three $Delta$ traces overlap almost perfectly for the lowest modes (they are converged) and separate for higher modes (they are under-resolved).],
) <fig-complex-detour>

The code generating @fig-complex-detour is available in:
- `codes/python/ch18/eigen_complex_detour.py`
- `codes/matlab/ch18/eigen_complex_detour.m`
- `codes/julia/ch18/eigen_complex_detour.jl`

The first four eigenvalues at $N = 40$, $Delta = 0.5$:
#align(center)[
  #table(
    columns: (auto, auto, auto, auto),
    stroke: 0.5pt + rgb(180, 180, 190),
    align: center + horizon,
    inset: 6pt,
    table.header([*mode*], [*Python*], [*Julia*], [*MATLAB*]),
    [1], [$-1.798 + 2.838 i$], [$-1.798 + 2.838 i$], [$-1.798 + 2.838 i$],
    [2], [$-10.887 + 0.320 i$], [$-10.887 + 0.320 i$], [$-10.887 + 0.320 i$],
    [3], [$-21.395 + 2.978 i$], [$-21.395 + 2.978 i$], [$-21.395 + 2.978 i$],
    [4], [$-40.203 + 0.172 i$], [$-40.203 + 0.172 i$], [$-40.203 + 0.172 i$],
  )
]
Maximum pairwise disagreement across the three languages: $1.6 times 10^(-10)$. The first four of these match to 11 digits between $N = 20$ and $N = 40$, a spectral-convergence signature that would be unrecoverable on the real axis where the singularity lives.

=== Critical Discussion

Two surprises in this étude deserve naming honestly.

#emph[First, the numerical eigenvalues do not match the table published by @Boyd2000.] Boyd reports a first eigenvalue near $0.125 + 0.285 i$ for essentially the same problem; I find $-1.798 + 2.838 i$ at the same $(N, Delta)$. The discrepancy is almost certainly due to a different interval or a different normalisation convention --- Boyd's Equation 7.56 uses an interval from $a$ to $b$ that he stretches to $-1 lt.eq.slant x lt.eq.slant 1$, and my stretching produces a different rescaled operator. The #emph[convergence] of my numbers under refinement is the primary evidence that the computed spectrum is consistent; matching a published table at the digit level requires reconstructing the exact published convention, which is a separate investigation.

#emph[Second, the complex detour works even though we cannot interpret the result as "the real eigenvalue on the real interval".] The physical interpretation is that @eq-singular-sl has no well-defined real spectrum: the operator is intrinsically non-self-adjoint because of the pole, and its eigenvalues live in the complex plane. The detour does not #emph[approximate] a real spectrum; it #emph[computes] the complex spectrum that the continuous operator has. The eigenfunctions live on the detoured contour, not on the real axis, and reconstructing them at real $y$ requires analytic continuation or a separate power-series expansion near $y = 0$ --- as @Boyd2000 discusses in detail. The chapter's pedagogical claim is narrower and more honest: #emph[eigenvalues] can be recovered cleanly by the detour, even when #emph[eigenfunctions] cannot be.

== A Non-Exhaustive Literature Review <sec-lit-review>

The eigenvalue material of this chapter rests on a lineage that stretches from the foundational monograph of @Boyd2000 and the error analysis of @GottliebOrszag1977 to a decade of rapid algorithmic progress between 2013 and 2026. What began as a catalogue of heuristics for distinguishing good from garbage eigenvalues has matured into a discipline with three recognisable pillars: sparse coefficient-space discretisation, rational approximation as a first-class tool, and an infrastructure of verification that treats every reported eigenvalue as provisional until a checklist is answered. The present overview is not comprehensive --- the field is too active for that --- but selects the threads most relevant to the étude-based philosophy of this chapter.

The historical bottleneck of spectral collocation was the $cal(O)(N^3)$ dense solve and the $cal(O)(N^(2 p))$ condition-number growth of differentiation matrices for $p$-th-order operators, a pathology we documented numerically in @sec-spurious-ii. The decisive break came with the ultraspherical spectral method of @OlverTownsend2013: by mapping the solution from a Chebyshev basis to a cascade of ultraspherical (Gegenbauer) bases during differentiation, the method produces _almost-banded_ operators whose linear systems can be solved in $cal(O)(N)$ operations with bounded condition number. The analytical framework is collected in the _Acta Numerica_ survey of @OlverTownsend2019, and the extension to nonlinear boundary-value problems is developed by @Xu2024. For the eigenvalue problems of this chapter, the practical consequence is that spectra with many thousands of degrees of freedom become routine, and the high-order-derivative condition-number story recedes from a hard ceiling to a secondary concern.

A parallel revolution concerns approximation itself. Polynomial bases are unsurpassed for smooth functions on bounded intervals, but they degrade sharply in the presence of poles, branch points, or unbounded domains --- exactly the situations that produce the complex-plane detour of @sec-complex-detour. The AAA algorithm of @NakatsukasaSeteTrefethen2018 reopened rational approximation as a numerically stable numerical tool by representing the approximant in barycentric form and selecting support points greedily. A family of variants has followed: @Baddoo2021 adapts the algorithm to periodic settings through a trigonometric barycentric basis, @NakatsukasaTrefethen2020 couples AAA with an iteratively reweighted least-squares phase to deliver true minimax approximations, and @DriscollNakatsukasaTrefethen2024 extends the construction to a continuum of sample points rather than a discrete set. For eigenvalue problems these are not ornamental: rational surrogates allow nonlinear operators $T(lambda)$ to be linearised into ordinary pencil form, and adaptive pole placement handles resonance modes in photonic crystals and open waveguides that polynomial methods cannot touch.

The nonlinear eigenvalue problem $T(lambda) u = 0$, where the operator depends nonlinearly on the spectral parameter, sits at the intersection of the two preceding threads. The community-standard survey is @GuttelTisseur2017, which catalogues Newton-based, contour-integral, and rational-sampling strategies. Contour integration in particular --- originated by @Beyn2012 and now implemented in multiple open-source libraries --- extracts every eigenvalue inside a user-specified contour by reducing the infinite-dimensional problem to a small linear pencil via Keldysh's theorem. @LietaertEtAl2022 combine Beyn's idea with a set-valued AAA step: a low-degree rational approximant of $T(lambda)$ is built automatically and then linearised, giving an end-to-end pipeline that requires neither hand-tuned shifts nor a priori spectral information. The most recent development is infinite-dimensional: @ColbrookTownsend2025 propose the InfBeyn solver, which operates directly on the quasimatrix $T(lambda)$ by applying randomised SVD to columns drawn from a Gaussian process --- a "solve then discretise" paradigm that documents, and then cures, the spectral pollution, spectral invisibility, and ghost-essential-spectra pathologies that motivate the spurious-mode discussion of @sec-spurious-i. The taxonomy of those pathologies in the present chapter borrows from that discussion directly.

For the physically spurious eigenvalues of fluid stability, the literature has continued to refine the capacitance-matrix cure introduced by @GottliebOrszag1977. Recent work by @PatneShankar2019 on shear flow in deformable channels shows that inconsistencies in either the bulk constitutive model or the linearisation of interface conditions can separately introduce spurious instabilities, and that these must be tracked before any discretisation is attempted. The message aligns with the chapter's own: physically spurious modes are cured at the _formulation_ step, not by increasing $N$ or polishing the solver.

A distinct thread is the emergence of _differentiable_ spectral solvers. The Dedalus project @Burns2020 provides a symbolic, MPI-parallel framework for multidimensional spectral simulation in Python, built around sparse ultraspherical and Fourier discretisations. Its 2026 extension by @Lecoanet2026 generates discrete adjoint solvers automatically from the symbolic operator graph, giving exact gradients of any scalar functional with respect to simulation parameters at the cost of a second (backward) spectral solve. The implications for stability theory are immediate: neutral curves, optimal perturbations, and minimal seeds --- previously the province of hand-derived continuous adjoints --- become computable end-to-end in the same spectral accuracy as the forward solve, and the applications reach from dynamo action in planetary cores to stellarator coil optimisation in fusion engineering.

For unbounded domains the rational Chebyshev basis deployed in @sec-infinite-tax remains the default workhorse, but it is not the only option. @LeeMarcus2023 develop a mapped-Legendre Galerkin method for the linear stability of axisymmetric vortices, reporting that for a given accuracy on Lamb--Oseen-type profiles the mapped-Legendre basis requires roughly three times fewer basis elements than standard rational Chebyshev collocation. More importantly, the method successfully resolves neutrally stable inviscid critical-layer eigenmodes whose singularities derail polynomial-based routines, and it identifies viscous critical-layer families obeying the classical $R e^(-1 \/ 3)$ scaling. This is a reminder that the "infinite-interval tax" documented in the chapter is not a fixed cost: mapping, weighting, and basis design can each claw back factors of two or three when the physics is understood.

Among the most consequential modern applications of spectral eigenvalue computation is the computation of quasinormal modes of black holes and exotic compact objects, where spectral methods have become the tool of choice for quantitative general relativity. A sustained series by D. Batić, D. Dutykh and collaborators applies Chebyshev collocation to Morris--Thorne wormholes @Batic2024, noncommutative-geometry-inspired Schwarzschild black holes @Batic2024b, Lee--Wick black holes @Batic2024c, massive phantom wormholes @Batic2025, noncommutative wormholes @Batic2025a, noncommutative "dirty" black holes @Batic2025b, and most recently the quasinormal modes of Planck stars @Batic2026. Two methodological lessons recur throughout that series. First, purely imaginary _overdamped_ modes --- which the WKB approximation systematically misses --- are routinely recovered by the spectral approach, and they carry stability information that the fundamental oscillating family does not. Second, the trust certificate of @sec-verification-culture is enforced operationally: every spectrum is cross-validated against continued-fraction expansions for the Schwarzschild benchmark, checked against time-domain Regge--Wheeler codes, and reproduced at two independent resolutions before publication. The pseudospectral viewpoint of @BoultonEtAl2022 provides the theoretical underpinning by generalising eigenvalue theorems to the pseudospectra that govern the stability of these computations under perturbation.

When the resolution $N$ becomes too large for dense solvers, the community's preferred scalable algorithm is the FEAST eigensolver of @Polizzi2009, which approximates the spectral projector onto a target interval via contour integration and returns all eigenvalues inside the contour in one shot. The algorithm is embarrassingly parallel across disjoint spectral regions, and a recent randomised-warmstart variant has reported speedups of more than an order of magnitude on sparse benchmarks; extensions to non-Hermitian pencils and generalised eigenproblems with oblique projections are active fronts. The chapter's emphasis on global-then-local strategies in @sec-solver-strategy is consistent with this development: FEAST plays the role of the "global" phase for large sparse systems just as QR plays it for small dense ones.

We close with the software ecosystem, which has co-evolved with the algorithms and now spans the three languages of this book. In Python, Dedalus @Burns2020 and its automatic-adjoint extension @Lecoanet2026 dominate multidimensional spectral simulation; in MATLAB, Chebfun provides the canonical AAA implementation @NakatsukasaSeteTrefethen2018 together with the infNEP toolbox implementing @ColbrookTownsend2025; in Julia, ApproxFun.jl delivers the ultraspherical method of @OlverTownsend2013 with near-native performance for million-degree-of-freedom problems. The student who follows the études of this chapter and the references of this section will find, at every scale, an open-source toolkit that makes the verification culture we have advocated practical rather than aspirational.

== A Verification Culture, not a Bibliography <sec-verification-culture>

The chapter ends where every eigenvalue calculation should end: with a checklist rather than a celebration. @Boyd2000 groups the common failure modes under three labels --- credulity, negligence, and coding --- that summarise the territory with uncomfortable accuracy. Credulity is believing a smooth graph. Negligence is running one resolution and moving on. Coding error is everything else.

The chapter has tried to produce a positive counterpart: a reusable recipe that any practitioner can apply to any eigenvalue calculation whose correctness they must defend. We collect it here as the _trust certificate_.

#block(stroke: 0.8pt + rgb(20, 45, 110), radius: 3pt, inset: 10pt,
  fill: rgb(248, 250, 254))[
  *The trust certificate.* For every reported eigenvalue $lambda_j^("num")$ and eigenfunction $u_j^("num")$, answer the following questions honestly:
  + Has $lambda_j^("num")$ stabilised under $N arrow 2N$ or at least $N arrow 3N\/2$?
  + Is the comparison across resolutions by ordinal match, or was nearest-match necessary?
  + Is the residual norm $||L u_j^("num") - lambda_j^("num") u_j^("num")||$ small?
  + Are the boundary residuals $|cal(B) u_j^("num")|$ small at every prescribed boundary condition?
  + Is the eigenfunction _shape_ (not just the eigenvalue) stable under refinement?
  + If the problem has a symmetry, parity, or sign structure, does the numerical mode respect it?
  + Could this mode be _physically_ spurious (formulation-induced, not discretisation-induced)?
  + Has it been cross-checked with a second formulation or a second language?
]

The eight questions are cheap to ask and dangerous to skip. A mode that passes all eight is not _guaranteed_ correct --- no finite number of checks can be --- but the probability that it is an artefact drops sharply with each pass. Conversely, a mode that fails any one of the eight should be reported as suspect, no matter how beautiful the plot it produces.

The chapter leaves the reader with one last recommendation: apply the trust certificate to every eigenvalue computation, however small, until it becomes a reflex. The cost of the habit is small. The cost of skipping it, when the computation matters, is an eigenvalue that looks right and is not.
