// textbook/chapters/classical_pdes.typ
// Chapter 2: Classical Second Order PDEs and Separation of Variables
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Email: denys.dutykh@ku.ac.ae
// Homepage: https://www.denys-dutykh.com/
// Last modified: February 2026
#import "../styles/template.typ": dropcap

= Classical Second Order PDEs and Separation of Variables

#dropcap[In this opening chapter we derive exact solutions for three classical linear partial differential equations: the _heat equation_ (parabolic), the _wave equation_ (hyperbolic), and the _Laplace equation_ (elliptic). These solutions are found by the _method of separation of variables_, which expresses the solution as an infinite series of eigenfunctions.]

Why begin a book on _numerical_ methods with _analytical_ solutions? Because separation of variables is the theoretical ancestor of spectral methods @Orszag1971 @GottliebOrszag1977 @Trefethen2000. When we later truncate these infinite series at some finite $N$ and compute with only the first $N$ modes, we are doing exactly what a spectral solver does, but with pen and paper first. This chapter thus serves as the conceptual bridge between classical analysis and modern computation. The methods presented here are classical and thoroughly developed in standard references such as @Tikhonov1963, @Boyd2000, and @Canuto2006.

We treat three model problems:

- heat equation with periodic boundary conditions in one spatial dimension,
- wave equation on a bounded interval,
- Laplace equation in a simple domain.

We begin with a complete analytic solution of the heat equation. The other two examples will follow the same pattern.

== Heat Equation with Periodic Boundary Conditions

We consider the one dimensional heat equation on the interval $[0, 2 pi]$ with periodic boundary conditions. The unknown $u(x,t)$ represents, for example, the temperature at point $x$ and time $t$.

The problem is
$ frac(partial u, partial t) (x,t) = frac(partial^2 u, partial x^2) (x,t), quad x in [0, 2 pi], space t > 0, $
with periodic boundary conditions
$ u(x + 2 pi, t) = u(x,t) $
for all real $x$ and all $t > 0$, and initial condition
$ u(x,0) = f(x), quad x in [0, 2 pi]. $

We assume that $f$ is smooth and $2 pi$ periodic:
$ f(x + 2 pi) = f(x). $

Our goal is to obtain an explicit representation of $u(x,t)$ as an infinite series and to see how separation of variables leads naturally to a Fourier series in space. The technique we employ here follows the classical approach presented in @Tikhonov1963.

=== Step 1: Separation Ansatz

We look for nontrivial solutions of the form
$ u(x,t) = X(x) dot T(t), $
where $X$ depends only on $x$ and $T$ depends only on $t$.

Substituting into the PDE gives
$ X(x) dot T'(t) = X''(x) dot T(t). $

We assume $X$ and $T$ are not identically zero, so we can divide both sides by $X(x) dot T(t)$:
$ frac(T'(t), T(t)) = frac(X''(x), X(x)). $

The left side depends only on $t$, the right side only on $x$. Therefore both sides must be equal to the same constant, which we denote by $-lambda$:
$ frac(T'(t), T(t)) = frac(X''(x), X(x)) = -lambda. $

We obtain two ordinary differential equations:
$ T'(t) + lambda T(t) = 0, $
$ X''(x) + lambda X(x) = 0. $

The periodic boundary conditions for $u$ imply periodic conditions for $X$:
$ X(0) = X(2 pi), quad X'(0) = X'(2 pi). $

We have arrived at a spatial eigenvalue problem for $X$. The systematic treatment of such eigenvalue problems is a cornerstone of the theory of partial differential equations, with roots in Fourier's original treatise @Fourier1822; see @Tikhonov1963 for a comprehensive exposition.

=== Step 2: Spatial Eigenvalue Problem with Periodic Boundary Conditions

We now solve
$ X''(x) + lambda X(x) = 0, $
with
$ X(0) = X(2 pi), quad X'(0) = X'(2 pi). $

We consider three cases: $lambda < 0$, $lambda = 0$, and $lambda > 0$.

*Case 1: $lambda < 0$*

Write $lambda = -mu^2$ with $mu > 0$. The equation becomes
$ X''(x) - mu^2 X(x) = 0. $

The general solution is
$ X(x) = A e^(mu x) + B e^(-mu x). $

Imposing periodicity $X(0) = X(2 pi)$ gives
$ A + B = A e^(2 mu pi) + B e^(-2 mu pi). $

Imposing $X'(0) = X'(2 pi)$ gives
$ mu(A - B) = mu(A e^(2 mu pi) - B e^(-2 mu pi)). $

Since $mu > 0$, we can divide by $mu$ and rewrite both conditions as a homogeneous linear system:
$ A(1 - e^(2 mu pi)) + B(1 - e^(-2 mu pi)) = 0, $
$ A(1 - e^(2 mu pi)) - B(1 - e^(-2 mu pi)) = 0. $

Adding these two equations gives
$ 2 A (1 - e^(2 mu pi)) = 0. $
Since $mu > 0$, we have $e^(2 mu pi) > 1$, so $1 - e^(2 mu pi) eq.not 0$. Therefore $A = 0$.

Subtracting the second equation from the first gives
$ 2 B (1 - e^(-2 mu pi)) = 0. $
Since $mu > 0$, we have $e^(-2 mu pi) < 1$, so $1 - e^(-2 mu pi) eq.not 0$. Therefore $B = 0$.

Hence the only solution is $A = B = 0$, the trivial solution. There are no nontrivial periodic eigenfunctions for $lambda < 0$, and we discard this case.

*Case 2: $lambda = 0$*

The equation reduces to
$ X''(x) = 0. $

Its general solution is
$ X(x) = A + B x. $

Periodicity $X(0) = X(2 pi)$ gives
$ A = A + 2 pi B $
so $B = 0$. Then $X(x) = A$ is constant.

The derivative is $X'(x) = 0$, so $X'(0) = X'(2 pi)$ is automatically satisfied.

Thus $lambda = 0$ gives one eigenfunction
$ X_0(x) = 1 $
(up to a multiplicative constant).

*Case 3: $lambda > 0$*

Write $lambda = k^2$ with $k > 0$. The equation becomes
$ X''(x) + k^2 X(x) = 0. $

The general solution is
$ X(x) = A cos(k x) + B sin(k x). $

We now impose periodicity. First,
$ X(0) = A, quad X(2 pi) = A cos(2 pi k) + B sin(2 pi k). $
The condition $X(0) = X(2 pi)$ gives
$ A = A cos(2 pi k) + B sin(2 pi k). $

Next,
$ X'(x) = -A k sin(k x) + B k cos(k x), $
so
$ X'(0) = B k, quad X'(2 pi) = -A k sin(2 pi k) + B k cos(2 pi k). $
The condition $X'(0) = X'(2 pi)$ gives
$ B k = -A k sin(2 pi k) + B k cos(2 pi k). $

We can divide by $k$ the last identity (for $k != 0$) and write the system as
$ A(1 - cos(2 pi k)) - B sin(2 pi k) = 0, $
$ A sin(2 pi k) + B(1 - cos(2 pi k)) = 0. $

For a nontrivial pair $(A,B)$ the determinant of the system must vanish:
$ (1 - cos(2 pi k))^2 + (sin(2 pi k))^2 = 0. $

The left side is a sum of squares, so it is zero if and only if
$ 1 - cos(2 pi k) = 0, quad sin(2 pi k) = 0. $

Hence
$ cos(2 pi k) = 1, quad sin(2 pi k) = 0. $

This happens exactly when $k$ is an integer:
$ k = n, quad n in ZZ. $

The case $n = 0$ corresponds to $lambda = 0$, which we have already treated. For $n gt.eq.slant 1$ we obtain eigenvalues
$ lambda_n = n^2, quad n = 1, 2, 3, dots $

For each $n gt.eq.slant 1$ the corresponding eigenfunctions can be chosen as
$ X_n^((c))(x) = cos(n x), quad X_n^((s))(x) = sin(n x). $

These functions are $2 pi$ periodic, and their derivatives are also $2 pi$ periodic, so the boundary conditions are satisfied.

We have therefore found a complete set of spatial eigenfunctions for the heat equation with periodic boundary conditions:

- a constant mode $X_0(x) = 1$ (eigenvalue $lambda_0 = 0$),
- cosine modes $X_n^((c))(x) = cos(n x)$,
- sine modes $X_n^((s))(x) = sin(n x)$,

with eigenvalues $lambda_n = n^2$ for $n gt.eq.slant 1$.

=== Step 3: Time Dependent Factors

For each eigenvalue $lambda$ the corresponding time factor satisfies
$ T'(t) + lambda T(t) = 0. $

The solution is
$ T(t) = C e^(-lambda t). $

For the constant mode $lambda_0 = 0$ we obtain
$ T_0(t) = C_0 $
(constant in time).

For $n gt.eq.slant 1$ we have
$ T_(n) (t) = C_n e^(-n^2 t). $

Combining space and time, we form products $u(x,t) = X(x) dot T(t)$ to obtain separated solutions. For each mode, the arbitrary constant from the time factor can be absorbed into a single new constant. Specifically, we write
$ u_0(x,t) = A_0, $
$ u_n^((c))(x,t) = A_n cos(n x) dot e^(-n^2 t), $
$ u_n^((s))(x,t) = B_n sin(n x) dot e^(-n^2 t), $
where we have renamed the constants: $A_0 = C_0$, and for $n gt.eq.slant 1$, the constant $A_n$ absorbs the coefficient from the cosine mode while $B_n$ absorbs the coefficient from the sine mode. Note that each separated solution carries its own independent arbitrary constant.

Since the heat equation is linear, any linear combination of these separated solutions is again a solution. Therefore the general solution that satisfies the periodic boundary conditions can be written as an infinite series
$ u(x,t) = A_0 + sum_(n=1)^infinity (A_n cos(n x) + B_n sin(n x)) e^(-n^2 t). $

The coefficients $A_0$, $A_n$, $B_n$ remain to be determined from the initial condition.

=== Step 4: Imposing the Initial Condition and Fourier Series

We impose the initial condition
$ u(x,0) = f(x). $

Setting $t = 0$ in the general solution we obtain
$ u(x,0) = A_0 + sum_(n=1)^infinity (A_n cos(n x) + B_n sin(n x)) = f(x). $

This is exactly the Fourier series expansion of the $2 pi$ periodic function $f$. Under mild regularity assumptions, $f$ has a Fourier series
$ f(x) = a_0 + sum_(n=1)^infinity (a_n cos(n x) + b_n sin(n x)) $
with Fourier coefficients
$ a_0 = frac(1, 2 pi) integral_0^(2 pi) f(x) dif x, $
$ a_n = frac(1, pi) integral_0^(2 pi) f(x) cos(n x) dif x, quad n gt.eq.slant 1, $
$ b_n = frac(1, pi) integral_0^(2 pi) f(x) sin(n x) dif x, quad n gt.eq.slant 1. $

By uniqueness of the Fourier expansion, we must have
$ A_0 = a_0, quad A_n = a_n, quad B_n = b_n. $

Thus the coefficients in the heat equation solution are exactly the Fourier coefficients of the initial data.

=== Step 5: Final Explicit Solution and Interpretation

Substituting these coefficients into the general solution we obtain the explicit formula
$ u(x,t) = a_0 + sum_(n=1)^infinity (a_n cos(n x) + b_n sin(n x)) e^(-n^2 t), quad t > 0. $

This infinite sum solves the heat equation with periodic boundary conditions and initial data $f$. Each Fourier mode decays exponentially in time at a rate proportional to its eigenvalue $n^2$. High frequency modes (large $n$) decay faster, which expresses the smoothing effect of the heat equation.

The constant term $a_0$ does not decay. It represents the average value of $f$ on $[0, 2 pi]$, which is preserved by the heat flow.

From a spectral viewpoint, the functions
$ 1, cos(x), sin(x), cos(2 x), sin(2 x), dots $
form an eigenbasis of the spatial operator
$ L u = frac(partial^2 u, partial x^2) $
with periodic boundary conditions. The evolution of each eigenmode is independent and given simply by multiplication by $e^(-n^2 t)$ in time.

Later, in the numerical part of this book, we will approximate $u(x,t)$ by truncating the sum to finitely many modes:
$ u_N(x,t) = a_0 + sum_(n=1)^N (a_n cos(n x) + b_n sin(n x)) e^(-n^2 t). $

This truncation is the essence of a Fourier spectral method; recent extensions treat both space and time spectrally @Kaur2025. The analytic solution derived here is the infinite dimensional limit of that numerical procedure.

== Computational Étude 2.1: Heat Equation on a Ring

To visualize the smoothing effect of heat diffusion, we compute the truncated Fourier series solution for a triangle wave initial condition:
$ f(x) = pi - |x - pi|, quad x in [0, 2 pi]. $

This function is continuous but has a corner (non-differentiable point) at $x = pi$. Its Fourier series contains only cosine terms with coefficients decaying as $1 \/ n^2$:
$ f(x) = frac(pi, 2) + frac(4, pi) sum_(k=1)^infinity frac(cos((2k-1) x), (2k-1)^2). $

The key portion of the implementation evaluates the truncated Fourier series at any point in space and time. In Python:

```python
def heat_solution(x, t, a0, a_n, b_n):
    u = np.full_like(x, a0, dtype=float)
    n_modes = len(a_n) - 1
    for n in range(1, n_modes + 1):
        decay = np.exp(-n**2 * t)
        u += (a_n[n] * np.cos(n * x) + b_n[n] * np.sin(n * x)) * decay
    return u
```

The equivalent MATLAB implementation:

```matlab
u = a0 * ones(size(x));
for n = 1:N_MODES
    decay = exp(-n^2 * t);
    u = u + (a_n(n+1) * cos(n*x) + b_n(n+1) * sin(n*x)) * decay;
end
```

The Julia implementation:

```julia
function heat_solution(x, t, a0, a_n, b_n)
    u = fill(a0, length(x))
    n_modes = length(a_n) - 1
    for n in 1:n_modes
        decay = exp(-n^2 * t)
        @. u += (a_n[n+1] * cos(n * x) + b_n[n+1] * sin(n * x)) * decay
    end
    return u
end
```

@fig-heat-evolution shows the evolution of $u_N (x, t)$ with $N = 50$ modes at several time values. At $t = 0$ the triangle wave is faithfully reproduced. As time increases, the higher frequency modes decay exponentially faster than the lower ones (the $n$-th mode decays as $e^(-n^2 t)$), and the solution rapidly smooths toward the constant equilibrium $u = pi \/ 2$.

#figure(
  image("../figures/ch02/python/heat_evolution.pdf", width: 95%),
  caption: [Evolution of the heat equation solution with a triangle wave initial condition. The higher frequency modes decay rapidly, smoothing the initial corner at $x = pi$.],
) <fig-heat-evolution>

The code generating @fig-heat-evolution is available in:
- `codes/python/ch02/heat_equation_evolution.py`
- `codes/matlab/ch02/heat_equation_evolution.m`
- `codes/julia/ch02/heat_equation_evolution.jl`

A complementary view of the solution is provided by the waterfall plot in @fig-heat-waterfall, which displays the entire space-time evolution as a three-dimensional surface. The smoothing effect of the heat equation is clearly visible: the initial sharp triangle wave rapidly flattens as time progresses, with the solution approaching the constant equilibrium state $u = pi\/2$.

#figure(
  image("../figures/ch02/python/heat_waterfall.pdf", width: 95%),
  caption: [Waterfall plot showing the complete space-time evolution of the heat equation solution. The initial triangle wave smooths rapidly as higher frequency modes decay exponentially.],
) <fig-heat-waterfall>

The code generating @fig-heat-waterfall is available in:
- `codes/python/ch02/heat_equation_waterfall.py`
- `codes/matlab/ch02/heat_equation_waterfall.m`
- `codes/julia/ch02/heat_equation_waterfall.jl`

=== Discussion

@fig-heat-evolution and @fig-heat-waterfall vividly illustrate the smoothing property of the heat equation: although the triangle-wave initial condition has a corner at $x = pi$, the solution becomes infinitely smooth for any $t > 0$. The mechanism is the exponential decay factor $e^(-n^2 t)$, which suppresses high-frequency modes far more aggressively than low-frequency ones. With only $N = 50$ modes the truncated series is visually indistinguishable from the exact solution at all displayed times, foreshadowing the rapid convergence that spectral methods will exploit.

The long-time limit $u arrow.r a_0 = pi\/2$ (the spatial mean of the initial data) is clearly visible in the waterfall plot: all oscillatory structure disappears, and only the zero-frequency mode survives. This étude thus demonstrates both the physics of diffusive equilibration and the mathematics of Fourier coefficient decay that will underpin the convergence theory in later chapters.

== Wave Equation with Dirichlet Boundary Conditions

We now consider the one dimensional wave equation on a finite interval with homogeneous Dirichlet boundary conditions. This is the classical model for a vibrating string of length $L$ with both ends fixed, first studied by d'Alembert @dAlembert1747 and Euler @Euler1748.

Let $u(x,t)$ denote the vertical displacement of the string at position $x in [0,L]$ and time $t > 0$. The equation of motion is
$ frac(partial^2 u, partial t^2) (x,t) = c^2 frac(partial^2 u, partial x^2) (x,t), quad 0 < x < L, space t > 0, $
where $c > 0$ is the wave speed.

The boundary conditions express that the endpoints of the string are clamped:
$ u(0,t) = 0, quad u(L,t) = 0, quad t > 0. $

We prescribe the initial displacement and initial velocity:
$ u(x,0) = f(x), quad frac(partial u, partial t) (x,0) = g(x), quad 0 < x < L, $
with suitable functions $f$ and $g$ that vanish at $x = 0$ and $x = L$.

As in the heat equation example, we use separation of variables and obtain a representation of the solution as an infinite series in spatial eigenfunctions. This time the temporal factors are oscillatory instead of decaying. The mathematical theory of the vibrating string is a classical topic covered extensively in @Tikhonov1963.

=== Step 1: Separation Ansatz

We look for nontrivial solutions of the form
$ u(x,t) = X(x) dot T(t). $

Substituting into the wave equation gives
$ X(x) dot T''(t) = c^2 X''(x) dot T(t). $

Assuming $X$ and $T$ are not identically zero, we divide both sides by $c^2 X(x) dot T(t)$:
$ frac(T''(t), c^2 T(t)) = frac(X''(x), X(x)). $

The left side depends only on $t$, the right side only on $x$. Therefore both sides must be equal to a constant, which we denote by $-lambda$:
$ frac(T''(t), c^2 T(t)) = frac(X''(x), X(x)) = -lambda. $

We obtain the pair of ordinary differential equations
$ T''(t) + c^2 lambda T(t) = 0, $
$ X''(x) + lambda X(x) = 0, $
with boundary conditions
$ X(0) = 0, quad X(L) = 0. $

As in the heat equation case, we have a spatial eigenvalue problem for $X$.

=== Step 2: Spatial Eigenvalue Problem with Dirichlet Boundary Conditions

We must solve
$ X''(x) + lambda X(x) = 0, quad 0 < x < L, $
with
$ X(0) = 0, quad X(L) = 0. $

We again consider three cases: $lambda < 0$, $lambda = 0$, and $lambda > 0$.

*Case 1: $lambda < 0$*

Write $lambda = -mu^2$ with $mu > 0$. The equation becomes
$ X''(x) - mu^2 X(x) = 0. $

The general solution is
$ X(x) = A e^(mu x) + B e^(-mu x). $

The boundary condition at $x = 0$ gives
$ X(0) = A + B = 0 quad arrow.r.double quad B = -A. $

Then
$ X(L) = A e^(mu L) - A e^(-mu L) = A (e^(mu L) - e^(-mu L)). $

The condition $X(L) = 0$ implies
$ A (e^(mu L) - e^(-mu L)) = 0. $

Since $e^(mu L) eq.not e^(-mu L)$ for $mu > 0$, we must have $A = 0$. Then $B = 0$ and the solution is trivial. Therefore there are no nontrivial eigenfunctions for $lambda < 0$.

*Case 2: $lambda = 0$*

The equation reduces to
$ X''(x) = 0, $
whose general solution is
$ X(x) = A + B x. $

The boundary conditions give
$ X(0) = A = 0, $
$ X(L) = A + B L = B L = 0. $

Hence $B = 0$ and $X$ is again trivial. There is no nontrivial eigenfunction for $lambda = 0$.

*Case 3: $lambda > 0$*

Write $lambda = k^2$ with $k > 0$. The equation becomes
$ X''(x) + k^2 X(x) = 0. $

The general solution is
$ X(x) = A cos(k x) + B sin(k x). $

The boundary condition at $x = 0$ gives
$ X(0) = A = 0. $

So $X(x) = B sin(k x)$. The boundary condition at $x = L$ gives
$ X(L) = B sin(k L) = 0. $

For a nontrivial solution we need $B eq.not 0$, so we must have
$ sin(k L) = 0. $

Therefore
$ k L = n pi, quad n = 1, 2, 3, dots $

The corresponding values of $k$ are
$ k_n = frac(n pi, L), quad n = 1, 2, 3, dots $

We conclude that the eigenvalues and eigenfunctions are
$ lambda_n = k_n^2 = (frac(n pi, L))^2, $
$ X_n (x) = sin(frac(n pi x, L)), quad n = 1, 2, 3, dots $

Each $X_n$ vanishes at $x = 0$ and $x = L$, as required by the Dirichlet boundary conditions. The family ${X_n}_(n gt.eq.slant 1)$ is orthogonal in $L^2 (0,L)$:
$ integral_0^L sin(frac(n pi x, L)) sin(frac(m pi x, L)) dif x = cases(0 & "if" n eq.not m, L\/2 & "if" n = m.) $

These eigenfunctions will form the spatial basis in our series solution.

=== Step 3: Time Dependent Factors

For each eigenvalue $lambda_n$ the time factor $T_n$ satisfies
$ T_(n)^('')(t) + c^2 lambda_n T_n (t) = 0. $

Using $lambda_n = (n pi \/ L)^2$ we can write
$ T_(n)^('')(t) + omega_n^2 T_n (t) = 0, $
where
$ omega_n = c frac(n pi, L), quad n = 1, 2, 3, dots $

The general solution of this second order linear ODE is
$ T_n (t) = A_n cos(omega_n t) + B_n sin(omega_n t), $
where $A_n$ and $B_n$ are constants.

Combining the space and time factors, we obtain separated solutions
$ u_n (x,t) = (A_n cos(omega_n t) + B_n sin(omega_n t)) sin(frac(n pi x, L)), quad n = 1, 2, 3, dots $

Because the wave equation is linear, any linear combination of these separated solutions is again a solution. Therefore the general solution satisfying the Dirichlet boundary conditions can be written as an infinite series
$ u(x,t) = sum_(n=1)^infinity (a_n cos(omega_n t) + b_n sin(omega_n t)) sin(frac(n pi x, L)), $
for suitable coefficients $a_n$ and $b_n$.

These coefficients will be determined from the initial conditions.

=== Step 4: Imposing the Initial Conditions and Sine Series

We now use the initial displacement and velocity.

At $t = 0$ we have
$ u(x,0) = sum_(n=1)^infinity (a_n cos(0) + b_n sin(0)) sin(frac(n pi x, L)) = sum_(n=1)^infinity a_n sin(frac(n pi x, L)). $

The initial condition $u(x,0) = f(x)$ becomes
$ f(x) = sum_(n=1)^infinity a_n sin(frac(n pi x, L)). $

This is the Fourier sine series of $f$ on the interval $(0,L)$.

Similarly, we differentiate $u$ with respect to $t$:
$ frac(partial u, partial t) (x,t) = sum_(n=1)^infinity (-a_n omega_n sin(omega_n t) + b_n omega_n cos(omega_n t)) sin(frac(n pi x, L)). $

Evaluating at $t = 0$ gives
$ frac(partial u, partial t) (x,0) = sum_(n=1)^infinity b_n omega_n sin(frac(n pi x, L)). $

The initial condition $frac(partial u, partial t) (x,0) = g(x)$ becomes
$ g(x) = sum_(n=1)^infinity b_n omega_n sin(frac(n pi x, L)). $

So the sequence ${a_n}$ consists of the Fourier sine coefficients of $f$, and the sequence ${b_n omega_n}$ consists of the Fourier sine coefficients of $g$.

Using the orthogonality relations, we obtain explicit formulas for the coefficients. For $n gt.eq.slant 1$,
$ a_n = frac(2, L) integral_0^L f(x) sin(frac(n pi x, L)) dif x, $
and
$ b_n omega_n = frac(2, L) integral_0^L g(x) sin(frac(n pi x, L)) dif x. $

Therefore
$ b_n = frac(2, L omega_n) integral_0^L g(x) sin(frac(n pi x, L)) dif x = frac(2, L) frac(1, c n pi \/ L) integral_0^L g(x) sin(frac(n pi x, L)) dif x. $

In summary,
$ a_n = frac(2, L) integral_0^L f(x) sin(frac(n pi x, L)) dif x, $
$ b_n = frac(2, L omega_n) integral_0^L g(x) sin(frac(n pi x, L)) dif x, quad omega_n = c frac(n pi, L). $

=== Step 5: Final Explicit Solution and Interpretation

Substituting these coefficients into the series, we obtain the explicit solution of the wave equation with Dirichlet boundary conditions:
$ u(x,t) = sum_(n=1)^infinity [a_n cos(omega_n t) + b_n sin(omega_n t)] sin(frac(n pi x, L)), $
where $omega_n = c n pi \/ L$ and the Fourier sine coefficients are
$ a_n = frac(2, L) integral_0^L f(y) sin(frac(n pi y, L)) dif y, quad b_n = frac(2, n pi c) integral_0^L g(y) sin(frac(n pi y, L)) dif y. $

Each term in the sum is a normal mode of vibration: a standing wave with spatial shape $sin(n pi x \/ L)$ and temporal oscillation at frequency $omega_n$. The coefficients of $cos(omega_n t)$ and $sin(omega_n t)$ are determined by the initial displacement $f$ and initial velocity $g$ through their Fourier sine coefficients.

Comparing with the heat equation:

- For the heat equation, each mode decayed like $e^(-n^2 t)$ and the solution became smoother in time.
- For the wave equation, each mode oscillates periodically in time with constant amplitude, reflecting conservation of energy in the undamped string.

From a spectral viewpoint, the functions
$ sin(frac(pi x, L)), sin(frac(2 pi x, L)), sin(frac(3 pi x, L)), dots $
form an eigenbasis of the spatial operator
$ L u = frac(partial^2 u, partial x^2) $
with Dirichlet boundary conditions. In this basis, the evolution is diagonal: each mode evolves independently according to a simple harmonic oscillator in time.

As in the heat equation example, a spectral method will approximate $u(x,t)$ by truncating the infinite sum. For some integer $N gt.eq.slant 1$ we consider the finite approximation
$ u_N (x,t) = sum_(n=1)^N (a_n cos(omega_n t) + b_n sin(omega_n t)) sin(frac(n pi x, L)). $

The analytic series above is the infinite dimensional limit of this spectral representation.

== Computational Étude 2.2: Vibrating String

To visualize the oscillatory behavior of the vibrating string, we compute the truncated Fourier sine series solution for a plucked string initial condition. The string is plucked at its center, forming a triangular initial displacement:
$ f(x) = cases(frac(2h,L) x & "for" 0 lt.eq.slant x lt.eq.slant L\/2, 2h (1 - x\/L) & "for" L\/2 lt.eq.slant x lt.eq.slant L) $
with zero initial velocity $g(x) = 0$. Here $h$ denotes the height of the pluck at the center.

The Fourier sine coefficients of this triangular shape are
$ a_n = frac(8h, n^2 pi^2) sin(frac(n pi, 2)), $
which gives nonzero values only for odd $n$, with alternating signs. The slow algebraic decay of these coefficients reflects the limited regularity of the triangular initial shape; see @Boyd2000 for a thorough discussion of convergence rates for non-smooth functions.

The key portion of the implementation computes the solution at any point in space and time. In Python:

```python
def wave_solution(x, t, a_n, b_n, L, c):
    u = np.zeros_like(x, dtype=float)
    n_modes = len(a_n) - 1
    for n in range(1, n_modes + 1):
        omega_n = c * n * np.pi / L
        spatial = np.sin(n * np.pi * x / L)
        temporal = a_n[n] * np.cos(omega_n * t) + b_n[n] * np.sin(omega_n * t)
        u += temporal * spatial
    return u
```

The equivalent MATLAB implementation:

```matlab
u = zeros(size(x));
for n = 1:N_MODES
    omega_n = C * n * pi / L;
    spatial = sin(n * pi * x / L);
    temporal = a_n(n+1) * cos(omega_n * t) + b_n(n+1) * sin(omega_n * t);
    u = u + temporal * spatial;
end
```

The Julia implementation:

```julia
function wave_solution(x, t, a_n, b_n, L, c)
    u = zeros(length(x))
    n_modes = length(a_n) - 1
    for n in 1:n_modes
        omega_n = c * n * pi / L
        spatial  = @. sin(n * pi * x / L)
        temporal = a_n[n+1] * cos(omega_n * t) + b_n[n+1] * sin(omega_n * t)
        @. u += temporal * spatial
    end
    return u
end
```

@fig-wave-evolution shows the evolution of $u_N (x, t)$ with $N = 50$ modes at several time values within half a period $T = 2 L \/ c$. The string oscillates back and forth, with the triangular shape inverting at $t = T\/2$. Unlike the heat equation, the wave equation preserves energy and the solution does not decay; it continues oscillating indefinitely.

#figure(
  image("../figures/ch02/python/wave_evolution.pdf", width: 95%),
  caption: [Evolution of the wave equation solution with a plucked string initial condition. The string oscillates with period $T = 2 pi$, inverting at $t = T\/2$.],
) <fig-wave-evolution>

The code generating @fig-wave-evolution is available in:
- `codes/python/ch02/wave_equation_evolution.py`
- `codes/matlab/ch02/wave_equation_evolution.m`
- `codes/julia/ch02/wave_equation_evolution.jl`

The waterfall plot in @fig-wave-waterfall provides a complete view of the oscillatory dynamics over one full period. Unlike the heat equation, the wave equation conserves energy: the solution oscillates indefinitely without decay, and the periodic nature of the motion is clearly visible in the three-dimensional representation.

#figure(
  image("../figures/ch02/python/wave_waterfall.pdf", width: 95%),
  caption: [Waterfall plot showing the complete space-time evolution of the wave equation solution over one period $T$. The plucked string oscillates between its initial shape and its mirror image.],
) <fig-wave-waterfall>

The code generating @fig-wave-waterfall is available in:
- `codes/python/ch02/wave_equation_waterfall.py`
- `codes/matlab/ch02/wave_equation_waterfall.m`
- `codes/julia/ch02/wave_equation_waterfall.jl`

=== Discussion

The contrast with the heat equation is striking. In @fig-wave-evolution and @fig-wave-waterfall the triangular profile oscillates indefinitely without losing amplitude --- energy is conserved, and no smoothing takes place. The Fourier sine coefficients $a_n tilde.op 1\/n^2$ reflect the corner in the initial data, and because there is no temporal decay factor, those high-frequency components persist for all time. With $N = 50$ modes, faint Gibbs-like ripples remain near the kink, a reminder that truncation errors are not erased by the physics when there is no dissipation.

This étude highlights a key distinction that will recur throughout the book: _parabolic_ equations are forgiving of spectral truncation because diffusion damps the unresolved modes, whereas _hyperbolic_ equations demand higher resolution to maintain accuracy over long integration times.

== Laplace Equation in a Periodic Strip

For the elliptic case we consider the Laplace equation in a simple two dimensional domain that is periodic in one direction and bounded in the other. This setting connects naturally with the periodic heat equation example and again leads to a Fourier series representation in the periodic direction.

Let
$ D = { (x,y) in RR^2 : 0 < x < 2 pi, space 0 < y < 1 }. $

We seek a harmonic function $u(x,y)$ solving
$ u_(x x) (x,y) + u_(y y) (x,y) = 0, quad (x,y) in D, $
with periodic boundary conditions in $x$
$ u(x + 2 pi, y) = u(x,y), quad "for all real" space x, space 0 < y < 1, $
and Dirichlet conditions in $y$
$ u(x,0) = f(x), quad u(x,1) = 0, quad 0 lt.eq.slant x lt.eq.slant 2 pi. $

We assume that $f$ is $2 pi$ periodic and smooth:
$ f(x + 2 pi) = f(x). $

As in the parabolic and hyperbolic examples, we apply separation of variables and obtain a representation of $u$ as an infinite Fourier series in $x$ with $y$ dependent coefficients. The theory of harmonic functions and Laplace's equation is presented in detail in @Tikhonov1963. Modern spectral approaches for elliptic problems in complex domains include multi-domain methods @Yao2025 and fictitious domain techniques with circular embedding @GuZhou2021.

=== Step 1: Separation Ansatz

We look for nontrivial separated solutions of the form
$ u(x,y) = X(x) dot Y(y). $

Substituting into the Laplace equation gives
$ X''(x) dot Y(y) + X(x) dot Y''(y) = 0. $

Assuming $X$ and $Y$ are not identically zero, we divide by $X(x) dot Y(y)$:
$ frac(X''(x), X(x)) + frac(Y''(y), Y(y)) = 0. $

The first term depends only on $x$, the second only on $y$. Therefore both must be equal to constants whose sum is zero. We introduce a separation constant $lambda$ and write
$ frac(X''(x), X(x)) = -lambda, quad frac(Y''(y), Y(y)) = lambda. $

This leads to the two ordinary differential equations
$ X''(x) + lambda X(x) = 0, $
$ Y''(y) - lambda Y(y) = 0. $

The periodic boundary conditions for $u$ in the $x$ direction imply the periodic conditions
$ X(x + 2 pi) = X(x), quad "for all real" space x. $

The Dirichlet conditions in $y$ will be enforced later on the full series. For the separated functions $Y$ we will impose appropriate conditions at $y = 1$, while the condition at $y = 0$ will be handled via the Fourier coefficients of $f$.

We have once more an eigenvalue problem in the periodic direction.

=== Step 2: Eigenfunctions in the Periodic Direction

The equation for $X$ with periodic boundary conditions is exactly the same as in the heat equation example:
$ X''(x) + lambda X(x) = 0, quad X(x + 2 pi) = X(x). $

We recall the result: there is a constant mode with eigenvalue $lambda_0 = 0$,
$ X_0 (x) = 1, $
and for each integer $n gt.eq.slant 1$ there are cosine and sine modes with eigenvalues
$ lambda_n = n^2, $
$ X_n^((c)) (x) = cos(n x), quad X_n^((s)) (x) = sin(n x). $

The family
$ 1, cos(x), sin(x), cos(2 x), sin(2 x), dots $
forms an orthogonal basis in $L^2 (0, 2 pi)$ for $2 pi$ periodic functions.

For each such eigenvalue we now solve the corresponding equation for $Y$.

=== Step 3: Equations for the $y$ Dependent Factors

For each eigenvalue $lambda$ we have
$ Y''(y) - lambda Y(y) = 0. $

We treat separately the constant mode $lambda_0 = 0$ and the nonzero modes $lambda_n = n^2$ for $n gt.eq.slant 1$.

*Constant mode $lambda_0 = 0$*

For $lambda_0 = 0$ the equation reduces to
$ Y''_(0)(y) = 0. $

The general solution is
$ Y_0 (y) = A_0 + B_0 y. $

We want separated solutions that satisfy the homogeneous boundary condition at $y = 1$:
$ u(x,1) = 0 quad "for all" space x. $

For the constant mode this means
$ X_0 (x) dot Y_0 (1) = Y_0 (1) = 0, $
hence
$ Y_0 (1) = A_0 + B_0 = 0. $

We choose a convenient normalization so that $Y_0 (0) = 1$. Then $A_0 = 1$ and the relation $A_0 + B_0 = 0$ gives $B_0 = -1$. Thus
$ Y_0 (y) = 1 - y. $

This separated mode
$ u_0 (x,y) = X_0 (x) dot Y_0 (y) = 1 - y $
is harmonic, periodic in $x$, and vanishes at $y = 1$, with value $1$ at $y = 0$.

*Higher modes $lambda_n = n^2$ for $n gt.eq.slant 1$*

For $lambda_n = n^2$ the equation is
$ Y_(n)^('')(y) - n^2 Y_n (y) = 0. $

The general solution can be written in hyperbolic form
$ Y_n (y) = alpha_n cosh(n y) + beta_n sinh(n y). $

We require that each separated mode vanish at $y = 1$:
$ Y_n (1) = 0. $

To match later the Fourier coefficients of $f$ at $y = 0$, it is convenient to normalize so that
$ Y_n (0) = 1. $

Imposing $Y_n (0) = 1$ gives
$ alpha_n = 1. $

Then
$ Y_n (1) = cosh(n) + beta_n sinh(n) = 0 $
so
$ beta_n = - frac(cosh(n), sinh(n)) = -coth(n). $

Thus
$ Y_n (y) = cosh(n y) - coth(n) dot sinh(n y). $

An alternative and simpler expression uses the hyperbolic sine function of the distance to the boundary $y = 1$. One checks that
$ Y_n (y) = frac(sinh(n (1 - y)), sinh(n)) $
satisfies
$ Y_(n)^('')(y) - n^2 Y_n (y) = 0, $
$ Y_n (1) = 0, $
$ Y_n (0) = 1. $

Indeed, $Y_n (1) = sinh(0) \/ sinh(n) = 0$, $Y_n (0) = sinh(n) \/ sinh(n) = 1$, and
$ Y_(n)^('')(y) = n^2 frac(sinh(n (1 - y)), sinh(n)) = n^2 Y_n (y). $

We will use the form
$ Y_n (y) = frac(sinh(n (1 - y)), sinh(n)), quad n gt.eq.slant 1. $

=== Step 4: Building the Series Solution

Each separated solution corresponding to the eigenfunctions in $x$ and the functions $Y_n$ in $y$ has the form
$ u_0 (x,y) = C_0 dot Y_0 (y) = C_0 (1 - y), $
$ u_n^((c)) (x,y) = C_n dot cos(n x) dot Y_n (y), $
$ u_n^((s)) (x,y) = D_n dot sin(n x) dot Y_n (y), quad n gt.eq.slant 1, $
for some constants $C_0$, $C_n$, $D_n$.

Since the Laplace equation is linear and the boundary condition at $y = 1$ is homogeneous, any linear combination of these separated solutions is again a solution that vanishes at $y = 1$. Therefore a general solution satisfying the periodic condition in $x$ and the Dirichlet condition $u(x,1) = 0$ can be written as
$ u(x,y) = C_0 (1 - y) + sum_(n=1)^infinity [C_n cos(n x) + D_n sin(n x)] frac(sinh(n (1 - y)), sinh(n)). $

It remains to impose the boundary condition at $y = 0$,
$ u(x,0) = f(x). $

At $y = 0$ we obtain
$ u(x,0) = C_0 + sum_(n=1)^infinity [C_n cos(n x) + D_n sin(n x)], $
because $Y_0 (0) = 1$ and $Y_n (0) = 1$ for $n gt.eq.slant 1$.

Hence the boundary condition $u(x,0) = f(x)$ becomes
$ f(x) = C_0 + sum_(n=1)^infinity [C_n cos(n x) + D_n sin(n x)]. $

This is exactly the Fourier series expansion of $f$. For a $2 pi$ periodic $f$ we have
$ f(x) = a_0 + sum_(n=1)^infinity [a_n cos(n x) + b_n sin(n x)], $
with coefficients
$ a_0 = frac(1, 2 pi) integral_0^(2 pi) f(x) dif x, $
$ a_n = frac(1, pi) integral_0^(2 pi) f(x) cos(n x) dif x, quad n gt.eq.slant 1, $
$ b_n = frac(1, pi) integral_0^(2 pi) f(x) sin(n x) dif x, quad n gt.eq.slant 1. $

By uniqueness of the Fourier expansion, we must have
$ C_0 = a_0, quad C_n = a_n, quad D_n = b_n. $

=== Step 5: Final Explicit Solution and Interpretation

Substituting these coefficients into the series we obtain the explicit representation
$ u(x,y) = a_0 (1 - y) + sum_(n=1)^infinity [a_n cos(n x) + b_n sin(n x)] frac(sinh(n (1 - y)), sinh(n)), quad 0 < y < 1. $

Here $a_0$, $a_n$, and $b_n$ are the Fourier coefficients of the boundary data $f$ as defined above.

This series converges (under mild assumptions on $f$) to the unique harmonic function that is periodic in $x$, equal to $f$ on $y = 0$, and zero on $y = 1$.

From a spectral viewpoint:

- The functions $1$, $cos(n x)$, $sin(n x)$ are eigenfunctions of the one dimensional Laplacian $X arrow.r.bar X''$ with periodic boundary conditions in $x$, with eigenvalues $lambda_0 = 0$ and $lambda_n = n^2$.

- For each spatial frequency $n$ in the periodic direction, the dependence in the transverse direction $y$ is determined by the simple ordinary differential equation $Y'' - lambda_n Y = 0$ with boundary condition $Y(1) = 0$ and normalization $Y(0) = 1$. This gives the hyperbolic profiles
$ Y_0 (y) = 1 - y, $
$ Y_n (y) = frac(sinh(n (1 - y)), sinh(n)), quad n gt.eq.slant 1. $

- The boundary data $f$ at $y = 0$ is expanded in the eigenbasis ${1, cos(n x), sin(n x)}$ and each Fourier mode is propagated into the interior of the strip with its own $y$ dependent factor $Y_n (y)$.

Analytically, the solution is an infinite sum of separated solutions. In a spectral method we will truncate this sum to finitely many modes in $x$,
$ u_N (x,y) = a_0 (1 - y) + sum_(n=1)^N [a_n cos(n x) + b_n sin(n x)] frac(sinh(n (1 - y)), sinh(n)), $
and approximate the harmonic function inside the strip by this finite Fourier representation. Related spectral techniques for Stokes and Laplace problems in confined periodic domains are developed in @Pelaez2025.

== Computational Étude 2.3: Steady-State Temperature in a Strip

To visualize the structure of harmonic functions in the strip, we compute the truncated Fourier series solution for a boundary condition that contains two modes:
$ f(x) = sin(x) + frac(1,2) sin(3 x). $

For this particular boundary data, the Fourier coefficients are simply $b_1 = 1$ and $b_3 = 1\/2$, with all other coefficients zero. The solution can be written explicitly as
$ u(x,y) = sin(x) frac(sinh(1 - y), sinh(1)) + frac(1,2) sin(3 x) frac(sinh(3(1-y)), sinh(3)). $

The key portion of the implementation evaluates the truncated series on a two-dimensional grid. In Python:

```python
def laplace_solution(x, y, a0, a_n, b_n):
    u = a0 * (1 - y)
    n_modes = len(a_n) - 1
    for n in range(1, n_modes + 1):
        if abs(a_n[n]) < 1e-15 and abs(b_n[n]) < 1e-15:
            continue
        y_factor = np.sinh(n * (1 - y)) / np.sinh(n)
        u += (a_n[n] * np.cos(n * x) + b_n[n] * np.sin(n * x)) * y_factor
    return u
```

The equivalent MATLAB implementation:

```matlab
U = a0 * (1 - Y);
for n = 1:N_MODES
    if abs(a_n(n+1)) < 1e-15 && abs(b_n(n+1)) < 1e-15
        continue;
    end
    y_factor = sinh(n * (1 - Y)) / sinh(n);
    U = U + (a_n(n+1) * cos(n*X) + b_n(n+1) * sin(n*X)) .* y_factor;
end
```

The Julia implementation:

```julia
function laplace_solution(X, Y, a0, a_n, b_n)
    U = a0 .* (1.0 .- Y)
    n_modes = length(a_n) - 1
    for n in 1:n_modes
        if abs(a_n[n+1]) < 1e-15 && abs(b_n[n+1]) < 1e-15
            continue
        end
        y_factor = sinh.(n .* (1.0 .- Y)) ./ sinh(n)
        @. U += (a_n[n+1] * cos(n * X) + b_n[n+1] * sin(n * X)) * y_factor
    end
    return U
end
```

@fig-laplace-solution shows the solution $u(x,y)$ in the strip $[0, 2 pi] times [0, 1]$. At the bottom boundary $y = 0$, the solution matches the prescribed boundary data $f(x)$. As $y$ increases toward the top boundary, the solution decays to zero. Crucially, the higher frequency mode ($n = 3$) decays much faster than the lower frequency mode ($n = 1$), as the hyperbolic factor $sinh(n(1-y))\/sinh(n)$ decreases more rapidly for larger $n$. This illustrates the smoothing effect of harmonic extension into the interior.

#figure(
  image("../figures/ch02/python/laplace_solution.pdf", width: 95%),
  caption: [Solution of the Laplace equation in the periodic strip with boundary data $f(x) = sin(x) + frac(1,2) sin(3x)$ at $y = 0$ and $u = 0$ at $y = 1$. Higher frequency modes decay faster toward the interior.],
) <fig-laplace-solution>

The code generating @fig-laplace-solution is available in:
- `codes/python/ch02/laplace_equation_2d.py`
- `codes/matlab/ch02/laplace_equation_2d.m`
- `codes/julia/ch02/laplace_equation_2d.jl`

=== Discussion

@fig-laplace-solution reveals the interior smoothing characteristic of elliptic equations. The boundary data $f(x) = sin(x) + frac(1,2) sin(3x)$ contains two modes, yet the higher-frequency component ($n = 3$) is almost entirely suppressed by the time $y approx 0.3$ because the hyperbolic factor $sinh(n(1-y))\/sinh(n)$ decays much faster for larger $n$. The Laplace equation thus acts as a natural low-pass filter in the transverse direction: the further one moves from the boundary, the smoother the solution becomes.

This behaviour previews a fundamental principle of spectral methods for elliptic problems: boundary data with many high-frequency components generates solutions that are smooth in the interior, so moderate truncation orders $N$ suffice for excellent accuracy away from the boundary. The analytical solution derived here will serve as a benchmark when we later solve the Laplace equation numerically using Chebyshev collocation in @ch-bvp and Kronecker product methods in @ch-advanced-bc.

== A non-exhaustive literature overview

The method of separation of variables and the spectral methods derived from it share a common intellectual lineage that spans three centuries. This overview traces that arc from the foundational debates of the Enlightenment to the high-performance computing algorithms of the present day.

*The Historical Foundations.* The techniques presented in this chapter originated in the "vibrating string controversy" of the mid-18th century, a debate that defined the very concept of a mathematical function. In 1747, Jean le Rond d'Alembert derived the wave equation and proposed a solution based on traveling waves, $u(x,t) = f(x - c t) + g(x + c t)$, but argued that the functions must be analytically given @dAlembert1747. Leonhard Euler expanded this view in 1748, arguing for "discontinuous" functions (in the sense of having corners, like the plucked string) @Euler1748. It was Daniel Bernoulli who, in 1753, proposed that the solution could be represented as an infinite sum of trigonometric modes --- the series solutions used in this chapter.

The legitimacy of these infinite series was firmly established by Joseph Fourier in his 1822 treatise _Théorie Analytique de la Chaleur_ @Fourier1822. Fourier demonstrated that even discontinuous functions could be represented by trigonometric series, applying the method to the heat equation. This work laid the groundwork for Sturm--Liouville theory @Sturm1836 @SturmLiouville1837, which generalized the concept of eigenfunctions to a broad class of second-order operators.

*The Computational Era.* The transition from analytical separation of variables to numerical spectral methods began with Cornelius Lanczos in 1938, who utilized Chebyshev polynomials for differential equations (the Tau method) @Lanczos1938. The modern era of spectral methods was launched by the rediscovery of the Fast Fourier Transform (FFT) by Cooley and Tukey in 1965 @Cooley1965, which reduced the computational cost of transformations from $O(N^(2))$ to $O(N log N)$. In the 1970s, Steven Orszag and David Gottlieb codified the theory of spectral methods, demonstrating their superior convergence properties for smooth flow problems compared to finite difference methods @Orszag1971 @GottliebOrszag1977.

*State of the Art (2024--2025).* Today, spectral methods remain at the forefront of computational research. Recent developments have focused on overcoming historical limitations regarding domain geometry and temporal discretization.

- *Complex Geometries:* Yao, Wang, and Wang @Yao2025 have developed multi-domain Legendre spectral methods capable of solving wave equations in exterior domains with complex obstacles. Similarly, fictitious domain methods with circular embedding have been refined to solve elliptic PDEs in irregular geometries @GuZhou2021.
- *Space-Time Methods:* Moving beyond the traditional method of lines, Kaur _et al._ @Kaur2025 have introduced space-time spectral methods that treat time as a fourth spectral dimension, achieving high-order accuracy in both space and time for linear and nonlinear PDEs.
- *Non-Local Operators:* Mustapha _et al._ @Mustapha2025 have extended Fourier spectral methods to non-local equations and fractional Laplacians on bounded domains, addressing the unique challenges of boundary conditions in non-local mechanics.
- *Fluid Dynamics:* In the _Journal of Fluid Mechanics_, Peláez _et al._ @Pelaez2025 recently applied spectral solvers to oscillatory Stokes flow in doubly periodic confined domains, demonstrating the method's continued vitality in microfluidics research.

These recent works illustrate that the fundamental principle of this chapter --- expanding a solution in a global basis of eigenfunctions --- remains a powerful and evolving tool in modern scientific computing.

== Summary

The three examples presented in this chapter (the heat equation, the wave equation, and the Laplace equation) share a common mathematical structure that will guide us throughout the rest of this book. In each case, separation of variables reduces a partial differential equation to a family of ordinary differential equations: an eigenvalue problem in the spatial variable and a simpler equation governing the behavior in the remaining variable (time or the transverse coordinate). The eigenfunctions form an orthogonal basis, and the solution is expressed as an infinite series
$ u(x,t) = sum_(k) hat(u)_k (t) phi_k (x) $
whose coefficients are determined by the initial or boundary data through Fourier projections. This is the DNA of spectral methods.

Let us distill the key ideas:

+ *Separation*. We decompose the solution into modes that evolve independently (or nearly so). Each mode satisfies a simpler equation than the original PDE.

+ *Basis*. The spatial part of each mode is an eigenfunction of a differential operator: Fourier exponentials for periodic problems, trigonometric functions for Dirichlet conditions, Chebyshev or Legendre polynomials for more general settings.

+ *Truncation*. In practice we cannot sum infinitely many terms due to the apparent finitude of our Universe. We retain only the first $N$ modes, and the accuracy of this approximation depends critically on how fast the coefficients $hat(u)_k$ decay.

When the solution is smooth, the coefficients decay _exponentially_, and a modest $N$ suffices for high accuracy. This is the source of spectral methods' legendary efficiency. For a comprehensive treatment of these classical methods and their mathematical foundations, the reader is referred to @Tikhonov1963.

The analytical solutions derived in this chapter are beautiful, but they are also fragile. They apply only to linear equations on simple domains with special boundary conditions. The moment we encounter a nonlinearity, a complicated geometry, or variable coefficients, we must turn to computation. The chapters that follow will develop the computational machinery needed to turn these analytical insights into practical algorithms:

- *FFT and its cousins*: how to move efficiently between physical space and coefficient space.
- *Differentiation matrices*: how to compute derivatives spectrally.
- *Quadrature rules*: how to compute inner products and projections.
- *Time stepping*: how to advance the ODE system for the coefficients.

Armed with these tools, we will be able to solve problems far beyond the reach of pen-and-paper analysis.

== Exercises <sec-classical-pdes-exercises>

*Exercise 2.1* (_Separation of Variables for the Wave Equation_). Consider the one-dimensional wave equation $u_(t t) = c^2 u_(x x)$ on $[0, pi]$ with homogeneous Dirichlet boundary conditions $u(0, t) = u(pi, t) = 0$ and initial conditions $u(x, 0) = sin(x) + (1\/2) sin(3x)$, $u_t (x, 0) = 0$. (a) Derive the analytical solution using separation of variables. (b) Verify that the solution is periodic in time and determine the recurrence period. (c) If $c = 1$ and you sample the solution at $N = 16$ equispaced interior points in $x$, what is the maximum time step $Delta t$ permitted by the CFL condition $c Delta t \/ Delta x lt.eq.slant 1$?

*Exercise 2.2* (_Fourier Convergence for Discontinuous Data_). The function $f(x) = 1$ for $0 < x < pi\/2$ and $f(x) = 0$ for $pi\/2 < x < pi$ has the Fourier sine series $f(x) = sum_(k=1)^infinity b_k sin(k x)$. (a) Compute the coefficients $b_k$ analytically. (b) Plot the partial sums $S_N (x) = sum_(k=1)^N b_k sin(k x)$ for $N = 5, 10, 20, 50$ and observe the Gibbs phenomenon near $x = pi\/2$. (c) Determine the rate at which $b_k$ decays and explain why algebraic (not exponential) convergence is expected for this non-smooth function.

*Exercise 2.3* (_Heat Equation with Non-homogeneous Boundary Conditions_). Solve the heat equation $u_t = u_(x x)$ on $[0, 1]$ with boundary conditions $u(0, t) = 0$, $u(1, t) = 1$ and initial condition $u(x, 0) = x + sin(pi x)$. (a) Introduce the steady-state solution $u_s (x) = x$ and write the PDE for the transient part $v(x, t) = u(x, t) - u_s (x)$. (b) Solve for $v(x, t)$ by separation of variables. (c) How many Fourier modes are needed to approximate the solution to within $10^(-6)$ at $t = 0.1$?

*Exercise 2.4* (_Comparing Analytical and Numerical Fourier Solutions_). For the heat equation $u_t = u_(x x)$ on $[0, pi]$ with $u(0, t) = u(pi, t) = 0$ and $u(x, 0) = sin(x)$: (a) Write down the exact solution. (b) Implement a numerical solver using the truncated Fourier series with $N$ modes and explicit Euler time stepping in the coefficient space. (c) Compare the numerical and analytical solutions at $t = 0.5$ for $N = 4, 8, 16$ and show that the spatial error is dominated by the truncation of the Fourier series. (d) Investigate the stability constraint on $Delta t$ as a function of $N$ and explain why explicit time stepping in spectral space becomes increasingly restrictive.
