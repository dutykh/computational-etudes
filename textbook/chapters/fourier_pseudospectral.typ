// textbook/chapters/fourier_pseudospectral.typ
// Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Email: denys.dutykh@ku.ac.ae
// Homepage: https://www.denys-dutykh.com/
// Last modified: February 2026

#import "../styles/template.typ": dropcap, num, format-table, etude-conclusion, idx

// Enable equation numbering for this chapter

= Fourier Pseudospectral Methods for Periodic PDEs <ch-fourier-pseudo>

#dropcap[The spectral PDE solvers of the previous chapter relied on Chebyshev grids, which excel on bounded, non-periodic domains but impose severe CFL constraints on explicit time stepping. When the geometry is periodic, there is a far more natural alternative: _Fourier pseudospectral method#idx("pseudospectral method")s_. The core idea is beguilingly simple: compute derivatives in Fourier space by multiplication with $i k$, evaluate nonlinear products pointwise in physical space, and shuttle between the two representations with the FFT. The result is a family of solvers that combine spectral accuracy in space with remarkably mild stability restrictions in time, and that can simulate soliton#idx("soliton") collisions, spatiotemporal chaos#idx("spatiotemporal chaos"), and two-dimensional turbulence in astonishingly compact code.]

By the end of this chapter, you should be able to:

1. Derive Fourier pseudospectral semidiscretisations for one- and two-dimensional periodic PDEs.
2. Implement spectral derivatives, nonlinear terms, and Poisson solves efficiently via FFT.
3. Understand CFL and stability restrictions for explicit schemes on equispaced periodic grids, and contrast them with the Chebyshev case.
4. Control aliasing in nonlinear terms using the two-thirds rule#idx("two-thirds rule") and zero padding.
5. Overcome stiffness in dispersive and high-order PDEs using integrating factor#idx("integrating factor") and exponential time differencing#idx("exponential time differencing") schemes.
6. Simulate, from scratch: advection, viscous Burgers fronts, KdV soliton collisions, nonlinear Schrödinger recurrence, Kuramoto--Sivashinsky chaos, and two-dimensional Navier--Stokes vortex dynamics.

== Fourier Pseudospectral Derivatives <sec-fourier-ps-deriv>

We begin by recalling, from @ch-fourier-grids, the essential machinery. Let $N in 2 NN$ and define the equispaced periodic grid
$ x_j = frac(2 pi j, N), quad j = 0, 1, dots, N - 1, $ <eq-periodic-grid>
with spacing $h = 2 pi \/ N$. The discrete Fourier transform (DFT) of a grid function $bold(u) = (u_0, dots, u_(N-1))^top$ is
$ hat(u)_k = frac(1, N) sum_(j = 0)^(N-1) u_j e^(-i k x_j), quad k in K_N, $ <eq-dft>
where $K_N = {0, 1, dots, N\/2, -N\/2+1, dots, -1}$ is the standard FFT wavenumber set.

For a smooth periodic function, the Fourier series $u(x) = sum_k hat(u)_k e^(i k x)$ satisfies
$ u_x (x) = sum_k i k hat(u)_k e^(i k x), quad u_(x x) (x) = sum_k (-k^2) hat(u)_k e^(i k x). $
On the grid, spectral differentiation is thus implemented by:

+ *Forward FFT*: compute $hat(u)_k$ from $u_j$.
+ *Multiply*: form $i k hat(u)_k$ (first derivative) or $-k^2 hat(u)_k$ (second derivative).
+ *Inverse FFT*: recover the derivative values at grid points.

This FFT-based approach is intimately connected to the _circulant differentiation matrix_ $D$ from @ch-differentiation. That matrix has eigenvalues $i k$ for $k in K_N$, and the FFT diagonalises it @Cooley1965. In other words, the FFT turns an $O(N^2)$ matrix-vector product into $O(N log N)$ operations @KreissOliger1972 @Fornberg1987. In this chapter we will _never_ build $D$ explicitly (except for pedagogical comparison at small $N$), but the eigenvalue picture is essential for stability analysis.

=== Computational Étude 11.1: Fourier Differentiation, Matrix versus FFT <etude-fourier-diff>

To make the connection concrete, we differentiate the test function $u(x) = e^(sin x)$ on $[0, 2 pi)$ using both the dense matrix $D$ and the FFT approach. The exact derivative is $u'(x) = cos x dot e^(sin x)$.

The core of the FFT approach is a three-line computation. In Python:

```python
def fourier_diff_fft(u):
    """Compute du/dx on a periodic grid via FFT."""
    N = len(u)
    u_hat = np.fft.fft(u)
    k = np.fft.fftfreq(N, d=1.0 / N)
    if N % 2 == 0:
        k[N // 2] = 0  # zero Nyquist for odd derivative
    return np.fft.ifft(1j * k * u_hat).real
```

The equivalent MATLAB implementation:

```matlab
function du = fourier_diff_fft(u)
    N = length(u);
    u_hat = fft(u);
    k = [0:N/2-1, 0, -N/2+1:-1]';  % zero Nyquist
    du = real(ifft(1i * k .* u_hat));
end
```

The Julia implementation:

```julia
using FFTW

function fourier_diff_fft(u::Vector{Float64})
    N = length(u)
    u_hat = fft(u)
    k = vcat(0:N÷2-1, 0, -N÷2+1:-1)
    real.(ifft(im .* k .* u_hat))
end
```

#figure(
  image("../figures/ch11/python/fourier_diff_accuracy.pdf", width: 95%),
  caption: [Spectral convergence of Fourier differentiation for the test function $u(x) = e^(sin x)$. Both the matrix approach ($D bold(u)$) and the FFT approach ($cal(F)^(-1)(i k hat(u))$) achieve machine precision for $N gt.eq.slant 30$, confirming that the two methods are algebraically equivalent.],
) <fig-fourier-diff-accuracy>

@fig-fourier-diff-accuracy shows that both methods converge to machine precision at the same rate: the error drops exponentially with $N$ and saturates near $10^(-15)$. @fig-fourier-diff-timing reveals the payoff: the matrix method scales as $O(N^2)$, while the FFT method scales as $O(N log N)$, making it orders of magnitude faster for large $N$.

#figure(
  image("../figures/ch11/python/fourier_diff_timing.pdf", width: 85%),
  caption: [Computational cost of Fourier differentiation. The dense matrix approach scales as $O(N^2)$; the FFT approach as $O(N log N)$. For $N = 2048$ the FFT is roughly 1000 times faster.],
) <fig-fourier-diff-timing>

#etude-conclusion[
  The convergence curves confirm the *algebraic equivalence* of the matrix and FFT routes to spectral differentiation: both reach machine precision at the same modest $N$. Yet the timing data expose a practical gap that widens dramatically with resolution: at $N = 2048$ the FFT is roughly *three orders of magnitude* faster than the dense matrix-vector product, and the asymptotic separation only grows for larger $N$. Every remaining étude in this chapter computes derivatives exclusively through the FFT. The exercise illustrates a recurring theme: two formulations that are mathematically identical can differ enormously in cost, and choosing the right representation --- here the eigenvalue decomposition furnished by the FFT --- is often the decisive step.
]

The code generating @fig-fourier-diff-accuracy and @fig-fourier-diff-timing is available in:
- `codes/python/ch11/fourier_diff_comparison.py`
- `codes/matlab/ch11/fourier_diff_comparison.m`
- `codes/julia/ch11/fourier_diff_comparison.jl`

== Method of Lines and Semidiscretisation <sec-mol-fourier>

With spectral differentiation in hand, we can discretise time-dependent PDEs on periodic domains by the _method of lines_ @Schiesser1991: replace spatial derivatives with their Fourier pseudospectral approximations, leaving time continuous. The result is a system of ordinary differential equations (ODEs) of the form
$ dot(bold(u))(t) = bold(F)(bold(u)(t)), $ <eq-mol>
where $bold(F)$ computes spatial derivatives via FFTs. This is precisely the framework introduced in @ch-spectral-pde, but with a crucial simplification: on a periodic, equispaced grid, the differentiation operator is _circulant_, and its eigenvalues are simply $i k$ (or $-k^2$, $-i k^3$, etc.). This makes the analysis and implementation far cleaner than the Chebyshev case.

=== Linear Constant-Coefficient Examples

For a linear PDE $u_t = cal(L) u$ with constant coefficients, the Fourier transform decouples the system completely. Each mode evolves independently:

- *Advection* ($u_t + c u_x = 0$): each Fourier mode rotates on the complex plane,
  $ hat(u)_k '(t) = -i c k hat(u)_k (t). $ <eq-advection-mode>
  The solution is $hat(u)_k (t) = hat(u)_k (0) e^(-i c k t)$: pure transport without distortion.

- *Diffusion* ($u_t = kappa u_(x x)$): each mode decays exponentially,
  $ hat(u)_k '(t) = -kappa k^2 hat(u)_k (t). $ <eq-diffusion-mode>
  High wavenumbers decay fastest, explaining the smoothing effect of diffusion.

This decoupling is the fundamental reason why _linear, constant-coefficient, periodic PDEs are exactly diagonalised by the Fourier transform_. Spatial discretisation introduces only spectral truncation error (the omission of modes with $|k| > N\/2$), not the additional approximation errors that arise from non-periodic grids.

== Time Stepping and Stability on Periodic Grids <sec-stability-periodic>

The semidiscrete system @eq-mol must be integrated in time. The choice of time stepper and the maximum allowable step size are governed by _stability regions_, a concept introduced in @ch-spectral-pde and treated in detail by Trefethen @Trefethen2000 and in the classical monograph of Gottlieb and Orszag @GottliebOrszag1977. The CFL analysis for spectral methods was further developed by Gottlieb and Turkel @GottliebTurkel1980.

=== Eigenvalues of Fourier Operators

For the first-derivative operator on the periodic grid, the eigenvalues are
$ lambda_k = -i c k, quad k in K_N. $
These lie on the imaginary axis. For the second-derivative operator,
$ lambda_k = -kappa k^2, $
which lie on the negative real axis. The maximum eigenvalue magnitude is set by the highest resolved wavenumber $k_max = N\/2$:

- Advection: $|lambda_max| = |c| N\/2$.
- Diffusion: $|lambda_max| = kappa (N\/2)^2$.

=== CFL Conditions for Fourier Pseudospectral Methods

For stability, the scaled eigenvalues $lambda_k Delta t$ must lie inside the stability region of the time stepper. For the classical fourth-order Runge--Kutta method (RK4), the stability region extends approximately 2.8 units along the imaginary axis and 2.8 units along the negative real axis.

This gives the following scaling for the maximum time step:

- *Advection (first derivative)*: $Delta t lt.eq.slant C_1 \/ N$, where $C_1$ depends on the time stepper.
- *Diffusion (second derivative)*: $Delta t lt.eq.slant C_2 \/ N^2$.

These scalings are dramatically better than the Chebyshev case, where the clustering of grid points near the boundaries forces $Delta t tilde.op N^(-2)$ even for advection, and $Delta t tilde.op N^(-4)$ for diffusion. *Fourier pseudospectral methods on periodic grids are far friendlier to explicit time stepping than their Chebyshev counterparts.*

=== Exact Integration for Linear Problems

For linear, constant-coefficient problems, one can bypass time stepping entirely. Since each Fourier mode satisfies $hat(u)_k (t) = e^(lambda_k t) hat(u)_k (0)$, one has
$ hat(u)_k (t + Delta t) = e^(lambda_k Delta t) hat(u)_k (t). $ <eq-exact-integration>
This "exponential time stepping" is unconditionally stable and introduces _no_ temporal discretisation error for the linear part. It is the prototype for the integrating factor and exponential time differencing methods that we develop in @sec-integrating-factor for nonlinear problems.

=== Computational Étude 11.2: Non-Dispersive Fourier Advection <etude-advection>

We illustrate these ideas with the simplest possible PDE: $u_t + c u_x = 0$ on $[0, 2 pi)$ with $c = 1$ and the smooth initial condition $u(x, 0) = exp(-10 (x - pi)^2)$. We solve with Fourier pseudospectral + RK4 and compare against a second-order upwind finite difference scheme on the same grid.

In Python, the right-hand side for the Fourier scheme is:

```python
def rhs_fourier(u, k, c):
    """Right-hand side: -c * u_x via FFT."""
    u_hat = np.fft.fft(u)
    ux = np.fft.ifft(1j * k * u_hat).real
    return -c * ux
```

The equivalent MATLAB implementation:

```matlab
function f = rhs_fourier(u, k, c)
    u_hat = fft(u);
    ux = real(ifft(1i * k .* u_hat));
    f = -c * ux;
end
```

The Julia implementation:

```julia
function rhs_fourier(u, k, c)
    u_hat = fft(u)
    ux = real.(ifft(im .* k .* u_hat))
    return -c .* ux
end
```

#figure(
  image("../figures/ch11/python/advection_comparison.pdf", width: 95%),
  caption: [Fourier pseudospectral versus second-order upwind finite differences for linear advection after one full period ($t = 2 pi$, $N = 128$). _Left_: The spectral solution (blue) is visually indistinguishable from the exact initial condition. The finite difference solution (red dashes) shows significant dispersion and dissipation. _Right_: Pointwise error on a logarithmic scale, confirming spectral accuracy (errors near machine precision) versus second-order accuracy.],
) <fig-advection>

The code generating @fig-advection is available in:
- `codes/python/ch11/advection_fourier.py`
- `codes/matlab/ch11/advection_fourier.m`
- `codes/julia/ch11/advection_fourier.jl`

#etude-conclusion[
  After one full period, the spectral solution returns to the initial condition with errors near machine precision, while the finite-difference solution has been visibly degraded by numerical dispersion and dissipation. This is the *spectral promise* in action. The comparison is especially instructive because both methods use the *same* number of grid points: what differs is not resolution but the quality of the spatial approximation. A second-order upwind scheme needs roughly $10$--$20$ points per wavelength; the Fourier method needs only $pi$ @KreissOliger1972. For smooth periodic data the gap between low-order and spectral methods is not merely quantitative but *qualitative* --- algebraic versus faster-than-any-power decay.
]

== Nonlinearities, Aliasing, and Dealiasing <sec-aliasing>

The true power of pseudospectral methods emerges with _nonlinear_ PDEs. The pseudospectral strategy for evaluating a nonlinear term such as $u u_x$ is:

+ Transform $u$ to physical space (if not already there).
+ Compute the derivative $u_x$ in Fourier space, then transform back.
+ Form the product $u dot u_x$ _pointwise_ in physical space.
+ Transform the result back to Fourier space.

This costs $O(N log N)$ per time step, compared with $O(N^2)$ for direct convolution in Fourier space. However, the pointwise product introduces _aliasing_: when two Fourier modes with wavenumbers $k_1$ and $k_2$ are multiplied, the result has wavenumber $k_1 + k_2$, which may lie outside the resolved range $[-N\/2, N\/2]$. Such out-of-range modes "fold back" into the resolved band, corrupting the low-frequency content.

=== The Two-Thirds Dealiasing Rule

For quadratic nonlinearities, the classic fix is the _Orszag two-thirds rule_ @Orszag1971: before transforming the nonlinear product back to Fourier space, zero out all modes with $|k| > N\/3$. This ensures that the product of any two resolved modes (each with $|k| lt.eq.slant N\/3$) produces modes that remain within the Nyquist limit.

In practice, this is implemented as a mask applied to the Fourier coefficients of the nonlinear term:
$
tilde(hat(f))_k = cases(hat(f)_k\, & "if" |k| lt.eq.slant N\/3, 0\, & "if" N\/3 < |k| lt.eq.slant N\/2.)
$ <eq-dealias-mask>

An alternative approach is _zero-padding_: embed the data in a grid of size $M gt.eq.slant 3 N\/2$, compute the nonlinearity on the larger grid, then truncate back to $N$ modes. Both approaches are equivalent for quadratic terms; see Patterson and Orszag @PattersonOrszag1971 for the original analysis.

=== Computational Étude 11.3: Viscous Burgers Equation <etude-burgers>

The viscous Burgers equation#idx("Burgers equation")
$ u_t + u u_x = nu u_(x x), quad x in [0, 2 pi), $ <eq-burgers>
with periodic boundary conditions and initial condition $u_0 (x) = sin x$, is an ideal test bed for the pseudospectral treatment of nonlinearities. The balance between nonlinear steepening ($u u_x$) and viscous smoothing ($nu u_(x x)$) produces sharp fronts that challenge numerical methods. Tadmor @Tadmor1989 analysed the spectral viscosity approach for conservation laws related to this problem.

The semidiscrete system in Fourier space is
$ hat(u)_k '(t) = -hat(u u_x)_k (t) - nu k^2 hat(u)_k (t), $
where $hat(u u_x)$ is computed pseudospectrally (with optional dealiasing#idx("dealiasing")) and the diffusive term is evaluated directly in Fourier space.

For the time step, we must accommodate both the advective restriction ($Delta t tilde.op 1\/N$) and the diffusive restriction ($Delta t tilde.op 1\/(nu N^2)$), taking the smaller of the two.

In Python, the right-hand side including dealiasing is:

```python
def rhs(u_phys):
    u_hat = np.fft.fft(u_phys)
    ux = np.fft.ifft(1j * k * u_hat).real
    f_nl = -u_phys * ux
    f_nl_hat = np.fft.fft(f_nl)
    f_nl_hat[~dealias_mask] = 0.0
    f_diff_hat = -nu * k**2 * u_hat
    return np.fft.ifft(f_nl_hat + f_diff_hat).real
```

The equivalent MATLAB implementation:

```matlab
function f = rhs(u_phys)
    u_hat = fft(u_phys);
    ux = real(ifft(1i * k .* u_hat));
    f_nl = -u_phys .* ux;
    f_nl_hat = fft(f_nl);
    f_nl_hat(~mask) = 0;
    f_diff_hat = -nu * k.^2 .* u_hat;
    f = real(ifft(f_nl_hat + f_diff_hat));
end
```

The Julia implementation:

```julia
function rhs(u_phys, k, nu, dealias_mask)
    u_hat = fft(u_phys)
    ux = real.(ifft(im .* k .* u_hat))
    f_nl = -u_phys .* ux
    f_nl_hat = fft(f_nl)
    f_nl_hat[.!dealias_mask] .= 0.0
    f_diff_hat = -nu .* k.^2 .* u_hat
    return real.(ifft(f_nl_hat .+ f_diff_hat))
end
```

#figure(
  image("../figures/ch11/python/burgers_evolution.pdf", width: 95%),
  caption: [Viscous Burgers equation @eq-burgers with $nu = 0.01$ and initial condition $u_0 (x) = sin x$. _Left_: Space-time diagram showing the steepening of the initial sine wave into a sharp front. _Right_: Snapshots at selected times illustrating the progressive formation of the viscous shock structure.],
) <fig-burgers-evolution>

#figure(
  image("../figures/ch11/python/burgers_dealiasing.pdf", width: 95%),
  caption: [Effect of dealiasing on the Burgers equation at $t = 1.0$. _Left_: Energy spectra $|hat(u)_k|^2$ with and without the two-thirds rule. Without dealiasing, energy piles up near the Nyquist frequency. _Right_: Physical space solutions; the differences are subtle here but grow for longer integration times or smaller viscosities.],
) <fig-burgers-dealiasing>

The code generating @fig-burgers-evolution and @fig-burgers-dealiasing is available in:
- `codes/python/ch11/burgers_fourier.py`
- `codes/matlab/ch11/burgers_fourier.m`
- `codes/julia/ch11/burgers_fourier.jl`

#etude-conclusion[
  The first figure shows the classic Burgers behaviour: the initial sine steepens nonlinearly into a narrow front near $x = pi$, regularised by viscosity. The second figure exposes the effect of aliasing on the energy spectrum: *without the two-thirds rule, energy accumulates near the highest resolved wavenumbers* --- a signature of aliasing errors that can destabilise long-time computations. The lesson is one that pure theory can only hint at: the *visual* difference between aliased and dealiased solutions may look small at moderate times, but the spectral diagnostic already shows a qualitative defect. The two-thirds rule is *not* a theoretical nicety; it is a practical necessity for any pseudospectral code intended for quantitative work.
]

== Overcoming Stiffness: Integrating Factors and ETD <sec-integrating-factor>

The PDEs studied so far (advection, Burgers) have relatively mild CFL restrictions because their highest-order spatial derivative is second order. Many physically important equations, however, contain third- or fourth-order derivatives (dispersion, biharmonic diffusion) whose eigenvalues grow as $k^3$ or $k^4$. For the Kuramoto--Sivashinsky equation ($u_(x x x x)$ term), explicit methods require $Delta t tilde.op N^(-4)$; for $N = 128$ this means $Delta t approx 10^(-9)$, rendering simulations impractical.

The _integrating factor_ (IF) method, originating with Lawson @Lawson1967, eliminates this stiffness by treating the linear part exactly. Consider a PDE of the form
$ hat(u)_k '(t) = lambda_k hat(u)_k (t) + hat(cal(N))_k (t), $ <eq-linear-nonlinear>
where $lambda_k$ is the linear symbol (e.g., $i k^3$ for KdV, $k^2 - k^4$ for KS) and $hat(cal(N))_k$ is the nonlinear part. Define the integrating factor variable
$ hat(v)_k (t) = e^(-lambda_k t) hat(u)_k (t). $
Then
$ hat(v)_k '(t) = e^(-lambda_k t) hat(cal(N))_k (t), $ <eq-if-ode>
which is a _non-stiff_ ODE that can be advanced with a standard explicit method such as RK4.

=== Per-Step IF-RK4 Algorithm

In practice, to avoid numerical overflow when $lambda_k$ has large real parts (as in the Kuramoto--Sivashinsky equation), we apply the integrating factor _per step_ rather than cumulatively from $t = 0$. Let $E = exp(lambda_k Delta t\/2)$ and $E_2 = E^2 = exp(lambda_k Delta t)$. Then a single IF-RK4 step advancing $hat(bold(u))$ from time $t_n$ to $t_n + Delta t$ is:

$
bold(N)_a &= Delta t dot hat(cal(N))(hat(bold(u))^n), \
bold(N)_b &= Delta t dot hat(cal(N))(E dot.o hat(bold(u))^n + bold(N)_a \/ 2), \
bold(N)_c &= Delta t dot hat(cal(N))(E dot.o hat(bold(u))^n + bold(N)_b \/ 2), \
bold(N)_d &= Delta t dot hat(cal(N))(E_2 dot.o hat(bold(u))^n + E dot.o bold(N)_c), \
hat(bold(u))^(n+1) &= E_2 dot.o hat(bold(u))^n + frac(1, 6)(E_2 dot.o bold(N)_a + 2 E dot.o (bold(N)_b + bold(N)_c) + bold(N)_d),
$ <eq-ifrk4>

where $dot.o$ denotes the Hadamard (componentwise) product. This is the approach used by Trefethen @Trefethen2000 and Kassam and Trefethen @KassamTrefethen2005. The key advantage is that $E$ and $E_2$ involve only a _single time step_ $Delta t$, so their magnitudes remain bounded even when $lambda_k$ has positive real parts.

For more demanding problems, the _exponential time differencing_ (ETDRK4#idx("ETDRK4")) scheme, proposed by Cox and Matthews @CoxMatthews2002 and refined by Kassam and Trefethen @KassamTrefethen2005, offers improved accuracy by incorporating the integrating factor into the RK weights themselves. Recent work by Fu, Shen, and Yang @FuShenYang2024 develops higher-order energy-decreasing ETD Runge--Kutta methods for gradient flows. The IF-RK4 scheme above is simpler and sufficient for all the études in this chapter.

== KdV Equation: Soliton Interactions <sec-kdv>

The Korteweg--de Vries (KdV) equation
$ u_t + 6 u u_x + u_(x x x) = 0 $ <eq-kdv>
is one of the most celebrated nonlinear PDEs. It combines a quadratic nonlinear steepening term with a third-order dispersive term, and admits the remarkable _soliton_ solutions
$ u(x, t) = 2 a^2 op("sech")^2 (a (x - x_0 - 4 a^2 t)), $ <eq-soliton>
where the parameter $a > 0$ controls both the amplitude ($2 a^2$) and the speed ($4 a^2$). Taller solitons travel faster.

=== Fourier Pseudospectral Formulation

Writing @eq-kdv in the form $u_t = cal(L) u + cal(N)(u)$ with $cal(L) u = -u_(x x x)$ and $cal(N)(u) = -6 u u_x$, the Fourier-space system is
$ hat(u)_k '(t) = i k^3 hat(u)_k (t) - hat(6 u u_x)_k (t). $ <eq-kdv-fourier>
The linear symbol $lambda_k = i k^3$ is purely imaginary, so $|e^(lambda_k t)| = 1$ for all $t$; there is no exponential growth. We apply the IF-RK4 scheme from @sec-integrating-factor, evaluating the nonlinear term $-6 u u_x$ pseudospectrally with two-thirds dealiasing.

=== Computational Étude 11.4: KdV Solitons <etude-kdv>

==== Part A: The Zabusky--Kruskal Experiment

In their seminal 1965 paper @ZabuskyKruskal1965, Zabusky and Kruskal studied the KdV equation#idx("KdV equation") (in a slightly different normalisation) with the cosine initial condition $u(x, 0) = cos(pi x)$ on $[0, 2)$. They observed that the initial wave steepens, breaks into a train of solitary pulses, and then (after a recurrence time) reassembles into an approximation of the original cosine. This was the numerical experiment that led to the discovery of _solitons_ and the connection to the Fermi--Pasta--Ulam recurrence phenomenon.

In Python, the core of the integrating factor solver is:

```python
lam = -1j * delta2 * k**3  # linear symbol

def G(t_now, v_hat_now):
    """RHS in integrating factor variables."""
    E = np.exp(lam * t_now)
    u_hat = E * v_hat_now
    u = np.fft.ifft(u_hat).real
    ux = np.fft.ifft(1j * k * u_hat).real
    f_nl_hat = np.fft.fft(-u * ux)
    f_nl_hat[~mask] = 0.0
    return np.exp(-lam * t_now) * f_nl_hat
```

The equivalent MATLAB implementation:

```matlab
lam = -1i * delta2 * k.^3;

function g = G(t_now, v_hat_now)
    E = exp(lam * t_now);
    u_hat = E .* v_hat_now;
    u = real(ifft(u_hat));
    ux = real(ifft(1i * k .* u_hat));
    f_nl_hat = fft(-u .* ux);
    f_nl_hat(~mask) = 0;
    g = exp(-lam * t_now) .* f_nl_hat;
end
```

The Julia implementation:

```julia
lam = -im * delta2 .* k.^3

function G(t_now, v_hat_now, lam, k, mask)
    E = exp.(lam .* t_now)
    u_hat = E .* v_hat_now
    u = real.(ifft(u_hat))
    ux = real.(ifft(im .* k .* u_hat))
    f_nl_hat = fft(-u .* ux)
    f_nl_hat[.!mask] .= 0.0
    return exp.(-lam .* t_now) .* f_nl_hat
end
```

#figure(
  image("../figures/ch11/python/kdv_zabusky_kruskal.pdf", width: 95%),
  caption: [The Zabusky--Kruskal experiment: KdV with initial condition $u(x, 0) = cos(pi x)$ on $[0, 2)$. _Left_: Space-time diagram showing the formation of a soliton train from the steepening cosine wave. _Right_: Snapshots at selected times. The initial cosine steepens, breaks into distinct solitons of different amplitudes and speeds, and exhibits the characteristic structure that led Zabusky and Kruskal to coin the term "soliton."],
) <fig-kdv-zk>

The code generating @fig-kdv-zk is available in:
- `codes/python/ch11/kdv_zabusky_kruskal.py`
- `codes/matlab/ch11/kdv_zabusky_kruskal.m`
- `codes/julia/ch11/kdv_zabusky_kruskal.jl`

==== Part B: Two-Soliton Elastic Collision

A hallmark of KdV integrability is the _elastic interaction_ of solitons. We place two solitons of the form @eq-soliton with parameters $a_1 = 1$ and $a_2 = 0.5$ on a large periodic domain $[0, 80)$. Since $a_1 > a_2$, the taller soliton (amplitude $2$, speed $4$) is faster than the shorter one (amplitude $0.5$, speed $1$). After a collision, both solitons re-emerge with their original shapes, differing only by _phase shifts_.

#figure(
  image("../figures/ch11/python/kdv_snapshots.pdf", width: 95%),
  caption: [Two-soliton collision for the KdV equation @eq-kdv with $a_1 = 1$, $a_2 = 0.5$ on $[0, 80)$. The four panels show the solution before, during, and after the interaction. The taller soliton overtakes the shorter one; after the collision, both re-emerge essentially unchanged, shifted in phase.],
) <fig-kdv-collision>

The code generating @fig-kdv-collision is available in:
- `codes/python/ch11/kdv_soliton.py`
- `codes/matlab/ch11/kdv_soliton.m`
- `codes/julia/ch11/kdv_soliton.jl`

#etude-conclusion[
  The two parts of this étude illustrate complementary facets of KdV. Part A reproduces the *landmark experiment of Zabusky and Kruskal* @ZabuskyKruskal1965: a smooth cosine initial condition steepens under $6 u u_x$ and, just before "breaking", is rescued by the third-order dispersion $u_(x x x)$ which splits the wave into a train of solitary pulses whose amplitudes and speeds are set by the inverse scattering transform. Part B demonstrates the hallmark of KdV *integrability*: two solitons of different amplitudes collide and re-emerge with their shapes perfectly intact, differing only by phase shifts. The numerical conservation of the soliton profiles over long integration times is a stringent test of both the spatial spectral accuracy and the temporal integrating-factor scheme. Fourier pseudospectral methods are ideally suited to dispersive wave equations on periodic domains: the linear symbol $lambda_k = i k^3$ is purely imaginary, avoiding exponential growth and allowing large, stable time steps.
]

== Nonlinear Schrödinger Equation: Recurrence <sec-nls>

The focusing nonlinear Schrödinger (NLS) equation
$ i u_t + u_(x x) + 2 |u|^2 u = 0, quad x in [0, 2 pi), $ <eq-nls>
arises in nonlinear optics, water waves, and Bose--Einstein condensation. It is Hamiltonian and preserves the $L^2$ norm $||u||_2^2 = integral |u|^2 dif x$.

=== Split-Step Fourier Method

For NLS, an elegant alternative to the integrating factor is the _split-step Fourier method_ (Strang splitting#idx("Strang splitting")), introduced by Hardin and Tappert @HardinTappert1973 and systematically studied by Taha and Ablowitz @Taha1984. We decompose the evolution into a linear step and a nonlinear step:

+ *Linear half-step* (in Fourier space): $hat(u)_k arrow.l e^(-i k^2 Delta t\/2) hat(u)_k$.
+ *Nonlinear full step* (in physical space): $u arrow.l e^(2 i |u|^2 Delta t) u$.
+ *Linear half-step*: $hat(u)_k arrow.l e^(-i k^2 Delta t\/2) hat(u)_k$.

Each substep is exact for its respective part, and Strang's symmetric ordering gives second-order accuracy in $Delta t$. Crucially, the split-step method#idx("split-step method") preserves the $L^2$ norm exactly (up to roundoff), making it a _structure-preserving integrator_ well suited to Hamiltonian PDEs. Li and Qin @LiQin2025 recently established uniform error estimates for energy-preserving exponential wave integrators applied to the NLS equation.

=== Computational Étude 11.5: NLS Modulation Instability <etude-nls>

We take the initial condition $u(x, 0) = 1 + 0.1 cos x$, a small perturbation of a plane wave. The focusing nonlinearity amplifies this perturbation through _modulation instability#idx("modulation instability")_, generating a pattern of recurrent modulation cycles.

#figure(
  image("../figures/ch11/python/nls_recurrence.pdf", width: 95%),
  caption: [NLS modulation instability and recurrence with initial condition $u(x, 0) = 1 + 0.1 cos x$. _Left_: Space-time diagram of $|u(x, t)|$ showing periodic modulation cycles driven by the focusing nonlinearity. _Right_: Relative deviation of the $L^2$ norm from its initial value, confirming that the split-step method preserves the norm to machine precision.],
) <fig-nls>

The code generating @fig-nls is available in:
- `codes/python/ch11/nls_recurrence.py`
- `codes/matlab/ch11/nls_recurrence.m`
- `codes/julia/ch11/nls_recurrence.jl`

#etude-conclusion[
  Two phenomena reward numerical exploration here. First, the small cosine perturbation of the plane wave is amplified by the focusing nonlinearity into a sequence of nearly periodic modulation cycles --- the *Fermi--Pasta--Ulam--Tsingou recurrence* in the NLS setting. The modulation depth grows, saturates, and returns close to the initial state before the cycle repeats; this quasi-periodic behaviour is a signature of NLS integrability, and would be extremely difficult to predict quantitatively from perturbation theory alone. Second, the split-step Fourier method *preserves $|| u ||_2$ to machine precision* over the entire integration, because each substep is individually norm-preserving. This structure-preserving property is essential for long-time simulations of Hamiltonian PDEs, where even tiny systematic drift in a conserved quantity will eventually corrupt the dynamics.
]

== Kuramoto--Sivashinsky Equation: Spatiotemporal Chaos <sec-ks>

The Kuramoto--Sivashinsky (KS) equation
$ u_t + u u_x + u_(x x) + u_(x x x x) = 0, quad x in [0, L), $ <eq-ks>
models flame-front instabilities @Kuramoto1978, thin-film flows, and other systems that exhibit spatiotemporal chaos @Sivashinsky1977. The second derivative $u_(x x)$ acts as _negative diffusion_ at low wavenumbers, injecting energy into the system, while the fourth derivative $u_(x x x x)$ provides damping at high wavenumbers. The balance between these two mechanisms, mediated by the quadratic nonlinearity, produces rich chaotic dynamics.

=== Fourier Pseudospectral Formulation

The linear symbol is $lambda_k = k^2 - k^4$: positive for $|k| < 1$ (unstable) and strongly negative for $|k| >> 1$ (stiff). An explicit method would require $Delta t tilde.op N^(-4)$. The per-step IF-RK4 from @eq-ifrk4 @KassamTrefethen2005 handles this stiffness exactly, with the nonlinear term $hat(cal(N))_k = -i k\/2 dot hat(u^2)_k$ (using the conservation form $u u_x = (u^2)_x \/ 2$) evaluated pseudospectrally with two-thirds dealiasing.

=== Computational Étude 11.6: Kuramoto--Sivashinsky Chaos <etude-ks>

We set $L = 32 pi$ (large enough for the chaotic regime) and start from a small random perturbation. @fig-ks-spacetime shows the characteristic space-time "fingerprint" of KS chaos: an initial transient followed by statistically steady cellular dynamics resembling one-dimensional turbulence.

In Python, the core per-step IF-RK4 loop is:

```python
E = np.exp(lam * dt / 2)  # half-step propagator
E2 = E**2                  # full-step propagator

def Nhat(u_hat):
    u = np.fft.ifft(u_hat).real
    f_hat = -0.5j * k * np.fft.fft(u**2)
    f_hat[~mask] = 0.0
    return f_hat

for n in range(n_steps):
    Na = dt * Nhat(u_hat)
    Nb = dt * Nhat(E * u_hat + Na / 2)
    Nc = dt * Nhat(E * u_hat + Nb / 2)
    Nd = dt * Nhat(E2 * u_hat + E * Nc)
    u_hat = E2 * u_hat + (E2*Na + 2*E*(Nb+Nc) + Nd) / 6
```

The equivalent MATLAB implementation:

```matlab
E = exp(lam * dt / 2);
E2 = E.^2;

for n = 1:nsteps
    Na = dt * Nhat(u_hat);
    Nb = dt * Nhat(E .* u_hat + Na/2);
    Nc = dt * Nhat(E .* u_hat + Nb/2);
    Nd = dt * Nhat(E2 .* u_hat + E .* Nc);
    u_hat = E2 .* u_hat + (E2.*Na + 2*E.*(Nb+Nc) + Nd)/6;
end
```

The Julia implementation:

```julia
E = exp.(lam .* dt / 2)
E2 = E.^2

for n in 1:n_steps
    Na = dt .* Nhat(u_hat)
    Nb = dt .* Nhat(E .* u_hat .+ Na ./ 2)
    Nc = dt .* Nhat(E .* u_hat .+ Nb ./ 2)
    Nd = dt .* Nhat(E2 .* u_hat .+ E .* Nc)
    u_hat .= E2 .* u_hat .+ (E2 .* Na .+ 2 .* E .* (Nb .+ Nc) .+ Nd) ./ 6
end
```

#figure(
  image("../figures/ch11/python/ks_spacetime.pdf", width: 85%),
  caption: [Kuramoto--Sivashinsky equation @eq-ks with $L = 32 pi$: the space-time "fingerprint" of spatiotemporal chaos. The initial random perturbation develops into a pattern of merging and splitting cellular structures, a hallmark of the KS dynamics. The per-step IF-RK4 handles the severe stiffness ($lambda_k tilde.op k^4$) without difficulty.],
) <fig-ks-spacetime>

#figure(
  image("../figures/ch11/python/ks_spectrum.pdf", width: 70%),
  caption: [Time-averaged energy spectrum of the KS equation in the statistically steady chaotic regime. Energy is concentrated at low wavenumbers (the linearly unstable range $|k| < 1$) and decays rapidly at high wavenumbers due to the fourth-derivative damping.],
) <fig-ks-spectrum>

The code generating @fig-ks-spacetime and @fig-ks-spectrum is available in:
- `codes/python/ch11/kuramoto_sivashinsky.py`
- `codes/matlab/ch11/kuramoto_sivashinsky.m`
- `codes/julia/ch11/kuramoto_sivashinsky.jl`

#etude-conclusion[
  Kuramoto--Sivashinsky is the first genuinely *stiff* problem in this chapter, and it showcases why the *integrating-factor approach is indispensable*. The linear symbol $lambda_k = k^2 - k^4$ is positive for $|k| < 1$ (energy injection) and grows as $-k^4$ for large $|k|$ (stiff damping). A naive explicit method would require $Delta t tilde.op N^(-4)$, prohibitive at large $N$; per-step IF-RK4 absorbs the stiffness into exponential propagators, allowing a time step governed only by the nonlinear CFL condition. The space-time figure displays the resulting *spatiotemporal chaos* --- cellular structures that merge, split, and interact indefinitely. The spectrum confirms that energy is concentrated in the linearly unstable band $|k| < 1$ and decays sharply at higher wavenumbers. KS is a *miniature model of turbulence*: chaotic, yet low-dimensional enough that a modest Fourier pseudospectral code captures its full dynamics.
]

== 2D Navier--Stokes: Vortex Dynamics and Decaying Turbulence <sec-ns2d>

As the capstone of this chapter, we simulate two-dimensional incompressible flow on the doubly periodic domain $[0, 2 pi)^2$. This is the problem that originally demonstrated the superiority of spectral methods in the 1970s @Orszag1971, and the formulation is a natural extension of the one-dimensional techniques developed above.

=== Vorticity--Streamfunction Formulation

For two-dimensional incompressible flow, the velocity field $bold(u) = (u, v)$ satisfies $nabla dot bold(u) = 0$, which allows us to introduce a streamfunction $psi$ such that $u = psi_y$ and $v = -psi_x$. The scalar vorticity $omega = v_x - u_y$ then satisfies
$ omega_t + J(psi, omega) = nu Delta omega, quad Delta psi = -omega, $ <eq-ns2d>
where $J(psi, omega) = psi_x omega_y - psi_y omega_x$ is the Jacobian (advection by the velocity field) and $nu$ is the kinematic viscosity.

=== Spectral Implementation

In Fourier space, the key operations become remarkably simple:

- *Poisson solve*: $hat(psi)_(k_x, k_y) = hat(omega)_(k_x, k_y) \/ (k_x^2 + k_y^2)$ for $(k_x, k_y) eq.not (0, 0)$. This is _trivial division_, the spectral analogue of a linear system solve.
- *Velocity*: $hat(u) = i k_y hat(psi)$, $hat(v) = -i k_x hat(psi)$.
- *Jacobian*: compute $u$, $v$, $omega_x$, $omega_y$ in physical space via inverse FFTs, then form $J = u omega_x + v omega_y$ pointwise and transform back. Apply two-thirds dealiasing in both directions.

In Python, the right-hand side is:

```python
def rhs(omega):
    omega_hat = np.fft.fft2(omega)
    psi_hat = omega_hat / K2   # Poisson solve
    u  = np.fft.ifft2( 1j * KY * psi_hat).real
    v  = np.fft.ifft2(-1j * KX * psi_hat).real
    ox = np.fft.ifft2( 1j * KX * omega_hat).real
    oy = np.fft.ifft2( 1j * KY * omega_hat).real
    jac_hat = np.fft.fft2(u * ox + v * oy)
    jac_hat[~mask] = 0.0
    return np.fft.ifft2(-jac_hat - nu * K2 * omega_hat).real
```

The equivalent MATLAB implementation:

```matlab
function f = rhs(omega)
    omega_hat = fft2(omega);
    psi_hat = omega_hat ./ K2;
    u  = real(ifft2( 1i * KY .* psi_hat));
    v  = real(ifft2(-1i * KX .* psi_hat));
    ox = real(ifft2( 1i * KX .* omega_hat));
    oy = real(ifft2( 1i * KY .* omega_hat));
    jac_hat = fft2(u .* ox + v .* oy);
    jac_hat(~mask) = 0;
    f = real(ifft2(-jac_hat - nu * K2 .* omega_hat));
end
```

The Julia implementation:

```julia
function rhs(omega, KX, KY, K2, nu, mask)
    omega_hat = fft(omega)
    psi_hat = omega_hat ./ K2
    u  = real.(ifft( im .* KY .* psi_hat))
    v  = real.(ifft(-im .* KX .* psi_hat))
    ox = real.(ifft( im .* KX .* omega_hat))
    oy = real.(ifft( im .* KY .* omega_hat))
    jac_hat = fft(u .* ox .+ v .* oy)
    jac_hat[.!mask] .= 0.0
    return real.(ifft(-jac_hat .- nu .* K2 .* omega_hat))
end
```

=== Computational Étude 11.7: 2D Decaying Turbulence <etude-ns2d>

We initialise the vorticity with random Fourier modes in the band $1 < |bold(k)| < 10$ and set $nu = 10^(-3)$. @fig-ns2d-vorticity shows the vorticity field at four times: the initially random small-scale vortices merge into progressively larger coherent structures, the hallmark of the _inverse energy cascade_ in two-dimensional turbulence.

#figure(
  image("../figures/ch11/python/ns2d_vorticity.pdf", width: 95%),
  caption: [2D Navier--Stokes: decaying turbulence on $[0, 2 pi)^2$ with $nu = 10^(-3)$ and $N = 128$. Vorticity snapshots at $t = 0$, $5$, $10$, and $20$. Small-scale structures merge into large coherent vortices, demonstrating the inverse energy cascade characteristic of two-dimensional turbulence.],
) <fig-ns2d-vorticity>

#figure(
  image("../figures/ch11/python/ns2d_energy.pdf", width: 95%),
  caption: [Diagnostics for 2D decaying turbulence. _Left_: Kinetic energy $E(t) = frac(1, 2) chevron.l |bold(u)|^2 chevron.r$ decays slowly as energy cascades to large scales and is eventually dissipated. _Right_: Enstrophy $Z(t) = frac(1, 2) chevron.l omega^2 chevron.r$ decays much faster, consistent with the selective decay principle of two-dimensional turbulence.],
) <fig-ns2d-energy>

The code generating @fig-ns2d-vorticity and @fig-ns2d-energy is available in:
- `codes/python/ch11/navier_stokes_2d.py`
- `codes/matlab/ch11/navier_stokes_2d.m`
- `codes/julia/ch11/navier_stokes_2d.jl`

#etude-conclusion[
  The vorticity snapshots illustrate the *inverse energy cascade*: small vortices from the random initial condition merge into progressively larger coherent structures, transferring kinetic energy to ever-larger scales. Energy $E(t)$ decays slowly (it would be conserved in the inviscid limit) while enstrophy $Z(t)$ decays much faster, consistent with the *selective decay principle* of 2D turbulence @Kraichnan1967. The entire 2D Navier--Stokes solver fits in roughly 50 lines of code per language --- one of the most compelling arguments for Fourier pseudospectral methods. Production-level turbulence simulations (weather, climate, astrophysical fluids) are built on precisely the same algorithmic skeleton. The Poisson solve, the most expensive step in an FD/FE code, reduces here to a *trivial division in Fourier space*, and the 2D dealiasing mask extends naturally from its 1D counterpart.
]

== Three Further Applications

The études above illustrate the main algorithmic ideas: Fourier pseudospectral differentiation, dealiasing, integrating factors, and split-step methods. We close the computational part of this chapter with three additional applications that broaden the palette of PDEs and time-stepping strategies.

=== Computational Étude 11.8: Variable-Coefficient Transport <etude-transport-variable>

In @etude-advection we solved constant-coefficient advection $u_t + c u_x = 0$ with $c = 1$. Real wave propagation, however, often occurs in media with spatially varying properties. Consider
$ u_t + c(x) u_x = 0, quad 0 < x < 2 pi, quad u "periodic", $ <eq-transport-variable>
with
$ c(x) = 0.3 + sin^2(x - 1). $

The same Fourier differentiation and RK4 time stepping apply without modification; only the right-hand side changes. In Python:

```python
def fourier_diff(v):
    """Derivative of a periodic function via FFT."""
    N = len(v)
    v_hat = np.fft.fft(v)
    k = np.concatenate([np.arange(N//2), [0], np.arange(1-N//2, 0)])
    return np.fft.ifft(1j * k * v_hat).real

# RHS and RK4 time stepping
rhs = lambda u: -c * fourier_diff(u)
for n in range(nsteps):
    k1 = rhs(u)
    k2 = rhs(u + 0.5*dt*k1)
    k3 = rhs(u + 0.5*dt*k2)
    k4 = rhs(u + dt*k3)
    u  = u + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
```

The equivalent MATLAB implementation:

```matlab
function w = fourier_diff(v)
% Derivative of a periodic function via FFT.
    N = length(v);
    v_hat = fft(v);
    k = [0:N/2-1, 0, -N/2+1:-1]';
    w = real(ifft(1i * k .* v_hat));
end

% RHS and RK4 time stepping
rhs = @(u) -c .* fourier_diff(u);
for n = 1:nsteps
    k1 = rhs(u);
    k2 = rhs(u + 0.5*dt*k1);
    k3 = rhs(u + 0.5*dt*k2);
    k4 = rhs(u + dt*k3);
    u  = u + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end
```

The Julia implementation:

```julia
function fourier_diff(v)
    N = length(v)
    v_hat = fft(v)
    k = [collect(0:N÷2-1); 0; collect(1-N÷2:-1)]
    real.(ifft(1im .* k .* v_hat))
end

# RHS and RK4 time stepping
rhs(u) = -c .* fourier_diff(u)
for n in 1:nsteps
    k1 = rhs(u)
    k2 = rhs(u .+ 0.5 .* dt .* k1)
    k3 = rhs(u .+ 0.5 .* dt .* k2)
    k4 = rhs(u .+ dt .* k3)
    u  = u .+ (dt / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
end
```

#figure(
  image("../figures/ch11/python/transport_variable.pdf", width: 85%),
  caption: [Solution of the variable-coefficient transport equation @eq-transport-variable. Left: space-time diagram showing the characteristic trajectory of a Gaussian pulse in a medium where $c(x) = 0.3 + sin^2(x - 1)$. The pulse accelerates in fast regions and decelerates in slow regions. Right: waterfall view of the same solution with vertical offset proportional to time.],
) <fig-transport-variable>

The code generating @fig-transport-variable is available in:
- `codes/python/ch11/transport_variable.py`
- `codes/matlab/ch11/transport_variable.m`
- `codes/julia/ch11/transport_variable.jl`

#etude-conclusion[
  The Gaussian pulse accelerates through the fast region of $c(x) = 0.3 + sin^2(x - 1)$ and decelerates through the slow region, tracing *curved characteristics* in the space-time plane. The Fourier pseudospectral derivative requires *no modification* for the variable coefficient: the product $c(x) u_x$ is formed pointwise in physical space, exactly as in the Burgers nonlinearity of @etude-burgers. The pseudospectral approach handles variable-coefficient linear problems and nonlinear problems with equal ease --- the only requirement is that the solution remain smooth and periodic. Because $c(x) > 0$ everywhere, the equation is purely advective with no diffusion, and the pulse maintains its shape (merely translating at variable speed) to the accuracy of the spectral discretisation; any visible distortion would be a numerical artefact, making this a sensitive solver test.
]

=== Computational Étude 11.9: Schrödinger Equation <etude-schrodinger>

The time-dependent Schrödinger equation with a harmonic potential is
$ i u_t = -u_(x x) + x^2 u, $ <eq-schrodinger>
where $u(x, t)$ is the complex wavefunction and $x^2$ is the harmonic oscillator potential.

A natural approach is _Strang splitting_, an operator splitting strategy introduced by Strang @Strang1968 that has become the method of choice for dispersive equations. Taha and Ablowitz @Taha1984 established the split-step Fourier method as the robust standard for the nonlinear Schrödinger equation#idx("nonlinear Schrödinger equation"), while Bao, Jin, and Markowich @BaoJinMarkowich2002 provided rigorous error analysis for time-splitting spectral approximations in the semiclassical regime. The algorithm proceeds in three sub-steps:
1. *Half potential step* (in physical space): $u arrow u dot e^(-i x^2 Delta t \/ 2)$
2. *Full kinetic step* (in Fourier space): $hat(u)_k arrow hat(u)_k dot e^(-i k^2 Delta t)$
3. *Half potential step*: $u arrow u dot e^(-i x^2 Delta t \/ 2)$

This splitting is second-order accurate in time and preserves the $L^2$ norm of the wavefunction (unitarity). The kinetic step uses exactly the same FFT machinery as the split-step method of @etude-nls, but here the nonlinear cubic is replaced by an external potential.

The Strang splitting loop in Python:

```python
# Precompute splitting operators
half_potential = np.exp(-0.5j * V * dt)
full_kinetic  = np.exp(-1j * k**2 * dt)

for n in range(nsteps):
    u     = half_potential * u              # half potential
    u_hat = full_kinetic * np.fft.fft(u)    # full kinetic
    u     = half_potential * np.fft.ifft(u_hat)  # half potential
```

The equivalent MATLAB implementation:

```matlab
% Precompute splitting operators
half_potential = exp(-0.5i * V * dt);
full_kinetic   = exp(-1i * k.^2 * dt);

for n = 1:nsteps
    u     = half_potential .* u;             % half potential
    u_hat = full_kinetic .* fft(u);          % full kinetic
    u     = half_potential .* ifft(u_hat);   % half potential
end
```

The Julia implementation:

```julia
# Precompute splitting operators
half_potential = exp.(-0.5im .* V .* dt)
full_kinetic   = exp.(-1im .* k.^2 .* dt)

for n in 1:nsteps
    u     .= half_potential .* u              # half potential
    u_hat  = full_kinetic .* fft(u)           # full kinetic
    u     .= half_potential .* ifft(u_hat)    # half potential
end
```

#figure(
  image("../figures/ch11/python/schrodinger.pdf", width: 85%),
  caption: [Solution of the Schrödinger equation @eq-schrodinger with harmonic potential using Strang splitting. Left: space-time evolution of the probability density $|u(x, t)|^2$ showing the oscillatory dynamics of a Gaussian wavepacket in the harmonic trap. Right: waterfall plot of $|u(x, t)|^2$ at selected times.],
) <fig-schrodinger>

The code generating @fig-schrodinger is available in:
- `codes/python/ch11/schrodinger.py`
- `codes/matlab/ch11/schrodinger.m`
- `codes/julia/ch11/schrodinger.jl`

#etude-conclusion[
  The figure shows a Gaussian wavepacket oscillating in the harmonic potential $V(x) = x^2$: the packet *breathes and translates periodically*, exactly as predicted by the analytical theory of coherent states. The Strang splitting used here is the kin of the split-step method in @etude-nls, with the external potential $x^2 u$ replacing the cubic nonlinearity, and with the splitting order reversed (potential halved here, kinetic halved there). Both substeps --- potential kick in physical space, kinetic propagation in Fourier space --- are *exact and individually unitary*, so the composite scheme preserves $||u||_2$ to machine precision. In quantum mechanics this is not aesthetic: $|| u ||_2$ is total probability, and any scheme that fails to conserve it will eventually produce unphysical results.
]

=== Computational Étude 11.10: Allen--Cahn Reaction-Diffusion <etude-allen-cahn>

The études above treated stiffness with integrating factors (@etude-kdv, @etude-ks) or split-step methods (@etude-nls, @etude-schrodinger). An alternative for reaction-diffusion equations is the _IMEX_ (implicit-explicit) strategy, which treats the stiff linear part implicitly and the nonlinear part explicitly.

The Allen--Cahn equation combines diffusion with a bistable nonlinearity:
$ u_t = epsilon^2 u_(x x) + u(1 - u^2), $ <eq-allen-cahn>
where the parameter $epsilon > 0$ controls the width of the transition layers between the two stable phases $u = plus.minus 1$. The nonlinear term $u(1 - u^2)$ has stable equilibria at $u = plus.minus 1$ and an unstable equilibrium at $u = 0$. Solutions tend to separate into regions where $u approx 1$ and $u approx -1$, connected by thin fronts of width $O(epsilon)$, which then slowly annihilate (coarsening). The Allen--Cahn equation, originally formulated by Allen and Cahn @AllenCahn1979 to model antiphase boundary motion in crystalline solids, demands numerical schemes that respect the thermodynamic structure of the system. Shen and Yang @ShenYang2010 introduced unconditionally energy-stable spectral schemes that preserve the dissipation of the Ginzburg--Landau free energy, building on the convex splitting framework of Eyre @Eyre1998.

The IMEX scheme treats the stiff linear diffusion implicitly and the nonlinear reaction explicitly:
$ (I - epsilon^2 Delta t D_(2,i)) bold(u)^(n+1) = bold(u)^n + Delta t (bold(u)^n - (bold(u)^n)^3), $
where the cubic is applied componentwise. In Fourier space, the implicit diffusion becomes a diagonal division, so no linear system solve is needed.

The IMEX time loop in Python:

```python
# IMEX denominator: 1 / (1 + eps^2 * k^2 * dt)
imex_denom = 1.0 + eps**2 * k**2 * dt

for n in range(nsteps):
    reaction_hat = np.fft.fft(u - u**3)
    u_hat = (np.fft.fft(u) + dt * reaction_hat) / imex_denom
    u = np.fft.ifft(u_hat).real
```

The equivalent MATLAB implementation:

```matlab
% IMEX denominator: 1 / (1 + eps^2 * k^2 * dt)
imex_denom = 1.0 + eps^2 * k.^2 * dt;

for n = 1:nsteps
    reaction_hat = fft(u - u.^3);
    u_hat = (fft(u) + dt * reaction_hat) ./ imex_denom;
    u = real(ifft(u_hat));
end
```

The Julia implementation:

```julia
# IMEX denominator: 1 / (1 + eps^2 * k^2 * dt)
imex_denom = 1.0 .+ eps^2 .* k.^2 .* dt

for n in 1:nsteps
    reaction_hat = fft(u .- u.^3)
    u_hat = (fft(u) .+ dt .* reaction_hat) ./ imex_denom
    u = real.(ifft(u_hat))
end
```

#figure(
  image("../figures/ch11/python/allen_cahn.pdf", width: 85%),
  caption: [Solution of the Allen--Cahn equation @eq-allen-cahn with $epsilon = 0.05$ using IMEX time stepping. Left: space-time evolution showing phase separation into regions where $u approx plus.minus 1$, followed by slow coarsening as fronts annihilate. Right: snapshots at selected times showing the initial perturbation, phase separation, and progressive front reduction.],
) <fig-allen-cahn>

The code generating @fig-allen-cahn is available in:
- `codes/python/ch11/allen_cahn.py`
- `codes/matlab/ch11/allen_cahn.m`
- `codes/julia/ch11/allen_cahn.jl`

#etude-conclusion[
  The figure illustrates the two-stage dynamics of Allen--Cahn: a rapid *phase separation* into regions where $u approx plus.minus 1$ separated by thin $cal(O)(epsilon)$ fronts, followed by a much slower *coarsening* in which neighbouring fronts annihilate in pairs. The *IMEX time stepper* treats the stiff diffusion $epsilon^2 u_(x x)$ implicitly and the bistable reaction $u(1 - u^2)$ explicitly. In Fourier space the implicit step reduces to a diagonal division --- *no linear system solve* --- so the scheme is unconditionally stable with respect to diffusion and trivial to implement. This is the third distinct time-stepping strategy in the chapter (after integrating factors for KdV/KS and split-step for NLS/Schrödinger), and the comparison highlights a general principle: *the best integrator depends on the structure of the PDE, not merely its order*. IMEX is simpler and more efficient than integrating factors when the linear part is purely dissipative.
]

== A Non-Exhaustive Literature Overview

The Fourier pseudospectral approach to periodic PDEs has a rich history that parallels the development of modern scientific computing. The efficient computation of discrete Fourier transforms was revolutionised by the fast Fourier transform (FFT) algorithm of Cooley and Tukey @Cooley1965, which reduced the cost from $O(N^2)$ to $O(N log N)$ and made spectral methods computationally viable. Orszag @Orszag1971 pioneered the pseudospectral treatment of fluid dynamics in the early 1970s, and the technique rapidly became the dominant method for homogeneous turbulence simulations. Kreiss and Oliger @KreissOliger1972 established that spectral methods require far fewer points per wavelength than finite differences ($pi$ versus 10--20), a result that remains central to the practical appeal of spectral computation. Fornberg @Fornberg1987 provided a comprehensive comparison of pseudospectral methods with finite differences, further clarifying the eigenvalue structure of circulant differentiation matrices.

The treatment of aliasing in nonlinear terms was addressed early by Orszag @Orszag1971 (the two-thirds rule) and Patterson and Orszag @PattersonOrszag1971, who showed that dealiasing is essential for energy-conserving simulations of turbulence. Tadmor @Tadmor1989 developed the spectral viscosity approach for conservation laws, providing a rigorous framework for convergence of spectral approximations to nonlinear problems with shocks. The connection between aliasing and energy conservation remains an active area of research in large-eddy simulation and climate modelling.

The theoretical foundations of spectral methods were laid out comprehensively by Gottlieb and Orszag @GottliebOrszag1977 in their influential monograph, and the stability analysis for spectral discretisations was developed by Gottlieb and Turkel @GottliebTurkel1980. The method of lines framework, which we use throughout this chapter, was systematically treated by Schiesser @Schiesser1991. For general references on Chebyshev and Fourier spectral methods, the reader may consult Fornberg @Fornberg1996, Trefethen @Trefethen2000, and Canuto _et al._ @Canuto2006, each of which covers the periodic pseudospectral framework from a complementary perspective.

For stiff problems, the integrating factor method dates back to Lawson @Lawson1967, who proposed treating the linear part of the ODE exactly via an exponential change of variables. The exponential time differencing approach was developed into a practical fourth-order Runge--Kutta scheme (ETDRK4) by Cox and Matthews @CoxMatthews2002, and Kassam and Trefethen @KassamTrefethen2005 identified and corrected numerical instabilities in the original formulation, producing a robust algorithm that is now standard for stiff pseudospectral PDE solvers. More recently, Fu, Shen, and Yang @FuShenYang2024 have developed higher-order energy-decreasing ETD Runge--Kutta methods for gradient flows, extending the integrating factor framework to structure-preserving settings. The split-step Fourier method for the NLS equation was introduced by Hardin and Tappert @HardinTappert1973, systematically studied by Taha and Ablowitz @Taha1984, and its symplectic structure was analysed in the broader context of geometric numerical integration. Li and Qin @LiQin2025 have recently established uniform error estimates for energy-preserving exponential wave integrators applied to the NLS equation.

The KdV equation and its soliton solutions have a fascinating history. The numerical discovery of solitons by Zabusky and Kruskal @ZabuskyKruskal1965, building on the earlier Fermi--Pasta--Ulam numerical experiments, led to the development of inverse scattering theory and the modern theory of integrable systems. Fornberg @Fornberg1996 provides an excellent account of the role of pseudospectral methods in simulating nonlinear wave phenomena. The Kuramoto--Sivashinsky equation, studied by Kuramoto @Kuramoto1978 and Sivashinsky @Sivashinsky1977 in the context of chemical waves and flame fronts, has become a canonical model for spatiotemporal chaos and a standard benchmark for stiff PDE solvers.

For two-dimensional turbulence, the inverse energy cascade was predicted theoretically by Kraichnan @Kraichnan1967 and confirmed by early pseudospectral simulations. Modern DNS (direct numerical simulation) of turbulence relies heavily on the algorithmic framework developed in this chapter, scaled to massively parallel architectures. Boffetta _et al._ @Boffetta2025 have recently used GPU-accelerated pseudospectral codes to study the direct enstrophy cascade at unprecedented resolutions. Canuto _et al._ @Canuto2006 provide a comprehensive treatment of spectral methods for fluid dynamics, including the vorticity--streamfunction formulation used here.

Operator splitting for Schrödinger-type equations was introduced by Strang @Strang1968; the resulting split-step Fourier method preserves unitarity and extends naturally from linear to nonlinear (NLS) settings. Bao, Jin, and Markowich @BaoJinMarkowich2002 provided rigorous error analysis for time-splitting spectral approximations in the semiclassical regime. For reaction-diffusion equations such as the Allen--Cahn equation @AllenCahn1979, implicit-explicit (IMEX) schemes treat diffusion implicitly in Fourier space (a diagonal operation) while keeping the nonlinear reaction explicit. Shen and Yang @ShenYang2010 developed unconditionally energy-stable spectral IMEX schemes, building on the convex splitting framework of Eyre @Eyre1998.

== Summary <sec-ch11-summary>

This chapter has developed Fourier pseudospectral methods for nonlinear time-dependent PDEs on periodic domains. The main themes are:

+ *Fourier pseudospectral derivatives* are computed via FFT in $O(N log N)$ operations, achieving the same spectral accuracy as the dense $O(N^2)$ matrix approach.

+ *Periodic grids are friendlier than Chebyshev grids* for explicit time stepping: the CFL constraint scales as $Delta t tilde.op 1\/N$ for first-order operators (versus $Delta t tilde.op 1\/N^2$ on Chebyshev grids) and $Delta t tilde.op 1\/N^2$ for second-order operators (versus $Delta t tilde.op 1\/N^4$).

+ *Integrating factor methods* eliminate the stiffness of high-order linear operators by treating them exactly via exponential propagators, leaving only the (non-stiff) nonlinear part for explicit time stepping.

+ *Dealiasing* (the two-thirds rule) prevents the corruption of resolved modes by aliased nonlinear interactions and is essential for reliable long-time simulations.

+ *Compact, powerful solvers*: with the techniques developed here, one can simulate:
  - Nondispersive transport (@etude-advection)
  - Viscous Burgers shock formation (@etude-burgers)
  - KdV soliton formation and elastic collisions (@etude-kdv)
  - NLS modulation instability and recurrence (@etude-nls)
  - Kuramoto--Sivashinsky spatiotemporal chaos (@etude-ks)
  - 2D Navier--Stokes decaying turbulence (@etude-ns2d)
  - Variable-coefficient transport (@etude-transport-variable)
  - Schrödinger wavepacket dynamics (@etude-schrodinger)
  - Allen--Cahn phase separation (@etude-allen-cahn)

  Each solver requires fewer than 100 lines of code per language.

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
        table.cell(fill: rgb("#142D6E").lighten(85%))[*PDE*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Key Technique*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Time Stepper*],
      ),
      table.hline(stroke: 0.5pt + luma(180)),
      [11.1], [Fourier differentiation], [Matrix vs FFT], [--],
      [11.2], [$u_t + c u_x = 0$], [Fourier PS], [RK4],
      [11.3], [Burgers: $u_t + u u_x = nu u_(x x)$], [Dealiasing (2/3 rule)], [RK4],
      [11.4], [KdV: $u_t + 6 u u_x + u_(x x x) = 0$], [Integrating factor], [IF-RK4],
      [11.5], [NLS: $i u_t + u_(x x) + 2|u|^2 u = 0$], [Split-step Fourier], [Strang],
      [11.6], [KS: $u_t + u u_x + u_(x x) + u_(x x x x) = 0$], [Per-step IF + dealiasing], [IF-RK4],
      [11.7], [2D NS: $omega_t + J(psi, omega) = nu Delta omega$], [2D FFT + Poisson solve], [RK4],
      [11.8], [$u_t + c(x) u_x = 0$], [Fourier PS (variable coeff.)], [RK4],
      [11.9], [Schrödinger: $i u_t = -u_(x x) + x^2 u$], [Split-step Fourier], [Strang],
      [11.10], [Allen--Cahn: $u_t = epsilon^2 u_(x x) + u(1-u^2)$], [IMEX (Fourier space)], [Euler],
    ),
  ),
  caption: [Summary of computational études in this chapter.],
) <tbl-ch11-summary>

== Exercises <sec-fourier-pseudospectral-exercises>

*Exercise 11.1* (_KdV Soliton Propagation_). Implement a Fourier pseudospectral solver for the KdV equation $u_t + 6 u u_x + u_(x x x) = 0$ on a periodic domain $[0, L]$ with $L = 40$. Use the initial condition $u(x, 0) = 2 "sech"^2(x - 10)$ (a single soliton of speed $c = 4$). (a) Solve using the integrating factor method with RK4. (b) Verify that the soliton propagates at speed $c = 4$ by plotting $u(x, t)$ at $t = 0, 2, 4, 6$. (c) Compute the conserved quantities $I_1 = integral u dif x$ and $I_2 = integral u^2 dif x$ at each time step and verify conservation to machine precision. (d) Repeat with two solitons ($c_1 = 16$ and $c_2 = 4$) and observe elastic collision.

*Exercise 11.2* (_Allen--Cahn Equation with IMEX Splitting_). Solve the Allen--Cahn equation $u_t = epsilon^2 u_(x x) + u(1 - u^2)$ on $[0, 2pi]$ with $epsilon = 0.1$ and random initial data $u(x, 0) = 0.1 xi(x)$ where $xi$ consists of $N\/4$ random Fourier modes. (a) Implement the IMEX scheme: treat the diffusion term implicitly in Fourier space and the nonlinear reaction term explicitly. (b) Run with $N = 256$ and $Delta t = 0.01$ until $t = 50$. (c) Plot the solution at several times and describe the coarsening dynamics. (d) Measure the discrete energy $E = integral (epsilon^2\/2 |u_x|^2 + (1 - u^2)^2\/4) dif x$ and verify that it decreases monotonically.

*Exercise 11.3* (_Energy Conservation in a Pseudospectral Schrödinger Solver_). Solve the cubic nonlinear Schrödinger equation $i u_t + u_(x x) + 2|u|^2 u = 0$ on $[0, 2pi]$ with initial condition $u(x, 0) = "sech"(x - pi) e^(2 i x)$ using the split-step Fourier method. (a) Run with $N = 256$ and $Delta t = 10^(-3)$ to $t = 10$. (b) Monitor three conserved quantities: mass $M = integral |u|^2 dif x$, momentum $P = "Im" integral overline(u) u_x dif x$, and energy $H = integral (|u_x|^2 - |u|^4) dif x$. (c) Plot the relative error in each invariant versus time. (d) Compare Strang splitting ($Delta t$ and $Delta t\/2$ half-steps) with Lie splitting (full step of each operator) and show that Strang splitting preserves the invariants more accurately.

*Exercise 11.4* (_Dealiasing Effect on Long-Time Burgers Simulation_). Solve the viscous Burgers equation $u_t + u u_x = nu u_(x x)$ on $[0, 2pi]$ with $nu = 0.01$ and $u(x, 0) = sin(x)$ using the Fourier pseudospectral method with RK4 time stepping. (a) Run with $N = 64$ and $Delta t = 10^(-3)$ to $t = 5$, both with and without the $2\/3$ dealiasing rule. (b) Plot the Fourier spectrum $|hat(u)_k|$ at $t = 1, 3, 5$ for both cases. (c) Identify the aliased energy that contaminates the low modes in the undealiased case. (d) Compute $||u_(Delta t\/2) - u_(Delta t)||_infinity$ at $t = 5$ for both dealiased and undealiased cases to assess the practical impact on accuracy.
