// textbook/chapters/spectral_pde_solvers.typ
// Chapter 10: Spectral PDE Solvers with Chebyshev and Fourier Grids
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Email: denys.dutykh@ku.ac.ae
// Homepage: https://www.denys-dutykh.com/
// Last modified: February 2026

#import "../styles/template.typ": dropcap, num, format-table

// Enable equation numbering for this chapter
#set math.equation(numbering: "(1)")

= Spectral PDE Solvers with Chebyshev and Fourier Grids <ch-spectral-pde>

#dropcap[The machinery of spectral differentiation, developed carefully over the preceding chapters, now stands ready for its primary purpose: solving partial differential equations. In this chapter, we bring together Chebyshev and Fourier methods to tackle the three fundamental classes of PDEs: hyperbolic equations (waves), parabolic equations (diffusion), and elliptic equations (equilibrium states). The unifying theme is the _method of lines_: discretise space with spectral accuracy, then advance in time (or solve directly for steady states) using appropriate numerical schemes.]

By the end of this chapter, you should be able to:

1. Build Chebyshev differentiation operators via FFT and apply them in 1D and 2D.
2. Discretise in space and time the main PDE classes: wave equations, heat equations, and elliptic problems (Laplace, Poisson, Helmholtz).
3. Understand CFL-type stability limits for explicit schemes, and why implicit schemes are often preferred for parabolic and elliptic problems.
4. Read space--time plots and contour plots, and relate numerical behaviour to physical intuition.
5. Implement spectral PDE solvers for vibrating strings, drum membranes, heat diffusion, and electrostatics.

== Chebyshev--Fourier Recap <sec-cheb-fourier-recap>

Before diving into PDEs, let us briefly recall the geometric connection between Chebyshev and Fourier methods that makes FFT-accelerated Chebyshev differentiation possible. This connection, developed in @ch-chebyshev and @ch-fourier-grids, is the foundation for everything that follows.

=== The Circle-to-Interval Map

The Chebyshev--Gauss--Lobatto points on $[-1, 1]$ are
$ x_j = cos(j pi / N), quad j = 0, 1, dots, N. $ <eq-cheb-points>
These $N + 1$ points arise naturally as projections: if we place $N + 1$ equally-spaced angles
$ theta_j = j pi / N, quad j = 0, 1, dots, N, $
on the upper semicircle and project them onto the $x$-axis, we obtain exactly the Chebyshev points. The map connecting the angle and physical variables is
$ x = cos theta. $ <eq-x-cos-theta>

This geometric interpretation explains the characteristic clustering of Chebyshev points near the boundaries $x = plus.minus 1$: points near $theta = 0$ and $theta = pi$ project to $x approx plus.minus 1$, but the horizontal projection compresses them together. Near the centre ($theta approx pi / 2$), the points spread out more evenly.

=== Chebyshev Polynomials as Wrapped Cosines

The Chebyshev polynomial of degree $n$ is defined by the remarkable formula
$ T_n (x) = cos(n theta), quad "where" x = cos theta. $ <eq-cheb-poly-cos>
This identity reveals that Chebyshev polynomials are simply cosine waves "in disguise": the $n$th Chebyshev polynomial oscillates exactly as $cos(n theta)$ would, but viewed through the nonlinear lens of $x = cos theta$.

The practical consequence is profound: data $v_j = f(x_j)$ sampled at Chebyshev points corresponds to an _even, $2 pi$-periodic function_ in the $theta$ variable:
$ F(theta) = f(cos theta). $
Since $F(theta) = F(-theta)$ (evenness) and $F(theta + 2 pi) = F(theta)$ (periodicity), we can apply Fourier methods, specifically the FFT, to work with Chebyshev expansions efficiently.

#figure(
  image("../figures/ch10/python/cheb_fourier_geometry.pdf", width: 95%),
  caption: [The Chebyshev--Fourier connection. _Left_: The unit semicircle with equally-spaced angles $theta_j = j pi / N$ projecting to Chebyshev points $x_j = cos theta_j$ on the interval $[-1, 1]$. The clustering near $x = plus.minus 1$ is evident. _Right_: A Chebyshev polynomial $T_n (x)$ can be visualised as a cosine wave $cos(n theta)$ "wrapped around a cylinder and viewed from the side."],
) <fig-cheb-fourier-geometry>

The code generating @fig-cheb-fourier-geometry is available in:
- `codes/python/ch10/cheb_fourier_geometry.py`
- `codes/matlab/ch10/cheb_fourier_geometry.m`
- `codes/julia/ch10/cheb_fourier_geometry.jl`

== Chebyshev Differentiation via FFT <sec-chebfft>

Having the geometric picture, we now build the computational machinery for FFT-based Chebyshev differentiation. Given function values $v_j$ at the Chebyshev points @eq-cheb-points, our goal is to compute accurate approximations to the derivatives $u'(x_j)$ and $u''(x_j)$.

=== The "x to $theta$ to Fourier" Pipeline

The algorithm proceeds through a sequence of transformations. The efficiency of this differentiation relies on the $O(N log N)$ Fast Fourier Transform algorithm introduced by Cooley and Tukey @Cooley1965. The specific implementation using the "extended vector" approach to enforce symmetries is a staple of modern spectral software, comprehensively detailed by Trefethen @Trefethen2000 and Boyd @Boyd2000.

1. *Extend to an even vector*: Given values $v_0, v_1, dots, v_N$ at Chebyshev points, form the extended vector $V$ of length $2 N$ by reflecting:
   $ V_j = v_j quad (j = 0, dots, N), quad V_(2 N - j) = v_j quad (j = 1, dots, N - 1). $
   This extended vector represents an even function sampled at $2 N$ equally-spaced angles in $[0, 2 pi)$.

2. *FFT in $theta$-space*: Compute the discrete Fourier transform $hat(V)_k$ using the FFT.

3. *Differentiate in Fourier space*: Multiply by $i k$ (for the first derivative) or $(i k)^2$ (for the second derivative). For odd-order derivatives, set the Nyquist mode $hat(V)_(N)$ to zero to avoid aliasing artefacts.

4. *Inverse FFT*: Transform back to obtain the derivative $W_j = Q'(theta_j)$ in $theta$-space.

5. *Convert from $theta$ to $x$*: Apply the chain rule to transform derivatives with respect to $theta$ into derivatives with respect to $x$. Since $x = cos theta$, we have $dif x / dif theta = -sin theta = -sqrt(1 - x^2)$, giving
   $ q'(x) = frac(Q'(theta), dif x \/ dif theta) = -frac(Q'(theta), sqrt(1 - x^2)). $ <eq-chain-rule-first>

For the second derivative, the relation is more involved:
$ q''(x) = frac(x, (1 - x^2)^(3\/2)) Q'(theta) + frac(1, 1 - x^2) Q''(theta). $ <eq-chain-rule-second>

6. *Boundary values*: At $x = plus.minus 1$, the formulas @eq-chain-rule-first and @eq-chain-rule-second involve division by zero. Special formulas, derived via l'Hôpital's rule, handle the endpoints:
   $ q'(1) = sum_(n=0)^N n^2 a_n, quad q'(-1) = sum_(n=0)^N (-1)^(n+1) n^2 a_n, $
   where $a_n$ are the Chebyshev coefficients.

=== Implementation: The `chebfft` Function

The following Python code implements the first derivative via FFT on the Chebyshev grid:

```python
def chebfft(v):
    """Chebyshev first derivative via FFT."""
    N = len(v) - 1
    if N == 0:
        return np.array([0.0])

    x = np.cos(np.pi * np.arange(N + 1) / N)

    # Extend to even periodic function in theta
    V = np.concatenate([v, v[-2:0:-1]])

    # FFT and multiply by ik
    U = np.fft.fft(V).real
    k = np.concatenate([np.arange(N), [0], np.arange(1-N, 0)])
    W = np.fft.ifft(1j * k * np.fft.fft(V)).real

    # Transform from theta to x derivative
    w = np.zeros(N + 1)
    w[1:N] = -W[1:N] / np.sqrt(1 - x[1:N]**2)

    # Boundary values via summation formulas
    ii = np.arange(N)
    w[0] = np.sum(ii**2 * U[ii]) / N + 0.5 * N * U[N]
    w[N] = np.sum((-1)**(ii+1) * ii**2 * U[ii]) / N + \
           0.5 * (-1)**(N+1) * N * U[N]

    return w
```

The equivalent MATLAB implementation:

```matlab
function w = chebfft(v)
% CHEBFFT  Chebyshev first derivative via FFT.
    N = length(v) - 1;
    if N == 0, w = 0; return, end

    x = cos((0:N)' * pi / N);
    V = [v(:); flipud(v(2:N))];
    U = real(fft(V));

    ii = 0:N-1;
    W = real(ifft(1i * [ii, 0, 1-N:-1]' .* fft(V)));

    w = zeros(N+1, 1);
    w(2:N) = -W(2:N) ./ sqrt(1 - x(2:N).^2);
    w(1) = sum(ii'.^2 .* U(ii'+1)) / N + 0.5*N*U(N+1);
    w(N+1) = sum((-1).^(ii+1)' .* ii'.^2 .* U(ii'+1)) / N + ...
             0.5*(-1)^(N+1)*N*U(N+1);
end
```

The Julia implementation:

```julia
function chebfft(v)
    N = length(v) - 1
    N == 0 && return [0.0]

    x = [cos(π * j / N) for j in 0:N]

    # Extend to even periodic function in theta
    V = [v; v[N:-1:2]]

    # FFT and multiply by ik
    U = real.(fft(V))
    k = [collect(0:N-1); 0; collect(1-N:-1)]
    W = real.(ifft(1im .* k .* fft(V)))

    # Transform from theta to x derivative
    w = zeros(N + 1)
    w[2:N] = -W[2:N] ./ sqrt.(1.0 .- x[2:N].^2)

    # Boundary values via summation formulas
    ii = 0:N-1
    w[1]   = sum(ii.^2 .* U[ii .+ 1]) / N + 0.5 * N * U[N+1]
    w[N+1] = sum((-1.0).^(ii .+ 1) .* ii.^2 .* U[ii .+ 1]) / N +
             0.5 * (-1.0)^(N + 1) * N * U[N+1]
    return w
end
```

The code generating the figures in this section is available in:
- `codes/python/ch10/chebfft.py`
- `codes/matlab/ch10/chebfft.m`
- `codes/julia/ch10/chebfft.jl`

=== Accuracy Demonstration

To verify spectral accuracy, consider the test function
$ f(x) = e^x cos(4 x), $
with exact derivative
$ f'(x) = e^x (cos(4 x) - 4 sin(4 x)). $

@fig-chebfft-accuracy shows the function, its spectral derivative, and the convergence of the maximum error as $N$ increases.

#figure(
  image("../figures/ch10/python/chebfft_accuracy.pdf", width: 95%),
  caption: [Chebyshev differentiation via FFT applied to $f(x) = e^x cos(4 x)$. _Left_: The function $f(x)$ on $[-1, 1]$ with Chebyshev nodes marked. _Centre_: Comparison of exact and FFT-computed derivative for $N = 16$. _Right_: Maximum error versus $N$ on a semilog scale, showing the characteristic straight-line descent of spectral convergence until machine precision is reached around $N = 24$.],
) <fig-chebfft-accuracy>

The code generating @fig-chebfft-accuracy is available in:
- `codes/python/ch10/chebfft_accuracy.py`
- `codes/matlab/ch10/chebfft_accuracy.m`
- `codes/julia/ch10/chebfft_accuracy.jl`

The straight-line descent in the error plot confirms exponential (spectral) convergence: each additional grid point reduces the error by a roughly constant factor. This is the hallmark of spectral methods applied to smooth functions.

== Method of Lines and Time Stepping <sec-method-of-lines>

Before tackling specific PDEs, we establish the general framework connecting spatial discretisation to time evolution. This _method of lines_ approach treats space and time separately: first discretise in space to obtain a system of ordinary differential equations (ODEs), then apply a suitable time integrator. The Method of Lines was formally systematised by Schiesser @Schiesser1991, providing a framework that decouples spatial accuracy (spectral) from temporal stability (A-stability of ODE solvers). The connection between this separation and the stiffness inherent in Chebyshev spectral derivatives has been explored in the comprehensive review by Hamdi _et al._ @Hamdi2007.

=== Semidiscretisation

Consider a PDE of the form
$ u_t = cal(L) u, $
where $cal(L)$ is a spatial differential operator (e.g., $cal(L) = kappa (partial^2) / (partial x^2)$ for the heat equation). After discretising in space with $N + 1$ grid points, we obtain the _semidiscrete_ system
$ dot(bold(u)) = bold(L) bold(u), $ <eq-semidiscrete-first>
where $bold(u)(t) = (u_0(t), u_1(t), dots, u_N(t))^top$ is the vector of grid values and $bold(L)$ is the matrix (or matrix-free operator) representing the discretised spatial operator.

For second-order in time PDEs like the wave equation $u_(t t) = c^2 u_(x x)$, the semidiscrete form is
$ dot.double(bold(u)) = bold(L) bold(u), $ <eq-semidiscrete-second>
which can be rewritten as a first-order system by introducing the velocity $bold(v) = dot(bold(u))$.

=== The Spatial Operator

The operator $bold(L)$ can be implemented in two ways:

1. *Matrix form*: Build the full Chebyshev differentiation matrix $D_N$ (or $D_N^2$ for second derivatives) as in @ch-chebyshev, then apply it via matrix-vector multiplication.

2. *Matrix-free form*: Use the FFT-based `chebfft` function, avoiding explicit matrix storage. For the second derivative, simply apply `chebfft` twice.

For 2D problems on tensor-product grids, the Laplacian takes the form
$ bold(L) = D_2 times.o I + I times.o D_2, $
where $times.o$ denotes the Kronecker product. In practice, we compute $D_2 U + U D_2^top$ using matrix-matrix multiplication, exploiting the separable structure.

=== Time Integrators by PDE Type

Different PDE types call for different time-stepping strategies. @tbl-time-integrators summarises the typical choices.

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    table(
      columns: (1fr, 1.5fr, 1.5fr, 2fr),
      align: (left, left, left, left),
      inset: (x: 0.8em, y: 0.6em),
      stroke: none,
      table.hline(stroke: 0.75pt + rgb("#142D6E")),
      table.header(
        table.cell(fill: rgb("#142D6E").lighten(85%))[*PDE Type*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Semidiscrete Form*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Typical Schemes*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Stability Notes*],
      ),
      table.hline(stroke: 0.5pt + luma(180)),
      [Wave], [$dot.double(bold(u)) = bold(L) bold(u)$], [Leapfrog, Newmark], [Explicit; $Delta t tilde N^(-2)$],
      [Diffusion], [$dot(bold(u)) = bold(L) bold(u)$], [Backward Euler, CN], [Implicit preferred],
      [Advection], [$dot(bold(u)) = -c bold(D) bold(u)$], [Leapfrog, RK3], [Explicit; $Delta t tilde N^(-2)$],
    ),
  ),
  caption: [Model time integrators for different PDE types. CN denotes Crank--Nicolson.],
) <tbl-time-integrators>

=== The Chebyshev CFL Condition

A crucial point for explicit time stepping: on Chebyshev grids, the time step must scale as
$ Delta t = O(N^(-2)), $
not $O(N^(-1))$ as for equispaced finite differences. This more restrictive stability limit arises from the $O(N^(-2))$ clustering of Chebyshev points near the boundaries. The explicit stability limit for Chebyshev spectral methods was rigorously analysed by Gottlieb and Turkel @GottliebTurkel1980, who demonstrated that the clustering of nodes necessitates time steps scaling as the square of the reciprocal of the number of modes. Trefethen and Trummer @TrefethenTrummer1987 further illuminated this by analysing the spectra of Chebyshev differentiation matrices, revealing that the eigenvalues of the first derivative operator grow as $O(N^2)$, creating the severe stability constraint known as the "Chebyshev penalty."

Physically, information propagates across the smallest grid spacing in each time step. Since the minimum spacing near the boundaries is $O(N^(-2))$, the CFL condition requires correspondingly small time steps. This is the price we pay for the superior spatial accuracy of spectral methods.

Implicit methods (backward Euler, Crank--Nicolson) avoid this restriction entirely, making them attractive for diffusion-dominated problems where time accuracy requirements are less stringent.

== One-Dimensional Wave Equation <sec-wave-1d>

We begin our tour of spectral PDE solvers with the classic vibrating string: the one-dimensional wave equation with fixed endpoints.

=== Problem Formulation

The _wave equation_ in one spatial dimension is
$ u_(t t) = c^2 u_(x x), quad -1 < x < 1, quad t > 0, $ <eq-wave-1d>
where $c > 0$ is the wave speed. The _Dirichlet boundary conditions_
$ u(-1, t) = u(1, t) = 0 $
model a string fixed at both ends. We also need initial conditions:
$ u(x, 0) = u_0(x), quad u_t (x, 0) = v_0(x). $

For our demonstration, we choose a half-sine initial displacement with zero initial velocity:
$ u(x, 0) = sin(frac(pi, 2)(1 + x)), quad u_t (x, 0) = 0. $ <eq-wave-1d-ic>
This initial profile satisfies the boundary conditions and has an explicit eigenfunction expansion, allowing comparison with the exact solution.

=== Spatial Discretisation

We discretise space using $N + 1$ Chebyshev points @eq-cheb-points. The second spatial derivative is approximated using either:
- The Chebyshev differentiation matrix $D_N^2 = D_N times D_N$, or
- Two applications of `chebfft`: $u_(x x) approx$ `chebfft(chebfft(u))`.

To enforce the Dirichlet boundary conditions, we work with the interior nodes $j = 1, 2, dots, N - 1$ and set $u_0 = u_N = 0$ after each time step.

=== Time Stepping: Leapfrog with Taylor Start

For the second-order-in-time wave equation, the _leapfrog_ (Störmer--Verlet) scheme is a natural choice. The seminal comparison of spectral and finite difference methods for hyperbolic equations by Kreiss and Oliger @KreissOliger1972 demonstrated that spectral methods require significantly fewer grid points per wavelength ($N approx pi$) compared to finite differences ($N approx 10$--$20$), making them ideal for long-range wave propagation. The leapfrog scheme is:
$ bold(u)^(n+1) = 2 bold(u)^n - bold(u)^(n-1) + (c Delta t)^2 D_2 bold(u)^n. $ <eq-leapfrog>

This scheme requires values at two previous time levels. For the first step ($n = 0$), we use a Taylor expansion:
$ bold(u)^1 = bold(u)^0 + Delta t dot(bold(u))^0 + frac((Delta t)^2, 2) dot.double(bold(u))^0 = bold(u)^0 + frac((c Delta t)^2, 2) D_2 bold(u)^0, $
where we used $dot(bold(u))^0 = bold(v)_0 = 0$ (zero initial velocity).

For stability, we choose
$ Delta t = frac(alpha, N^2), $
with $alpha approx 4$ to $8$. This $O(N^(-2))$ scaling reflects the Chebyshev CFL condition.

=== Implementation

The following Python code solves the 1D wave equation:

```python
def wave1d_cheb(N=80, c=1.0, tmax=6.0, alpha=4.0):
    """Solve 1D wave equation on Chebyshev grid."""
    # Chebyshev grid
    x = np.cos(np.pi * np.arange(N + 1) / N)

    # Initial data: half-sine
    u = np.sin(0.5 * np.pi * (1 + x))

    # Time stepping
    dt = alpha / N**2
    nsteps = int(np.ceil(tmax / dt))

    # Taylor start (zero initial velocity)
    u_xx = chebfft(chebfft(u))
    u_xx[0] = u_xx[N] = 0  # enforce BCs
    u_prev = u - 0.5 * (c * dt)**2 * u_xx

    for n in range(nsteps):
        u_xx = chebfft(chebfft(u))
        u_xx[0] = u_xx[N] = 0
        u_new = 2*u - u_prev + (c * dt)**2 * u_xx
        u_prev, u = u, u_new

    return x, u
```

The equivalent MATLAB implementation:

```matlab
function [x, u] = wave1d_cheb(N, c, tmax, alpha)
% Solve 1D wave equation on Chebyshev grid.
    x = cos(pi * (0:N)' / N);

    % Initial data: half-sine
    u = sin(0.5 * pi * (1 + x));

    dt = alpha / N^2;
    nsteps = ceil(tmax / dt);

    % Taylor start
    u_xx = chebfft(chebfft(u));
    u_xx(1) = 0; u_xx(N+1) = 0;
    u_prev = u - 0.5 * (c*dt)^2 * u_xx;

    for n = 1:nsteps
        u_xx = chebfft(chebfft(u));
        u_xx(1) = 0; u_xx(N+1) = 0;
        u_new = 2*u - u_prev + (c*dt)^2 * u_xx;
        u_prev = u; u = u_new;
    end
end
```

The Julia implementation:

```julia
function wave1d_cheb(; N=80, c=1.0, tmax=6.0, alpha=4.0)
    x = [cos(pi * j / N) for j in 0:N]

    # Initial data: half-sine
    u = sin.(0.5 * pi .* (1.0 .+ x))

    dt = alpha / N^2
    nsteps = Int(ceil(tmax / dt))

    # Taylor start (zero initial velocity)
    u_xx = chebfft(chebfft(u))
    u_xx[1] = 0.0; u_xx[N+1] = 0.0
    u_prev = u .- 0.5 * (c * dt)^2 .* u_xx

    for n in 1:nsteps
        u_xx = chebfft(chebfft(u))
        u_xx[1] = 0.0; u_xx[N+1] = 0.0
        u_new = 2.0 .* u .- u_prev .+ (c * dt)^2 .* u_xx
        u_prev = u; u = u_new
    end
    return x, u
end
```

=== Results and Physical Interpretation

@fig-wave-1d-waterfall shows the solution as a waterfall plot in space and time.

#figure(
  image("../figures/ch10/python/wave1d_waterfall.pdf", width: 90%),
  caption: [Solution of the 1D wave equation @eq-wave-1d with initial half-sine profile. The waterfall plot shows $u(x, t)$ for $t in [0, 6]$. The initial shape splits into two counter-propagating waves, which reflect at the fixed boundaries ($x = plus.minus 1$) with a sign change. The waves pass through each other and interfere, creating the characteristic standing wave pattern.],
) <fig-wave-1d-waterfall>

The physics is clearly visible: the initial displacement splits into left- and right-moving pulses (d'Alembert's solution). Upon reaching the boundaries, each pulse reflects with a sign change (required by the Dirichlet condition). As the waves pass through each other, they superpose linearly, creating the complex interference patterns visible in the waterfall plot.

The code generating @fig-wave-1d-waterfall is available in:
- `codes/python/ch10/wave1d_cheb.py`
- `codes/matlab/ch10/wave1d_cheb.m`
- `codes/julia/ch10/wave1d_cheb.jl`

== Two-Dimensional Wave Equation <sec-wave-2d>

We now extend the wave equation to two spatial dimensions: the vibrating drum membrane.

=== Problem Formulation

The 2D wave equation on the square $[-1, 1]^2$ is
$ u_(t t) = c^2 (u_(x x) + u_(y y)), quad -1 < x, y < 1, quad t > 0, $ <eq-wave-2d>
with homogeneous Dirichlet boundary conditions $u = 0$ on $partial Omega$ (the membrane is fixed along its entire boundary).

For initial data, we choose an offset Gaussian bump:
$ u(x, y, 0) = exp(-30[(x - 0.2)^2 + (y + 0.3)^2]), quad u_t (x, y, 0) = 0. $ <eq-wave-2d-ic>
The asymmetric position $(0.2, -0.3)$ breaks the symmetry of the domain, producing richer interference patterns.

=== Spatial Discretisation: Tensor Products

We use a _tensor product grid_: Chebyshev points in both $x$ and $y$ directions:
$ x_i = cos(i pi / N), quad y_j = cos(j pi / N), quad i, j = 0, 1, dots, N. $

The solution is represented as an $(N + 1) times (N + 1)$ matrix $U$ with entries $U_(i j) approx u(x_i, y_j)$.

The key insight is that the 2D Laplacian separates into 1D operators:
$ Delta u approx D_2 U + U D_2^top, $ <eq-tensor-laplacian>
where $D_2 = D_N^2$ is the 1D second-derivative matrix. The first term applies $D_2$ to each row (differentiation in $x$); the second applies $D_2^top$ to each column (differentiation in $y$).

This tensor-product structure is computationally efficient: applying @eq-tensor-laplacian costs $O(N^3)$ operations (two matrix-matrix products), compared to $O(N^4)$ for a naive implementation using the full $(N + 1)^2 times (N + 1)^2$ Laplacian matrix.

=== Time Stepping

Leapfrog extends naturally to 2D:
$ U^(n+1) = 2 U^n - U^(n-1) + (c Delta t)^2 (D_2 U^n + U^n D_2^top). $ <eq-leapfrog-2d>

After each time step, we enforce the boundary conditions by setting the boundary rows and columns of $U$ to zero.

For stability, we take $Delta t = alpha / N^2$ with $alpha approx 3$ to $6$.

=== Implementation

```python
def wave2d_cheb(N=32, c=1.0, tmax=1.0, alpha=3.0):
    """Solve 2D wave equation on Chebyshev tensor grid."""
    # Grid
    x = np.cos(np.pi * np.arange(N + 1) / N)
    xx, yy = np.meshgrid(x, x)

    # Build D2 matrix
    D, _ = cheb_matrix(N)
    D2 = D @ D

    # Initial data: offset Gaussian
    U = np.exp(-30 * ((xx - 0.2)**2 + (yy + 0.3)**2))

    dt = alpha / N**2
    nsteps = int(np.ceil(tmax / dt))

    # Taylor start
    Lap_U = D2 @ U + U @ D2.T
    U_prev = U - 0.5 * (c * dt)**2 * Lap_U

    for n in range(nsteps):
        Lap_U = D2 @ U + U @ D2.T
        U_new = 2*U - U_prev + (c * dt)**2 * Lap_U

        # Enforce boundary conditions
        U_new[0, :] = U_new[-1, :] = 0
        U_new[:, 0] = U_new[:, -1] = 0

        U_prev, U = U, U_new

    return xx, yy, U
```

The equivalent MATLAB implementation:

```matlab
function [xx, yy, U] = wave2d_cheb(N, c, tmax, alpha)
% Solve 2D wave equation on Chebyshev tensor grid.
    x = cos(pi * (0:N)' / N);
    [xx, yy] = meshgrid(x, x);

    [D, ~] = cheb_matrix(N);
    D2 = D * D;

    % Initial data
    U = exp(-30 * ((xx - 0.2).^2 + (yy + 0.3).^2));

    dt = alpha / N^2;
    nsteps = ceil(tmax / dt);

    Lap_U = D2 * U + U * D2';
    U_prev = U - 0.5 * (c*dt)^2 * Lap_U;

    for n = 1:nsteps
        Lap_U = D2 * U + U * D2';
        U_new = 2*U - U_prev + (c*dt)^2 * Lap_U;
        U_new(1,:) = 0; U_new(end,:) = 0;
        U_new(:,1) = 0; U_new(:,end) = 0;
        U_prev = U; U = U_new;
    end
end
```

The Julia implementation:

```julia
function wave2d_cheb(; N=32, c=1.0, tmax=1.0, alpha=3.0)
    x = [cos(pi * j / N) for j in 0:N]
    xx = [x[i] for i in 1:N+1, j in 1:N+1]
    yy = [x[j] for i in 1:N+1, j in 1:N+1]

    D, _ = cheb_matrix(N)
    D2 = D * D

    # Initial data: offset Gaussian
    U = exp.(-30.0 .* ((xx .- 0.2).^2 .+ (yy .+ 0.3).^2))

    dt = alpha / N^2
    nsteps = Int(ceil(tmax / dt))

    # Taylor start
    Lap_U = D2 * U + U * D2'
    U_prev = U .- 0.5 * (c * dt)^2 .* Lap_U

    for n in 1:nsteps
        Lap_U = D2 * U + U * D2'
        U_new = 2.0 .* U .- U_prev .+ (c * dt)^2 .* Lap_U
        U_new[1, :] .= 0.0;   U_new[end, :] .= 0.0
        U_new[:, 1] .= 0.0;   U_new[:, end] .= 0.0
        U_prev = U; U = U_new
    end
    return xx, yy, U
end
```

=== Results

@fig-wave-2d-snapshots shows four snapshots of the drum membrane at different times.

#figure(
  image("../figures/ch10/python/wave2d_snapshots.pdf", width: 95%),
  caption: [Solution of the 2D wave equation @eq-wave-2d showing a vibrating drum membrane. Four snapshots at times $t = 0$, $0.33$, $0.67$, and $1.0$. The initial Gaussian bump spreads into a circular wavefront, which reflects off the fixed boundaries and creates complex interference patterns. The tensor-product Chebyshev grid clusters points near all four edges.],
) <fig-wave-2d-snapshots>

The initial Gaussian bump spreads radially as a circular wavefront. Upon reaching the boundaries, the wave reflects and interferes with itself, creating intricate patterns. Note that despite using only $N = 24$ points in each direction (576 total grid points), the solution captures the wave dynamics accurately, a testament to spectral accuracy.

The code generating @fig-wave-2d-snapshots is available in:
- `codes/python/ch10/wave2d_cheb.py`
- `codes/matlab/ch10/wave2d_cheb.m`
- `codes/julia/ch10/wave2d_cheb.jl`

== One-Dimensional Heat Equation <sec-heat-1d>

We now turn to parabolic equations, where diffusion rather than wave propagation dominates the physics.

=== Problem Formulation

The _heat equation_ (or diffusion equation) in one dimension is
$ u_t = kappa u_(x x), quad -1 < x < 1, quad t > 0, $ <eq-heat-1d>
where $kappa > 0$ is the thermal diffusivity. With homogeneous Dirichlet conditions $u(plus.minus 1, t) = 0$, this models a rod with endpoints held at zero temperature.

For initial data, we choose a two-bump profile:
$ u(x, 0) = exp(-20(x + 0.5)^2) + 0.5 exp(-30(x - 0.4)^2). $ <eq-heat-1d-ic>
This allows us to observe how two separate hot spots diffuse and eventually merge.

=== Time Stepping: Crank--Nicolson

Unlike the wave equation, the heat equation is _dissipative_: solutions decay over time as heat diffuses away. Explicit methods face severe stability restrictions ($Delta t lt.eq.slant O(N^(-4))$ for spectral methods!), making implicit schemes essential. The Implicit-Explicit (IMEX) strategy, formalized by Ascher, Ruuth, and Wetton @Ascher1995, treats the stiff linear diffusion implicitly while keeping nonlinear terms explicit. For highly stiff systems with diagonal linear operators, Kassam and Trefethen @KassamTrefethen2005 demonstrated the superiority of Exponential Time Differencing (ETD) methods, which analytically incorporate the stiff linear part.

The _Crank--Nicolson_ scheme averages explicit and implicit Euler:
$ frac(bold(u)^(n+1) - bold(u)^n, Delta t) = frac(kappa, 2) (D_(2,i) bold(u)^(n+1) + D_(2,i) bold(u)^n), $ <eq-crank-nicolson>
where $D_(2,i)$ is the interior block of $D_2$ (rows and columns 1 through $N - 1$).

Rearranging gives the linear system:
$ (I - frac(kappa Delta t, 2) D_(2,i)) bold(u)^(n+1) = (I + frac(kappa Delta t, 2) D_(2,i)) bold(u)^n. $

Crank--Nicolson is:
- *Unconditionally stable*: no CFL restriction on $Delta t$
- *Second-order accurate* in time: $O((Delta t)^2)$ local truncation error
- *Implicit*: requires solving a linear system at each step

The linear system is dense (because $D_(2,i)$ is dense), but for moderate $N$, direct solution via LU factorisation is efficient. For very large $N$, iterative methods become preferable.

=== Implementation

```python
def heat1d_cheb(N=64, kappa=1.0, tmax=1.0, dt=0.01):
    """Solve 1D heat equation with Crank-Nicolson."""
    D, x = cheb_matrix(N)
    D2 = D @ D

    # Interior block
    ii = slice(1, N)
    D2i = D2[np.ix_(np.arange(1, N), np.arange(1, N))]
    I = np.eye(N - 1)

    # Crank-Nicolson matrices
    A = I - 0.5 * kappa * dt * D2i
    B = I + 0.5 * kappa * dt * D2i

    # Initial condition (interior points only)
    u_full = np.exp(-20*(x + 0.5)**2) + 0.5*np.exp(-30*(x - 0.4)**2)
    u = u_full[ii]

    nsteps = int(np.ceil(tmax / dt))
    for n in range(nsteps):
        u = np.linalg.solve(A, B @ u)

    # Reconstruct full solution with BCs
    u_full = np.zeros(N + 1)
    u_full[ii] = u
    return x, u_full
```

The equivalent MATLAB implementation:

```matlab
function [x, u_full] = heat1d_cheb(N, kappa, tmax, dt)
% Solve 1D heat equation with Crank-Nicolson.
    [D, x] = cheb_matrix(N);
    D2 = D * D;

    % Interior block
    D2i = D2(2:N, 2:N);
    I = eye(N - 1);

    A = I - 0.5 * kappa * dt * D2i;
    B = I + 0.5 * kappa * dt * D2i;

    u_full = exp(-20*(x + 0.5).^2) + 0.5*exp(-30*(x - 0.4).^2);
    u = u_full(2:N);

    nsteps = ceil(tmax / dt);
    for n = 1:nsteps
        u = A \ (B * u);
    end

    u_full = zeros(N+1, 1);
    u_full(2:N) = u;
end
```

The Julia implementation:

```julia
function heat1d_cheb(; N=64, kappa=1.0, tmax=1.0, dt=0.01)
    D, x = cheb_matrix(N)
    D2 = D * D

    # Interior block
    ii = 2:N
    D2i = D2[ii, ii]
    Ie = Matrix{Float64}(I, N - 1, N - 1)

    # Crank-Nicolson matrices
    A = Ie - 0.5 * kappa * dt * D2i
    B = Ie + 0.5 * kappa * dt * D2i

    # Initial condition (interior points only)
    u_full = exp.(-20.0 .* (x .+ 0.5).^2) .+
             0.5 .* exp.(-30.0 .* (x .- 0.4).^2)
    u = u_full[ii]

    nsteps = Int(ceil(tmax / dt))
    for n in 1:nsteps
        u = A \ (B * u)
    end

    # Reconstruct full solution with BCs
    u_full = zeros(N + 1)
    u_full[ii] = u
    return x, u_full
end
```

=== Results

@fig-heat-1d-evolution shows the temperature evolution as a space--time surface and as four snapshots.

#figure(
  image("../figures/ch10/python/heat1d_evolution.pdf", width: 95%),
  caption: [Solution of the 1D heat equation @eq-heat-1d with two-bump initial data. _Main panel_: Space--time surface showing diffusion and decay. _Insets_: Four snapshots at $t = 0$, $0.1$, $0.3$, and $1.0$. The two initial hot spots spread, merge, and eventually decay towards the zero steady state required by the boundary conditions.],
) <fig-heat-1d-evolution>

The physics of diffusion is evident: sharp features smooth out, amplitudes decrease, and the solution approaches the steady state $u equiv 0$. The two initial bumps spread and eventually merge into a single, broader profile before decaying completely.

The code generating @fig-heat-1d-evolution is available in:
- `codes/python/ch10/heat1d_cheb.py`
- `codes/matlab/ch10/heat1d_cheb.m`
- `codes/julia/ch10/heat1d_cheb.jl`

== Two-Dimensional Heat Equation <sec-heat-2d>

The extension to 2D follows the same pattern as the wave equation: tensor-product grids and Kronecker structure.

=== Problem Formulation

The 2D heat equation on the square $[-1, 1]^2$ is
$ u_t = u_(x x) + u_(y y), quad -1 < x, y < 1, quad t > 0, $ <eq-heat-2d>
with homogeneous Dirichlet conditions $u = 0$ on $partial Omega$.

Initial data:
$ u(x, y, 0) = exp(-25[(x + 0.3)^2 + (y - 0.1)^2]). $ <eq-heat-2d-ic>

=== Time Stepping: Backward Euler with Kronecker Structure

The semidiscrete system is
$ U_t = D_2 U + U D_2^top. $

For implicit time stepping (e.g., backward Euler), we must solve
$ (I - Delta t bold(L)) op("vec")(U^(n+1)) = op("vec")(U^n), $
where $bold(L) = D_(2,i) times.o I + I times.o D_(2,i)$ is the _Kronecker sum_ of the interior blocks.

The Kronecker sum $bold(L)$ is an $(N - 1)^2 times (N - 1)^2$ matrix. For moderate $N$ (say, $N lt.eq.slant 64$), we can form it explicitly and solve with a direct method. For larger $N$, iterative methods (conjugate gradient, GMRES) or specialised tensor solvers are more efficient.

=== Implementation

```python
def heat2d_cheb(N=32, tmax=0.5, dt=0.01):
    """Solve 2D heat equation with backward Euler."""
    D, x = cheb_matrix(N)
    D2 = D @ D

    xx, yy = np.meshgrid(x, x)

    # Interior indices
    ii = np.arange(1, N)
    D2i = D2[np.ix_(ii, ii)]
    I = np.eye(N - 1)

    # Kronecker sum Laplacian
    L = np.kron(D2i, I) + np.kron(I, D2i)
    A = np.eye((N-1)**2) - dt * L

    # Initial condition
    U_full = np.exp(-25 * ((xx + 0.3)**2 + (yy - 0.1)**2))
    u = U_full[np.ix_(ii, ii)].flatten()

    nsteps = int(np.ceil(tmax / dt))
    for n in range(nsteps):
        u = np.linalg.solve(A, u)

    # Reconstruct full solution
    U_full = np.zeros((N+1, N+1))
    U_full[np.ix_(ii, ii)] = u.reshape(N-1, N-1)
    return xx, yy, U_full
```

The equivalent MATLAB implementation:

```matlab
function [xx, yy, U] = heat2d_cheb(N, tmax, dt)
    [D, x] = cheb_matrix(N);
    D2 = D * D;
    [xx, yy] = meshgrid(x, x);

    % Interior block and Kronecker sum
    D2i = D2(2:N, 2:N);
    I = eye(N - 1);
    L = kron(D2i, I) + kron(I, D2i);
    A = eye((N-1)^2) - dt * L;

    % Initial condition
    U = exp(-25 * ((xx+0.3).^2 + (yy-0.1).^2));
    U(1,:) = 0; U(N+1,:) = 0;
    U(:,1) = 0; U(:,N+1) = 0;

    nsteps = ceil(tmax / dt);
    for n = 1:nsteps
        u = reshape(U(2:N, 2:N), [], 1);
        U = zeros(N+1, N+1);
        U(2:N, 2:N) = reshape(A \ u, N-1, N-1);
    end
end
```

The Julia implementation:

```julia
function heat2d_cheb(; N=32, tmax=0.5, dt=0.01)
    D, x = cheb_matrix(N)
    D2 = D * D
    xx = [x[i] for i in 1:N+1, j in 1:N+1]
    yy = [x[j] for i in 1:N+1, j in 1:N+1]

    # Interior block and Kronecker sum
    D2i = D2[2:N, 2:N]
    Imat = Matrix{Float64}(I, N-1, N-1)
    L = kron(D2i, Imat) + kron(Imat, D2i)
    A = Matrix{Float64}(I, (N-1)^2, (N-1)^2) - dt * L

    # Initial condition
    U = exp.(-25 .* ((xx .+ 0.3).^2 .+ (yy .- 0.1).^2))
    u = vec(U[2:N, 2:N])

    nsteps = ceil(Int, tmax / dt)
    for n in 1:nsteps
        u = A \ u
    end

    U = zeros(N+1, N+1)
    U[2:N, 2:N] = reshape(u, N-1, N-1)
    return xx, yy, U
end
```

=== Results

@fig-heat-2d-snapshots shows four snapshots of the diffusing hot spot.

#figure(
  image("../figures/ch10/python/heat2d_snapshots.pdf", width: 95%),
  caption: [Solution of the 2D heat equation @eq-heat-2d. Four snapshots at times $t = 0$, $0.05$, $0.15$, and $0.5$. The initial hot spot diffuses isotropically, spreading and decaying towards the zero steady state.],
) <fig-heat-2d-snapshots>

The code generating @fig-heat-2d-snapshots is available in:
- `codes/python/ch10/heat2d_cheb.py`
- `codes/matlab/ch10/heat2d_cheb.m`
- `codes/julia/ch10/heat2d_cheb.jl`

== Two-Dimensional Poisson Equation <sec-poisson-2d>

We now turn to _elliptic_ equations, which describe equilibrium (steady-state) phenomena rather than time evolution.

=== Problem Formulation

The _Poisson equation_ on the square $[-1, 1]^2$ is
$ -Delta u(x, y) = f(x, y), quad -1 < x, y < 1, $ <eq-poisson-spectral>
with Dirichlet boundary conditions $u = g$ on $partial Omega$.

Physically, @eq-poisson-spectral describes:
- *Electrostatics*: $u$ is the electric potential, $f$ is the charge density
- *Steady heat conduction*: $u$ is the temperature, $f$ is the heat source
- *Membrane deflection*: $u$ is the vertical displacement, $f$ is the applied load

=== Manufactured Solution

To test our solver, we use the _method of manufactured solutions_: choose an exact solution, compute the corresponding right-hand side $f$, and verify that the numerical solution matches.

We take
$ u(x, y) = (1 - x^2)(1 - y^2) cos(frac(pi x, 2)) cos(frac(pi y, 2)). $ <eq-poisson-manufactured>
The factors $(1 - x^2)$ and $(1 - y^2)$ ensure that $u = 0$ on the boundary (homogeneous Dirichlet conditions). The right-hand side $f = -Delta u$ can be computed symbolically or numerically by applying the discrete Laplacian to the exact solution values.

=== Discrete System

The discrete Poisson equation is a linear system:
$ bold(A) bold(u) = bold(f), $
where:
- $bold(u)$ contains the $(N - 1)^2$ interior unknowns, stacked column-by-column
- $bold(f)$ is the corresponding right-hand side
- $bold(A) = D_(2,i) times.o I + I times.o D_(2,i)$ is the Kronecker sum of interior second-derivative blocks

The efficiency of spectral elliptic solvers relies on exploiting the structure of the operator. Haidvogel and Zang @HaidvogelZang1979 introduced the matrix diagonalization method, which decouples the multidimensional system into a sequence of independent one-dimensional problems, reducing the complexity from $O(N^4)$ to $O(N^3)$. Shen @Shen1995 later developed spectral-Galerkin methods using basis functions that satisfy the boundary conditions automatically, leading to sparse matrices and optimal $O(N^2)$ complexity.

=== Implementation

```python
def poisson2d_cheb(N=32):
    """Solve 2D Poisson equation with manufactured solution."""
    D, x = cheb_matrix(N)
    D2 = D @ D

    xx, yy = np.meshgrid(x, x)

    # Interior indices
    ii = np.arange(1, N)
    D2i = D2[np.ix_(ii, ii)]
    I = np.eye(N - 1)

    # Kronecker sum Laplacian
    A = np.kron(D2i, I) + np.kron(I, D2i)

    # Manufactured exact solution
    U_exact = (1 - xx**2) * (1 - yy**2) * \
              np.cos(np.pi*xx/2) * np.cos(np.pi*yy/2)

    # RHS from discrete Laplacian of exact solution
    Lap_U = D2 @ U_exact + U_exact @ D2.T
    f = -Lap_U[np.ix_(ii, ii)]

    # Solve
    u_vec = np.linalg.solve(A, f.flatten())

    # Reconstruct full solution
    U = np.zeros((N+1, N+1))
    U[np.ix_(ii, ii)] = u_vec.reshape(N-1, N-1)

    return xx, yy, U, U_exact
```

The equivalent MATLAB implementation:

```matlab
function [xx, yy, U, U_exact] = poisson2d_cheb(N)
% Solve 2D Poisson equation with manufactured solution.
    [D, x] = cheb_matrix(N);
    D2 = D * D;

    [xx, yy] = meshgrid(x, x);

    ii = 2:N;
    D2i = D2(ii, ii);
    I = eye(N - 1);

    A = kron(D2i, I) + kron(I, D2i);

    U_exact = (1 - xx.^2) .* (1 - yy.^2) .* ...
              cos(pi*xx/2) .* cos(pi*yy/2);

    Lap_U = D2 * U_exact + U_exact * D2';
    f = -Lap_U(ii, ii);

    u_vec = A \ f(:);

    U = zeros(N+1, N+1);
    U(ii, ii) = reshape(u_vec, N-1, N-1);
end
```

The Julia implementation:

```julia
function poisson2d_cheb(; N=32)
    D, x = cheb_matrix(N)
    D2 = D * D

    xx = [x[i] for i in 1:N+1, j in 1:N+1]
    yy = [x[j] for i in 1:N+1, j in 1:N+1]

    # Interior indices
    ii = 2:N
    D2i = D2[ii, ii]
    Ie = Matrix{Float64}(I, N - 1, N - 1)

    # Kronecker sum Laplacian
    L = kron(D2i, Ie) + kron(Ie, D2i)

    # Manufactured exact solution
    U_exact = (1 .- xx.^2) .* (1 .- yy.^2) .*
              cos.(pi .* xx ./ 2) .* cos.(pi .* yy ./ 2)

    # RHS from discrete Laplacian of exact solution
    Lap_U = D2 * U_exact + U_exact * D2'
    f = -Lap_U[ii, ii]

    # Solve -L * u = f
    u_vec = (-L) \ vec(f)

    # Reconstruct full solution
    U = zeros(N + 1, N + 1)
    U[ii, ii] = reshape(u_vec, N - 1, N - 1)
    return xx, yy, U, U_exact
end
```

=== Results

@fig-poisson-2d-solution compares the exact and numerical solutions.

#figure(
  image("../figures/ch10/python/poisson2d_solution.pdf", width: 95%),
  caption: [Solution of the 2D Poisson equation @eq-poisson-spectral with manufactured solution @eq-poisson-manufactured. _Left_: Exact solution. _Centre_: Numerical solution (visually indistinguishable from exact). _Right_: Error on a logarithmic scale, showing machine-precision accuracy throughout the interior. The spectral method achieves essentially exact results for this smooth problem.],
) <fig-poisson-2d-solution>

The error is at machine precision (roughly $10^(-14)$ to $10^(-15)$), confirming spectral accuracy. This is the ideal case: a smooth manufactured solution where the only errors come from finite-precision arithmetic.

The code generating @fig-poisson-2d-solution is available in:
- `codes/python/ch10/poisson2d_cheb.py`
- `codes/matlab/ch10/poisson2d_cheb.m`
- `codes/julia/ch10/poisson2d_cheb.jl`

=== Laplace Equation

The _Laplace equation_ is the special case $f = 0$:
$ Delta u = 0. $
Solutions are called _harmonic functions_. Physically, this describes equilibrium with no sources: the electrostatic potential in a charge-free region, or the steady temperature in a plate with no internal heat sources.

For Laplace's equation, the solution is entirely determined by the boundary conditions. The code adapts trivially: set the interior right-hand side to zero and include any nonzero boundary data in the system.

=== Helmholtz Equation

The _Helmholtz equation_ adds a "mass" term:
$ -Delta u - k^2 u = f(x, y), $ <eq-helmholtz-spectral>
where $k > 0$ is the wavenumber. This equation arises in time-harmonic wave problems (acoustics, electromagnetics) and quantum mechanics.

The discrete system becomes
$ (-bold(L) - k^2 I) bold(u) = bold(f). $

As $k$ increases, the matrix $-bold(L) - k^2 I$ can become ill-conditioned, especially when $k^2$ approaches an eigenvalue of $-bold(L)$. Near such _resonance_ values, the problem becomes nearly singular, and iterative solvers may struggle. This is a fundamental physical phenomenon: near resonance, small forcing produces large response.

== Additional Physics Études <sec-physics-etudes>

This section presents three additional examples demonstrating the versatility of spectral methods. These études are optional but illustrate how the same FFT-based machinery applies across different physical contexts.

=== Variable-Coefficient Transport

Consider advection with spatially varying velocity:
$ u_t + c(x) u_x = 0, quad 0 < x < 2 pi, quad u "periodic", $ <eq-transport-variable>
with
$ c(x) = 0.3 + sin^2(x - 1). $

This models transport in a medium where the wave speed depends on position. Using a _periodic Fourier grid_ with FFT-based differentiation and RK4 time stepping, we can track a pulse as it moves through the variable medium.

#figure(
  image("../figures/ch10/python/transport_variable.pdf", width: 85%),
  caption: [Solution of the variable-coefficient transport equation @eq-transport-variable. Left: space-time diagram showing the characteristic trajectory of a Gaussian pulse in a medium where $c(x) = 0.3 + sin^2(x - 1)$. The pulse accelerates in fast regions and decelerates in slow regions. Right: waterfall view of the same solution with vertical offset proportional to time.],
) <fig-transport-variable>

The code generating @fig-transport-variable is available in:
- `codes/python/ch10/transport_variable.py`
- `codes/matlab/ch10/transport_variable.m`
- `codes/julia/ch10/transport_variable.jl`

=== Schrödinger Equation

The time-dependent Schrödinger equation with a harmonic potential is
$ i u_t = -u_(x x) + x^2 u, $ <eq-schrodinger>
where $u(x, t)$ is the complex wavefunction and $x^2$ is the harmonic oscillator potential.

A natural approach is _Strang splitting_, an operator splitting strategy introduced by Strang @Strang1968 that has become the method of choice for dispersive equations. Taha and Ablowitz @Taha1984 established the split-step Fourier method as the robust standard for the nonlinear Schrödinger equation, while Bao, Jin, and Markowich @BaoJinMarkowich2002 provided rigorous error analysis for time-splitting spectral approximations in the semiclassical regime. The algorithm proceeds in three sub-steps:
1. *Half potential step* (in physical space): $u arrow u dot e^(-i x^2 Delta t \/ 2)$
2. *Full kinetic step* (in Fourier space): $hat(u)_k arrow hat(u)_k dot e^(-i k^2 Delta t)$
3. *Half potential step*: $u arrow u dot e^(-i x^2 Delta t \/ 2)$

This splitting is second-order accurate in time and preserves the $L^2$ norm of the wavefunction (unitarity). The kinetic step uses exactly the same FFT machinery as our other spectral methods.

#figure(
  image("../figures/ch10/python/schrodinger.pdf", width: 85%),
  caption: [Solution of the Schrödinger equation @eq-schrodinger with harmonic potential using Strang splitting. Left: space-time evolution of the probability density $|u(x, t)|^2$ showing the oscillatory dynamics of a Gaussian wavepacket in the harmonic trap. Right: waterfall plot of $|u(x, t)|^2$ at selected times.],
) <fig-schrodinger>

The code generating @fig-schrodinger is available in:
- `codes/python/ch10/schrodinger.py`
- `codes/matlab/ch10/schrodinger.m`
- `codes/julia/ch10/schrodinger.jl`

=== Nonlinear Reaction-Diffusion

The Allen--Cahn equation combines diffusion with a bistable nonlinearity:
$ u_t = u_(x x) + u(1 - u^2). $ <eq-allen-cahn>

The nonlinear term $u(1 - u^2)$ has stable equilibria at $u = plus.minus 1$ and an unstable equilibrium at $u = 0$. Solutions tend to separate into regions where $u approx 1$ and $u approx -1$, with thin transition layers (fronts) between them. The Allen--Cahn equation, originally formulated by Allen and Cahn @AllenCahn1979 to model antiphase boundary motion in crystalline solids, demands numerical schemes that respect the thermodynamic structure of the system. Shen and Yang @ShenYang2010 introduced unconditionally energy-stable spectral schemes that preserve the dissipation of the Ginzburg--Landau free energy, building on the convex splitting framework of Eyre @Eyre1998.

An _IMEX_ (implicit-explicit) scheme treats the stiff linear diffusion implicitly and the nonlinear reaction explicitly:
$ (I - Delta t D_(2,i)) bold(u)^(n+1) = bold(u)^n + Delta t (bold(u)^n - (bold(u)^n)^3), $
where the cubic is applied componentwise.

#figure(
  image("../figures/ch10/python/allen_cahn.pdf", width: 85%),
  caption: [Solution of the Allen--Cahn equation @eq-allen-cahn using IMEX time stepping. Left: space-time evolution showing phase separation into regions where $u approx plus.minus 1$. Right: snapshots at selected times showing the formation and coarsening of fronts between the two stable phases.],
) <fig-allen-cahn>

The code generating @fig-allen-cahn is available in:
- `codes/python/ch10/allen_cahn.py`
- `codes/matlab/ch10/allen_cahn.m`
- `codes/julia/ch10/allen_cahn.jl`

== Efficiency: Matrices versus FFT <sec-efficiency>

Throughout this chapter, we have presented both matrix-based and FFT-based approaches. @tbl-complexity summarises the computational complexity.

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    table(
      columns: (1.5fr, 1fr, 1fr),
      align: (left, center, center),
      inset: (x: 1em, y: 0.6em),
      stroke: none,
      table.hline(stroke: 0.75pt + rgb("#142D6E")),
      table.header(
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Operation*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Matrix-Based*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*FFT-Based*],
      ),
      table.hline(stroke: 0.5pt + luma(180)),
      [1D first derivative], [$O(N^2)$], [$O(N log N)$],
      [1D second derivative], [$O(N^2)$], [$O(N log N)$],
      [2D Laplacian (general)], [$O(N^4)$], [$O(N^2 log N)$],
      [2D Laplacian (tensor)], [$O(N^3)$], [$O(N^2 log N)$],
    ),
  ),
  caption: [Computational complexity of spectral differentiation by method. The matrix-based approach uses dense matrix-vector or matrix-matrix multiplication. The FFT-based approach uses the Chebyshev--Fourier connection. The tensor-product form $D_2 U + U D_2^top$ achieves $O(N^3)$ without explicit Kronecker products.],
) <tbl-complexity>

For small to moderate $N$ (say, $N lt.eq.slant 128$), the matrix approach is often simpler and fast enough. For larger $N$, the FFT approach becomes essential.

A further improvement, not explored here, uses the _Discrete Cosine Transform (DCT)_ instead of the full FFT. Since the extended Chebyshev data is both real and even, the DCT can exploit both symmetries for an additional factor of 2--4 speedup. See the references for details.

== A non-exhaustive literature overview

The development of spectral methods for partial differential equations is a rich narrative that intertwines approximation theory, linear algebra, and computational physics. What began as a technique for high-accuracy meteorological forecasting has evolved into a dominant paradigm for simulating turbulence, quantum dynamics, and general relativity. This overview highlights the seminal contributions that shaped the field and points toward contemporary advancements.

While the roots of spectral expansions can be traced to Fourier and Chebyshev, the computational era of spectral PDEs properly began in the early 1970s. The seminal work of Kreiss and Oliger @KreissOliger1972 compared Fourier methods with finite difference schemes for hyperbolic equations, famously demonstrating that for a fixed error tolerance, spectral methods require significantly fewer grid points per wavelength ($N approx pi$) compared to finite differences ($N approx 10$--$20$). This efficiency was critical for problems where memory was a bottleneck. Simultaneously, Orszag @Orszag1971 applied spectral methods to the Orr--Sommerfeld stability equation, solving problems in fluid stability that were previously intractable. This work, along with his development of transform methods to evaluate convolution sums efficiently, marked the transition from "spectral Galerkin" methods (working entirely in coefficient space) to "pseudospectral" or "collocation" methods (working on physical grids), which we have utilised throughout this chapter. The theoretical bedrock for these methods was laid in the monograph by Gottlieb and Orszag @GottliebOrszag1977, which provided the first rigorous error analysis, establishing the principle of "spectral accuracy" --- that for analytic functions, the error decays exponentially with the number of modes.

A central challenge in spectral PDE solvers is the stiffness introduced by the high-order spatial discretisation. As derived in this chapter, the second-derivative operator on a Chebyshev grid scales as $O(N^4)$, imposing severe time-step restrictions for explicit schemes. Gottlieb and Turkel @GottliebTurkel1980 and Trefethen and Trummer @TrefethenTrummer1987 analysed these stability limits, revealing that the eigenvalues of spectral differentiation matrices are extremely sensitive to boundary conditions and grid distribution. Trefethen's analysis of the "pseudo-spectra" of these non-normal matrices provided deep insight into why spectral schemes could be sensitive to rounding errors and initial transients. To mitigate these restrictions, the community moved toward implicit and semi-implicit schemes. The Implicit-Explicit (IMEX) strategies formalised by Ascher, Ruuth, and Wetton @Ascher1995 became standard for convection-diffusion problems, allowing the stiff diffusive terms to be treated implicitly while the nonlinear convective terms remain explicit. For highly stiff PDEs like the Kuramoto--Sivashinsky equation, Kassam and Trefethen @KassamTrefethen2005 demonstrated the superiority of Exponential Time Differencing (ETD) and integrating factor methods, sparking a renewed interest in fourth-order time-stepping for stiff PDEs.

For equilibrium problems (Poisson, Helmholtz), the "matrix surgery" approach presented in this chapter is effective for moderate $N$. However, for large-scale simulations, efficiency becomes paramount. Haidvogel and Zang @HaidvogelZang1979 introduced a matrix diagonalisation technique that decoupled the system into independent one-dimensional problems, reducing the complexity of two-dimensional Poisson solves from $O(N^4)$ to $O(N^3)$. This was refined by Shen @Shen1995, who constructed basis functions satisfying the boundary conditions automatically, leading to sparse matrices and $O(N^2)$ complexity. The tensor product structure emphasised in this chapter --- exploiting Kronecker sums --- is a modern realisation of these early insights. The modern canonical treatment of these algorithms is found in the texts by Shen, Tang, and Wang @ShenTangWang2011 and Boyd @Boyd2000, which detail how to solve Poisson and Helmholtz equations in complex geometries using spectral-element and domain decomposition techniques.

The versatility of spectral methods is best seen in their application to specific physical systems, where preserving physical invariants is often as important as raw accuracy. For the Schrödinger equation, the preservation of unitarity (probability density) is critical in quantum mechanics. Taha and Ablowitz @Taha1984 compared various schemes for the nonlinear Schrödinger (NLS) equation, identifying the split-step Fourier method as a robust performer. Bao, Jin, and Markowich @BaoJinMarkowich2002 later provided a rigorous convergence analysis for time-splitting spectral approximations in the semiclassical regime, where highly oscillatory wavefunctions require the high resolution that only spectral methods can provide. For the Allen--Cahn and Cahn--Hilliard equations, energy stability is the primary concern. Eyre @Eyre1998 proposed the convex splitting idea, which ensures unconditional gradient stability. Shen and Yang @ShenYang2010 extended these ideas to spectral discretisations, proving that high-order spectral accuracy can be combined with thermodynamic consistency to accurately model phase separation dynamics.

The field continues to evolve, increasingly intersecting with machine learning and high-performance computing. The integration of spectral methods with deep learning has led to Fourier Neural Operators (FNO). Li _et al._ @LiFNO2021 introduced the FNO architecture, and subsequent works such as Cao _et al._ @Cao2025 use spectral layers to learn PDE solution operators that generalise across different discretisations. These "Spectral-Refiner" models fine-tune FNO predictions using physical constraints, merging the data-driven power of AI with the rigorous accuracy of spectral methods. Moving beyond the method of lines, Kaur _et al._ @Kaur2025 have proposed discretising time spectrally as well, treating the entire spacetime domain as a single tensor product grid, yielding exponential convergence in both space and time. Spectral methods are also uniquely suited for nonlocal operators (like the fractional Laplacian) which appear in anomalous diffusion models. Mustapha _et al._ @Mustapha2025 and Zhang and Wang @ZhangWang2025 have developed fast spectral solvers for these equations, exploiting the fact that nonlocal convolution operators become diagonal multiplications in the spectral domain. Recent implementations, such as those described by Melia _et al._ @Melia2026, leverage graphical processing units (GPUs) to accelerate hierarchical spectral solvers, enabling the simulation of turbulent flows and large-scale elliptic problems at resolutions previously thought impossible. This brief overview underscores that spectral methods are not a static collection of algorithms but a dynamic field adapting to new computational architectures and physical models.

== Summary <sec-spectral-pde-summary>

This chapter has developed spectral methods for the three fundamental classes of PDEs:

1. *Hyperbolic equations* (wave type):
   - Use explicit time stepping (leapfrog, Verlet)
   - Stability requires $Delta t = O(N^(-2))$ on Chebyshev grids
   - The wave equation demonstrates propagation, reflection, and interference

2. *Parabolic equations* (diffusion type):
   - Use implicit time stepping (Crank--Nicolson, backward Euler)
   - Avoid severe CFL restrictions; take time steps dictated by accuracy
   - The heat equation demonstrates smoothing and decay

3. *Elliptic equations* (equilibrium type):
   - No time stepping; solve a linear system directly
   - Kronecker structure enables efficient tensor-product methods
   - Poisson and Helmholtz equations arise throughout physics

Key computational tools:
- *Chebyshev differentiation via FFT* (`chebfft`): The $x arrow theta arrow$ Fourier pipeline transforms Chebyshev differentiation into FFT operations
- *Tensor products*: The 2D Laplacian $D_2 U + U D_2^top$ exploits separability for efficiency
- *Kronecker sums*: For implicit solvers, $D_(2,i) times.o I + I times.o D_(2,i)$ gives the discrete Laplacian

@tbl-chapter-summary provides a quick reference for all PDEs treated in this chapter.

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    table(
      columns: (0.8fr, 1.2fr, 1fr, 1fr, 1.2fr),
      align: (center, left, left, left, left),
      inset: (x: 0.6em, y: 0.5em),
      stroke: none,
      table.hline(stroke: 0.75pt + rgb("#142D6E")),
      table.header(
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Section*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*PDE*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Grid Type*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Time Scheme*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Key Operation*],
      ),
      table.hline(stroke: 0.5pt + luma(180)),
      [10.4], [1D wave], [Chebyshev], [Leapfrog], [`chebfft` twice],
      [10.5], [2D wave], [Cheb. tensor], [Leapfrog], [$D_2 U + U D_2^top$],
      [10.6], [1D heat], [Chebyshev], [Crank--Nicolson], [Implicit solve],
      [10.7], [2D heat], [Cheb. tensor], [Backward Euler], [Kronecker system],
      [10.8], [2D Poisson], [Cheb. tensor], [Direct solve], [Kronecker sum],
      [10.8], [2D Helmholtz], [Cheb. tensor], [Direct solve], [$bold(L) + k^2 I$],
      [10.9], [1D transport], [Fourier], [RK4], [FFT derivative],
    ),
  ),
  caption: [Summary of PDEs and methods presented in this chapter.],
) <tbl-chapter-summary>
