// textbook/chapters/time_stability.typ
// Chapter 17: Time Stepping, Stability, and the CFL Constraint
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Email: denys.dutykh@ku.ac.ae
// Homepage: https://www.denys-dutykh.com/
// Last modified: April 2026

#import "../styles/template.typ": dropcap, num, format-table

// Enable equation numbering for this chapter
#set math.equation(numbering: "(1)")

= Time Stepping, Stability, and the CFL Constraint <ch-time-stability>

#dropcap[Every spectral method we have developed so far excels in one dimension: space. Fourier and Chebyshev discretisations deliver exponential accuracy with remarkably few degrees of freedom, and the preceding chapters have assembled a formidable toolkit for boundary value problems, eigenvalue computations, and spatial differentiation. Yet the moment we couple a spectral spatial discretisation with a time-stepping scheme, a new and often devastating failure mode appears. A computation that runs beautifully at one time step can explode catastrophically when that step is increased by as little as ten percent. The culprit is the Courant--Friedrichs--Lewy (CFL) condition @Courant1928 @Courant1967, a stability constraint that ties the admissible time step to the spatial resolution. For spectral methods, the consequences are especially severe: the dense eigenvalue spectra of Fourier and Chebyshev differentiation matrices impose time-step restrictions that scale as $cal(O)(N^(-1))$, $cal(O)(N^(-2))$, or even $cal(O)(N^(-4))$ depending on the PDE and the discretisation. This chapter develops the geometric framework that explains these restrictions, makes them visible through eigenvalue plots and stability regions, and surveys the principal strategies for overcoming them: implicit methods, integrating factors, IMEX splitting, and several instructive alternatives from the finite-difference tradition. The central thesis is that success or failure in time stepping is governed by the relationship between the spectrum of the spatial operator and the stability region of the time integrator; once the reader internalises this geometric picture, the CFL condition ceases to be a mysterious formula and becomes a predictable, diagnosable, and often curable constraint.]

By the end of this chapter, you should be able to:

1. Recognise time-stepping instability in spectral computations and trace it to a violation of the CFL condition.
2. Write any spectral PDE discretisation in the method-of-lines form $dif bold(u) / dif t = L_N bold(u)$ and identify the relevant eigenvalues of $L_N$.
3. State the _stability rule of thumb_: a time-stepping scheme is stable when the scaled spectrum $Delta t dot sigma(L_N)$ lies inside the stability region of the time integrator.
4. Plot stability regions for forward Euler, classical RK4, leapfrog, Adams--Bashforth, Crank--Nicolson, and BDF methods, and superimpose spectral eigenvalues to predict stable time steps.
5. Explain why Chebyshev discretisations are stiffer than Fourier discretisations, identify the outlier eigenvalues that scale as $cal(O)(N^4)$ for second-derivative operators, and distinguish physical modes from boundary-localised spurious modes.
6. Describe when eigenvalue analysis fails for nonnormal operators and explain the role of pseudospectra.
7. Compare six strategies for overcoming the CFL barrier: implicit methods, integrating factors, IMEX splitting, Dufort--Frankel, Saulyev, and hyperbolisation; articulate the trade-offs of each.
8. Apply integrating-factor Runge--Kutta methods to stiff nonlinear PDEs and quantify the gain over naive explicit integration.

// ============================================================================
== When a Beautiful Spectral Computation Explodes <sec-catastrophe>
// ============================================================================

We begin not with theory but with catastrophe. The reader who has implemented the spectral PDE solvers from @ch-spectral-pde and @ch-fourier-pseudo already has working codes for wave propagation, heat diffusion, and nonlinear dispersive dynamics. In this section, we deliberately sabotage these codes by increasing the time step just past its stability threshold. The results are dramatic, instructive, and entirely predictable once the theory is in place.

=== Fourier advection gone wrong

Consider the periodic advection equation $u_t + u_x = 0$ on $[0, 2 pi)$ discretised with $N$ Fourier modes, as in @ch-fourier-pseudo. With RK4 time stepping, the computation is stable for $Delta t$ below a threshold that scales as $cal(O)(N^(-1))$. Increase $Delta t$ by 20 percent beyond this threshold and the solution develops sawtooth oscillations within a few time steps: the highest-frequency Fourier modes, those with wavenumber near $N\/2$, are the first to blow up. Within ten to twenty steps the oscillations have grown to fill the entire domain, and the computation is irretrievable.

=== Chebyshev wave equation gone wrong

The 1D wave equation $u_(t t) = u_(x x)$ on $[-1, 1]$ with Chebyshev collocation and leapfrog time stepping (@ch-spectral-pde) is stable for $Delta t = cal(O)(N^(-2))$. Increase $Delta t$ to $cal(O)(N^(-1))$ and the instability appears not as sawtooth oscillations across the domain, but as violent oscillations localised near the boundaries $x = plus.minus 1$. This boundary localisation is a signature of the _outlier eigenvalues_ of the Chebyshev second-derivative matrix, which we shall analyse in detail in @sec-chebyshev-stiffness.

=== A parabolic example

The heat equation $u_t = u_(x x)$ with an explicit scheme such as forward Euler and Chebyshev spatial discretisation faces the most severe restriction of all: $Delta t = cal(O)(N^(-4))$. For $N = 32$ this already forces time steps of order $10^(-6)$. Doubling $Delta t$ past the threshold produces an immediate blow-up, typically within one or two steps.

These three examples illustrate a universal pattern: _the instability always appears first in the modes associated with the largest eigenvalues of the spatial discretisation_, whether those are the high-frequency Fourier modes or the boundary-localised Chebyshev outliers. The question "Where did this time-step restriction come from?" has a precise geometric answer, which we now develop.

// ============================================================================
== Computational Étude 17.1: The Catastrophe Gallery <sec-etude-catastrophe>
// ============================================================================

This étude asks the reader to reproduce three deliberate blow-ups. Take the Fourier advection code from @ch-fourier-pseudo (Étude 11.2), the Chebyshev wave equation from @ch-spectral-pde (Étude 10.3), and a forward Euler discretisation of the Chebyshev heat equation, and in each case increase $Delta t$ by 20 percent beyond the stable threshold. Plot the onset of instability and identify where it first appears.

```python
import numpy as np
import matplotlib.pyplot as plt

def fourier_advection_blowup(N=64, dt_factor=1.2):
    """Fourier advection with deliberately unstable time step."""
    x = 2 * np.pi * np.arange(N) / N
    u = np.exp(-10 * (x - np.pi)**2)
    k = np.fft.fftfreq(N, d=1.0/N)
    dt_stable = 2.0 / N       # approximate CFL limit for RK4
    dt = dt_factor * dt_stable
    for step in range(80):
        u_hat = np.fft.fft(u)
        du = np.real(np.fft.ifft(1j * k * u_hat))
        u = u - dt * du        # forward Euler for simplicity
    return x, u

def cheb_wave_blowup(N=40, dt_factor=1.5):
    """Chebyshev wave equation with too-large time step."""
    j = np.arange(N + 1)
    x = np.cos(np.pi * j / N)
    u = np.exp(-30 * (x - 0.2)**2)
    dt_stable = 4.0 / N**2    # approximate CFL
    dt = dt_factor * dt_stable
    # Build D^2 via explicit matrix
    D = cheb_diff_matrix(N)
    D2 = D @ D
    D2 = D2[1:N, 1:N]         # strip boundaries
    u_int = u[1:N]
    u_prev = u_int.copy()
    for step in range(200):
        u_new = 2*u_int - u_prev + dt**2 * D2 @ u_int
        u_prev, u_int = u_int, u_new
    return x, u_int

def cheb_heat_blowup(N=20, dt_factor=2.0):
    """Chebyshev heat equation with forward Euler, unstable."""
    j = np.arange(N + 1)
    x = np.cos(np.pi * j / N)
    u = np.sin(np.pi * (1 + x) / 2)
    D = cheb_diff_matrix(N)
    D2 = D @ D
    D2 = D2[1:N, 1:N]
    dt_stable = 9.0 / N**4    # approximate CFL
    dt = dt_factor * dt_stable
    u_int = u[1:N]
    for step in range(20):
        u_int = u_int + dt * D2 @ u_int
    return x, u_int
```

The equivalent MATLAB implementation:

```matlab
function [x, u] = fourier_advection_blowup(N, dt_factor)
    x = 2*pi*(0:N-1)/N;
    u = exp(-10*(x - pi).^2);
    k = [0:N/2-1, 0, -N/2+1:-1];
    dt = dt_factor * 2.0 / N;
    for step = 1:80
        u_hat = fft(u);
        du = real(ifft(1i * k .* u_hat));
        u = u - dt * du;
    end
end
```

The Julia implementation:

```julia
function fourier_advection_blowup(; N=64, dt_factor=1.2)
    x = 2π * (0:N-1) / N
    u = exp.(-10 .* (x .- π).^2)
    k = vcat(0:N÷2-1, 0, -N÷2+1:-1)
    dt = dt_factor * 2.0 / N
    for step in 1:80
        u_hat = fft(u)
        du = real(ifft(1im .* k .* u_hat))
        u .= u .- dt .* du
    end
    return x, u
end
```

#figure(
  image("../figures/ch17/python/catastrophe_gallery.pdf", width: 95%),
  caption: [Three deliberate blow-ups. _Left_: Fourier advection with forward Euler and $Delta t$ 20% above the stability threshold; sawtooth oscillations appear at the highest wavenumbers and grow to fill the domain. _Centre_: Chebyshev wave equation with leapfrog and $Delta t$ 50% too large; the instability is localised near the boundaries $x = plus.minus 1$. _Right_: Chebyshev heat equation with forward Euler and $Delta t$ tripled beyond the stability limit; boundary-localised oscillations grow explosively within a few steps. In each case, the instability first appears in the modes associated with the largest eigenvalues of the spatial operator.],
) <fig-catastrophe-gallery>

The code generating @fig-catastrophe-gallery is available in:
- `codes/python/ch17/catastrophe_gallery.py`
- `codes/matlab/ch17/catastrophe_gallery.m`
- `codes/julia/ch17/catastrophe_gallery.jl`

=== Discussion

@fig-catastrophe-gallery makes three points that the rest of the chapter will formalise. First, instability is _sudden_: a small increase in $Delta t$ triggers explosive growth. Second, the _spatial location_ of the instability is diagnostic: sawtooth patterns across the domain for Fourier methods, boundary localisation for Chebyshev methods. Third, the severity of the restriction varies dramatically: $Delta t = cal(O)(N^(-1))$ for hyperbolic Fourier problems, $cal(O)(N^(-2))$ for Chebyshev wave equations, and $cal(O)(N^(-4))$ for Chebyshev diffusion. Understanding these scalings requires the method of lines and eigenvalue analysis, which we develop next.

// ============================================================================
== The Method of Lines as Master Viewpoint <sec-method-of-lines-ch17>
// ============================================================================

Every time-dependent PDE discretised spectrally in space becomes an ordinary differential equation (ODE) system. This perspective, the _method of lines_, is the key to understanding time-stepping stability.

=== From PDE to ODE system

Consider a general PDE of the form $u_t = cal(L) u + cal(N)(u)$, where $cal(L)$ is a linear spatial operator (such as $partial_(x x)$ or $c partial_x$) and $cal(N)$ is a nonlinear term. After spectral discretisation on $N$ grid points, we obtain the semidiscrete system
$ frac(dif bold(u), dif t) = L_N bold(u) + bold(N)(bold(u)), $ <eq-mol-ch17>
where $bold(u)(t) in bb(R)^N$ (or $bb(C)^N$) is the vector of nodal values, $L_N$ is the spectral discretisation of $cal(L)$, and $bold(N)$ is the discrete nonlinear operator. For second-order-in-time problems such as the wave equation, the system takes the form
$ frac(dif^2 bold(u), dif t^2) = L_N bold(u), $ <eq-mol-ch17-second>
which can be rewritten as a first-order system of twice the dimension.

The eigenvalues of $L_N$ (or, more precisely, of $Delta t dot L_N$) govern whether a given time-stepping scheme will be stable. This is the content of the next subsection.

=== The stability rule of thumb

#block(
  fill: rgb("#142D6E").lighten(92%),
  stroke: (left: 3pt + rgb("#142D6E")),
  inset: (left: 12pt, y: 10pt, right: 10pt),
  width: 100%,
)[
*Rule of Thumb* (Stability of the Method of Lines). A time-stepping method applied to the semidiscrete system @eq-mol-ch17 is stable if and only if every scaled eigenvalue $Delta t dot lambda_j$ of $Delta t dot L_N$ lies within the _stability region_ of the time integrator. Formally, if $sigma(L_N) = {lambda_1, dots, lambda_N}$ denotes the spectrum of $L_N$, then stability requires
$ Delta t dot lambda_j in cal(S) quad "for all" j = 1, dots, N, $ <eq-rule-of-thumb>
where $cal(S) subset.eq bb(C)$ is the stability region of the time-stepping formula.
]

This rule is _exact_ for linear, constant-coefficient problems and serves as an excellent guide for nonlinear and variable-coefficient problems via frozen-coefficient (local eigenvalue) analysis.

=== Stability regions of classical time integrators

The _stability region_ $cal(S)$ of a one-step method is the set of $z = Delta t dot lambda in bb(C)$ for which the numerical solution of the scalar test equation $u' = lambda u$ remains bounded. For the principal methods used in spectral computations:

*Forward Euler.* The scheme $u^(n+1) = u^n + Delta t lambda u^n = (1 + z) u^n$ has stability region
$ cal(S)_("FE") = {z in bb(C) : |1 + z| lt.eq.slant 1}, $ <eq-sr-fe>
a disk of radius 1 centred at $z = -1$.

*Classical Runge--Kutta (RK4).* The fourth-order Runge--Kutta method has the amplification factor
$ g(z) = 1 + z + frac(z^2, 2) + frac(z^3, 6) + frac(z^4, 24), $ <eq-sr-rk4>
and $cal(S)_("RK4") = {z : |g(z)| lt.eq.slant 1}$. This region extends approximately $2.8$ units along the negative real axis and about $2.8$ units along the imaginary axis, making RK4 suitable for both diffusive and advective problems.

*Leapfrog (for first-order systems).* The scheme $v^(n+1) - v^(n-1) = 2 z v^n$ is stable only on the open imaginary interval
$ cal(S)_("LF") = {z in bb(C) : z in (-2 i, 2 i)}, $ <eq-sr-leapfrog>
an infinitesimally thin segment on the imaginary axis. This makes leapfrog appropriate for purely advective problems (imaginary eigenvalues) but _never_ for diffusion (real negative eigenvalues).

*Leapfrog (for second-order systems).* The scheme $(v^(n+1) - 2 v^n + v^(n-1)) / (Delta t)^2 = lambda v^n$ gives the characteristic equation $g + g^(-1) = z + 2$, which requires $z in (-4, 0)$, i.e., the stability region is the open real interval from $-4$ to $0$ in the $z = lambda (Delta t)^2$ plane.

*Crank--Nicolson.* The scheme $u^(n+1) = u^n + frac(Delta t, 2)(lambda u^n + lambda u^(n+1))$ has amplification factor $g(z) = (1 + z\/2)\/(1 - z\/2)$. The stability region is the entire left half-plane:
$ cal(S)_("CN") = {z in bb(C) : "Re"(z) lt.eq.slant 0}, $ <eq-sr-cn>
making Crank--Nicolson unconditionally stable for problems with eigenvalues in the left half-plane, which includes all parabolic problems.

*Backward Euler.* The amplification factor is $g(z) = 1\/(1 - z)$, with stability region $cal(S)_("BE") = {z : |1 - z| gt.eq.slant 1}$, also containing the entire left half-plane and more. This strong stability comes at the cost of only first-order accuracy in time.

*Adams--Bashforth (second order).* The two-step explicit method $u^(n+1) = u^n + frac(Delta t, 2)(3 f^n - f^(n-1))$ has a stability region that is significantly smaller than RK4, extending only about $1$ unit along the negative real axis and $0.5$ units along the imaginary axis.

=== A starred remark on pseudospectra $star$

The rule of thumb @eq-rule-of-thumb is based on eigenvalue analysis and is exact for normal operators (those satisfying $L_N^* L_N = L_N L_N^*$). However, spectral differentiation matrices are typically _nonnormal_: the eigenvectors are far from orthogonal, and the resolvent $(z I - L_N)^(-1)$ can be very large even when $z$ is far from any eigenvalue. In such cases, even though $Delta t dot lambda_j in cal(S)$ for all $j$, the computation may exhibit large _transient growth_ before eventually decaying @TrefethenEmbree2005. The correct diagnostic for nonnormal operators is the _pseudospectrum_:
$ sigma_epsilon (L_N) = {z in bb(C) : ||(z I - L_N)^(-1)|| gt.eq.slant epsilon^(-1)}, $ <eq-pseudospectrum>
which enlarges each eigenvalue into a region whose size reflects the nonnormality. We explore this in the starred Étude 17.5. For most of this chapter, the eigenvalue rule of thumb is sufficient.

// ============================================================================
== Computational Étude 17.2: Stability Regions Meet Spectra <sec-etude-stability-regions>
// ============================================================================

This étude constructs the stability regions of six classical time integrators and superimposes the scaled eigenvalues $Delta t dot lambda_j$ of specific spectral operators. The visual overlap (or lack thereof) immediately predicts stability or instability.

```python
import numpy as np
import matplotlib.pyplot as plt

def stability_region_fe(x, y):
    """Forward Euler: |1 + z| <= 1."""
    z = x + 1j * y
    return np.abs(1 + z)

def stability_region_rk4(x, y):
    """RK4: |g(z)| <= 1 where g = 1 + z + z^2/2 + z^3/6 + z^4/24."""
    z = x + 1j * y
    g = 1 + z + z**2/2 + z**3/6 + z**4/24
    return np.abs(g)

def plot_stability_regions():
    """Plot stability regions for FE, RK4, AB2, and leapfrog."""
    x = np.linspace(-5, 2, 400)
    y = np.linspace(-4, 4, 400)
    X, Y = np.meshgrid(x, y)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # Forward Euler
    axes[0].contour(X, Y, stability_region_fe(X, Y), [1.0], colors='navy')
    axes[0].contourf(X, Y, stability_region_fe(X, Y), [0, 1.0], alpha=0.15)
    axes[0].set_title('Forward Euler')

    # RK4
    axes[1].contour(X, Y, stability_region_rk4(X, Y), [1.0], colors='navy')
    axes[1].contourf(X, Y, stability_region_rk4(X, Y), [0, 1.0], alpha=0.15)
    axes[1].set_title('Classical RK4')

    # Leapfrog (first-order): stability on imaginary axis only
    axes[2].plot([0, 0], [-2, 2], 'navy', lw=2.5)
    axes[2].set_title('Leapfrog (1st order)')

    for ax in axes:
        ax.set_xlabel(r'Re($z$)')
        ax.set_ylabel(r'Im($z$)')
        ax.set_aspect('equal')
        ax.axhline(0, color='k', lw=0.5)
        ax.axvline(0, color='k', lw=0.5)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig
```

The equivalent MATLAB implementation:

```matlab
function plot_stability_regions()
    [X, Y] = meshgrid(linspace(-5, 2, 400), linspace(-4, 4, 400));
    Z = X + 1i*Y;

    figure('Position', [100 100 1200 400]);

    % Forward Euler
    subplot(1,3,1);
    contour(X, Y, abs(1 + Z), [1 1], 'Color', [0.08 0.17 0.43], 'LineWidth', 1.5);
    title('Forward Euler'); axis equal; grid on;

    % RK4
    subplot(1,3,2);
    g = 1 + Z + Z.^2/2 + Z.^3/6 + Z.^4/24;
    contour(X, Y, abs(g), [1 1], 'Color', [0.08 0.17 0.43], 'LineWidth', 1.5);
    title('Classical RK4'); axis equal; grid on;

    % Leapfrog
    subplot(1,3,3);
    plot([0 0], [-2 2], 'Color', [0.08 0.17 0.43], 'LineWidth', 2.5);
    title('Leapfrog (1st order)'); axis equal; grid on;
end
```

The Julia implementation:

```julia
using LinearAlgebra, Plots

function stability_region_rk4(z)
    abs(1 + z + z^2/2 + z^3/6 + z^4/24)
end

function plot_stability_regions()
    x = range(-5, 2, length=400)
    y = range(-4, 4, length=400)
    Z = [xi + yi*im for yi in y, xi in x]

    p1 = contour(x, y, abs.(1 .+ Z), levels=[1.0], title="Forward Euler")
    p2 = contour(x, y, stability_region_rk4.(Z), levels=[1.0], title="RK4")
    plot(p1, p2, layout=(1,2), size=(900, 400), aspect_ratio=:equal)
end
```

#figure(
  image("../figures/ch17/python/stability_regions.pdf", width: 95%),
  caption: [Stability regions of six classical time integrators. _Top left_: Forward Euler (disk of radius 1 centred at $-1$). _Top centre_: Classical RK4 (extends $approx 2.8$ units along both the real and imaginary axes). _Top right_: Leapfrog for first-order systems (the open segment $(-2 i, 2 i)$ on the imaginary axis). _Bottom left_: Adams--Bashforth 2 (small region near the origin). _Bottom centre_: Crank--Nicolson (the entire left half-plane, shaded). _Bottom right_: Backward Euler (complement of the unit disk centred at $1$; contains the left half-plane and more). The shaded interiors are the stable regions.],
) <fig-stability-regions>

The code generating @fig-stability-regions is available in:
- `codes/python/ch17/stability_regions.py`
- `codes/matlab/ch17/stability_regions.m`
- `codes/julia/ch17/stability_regions.jl`

#figure(
  image("../figures/ch17/python/spectra_overlay.pdf", width: 95%),
  caption: [Stability regions with superimposed spectral eigenvalues. _Left_: The RK4 stability region with the purely imaginary eigenvalues $Delta t dot i k$ of the Fourier first-derivative matrix for $N = 32$ and two values of $Delta t$: stable (dots inside the region) and unstable (dots outside). _Right_: The forward Euler stability region with the real negative eigenvalues $Delta t dot lambda_j$ of the Chebyshev second-derivative matrix for $N = 16$. The outlier eigenvalues near $-0.048 N^4$ force an extremely small $Delta t$ to keep all dots inside the unit disk.],
) <fig-spectra-overlay>

The code generating @fig-spectra-overlay is available in:
- `codes/python/ch17/stability_regions.py`
- `codes/matlab/ch17/stability_regions.m`
- `codes/julia/ch17/stability_regions.jl`

=== Discussion

@fig-stability-regions makes the geometry of stability concrete. The reader should note two qualitative features. First, explicit methods have _bounded_ stability regions: there is always a maximum $|z|$ beyond which the method is unstable. This means that large eigenvalues $lambda_j$ force small $Delta t$. Second, implicit methods such as Crank--Nicolson and backward Euler have stability regions that contain the entire left half-plane (or more), so they can accommodate arbitrarily large negative real eigenvalues without any $Delta t$ restriction. The price, of course, is that each time step requires the solution of a linear system. @fig-spectra-overlay makes the connection to spectral methods explicit: the eigenvalue distribution of the spatial operator, scaled by $Delta t$, must fit inside the stability region. We now quantify this for Fourier and Chebyshev discretisations.

// ============================================================================
== Hyperbolic Problems and the Spectral CFL in Fourier Space <sec-fourier-cfl>
// ============================================================================

The cleanest setting for understanding the spectral CFL condition is periodic advection on $[0, 2 pi)$, discretised with Fourier collocation.

=== Eigenvalues of the Fourier differentiation matrix

The Fourier spectral derivative approximation of $partial_x$ on the equispaced grid $x_j = 2 pi j / N$ is equivalent to multiplication by $i k$ in Fourier space, where $k = 0, plus.minus 1, dots, plus.minus (N\/2 - 1), N\/2$. The eigenvalues of the $N times N$ Fourier differentiation matrix $D_F$ are therefore
$ lambda_k = i k, quad k = -N\/2 + 1, dots, N\/2 - 1, $ <eq-fourier-eigenvalues>
lying on the imaginary axis with $|lambda_k| lt.eq.slant N\/2$. The spectral radius is
$ rho(D_F) = N\/2 - 1 approx N\/2. $ <eq-fourier-spectral-radius>
This linear growth with $N$ is the mildest possible scaling among spectral differentiation operators.

=== The CFL condition for Fourier advection

For the advection equation $u_t + c u_x = 0$ with constant wave speed $c > 0$, the semidiscrete eigenvalues are $lambda_k = -i c k$. For the leapfrog scheme, whose stability region is the imaginary interval $(-2 i, 2 i)$, stability requires
$ |Delta t dot c dot k| lt.eq.slant 2 quad "for all" k, $
which gives
$ Delta t lt.eq.slant frac(2, c (N\/2 - 1)) approx frac(4, c N). $ <eq-fourier-cfl>
For RK4, the stability region extends approximately $2.83$ units along the imaginary axis, so the condition becomes $Delta t lt.eq.slant 2.83 / (c dot N\/2)$, i.e., $Delta t approx 5.66 / (c N)$. In either case, the critical time step scales as $cal(O)(N^(-1))$.

Trefethen @Trefethen2000 reports the empirical stability limit $Delta t < 1.9 / N$ for the Fourier advection code with RK4, which is consistent with this estimate once the specific stability region boundary is accounted for.

=== Variable coefficients and frozen-coefficient estimates

For variable-coefficient advection $u_t + c(x) u_x = 0$, the spatial operator is no longer diagonal in Fourier space. However, the _frozen-coefficient heuristic_ predicts that the maximum eigenvalue scales as $c_max dot N\/2$, where $c_max = max_x |c(x)|$. Trefethen @Trefethen2000 shows that for $c(x) in [6\/5, 1 + 3 pi \/ 2]$, the predicted threshold $Delta t lt.eq.slant (5\/6) N^(-1)$ agrees well with experiment for large $N$, though the observed limit $Delta t < 1.9 / N$ is slightly more permissive for small $N$.

// ============================================================================
== Computational Étude 17.3: The First Unstable Fourier Mode <sec-etude-fourier-cfl>
// ============================================================================

This étude computes the eigenvalues of the Fourier differentiation matrix, identifies the extreme modes, and compares the predicted CFL threshold with experimental observation.

```python
import numpy as np

def fourier_diff_eigenvalues(N):
    """Eigenvalues of the N-point Fourier differentiation matrix."""
    if N % 2 == 0:
        k = np.concatenate([np.arange(N//2), [-N//2+1+j for j in range(N//2)]])
        k = np.arange(-N//2+1, N//2+1)  # symmetric ordering
    return 1j * k

def find_cfl_threshold(N, method='rk4', c=1.0, tol=1e-3):
    """Binary search for the critical dt of Fourier advection."""
    dt_lo, dt_hi = 0.001 / N, 10.0 / N
    x = 2 * np.pi * np.arange(N) / N
    u0 = np.exp(np.cos(x))
    k = np.fft.fftfreq(N, d=1.0/N)

    for _ in range(50):
        dt = 0.5 * (dt_lo + dt_hi)
        u = u0.copy()
        stable = True
        for step in range(500):
            u_hat = np.fft.fft(u)
            rhs = lambda v: -c * np.real(np.fft.ifft(1j*k*np.fft.fft(v)))
            # RK4 step
            k1 = rhs(u)
            k2 = rhs(u + 0.5*dt*k1)
            k3 = rhs(u + 0.5*dt*k2)
            k4 = rhs(u + dt*k3)
            u = u + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
            if np.max(np.abs(u)) > 1e6:
                stable = False
                break
        if stable:
            dt_lo = dt
        else:
            dt_hi = dt
        if dt_hi - dt_lo < tol / N:
            break
    return 0.5 * (dt_lo + dt_hi)

# Verify O(N^{-1}) scaling
for N in [16, 32, 64, 128, 256]:
    dt_crit = find_cfl_threshold(N)
    print(f"N={N:4d}  dt_crit={dt_crit:.6f}  N*dt={N*dt_crit:.4f}")
```

The equivalent MATLAB implementation:

```matlab
function fourier_cfl_scaling()
    for N = [16, 32, 64, 128, 256]
        dt_crit = find_cfl_threshold(N);
        fprintf('N=%4d  dt_crit=%.6f  N*dt=%.4f\n', N, dt_crit, N*dt_crit);
    end
end

function dt_crit = find_cfl_threshold(N)
    x = 2*pi*(0:N-1)/N;
    u0 = exp(cos(x))';
    k = [0:N/2-1, 0, -N/2+1:-1]';
    dt_lo = 0.001/N;  dt_hi = 10/N;
    for iter = 1:50
        dt = 0.5*(dt_lo + dt_hi);
        u = u0;  stable = true;
        for step = 1:500
            rhs = @(v) -real(ifft(1i*k.*fft(v)));
            k1 = rhs(u);  k2 = rhs(u+0.5*dt*k1);
            k3 = rhs(u+0.5*dt*k2);  k4 = rhs(u+dt*k3);
            u = u + dt/6*(k1+2*k2+2*k3+k4);
            if max(abs(u)) > 1e6, stable = false; break; end
        end
        if stable, dt_lo = dt; else, dt_hi = dt; end
    end
    dt_crit = 0.5*(dt_lo + dt_hi);
end
```

The Julia implementation:

```julia
using FFTW

function find_cfl_threshold(N; c=1.0)
    x = 2π * (0:N-1) / N
    u0 = exp.(cos.(x))
    k = vcat(0:N÷2-1, 0, -N÷2+1:-1)
    dt_lo, dt_hi = 0.001/N, 10.0/N
    rhs(v) = -c * real(ifft(1im .* k .* fft(v)))
    for _ in 1:50
        dt = 0.5(dt_lo + dt_hi)
        u = copy(u0); stable = true
        for step in 1:500
            k1 = rhs(u); k2 = rhs(u .+ 0.5dt.*k1)
            k3 = rhs(u .+ 0.5dt.*k2); k4 = rhs(u .+ dt.*k3)
            u .+= dt/6 .* (k1 .+ 2k2 .+ 2k3 .+ k4)
            if maximum(abs.(u)) > 1e6; stable = false; break; end
        end
        stable ? (dt_lo = dt) : (dt_hi = dt)
    end
    return 0.5(dt_lo + dt_hi)
end
```

#figure(
  image("../figures/ch17/python/fourier_cfl_scaling.pdf", width: 85%),
  caption: [_Left_: Eigenvalues of the Fourier differentiation matrix for $N = 32$, lying on the imaginary axis from $-16 i$ to $16 i$. _Right_: Critical time step $Delta t_"crit"$ versus $N$ on a log-log scale for RK4 applied to Fourier advection with $c = 1$. The slope of $-1$ confirms the $cal(O)(N^(-1))$ scaling predicted by @eq-fourier-cfl. The dashed line shows $Delta t = 2.83 / (N\/2)$.],
) <fig-fourier-cfl-scaling>

The code generating @fig-fourier-cfl-scaling is available in:
- `codes/python/ch17/fourier_cfl.py`
- `codes/matlab/ch17/fourier_cfl.m`
- `codes/julia/ch17/fourier_cfl.jl`

=== Discussion

The log-log plot in @fig-fourier-cfl-scaling confirms that for Fourier advection, the critical time step scales exactly as $N^(-1)$, with a proportionality constant that matches the stability region of RK4 along the imaginary axis. The product $N dot Delta t_"crit"$ is approximately constant, which is the hallmark of a CFL condition. This $cal(O)(N^(-1))$ restriction is typically acceptable for hyperbolic problems: doubling $N$ merely halves the allowed $Delta t$, and the total work to reach a fixed time $T$ scales as $N^2$ (or $N^2 log N$ with FFT-based differentiation).

=== Variable-coefficient extension: the frozen-coefficient heuristic in practice

The constant-coefficient analysis above relied on the explicit eigenvalue formula $lambda_k = -i c k$. When the wave speed varies in space, the multiplication operator $c(x) partial_x$ is no longer diagonal in Fourier space, and exact eigenvalue formulas are unavailable. Instead, we must determine the critical time step _experimentally_ by running the time-stepping scheme and checking whether the solution remains bounded. The frozen-coefficient heuristic predicts that the CFL threshold is controlled by $c_max = max_x |c(x)|$:
$ Delta t_max approx frac(2.83, c_max dot N \/ 2). $ <eq-frozen-cfl>
We now test this prediction for the smooth speed profile $c(x) = 1 + 1\/2 sin(x)$, which has $c_max = 3\/2$ and mean $chevron.l c chevron.r = 1$.

The Python implementation performs a binary search over $Delta t$, declaring instability when $max |u| > 10^6$ after 500 RK4 steps:

```python
def variable_rhs(u, c, k):
    """RHS of u_t + c(x)*u_x = 0."""
    u_hat = np.fft.fft(u)
    du = np.real(np.fft.ifft(1j * k * u_hat))
    return -c * du

def is_stable(dt, u0, c, k, n_steps=500):
    """Run n_steps of RK4; return True if solution stays bounded."""
    u = u0.copy()
    for _ in range(n_steps):
        k1 = variable_rhs(u, c, k)
        k2 = variable_rhs(u + 0.5*dt*k1, c, k)
        k3 = variable_rhs(u + 0.5*dt*k2, c, k)
        k4 = variable_rhs(u + dt*k3, c, k)
        u = u + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        if np.max(np.abs(u)) > 1e6:
            return False
    return True

def find_dt_crit(N, c_func, tol=1e-6):
    """Binary search for critical dt."""
    x = 2*np.pi * np.arange(N) / N
    c = c_func(x);  u0 = np.exp(np.cos(x))
    k = np.fft.fftfreq(N, d=1.0/N)
    dt_lo, dt_hi = 0.1/(np.max(np.abs(c))*N), 10.0/N
    for _ in range(60):
        dt_mid = 0.5 * (dt_lo + dt_hi)
        if is_stable(dt_mid, u0, c, k):
            dt_lo = dt_mid
        else:
            dt_hi = dt_mid
        if dt_hi - dt_lo < tol * dt_lo:
            break
    return 0.5 * (dt_lo + dt_hi)
```

The equivalent MATLAB implementation:

```matlab
function rhs_val = variable_rhs(u, c, k)
    rhs_val = -c .* real(ifft(1i*k.*fft(u)));
end

function stable = is_stable_rk4(dt, u0, c, k, n_steps)
    u = u0; stable = true;
    for step = 1:n_steps
        k1 = variable_rhs(u,c,k); k2 = variable_rhs(u+0.5*dt*k1,c,k);
        k3 = variable_rhs(u+0.5*dt*k2,c,k); k4 = variable_rhs(u+dt*k3,c,k);
        u = u + dt/6*(k1+2*k2+2*k3+k4);
        if max(abs(u)) > 1e6, stable = false; return; end
    end
end
```

The Julia implementation:

```julia
function variable_rhs(u, c, k)
    u_hat = fft(u)
    du = real.(ifft(1im .* k .* u_hat))
    return -c .* du
end

function is_stable(dt, u0, c, k; n_steps=500)
    u = copy(u0)
    for _ in 1:n_steps
        k1 = variable_rhs(u, c, k)
        k2 = variable_rhs(u .+ 0.5dt.*k1, c, k)
        k3 = variable_rhs(u .+ 0.5dt.*k2, c, k)
        k4 = variable_rhs(u .+ dt.*k3, c, k)
        u .+= dt/6 .* (k1 .+ 2k2 .+ 2k3 .+ k4)
        maximum(abs.(u)) > 1e6 && return false
    end
    return true
end
```

#figure(
  image("../figures/ch17/python/fourier_cfl_variable.pdf", width: 85%),
  caption: [_Left_: Variable speed profile $c(x) = 1 + 1\/2 sin(x)$ with $c_max = 3\/2$ marked. _Right_: Critical time step $Delta t_"crit"$ versus $N$ on a log-log scale for RK4 applied to variable-coefficient Fourier advection. The dashed line shows the frozen-coefficient prediction $Delta t = 2.83 / (c_max dot N\/2)$; the dotted line shows the constant-coefficient reference $2.83 / (N\/2)$. The experimental thresholds converge to the frozen-coefficient prediction for large $N$, confirming the accuracy of the heuristic for this smooth speed profile.],
) <fig-fourier-cfl-variable>

The code generating @fig-fourier-cfl-variable is available in:
- `codes/python/ch17/fourier_cfl_variable.py`
- `codes/matlab/ch17/fourier_cfl_variable.m`
- `codes/julia/ch17/fourier_cfl_variable.jl`

=== Discussion

The frozen-coefficient prediction @eq-frozen-cfl matches the experimental threshold to within a few percent for $N gt.eq.slant 64$, confirming that for smooth, strictly positive $c(x)$, the heuristic $Delta t lt.eq.slant 2.83 / (c_max dot N\/2)$ is reliable. For small $N$, the experimental threshold is slightly more permissive than predicted: the frozen-coefficient estimate is conservative because the non-commutativity of multiplication by $c(x)$ and Fourier differentiation introduces coupling that slightly reduces the effective spectral radius below $c_max dot N\/2$. The key practical takeaway is that for variable-coefficient hyperbolic problems with Fourier discretisation, the frozen-coefficient CFL estimate provides a safe and sharp bound. The situation becomes dramatically worse for Chebyshev discretisations, as we now show.

// ============================================================================
== Chebyshev Discretisation, Stiffness, and Outliers <sec-chebyshev-stiffness>
// ============================================================================

This is the conceptual centre of the chapter. The Chebyshev differentiation matrices developed in @ch-chebyshev and used throughout @ch-spectral-pde have eigenvalue spectra that are qualitatively different from, and far more severe than, the Fourier case.

=== Eigenvalues of the Chebyshev second-derivative matrix

Consider the Chebyshev second-derivative operator $D_N^2$ on $[-1, 1]$ with homogeneous Dirichlet boundary conditions (i.e., the interior submatrix $tilde(D)_2$ obtained by stripping the first and last rows and columns, as in @ch-spectral-pde). The eigenvalues of this $(N-1) times (N-1)$ matrix are all real and negative, as expected for an approximation to $partial_(x x)$ with Dirichlet conditions.

However, the eigenvalue distribution is far from uniform. Trefethen @Trefethen2000 shows that approximately $2 \/ pi$ of the eigenvalues correspond to good approximations of the exact eigenvalues $-n^2 pi^2 \/ 4$ of $partial_(x x)$ on $[-1, 1]$ (these are the _physical modes_). The remaining eigenvalues are _outliers_: they are much larger in magnitude, growing as $cal(O)(N^4)$, and their corresponding eigenvectors are localised near the boundaries $x = plus.minus 1$.

More precisely, the largest eigenvalue in magnitude satisfies
$ |lambda_max| approx 0.048 N^4 $ <eq-cheb-outlier>
for large $N$. For $N = 60$, Trefethen reports $|lambda_max| approx 0.047438 N^4 approx 6.2 times 10^5$.

=== Physical modes versus boundary-localised outliers

The physical eigenvectors of $tilde(D)_2$ resemble $sin(n pi (x + 1) \/ 2)$ for small $n$, oscillating smoothly across $[-1, 1]$. The outlier eigenvectors, by contrast, are concentrated near $x = plus.minus 1$ and decay rapidly into the interior. This boundary localisation is a consequence of the Chebyshev point clustering: the grid spacing near $x = plus.minus 1$ is $cal(O)(N^(-2))$, so the effective spatial resolution there is much finer than in the interior, and the discrete operator "sees" modes that oscillate on the scale of the boundary grid spacing.

This distinction is crucial: the outlier eigenvalues do not correspond to any physically meaningful mode of $partial_(x x)$. They are artefacts of the spatial discretisation, yet they _completely determine the stability constraint_ for explicit time stepping.

=== Three scaling regimes

The eigenvalue growth rate of the spectral operator determines the CFL restriction. For the three principal PDE types:

#block(
  fill: rgb("#142D6E").lighten(92%),
  stroke: (left: 3pt + rgb("#142D6E")),
  inset: (left: 12pt, y: 10pt, right: 10pt),
  width: 100%,
)[
*Scaling regimes for explicit time stepping:*
+ *Fourier advection* ($partial_x$ on $[0, 2 pi)$): $rho(D_F) = cal(O)(N)$, so $Delta t = cal(O)(N^(-1))$.
+ *Chebyshev wave equation* ($partial_(x x)$ in the second-order form @eq-mol-ch17-second): the stability region is $(-4, 0)$ in the $lambda (Delta t)^2$ plane, with $|lambda_max| = cal(O)(N^4)$, giving $Delta t = cal(O)(N^(-2))$.
+ *Chebyshev heat equation* ($partial_(x x)$ in the first-order form @eq-mol-ch17): for forward Euler, $Delta t dot |lambda_max| lt.eq.slant 2$ requires $Delta t = cal(O)(N^(-4))$. This is the most severe restriction and makes naive explicit time stepping impractical for diffusion on Chebyshev grids with moderate to large $N$.
]

The $N^(-4)$ scaling for Chebyshev heat equations is especially striking: at $N = 32$, forward Euler requires $Delta t approx 9 times 10^(-6)$, while at $N = 64$ the constraint tightens to $Delta t approx 5.6 times 10^(-7)$. This is the primary motivation for the implicit and semi-implicit methods discussed in @sec-overcoming-cfl.

// ============================================================================
== Computational Étude 17.4: Physical Modes Versus Spurious Outliers <sec-etude-cheb-stiffness>
// ============================================================================

This étude computes the eigenvalues and eigenvectors of the Chebyshev second-derivative matrix $tilde(D)_2$ (with Dirichlet boundary conditions), plots the eigenvalue spectrum, compares physical and outlier eigenvectors, and verifies the $N^4$ scaling of the spectral radius.

```python
import numpy as np
import matplotlib.pyplot as plt

def cheb_diff_matrix(N):
    """Chebyshev differentiation matrix of size (N+1) x (N+1)."""
    if N == 0:
        return np.array([[0.0]])
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.ones(N + 1)
    c[0] = 2.0; c[N] = 2.0
    c *= (-1.0)**np.arange(N + 1)
    X = np.outer(x, np.ones(N + 1))
    dX = X - X.T + np.eye(N + 1)
    D = np.outer(c, 1.0 / c) / dX
    D -= np.diag(D.sum(axis=1))
    return D

def analyse_cheb_stiffness(N=32):
    """Eigenvalue analysis of Chebyshev D^2 with Dirichlet BCs."""
    D = cheb_diff_matrix(N)
    D2 = D @ D
    # Strip boundary rows/columns
    D2_int = D2[1:N, 1:N]

    eigenvalues, eigenvectors = np.linalg.eig(D2_int)
    idx = np.argsort(eigenvalues.real)  # sort by real part
    eigenvalues = eigenvalues[idx].real
    eigenvectors = eigenvectors[:, idx].real

    x = np.cos(np.pi * np.arange(N + 1) / N)
    x_int = x[1:N]

    # Physical mode (e.g., mode N//4) and outlier (last mode)
    phys_idx = N // 4
    outlier_idx = -1  # largest eigenvalue

    return eigenvalues, eigenvectors, x_int, phys_idx, outlier_idx

# Verify N^4 scaling
for N in [16, 24, 32, 48, 64]:
    D = cheb_diff_matrix(N)
    D2 = D @ D
    D2_int = D2[1:N, 1:N]
    lam_max = np.max(np.abs(np.linalg.eigvals(D2_int)))
    print(f"N={N:3d}  |lam_max|={lam_max:.2f}  ratio={lam_max/N**4:.6f}")
```

The equivalent MATLAB implementation:

```matlab
function analyse_cheb_stiffness(N)
    D = cheb_diff_matrix(N);
    D2 = D * D;
    D2_int = D2(2:N, 2:N);
    [V, Lambda] = eig(D2_int);
    lam = diag(Lambda);
    [lam, idx] = sort(real(lam));
    V = V(:, idx);
    fprintf('N=%3d  |lam_max|=%.2f  ratio=%.6f\n', ...
            N, abs(lam(end)), abs(lam(end))/N^4);
end
```

The Julia implementation:

```julia
using LinearAlgebra

function analyse_cheb_stiffness(N)
    D = cheb_diff_matrix(N)
    D2 = D * D
    D2_int = D2[2:N, 2:N]
    λ = eigvals(D2_int)
    λ_max = maximum(abs.(λ))
    println("N=$N  |λ_max|=$(round(λ_max, digits=2))  ratio=$(λ_max/N^4)")
    return λ
end
```

#figure(
  image("../figures/ch17/python/cheb_eigenvalues.pdf", width: 85%),
  caption: [Eigenvalues of the Chebyshev second-derivative matrix $tilde(D)_2$ (with Dirichlet boundary conditions) for $N = 32$. All eigenvalues are real and negative. The physical eigenvalues (approximately $2\/pi$ of the total) cluster near the exact values $-n^2 pi^2 \/ 4$ (shown as circles). The outlier eigenvalues extend to $|lambda| approx 0.048 N^4 approx 5 times 10^4$, far beyond the physical range.],
) <fig-cheb-eigenvalues>

#figure(
  image("../figures/ch17/python/cheb_eigenvectors.pdf", width: 90%),
  caption: [_Left_: A physical eigenvector of $tilde(D)_2$ (mode $N\/4 = 8$ for $N = 32$), oscillating smoothly across the interval and closely approximating $sin(8 pi (x+1)\/2)$. _Right_: The outlier eigenvector corresponding to the largest eigenvalue, localised near the boundaries $x = plus.minus 1$ and decaying rapidly into the interior. This boundary-localised mode has no physical counterpart yet completely determines the CFL constraint for explicit time stepping.],
) <fig-cheb-eigenvectors>

#figure(
  image("../figures/ch17/python/cheb_cfl_scaling.pdf", width: 70%),
  caption: [Spectral radius $rho(tilde(D)_2)$ versus $N$ on a log-log scale. The slope of $4$ confirms the $cal(O)(N^4)$ scaling @eq-cheb-outlier. The dashed line shows the fit $0.048 N^4$.],
) <fig-cheb-cfl-scaling>

The code generating @fig-cheb-eigenvalues, @fig-cheb-eigenvectors, and @fig-cheb-cfl-scaling is available in:
- `codes/python/ch17/chebyshev_stiffness.py`
- `codes/matlab/ch17/chebyshev_stiffness.m`
- `codes/julia/ch17/chebyshev_stiffness.jl`

=== Discussion

The eigenvalue analysis in @fig-cheb-eigenvalues reveals the structural origin of stiffness in Chebyshev spectral methods. The physical eigenvalues, those that approximate the true eigenvalues of $partial_(x x)$, are well-behaved and clustered. But the outliers, which grow as $N^4$, are entirely responsible for the devastating CFL constraint. @fig-cheb-eigenvectors shows that these outlier modes are boundary-localised and physically meaningless; they are artefacts of the $cal(O)(N^(-2))$ grid spacing near $x = plus.minus 1$. This explains why the instability in the Chebyshev blow-up (@fig-catastrophe-gallery, centre panel) appears at the boundaries: it is precisely these outlier modes that grow first when $Delta t$ exceeds the threshold. The $N^4$ scaling verified in @fig-cheb-cfl-scaling confirms that for explicit methods, doubling $N$ tightens the time-step restriction by a factor of $16$. This is the fundamental motivation for implicit, semi-implicit, and integrating-factor methods on Chebyshev grids.

// ============================================================================
== When Eigenvalues Are Not Enough $star$ <sec-pseudospectra-ch17>
// ============================================================================

This starred section addresses a subtlety that the reader may safely skip on first reading: the eigenvalue rule of thumb @eq-rule-of-thumb can fail for _nonnormal_ operators.

=== Nonnormal operators and transient growth

A matrix $A$ is normal if $A^* A = A A^*$; equivalently, if it has an orthonormal basis of eigenvectors. The Fourier differentiation matrix (a circulant) is normal, and so the eigenvalue analysis of @sec-fourier-cfl is exact. However, Chebyshev differentiation matrices with boundary conditions are typically nonnormal, sometimes highly so @TrefethenTrummer1987.

For nonnormal operators, the norm of the solution $bold(u)(t) = e^(L_N t) bold(u)(0)$ can grow substantially before eventually decaying, even when all eigenvalues of $L_N$ have negative real parts. This _transient growth_ is bounded by the _Kreiss constant_ and quantified by the pseudospectrum @TrefethenEmbree2005.

=== Pseudospectral contours as the right diagnostic

The $epsilon$-pseudospectrum $sigma_epsilon(L_N)$ is the set of all $z in bb(C)$ where the resolvent norm $||(z I - L_N)^(-1)||$ exceeds $epsilon^(-1)$, as defined in @eq-pseudospectrum. When the pseudospectral contours extend significantly beyond the eigenvalues, transient growth is possible, and the eigenvalue rule of thumb may be optimistic. For the Chebyshev first-derivative matrix with boundary conditions, the pseudospectra can bulge substantially into the right half-plane even though all eigenvalues have zero or negative real parts.

In practice, the pseudospectral warning is most relevant for _advection-dominated_ problems on Chebyshev grids, where the first-derivative operator is highly nonnormal. For second-derivative operators (diffusion), the matrix is closer to normal and the eigenvalue analysis is more reliable.

// ============================================================================
== Computational Étude 17.5: Same Eigenvalues, Different Behaviour $star$ <sec-etude-pseudospectra>
// ============================================================================

This étude compares the transient behaviour of two operators with identical eigenvalues but different degrees of nonnormality.

```python
import numpy as np
from scipy.linalg import expm, svd

def pseudospectra_demo(N=32):
    """Compare normal and nonnormal operators with same eigenvalues."""
    # Normal operator: diagonal with eigenvalues -1, -2, ..., -N
    Lambda = np.diag(-np.arange(1, N + 1, dtype=float))

    # Nonnormal operator: same eigenvalues, upper triangular perturbation
    A = Lambda.copy()
    for i in range(N - 1):
        A[i, i + 1] = N  # strong nonnormality

    # Compute resolvent norm on a grid
    x = np.linspace(-N - 5, 5, 200)
    y = np.linspace(-N/2, N/2, 200)
    X, Y = np.meshgrid(x, y)
    Z = X + 1j * Y

    sigma_normal = np.zeros_like(X)
    sigma_nonnormal = np.zeros_like(X)
    I = np.eye(N)
    for i in range(len(x)):
        for j in range(len(y)):
            z = Z[j, i]
            s_n = svd(z * I - Lambda, compute_uv=False)
            sigma_normal[j, i] = 1.0 / s_n[-1]
            s_a = svd(z * I - A, compute_uv=False)
            sigma_nonnormal[j, i] = 1.0 / s_a[-1]

    return X, Y, sigma_normal, sigma_nonnormal
```

The equivalent MATLAB implementation:

```matlab
function pseudospectra_demo(N)
    Lambda = diag(-(1:N));
    A = Lambda;
    for i = 1:N-1
        A(i, i+1) = N;
    end
    x = linspace(-N-5, 5, 200);
    y = linspace(-N/2, N/2, 200);
    [X, Y] = meshgrid(x, y);
    sig = zeros(size(X));
    I_mat = eye(N);
    for ii = 1:numel(X)
        z = X(ii) + 1i*Y(ii);
        s = svd(z*I_mat - A);
        sig(ii) = 1 / s(end);
    end
    contour(X, Y, log10(sig), -1:0.5:6);
end
```

The Julia implementation:

```julia
using LinearAlgebra

function pseudospectra_demo(; N=32)
    Λ = Diagonal(Float64.(-1:-1:-N))
    A = Matrix(Λ)
    for i in 1:N-1
        A[i, i+1] = N
    end
    x = range(-N-5, 5, length=200)
    y = range(-N÷2, N÷2, length=200)
    σ = [1.0 / svdvals(complex(xi, yi)*I - A)[end]
         for yi in y, xi in x]
    return x, y, σ
end
```

#figure(
  image("../figures/ch17/python/pseudospectra_comparison.pdf", width: 90%),
  caption: [Pseudospectra of two operators with identical eigenvalues $-1, -2, dots, -N$ (shown as dots). _Left_: A normal (diagonal) operator; the pseudospectral contours are concentric disks around each eigenvalue. _Right_: A highly nonnormal (upper triangular) operator; the pseudospectral contours bulge far beyond the eigenvalues into the right half-plane, indicating potential for large transient growth despite all eigenvalues being stable. The contour levels are $log_(10) epsilon^(-1)$ for $epsilon = 10^(-1), 10^(-2), dots$],
) <fig-pseudospectra-ch17>

The code generating @fig-pseudospectra-ch17 is available in:
- `codes/python/ch17/pseudospectra_demo.py`
- `codes/matlab/ch17/pseudospectra_demo.m`
- `codes/julia/ch17/pseudospectra_demo.jl`

=== Discussion

@fig-pseudospectra-ch17 delivers a visual warning: two operators with identical spectra can have vastly different pseudospectra and hence vastly different transient behaviour. For the normal operator, the pseudospectral contours hug the eigenvalues tightly; the eigenvalue rule of thumb is exact. For the nonnormal operator, the contours extend far into the right half-plane, meaning that the resolvent is large (and transient growth possible) even at points well outside the spectrum. In practice, this warns us that for advection-dominated problems on Chebyshev grids, the eigenvalue-based CFL prediction may be optimistic. The reader who wishes to pursue this topic in depth should consult Trefethen and Embree @TrefethenEmbree2005.

// ============================================================================
== How to Overcome the CFL Condition <sec-overcoming-cfl>
// ============================================================================

We now survey the principal strategies for relaxing or eliminating the CFL barrier. Each strategy has a clear "price to pay," and the goal of this section is to make those trade-offs explicit.

=== Implicit and semi-implicit methods <sec-implicit>

The most classical cure for CFL restrictions is to treat the stiff spatial operator _implicitly_. The backward Euler scheme applied to the semidiscrete heat equation $(dif bold(u)) / (dif t) = L_N bold(u)$ reads
$ frac(bold(u)^(n+1) - bold(u)^n, Delta t) = L_N bold(u)^(n+1), $
which gives
$ (I - Delta t dot L_N) bold(u)^(n+1) = bold(u)^n. $ <eq-backward-euler>
Since the stability region of backward Euler contains the entire left half-plane (and more), this scheme is unconditionally stable: _any_ $Delta t$ yields a bounded solution, regardless of the spectral radius of $L_N$.

The Crank--Nicolson scheme,
$ (I - frac(Delta t, 2) L_N) bold(u)^(n+1) = (I + frac(Delta t, 2) L_N) bold(u)^n, $ <eq-crank-nicolson-ch17>
combines unconditional stability with second-order accuracy in time. Both schemes require solving a linear system at each step. For one-dimensional Chebyshev problems, this system is dense (the Chebyshev $D^2$ matrix is full), but for Fourier problems the system is diagonal in Fourier space, making implicit methods nearly as fast as explicit ones.

Higher-order implicit methods include the _BDF_ (backward differentiation formulas) family, which achieve up to sixth-order accuracy while maintaining A-stability (for orders $lt.eq.slant 2$) or A($alpha$)-stability (for orders 3--6).

*The price:* Each time step requires a linear system solve ($cal(O)(N^2)$ for dense Chebyshev matrices, $cal(O)(N)$ for tridiagonal finite-difference or $cal(O)(N log N)$ for Fourier-diagonal systems). For nonlinear problems, Newton iteration may be needed.

=== Integrating factors and exponential time differencing <sec-integrating-factor-ch17>

For problems with a stiff linear part and a mild nonlinear part, the _integrating factor_ (IF) method @Lawson1967 removes the stiffness entirely by treating the linear term exactly. Consider the Fourier-space ODE for mode $k$:
$ hat(u)_k '(t) = lambda_k hat(u)_k (t) + hat(cal(N))_k (t), $
where $lambda_k$ is the (potentially stiff) linear symbol and $hat(cal(N))_k$ is the nonlinear term. Define the integrating factor variable $hat(v)_k = e^(-lambda_k t) hat(u)_k$; then
$ hat(v)_k '(t) = e^(-lambda_k t) hat(cal(N))_k (t), $ <eq-if-variable>
which is a _nonstiff_ equation that can be integrated with a standard explicit method such as RK4.

The per-step IF-RK4 algorithm, developed in detail in @ch-fourier-pseudo, applies the integrating factor over a single time step rather than cumulatively, avoiding numerical overflow. This approach was used for the KdV and Kuramoto--Sivashinsky computations in @ch-fourier-pseudo and is revisited in Étude 17.7.

_Exponential time differencing_ (ETD) methods @Cox2002 @KassamTrefethen2005 go further by incorporating the matrix exponential directly into the time-stepping formula, achieving high-order accuracy for the linear part while treating the nonlinear part explicitly. Kassam and Trefethen @KassamTrefethen2005 showed that ETD schemes of fourth order can dramatically outperform both standard explicit and IMEX methods for stiff PDEs such as the Kuramoto--Sivashinsky equation.

*The price:* The method is most natural when the linear operator is diagonal in the chosen spectral basis (e.g., Fourier modes for periodic problems). For Chebyshev problems, where $L_N$ is dense, computing the matrix exponential is more expensive.

=== IMEX splitting for semilinear PDEs <sec-imex>

For semilinear PDEs of the form
$ u_t = cal(L) u + cal(N)(u), $ <eq-imex-splitting>
the _implicit-explicit_ (IMEX) strategy @Ascher1995 @Ascher1997 treats the linear operator $cal(L)$ implicitly and the nonlinear operator $cal(N)$ explicitly:
$ (I - Delta t dot L_N) bold(u)^(n+1) = bold(u)^n + Delta t dot bold(N)(bold(u)^n). $ <eq-imex-scheme>
This removes the CFL restriction associated with the stiff linear part while avoiding the cost of solving a nonlinear system. In Fourier space, the implicit solve is diagonal and therefore trivial; in Chebyshev space, it requires a dense linear solve (though this can be precomputed and factored once).

IMEX schemes have become the workhorse of spectral PDE simulation. Higher-order IMEX Runge--Kutta methods @Ascher1997 combine implicit treatment of the linear part with explicit multi-stage treatment of the nonlinear part, achieving third- or fourth-order accuracy overall. The Allen--Cahn computation in @ch-fourier-pseudo (Étude 11.10) is a concrete example of IMEX in action.

*The price:* The implicit linear solve at each step, and the splitting error (the nonlinear part may limit accuracy even though it is not stiff).

=== A comparative interlude: Dufort--Frankel, Saulyev, and hyperbolisation <sec-comparative-interlude>

The methods above (implicit, IF, IMEX) are the most natural tools for global spectral discretisations. We now discuss three alternative approaches that arise most naturally from _local_ finite-difference stencils @Dutykh2016. They are presented here as _instructive neighbouring ideas_ rather than as standard spectral technology, because they teach valuable lessons about modified equations, stability-accuracy trade-offs, and the principle that stability is never obtained for free.

==== Dufort--Frankel

The classical leapfrog scheme for the heat equation $u_t = nu u_(x x)$ is _unconditionally unstable_ (as shown in @Dutykh2016). The Dufort--Frankel modification replaces $2 u_j^n$ in the centred-difference Laplacian by $u_j^(n-1) + u_j^(n+1)$:
$ frac(u_j^(n+1) - u_j^(n-1), 2 Delta t) = nu frac(u_(j-1)^n - (u_j^(n-1) + u_j^(n+1)) + u_(j+1)^n, (Delta x)^2). $ <eq-dufort-frankel>
This is _explicit_ (it can be solved directly for $u_j^(n+1)$) and _unconditionally stable_ by von Neumann analysis. However, the consistency error reveals a hidden hyperbolic perturbation: the scheme is consistent not with $u_t = nu u_(x x)$ but with
$ tau u_(t t) + u_t = nu u_(x x), quad "where" tau = nu (Delta t \/ Delta x)^2. $ <eq-df-modified>
This means that Dufort--Frankel achieves unconditional stability by implicitly adding a small second-derivative-in-time term that changes the character of the PDE from parabolic to hyperbolic. The scheme is first-order accurate in space if $Delta t = cal(O)((Delta x)^(3\/2))$ and second-order if $Delta t = cal(O)((Delta x)^2)$ @DuFort1953.

==== Saulyev

Saulyev's method @Saulyev1960 uses _asymmetric_ spatial differences in two alternating stages:

*Stage 1* (left to right): at each node $j$, the unknown $u_j^(n+1)$ appears on both sides of the equation, but the left neighbour $u_(j-1)^(n+1)$ is _already known_ from the previous $j-1$ step. Solved from left to right.

*Stage 2* (right to left): the sweep direction is reversed, using the right neighbour $u_(j+1)^(n+2)$ from the previous $j+1$ step. Solved from right to left.

Each individual stage introduces an $cal(O)(Delta t)$ directional bias, but after _both_ stages are completed the biases cancel, yielding a scheme that is unconditionally stable and second-order accurate (provided $Delta t = cal(O)((Delta x)^(3\/2))$). The scheme is fully explicit; no matrix inversions are needed.

==== Hyperbolisation

The _hyperbolisation_ approach @Chetverushkin2012 @Myshetskaya2015 modifies the PDE itself by adding a small $u_(t t)$ term:
$ tau u_(t t) + u_t = nu u_(x x), $ <eq-hyperbolised-heat>
where $tau > 0$ is a small parameter. The effect is to change the character of the equation from parabolic (infinite propagation speed) to hyperbolic (finite propagation speed $sqrt(nu \/ tau)$). For the resulting hyperbolic equation, an explicit leapfrog discretisation is stable under the relaxed CFL condition
$ frac(Delta t, Delta x) lt.eq.slant sqrt(frac(tau, nu)). $ <eq-hyp-cfl>
With the choice $tau = nu Delta x$, this gives $Delta t lt.eq.slant (Delta x)^(3\/2)$, which is weaker than the parabolic restriction $Delta t = cal(O)((Delta x)^2)$. The error introduced by the hyperbolic perturbation is $cal(O)(tau)$ and can be made arbitrarily small.

==== Editorial remark: instructive neighbours, not native spectral tools

These three methods share a common lesson: _stability is never obtained for free_. Dufort--Frankel trades consistency (adding a hidden $u_(t t)$ perturbation), Saulyev trades isotropy and temporal granularity (consistent only every two time steps), and hyperbolisation trades fidelity to the original PDE (adding an explicit $u_(t t)$ term). All three arise most naturally from local stencil manipulations on uniform grids, whereas global Fourier and Chebyshev collocation leads to dense spatial operators for which directional sweeps and local stencil tricks are less natural. Nevertheless, they provide valuable conceptual insight and serve as useful comparison points for the spectral methods that follow.

=== Stabilised explicit methods <sec-stabilised-explicit>

A more spectral-compatible approach to extending the stability of explicit methods is the _Runge--Kutta--Chebyshev_ (RKC) family @VanDerHouwen1980 @Verwer2004. The key idea is to construct an $s$-stage explicit Runge--Kutta method whose stability region extends approximately $0.65 s^2$ units along the negative real axis (compared to $approx 2$ for forward Euler and $approx 2.8$ for RK4). The number of stages $s$ is chosen to match the spectral radius: for a problem with $rho(L_N) = R$, one takes $s approx sqrt(R Delta t \/ 0.65)$, which gives an overall cost of $s$ function evaluations per step.

For the Chebyshev heat equation with $rho approx 0.048 N^4$, forward Euler requires $Delta t = cal(O)(N^(-4))$, while an $s$-stage RKC method with $s = cal(O)(N^2)$ stages allows $Delta t = cal(O)(1)$. The total cost per step is then $cal(O)(N^2)$ function evaluations, compared to the $cal(O)(N^4)$ steps (each costing $cal(O)(N)$) required by forward Euler to cover the same time interval. The net effect is a substantial speedup for problems where the stiffness is dominated by a few large negative real eigenvalues (precisely the Chebyshev outlier scenario).

Abdulle and Medovikov @Abdulle2002 further developed _second-order_ Chebyshev methods (ROCK2), achieving both extended stability and second-order accuracy.

*The price:* More stages per step (but each stage is just a matrix-vector product, not a linear system solve). The method is most effective when eigenvalues lie on or near the negative real axis; for problems with significant imaginary eigenvalue components (advection), the advantage is reduced.

// ============================================================================
== Computational Étude 17.6: A Fair Comparison of Cures <sec-etude-fair-comparison>
// ============================================================================

This is the signature étude of the chapter. We solve the 1D heat equation $u_t = nu u_(x x)$ on $[0, 2 pi)$ with periodic boundary conditions and smooth initial data, using six different time-stepping strategies. We compare: (1) the maximum stable $Delta t$, (2) the global error at a fixed time, (3) the wall-clock time, and (4) the coding complexity.

We use a _Fourier_ spatial discretisation (rather than Chebyshev) to isolate the time-stepping effects. In Fourier space, $u_t = nu u_(x x)$ becomes $hat(u)_k ' = -nu k^2 hat(u)_k$, with real negative eigenvalues $lambda_k = -nu k^2$ ranging from $0$ to $-nu (N\/2)^2$.

The six methods are:
+ *Forward Euler* (FE): explicit, first-order, CFL-limited.
+ *Backward Euler* (BE): implicit, first-order, unconditionally stable.
+ *Crank--Nicolson* (CN): implicit, second-order, unconditionally stable.
+ *Integrating factor + RK4* (IF-RK4): exact linear part, fourth-order nonlinear.
+ *Dufort--Frankel* (DF): explicit, unconditionally stable, modified equation.
+ *RKC* ($s$-stage): explicit, extended stability, $s$ stages per step.

```python
import numpy as np
import time

def heat_equation_comparison(N=64, nu=0.1, T=1.0):
    """Compare six time-stepping methods for Fourier heat equation."""
    k = np.fft.fftfreq(N, d=1.0/(2*np.pi*N/(2*np.pi)))
    k = np.fft.fftfreq(N, d=1.0/N)  # wavenumbers
    lam = -nu * k**2                  # eigenvalues
    x = 2 * np.pi * np.arange(N) / N
    u0 = np.sin(x) + 0.5 * np.cos(3*x)
    u0_hat = np.fft.fft(u0)

    # Exact solution in Fourier space
    u_exact_hat = u0_hat * np.exp(lam * T)
    u_exact = np.real(np.fft.ifft(u_exact_hat))

    results = {}

    # 1. Forward Euler
    dt_fe = 1.8 / (nu * (N/2)**2)     # near CFL limit
    nsteps = int(np.ceil(T / dt_fe))
    dt_fe = T / nsteps
    t0 = time.time()
    u_hat = u0_hat.copy()
    for n in range(nsteps):
        u_hat = u_hat + dt_fe * lam * u_hat
    elapsed = time.time() - t0
    err = np.max(np.abs(np.real(np.fft.ifft(u_hat)) - u_exact))
    results['FE'] = (dt_fe, nsteps, err, elapsed)

    # 2. Backward Euler
    dt_be = 0.01
    nsteps = int(np.ceil(T / dt_be))
    dt_be = T / nsteps
    t0 = time.time()
    u_hat = u0_hat.copy()
    for n in range(nsteps):
        u_hat = u_hat / (1 - dt_be * lam)
    elapsed = time.time() - t0
    err = np.max(np.abs(np.real(np.fft.ifft(u_hat)) - u_exact))
    results['BE'] = (dt_be, nsteps, err, elapsed)

    # 3. Crank-Nicolson
    dt_cn = 0.01
    nsteps = int(np.ceil(T / dt_cn))
    dt_cn = T / nsteps
    t0 = time.time()
    u_hat = u0_hat.copy()
    for n in range(nsteps):
        u_hat = u_hat * (1 + 0.5*dt_cn*lam) / (1 - 0.5*dt_cn*lam)
    elapsed = time.time() - t0
    err = np.max(np.abs(np.real(np.fft.ifft(u_hat)) - u_exact))
    results['CN'] = (dt_cn, nsteps, err, elapsed)

    # 4. Integrating Factor + RK4 (for linear problem, this is exact)
    dt_if = 0.01
    nsteps = int(np.ceil(T / dt_if))
    dt_if = T / nsteps
    t0 = time.time()
    u_hat = u0_hat.copy()
    E = np.exp(lam * dt_if)
    for n in range(nsteps):
        u_hat = E * u_hat  # exact for linear problem
    elapsed = time.time() - t0
    err = np.max(np.abs(np.real(np.fft.ifft(u_hat)) - u_exact))
    results['IF'] = (dt_if, nsteps, err, elapsed)

    return results, u_exact, x
```

The equivalent MATLAB implementation:

```matlab
function results = heat_equation_comparison(N, nu, T)
    k = [0:N/2-1, 0, -N/2+1:-1]';
    lam = -nu * k.^2;
    x = 2*pi*(0:N-1)'/N;
    u0 = sin(x) + 0.5*cos(3*x);
    u0_hat = fft(u0);
    u_exact_hat = u0_hat .* exp(lam * T);

    % Forward Euler
    dt = 1.8 / (nu * (N/2)^2);
    nsteps = ceil(T / dt);  dt = T / nsteps;
    u_hat = u0_hat;
    for n = 1:nsteps
        u_hat = u_hat + dt * lam .* u_hat;
    end
    err_fe = max(abs(real(ifft(u_hat - u_exact_hat))));

    % Crank-Nicolson
    dt = 0.01;  nsteps = ceil(T/dt);  dt = T/nsteps;
    u_hat = u0_hat;
    for n = 1:nsteps
        u_hat = u_hat .* (1 + 0.5*dt*lam) ./ (1 - 0.5*dt*lam);
    end
    err_cn = max(abs(real(ifft(u_hat - u_exact_hat))));
end
```

The Julia implementation:

```julia
using FFTW, BenchmarkTools

function heat_comparison(; N=64, ν=0.1, T=1.0)
    k = vcat(0:N÷2-1, 0, -N÷2+1:-1)
    λ = -ν .* k.^2
    x = 2π * (0:N-1) / N
    u0 = sin.(x) .+ 0.5cos.(3x)
    u0_hat = fft(u0)
    u_exact_hat = u0_hat .* exp.(λ .* T)

    # Forward Euler
    dt = 1.8 / (ν * (N÷2)^2)
    nsteps = ceil(Int, T / dt); dt = T / nsteps
    u_hat = copy(u0_hat)
    for n in 1:nsteps
        u_hat .+= dt .* λ .* u_hat
    end
    err_fe = maximum(abs.(real(ifft(u_hat .- u_exact_hat))))

    # Crank-Nicolson
    dt = 0.01; nsteps = ceil(Int, T/dt); dt = T/nsteps
    u_hat = copy(u0_hat)
    for n in 1:nsteps
        u_hat .= u_hat .* (1 .+ 0.5dt.*λ) ./ (1 .- 0.5dt.*λ)
    end
    err_cn = maximum(abs.(real(ifft(u_hat .- u_exact_hat))))

    return err_fe, err_cn
end
```

#figure(
  image("../figures/ch17/python/fair_comparison.pdf", width: 95%),
  caption: [A fair comparison of six time-stepping strategies for the Fourier heat equation with $N = 64$ and $nu = 0.1$. _Left_: Error versus $Delta t$ on a log-log scale. Forward Euler is limited to tiny $Delta t$ by the CFL condition, while CN and IF achieve much smaller errors at larger $Delta t$. _Right_: Error versus wall-clock time, showing the efficiency trade-off. IF-RK4 and CN dominate for this problem, achieving machine-precision-level accuracy at modest cost.],
) <fig-fair-comparison>

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    table(
      columns: (1.4fr, 0.9fr, 0.9fr, 0.9fr, 1.2fr),
      align: (left, center, center, center, left),
      inset: (x: 0.8em, y: 0.5em),
      stroke: none,
      table.hline(stroke: 0.75pt + rgb("#142D6E")),
      table.header(
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Method*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Order*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Stable $Delta t$*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Implicit?*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Price*],
      ),
      table.hline(stroke: 0.5pt + luma(180)),
      [Forward Euler], [1], [$cal(O)(N^(-2))$], [No], [Tiny $Delta t$],
      [Backward Euler], [1], [Any], [Yes], [Linear solve],
      [Crank--Nicolson], [2], [Any], [Yes], [Linear solve],
      [IF-RK4], [4], [Any], [No], [Matrix exp.],
      [Dufort--Frankel], [1--2], [Any], [No], [Modified PDE],
      [RKC ($s$ stages)], [1--2], [$cal(O)(s^2 / rho)$], [No], [$s$ evaluations],
    ),
  ),
  caption: [Summary of time-stepping strategies for diffusive problems. The "stable $Delta t$" column gives the CFL constraint for the Fourier heat equation with $N$ modes. "Any" means unconditionally stable (or exact). The "Price" column lists the main cost or compromise.],
) <tbl-method-comparison>

The code generating @fig-fair-comparison is available in:
- `codes/python/ch17/fair_comparison.py`
- `codes/matlab/ch17/fair_comparison.m`
- `codes/julia/ch17/fair_comparison.jl`

=== Discussion

@fig-fair-comparison deserves a careful reading because it encodes several lessons that are easy to miss at a glance.

*The left panel* (error versus $Delta t$) reveals the convergence order and the stability barrier of each method. Forward Euler (blue circles) is first-order: the error decreases linearly with $Delta t$, but the curve _terminates abruptly_ at $Delta t approx 2 times 10^(-2)$ because larger steps violate the CFL condition $Delta t lt.eq.slant 2 \/ (nu (N\/2)^2)$ and produce blow-up. No amount of patience can make forward Euler accurate at large $Delta t$; it simply ceases to exist. Backward Euler (green squares) has no stability barrier; the curve extends to large $Delta t$ without interruption. However, as a first-order method, its error decreases only as $cal(O)(Delta t)$, so achieving high accuracy requires small time steps regardless. At $Delta t = 10^(-4)$, backward Euler is roughly four orders of magnitude less accurate than Crank--Nicolson at the same step size. Crank--Nicolson (red triangles) is the first method to show the power of higher order: the error decreases as $cal(O)((Delta t)^2)$, a steeper slope on the log-log plot, and the method is unconditionally stable. This combination means that CN reaches machine precision at moderate $Delta t$ values (around $10^(-3)$), far larger than what forward Euler can tolerate. The integrating factor (orange diamonds) is exact for the linear heat equation: since all eigenvalues are treated via the matrix exponential $e^(lambda_k Delta t)$, the only error comes from floating-point arithmetic. Its error sits near $10^(-13)$ to $10^(-14)$ regardless of $Delta t$, a flat line at the bottom of the plot. For a _nonlinear_ problem, the IF-RK4 error would depend on $Delta t$ through the explicit treatment of the nonlinear term, but for the linear heat equation it is essentially exact.

*The right panel* (error versus wall-clock time) is the practitioner's plot: it answers the question "which method gives me the best accuracy for a given computational budget?" Here the ordering reshuffles. Forward Euler, despite being the cheapest per step (no linear solve), is the _least_ efficient overall: it spends enormous time taking tiny CFL-limited steps and still delivers only modest accuracy. Backward Euler is better (it can take larger steps) but its first-order accuracy limits the return on investment. Crank--Nicolson occupies a sweet spot: its second-order accuracy means that doubling the computational effort (halving $Delta t$) reduces the error by a factor of four. The integrating factor dominates completely for this linear problem, reaching machine precision in negligible time. The lesson is clear: _the cheapest method per step is not the cheapest method per digit of accuracy_. The implicit overhead of Crank--Nicolson (a diagonal division in Fourier space) is trivial compared to the thousands of extra time steps that forward Euler must take.

@tbl-method-comparison summarises these trade-offs. The reader should note that the comparison depends on the spatial discretisation: in Fourier space, implicit methods are diagonal and essentially free, which is why CN and IF dominate. On Chebyshev grids, where the spatial operator is dense, the cost of implicit solves is higher, and the balance shifts toward IMEX or stabilised explicit methods (RKC). The choice of time-stepping strategy should always be informed by the _structure_ of the spatial operator, not just its spectral radius.

// ============================================================================
== Nonlinear Capstone: The KdV Equation Revisited <sec-nonlinear-capstone>
// ============================================================================

We close the computational part of the chapter with a nonlinear problem that demonstrates the practical importance of the methods developed above. The Korteweg--de Vries (KdV) equation
$ u_t + 6 u u_x + u_(x x x) = 0 $ <eq-kdv-ch17>
was solved in @ch-fourier-pseudo (Étude 11.4) using integrating-factor RK4. Here we compare plain RK4 with IF-RK4 to quantify the gain.

=== Stiffness in nonlinear dispersive PDEs

The linear dispersive symbol $lambda_k = i k^3$ is purely imaginary, with $|lambda_k|$ growing as $k^3 = cal(O)(N^3)$ for the highest wavenumber. For plain RK4, the stability region extends approximately $2.83$ units along the imaginary axis, so the CFL condition is
$ Delta t lt.eq.slant frac(2.83, (N\/2)^3) = cal(O)(N^(-3)). $
At $N = 256$ this gives $Delta t approx 1.7 times 10^(-6)$, making long-time simulation prohibitively expensive.

The IF-RK4 method removes this restriction entirely by treating the linear dispersive term exactly. The remaining nonlinear term $-6 u u_x$ has mild derivatives and does not impose a severe time-step restriction. In practice, IF-RK4 allows time steps five to ten times larger than plain RK4, with the gain increasing as $N$ grows.

// ============================================================================
== Computational Étude 17.7: Plain RK Versus Integrating-Factor RK <sec-etude-kdv-comparison>
// ============================================================================

This étude runs the KdV equation with both plain RK4 and IF-RK4, comparing the maximum admissible time step and the accuracy of the solution.

```python
import numpy as np
from numpy.fft import fft, ifft

def kdv_plain_rk4(N=128, dt=None, T=1.0):
    """Solve KdV with plain RK4 (no integrating factor)."""
    x = 2 * np.pi * np.arange(N) / N
    k = np.fft.fftfreq(N, d=1.0/N)
    lam = 1j * k**3  # linear dispersion
    u = np.cos(np.pi * x)
    nsteps = int(np.ceil(T / dt))
    dt = T / nsteps

    def rhs(u_now):
        u_hat = fft(u_now)
        u_x = np.real(ifft(1j * k * u_hat))
        u_xxx = np.real(ifft(lam * u_hat))
        return -6 * u_now * u_x - u_xxx

    for n in range(nsteps):
        k1 = rhs(u)
        k2 = rhs(u + 0.5*dt*k1)
        k3 = rhs(u + 0.5*dt*k2)
        k4 = rhs(u + dt*k3)
        u = u + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        if np.max(np.abs(u)) > 1e10:
            return None  # blow-up
    return u

def kdv_ifrk4(N=128, dt=0.01, T=1.0):
    """Solve KdV with integrating-factor RK4."""
    x = 2 * np.pi * np.arange(N) / N
    k = np.fft.fftfreq(N, d=1.0/N)
    lam = 1j * k**3
    u = np.cos(np.pi * x)
    nsteps = int(np.ceil(T / dt))
    dt = T / nsteps

    E = np.exp(lam * dt / 2)
    E2 = E**2

    def nonlinear(u_hat):
        u_phys = np.real(ifft(u_hat))
        u_x = np.real(ifft(1j * k * u_hat))
        return fft(-6 * u_phys * u_x)

    u_hat = fft(u)
    for n in range(nsteps):
        Na = nonlinear(u_hat)
        a = E * u_hat + 0.5 * dt * E * Na
        Nb = nonlinear(a)
        b = E * u_hat + 0.5 * dt * E * Nb
        Nc = nonlinear(b)
        c = E2 * u_hat + dt * E * Nc
        Nd = nonlinear(c)
        u_hat = E2 * u_hat + dt/6 * (E2*Na + 2*E*(Nb+Nc) + Nd)

    return np.real(ifft(u_hat))
```

The equivalent MATLAB implementation:

```matlab
function u = kdv_ifrk4(N, dt, T)
    x = 2*pi*(0:N-1)'/N;
    k = [0:N/2-1, 0, -N/2+1:-1]';
    lam = 1i*k.^3;
    u_hat = fft(cos(pi*x));
    E = exp(lam*dt/2);  E2 = E.^2;
    nsteps = round(T / dt);
    NL = @(v) fft(-6*real(ifft(v)).*real(ifft(1i*k.*v)));
    for n = 1:nsteps
        Na = NL(u_hat);
        a = E.*u_hat + 0.5*dt*E.*Na;
        Nb = NL(a);
        b = E.*u_hat + 0.5*dt*E.*Nb;
        Nc = NL(b);
        c = E2.*u_hat + dt*E.*Nc;
        Nd = NL(c);
        u_hat = E2.*u_hat + dt/6*(E2.*Na + 2*E.*(Nb+Nc) + Nd);
    end
    u = real(ifft(u_hat));
end
```

The Julia implementation:

```julia
using FFTW

function kdv_ifrk4(; N=128, dt=0.01, T=1.0)
    x = 2π * (0:N-1) / N
    k = vcat(0:N÷2-1, 0, -N÷2+1:-1)
    λ = 1im .* k.^3
    u_hat = fft(cos.(π .* x))
    E = exp.(λ .* dt/2); E2 = E.^2
    nsteps = round(Int, T / dt)
    NL(v) = fft(-6 .* real(ifft(v)) .* real(ifft(1im .* k .* v)))
    for n in 1:nsteps
        Na = NL(u_hat)
        a = E .* u_hat .+ 0.5dt .* E .* Na
        Nb = NL(a)
        b = E .* u_hat .+ 0.5dt .* E .* Nb
        Nc = NL(b)
        c = E2 .* u_hat .+ dt .* E .* Nc
        Nd = NL(c)
        u_hat .= E2 .* u_hat .+ dt/6 .* (E2.*Na .+ 2E.*(Nb.+Nc) .+ Nd)
    end
    return real(ifft(u_hat))
end
```

#figure(
  image("../figures/ch17/python/kdv_rk_comparison.pdf", width: 95%),
  caption: [KdV equation: plain RK4 versus IF-RK4. _Left_: Maximum admissible $Delta t$ as a function of $N$. Plain RK4 (circles) scales as $N^(-3)$ due to the $i k^3$ dispersion; IF-RK4 (squares) is limited only by the nonlinear term and scales much more gently. _Right_: Solution at $t = 1$ for $N = 128$; both methods agree to graphical accuracy when run within their respective CFL limits, but IF-RK4 reaches $t = 1$ in approximately 100 steps versus $approx 10^4$ for plain RK4.],
) <fig-kdv-comparison>

The code generating @fig-kdv-comparison is available in:
- `codes/python/ch17/kdv_rk_comparison.py`
- `codes/matlab/ch17/kdv_rk_comparison.m`
- `codes/julia/ch17/kdv_rk_comparison.jl`

=== Discussion

@fig-kdv-comparison provides a compelling demonstration of the integrating-factor strategy in a nonlinear setting. By treating the stiff dispersive term $u_(x x x)$ exactly, IF-RK4 removes the $cal(O)(N^(-3))$ CFL restriction and allows time steps governed only by the mild nonlinearity. The computational savings are dramatic: at $N = 256$, IF-RK4 is roughly $100$ times faster than plain RK4 for the same accuracy. This is the genuinely spectral cure for stiffness in dispersive PDEs and illustrates why integrating factors and ETD methods have become the default time-stepping strategy for Fourier spectral simulations of nonlinear waves.

// ============================================================================
== Closing Synthesis: A Decision Guide <sec-decision-guide>
// ============================================================================

We close with a concise guide for choosing a time-stepping strategy in spectral computations:

- *Periodic transport* ($u_t + c u_x = 0$ on $[0, 2 pi)$): the CFL condition $Delta t = cal(O)(N^(-1))$ is mild. Explicit methods (RK4, leapfrog for advection) are perfectly reasonable. No implicit or integrating-factor treatment is needed.

- *Chebyshev discretisations of diffusive operators*: stiffness is dominated by the $cal(O)(N^4)$ outlier eigenvalues, and naive explicit stepping forces $Delta t = cal(O)(N^(-4))$. Use Crank--Nicolson, backward Euler, or BDF for linear problems. For nonlinear problems, IMEX is the standard choice.

- *Semilinear PDEs with stiff linear part* ($u_t = cal(L) u + cal(N)(u)$, with $cal(L)$ stiff): integrating-factor (IF-RK4), exponential time differencing (ETD), or IMEX methods are the most natural spectral tools. IF-RK4 is especially effective when $cal(L)$ is diagonal in the spectral basis (Fourier problems).

- *Stabilised explicit alternatives*: when implicit solves are prohibitive (e.g., multidimensional Chebyshev problems with non-separable operators), Runge--Kutta--Chebyshev (RKC) methods offer a viable explicit alternative with extended stability along the negative real axis.

- *Finite-difference alternatives* (Dufort--Frankel, Saulyev, hyperbolisation): these teach valuable lessons about modified equations and stability-accuracy trade-offs, but are not the native tools for global spectral discretisations. They are most useful as conceptual comparisons and in settings where local stencil methods are appropriate.

// ============================================================================
== A non-exhaustive literature overview
// ============================================================================

The study of time-stepping stability for partial differential equations begins with the landmark 1928 paper of Courant, Friedrichs, and Lewy @Courant1928, which established the principle that the numerical domain of dependence must contain the analytical domain of dependence for convergence. The English translation @Courant1967 brought this insight to a wider audience, and Thomée @Thomee2001 provides a comprehensive historical account tracing the development from the CFL condition through the modern theory of finite difference and finite element methods. The geometric framework of stability regions for ODE methods was formalised in the textbooks by Hairer, Nørsett, and Wanner @Hairer2009 and Hundsdorfer and Verwer @HundsdorferVerwer2003. The connection between stability regions and the spectra of spatial discretisation operators was made explicit by Gottlieb and Turkel @GottliebTurkel1980, who showed that the $cal(O)(N^4)$ outlier eigenvalues in Chebyshev second-derivative matrices impose severe time-step restrictions that require careful treatment.

Trefethen @Trefethen2000 elevated this analysis into a geometric art in Chapter 10 of _Spectral Methods in MATLAB_, which is the primary inspiration for much of this chapter. His presentation of stability regions, eigenvalue scaling, and the method-of-lines perspective set the standard for teaching spectral time stepping. The instability phenomenon in spectral methods was further analysed by Trefethen and Trummer @TrefethenTrummer1987, who demonstrated that the eigenvalues of spectral differentiation matrices can be extremely sensitive to boundary conditions and grid distribution. The theory of pseudospectra, which extends eigenvalue analysis to nonnormal operators, is developed in the monograph by Trefethen and Embree @TrefethenEmbree2005; this work provides the mathematical foundation for understanding why eigenvalue-based stability predictions can fail for advection operators on Chebyshev grids. Recent literature has further underscored that Chebyshev and Legendre operators are exceptionally sensitive to machine-precision perturbations: an infinitesimal rounding error on a grid of $N$ points can induce chaotic eigenvalue deviations scaling as $cal(O)(N^2)$, explaining the long-standing paradox where Legendre collocation, which theoretically enjoys an $cal(O)(N^(-1))$ time-step restriction for advection, frequently manifests a far more restrictive $cal(O)(N^(-5))$ constraint in double-precision arithmetic. The development of randomised algorithms for large-scale pseudospectral estimation, using randomised SVD and sketching techniques, has transformed pseudospectral analysis from a purely offline diagnostic into a scalable tool for assessing transient dynamics in high-dimensional systems.

The classical implicit countermeasure to spectral stiffness begins with the Crank--Nicolson scheme @Crank1947 and the backward differentiation formulas (BDF), which were systematised in the subsequent decades. For spectral methods in Fourier space, the implicit solve is diagonal and therefore trivial; on Chebyshev grids, the matrices are dense but can be pre-factored. The implicit-explicit (IMEX) approach, treating the stiff linear part implicitly and the nonlinear part explicitly, was formalised by Ascher, Ruuth, and Wetton @Ascher1995 and extended to Runge--Kutta frameworks by Ascher, Ruuth, and Spiteri @Ascher1997. IMEX methods have become the workhorse of spectral PDE simulation, as documented comprehensively by Hundsdorfer and Verwer @HundsdorferVerwer2003.

A significant recent advance in implicit time stepping is the development of _tunable_ BDF and IMEX schemes by Akrivis, Li, and Li @AkrivisLiLi2024. Classical high-order BDF methods (orders three to six) are constrained by the Dahlquist barrier: their stability domains shrink near the imaginary axis, limiting applicability to highly stiff nonlinear problems. The new schemes expand the solution at a shifted temporal evaluation point $t^(n + beta)$ with $beta > 1$ as an independently tunable parameter; as $beta$ increases, the stability domain expands significantly. For example, a fourth-order scheme with $beta = 5$ achieves the same expansive stability as the classical second-order BDF. This means that spectral practitioners are no longer forced into the detrimental compromise between severe time-step constraints of explicit high-order methods and the low accuracy of A-stable low-order implicit methods.

In parallel with linear BDF enhancements, the treatment of complex nonlinear gradient flows (Cahn--Hilliard equations, Navier--Stokes systems, phase-field crystal models) has been transformed by the _Scalar Auxiliary Variable_ (SAV) methodology @ShenXuYang2018. Standard IMEX schemes can violate the energy dissipation laws inherent to these PDEs, leading to nonlinear blow-up even when the linear CFL condition is satisfied. The SAV approach introduces an artificial scalar variable tracking the square root of the free energy functional, thereby reformulating the system so that a discrete analogue of the continuous energy dissipation law is preserved unconditionally. Adaptive IMEX BDF2 schemes have been rigorously shown to preserve modified energy dissipation under remarkably mild restrictions on the ratio of adjacent time steps, validating the use of aggressive temporal adaptivity for spectral computations of phase-field dynamics.

The integrating factor method originates with Lawson @Lawson1967, who proposed using the matrix exponential to transform away the stiff linear part. Cox and Matthews @Cox2002 developed exponential time differencing (ETD) methods that incorporate the matrix exponential directly into the time-stepping formula, and Kassam and Trefethen @KassamTrefethen2005 demonstrated that fourth-order ETD schemes dramatically outperform both standard explicit and IMEX methods for stiff PDEs such as the Kuramoto--Sivashinsky equation. Krogstad @KrogstadEtAl2005 further developed generalised integrating factor methods. Historically, ETD methods were largely confined to Fourier spectral discretisations on periodic domains, where the linear operator is perfectly diagonal. Recent parallelised ETD architectures, using hybrid MPI and OpenMP strategies, have extended exact linear integration to dense Chebyshev operators on non-periodic domains, unlocking ETD for multi-dimensional geometries with strict boundary condition enforcement. A further major advance is the formulation of _Exponential Spectral Deferred Correction_ (ESDC) methods @CaliariESDC2025, which bridge polynomial spectral deferred correction with exact exponential integration: iterative correction sweeps using exponential Runge--Kutta integrators achieve arbitrarily high order without complex multi-stage Butcher tableaus, and display significantly extended stability domains compared to existing exponential-linear multistep methods.

Stabilised explicit methods address the CFL bottleneck from a purely geometric perspective. The standard Runge--Kutta--Chebyshev (RKC) method, pioneered by van der Houwen and Sommeijer @VanDerHouwen1980 and developed by Verwer, Sommeijer, and Hundsdorfer @Verwer2004, stretches the stability domain to approximately $0.65 s^2$ units along the negative real axis using $s$-stage recursive relations based on shifted Chebyshev polynomials. Abdulle and Medovikov @Abdulle2002 extended this to second-order methods (ROCK2). However, classical RKC methods possess a critical vulnerability: their stability region remains infinitesimally thin along the imaginary axis, making them unsuitable for problems with significant advective or wave-like components. The state-of-the-art _Partitioned Runge--Kutta--Chebyshev_ (PRKC) methods @Zbinden2025 address this by independently partitioning the semi-discrete system into a stiff diffusion term and a non-stiff advection term, processing each with separately tuned RKC stage sequences. The resulting stability domain forms a customised rectangular region covering up to $0.65 s^2$ along the negative real axis and $2.15 m$ along the imaginary axis, simultaneously neutralising both Chebyshev boundary outliers and Fourier high-frequency modes within a single explicit step. For problems where severe stiffness is localised (geometric singularities, adaptive mesh refinement), multirate stabilised explicit schemes @AbdulleGrote2022 integrate stiff local terms with micro-steps while using stabilised macro-steps for the global dynamics, bypassing the temporal bottleneck without implicit solves.

The Dufort--Frankel scheme @DuFort1953, the Saulyev method @Saulyev1960, and the hyperbolisation approach @Chetverushkin2012 @Myshetskaya2015 @Dutykh2016 provide instructive alternatives from the finite-difference tradition. Dufort--Frankel achieves unconditional stability at the cost of a hidden hyperbolic perturbation; Saulyev provides unconditional stability through alternating directional sweeps; and hyperbolisation relaxes the CFL by adding a small $u_(t t)$ term to convert the parabolic equation into a hyperbolic system with finite propagation speed. The perturbation error analysis of the hyperbolisation approach was carried out rigorously by Dutykh @Dutykh2016, who showed that the error is $cal(O)(tau)$ where $tau$ is the hyperbolisation parameter.

The intense theoretical advancements in spectral time stepping are increasingly being embedded in sophisticated open-source frameworks that shift the burden of stability management from the practitioner to the solver architecture. The Dedalus framework @Burns2020, operating on symbolic graph representations of user-defined PDEs, automatically maps governing equations into optimal sparse spectral bases across Cartesian, cylindrical, and spherical geometries, and parses and partitions equations to treat stiff linear terms implicitly via sparse direct solvers while integrating nonlinear terms explicitly using multi-stage IMEX schemes. A landmark 2026 extension @Lecoanet2026 adds fast automated discrete adjoints by applying reverse-mode automatic differentiation directly to the discrete sparse computational graphs, eliminating the need for manual continuous adjoint derivations. The spectral/hp element framework Nektar++ @MoxeyNektar2020 provides a complementary approach for geometrically complex domains, combining domain decomposition with high-order polynomial expansions and employing sub-stepping velocity splitting schemes that decouple advection from diffusion-pressure operators.

The quest to overcome the strictly sequential nature of time integration has driven aggressive refinement of _Parallel-in-Time_ (PinT) algorithms. Gander @Gander2015 provides a comprehensive survey of 50 years of development in this area. The Parareal algorithm and the Parallel Full Approximation Scheme in Space and Time (PFASST) have seen significant adaptations for spectral discretisations: recent hybridisations deploy semi-implicit methods for a cheap coarse propagator and highly accurate ESDC methods for the fine corrector, yielding high-order, globally stable parallel-in-time integration uniquely equipped for stiff multiscale nonlinearities. Convergence analyses confirm that the parallel speedup scales almost linearly with the fine integrator's degrees of freedom. In a related but distinct direction, Deeb and Dutykh @Deeb2024 proposed a radical departure from standard time marching by integrating the Navier--Stokes equations using a stabilised Time Series Expansion (TSE) within a finite element framework: artificial diffusion terms at each temporal power damp spurious high-frequency oscillations while preserving consistency, and when coupled with Borel--Padé--Laplace resummation techniques, the operational radius of temporal integration is exponentially enlarged.

Finally, the precision and stability enhancements in spectral time stepping are actively advancing research in relativistic physics. The calculation of quasinormal modes (QNMs) for black holes and wormholes relies on resolving highly damped, complex eigenvalue spectra that describe gravitational ringdown. Batic and Dutykh @Batic2024 deployed fully resolved spectral methods to analyse electromagnetic and gravitational perturbations of noncommutative-geometry-inspired Schwarzschild black holes, detecting purely imaginary overdamped modes that were entirely missed by decades of WKB approximations. Subsequent studies @Batic2025 @Batic2025a revealed a previously unknown "Martini glass" morphology in the oscillatory spectrum of massive phantom wormholes and identified critical instability thresholds corresponding to the ratio of the Schwarzschild radius to the wormhole geometry; these findings powerfully underscore how modern spectral PDE frameworks are essential for extracting the fine-grained observational signatures required by quantum-gravity-inspired models.

The broader context of time-stepping strategies for spectral methods is covered in the monographs by Boyd @Boyd2000, Canuto _et al._ @Canuto2006, and Shen, Tang, and Wang @ShenTangWang2011. The overarching narrative of recent research is clear: the $cal(O)(N^4)$ Chebyshev outliers and the nonnormal pseudospectra are no longer insurmountable barriers. Modern numerical stability is achieved through hybridised, structurally aware, and increasingly automated frameworks; the geometric expansion of stability domains via tunable BDF and PRKC schemes, the energy dissipation guarantees of SAV methods, the analytical precision of exponential integrators, and the parallelism of PinT algorithms together ensure that for the modern practitioner, the CFL condition has evolved from an inescapable bottleneck into a manageable and well-understood geometric constraint.

// ============================================================================
== Summary <sec-stability-summary>
// ============================================================================

This chapter has developed the geometric framework that connects spectral spatial discretisations to time-stepping stability:

- *The method of lines* reduces every spectral PDE discretisation to an ODE system $dif bold(u) / dif t = L_N bold(u)$, whose stability is governed by the eigenvalues of $L_N$.

- *The stability rule of thumb*: a time-stepping method is stable when $Delta t dot sigma(L_N) subset cal(S)$, the stability region of the time integrator.

- *Fourier advection*: eigenvalues of $D_F$ lie on the imaginary axis with $|lambda_max| = cal(O)(N)$, giving $Delta t = cal(O)(N^(-1))$ for explicit methods. This is the mildest spectral CFL condition.

- *Chebyshev diffusion*: the second-derivative matrix $tilde(D)_2$ has _outlier eigenvalues_ scaling as $0.048 N^4$, localised near the boundaries. For forward Euler, $Delta t = cal(O)(N^(-4))$; for the second-order wave equation with leapfrog, $Delta t = cal(O)(N^(-2))$.

- *Pseudospectra* (starred): for nonnormal operators, eigenvalue analysis can be optimistic. The pseudospectrum $sigma_epsilon(L_N)$ provides a more reliable diagnostic, especially for advection on Chebyshev grids.

- *Overcoming the CFL barrier*: six strategies were compared:
  + _Implicit methods_ (BE, CN, BDF): unconditionally stable; require a linear solve per step.
  + _Integrating factors_ (IF-RK4): exact treatment of the linear part; most effective in Fourier space.
  + _IMEX splitting_: implicit for the linear part, explicit for the nonlinear part; the workhorse of spectral PDE simulation.
  + _Dufort--Frankel_: explicit and unconditionally stable, but modifies the PDE.
  + _Saulyev_: explicit and unconditionally stable via alternating sweeps; consistent every two steps.
  + _Hyperbolisation_: relaxes the CFL by adding a small $u_(t t)$ perturbation.
  + _RKC methods_: explicit with extended stability along the negative real axis.

- *Stability is never obtained for free*: every cure has a price, whether it is a linear solve, a modified equation, extra stages, or a restriction to diagonal operators.

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    table(
      columns: (0.6fr, 1.3fr, 1.1fr, 1.2fr),
      align: (center, left, left, left),
      inset: (x: 0.8em, y: 0.5em),
      stroke: none,
      table.hline(stroke: 0.75pt + rgb("#142D6E")),
      table.header(
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Étude*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Topic*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Key Insight*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Scaling*],
      ),
      table.hline(stroke: 0.5pt + luma(180)),
      [17.1], [Catastrophe gallery], [Instability appears at largest eigenvalue modes], [Visual],
      [17.2], [Stability regions], [Geometric view: spectrum must fit inside $cal(S)$], [$Delta t dot lambda_j in cal(S)$],
      [17.3], [Fourier CFL], [Imaginary eigenvalues $i k$, spectral radius $cal(O)(N)$], [$Delta t = cal(O)(N^(-1))$],
      [17.4], [Chebyshev stiffness], [Outlier eigenvalues $cal(O)(N^4)$, boundary modes], [$Delta t = cal(O)(N^(-4))$],
      [17.5], [Pseudospectra $star$], [Nonnormal operators: eigenvalues can mislead], [Resolvent norm],
      [17.6], [Fair comparison], [Six cures, each with a price], [Method-dependent],
      [17.7], [KdV: RK vs IF-RK], [IF removes $cal(O)(N^(-3))$ CFL restriction], [$approx 100 times$ speedup],
    ),
  ),
  caption: [Summary of computational études in this chapter.],
) <tbl-ch17-summary>

// ============================================================================
== Exercises <sec-time-stability-exercises>
// ============================================================================

*Exercise 17.1* (_Stability Region of the $theta$-Method_). The $theta$-method for $u' = lambda u$ is $u^(n+1) = u^n + Delta t [(1 - theta) lambda u^n + theta lambda u^(n+1)]$, with $theta in [0, 1]$. (a) Derive the amplification factor $g(z) = (1 + (1-theta) z) / (1 - theta z)$. (b) Show that the stability region is $|g(z)| lt.eq.slant 1$ and plot it for $theta = 0$ (forward Euler), $theta = 1\/2$ (Crank--Nicolson), $theta = 1$ (backward Euler), and $theta = 3\/4$. (c) Prove that the method is A-stable if and only if $theta gt.eq.slant 1\/2$.

*Exercise 17.2* (_Eigenvalue Scaling for Chebyshev First Derivative_). The Chebyshev first-derivative matrix $D_N$ (with the first and last rows and columns stripped for Dirichlet boundary conditions) has complex eigenvalues. (a) Compute the eigenvalues of $tilde(D)_1$ for $N = 16, 32, 64, 128$ and plot them in the complex plane. (b) Determine the spectral radius $rho(tilde(D)_1)$ and verify the scaling $rho = cal(O)(N^2)$. (c) Compare with the Fourier case, where $rho = cal(O)(N)$, and explain the difference in terms of the Chebyshev grid spacing near $x = plus.minus 1$.

*Exercise 17.3* (_Frozen-Coefficient Breakdown for Non-Smooth Speeds_). The étude accompanying @fig-fourier-cfl-variable demonstrated the frozen-coefficient heuristic for the smooth speed profile $c(x) = 1 + 1\/2 sin(x)$. This exercise explores its limitations. (a) Replace the speed by $c(x) = 1 + 1\/2 |sin(x)|$ (Lipschitz but not $C^1$) and repeat the CFL experiment for $N = 32, 64, 128, 256$. Does the frozen-coefficient prediction @eq-frozen-cfl still agree with the experimental threshold? (b) Try the multi-scale speed $c(x) = 2 + sin(x) + 0.8 sin(3x)$ and compare. (c) Find a smooth $c(x) > 0$ for which the frozen-coefficient prediction is _most_ conservative, i.e., the ratio $Delta t_("crit, exp") / Delta t_("crit, frozen")$ is maximised. Hint: consider large-amplitude oscillations in $c$.

*Exercise 17.4* (_Dufort--Frankel Modified Equation_). (a) Starting from the Dufort--Frankel scheme @eq-dufort-frankel, perform a Taylor expansion about $(x_j, t^n)$ and show that the leading truncation error includes the term $nu (Delta t \/ Delta x)^2 u_(t t)$. (b) Explain why this means the scheme is consistent with the telegraph equation @eq-df-modified rather than the heat equation. (c) Under what condition on $Delta t / Delta x$ does the hyperbolic perturbation become negligible?

*Exercise 17.5* (_Chebyshev Heat Equation: Explicit versus Implicit_). Solve $u_t = u_(x x)$ on $[-1, 1]$ with $u(plus.minus 1, t) = 0$ and $u(x, 0) = sin(pi(x+1)\/2)$ using: (a) forward Euler with $Delta t$ at the CFL limit, (b) Crank--Nicolson with $Delta t = 0.01$. For $N = 16, 32, 64$: (i) compare the number of time steps needed to reach $t = 1$, (ii) compare the error at $t = 1$ against the exact solution, (iii) measure the wall-clock time. Explain why Crank--Nicolson is vastly more efficient despite requiring a linear solve at each step.

*Exercise 17.6* (_RKC for a 2D Heat Equation_). Implement a Runge--Kutta--Chebyshev method with $s$ stages for the 2D heat equation on $[0, 2 pi)^2$ with Fourier spatial discretisation. (a) For $N = 32$, compute the spectral radius $rho = nu (N\/2)^2$ and determine the number of stages $s$ needed to use $Delta t = 0.01$. (b) Compare the cost per time step (number of FFTs) with forward Euler at the CFL limit and with Crank--Nicolson. (c) For which range of $N$ is RKC more efficient than Crank--Nicolson?

*Exercise 17.7* (_Hyperbolisation Parameter Selection_). For the hyperbolised heat equation @eq-hyperbolised-heat with $nu = 1$, $Delta x = 0.1$: (a) compute the CFL condition @eq-hyp-cfl for $tau = nu Delta x$, $tau = nu (Delta x)^2$, and $tau = 0.01$. (b) Solve the original and hyperbolised equations numerically with these three values of $tau$ and compare the solutions at $t = 0.5$. (c) Plot the error $||u_"hyp" - u_"heat"||_infinity$ as a function of $tau$ and verify the linear convergence predicted by the perturbation theory.

*Exercise 17.8* (_KdV Integrating Factor: Phase Accuracy_). Solve the KdV equation with a single-soliton initial condition $u(x, 0) = 2 "sech"^2(x - 10)$ on $[0, 40)$ using both plain RK4 and IF-RK4 with $N = 256$. (a) Run both methods to $t = 20$ and measure the phase error (shift in the soliton position relative to the exact travelling-wave solution). (b) Plot the phase error versus $Delta t$ for both methods. (c) Show that IF-RK4 achieves a given phase accuracy at roughly $1\/100$th the computational cost of plain RK4.
