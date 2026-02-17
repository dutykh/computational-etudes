// textbook/chapters/fourier_grids.typ
// Chapter 9: Physical and Fourier Space on Grids
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Email: denys.dutykh@ku.ac.ae
// Homepage: https://www.denys-dutykh.com/
// Last modified: February 2026

#import "../styles/template.typ": dropcap, num, format-table

// Enable equation numbering for this chapter
#set math.equation(numbering: "(1)")

= Physical and Fourier Space on Grids <ch-fourier-grids>

#dropcap[Every function tells two stories. In _physical space_, we see it as a curve: values $u(x)$ plotted against position $x$. In _Fourier space_, we see it decomposed into waves: amplitude coefficients $hat(u)(k)$ plotted against wavenumber $k$. These two representations are equivalent; each contains complete information about the function. The art of spectral methods lies in moving fluently between these perspectives, exploiting whichever view makes a problem simpler @Trefethen2000.]

This chapter develops the mathematical machinery connecting physical and Fourier space across three increasingly practical settings:

1. *Continuous functions on $RR$*: The classical Fourier transform, where both $x$ and $k$ range over the entire real line.

2. *Functions on an infinite grid $h ZZ$*: The semidiscrete Fourier transform, where sampling at spacing $h$ restricts wavenumbers to a bounded interval.

3. *Functions on a periodic grid*: The discrete Fourier transform (DFT), computable via the Fast Fourier Transform (FFT), which is the workhorse of practical spectral computation.

Along the way, we encounter the phenomenon of _aliasing_, which explains why high frequencies masquerade as low frequencies on discrete grids, and the _sinc function_, which provides the bridge between discrete samples and continuous band-limited interpolants. These concepts form the foundation for everything that follows in spectral methods.

== Motivation: Two Views of the Same Function <sec-two-views>

Consider a smooth function such as $u(x) = e^(-x^2) cos(3x)$, a Gaussian-modulated cosine. @fig-two-views shows this function from both perspectives.

#figure(
  image("../figures/ch09/python/two_views_function.pdf", width: 95%),
  caption: [Two views of the same function $u(x) = e^(-x^2) cos(3 x)$. _Left_: Physical space representation showing the oscillating, localized wavepacket. _Right_: Fourier space representation showing the magnitude of the Fourier transform on a logarithmic scale. The rapid decay of $|hat(u)(k)|$ reflects the smoothness of $u(x)$.],
) <fig-two-views>

In physical space (left panel), we see a wavepacket: oscillations modulated by a Gaussian envelope. The function is smooth, localized, and decays rapidly as $|x| arrow infinity$.

In Fourier space (right panel), we see the same information differently. The Fourier transform $hat(u)(k)$ shows two peaks near $k = plus.minus 3$, corresponding to the carrier frequency $cos(3x)$. The transform decays rapidly as $|k| arrow infinity$, a hallmark of smooth functions @Trefethen2000 @Gegenbauer2025.

The following Python code approximates the Fourier transform using the FFT on a large computational domain:

```python
import numpy as np

def compute_spectrum(u_func, L=20.0, N=1024):
    """Approximate Fourier transform via FFT on [-L/2, L/2]."""
    x = np.linspace(-L/2, L/2, N, endpoint=False)
    dx = x[1] - x[0]
    u = u_func(x)

    # Approximate continuous FT via FFT
    k = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    u_hat = dx * np.fft.fft(u)
    return k, u_hat
```

The equivalent MATLAB implementation:

```matlab
function [k, u_hat] = compute_spectrum(u_func, L, N)
    % Approximate Fourier transform via FFT on [-L/2, L/2]
    x = linspace(-L/2, L/2, N+1); x(end) = [];
    dx = x(2) - x(1);
    u = u_func(x);

    % Wavenumber vector (MATLAB ordering)
    k = 2*pi * [0:N/2-1, -N/2:-1] / (N*dx);
    u_hat = dx * fft(u);
end
```

The Julia implementation:

```julia
using FFTW

function compute_spectrum(u_func; L=20.0, N=1024)
    x  = range(-L/2, stop=L/2, length=N+1)[1:N]
    dx = x[2] - x[1]
    u  = u_func.(x)

    # Approximate continuous FT via FFT
    k     = 2π .* fftfreq(N, 1/dx)
    u_hat = dx .* fft(u)
    return collect(x), u, k, u_hat
end
```

The code generating @fig-two-views is available in:
- `codes/python/ch09/two_views_function.py`
- `codes/matlab/ch09/two_views_function.m`
- `codes/julia/ch09/two_views_function.jl`

This visual preview illustrates the central theme: _smooth functions have rapidly decaying Fourier transforms_. This connection between physical smoothness and spectral decay is what makes spectral methods so powerful.

== The Continuous Fourier Transform on $RR$ <sec-continuous-ft>

=== Definitions and Conventions

The _Fourier transform_ of a function $u(x)$ defined on $RR$ is
$ hat(u)(k) = integral_(-infinity)^(infinity) e^(-i k x) u(x) dif x, quad k in RR. $ <eq-fourier-transform>
The quantity $hat(u)(k)$ represents the _amplitude density_ of $u$ at wavenumber $k$. This definition goes back to Fourier's foundational treatise on heat conduction @Fourier1822. The process of decomposing a function into its constituent waves is called _Fourier analysis_.

Conversely, we can reconstruct $u$ from $hat(u)$ by the _inverse Fourier transform_:
$ u(x) = frac(1, 2 pi) integral_(-infinity)^(infinity) e^(i k x) hat(u)(k) dif k, quad x in RR. $ <eq-inverse-ft>
This reconstruction, building $u$ from a superposition of complex exponentials $e^(i k x)$, is called _Fourier synthesis_.

The variable $x$ is the _physical variable_ (position), and $k$ is the _Fourier variable_ (wavenumber). For a wave $e^(i k x)$, the wavenumber $k$ determines the spatial frequency: the number of oscillations per unit length is $|k| \/ (2 pi)$.

=== Fundamental Properties

The Fourier transform satisfies several important properties that make it a powerful tool for analysis. @tbl-ft-properties summarizes the key relationships.

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    table(
      columns: (1.2fr, 2fr, 2fr),
      align: (left, center, center),
      inset: (x: 1em, y: 0.6em),
      stroke: none,
      table.hline(stroke: 0.75pt + rgb("#142D6E")),
      table.header(
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Property*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Physical Space*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Fourier Space*],
      ),
      table.hline(stroke: 0.5pt + luma(180)),
      [Linearity], [$a u + b v$], [$a hat(u) + b hat(v)$],
      [Translation], [$u(x - x_0)$], [$e^(-i k x_0) hat(u)(k)$],
      [Modulation], [$e^(i k_0 x) u(x)$], [$hat(u)(k - k_0)$],
      [Differentiation], [$u'(x)$], [$i k hat(u)(k)$],
      [Convolution], [$u * v$], [$hat(u) dot hat(v)$],
      [Multiplication], [$u dot v$], [$hat(u) * hat(v)$],
    ),
  ),
  caption: [Fundamental properties of the Fourier transform. Operations in physical space correspond to simple algebraic operations in Fourier space.],
) <tbl-ft-properties>

The _differentiation property_ is particularly important for spectral methods: differentiation in physical space becomes multiplication by $i k$ in Fourier space. This simple observation underlies the efficiency of spectral methods for differential equations @Trefethen2000 @Fourier1822.

=== Physical and Fourier Variables: A Summary

@tbl-physical-fourier summarizes the relationship between physical and Fourier variables across the three settings we develop in this chapter. The key insight is that discreteness in one domain implies boundedness in the other.

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    table(
      columns: (1.5fr, 1.5fr, 1.5fr),
      align: (left, center, center),
      inset: (x: 1em, y: 0.6em),
      stroke: none,
      table.hline(stroke: 0.75pt + rgb("#142D6E")),
      table.header(
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Setting*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Physical Variable*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Fourier Variable*],
      ),
      table.hline(stroke: 0.5pt + luma(180)),
      [Continuous functions], [$x in RR$], [$k in RR$],
      [Infinite grid $h ZZ$], [$x_j = j h$], [$k in [-pi\/h, pi\/h]$],
      [Periodic grid, $N$ points], [$x_j in [0, 2 pi)$], [integer $k in ZZ$],
    ),
  ),
  caption: [Physical versus Fourier variables for three settings. Discreteness in physical space implies boundedness in Fourier space; periodicity in physical space implies discreteness in Fourier space.],
) <tbl-physical-fourier>

== Sampling on an Infinite Grid: The Semidiscrete Fourier Transform <sec-semidiscrete>

=== The Infinite Grid $h ZZ$

We now consider sampling a continuous function at equally spaced points. The _infinite grid_ $h ZZ$ consists of points
$ x_j = j h, quad j in ZZ, $
where $h > 0$ is the _grid spacing_ or _mesh width_.

#figure(
  align(center)[
    #box(width: 80%)[
      #align(center)[
        #box[
          #place(dy: 3pt)[#line(length: 100%, stroke: 0.5pt)]
          #grid(
            columns: (1fr,) * 11,
            align: center,
            ..range(11).map(i => [
              #circle(radius: 3pt, fill: black)
            ])
          )
        ]
        #v(0.3em)
        #grid(
          columns: (1fr,) * 11,
          align: center,
          ..range(11).map(i => [
            #text(size: 0.8em)[$x_(#(i - 5))$]
          ])
        )
        #v(0.5em)
        #text(size: 0.9em)[Grid spacing $h$]
      ]
    ]
  ],
  caption: [The infinite grid $h ZZ$ with grid points $x_j = j h$ for $j in ZZ$.],
) <fig-infinite-grid>

A grid function $v$ assigns values $v_j$ to each grid point $x_j$. We want to define a Fourier transform for such grid functions.

=== Aliasing: Why Wavenumbers Are Bounded

The crucial observation is that on the discrete grid, not all wavenumbers are distinguishable. Consider two complex exponentials $f(x) = e^(i k_1 x)$ and $g(x) = e^(i k_2 x)$. On the continuum $RR$, these are different functions whenever $k_1 eq.not k_2$. But when restricted to the grid $h ZZ$:
$ f_j = e^(i k_1 x_j) = e^(i k_1 j h), quad g_j = e^(i k_2 x_j) = e^(i k_2 j h). $
If $k_1 - k_2$ is an integer multiple of $2 pi \/ h$, then $f_j = g_j$ for all $j$!

This phenomenon is called _aliasing_: waves with wavenumbers differing by multiples of $2 pi \/ h$ are indistinguishable on the grid. They are _aliases_ of each other @Shannon1949 @Jerri1977.

#block(
  fill: rgb("#142D6E").lighten(92%),
  stroke: (left: 3pt + rgb("#142D6E")),
  inset: (left: 12pt, y: 10pt, right: 10pt),
  width: 100%,
)[
*Key Insight*: On a grid with spacing $h$, wavenumbers differing by $2 pi \/ h$ are indistinguishable:
$ e^(i k x_j) = e^(i (k + 2 pi m \/ h) x_j) quad "for all" j in ZZ, m in ZZ. $
Consequently, it suffices to measure wavenumbers in an interval of length $2 pi \/ h$. By convention, we choose the symmetric interval $[-pi\/h, pi\/h]$, called the _Nyquist interval_ @Shannon1949 @Jerri1977.
]

=== The Semidiscrete Fourier Transform

For a grid function $v$ defined on $h ZZ$, the _semidiscrete Fourier transform_ is
$ hat(v)(k) = h sum_(j = -infinity)^(infinity) e^(-i k x_j) v_j, quad k in [-pi\/h, pi\/h]. $ <eq-semidiscrete-ft>
The inverse transform reconstructs the grid values:
$ v_j = frac(1, 2 pi) integral_(-pi\/h)^(pi\/h) e^(i k x_j) hat(v)(k) dif k, quad j in ZZ. $ <eq-semidiscrete-inverse>

Note the duality: the forward transform @eq-semidiscrete-ft is a sum (discrete in physical space), while the inverse @eq-semidiscrete-inverse is an integral (continuous in Fourier space, but bounded).

These formulas approximate the continuous Fourier transform and its inverse: @eq-semidiscrete-ft is a trapezoid rule approximation to @eq-fourier-transform, and @eq-semidiscrete-inverse truncates the integration domain of @eq-inverse-ft. As $h arrow 0$, both pairs of formulas converge.

=== Computational Étude 1: Aliasing of $sin(pi x)$ and $sin(9 pi x)$ <sec-etude-aliasing>

@fig-aliasing demonstrates aliasing concretely. Consider the two functions $sin(pi x)$ and $sin(9 pi x)$ sampled on the grid $h = 1\/4$ (i.e., $1\/4 ZZ$). Despite being completely different continuous functions, they produce _identical_ samples at every grid point!

#figure(
  image("../figures/ch09/python/aliasing_demo.pdf", width: 85%),
  caption: [Aliasing on the grid $h = 1\/4$. The functions $sin(pi x)$ and $sin(9 pi x)$ are distinct continuous functions (solid and dashed curves), but they coincide at every grid point (black dots). The wavenumber $9 pi$ aliases to $pi$ because $9 pi = pi + 2 pi dot 4 = pi + 2 pi \/ h$.],
) <fig-aliasing>

The aliasing relation is straightforward to verify. The Nyquist interval for $h = 1\/4$ is $[-4 pi, 4 pi]$. The wavenumber $9 pi$ lies outside this interval. To find its alias, we compute:
$ 9 pi = pi + 8 pi = pi + 2 pi dot frac(1, h) = pi + frac(2 pi, 1\/4). $
Thus $9 pi$ aliases to $pi$: both waves oscillate identically at the grid points.

The following Python code demonstrates this aliasing:

```python
import numpy as np

h = 0.25
x_grid = np.arange(-2, 2 + h, h)

# Both functions give identical samples!
f_samples = np.sin(np.pi * x_grid)
g_samples = np.sin(9 * np.pi * x_grid)

print(f"Max difference: {np.max(np.abs(f_samples - g_samples)):.2e}")
# Output: Max difference: 1.11e-16 (machine precision)
```

The equivalent MATLAB implementation:

```matlab
h = 0.25;
x_grid = -2:h:2;

% Both functions give identical samples
f_samples = sin(pi * x_grid);
g_samples = sin(9 * pi * x_grid);

fprintf('Max difference: %.2e\n', max(abs(f_samples - g_samples)));
% Output: Max difference: 1.11e-16
```

The Julia implementation:

```julia
h = 0.25
x_grid = collect(-2:h:2)

# Both functions give identical samples!
f_samples = sin.(π .* x_grid)
g_samples = sin.(9π .* x_grid)

println("Max difference: $(maximum(abs.(f_samples .- g_samples)))")
# Output: Max difference: 1.1102230246251565e-16
```

*General aliasing formula*: Given a wavenumber $kappa in RR$, its alias $k in [-pi\/h, pi\/h]$ is determined by
$ kappa = k + frac(2 pi m, h) $
for some integer $m$. Explicitly, $k = kappa - frac(2 pi, h) op("round")(frac(kappa h, 2 pi))$ @Trefethen2000.

The code generating @fig-aliasing is available in:
- `codes/python/ch09/aliasing_demo.py`
- `codes/matlab/ch09/aliasing_demo.m`
- `codes/julia/ch09/aliasing_demo.jl`

== Band-Limited Interpolation and the Sinc Kernel <sec-sinc>

=== From Discrete Samples to Continuous Functions

Given samples $v_j$ on the grid $h ZZ$, how do we construct a continuous function that interpolates these values? The semidiscrete Fourier transform provides a canonical answer: construct the unique _band-limited interpolant_ whose Fourier transform is supported in the Nyquist interval $[-pi\/h, pi\/h]$. The mathematical foundations of this construction trace back to Whittaker @Whittaker1915 and have a rich history documented by Jerri @Jerri1977 and Meijering @Meijering2002.

=== The Discrete Delta Function

The simplest grid function is the _Kronecker delta_:
$ delta_j = cases(1 quad& "if" j = 0, 0 quad& "if" j eq.not 0.) $ <eq-kronecker-delta>
Its semidiscrete Fourier transform is constant:
$ hat(delta)(k) = h sum_(j = -infinity)^(infinity) e^(-i k x_j) delta_j = h dot e^0 = h, quad k in [-pi\/h, pi\/h]. $

=== The Sinc Function

Since $hat(delta)(k) = h$ on $[-pi\/h, pi\/h]$ and vanishes outside this interval, applying the inverse transform @eq-semidiscrete-inverse yields the band-limited interpolant of $delta$:
$ p(x) = frac(h, 2 pi) integral_(-pi\/h)^(pi\/h) e^(i k x) dif k = frac(sin(pi x \/ h), pi x \/ h). $
This famous function is called the _sinc function_:
$ S_h (x) = frac(sin(pi x \/ h), pi x \/ h). $ <eq-sinc>
By convention, $S_h (0) = 1$ (the limit as $x arrow 0$).

Sir Edmund Whittaker @Whittaker1915 called $S_1$ "a function of royal blood in the family of entire functions, whose distinguished properties separate it from its bourgeois brethren" @McNameeStenger1971 @Trefethen2000.

=== Sinc Interpolation

The sinc function is the fundamental interpolation kernel for band-limited functions. Any grid function $v$ can be written as
$ v_j = sum_(m = -infinity)^(infinity) v_m delta_(j - m), $
where $delta_(j-m)$ equals $1$ if $j = m$ and $0$ otherwise. Since band-limited interpolation is linear and translation-invariant, the interpolant of $v$ is:
$ p(x) = sum_(m = -infinity)^(infinity) v_m S_h (x - x_m). $ <eq-sinc-interpolation>

This is the _Whittaker--Shannon interpolation formula_, the foundation of the sampling theorem.

#block(
  fill: rgb("#142D6E").lighten(92%),
  stroke: (left: 3pt + rgb("#142D6E")),
  inset: (left: 12pt, y: 10pt, right: 10pt),
  width: 100%,
)[
*Sampling Theorem* (Whittaker--Shannon--Nyquist): A band-limited function $u$ with $hat(u)(k) = 0$ for $|k| > pi\/h$ is completely determined by its samples ${u(x_j)}$ on $h ZZ$. The reconstruction is given by @eq-sinc-interpolation.
]

For a comprehensive tutorial on the history and various extensions of this theorem, see the survey by Jerri @Jerri1977, which traces the result from its origins with Whittaker and Nyquist to Shannon's @Shannon1949 information-theoretic formulation.

=== Computational Étude 2: Sinc Interpolation of Three Signals <sec-etude-sinc>

@fig-sinc-interpolation shows the band-limited interpolants of three grid functions: a discrete delta, a discrete square wave, and a discrete triangular (hat) function.

#figure(
  image("../figures/ch09/python/sinc_interpolation.pdf", width: 85%),
  caption: [Band-limited interpolation of three grid functions ($h = 1$). _Top_: The Kronecker delta $delta_j$ and its interpolant, the sinc function $S_h (x)$. _Middle_: A discrete square wave and its interpolant, exhibiting the Gibbs phenomenon near discontinuities. _Bottom_: A discrete triangular (hat) function and its interpolant, showing moderate oscillations.],
) <fig-sinc-interpolation>

Several observations:
- *Top panel*: The sinc function is smooth and analytic, decaying as $|x|^(-1)$.
- *Middle panel*: The square wave interpolant exhibits _Gibbs phenomenon_: oscillations near discontinuities that do not diminish as $h arrow 0$ @Gegenbauer2025 @Trefethen2000. Band-limited interpolation cannot approximate discontinuous functions well.
- *Bottom panel*: The hat function is continuous but not differentiable at its corners. Its interpolant is smoother than the square wave case but still oscillatory.

The core of the sinc interpolation algorithm is simple:

```python
def sinc_kernel(z):
    """Compute sinc(z) = sin(pi*z)/(pi*z), with sinc(0) = 1."""
    out = np.ones_like(z, dtype=float)
    mask = z != 0
    out[mask] = np.sin(np.pi * z[mask]) / (np.pi * z[mask])
    return out

def sinc_interpolate(x_grid, v, x_fine, h=1.0):
    """Band-limited interpolation via sinc functions."""
    p = np.zeros_like(x_fine)
    for xi, vi in zip(x_grid, v):
        p += vi * sinc_kernel((x_fine - xi) / h)
    return p
```

The equivalent MATLAB implementation:

```matlab
function p = sinc_interpolate(x_grid, v, x_fine, h)
    % Band-limited interpolation via sinc functions
    p = zeros(size(x_fine));
    for i = 1:length(x_grid)
        z = (x_fine - x_grid(i)) / h;
        s = ones(size(z));
        mask = z ~= 0;
        s(mask) = sin(pi*z(mask)) ./ (pi*z(mask));
        p = p + v(i) * s;
    end
end
```

The Julia implementation:

```julia
function sinc_kernel(z)
    out = ones(length(z))
    for i in eachindex(z)
        z[i] != 0 && (out[i] = sin(π * z[i]) / (π * z[i]))
    end
    return out
end

function sinc_interpolate(x_grid, v, x_fine; h=1.0)
    p = zeros(length(x_fine))
    for (xi, vi) in zip(x_grid, v)
        p .+= vi .* sinc_kernel((x_fine .- xi) ./ h)
    end
    return p
end
```

The code generating @fig-sinc-interpolation is available in:
- `codes/python/ch09/sinc_interpolation.py`
- `codes/matlab/ch09/sinc_interpolation.m`
- `codes/julia/ch09/sinc_interpolation.jl`

== Periodic Grids: The Discrete Fourier Transform and FFT <sec-dft>

=== The Periodic Setting

We now turn to the practical setting: a bounded, periodic domain. Consider the interval $[0, 2 pi)$ with periodic boundary conditions. We place $N$ equispaced grid points:
$ x_j = frac(2 pi j, N), quad j = 0, 1, dots, N - 1. $ <eq-periodic-grid>
The grid spacing is $h = 2 pi \/ N$, which gives the fundamental relation:
$ frac(pi, h) = frac(N, 2). $ <eq-pi-over-h>
This equation appears repeatedly in periodic spectral methods.

A _periodic grid function_ $v$ satisfies $v_(j + N) = v_j$ for all $j$. Equivalently, we can view it as one period of length $N$ extracted from an infinite periodic sequence.

=== Fourier Space for Periodic Functions

Two constraints shape the Fourier space for periodic grids:

1. *Periodicity in physical space*: The function must have period $2 pi$, so only waves $e^(i k x)$ with _integer_ wavenumbers $k$ have the correct period.

2. *Discreteness in physical space*: Aliasing restricts distinguishable wavenumbers to an interval of length $2 pi \/ h = N$.

Combining these constraints:
#align(center)[
  #box(stroke: 0.5pt + rgb("#142D6E"), inset: 1em, radius: 3pt)[
    #grid(
      columns: (auto, auto, auto, auto, auto),
      align: center + horizon,
      gutter: 0.5em,
      [Physical space:], [discrete,], [bounded], [$arrow.r.double$], [$x_j in {0, h, 2h, dots, (N-1)h}$],
      [Fourier space:], [bounded,], [discrete], [$arrow.r.double$], [$k in {-N\/2 + 1, dots, N\/2}$],
    )
  ]
]

=== The Discrete Fourier Transform

For $N$ even, the _discrete Fourier transform (DFT)_ is
$ hat(v)_k = h sum_(j = 0)^(N - 1) e^(-i k x_j) v_j, quad k = -N\/2 + 1, dots, N\/2. $ <eq-dft>
The _inverse DFT_ is
$ v_j = frac(1, 2 pi) sum_(k = -N\/2 + 1)^(N\/2) e^(i k x_j) hat(v)_k, quad j = 0, 1, dots, N - 1. $ <eq-inverse-dft>

These are exact inverses: no approximation is involved. The DFT maps $N$ function values to $N$ Fourier coefficients and back.

=== The Nyquist Mode and Symmetry

The wavenumber $k = N\/2$ requires special treatment. The wave $e^(i (N\/2) x)$ equals $cos(N x \/ 2)$ on the grid (a "sawtooth" pattern alternating $plus.minus 1$). For differentiation, we symmetrize the highest mode by setting $hat(v)_(-N\/2) = hat(v)_(N\/2)$, giving:
$ v_j = frac(1, 2 pi) sum_(k = -N\/2)^(N\/2) {}' e^(i k x_j) hat(v)_k, $
where the prime indicates that the $k = plus.minus N\/2$ terms are multiplied by $1\/2$.

=== The Fast Fourier Transform

The DFT can be computed naively in $O(N^2)$ operations, but the _Fast Fourier Transform (FFT)_ reduces this to $O(N log N)$. Discovered by Cooley and Tukey @Cooley1965 in 1965 (though Gauss had the idea in 1805, as documented by Heideman, Johnson, and Burrus @HeidmanJohnson1984), the FFT revolutionized scientific computing.

The key insight is that a DFT of size $N$ can be decomposed into two DFTs of size $N\/2$: one for even-indexed entries, one for odd. This divide-and-conquer approach yields the $O(N log N)$ complexity.

=== FFT Indexing Conventions

MATLAB and Python store FFT output with different conventions than our mathematical formulas. @tbl-fft-indexing clarifies the correspondence.

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    table(
      columns: (1fr, 1fr, 1fr, 1.5fr),
      align: (center, center, center, left),
      inset: (x: 0.8em, y: 0.6em),
      stroke: none,
      table.hline(stroke: 0.75pt + rgb("#142D6E")),
      table.header(
        table.cell(fill: rgb("#142D6E").lighten(85%))[*MATLAB (1-based)*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Python (0-based)*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Math mode $k$*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Description*],
      ),
      table.hline(stroke: 0.5pt + luma(180)),
      [1], [0], [0], [Constant (mean)],
      [2 to $N\/2$], [1 to $N\/2 - 1$], [1 to $N\/2 - 1$], [Positive modes],
      [$N\/2 + 1$], [$N\/2$], [$N\/2$], [Nyquist mode],
      [$N\/2 + 2$ to $N$], [$N\/2 + 1$ to $N - 1$], [$-N\/2 + 1$ to $-1$], [Negative modes],
    ),
  ),
  caption: [FFT indexing conventions for $N$ even. The wavenumber vector in Python is `np.fft.fftfreq(N)*N`, and in MATLAB it is `[0:N/2, -N/2+1:-1]`.],
) <tbl-fft-indexing>

=== Spectral Differentiation via FFT

The deepest reason to move between physical and Fourier space is that _differentiation becomes multiplication_. In physical space, computing a derivative requires approximating a limit --- finite differences, polynomial interpolation, or some other local procedure that inevitably introduces truncation error. In Fourier space, the derivative is _exact_: each mode $e^(i k x)$ satisfies $dif / (dif x) e^(i k x) = i k e^(i k x)$, so differentiating a trigonometric polynomial amounts to multiplying each coefficient by $i k$. No approximation is needed, no stencil must be designed, and the only error comes from the truncation of the Fourier series itself, which, for smooth functions, is exponentially small.

This observation yields a remarkably simple algorithm. Given $N$ samples $v_j = v(x_j)$ of a periodic function on the grid $x_j = 2 pi j \/ N$, the spectral derivative is computed in three steps @Trefethen2000 @Frigo2005:
1. *Transform*: Compute $hat(v) = "FFT"(v)$, moving from physical to Fourier space.
2. *Multiply*: Set $hat(w)_k = i k hat(v)_k$ for each wavenumber $k$. This is the Fourier-space equivalent of differentiation.
3. *Invert*: Compute $w = "IFFT"(hat(w))$, returning to physical space.

One subtlety deserves attention: the Nyquist mode $k = N\/2$. This highest-frequency mode corresponds to the pattern $(-1)^j$ on the grid, which has no well-defined derivative direction --- it oscillates maximally but could be interpreted as either a positive or a negative frequency. Following the symmetrisation convention discussed above, we set $hat(w)_(N\/2) = 0$, effectively excluding this ambiguous mode from the derivative. For well-resolved functions (those whose Fourier coefficients have decayed to negligible values well before the Nyquist frequency), this choice has no practical effect.

The total cost is dominated by the two FFTs: $O(N log N)$ operations, compared to $O(N^2)$ for a dense differentiation matrix or $O(N)$ for a finite difference stencil. But unlike finite differences, the spectral derivative achieves _exponential_ accuracy for analytic functions --- the best of both worlds.

The following Python code implements this:

```python
def spectral_derivative(v):
    """Compute spectral derivative of periodic function."""
    N = len(v)
    v_hat = np.fft.fft(v)
    k = np.fft.fftfreq(N) * N  # Wavenumbers: 0,1,...,N/2,-N/2+1,...,-1
    k = 1j * k
    k[N//2] = 0  # Zero out Nyquist mode for odd derivatives
    w_hat = k * v_hat
    return np.real(np.fft.ifft(w_hat))
```

The equivalent MATLAB implementation:

```matlab
function w = spectral_derivative(v)
    N = length(v);
    v_hat = fft(v);
    k = [0:N/2-1, 0, -N/2+1:-1];  % Wavenumbers with zeroed Nyquist
    w_hat = 1i * k .* v_hat;
    w = real(ifft(w_hat));
end
```

The Julia implementation:

```julia
using FFTW

function spectral_derivative(v)
    N = length(v)
    v_hat = fft(v)
    k = [0:N÷2-1; 0; -N÷2+1:-1]  # Wavenumbers with zeroed Nyquist
    w_hat = 1im .* k .* v_hat
    return real.(ifft(w_hat))
end
```

=== Computational Étude 3: Spectral Differentiation Accuracy <sec-etude-spectral-diff>

To verify the spectral accuracy of FFT-based differentiation, we test it on the analytic periodic function $f(x) = exp(sin x)$, whose derivative is $f'(x) = cos(x) exp(sin x)$. For each grid size $N$, we sample $f$ at $x_j = 2 pi j \/ N$ for $j = 0, dots, N - 1$, compute the spectral derivative, and measure the maximum absolute error against the exact derivative.

@tab-spectral-diff-errors shows the results. The error decreases _exponentially_ with $N$, reaching machine precision ($approx 10^(-15)$) by $N = 28$. This dramatic convergence --- a hallmark of spectral methods for analytic functions --- should be compared with the algebraic convergence of finite difference methods, which would require hundreds or thousands of points to achieve comparable accuracy.

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    {
      show table: format-table(auto, auto, auto, auto)
      table(
        columns: 4,
        align: (center, center, center, center),
        inset: (x: 1em, y: 0.6em),
        stroke: none,
        table.hline(stroke: 0.75pt + rgb("#142D6E")),
        table.header(
          table.cell(fill: rgb("#142D6E").lighten(85%))[*$N$*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*Max error*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*$N$*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*Max error*],
        ),
        table.hline(stroke: 0.5pt + luma(180)),
        [4], num[1.75e-1], [20], num[4.99e-10],
        [8], num[4.32e-3], [24], num[9.55e-13],
        [12], num[3.82e-5], [28], num[5.55e-15],
        [16], num[1.76e-7], [32], num[2.89e-15],
      )
    },
  ),
  caption: [Maximum error in the spectral derivative of $f(x) = exp(sin x)$ for various grid sizes $N$. The exponential convergence to machine precision is characteristic of spectral methods applied to analytic functions.],
) <tab-spectral-diff-errors>

The following Python code computes this error table:

```python
f = lambda x: np.exp(np.sin(x))
f_prime = lambda x: np.cos(x) * np.exp(np.sin(x))

for N in [4, 8, 12, 16, 20, 24, 28, 32]:
    x = 2 * np.pi * np.arange(N) / N
    w = spectral_derivative(f(x))
    error = np.max(np.abs(w - f_prime(x)))
    print(f"N = {N:2d}, error = {error:.2e}")
```

The equivalent MATLAB implementation:

```matlab
for N = [4, 8, 12, 16, 20, 24, 28, 32]
    x = 2*pi*(0:N-1)/N;
    v = exp(sin(x));
    w = spectral_derivative(v);
    w_exact = cos(x) .* exp(sin(x));
    fprintf('N = %2d, error = %.2e\n', N, max(abs(w - w_exact)));
end
```

The Julia implementation:

```julia
for N in [4, 8, 12, 16, 20, 24, 28, 32]
    x = 2π .* collect(0:N-1) ./ N
    w = spectral_derivative(exp.(sin.(x)))
    w_exact = cos.(x) .* exp.(sin.(x))
    @printf("N = %2d, error = %.2e\n", N, maximum(abs.(w .- w_exact)))
end
```

The code generating @tab-spectral-diff-errors is available in:
- `codes/python/ch09/spectral_derivative_accuracy.py`
- `codes/matlab/ch09/spectral_derivative_accuracy.m`
- `codes/julia/ch09/spectral_derivative_accuracy.jl`

== Aliasing and Spectra on Periodic Grids <sec-periodic-aliasing>

=== Aliasing in FFT Computations

When we compute the FFT of a sampled function, any frequency content above the Nyquist frequency $N\/2$ folds back into the resolved range. This is the periodic analog of the aliasing we saw in @sec-semidiscrete. The foundational work on mitigating this error in nonlinear computations is due to Orszag @OrszagDealiasing1971, who introduced the famous "2/3 rule" for dealiasing; see also Bowman and Roberts @Bowman2011 for efficient modern implementations.

=== Computational Étude 4: Aliasing in FFT of High-Frequency Data <sec-etude-fft-aliasing>

Consider $u(x) = sin(17 x)$ sampled at $N = 32$ points. Since $17 > N\/2 = 16$, this frequency is above the Nyquist limit and will alias.

#figure(
  image("../figures/ch09/python/fft_aliasing.pdf", width: 85%),
  caption: [Aliasing in the FFT of $sin(17 x)$ with $N = 32$ points. The true frequency $k = 17$ lies above the Nyquist limit $N\/2 = 16$. The FFT shows a spike at $k = -15$ because $17 = -15 + 32$: the frequency "wraps around" the Nyquist boundary.],
) <fig-fft-aliasing>

The alias is computed as: $17 = -15 + 32$, so wavenumber 17 appears as wavenumber $-15$ in the FFT output.

```python
N = 32
x = 2 * np.pi * np.arange(N) / N
u = np.sin(17 * x)

u_hat = np.fft.fft(u) / N
k = np.fft.fftshift(np.fft.fftfreq(N) * N)

# The spike appears at k = -15, not k = 17
```

The equivalent MATLAB implementation:

```matlab
N = 32;
x = 2*pi*(0:N-1)/N;
u = sin(17*x);

u_hat = fft(u) / N;
k = fftshift([0:N/2-1, -N/2:-1]);

% The spike appears at k = -15, not k = 17
```

The Julia implementation:

```julia
N = 32
x = 2π * (0:N-1) / N
u = sin.(17 .* x)

u_hat = fft(u) / N
k = fftshift(fftfreq(N) .* N)

# The spike appears at k = -15, not k = 17
```

The code generating @fig-fft-aliasing is available in:
- `codes/python/ch09/fft_aliasing.py`
- `codes/matlab/ch09/fft_aliasing.m`
- `codes/julia/ch09/fft_aliasing.jl`

The concept of aliasing has found renewed significance beyond numerical analysis. In the field of deep learning, aliasing in neural radiance fields (NeRF) has been identified as a critical bottleneck; see Mildenhall _et al._ @Mildenhall2021 for a comprehensive treatment. The classical spectral perspective on aliasing developed by Trefethen @Trefethen2000 thus remains fundamental across a wide range of modern applications.

=== What Your Eye Aliases

Aliasing occurs not just in computation but in perception. Consider plotting $sin(n)$ for $n = 1, 2, dots, 3000$. The function $sin(n)$ oscillates rapidly --- its argument increases by 1 radian per step, which is roughly $1 \/ (2 pi) approx 16%$ of a full cycle. Yet when we scatter 3000 such points, the eye perceives a slowly undulating envelope, as shown in @fig-eye-aliasing. The explanation is pure aliasing: the true frequency ($1$ radian per sample) lies far above the "Nyquist limit" of our visual system, and folds back to a low frequency that our perception resolves as a smooth modulation.

#figure(
  image("../figures/ch09/python/eye_aliasing.pdf", width: 85%),
  caption: [Visual aliasing: $sin(n)$ for integer $n = 1, dots, 3000$. Although the function oscillates rapidly (nearly one radian per step), the plotted points appear to trace a slowly varying envelope. The eye, acting as a sampler, aliases the high frequency to a low-frequency pattern.],
) <fig-eye-aliasing>

The code generating @fig-eye-aliasing is available in:
- `codes/python/ch09/eye_aliasing.py`
- `codes/matlab/ch09/eye_aliasing.m`
- `codes/julia/ch09/eye_aliasing.jl`

== Spectra and Smoothness <sec-spectra-smoothness>

=== The Smoothness-Decay Connection

One of the most important results in Fourier analysis is the connection between smoothness in physical space and decay in Fourier space:

- *Discontinuous functions*: Fourier coefficients decay as $O(|k|^(-1))$
- *$p$ continuous derivatives*: Decay as $O(|k|^(-p-1))$
- *Infinitely smooth ($C^infinity$)*: Faster than any polynomial
- *Analytic functions*: Exponential decay $O(e^(-a |k|))$ @Boyd2000 @Trefethen2013

This hierarchy is a classical result in harmonic analysis; see, e.g., Trefethen @Trefethen2000, Fourier @Fourier1822, and the general theory of $(p,q)$-Fourier coefficient decay by Edmunds, Gurka, and Lang @FourierDecay2014. It explains why spectral methods are so accurate for smooth problems: smooth functions have negligible high-frequency content, so the aliasing error from discretization is tiny.

=== Computational Demonstration

@fig-smoothness-spectra shows the spectra of three periodic functions with different smoothness:

1. $exp(sin x)$: Analytic (infinitely differentiable and more), with exponential spectral decay
2. Periodic hat function: $C^0$ but not $C^1$ at corners, with algebraic decay $approx |k|^(-2)$
3. Square wave: Jump discontinuities, with slow decay $approx |k|^(-1)$

#figure(
  image("../figures/ch09/python/smoothness_spectra.pdf", width: 95%),
  caption: [Fourier spectra of three functions with different smoothness. _Left_: $exp(sin x)$ (analytic) shows exponential decay on a semilog plot. _Center_: Periodic hat function ($C^0$) shows algebraic decay $approx |k|^(-2)$ on a log-log plot. _Right_: Square wave (discontinuous) shows slow decay $approx |k|^(-1)$.],
) <fig-smoothness-spectra>

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    table(
      columns: (1.5fr, 1.5fr, 1.5fr),
      align: (left, center, center),
      inset: (x: 1em, y: 0.6em),
      stroke: none,
      table.hline(stroke: 0.75pt + rgb("#142D6E")),
      table.header(
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Function*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Regularity*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Decay of $|hat(f)_k|$*],
      ),
      table.hline(stroke: 0.5pt + luma(180)),
      [$exp(sin x)$], [Analytic], [Faster than any power],
      [Hat function], [$C^0$, piecewise $C^1$], [$approx |k|^(-2)$],
      [Square wave], [Jump discontinuities], [$approx |k|^(-1)$],
    ),
  ),
  caption: [Smoothness classes and their characteristic Fourier decay rates.],
) <tbl-smoothness-decay>

The code generating @fig-smoothness-spectra is available in:
- `codes/python/ch09/smoothness_spectra.py`
- `codes/matlab/ch09/smoothness_spectra.m`
- `codes/julia/ch09/smoothness_spectra.jl`

== Periodic Band-Limited Interpolation <sec-periodic-sinc>

=== The Periodic Sinc Function

Just as the sinc function is the band-limited interpolant of the Kronecker delta on $h ZZ$, the _periodic sinc_ is the band-limited interpolant of the periodic delta on the $N$-point periodic grid.

The periodic delta is defined on the grid $x_j = j h$ (where $h = 2 pi \/ N$) as:
$ delta_j = cases(1 quad& "if" j equiv 0 space (mod N), 0 quad& "otherwise.") $

To derive the continuous interpolant, we first move to Fourier space. The Discrete Fourier Transform (DFT) of the periodic delta is computed as:
$ hat(delta)_k = h sum_(j=1)^N e^(-i k x_j) delta_j. $
Since $delta_j$ is non-zero only when $j equiv 0 (mod N)$ (which corresponds to index $N$ in the sum), and $e^(-i k x_N) = e^(-i k 2 pi) = 1$, the sum collapses to a single term:
$ hat(delta)_k = h dot 1 = h. $
Thus, the Fourier coefficients are constant, $hat(delta)_k = h$, for all wavenumbers $k$.

The periodic sinc function $S_N (x)$ is obtained via the inverse DFT. However, because we are interpolating on an even number of points $N$, we must handle the highest wavenumber $k = N\/2$ carefully to preserve symmetry and ensure the result is real-valued. We use the symmetric inverse transform, where the terms at $k = plus.minus N\/2$ are weighted by $1\/2$:

$ S_N (x) &= 1 / (2 pi) sum_(k=-N\/2)^(N\/2) {}' hat(delta)_k e^(i k x) \
          &= h / (2 pi) sum_(k=-N\/2)^(N\/2) {}' e^(i k x) \
          &= 1 / N sum_(k=-N\/2)^(N\/2) {}' e^(i k x), $
where the prime indicates that the first and last terms of the summation are halved. Expanding the summation limits, we write:
$ S_N (x) = 1 / N [ 1 / 2 e^(-i N x \/ 2) + sum_(k=-N\/2 + 1)^(N\/2 - 1) e^(i k x) + 1 / 2 e^(i N x \/ 2) ]. $
Grouping the boundary terms and the internal geometric series:
$ S_N (x) = 1 / N [ cos(N x \/ 2) + sum_(k=- (N\/2 - 1))^(N\/2 - 1) e^(i k x) ]. $
The geometric sum in the brackets is the Dirichlet kernel of order $M = N\/2 - 1$:
$ sum_(k=-M)^(M) e^(i k x) = frac(sin((M + 1\/2)x), sin(x\/2)) = frac(sin((N-1)x \/ 2), sin(x\/2)). $
Substituting this back into the expression for $S_N (x)$:
$ S_N (x) = 1 / N [ cos(N x \/ 2) + frac(sin(N x \/ 2 - x \/ 2), sin(x\/2)) ]. $
Using the subtraction formula $sin(A - B) = sin A cos B - cos A sin B$, the numerator becomes:
$ sin(N x \/ 2) cos(x\/2) - cos(N x \/ 2) sin(x\/2). $
When divided by $sin(x\/2)$, the second term cancels the isolated cosine term in our expression for $S_N (x)$, leaving:
$ S_N (x) &= 1 / N [ cos(N x \/ 2) + frac(sin(N x \/ 2) cos(x\/2), sin(x\/2)) - cos(N x \/ 2) ] \
          &= 1 / N [ frac(sin(N x \/ 2) cos(x\/2), sin(x\/2)) ] \
          &= 1 / N sin(N x \/ 2) cot(x\/2). $
Finally, substituting $N = 2 pi \/ h$, we obtain the standard form. The term $sin(N x \/ 2)$ becomes $sin(pi x \/ h)$, and the factor $1\/N cot(x\/2)$ becomes $h \/ (2 pi) dot 1 \/ tan(x\/2)$.

This yields the explicit formula for the periodic sinc:
$ S_N (x) = frac(sin(pi x \/ h), (2 pi \/ h) tan(x \/ 2)), $ <eq-periodic-sinc>
with $h = 2 pi \/ N$. This derivation is presented in detail by Trefethen @Trefethen2000 and follows naturally from the discrete orthogonality of complex exponentials on the periodic grid @GottliebOrszag1977.

For small $x$, the periodic sinc behaves like the nonperiodic sinc: $S_N (x) approx sin(pi x \/ h) \/ (pi x \/ h)$. The difference appears for larger $|x|$, where the periodic sinc has period $2 pi$.

=== Computational Étude 5: Zero-Padding Interpolation via FFT <sec-etude-zero-padding>

Suppose we have $N$ samples of a periodic function and we want to evaluate the trigonometric interpolant at a finer set of $M = q N$ equally spaced points. A naive approach would build the $N$-term trigonometric polynomial and evaluate it at each of the $M$ points, costing $O(M N)$ operations. _Zero-padding_ in Fourier space achieves the same result in $O(M log M)$ operations, by exploiting a simple observation: if a function has no energy at high frequencies, then adding zero Fourier coefficients at those frequencies does not change it. The algorithm proceeds as follows @Trefethen2000 @Frigo2005:

1. Compute the DFT $hat(v)_k$ from the $N$ samples.
2. Create a larger array of size $M = q N$, placing the $N$ computed coefficients $hat(v)_k$ in the low-frequency positions and zeros in the high-frequency positions.
3. Compute the inverse DFT of the padded array to obtain $M$ interpolated values.

The zeros we insert correspond to wavenumbers that the coarse grid cannot resolve. By setting them to zero, we are stating: "the trigonometric interpolant has no content at those frequencies." The inverse DFT on the finer grid then evaluates this same trigonometric polynomial at $M$ points instead of $N$ --- an exact operation, not an approximation.

@fig-zero-padding illustrates this for $f(x) = exp(sin x)$ with $N = 32$ coarse samples and $M = 128$ interpolated points. The top panel shows that the interpolant (solid line) passes exactly through the original samples (circles) and lies on top of the true function (dashed line). The bottom panel confirms that the pointwise error is at machine precision ($approx 10^(-14)$), since $exp(sin x)$ is analytic and $N = 32$ is sufficient to resolve all its significant Fourier modes.

#figure(
  image("../figures/ch09/python/zero_padding_interpolation.pdf", width: 85%),
  caption: [Band-limited interpolation via zero-padding. _Top_: the true function $exp(sin x)$ (dashed), the interpolant from zero-padding to $M = 128$ points (solid), and the $N = 32$ coarse samples (circles). The curves are visually indistinguishable. _Bottom_: pointwise error $|p(x) - f(x)|$ on a logarithmic scale, confirming machine-precision accuracy for this analytic function.],
) <fig-zero-padding>

```python
def zero_pad_interpolate(v, q=4):
    """Interpolate by zero-padding in Fourier space."""
    N = len(v)
    M = q * N

    v_hat = np.fft.fft(v)
    v_hat_padded = np.zeros(M, dtype=complex)

    # Copy low frequencies (0 to N/2-1 and -N/2+1 to -1)
    v_hat_padded[:N//2] = v_hat[:N//2]
    v_hat_padded[-N//2+1:] = v_hat[-N//2+1:]
    # Handle Nyquist: split between +N/2 and -N/2
    v_hat_padded[N//2] = v_hat[N//2] / 2
    v_hat_padded[-N//2] = v_hat[N//2] / 2

    return np.real(np.fft.ifft(v_hat_padded)) * q
```

The equivalent MATLAB implementation:

```matlab
function v_fine = zero_pad_interpolate(v, q)
    N = length(v);
    M = q * N;

    v_hat = fft(v);
    v_hat_padded = zeros(1, M);

    % Copy low frequencies
    v_hat_padded(1:N/2) = v_hat(1:N/2);
    v_hat_padded(end-N/2+2:end) = v_hat(N/2+2:end);
    % Handle Nyquist mode
    v_hat_padded(N/2+1) = v_hat(N/2+1) / 2;
    v_hat_padded(M-N/2+1) = v_hat(N/2+1) / 2;

    v_fine = real(ifft(v_hat_padded)) * q;
end
```

The Julia implementation:

```julia
using FFTW

function zero_pad_interpolate(v; q=4)
    N = length(v)
    M = q * N
    v_hat = fft(v)
    v_hat_padded = zeros(ComplexF64, M)

    # Copy low frequencies
    v_hat_padded[1:N÷2] = v_hat[1:N÷2]
    v_hat_padded[end-N÷2+2:end] = v_hat[end-N÷2+2:end]
    # Handle Nyquist mode: split between +N/2 and -N/2
    v_hat_padded[N÷2+1]   = v_hat[N÷2+1] / 2
    v_hat_padded[M-N÷2+1] = v_hat[N÷2+1] / 2

    return real.(ifft(v_hat_padded)) .* q
end
```

The code generating @fig-zero-padding is available in:
- `codes/python/ch09/zero_padding_interpolation.py`
- `codes/matlab/ch09/zero_padding_interpolation.m`
- `codes/julia/ch09/zero_padding_interpolation.jl`

== A non-exhaustive literature overview

The mathematical machinery codified in this chapter --- the Fourier transform, sampling, aliasing, and the FFT --- is the product of three centuries of analytical progress. While the name "Fast Fourier Transform" is inextricably linked to the 1965 publication by Cooley and Tukey @Cooley1965, the intellectual history runs far deeper. Historical analysis of Carl Friedrich Gauss's _Nachlass_ (unpublished notes) reveals that in 1805, while seeking to interpolate the orbits of asteroids Pallas and Juno from limited observational data, Gauss derived an algorithm functionally equivalent to the FFT @HeidmanJohnson1984. Gauss's method, designed for hand calculation, recognized that trigonometric sums could be split based on the prime factors of the number of data points, effectively reducing the arithmetic load from $O(N^2)$ to $O(N log N)$ long before the "Big O" notation existed. It is a striking historical irony that this method remained largely dormant, overshadowed by Fourier's analytical work on heat conduction published in 1822 @Fourier1822. Fourier's insistence that arbitrary functions --- even those with discontinuities --- could be represented by trigonometric series was revolutionary and controversial, sparking decades of debate regarding convergence that eventually birthed the field of harmonic analysis. The _discrete_ algorithmic efficiency discovered by Gauss was lost to the broader scientific community until the digital age necessitated it. The 1965 Cooley--Tukey paper was the catalyst for the modern spectral revolution, not because it was the first discovery, but because it appeared at the precise moment when digital computers became capable of executing these transforms at scale.

The transition from continuous signals to discrete samples, treated in @sec-semidiscrete and @sec-sinc, rests on the Whittaker--Shannon--Kotelnikov (WKS) sampling theorem. While often attributed simply to Shannon (1949) or Nyquist (1928) in engineering contexts, the mathematical rigorousness of the theorem traces back to the interpolation theory of the early 20th century. E. T. Whittaker's seminal 1915 paper @Whittaker1915 is of particular importance to the narrative of this textbook. Whittaker investigated the properties of the "cardinal function" --- what we now call the sinc function expansion --- and established that it provides the unique consistent interpolation for band-limited functions. Whittaker famously described this function as "a function of royal blood in the family of entire functions, whose distinguished properties separate it from its bourgeois brethren" @McNameeStenger1971. This evocative phrasing highlights the exceptional nature of the sinc kernel: it is an entire function (analytic everywhere in the complex plane) that perfectly reconstructs band-limited data, satisfying the "spectral promise" of exponential accuracy. The genealogy extends even further back. Historical scholarship indicates that Augustin-Louis Cauchy was aware of the mechanisms of band-limited sampling as early as 1841, and Émile Borel stated the essential features of the reconstruction theorem in 1897 @Meijering2002. Shannon's contribution in 1949 @Shannon1949 was to transplant these mathematical identities into the fertile soil of information theory, proving that the sampling limit (the Nyquist rate) was a fundamental bound on information capacity, not just a mathematical curiosity. For a comprehensive tutorial review, see Jerri @Jerri1977. For the student of spectral methods, this distinction is crucial: the sampling theorem validates the use of discrete grids to solve continuous PDEs, guaranteeing that no information is lost as long as the solution remains sufficiently smooth (band-limited).

Aliasing, introduced in @sec-periodic-aliasing as the phenomenon where high frequencies "masquerade" as low frequencies, is often presented solely as a source of error to be eliminated. The literature, however, presents a more nuanced view, evolving from strict mitigation strategies in the 1970s to sophisticated exploitation in modern machine learning and optics. The foundational work on mitigating this error was established by Orszag @OrszagDealiasing1971, who introduced the famous "2/3 rule" (often called the 3/2 rule in terms of dealiasing padding), which proves that aliasing can be completely eliminated for quadratic nonlinearities by padding the spectrum with zeros to extend the grid by a factor of $3\/2$ before performing the nonlinear multiplication in physical space. This technique remains the gold standard in Direct Numerical Simulation (DNS) of turbulence, where preserving the fidelity of the energy cascade is paramount. Bowman and Roberts @Bowman2011 developed efficient dealiased convolution algorithms that avoid the overhead of explicit zero-padding. In the period 2024--2026, the literature reflects a paradigm shift where aliasing is viewed through new lenses. In the field of deep learning, specifically in Physics-Informed Neural Networks (PINNs) and Neural Operators (like the Fourier Neural Operator, FNO), aliasing has emerged as a critical bottleneck @Mildenhall2021. Standard activation functions introduce infinite spectral content, leading to aliasing errors that degrade the convergence of neural PDE solvers. This has led to the development of "aliasing-free" architectures and spectral activation functions designed to respect the Nyquist limit of the underlying grid @AntiAliasNet2025, effectively reinventing Orszag's principles for the age of AI. Furthermore, in the domain of optics and nanophotonics, recent research describes "anti-aliasing metasurfaces" @Metasurfaces2025. Here, aliasing is not a computational artifact but a physical diffraction phenomenon. By developing multidimensional sampling theories that transcend the traditional Nyquist limit, researchers have created flat optical devices that suppress diffraction noise, enabling ultra-compact high-resolution imaging systems. This illustrates a profound "third-order insight": the mathematical concept of aliasing, once confined to time-series analysis, now governs the design of physical materials and neural architectures.

This chapter primarily focuses on uniform grids, which allow for the use of the standard FFT. However, the "Spectral Promise" is increasingly being applied to problems where data cannot be sampled uniformly --- from MRI scans to astrophysical observations. The Non-Uniform Fast Fourier Transform (NUFFT) has become a critical area of algorithmic research, exploding in importance between 2023 and 2026. The NUFFT generalizes the FFT to off-grid data points, usually by convolving the non-uniform data onto a regular grid using a carefully chosen kernel (like the Kaiser--Bessel or the "exponential of semicircle" kernel). Barnett, Magland, and af Klinteberg @Barnett2019 developed a high-performance parallel NUFFT library (FINUFFT) based on this exponential of semicircle kernel, achieving order-of-magnitude speedups on modern hardware. In radio astronomy, the detection of the faint 21-cm hydrogen signal from the Epoch of Reionization requires processing visibility data from massive interferometer arrays; Cox _et al._ @FFTvis2025 demonstrated that high-performance NUFFT algorithms are essential for simulating these visibilities at the scale required for the Square Kilometer Array (SKA). Parallel to the NUFFT, the Sparse Fast Fourier Transform (sFFT) breaks the $O(N log N)$ barrier by exploiting the observation that many real-world signals have only $K$ non-zero coefficients in the frequency domain, where $K lt.double N$. The sFFT algorithms compute the transform in sub-linear time, often $O(K log N)$ or even $O(K)$ @Hassanieh2012. This literature suggests a divergence in spectral methods: one branch (NUFFT) focuses on geometric flexibility, while the other (sFFT) focuses on algorithmic efficiency for massive, sparse datasets.

Perhaps the most visionary development in the recent literature (2025--2026) is the intersection of Fourier analysis with quantum computing. While the Quantum Fourier Transform (QFT) has been a staple of quantum algorithms (e.g., Shor's algorithm) for decades, recent work has generalized Fourier analysis to higher-order structures to characterize quantum complexity. A series of breakthrough papers in 2025, including work published in _PNAS_, introduces "Quantum Higher-Order Fourier Analysis" (q-HOFA) @QHOFA2025. This theoretical framework generalizes classical higher-order Fourier analysis (used in additive combinatorics) to the quantum realm. The researchers define "quantum uniformity norms" (analogs of Gowers norms) that characterize the "Clifford hierarchy" --- a classification of quantum gates based on their complexity and fault-tolerance properties. This development is profound for the future of spectral methods. It suggests that the tools of Fourier analysis are not limited to classical scalar functions but can be lifted to operate on operators in Hilbert space. While classical spectral methods scale as $O(N log N)$, quantum spectral methods scale polylogarithmically with the dimension of the vector space, potentially offering exponential speedups for specific classes of smooth, high-dimensional PDEs.

Despite these advances, the core theoretical engine of spectral methods remains the connection between smoothness in physical space and decay in Fourier space. This relationship, formalized in @sec-spectra-smoothness, guarantees that for analytic functions, the Fourier coefficients decay exponentially ($O(e^(-a |k|))$), leading to the famed "spectral accuracy." The recent literature reaffirms and refines this principle. Work in _SIAM Journal on Numerical Analysis_ and _Mathematics of Computation_ (2024--2025) has derived sharp pointwise error estimates for functions with algebraic singularities, showing exactly how convergence deteriorates near non-smooth points @FourierExtension2024. Furthermore, the theory of "Fourier extensions" or "Fourier continuation" has matured, providing robust methods to map non-periodic functions onto periodic domains while preserving spectral convergence rates, effectively bypassing the Gibbs phenomenon for non-periodic BVPs @FourierExtLocal2026. The comprehensive treatments by Trefethen @Trefethen2013 and Boyd @Boyd2000 provide the rigorous justification for the computational heuristics used in this textbook. They explain _why_ the method works so well for the test cases in this and subsequent chapters and provide the diagnostic tools to understand failure modes when smoothness is lost.

== Summary <sec-fourier-summary>

This chapter has developed the mathematical framework connecting physical and Fourier space across three settings:

+ *Continuous Fourier transform*: Functions on $RR$ have Fourier transforms on $RR$. Differentiation becomes multiplication by $i k$.

+ *Semidiscrete Fourier transform*: Sampling at spacing $h$ restricts wavenumbers to the Nyquist interval $[-pi\/h, pi\/h]$. Wavenumbers differing by $2 pi \/ h$ are aliases of each other.

+ *Discrete Fourier transform*: On a periodic grid with $N$ points, both physical and Fourier variables are discrete. The FFT computes the DFT in $O(N log N)$ operations.

Key concepts:
- *Aliasing*: High frequencies masquerade as low frequencies on discrete grids
- *Sinc interpolation*: Band-limited reconstruction from samples via the Whittaker--Shannon formula
- *Smoothness implies decay*: Smooth functions have rapidly decaying Fourier coefficients, making aliasing negligible

These concepts form the foundation for spectral methods. In the chapters that follow, we apply them to solve differential equations, exploiting the connection between smoothness and spectral accuracy.

=== Quick Reference

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    table(
      columns: (1fr, 1fr, 1.5fr, 1.5fr),
      align: (left, left, left, left),
      inset: (x: 0.8em, y: 0.6em),
      stroke: none,
      table.hline(stroke: 0.75pt + rgb("#142D6E")),
      table.header(
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Setting*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Variables*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Forward Transform*],
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Inverse Transform*],
      ),
      table.hline(stroke: 0.5pt + luma(180)),
      [Continuous], [$x, k in RR$], [$hat(u) = integral e^(-i k x) u dif x$], [$u = frac(1, 2 pi) integral e^(i k x) hat(u) dif k$],
      [Semidiscrete], [$x_j = j h$, $k in [-frac(pi,h), frac(pi,h)]$], [$hat(v) = h sum_j e^(-i k x_j) v_j$], [$v_j = frac(1, 2 pi) integral e^(i k x_j) hat(v) dif k$],
      [Discrete (DFT)], [$x_j = frac(2 pi j, N)$, $k in ZZ$], [$hat(v)_k = h sum_j e^(-i k x_j) v_j$], [$v_j = frac(1, 2 pi) sum_k e^(i k x_j) hat(v)_k$],
    ),
  ),
  caption: [Quick reference for Fourier transforms in three settings.],
) <tbl-fourier-summary>

=== FFT Idioms

*Python (NumPy):*
```python
import numpy as np
k = np.fft.fftfreq(N) * N      # Wavenumbers: 0,1,...,N/2,-N/2+1,...,-1
u_hat = np.fft.fft(u)          # Forward FFT
u = np.fft.ifft(u_hat)         # Inverse FFT
k_shifted = np.fft.fftshift(k) # Reorder to -N/2,...,N/2-1 for plotting
```

*MATLAB:*
```matlab
k = [0:N/2-1, -N/2:-1];        % Wavenumbers
u_hat = fft(u);                % Forward FFT
u = ifft(u_hat);               % Inverse FFT
k_shifted = fftshift(k);       % Reorder for plotting
```

*Julia (FFTW.jl):*
```julia
using FFTW
k = fftfreq(N) .* N             # Wavenumbers: 0,1,...,N/2,-N/2+1,...,-1
u_hat = fft(u)                  # Forward FFT
u = ifft(u_hat)                 # Inverse FFT
k_shifted = fftshift(k)         # Reorder for plotting
```
