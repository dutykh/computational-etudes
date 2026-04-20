// textbook/chapters/coordinate_transformations.typ
// Chapter 19: Coordinate Transformations and Mapped Spectral Methods
// Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
// Email: denys.dutykh@ku.ac.ae
// Homepage: https://www.denys-dutykh.com/
// Last modified: April 2026

#import "../styles/template.typ": dropcap


= Coordinate Transformations and Mapped Spectral Methods <ch-coord-transforms>

#dropcap[A coordinate transformation is the cheapest trick in the spectral numericist's bag and, very often, the deepest. It costs a line of algebra, delivers in exchange a differently-distributed grid, a new Jacobian weight, a variable-coefficient operator, and sometimes a whole new convergence rate. The pedagogical theme of the chapter is deliberately one question applied over and over to each map we meet: _what defect in the problem is the map curing, and what new difficulty does the map introduce?_ We will build a small reusable toolkit, exercise it on eight computational études, and finish with a decision guide that maps problem pathologies onto mapping strategies. The underlying philosophy follows @Boyd2000 Chapter 16: a map is a _resolution design tool_, to be judged by whether the transformed solution is easier to approximate than the original one, not by whether the formula is elegant.]

By the end of this chapter, you should be able to:

1. Derive first and higher derivatives under a general one-dimensional map $y = f(x)$, and understand how orthogonality weights and quadrature nodes transform.
2. Implement Chebyshev methods either directly in the physical variable $x$ or through the computational variable $t = arccos(x)$, and recognise that these are two arithmetic paths to the same mathematics.
3. Choose and tune maps for semi-infinite and infinite intervals, including the algebraic and logarithmic families.
4. Recognise when exponential boundary clustering (tanh map) can heal a weak endpoint singularity, and equally important, when it cannot.
5. Distinguish practical tensor-product mapping from analytically demanding two-dimensional conformal mapping, and apply the simpler tool first.
6. Deploy the arctan/tan map for sharply localised periodic structures and tune its width parameter $L$ by parameter sweep rather than intuition.
7. Describe the logic of adaptive mappings for moving fronts, and account honestly for their cost.
8. Explain why the almost-equispaced Kosloff--Tal--Ezer grid can improve timestep restrictions but destroy spectral accuracy if the map parameter is chosen badly.

== Prelude: Where Should the Points Go? <sec-ct-prelude>

The opening shock of the chapter must arrive before any theorem. Consider the sharply concentrated, $2 pi$-periodic profile
$ f(y) = exp(-kappa (1 - cos y)), quad y in [-pi, pi], $ <eq-ct-pulse>
with $kappa = 80$. The function is smooth and periodic; the standard theorem on trigonometric interpolation promises spectral convergence, and the standard theorem delivers it. What the theorem does not promise is that the promised rate will arrive at the resolutions a student can actually afford. The width of the pulse is of order $1 \/ sqrt(kappa) approx 0.11$ radians; outside a $plus.minus 0.5$ radian window the function is essentially zero. An equispaced Fourier grid on $[-pi, pi]$ therefore spends the vast majority of its degrees of freedom resolving silence. This is the bad trade that motivates the whole chapter.

An alternative is to _stretch the computational coordinate_ so that it clusters near the pulse. For $2 pi$-periodic data the natural choice is the arctan/tan map, originally due to @Boyd2000 in his period-$pi$ form. We use the period-$2 pi$ variant
$ y = 2 arctan(L tan(x \/ 2)), quad x, y in [-pi, pi], $ <eq-ct-arctantan>
with map parameter $L$. For $L < 1$ the map clusters the grid near $y = 0$; for $L > 1$ it clusters near $y = plus.minus pi$; at $L = 1$ it is the identity. The inverse is explicit: $x = 2 arctan(tan(y \/ 2) \/ L)$.

=== Computational Étude 19.1: A Periodic Pulse on Two Grids <etude-ct-prelude>

We compute the trigonometric interpolant of @eq-ct-pulse on two grids: the uniform Fourier grid $y_j = -pi + 2 pi j \/ N$ and the mapped grid $y_j = 2 arctan(L tan(x_j \/ 2))$ with $L = 0.3$ and the same uniform $x$-spacing. For each $N in {8, 12, ..., 128}$ we evaluate both interpolants on a fine reference grid and record the $L^infinity$ error against the analytic target.

In Python:

```python
import numpy as np
kappa, L = 80.0, 0.3
def target(y):   return np.exp(-kappa * (1 - np.cos(y)))
def phys_grid(N, L):
    x = -np.pi + 2*np.pi*np.arange(N)/N
    return 2*np.arctan(L*np.tan(x/2)), x
def mapped_interp(x_nodes, f_nodes, y_eval, L):
    x_eval = 2*np.arctan(np.tan(y_eval/2)/L)
    c = np.fft.fft(f_nodes) / len(x_nodes)
    k = np.fft.fftfreq(len(x_nodes), d=1.0/len(x_nodes))
    return np.real(np.sum(c[:,None]*np.exp(1j*k[:,None]*(x_eval[None,:]+np.pi)), axis=0))
```

In MATLAB:

```matlab
kappa = 80; L = 0.3;
target = @(y) exp(-kappa*(1 - cos(y)));
N = 96; x = -pi + 2*pi*(0:N-1)/N;
y = 2*atan(L*tan(x/2));
c = fft(target(y)) / N;
k = [0:N/2-1, -N/2:-1];
```

In Julia:

```julia
const KAPPA, L = 80.0, 0.3
target(y) = exp(-KAPPA * (1 - cos(y)))
x = [-pi + 2pi*k/N for k in 0:N-1]
y = 2.0 .* atan.(L .* tan.(x ./ 2))
coeffs = fft(target.(y)) ./ N
```

#figure(
  image("../figures/ch19/python/periodic_pulse_two_grids.pdf", width: 98%),
  caption: [Étude 19.1: a periodic pulse on two grids. Left: the profile $f(y)$ (solid line) and the two grids at $N = 32$ --- uniform Fourier (crosses, offset for visibility) and arctan/tan with $L = 0.3$ (circles). Middle: $L^infinity$ error against $N$; the uniform grid needs $N approx 128$ to reach machine precision, while the mapped grid reaches it at $N approx 64$. Right: Fourier-coefficient magnitudes at $N = 96$; both spectra are geometric, but the mapped spectrum is geometric with a markedly larger decay rate, because in the $x$-coordinate the pulse is broader and smoother.],
) <fig-ct-prelude>

Source files: `codes/python/ch19/periodic_pulse_two_grids.py`, `codes/matlab/ch19/periodic_pulse_two_grids.m`, and `codes/julia/ch19/periodic_pulse_two_grids.jl`.

=== The Three Questions

Every mapped method in this chapter will be subjected to the same three-question interrogation:

#block(stroke: 0.8pt + rgb(20, 45, 110), radius: 3pt, inset: 10pt,
  fill: rgb(248, 250, 254))[
  *The three questions for any map.*
  1. Does the transformed solution $u(f(x))$ look smoother, broader, or less singular than the original $u(y)$?
  2. What variable coefficients, weights, and higher-derivative penalties has the map introduced into the operator?
  3. Do the map's own singularities (zeros of $f'$, poles of $f''$, branch points of $f$) sit safely away from the approximation interval?
]

Étude 19.1 scores the arctan/tan map near-perfectly on question 1 (the pulse is vastly broader in $x$ than in $y$), gently on question 2 (the mapped Fourier derivative is still trigonometric), and trivially on question 3 (the map is entire). By contrast, the very different map of @sec-ct-kte will score well on a different criterion but fail question 3 dramatically, and the consequences will be the chapter's closing warning.

== Chebyshev Polynomials are Cosine Functions in Disguise <sec-ct-cheby-cosine>

The most important map in numerical analysis is hidden in the definition of the Chebyshev polynomial:
$ x = cos(t), quad T_n (x) equiv cos(n t). $ <eq-ct-cheby-cosine>
This identity, which appears to be a curiosity, is the pedagogical entry point of the whole chapter. It says that the Chebyshev machinery developed in @ch-chebyshev is already a _mapped_ Fourier cosine method: a pseudospectral computation on the Chebyshev--Gauss--Lobatto grid is the same as a Fourier cosine pseudospectral computation on a uniform $t$-grid, with the variable change buried in the arithmetic.

=== Two Arithmetic Paths

For a model boundary-value problem
$ u_(x x) (x) - q u(x) = f(x), quad x in [-1, 1], quad u(plus.minus 1) = 0, $ <eq-ct-bvp>
there are two equivalent strategies:

- *The X-path.* Discretise $d^2 \/ d x^2$ on the Chebyshev--Gauss--Lobatto grid using the dense differentiation matrix $bold(D)^2$, remove the boundary rows and columns, solve the $(N-1) times (N-1)$ linear system. Cost: $cal(O)(N^3)$ per solve (or $cal(O)(N^2)$ for diagonal preconditioners).
- *The T-path.* Transform @eq-ct-bvp through $x = cos t$ to a differential equation on $t in [0, pi]$ @Boyd2000:
$ sin(t) u_(t t) - cos(t) u_t - q sin^3 (t) u = sin^3 (t) f(cos t). $ <eq-ct-bvp-t>
Evaluate derivatives via an $cal(O)(N log N)$ FFT on the cosine-sampled function. Cost: $cal(O)(N log N)$ per derivative.

Both paths give the _same_ numerical answer up to rounding, and both converge at the same spectral rate. The choice is a question of arithmetic economy, not of mathematical content. For small to medium $N$ the dense matrix is simpler; for very large $N$ or for many right-hand sides the FFT path wins.

=== Computational Étude 19.2: One Problem, Two Coordinates <etude-ct-cheby-cosine>

We solve @eq-ct-bvp with $q = 4$ and a manufactured source chosen so that $u_("ex") (x) = sin(pi x)$ is an exact solution, once directly in $x$ and once through the cosine FFT. The two codes produce numerically identical spectra of errors.

In Python, the X-path is the short routine we met in @sec-workflow; the T-path wraps an FFT-based Chebyshev differentiation twice:

```python
def chebfft(v):
    N = len(v) - 1
    V = np.concatenate([v, v[N-1:0:-1]])
    U = np.real(np.fft.fft(V))
    k = np.arange(N + 1)
    w_hat = 1j * np.concatenate([k[:N], [0], k[1:N] - N]) * U
    W = np.real(np.fft.ifft(w_hat))
    w = np.zeros(N + 1)
    x = np.cos(np.pi * np.arange(N + 1) / N)
    w[1:N] = -W[1:N] / np.sqrt(1 - x[1:N]**2)
    # endpoint formulas as in Trefethen, Program 21
    return w
```

In MATLAB:

```matlab
function w = chebfft(v)
    N = length(v) - 1;
    V = [v; v(N:-1:2)];
    U = real(fft(V));
    k = (0:N)';
    w_hat = 1i * [k(1:N); 0; k(2:N) - N] .* U;
    W = real(ifft(w_hat));
    % ...endpoint formulas as in Trefethen (2000)
end
```

In Julia:

```julia
function chebfft(v)
    N = length(v) - 1
    V = vcat(v, reverse(v[2:N]))
    U = real.(fft(V))
    k = 0:N
    w_hat = im .* vcat(collect(k[1:N]), 0, collect(k[2:N]) .- N) .* U
    W = real.(ifft(w_hat))
    # ...endpoint formulas
end
```

#figure(
  image("../figures/ch19/python/chebyshev_as_cosine.pdf", width: 78%),
  caption: [Étude 19.2: one problem, two coordinates. The X-path (dense differentiation matrix) and the T-path (FFT-based cosine differentiation) give identical spectral convergence of $u_("ex") (x) = sin(pi x)$ for the boundary-value problem @eq-ct-bvp. The left panel shows the interpolant at $N = 24$ overlaid on the exact solution; the right panel confirms that both paths reach machine precision at $N approx 20$, differing only by a factor close to unity at the floating-point noise floor.],
) <fig-ct-cheby-cosine>

Source files: `codes/python/ch19/chebyshev_as_cosine.py`, `codes/matlab/ch19/chebyshev_as_cosine.m`, and `codes/julia/ch19/chebyshev_as_cosine.jl`.

=== Verdict

The map $x = cos t$ has cured nothing, because there was nothing to cure; it is a _representational_ map, not a resolution map. Its role in this chapter is foundational: it normalises the idea that _working in a different coordinate_ is a legitimate, everyday computational move. With this preparation, we can now treat the more ambitious maps of the remaining sections without conceptual anxiety.

== The Calculus of One-Dimensional Mapped Methods <sec-ct-calculus>

Let $y = f(x)$ be a smooth, invertible map from the computational interval to the physical interval. The chain rule gives
$ dif / dif y = 1 / (f'(x)) dif / dif x, $ <eq-ct-chain>
and by iteration
$ dif^2 / dif y^2 = 1 / ((f'(x))^2) dif^2 / dif x^2 - (f''(x)) / ((f'(x))^3) dif / dif x. $ <eq-ct-chain-2>
When a Fourier cosine series on $x in [0, pi]$ is mapped so as to create a new basis, the orthogonality relation is preserved with the Jacobian weight @Boyd2000,
$ integral_0^pi cos(m x) cos(n x) dif x = integral_(f(0))^(f(pi)) (phi_m (y) phi_n (y)) / (f'(f^(-1) (y))) dif y, $ <eq-ct-ortho>
where $phi_n (y) equiv cos(n f^(-1) (y))$ is the mapped basis. The mapped quadrature rule is the image of the original equispaced rule under $f$: the nodes are unevenly distributed in $y$, but the _weights are all equal_ (in the midpoint-or-trapezoid sense). This last observation is surprisingly consequential and often underappreciated: much of the "complexity" of mapped spectral methods lives in the Jacobian weight, but the quadrature itself is trivial.

=== A Reusable Abstraction

To make the rest of the chapter cumulative rather than episodic, we introduce a small data structure that carries a map and its first two derivatives together, and builds the mapped derivative matrices from the underlying Chebyshev differentiation matrix.

In Python, a lightweight dataclass suffices:

```python
@dataclass
class Map1D:
    forward:      callable   # x -> y
    inverse:      callable   # y -> x
    fprime:       callable   # dy/dx
    fdoubleprime: callable   # d2y/dx2

    def derivative_matrices(self, Dx):
        x = np.cos(np.pi * np.arange(Dx.shape[0]) / (Dx.shape[0] - 1))
        fp, fpp = self.fprime(x), self.fdoubleprime(x)
        Dy = np.diag(1.0 / fp) @ Dx
        Dy2 = np.diag(1.0 / fp**2) @ (Dx @ Dx) - np.diag(fpp / fp**3) @ Dx
        return Dy, Dy2
```

In MATLAB we encapsulate the map into a struct with function handles; in Julia we use a `struct Map1D` with field-of-type-`Function`. The three-language codebook accompanying this chapter makes the same pattern available everywhere.

=== Computational Étude 19.3: Build a Reusable Map Toolkit <etude-ct-toolkit>

We validate the toolkit on two manufactured cases: the algebraic semi-infinite map $y = L(1 + x) \/ (1 - x)$ applied to $u(y) = exp(-y)$, and the tanh map $y = tanh(x)$ applied to $u(y) = 1 \/ (1 + y^2)$. In each case, first and second derivatives of $u$ are known analytically, and the toolkit should reproduce them to spectral accuracy.

#figure(
  image("../figures/ch19/python/map1d_toolkit.pdf", width: 90%),
  caption: [Étude 19.3: reusable map toolkit. Left: physical grids at $N = 24$ for the standard Chebyshev--Gauss--Lobatto nodes $x_j$ (navy), the algebraic semi-infinite image $y_j = 2 (1 + x_j) \/ (1 - x_j)$ (coral, clipped at $y = 12$), and the tanh image $y_j = tanh(x_j)$ (teal). Right: max-norm error of the first (solid) and second (dashed) mapped derivatives for both maps, validated against known analytic derivatives of $exp(-y)$ and $1 \/ (1 + y^2)$. Both first derivatives reach machine precision at $N approx 64$; the second derivatives converge with a roughly $N$-squared lag, as predicted by the higher condition number of $bold(D)^2$.],
) <fig-ct-toolkit>

Source files: `codes/python/ch19/map1d_toolkit.py`, `codes/matlab/ch19/map1d_toolkit.m`, and `codes/julia/ch19/map1d_toolkit.jl`. Every subsequent étude reuses this abstraction.

== Mapping Infinity to a Finite Interval <sec-ct-semi-inf>

On an unbounded domain, three strategies compete: _native unbounded bases_ (Hermite, Laguerre, sinc), _domain truncation_ to a large but finite $[0, L]$, and _coordinate mapping_ back to $[-1, 1]$. The chapter's focus is the third, but the student should understand the trade-offs of all three before choosing.

=== Why Truncation is Not Free

If the target solution decays exponentially as $y arrow infinity$, then truncation to $[0, L]$ introduces two distinct error sources @Boyd2000: the _spectral series error_ $E_S$ of approximating the solution on $[0, L]$, which decays geometrically with $N$, and the _domain truncation error_ $E_(D T) (L)$, which is exponentially small in $L$ but does not vanish at fixed $L$ however large $N$ becomes. Because both must shrink together, the total error is _subgeometric_: to drive the error below a prescribed tolerance one must grow $L$ and $N$ in tandem, and the effective convergence rate is the slower of the two.

#block(stroke: 0.6pt + rgb(180, 0, 0), radius: 3pt, inset: 10pt,
  fill: rgb(254, 248, 248))[
  *Warning (truncation creates two errors).* On an unbounded domain, truncation to a finite interval yields _two_ error components that must shrink together: the spectral-series error (geometric in $N$) and the domain-truncation error (exponential in $L$). Treating the two as independent is a standard beginner's mistake.
]

=== Algebraic Maps and Rational Chebyshev Functions

For a function decaying exponentially at infinity, @Boyd2000 argues that the _algebraic_ map
$ y = L (1 + x) / (1 - x), quad x in [-1, 1], quad y in [0, infinity), $ <eq-ct-algebraic>
is asymptotically superior to the logarithmic alternative $y = -L log(1 - x)$, even though the logarithmic map can look competitive at moderate $N$. The reason lies in the singularity structure of the mapped problem in the complex plane: the algebraic map places the nearest map-induced singularity at $x = 1$ as a simple pole, whereas the logarithmic map places a branch point at $x = 1$ and a slower-decaying tail of coefficients. The Chebyshev polynomials $T_n (x)$ pulled back through @eq-ct-algebraic are called _rational Chebyshev functions_ $T L_n (y) equiv T_n ((y - L) \/ (y + L))$, and their pseudospectral grid clusters near $y = 0$ as well as fanning out toward infinity.

The parameter $L$ controls the clustering. Small $L$ crowds the grid near the origin; large $L$ starves the origin and spreads the points out. Tuning $L$ is unavoidable, and the right answer depends on the characteristic decay length of the solution.

=== Computational Étude 19.4: Infinity Without Truncation <etude-ct-semi-inf>

We solve the semi-infinite benchmark
$ u''(y) - u(y) = 0, quad u(0) = 1, quad u(y) arrow 0 "as" y arrow infinity, $ <eq-ct-semi-inf-bvp>
whose exact solution is $u_("ex")(y) = exp(-y)$. The two methods are truncation to $[0, L]$ with $u(L) = 0$, and Chebyshev collocation under @eq-ct-algebraic. The map parameter $L$ is swept independently for each method.

In Python the solver is a one-liner over the mapped derivative matrices of Étude 19.3:

```python
def solve_algebraic_map(N, L):
    D, x = cheb_matrix(N)
    fp  = 2.0 * L / (1.0 - x)**2
    fpp = 4.0 * L / (1.0 - x)**3
    Dy  = np.diag(1.0/fp) @ D
    Dy2 = np.diag(1.0/fp**2) @ (D @ D) - np.diag(fpp/fp**3) @ D
    y = L * (1.0 + x) / (1.0 - x)
    A = Dy2[1:N, 1:N] - np.eye(N-1)          # interior operator
    rhs = -Dy2[1:N, N] * 1.0                  # Dirichlet data at x=-1 (y=0)
    u = np.zeros(N + 1); u[N] = 1.0
    u[1:N] = np.linalg.solve(A, rhs)
    return y, u
```

The equivalent MATLAB and Julia routines appear in the codebook.

#figure(
  image("../figures/ch19/python/semi_infinite_compare.pdf", width: 100%),
  caption: [Étude 19.4: truncation versus algebraic mapping for $u'' - u = 0$ on $[0, infinity)$. Left: solution at $N = 24$ for truncation at $L = 20$ (coral circles) and algebraic mapping with $L = 2$ (teal squares, clipped at $y = 12$). Middle: truncation error vs $N$ for three truncation lengths; each curve plateaus at the domain-truncation-error level $exp(-L)$. Right: algebraic-map error vs $N$ for four map parameters; each curve descends geometrically to machine precision with no accuracy floor. At fixed $N = 48$, the best mapped parameter achieves error $tilde.op 10^(-12)$ while the best truncation achieves error $tilde.op 10^(-8)$.],
) <fig-ct-semi-inf>

Source files: `codes/python/ch19/semi_infinite_compare.py`, `codes/matlab/ch19/semi_infinite_compare.m`, and `codes/julia/ch19/semi_infinite_compare.jl`.

=== Verdict

The algebraic map eliminates the domain-truncation error entirely; there is no truncation to commit. The map parameter $L$ must be tuned, but Figure 19.4 shows that the valley of good $L$ is broad, and a reasonable a-priori guess based on the decay length of the solution suffices. For this problem the optimal $L$ is roughly the decay length of the solution itself. Grown-up terminology: one pays for the map with a pair of extra diagonal matrix products (the terms $1 \/ f'$ and $f'' \/ (f')^3$) and gains in exchange a genuinely unbounded computation that converges geometrically in $N$.

== Weak Endpoint Singularities and Exponential Clustering <sec-ct-endpoint-sing>

Now we turn from _where_ the function lives to _how_ it behaves. Consider the elementary target
$ g(X) = sqrt(1 - X^2), quad X in [-1, 1]. $ <eq-ct-sqrt>
The function is bounded, smooth in the open interval, and innocuously symmetric. Its derivatives at $X = plus.minus 1$, however, blow up; the function has square-root branch points at the endpoints. The Chebyshev series for $g(X)$ is known in closed form @Boyd2000,
$ sqrt(1 - X^2) = 2 / pi { 1 - sum_(n=1)^infinity 2 / (4 n^2 - 1) T_(2 n) (X) }, $
with coefficients decaying only as $cal(O)(1 \/ n^2)$ --- _algebraically_. The max-norm error of an $N$-term Chebyshev truncation is therefore $cal(O)(1 \/ N)$: to achieve $10^(-12)$ accuracy one needs $N tilde.op 10^(12)$, which is hopeless.

=== The Tanh Map Changes Everything

The essential observation is due to Stenger (1981) and sharpened by @Boyd2000: if a mapping has $dif X \/ dif y$ decaying exponentially fast as $|y| arrow infinity$, then the image of a boundary layer in $X$ near the endpoints becomes a region of exponentially fine spacing in $X$, which is exactly what a weak endpoint singularity demands. The canonical such map is
$ X = tanh(y), quad y in (-infinity, infinity). $ <eq-ct-tanh>
Under this map,
$ g(tanh y) = sqrt(1 - tanh^2 (y)) = "sech" (y), $
which is an _entire_ function of $y$ with exponential decay at both infinities. A rational Chebyshev basis on $y in (-infinity, infinity)$ (or equivalently a Chebyshev basis on $y$ truncated to a wide window) recovers _geometric_ convergence. The numerical improvement is dramatic.

=== A Crossover Warning

Before the student embraces the tanh map universally, @Boyd2000 warns of a _crossover phenomenon_: for sufficiently weak singularities, the map's induced variable coefficients slow convergence at modest $N$, and only for $N > N_("cross")$ does the map begin to win. The practical advice is: if the singularity is weak and only modest accuracy is needed, try the unmapped method first. If high accuracy is required, the map is indispensable.

#block(stroke: 0.6pt + rgb(180, 0, 0), radius: 3pt, inset: 10pt,
  fill: rgb(254, 248, 248))[
  *Warning (asymptotics begin late).* For weak endpoint singularities there can be a resolution $N_("cross")$ below which the mapped method is worse than the unmapped one. The exponentially fine clustering that heals the singularity is asymptotic; its benefits emerge only when $N$ is large enough to see the exponential tail of $dif X \/ dif y$.
]

=== Computational Étude 19.5: Healing the Square-Root Branch Point <etude-ct-endpoint>

We compare two expansions of $g(X) = sqrt(1 - X^2)$:

- *Direct*: Chebyshev coefficients sampled at the Chebyshev--Gauss--Lobatto nodes of $X in [-1, 1]$.
- *Tanh-mapped*: Chebyshev coefficients sampled at a wide window $y in [-Y, Y]$ of $"sech"(y)$, with $Y = 10$ (so that $"sech"(plus.minus Y) < 10^(-8)$ and the truncation error is negligible).

In Python:

```python
def chebyshev_coeffs(v):
    N = len(v) - 1
    V = np.concatenate([v, v[N-1:0:-1]])
    a = np.real(np.fft.fft(V)) / N
    a[0] *= 0.5; a[N] *= 0.5
    return a[:N+1]

_, x = cheb_matrix(N); a_direct = chebyshev_coeffs(np.sqrt(1 - x**2))
_, xi = cheb_matrix(N); a_mapped = chebyshev_coeffs(1 / np.cosh(10 * xi))
```

The MATLAB and Julia versions follow the same template.

#figure(
  image("../figures/ch19/python/heal_branch_point.pdf", width: 100%),
  caption: [Étude 19.5: healing the square-root branch point. Left: the target function in physical coordinates $sqrt(1 - X^2)$ (coral) and in mapped coordinates $"sech" y$ (teal). Middle: max-norm error vs $N$ on a log-log scale; the direct expansion tracks $1 \/ N$ (dotted guide line), while the tanh-mapped expansion descends geometrically, reaching $10^(-9)$ at $N = 128$. Right: Chebyshev coefficients at $N = 64$; the direct expansion shows the algebraic envelope $|a_n| tilde.op 1 \/ n^2$ of the square-root series, while the tanh-mapped expansion is geometric with a decay rate of roughly $exp(-n \/ 3)$.],
) <fig-ct-endpoint>

Source files: `codes/python/ch19/heal_branch_point.py`, `codes/matlab/ch19/heal_branch_point.m`, and `codes/julia/ch19/heal_branch_point.jl`.

=== Verdict

The map $X = tanh y$ converts an algebraically-convergent expansion into a geometrically-convergent one, for exactly the reason Boyd identifies: the function $"sech" y$ is analytic on the entire real line (and in a horizontal strip of the complex $y$-plane), whereas $sqrt(1 - X^2)$ is not analytic at $X = plus.minus 1$. The map did not make the singularities vanish; it moved them to $y = plus.minus i pi \/ 2$, where they no longer obstruct the convergence of the Chebyshev basis.

== Corner Singularities: Tensor-Product Clustering First <sec-ct-corners>

Weak branch points in one dimension invite a one-dimensional map; in two dimensions they invite a two-dimensional analogue. The standard test is the Poisson problem
$ - Delta u = 1 "on" Omega = [-1, 1]^2, quad u = 0 "on" partial Omega, $ <eq-ct-poisson-sq>
whose solution has logarithmic-type branch points of the form $r^2 log r$ at the four corners @Boyd2000. Two response strategies compete: full two-dimensional _conformal_ mapping, which couples $x$ and $y$ into a single analytical transformation, and _tensor-product_ one-dimensional mapping, which applies independent one-dimensional maps in each coordinate. The former is intellectually satisfying, sometimes necessary, and almost always painful; the latter is simple, general, and frequently enough.

#block(stroke: 0.8pt + rgb(20, 45, 110), radius: 3pt, inset: 10pt,
  fill: rgb(248, 250, 254))[
  *Principle (practical first, analytic second).* For weak corner branch points, tensor-product one-dimensional clustering is almost always the right first attempt. Reach for full two-dimensional conformal maps only if (a) the singularity is strong, (b) the singularity structure is known analytically, and (c) the extra accuracy is genuinely worth the implementation cost.
]

=== Computational Étude 19.6: A Square Poisson Problem with Corner Stress <etude-ct-corners>

We solve @eq-ct-poisson-sq on an $(N+1) times (N+1)$ tensor-product Chebyshev grid with two formulations:

- *Unmapped*: standard Chebyshev--Gauss--Lobatto grid in both coordinates.
- *Tanh-clustered*: each coordinate is transformed by $X = tanh(alpha xi) \/ tanh(alpha)$ so that the grid clusters exponentially toward the walls; the mapped Laplacian uses the mapped first- and second-derivative matrices of Étude 19.3 applied in Kronecker-sum form.

Because the exact solution has no closed form, we benchmark both methods against a highly refined reference computed at $N_("ref") = 96$ on the unmapped grid. The Python assembly uses the Kronecker-sum trick of @ch-spectral-pde:

```python
D, x = cheb_matrix(N); D2 = D @ D
I = np.eye(N + 1)
L = np.kron(D2, I) + np.kron(I, D2)
# boundary/interior decomposition as in Chapter 10
```

for the unmapped problem, and the mapped version replaces `D2` by the transformed operator from the toolkit.

#figure(
  image("../figures/ch19/python/corner_tensor_clustering.pdf", width: 100%),
  caption: [Étude 19.6: square Poisson with corner stress. Left: contour plot of the solution at $N = 32$, showing the roughly circular level sets and the mild but unmistakable corner flattening. Middle: one-dimensional view of the physical grids at $N = 24$; the standard Chebyshev grid (navy) and the tanh-clustered grid with $alpha = 2$ (coral) differ only near the endpoints but differ sharply there. Right: max-norm error against the reference, for the unmapped method and for tanh clustering with $alpha in {1, 2, 3}$. The unmapped method _wins_ at every resolution tested up to $N = 64$, a striking confirmation of Boyd's crossover warning.],
) <fig-ct-corners>

=== Verdict: A Crossover Result in Two Dimensions

The Figure 19.6 right panel deserves careful reading. The unmapped Chebyshev method outperforms all three tanh-clustered variants across the entire resolution range $12 lt.eq.slant N lt.eq.slant 64$. This is not a failure of the mapping; it is the crossover phenomenon of @sec-ct-endpoint-sing in action. The corner singularity here is of type $r^2 log r$, giving an algebraic convergence index of roughly six for the unmapped Chebyshev expansion, and the clustering's exponential payoff only begins to dominate beyond $N approx 100$. At the resolutions of practical interest, the added variable coefficients of the mapped operator slow the method down more than the improved clustering speeds it up. @Boyd2000 anticipated this result precisely: "it is important to apply asymptotic concepts _asymptotically_."

The pedagogical lesson is the one we labelled a principle above: _try the unmapped method first_. In this particular corner problem, high accuracy can be reached without any mapping at all, and the extra algorithmic complexity of the mapping is pure cost. The mapping becomes useful only when the unmapped method has genuinely plateaued --- which, for weak singularities, is a long way off.

Source files: `codes/python/ch19/corner_tensor_clustering.py`, `codes/matlab/ch19/corner_tensor_clustering.m`, and `codes/julia/ch19/corner_tensor_clustering.jl`.

=== An Aside on Conformal Maps

For _strong_ corner singularities --- the L-shaped membrane being the classical example --- one-dimensional clustering is not enough, and the full two-dimensional conformal map $z mapsto z^(2 \/ 3)$ that straightens the $270 degree$ interior corner @Boyd2000 becomes unavoidable. The implementation introduces a curved mapped boundary and a "boundary factor" $Phi(u, v)$ that carries the known power-law behaviour; the resulting mapped problem is polynomial in a double Chebyshev series on the unit square, and Mason (1967) obtained five-decimal-place accuracy on an $81 times 81$ grid. An alternative is _singular basis enrichment_: augment the polynomial basis with a small number of $r^(2 \/ 3) sin(2 theta \/ 3)$-type functions. @Boyd2000 warns that too many weak singular basis functions destroy conditioning, because they become numerically indistinguishable from combinations of ordinary polynomials. We leave this as a starred exercise; the main text of this chapter does not dwell on it.

== Periodic Concentration: The Arctan/Tan Map <sec-ct-arctantan>

We have already met the arctan/tan map in the opening étude. Here we promote it from a demonstration trick to a well-characterised tool with a parameter that must be tuned deliberately. Boyd's original form,
$ y = arctan(L tan x), quad x, y in [0, pi], $ <eq-ct-arctantan-pi>
is period-$pi$; the period-$2 pi$ variant of @eq-ct-arctantan is pedagogically equivalent and slightly more convenient for Fourier practitioners.

=== Why the Map is so Convenient

Five properties of the arctan/tan map together account for its popularity @Boyd2000:

1. The forward and inverse maps are both explicit elementary functions, so evaluation back to the physical coordinate is free.
2. The derivative $dif y \/ dif x$ is a rational trigonometric polynomial in $x$, so differential equations with trigonometric-polynomial coefficients in $y$ transform into differential equations with trigonometric-polynomial coefficients in $x$. The map is _polynomial-coefficient preserving_.
3. The map is smooth (infinitely differentiable) everywhere on the real axis, so no new boundary layers or rapid-variation regions are introduced.
4. The width parameter $L$ has a clear geometric meaning: it is the slope of the map at the centre. Small $L$ concentrates the grid near the centre; large $L$ spreads it out.
5. The map is periodic in the computational coordinate with the same period as the physical coordinate, so the Fourier basis is preserved.

=== Computational Étude 19.7: A Localised Periodic Pulse in the Right Coordinate <etude-ct-arctan-sweep>

We approximate the same pulse @eq-ct-pulse as in Étude 19.1, now sweeping the map parameter $L$ over $[0.08, 1.5]$ at each $N in {12, ..., 96}$ to produce an error landscape. The practical question for the student is: how sensitive is the method to the choice of $L$? If the error depends delicately on $L$, parameter tuning becomes a research problem in itself; if the error is insensitive across a broad valley, the map is ergonomically friendly.

The Python evaluation loop is compact:

```python
for N in Ns:
    for L in Ls:
        x = -np.pi + 2*np.pi*np.arange(N)/N
        y = 2*np.arctan(L*np.tan(x/2))
        c = np.fft.fft(target(y)) / N
        x_eval = 2*np.arctan(np.tan(y_eval/2) / L)
        # synthesis via FFT coefficients at the inverse-mapped points
```

#figure(
  image("../figures/ch19/python/arctan_tan_sweep.pdf", width: 100%),
  caption: [Étude 19.7: parameter sweep for the arctan/tan map. Left: three mapped grids at $N = 32$ for $L = 0.1$ (coral), $L = 0.3$ (teal), $L = 1.0$ (orange); the profile $f(y)$ is shown for reference. Middle: error landscape in $(N, L)$-space; a broad dark-blue valley of good parameter choices is clearly visible. Right: slices at fixed $N$; the optimal $L$ ranges from roughly $0.15$ at $N = 16$ to roughly $0.45$ at $N = 64$, and the valley is decades wide at every resolution. The method is not sensitive to fine tuning of $L$.],
) <fig-ct-arctan-sweep>

=== Verdict

The optimal $L$ is a weak function of $N$, drifting slowly from $approx 0.15$ at $N = 16$ to $approx 0.45$ at $N = 64$, and the error valley is broad enough that a single fixed choice of $L$ (say $L = 0.3$) performs well across the entire resolution range. The method reaches machine precision at $N approx 64$, a factor-of-two saving over the unmapped Fourier grid of Étude 19.1. In a cubic-cost solve, that is a factor-of-eight saving in work.

Source files: `codes/python/ch19/arctan_tan_sweep.py`, `codes/matlab/ch19/arctan_tan_sweep.m`, and `codes/julia/ch19/arctan_tan_sweep.jl`.

== Adaptive Mappings for Moving Fronts <sec-ct-adaptive>

All the maps so far have been static: the map is chosen once and the grid stays put. For problems in which the region of rapid change _moves_ (a propagating front in a reaction-diffusion equation, a steepening shock in viscous Burgers, a coherent structure in turbulent convection), a static map can only do so much. The ambition of _adaptive_ mapping is to let the map's parameters evolve with the solution.

@Boyd2000 discusses two strategies. The first, arclength adaptation, re-parametrises along the solution curve; the second, parametric adaptation, fixes a functional form (typically arctan/tan) and lets its width $L(t)$ and centre $y_f(t)$ evolve in time, chosen at each update by minimising a smoothness functional such as
$ I(y_f, L) = integral_0^pi |u_(x x)|^2 + |u_x|^2 dif x. $ <eq-ct-smoothness>

The cost of adaptivity is substantial: the smoothness minimisation is a small nonlinear optimisation problem, and the transfer of the solution from the old grid to the new grid requires _off-grid_ interpolation, which cannot use the FFT (@Boyd2000, Rule of Thumb: off-grid interpolation costs roughly four FFTs' worth of work). A useful practical compromise is to update the map only every ten or so timesteps, not at every step. Within those ten steps, the front does not drift far enough to demand a fresh grid.

This section is presented as a capstone rather than as part of the routine toolkit. The student should understand the idea and the cost, and should _not_ reach for adaptivity as a first resort. A carefully chosen static map is almost always cheaper and, in a well-posed problem with a slowly-varying front location, almost as accurate.

#block(stroke: 0.8pt + rgb(20, 45, 110), radius: 3pt, inset: 10pt,
  fill: rgb(248, 250, 254))[
  *Principle (update the map only when necessary).* An adaptive map should be treated as an _occasional_ correction to a static one, not as a compulsory action at every timestep. In a well-designed code, the map is re-computed only when a smoothness diagnostic exceeds a preset threshold.
]

A full computational étude for adaptive Burgers integration is a substantial project and is left as an end-of-chapter exercise. The infrastructure --- the `Map1D` toolkit from Étude 19.3 and the parameter-sweep logic from Étude 19.7 --- is already in place; the student need only wrap them in a time loop.

== The Almost-Equispaced Grid: The Kosloff--Tal-Ezer Map <sec-ct-kte>

We close the chapter with the map most likely to mislead the unwary practitioner: the Kosloff--Tal-Ezer (KTE) arcsine map, proposed to cure the notorious $cal(O)(N^2)$ time-step restriction of Chebyshev pseudospectral methods. The map is
$ y = arcsin((1 - beta) x) / arcsin(1 - beta), quad x in [-1, 1], $ <eq-ct-kte>
which makes the image of the Chebyshev grid far more uniform in $y$ than the original grid. When $beta$ is chosen to depend on $N$ like $beta = C \/ N^2$, the minimum inter-grid spacing becomes $cal(O)(1 \/ N)$ (Boyd's Theorem 30), matching the spacing of an equispaced finite-difference grid. The explicit time-step restriction of an $cal(O)(1 \/ N^2)$ problem therefore becomes $cal(O)(1 \/ N)$, a dramatic improvement.

=== The Hidden Price

The price is paid in spectral accuracy. The map @eq-ct-kte has branch-point singularities at
$ t_s = m pi plus.minus i op("arccosh")(1 - beta), quad m in ZZ, $ <eq-ct-kte-branch>
in the trigonometric coordinate $t = arccos x$. For small $beta$, the imaginary part of the nearest branch point is approximately
$ op("Im")(t_s) approx sqrt(2 beta), $
so when $beta = C \/ N^2$ the nearest singularity sits at imaginary distance $approx sqrt(2 C) \/ N$ from the real axis. The geometric convergence factor for a Chebyshev (equivalently, Fourier cosine) series is $exp(-N op("Im")(t_s))$, so with this scaling
$ exp(-N op("Im")(t_s)) approx exp(-sqrt(2 C)) + cal(O)(1 \/ N). $
In words: the exponential convergence factor _does not decrease with $N$_. Spectral accuracy has been sacrificed to obtain the larger time-step. @Boyd2000 is emphatic on this point.

#block(stroke: 0.6pt + rgb(180, 0, 0), radius: 3pt, inset: 10pt,
  fill: rgb(254, 248, 248))[
  *Warning (do not buy timesteps by selling spectral accuracy).* A map that improves stiffness but destroys exponential coefficient decay is not automatically a good spectral map. The only question worth asking of the KTE map is: what conservative choice of $beta$ gives a useful time-step improvement _without_ sacrificing spectral accuracy?
]

Hesthaven, Dinesen, and Lynov (1999) chose $beta = 1 - cos(1 \/ 2)$, an $N$-_independent_ constant. This choice gives only a factor of about two improvement in the time-step, but it preserves geometric convergence in full.

=== Computational Étude 19.8: Accuracy versus Timestep <etude-ct-kte>

We reproduce Boyd's Figure 16.3 in full, plotting three things as functions of $N$:

1. the minimum grid spacing $min_j |y_(j+1) - y_j|$ under the standard, aggressive KTE ($beta = 1 - cos(1 \/ N)$), and conservative KTE ($beta = 1 - cos(1 \/ 2)$);
2. the spectral radius $rho(bold(D))$ of the first-derivative matrix, which controls the explicit time-step limit;
3. the Chebyshev coefficient magnitude $|a_N|$ of the simplest test function $f(y) = y$, which in the absence of any mapping is identically $T_1$ and has $|a_1| = 1$ and all other coefficients zero.

In Python:

```python
def kte_grid(N, beta):
    D, xi = cheb_matrix(N)
    denom = np.arcsin(1.0 - beta)
    y = np.arcsin((1.0 - beta) * xi) / denom
    fp = (1.0 - beta) / (np.sqrt(1 - (1-beta)**2 * xi**2) * denom)
    Dy = np.diag(1.0 / fp) @ D
    return y, xi, D, Dy
```

#figure(
  image("../figures/ch19/python/kosloff_tal_ezer.pdf", width: 100%),
  caption: [Étude 19.8: accuracy versus timestep under the Kosloff--Tal-Ezer map. Left: minimum grid spacing versus $N$ for the standard Chebyshev grid (navy, $cal(O)(1 \/ N^2)$ scaling), aggressive KTE with $beta tilde.op 1 \/ N^2$ (coral, $cal(O)(1 \/ N)$ scaling), and conservative KTE with $beta = 1 - cos(1 \/ 2)$ (teal). Middle: spectral radius of the first-derivative matrix; aggressive KTE reduces stiffness by nearly an order of magnitude at $N = 96$. Right: Chebyshev coefficient $|a_N|$ of $f(y) = y$ sampled on each grid, with the theoretical asymptote $0.488 \/ N^2$ (navy dotted) confirming Boyd's analysis: aggressive KTE destroys geometric convergence for the simplest polynomial, while conservative KTE preserves it down to the floating-point floor.],
) <fig-ct-kte>

Source files: `codes/python/ch19/kosloff_tal_ezer.py`, `codes/matlab/ch19/kosloff_tal_ezer.m`, and `codes/julia/ch19/kosloff_tal_ezer.jl`.

=== Verdict: the Chapter's Final Lesson

The Kosloff--Tal-Ezer map is the chapter's _cautionary climax_. A student who writes a finite-difference code with a Chebyshev grid and sees a dramatic improvement in timestep at modest $beta$ may be tempted to push $beta$ aggressively; the étude shows that doing so destroys spectral accuracy at the simplest test, the linear polynomial. The only defensible choice is the conservative one: fix $beta$ independently of $N$, accept a modest (factor-of-two to factor-of-four) timestep improvement, and retain geometric convergence. The aggressive $beta = cal(O)(1 \/ N^2)$ choice is correct only if one is willing to be a finite-difference method in disguise.

The deeper lesson generalises beyond KTE: _a more uniform grid is not automatically a more accurate grid_. Spectral accuracy lives in the analytic extension of the basis into the complex plane, and any map that pulls singularities toward the approximation interval --- no matter how innocent the grid picture --- pays a corresponding price in convergence. Every mapped spectral method in this chapter has been scored against the three questions of @sec-ct-prelude; the KTE map scores well on questions 1 and 2 but fails question 3 dramatically when $beta$ scales with $N$.

== A Decision Guide for Coordinate Transformations <sec-ct-decision>

We close the chapter with a compact decision tree that summarises its contents operationally.

#block(stroke: 0.8pt + rgb(20, 45, 110), radius: 4pt, inset: 12pt,
  fill: rgb(248, 250, 254))[
  *Mapping decision guide.* For any computational problem requiring a spectral discretisation, ask, in order:

  1. *What kind of domain?*  Bounded, semi-infinite, infinite, periodic.  Each class has a default basis: Chebyshev, rational Chebyshev, rational Chebyshev or sinc, Fourier.
  2. *What is the pathology?*  Is the difficulty at infinity, at an endpoint, at a corner, inside a narrow pulse, or in the timestep?  Each pathology has a preferred cure: algebraic mapping, exponential clustering, tensor-product clustering (or conformal map for strong singularities), arctan/tan, KTE with conservative parameter.
  3. *Is the singularity weak or strong?*  If weak, _try the unmapped method first_; a crossover $N$ may exist only out of reach.  If strong, the mapping is indispensable.
  4. *What parameter must be tuned?*  Every useful map has at least one: the algebraic $L$, the tanh clustering strength $alpha$, the arctan/tan width $L$, the KTE conservativeness $beta$.  Sweep the parameter and verify that the valley of good choices is broad.
  5. *What diagnostic certifies success?*  Coefficient decay, max-norm error against an independent reference, condition number of the mapped operator, or observed timestep?  The choice depends on the pathology; write it down _before_ running the code.
  6. *Are the map's singularities safely away from the approximation interval?*  This is the KTE question.  Compute the nearest map-induced singularity in the complex plane and compare its imaginary distance to $1 \/ N$.
]

The guide is not a substitute for judgement. Problems that mix pathologies --- an infinite domain with a sharp front moving across it, or a strong corner singularity in an unbounded region --- require chaining maps together, and the order of composition matters. The habit cultivated by this chapter is, at every stage, to stop and ask the three questions.

== A Non-Exhaustive Literature Review <sec-ct-lit-review>

Coordinate transformations are one of the oldest topics in numerical analysis of differential equations, and one of the most recently renovated. The classical reference for the whole subject as practiced in spectral methods is @Boyd2000 Chapter 16, from which the entire pedagogical arc of this chapter descends. @Trefethen2000 makes the Chebyshev-to-cosine map the organising principle of spectral-methods pedagogy, and @Trefethen2013 deepens the picture into the approximation-theory literature on analytic extension and nearest singularities.

The algebraic mapping of semi-infinite and infinite intervals is elaborated in the classical paper of Grosch and Orszag (1977), with the rational Chebyshev functions $T L_n (y)$ promoted to first-class status in @Boyd2000. For the doubly-infinite case, @WeidemanTrefethen2007 analyse the sinc- and rational-Chebyshev families side by side and settle several practical questions of parameter choice.

The tanh-map analysis of weak endpoint singularities originates with Stenger (1981) and is sharpened in the spectral-method context by @Boyd2000. For L-shaped membrane and other corner-singularity problems the landmark paper is Mason (1967), with a modern spectral-element treatment in Nektar++ @MoxeyNektar2020 and a tangential but illuminating application in the "lightning solver" series @GopalTrefethen2019 @HerremansHuybrechsTrefethen2023, which use _rational_ basis functions (tight clusters of poles exterior to the domain) to achieve spectral accuracy at corners without any mapping whatsoever. The "lightning" viewpoint is in some sense the complement of the mapped one, and students who complete this chapter should read those papers to see the same problem solved from the other end.

The arctan/tan map for periodic concentration appears in @Boyd2000 as one of the earliest success stories of pedagogically-motivated mapping; it has since been adapted for Laplace problems in periodic domains by @Baddoo2021 via his AAAtrig algorithm (a topic treated in the literature review of @ch-linear-eigen). Adaptive mapping for moving fronts has an extensive literature surveyed in Boyd's Table 16.2; for modern differentiable-solver treatments the reader is referred to the Dedalus ecosystem @Burns2020 @Lecoanet2026.

The Kosloff--Tal-Ezer map was introduced in Kosloff and Tal-Ezer (1993), analysed in detail by Tal-Ezer (1994), and exhaustively compared with alternatives in the papers catalogued in Boyd's Table 16.4. The conservative choice $beta = 1 - cos(1 \/ 2)$ advocated in this chapter is due to Hesthaven, Dinesen, and Lynov (1999); the subsequent literature has largely confirmed this as the sweet spot for applications that cannot afford to sacrifice spectral accuracy.

Finally, the ultraspherical spectral method @OlverTownsend2013 @OlverTownsend2019 @Xu2024 deserves mention as a modern competitor to the mapping approach: rather than transforming the grid to match the solution, the ultraspherical method transforms the _basis_ to give almost-banded differential operators with bounded condition number. The two strategies are not in opposition; for many problems the best approach is to apply a coordinate map followed by an ultraspherical discretisation of the mapped operator. A recent extension by @OlverTownsendVasil2019 carries sparse ultraspherical discretisation onto the _triangle_, an important geometric primitive that previously demanded a body-fitted conformal map; the basis-centric approach delivers linear-complexity solvers for bi-Laplacian and similar high-order operators on non-rectangular polygons without any physical-domain mapping whatsoever.

The period 2021--2026 has also produced a set of developments that refine, and in some cases overturn, the cautionary story we have been telling about the Kosloff--Tal-Ezer map. @DeMarchiMarchettiPerracchione2023 introduce a _Kosloff--Tal-Ezer least-squares_ (KTL) quadrature framework that decouples the polynomial-approximation degree from the number of spatial nodes: by solving an overdetermined least-squares problem evaluated on the mapped "fake nodes", one can use an aggressive $N$-dependent $beta$ to secure the $cal(O)(1 \/ N)$ timestep scaling _without_ destroying geometric convergence. The aggressive KTE map is thereby rehabilitated from a cautionary tale into a stable, time-dependent-PDE workhorse on quasi-uniform spatial grids. This is a reminder that a pedagogical caveat can have a genuine algorithmic solution that the textbook must acknowledge.

An equally consequential development for corner singularities is the _lightning solver_ family, which bypasses coordinate mapping altogether by placing poles outside the computational domain in exponentially clustered patterns. After the Laplace and Helmholtz origins of @GopalTrefethen2019, @BrubeckTrefethen2022 extend the approach to the biharmonic equation (Stokes flow) via the Goursat-function partition and rational approximation, reaching ten-digit accuracy on polygonal Stokes flows in under a second on a laptop. @XueWatersTrefethen2024 push the methodology further with the Lightning-AAA Rational Stokes (LARS) algorithm, resolving intricate Moffatt eddies in microfluidic bifurcations with the AAA algorithm providing automated pole placement. The "log-lightning" or reciprocal-log family of @NakatsukasaTrefethen2021Recip takes the idea further still, elevating the baseline root-exponential convergence of lightning methods to near-exponential $cal(O)(exp(-C N \/ log N))$ for analytic functions with branch-point singularities. The relationship between mapping and lightning methods can be summarised as follows: mapping warps the _physical_ coordinate so the basis sees a smoother function; lightning enriches the _basis_ with rational elements that absorb the singularities where they are, in the physical coordinate, leaving the domain untouched. Both are legitimate; both work; the practitioner now has a genuine choice.

Among recent developments at the boundary between unbounded-interval methods (the subject of @ch-unbounded) and coordinate transformations, @ChenShen2022LOF introduce _log-orthogonal functions_ built by applying a logarithmic map to standard Laguerre polynomials; these functions resolve weak endpoint singularities with exponential accuracy where ordinary Laguerre expansions yield only algebraic rates. @HuangShen2024MHF develop _mapped Hermite functions_ for the multi-point weakly singular Fredholm--Hammerstein integral equations that arise in two-dimensional nonlocal problems; in both cases the map is applied to the generating polynomials rather than to the macroscopic grid, which keeps the orthogonality structure of the basis intact. @ShenWen2026 propose a _fictitious-domain spectral method_ for fourth-order PDEs in complex geometries: rather than attempting to map an irregular physical domain to a disc or square via a conformal transformation (which would produce pathologically ill-conditioned bi-Laplacian matrices), the domain is _embedded_ in a larger disc and a polar transform is applied to the fictitious domain. The idea is the same one we pushed in @sec-ct-corners --- _prefer the simpler transformation when you can_ --- taken to its logical conclusion.

Finally, the Dedalus project @Burns2020 @Lecoanet2026 has continued its rapid expansion, introducing sparse tensorial bases built on weighted Jacobi polynomials that carry the correct geometric-regularity behaviour at polar and spherical coordinate singularities automatically, without parity filtering or manual basis modification. The resulting capability to simulate stars, planetary interiors, and moist convection end-to-end on automatic adjoint-aware spectral discretisations moves the "mapping" question from a numerical-analysis problem to a software-framework one: the user specifies the geometry, and the framework picks the right mapped basis.

== Summary and Exercises <sec-ct-summary>

The chapter's slogan, to borrow and amplify Boyd's own phrasing: _choose the coordinate that makes the solution smoother, but never forget to inspect what the coordinate has done to the operator, the singularities, and the cost._

=== Conceptual Exercises

1. Prove the identity $T_n (cos t) = cos(n t)$ by induction on $n$ using the three-term recurrence.
2. Derive @eq-ct-chain-2 from @eq-ct-chain by the chain rule, and verify it on $y = tanh(x)$ applied to $u(y) = 1 \/ (1 + y^2)$.
3. Show how the $L^2$ inner product on $y in [a, b]$ transforms into an inner product on $x in [f^(-1) (a), f^(-1) (b)]$ with weight $f'(x)$.
4. Explain why truncation of an exponentially-decaying function to $[0, L]$ creates a domain-truncation error of order $exp(-L)$ in addition to the spectral-series error.
5. Show that the algebraic map @eq-ct-algebraic places its only map-induced singularity at $x = 1$, and argue why this makes it asymptotically preferable to the logarithmic alternative.
6. For $g(X) = sqrt(1 - X^2)$, identify the branch points in the complex $X$-plane and show they become horizontal lines at $op("Im")(y) = plus.minus pi \/ 2$ under $X = tanh(y)$.
7. Explain in words why exponential boundary clustering can heal a weak endpoint singularity but does not cure a strong one (e.g.\ a simple pole).
8. Derive the inverse of the arctan/tan map @eq-ct-arctantan and compute its Jacobian.
9. Show that the KTE map @eq-ct-kte has branch-point singularities at @eq-ct-kte-branch, and deduce the asymptotic relation $op("Im")(t_s) approx sqrt(2 beta)$ for small $beta$.
10. Explain why a more uniform grid can still give worse spectral convergence than a clustered grid.

=== Computational Exercises

11. Implement the `Map1D` abstraction in the language of your choice and validate it on the tanh map applied to $u(y) = "sech"^2 (y)$, with analytic first and second derivatives.
12. Repeat Étude 19.4 for the semi-infinite problem $u'' + u = 0$, $u(0) = 0$, $u(y) arrow 0$, whose exact solution is $u = 0$ (the only decaying solution). Observe that the mapping produces the zero solution to machine precision, but the unmapped truncated problem has a subtle contamination from the lower endpoint. Explain why.
13. For Étude 19.6, find the $N$ at which the $alpha = 2$ tanh-clustered method catches up with the unmapped method. How far out must you go?
14. Re-run Étude 19.7 with the even narrower pulse $kappa = 200$. At what $N$ does the mapped method with best $L$ reach machine precision? How does that compare with the unmapped method at the same $N$?
15. Construct a static arctan/tan map for the 1D viscous Burgers equation with a single steep front, and compare the error at fixed time against an unmapped Fourier method. Report the factor of work saved.
16. Implement the KTE map with $beta = 1 - cos(1 \/ 2)$ in the 1D advection equation and compare the largest stable explicit time-step against the unmapped Chebyshev method at $N = 64$. Confirm the factor-of-two improvement.

=== Project-Style Exercises

17. Build a fully adaptive solver for viscous Burgers with a moving quasi-shock, using the arctan/tan map with time-varying centre $y_f (t)$ and width $L(t)$ chosen by minimising @eq-ct-smoothness every ten timesteps. Compare total cost against a carefully-chosen static map.
18. Study the L-shaped-membrane eigenvalue problem with two approaches in sequence: first a tensor-product tanh clustering, and then a $z^(2 \/ 3)$ conformal map with boundary factor. Report which method reaches five-decimal-place accuracy with fewer degrees of freedom.
19. Design a test suite of four test problems, each with a different pathology (unbounded, weak endpoint, periodic concentration, timestep-dominated), and for each one justify, implement, and evaluate a coordinate transformation. Prepare a three-page report in the style of the chapter's décision guide.
20. _When mapping hurts._ Write a short computational note giving two examples from the chapter where the unmapped method is preferable at moderate $N$. Explain precisely why.
