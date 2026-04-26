// textbook/chapters/notation.typ
// Notation and Symbols Glossary
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Last modified: April 2026

#import "../styles/template.typ": dropcap, ii

#heading(level: 1, numbering: none)[Notation and Symbols]

#dropcap[The book uses a small but consistent body of notation. The conventions below are introduced where they are first needed, then carried throughout: a reader who dips into Ch 13 looking for advanced boundary conditions, or into Ch 18 for an eigenproblem technique, can use this glossary as a one-stop reference for the basis families, derivative matrices, and convergence quantities that the surrounding text takes for granted.] Symbols whose meaning shifts between chapters (the most common case is $N$, the resolution, which is sometimes the polynomial degree, sometimes the number of grid points, and sometimes the number of Fourier modes) are explicitly redefined wherever they appear; the entries below give the most frequent reading.

#let navy = rgb(20, 45, 110)

#let glossary(rows) = block(
  stroke: (top: 1.2pt + navy, bottom: 1.2pt + navy),
  inset: 0pt,
  table(
    columns: (auto, 1fr),
    align: (right + horizon, left + horizon),
    inset: (x: 1em, y: 0.4em),
    stroke: none,
    table.hline(stroke: 0.5pt + navy),
    ..rows.map(((sym, desc)) => (sym, desc)).flatten(),
    table.hline(stroke: 0.5pt + navy),
  ),
)

#heading(level: 2, numbering: none)[Sets, constants, and basic conventions]

#glossary((
  ($bb(R), bb(C), bb(Z), bb(N)$,           [real numbers, complex numbers, integers, non-negative integers]),
  ($ii$,                                    [upright imaginary unit, $ii^2 = -1$ (the Typst keyword `ii`; equivalent to `\mathrm{i}` in LaTeX)]),
  ($delta_(j k)$,                           [Kronecker delta, $1$ if $j = k$ else $0$]),
  ($delta(x)$,                              [Dirac delta distribution]),
  ($f tilde.op g$,                          [asymptotic equivalence, $f \/ g arrow.r 1$ in the relevant limit]),
  ($cal(O), o$,                             [Landau big-O and little-o symbols]),
  ($epsilon_("mach")$,                      [machine epsilon, $approx 2.22 times 10^(-16)$ in IEEE 754 double precision]),
  ($epsilon$,                               [a small parameter (singular perturbation, regularisation, or numerical tolerance, context-dependent)]),
  ($lambda$,                                [eigenvalue (real or complex)]),
))

#heading(level: 2, numbering: none)[Polynomials, basis functions, and special functions]

#glossary((
  ($T_n (x)$,                               [Chebyshev polynomial of the first kind, $T_n (cos t) = cos(n t)$, $x in [-1, 1]$]),
  ($U_n (x)$,                               [Chebyshev polynomial of the second kind]),
  ($P_n (x)$,                               [Legendre polynomial]),
  ($H_n (y)$,                               [Hermite polynomial (physicist's normalisation)]),
  ($L_n (y)$,                               [Laguerre polynomial]),
  ($psi_n (y)$,                             [normalised Hermite function, $psi_n = (2^n n! sqrt(pi))^(-1\/2) H_n (y) e^(-y^2 \/ 2)$]),
  ($phi_n (y)$,                             [Laguerre function, $phi_n = e^(-y \/ 2) L_n (y)$, $y in [0, infinity)$]),
  ($T B_n (y)$,                             [rational Chebyshev on $bb(R)$, $T B_n (y) = T_n (y \/ sqrt(ell^2 + y^2))$]),
  ($T L_n (y)$,                             [rational Chebyshev on $[0, infinity)$, $T L_n (y) = T_n ((y - ell) \/ (y + ell))$]),
  ($S B_n (y)$,                             [sine-rational basis on $bb(R)$, $S B_n (y) = sin{(n + 1) op("arccot")(y \/ ell)}$]),
  ($J_m (r), K_m (r)$,                      [Bessel function of the first kind; modified Bessel function of the second kind]),
  ($"erf"(x), "erfc"(x)$,                   [error function and complementary error function]),
  ($op("sinc")(z)$,                         [normalised sinc, $op("sinc")(z) = sin(pi z) \/ (pi z)$]),
))

#heading(level: 2, numbering: none)[Grids, parameters, and coordinate maps]

#glossary((
  ($N$,                                     [resolution: number of modes, polynomial degree, or grid intervals (chapter-dependent; redefined where needed)]),
  ($h$,                                     [grid spacing; $h = 2 pi \/ N$ on the standard Fourier grid]),
  ($x_j, y_j$,                              [grid points in computational ($x$) and physical ($y$) coordinate]),
  ($x_j = cos(j pi \/ N)$,                  [Chebyshev--Gauss--Lobatto nodes, $j = 0, dots, N$]),
  ($x_j = 2 pi j \/ N$,                     [equispaced periodic Fourier nodes, $j = 0, dots, N - 1$]),
  ($L$,                                     [truncation half-width on $[-L, L]$ for unbounded-domain methods]),
  ($ell$,                                   [map parameter for rational-Chebyshev / arctan-tan / sinh maps]),
  ($alpha$,                                 [Hermite scale, $alpha = sqrt(2 A)$ matches a target $e^(-A y^2)$]),
  ($beta$,                                  [Kosloff--Tal-Ezer map parameter; $beta = 1 - cos(1 \/ 2)$ is the conservative choice]),
  ($x = cos t$,                             [Chebyshev cosine map, $t in [0, pi]$]),
  ($y = ell x \/ sqrt(1 - x^2)$,            [algebraic infinite map, $bb(R) arrow.l.r [-1, 1]$]),
  ($y = ell (1 + x) \/ (1 - x)$,            [algebraic semi-infinite map, $[0, infinity) arrow.l.r [-1, 1]$]),
  ($y = tanh(x)$,                           [tanh endpoint-clustering map]),
  ($y = 2 arctan(ell tan(x \/ 2))$,         [arctan/tan periodic map (period $2 pi$)]),
  ($y = sinh(ell t)$,                       [Weideman--Cloot sinh hybrid map]),
))

#heading(level: 2, numbering: none)[Discretisation: matrices and operators]

#glossary((
  ($bold(D), D$,                            [first-derivative differentiation matrix in the working basis]),
  ($bold(D)^2, D_N^((2))$,                  [second-derivative matrix (often $bold(D) dot bold(D)$ in practice)]),
  ($bold(D)^((m))$,                         [$m$-th derivative matrix]),
  ($D_y, D_y^2$,                            [derivative matrices in physical coordinate after a coordinate map]),
  ($bold(A)$,                               [system / left-hand-side matrix of a discrete BVP, eigenproblem, or pencil]),
  ($bold(B)$,                               [right-hand-side mass matrix of a generalised eigenvalue pencil $bold(A) bold(u) = lambda bold(B) bold(u)$]),
  ($bold(M)$,                               [mass matrix in a Galerkin / weighted-residual formulation]),
  ($bold(L)$,                               [linear operator (continuous or discrete, context-dependent)]),
  ($bold(S)$,                               [scaling/Heinrichs basis-recombination factor, e.g.\ $bold(S) = "diag"(1 - x_j^2)$]),
  ($bold(u), bold(v)$,                      [discrete solution / grid-function vectors]),
  ($bold(I), I$,                            [identity matrix]),
  ($times.o$,                          [Kronecker product (used in tensor-product 2D discretisations $bold(I) times.o bold(D)^2 + bold(D)^2 times.o bold(I)$)]),
  ($cal(F), cal(F)^(-1)$,                   [continuous Fourier transform and its inverse]),
  ($"DFT", "FFT", "DCT"$,                   [discrete Fourier transform, fast Fourier transform, discrete cosine transform]),
  ($cal(B)$,                                [boundary-condition operator, $cal(B) u = 0$]),
))

#heading(level: 2, numbering: none)[Spectral coefficients, errors, and convergence]

#glossary((
  ($hat(f)_k$,                              [continuous Fourier or Chebyshev coefficient]),
  ($tilde(f)_k$,                            [discrete (DFT/aliased) coefficient computed from $N$ samples]),
  ($a_n, c_n$,                              [series coefficients in a polynomial or trigonometric expansion]),
  ($p_N (x), u_N (x)$,                      [degree-$N$ truncation of the series; spectral approximation to $u(x)$]),
  ($Lambda_N$,                              [Lebesgue constant of an interpolation node set]),
  ($rho$,                                   [Bernstein-ellipse parameter; convergence rate $cal(O)(rho^(-N))$ for analytic functions]),
  ($a$,                                     [half-width of an analyticity strip in the complex plane]),
  ($kappa(bold(A))$,                        [spectral condition number of the matrix $bold(A)$, $kappa = sigma_max \/ sigma_min$]),
  ($sigma_max, sigma_min$,                  [largest and smallest singular values]),
  ($epsilon_N, ||f - p_N||_infinity$,       [max-norm (sup-norm) approximation error at resolution $N$]),
  ($||u||_2$,                               [$L^2$ or $ell^2$ norm of $u$]),
  ($||u||_infinity$,                        [sup-norm of $u$]),
  ($cal(R)$,                                [reflection coefficient (in scattering études); residual operator (in Galerkin études)]),
  ($delta(q)$,                              [eigenvalue correction $delta(q) = lambda(q) - lambda_0$ in a parameter sweep]),
))

#heading(level: 2, numbering: none)[Time integration and stability]

#glossary((
  ($Delta t$,                               [time step]),
  ($Delta t_max$,                           [maximum stable time step (CFL limit)]),
  ($nu_k = lambda_k Delta t$,               [scaled eigenvalue used to test if a step lies inside a stability region]),
  ($cal(R)(z)$,                             [stability function of an explicit or implicit Runge--Kutta scheme]),
  ($"CFL"$,                                 [Courant--Friedrichs--Lewy condition]),
  ($"IF", "ETD", "IMEX", "RKC"$,            [integrating-factor, exponential-time-differencing, implicit-explicit, and Runge--Kutta--Chebyshev families]),
))

#heading(level: 2, numbering: none)[Differential operators and PDE acronyms]

#glossary((
  ($partial_x, u_x$,                        [partial derivative; abbreviation for $partial u \/ partial x$]),
  ($nabla, Delta$,                          [gradient and Laplacian]),
  ($nabla^4$,                               [biharmonic operator]),
  ($"BVP", "ODE", "PDE"$,                   [boundary value problem, ordinary / partial differential equation]),
  ($"KdV"$,                                 [Korteweg--de Vries equation, $u_t + 6 u u_x + u_(x x x) = 0$]),
  ($"NLS"$,                                 [nonlinear Schrödinger equation, $i u_t + u_(x x) + 2 |u|^2 u = 0$]),
  ($"KS"$,                                  [Kuramoto--Sivashinsky equation, $u_t + u u_x + u_(x x) + u_(x x x x) = 0$]),
  ($"NS"$,                                  [Navier--Stokes equations]),
  ($"QNM"$,                                 [quasinormal mode (complex eigenvalue with outgoing-wave boundary conditions)]),
))

#v(0.6em)
A handful of secondary acronyms appear once or twice and are defined inline at first use: *AAA* (adaptive Antoulas--Anderson rational approximation), *NUFFT* (non-uniform FFT), *sFFT* (sparse FFT), *PINN* (physics-informed neural network), *FNO* (Fourier neural operator), *NeRF* (neural radiance field), *KTE/KTL* (Kosloff--Tal-Ezer / Kosloff--Tal-Ezer least-squares), *LARS* (Lightning--AAA Rational Stokes), *DNS* (direct numerical simulation), *MSC* (Mathematics Subject Classification).
