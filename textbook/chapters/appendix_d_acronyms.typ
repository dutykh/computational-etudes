// textbook/chapters/appendix_d_acronyms.typ
// Appendix D: Glossary of Acronyms and Initialisms
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Last modified: April 2026

#import "../styles/template.typ": dropcap, chapter-abstract

= Glossary of Acronyms and Initialisms <appendix-acronyms>

#chapter-abstract(keywords: [Acronyms · Glossary · Notation])[
This appendix collects, in one alphabetised place, the expansions of every acronym and initialism used in the book. Drawn from numerical analysis, signal processing, fluid mechanics, quantum physics, and scientific machine learning, each term is also defined inline at its first appearance in the text; the glossary serves the non-linear reader who arrives in the middle of a chapter.
]

#dropcap[Spectral methods sit at the intersection of numerical analysis, signal processing, fluid mechanics, quantum physics, and modern scientific machine learning, and the volume's prose accordingly inherits acronyms from all five.] Every acronym used in the body is defined inline at first appearance; this glossary collects the expansions in one alphabetised place for the benefit of non-linear readers.

#let navy = rgb(20, 45, 110)

#let acronyms(rows) = block(
  stroke: (top: 1.2pt + navy, bottom: 1.2pt + navy),
  inset: 0pt,
  table(
    columns: (auto, 1fr),
    align: (right + horizon, left + horizon),
    inset: (x: 1em, y: 0.45em),
    stroke: none,
    table.hline(stroke: 0.5pt + navy),
    ..rows.map(((sym, desc)) => (sym, desc)).flatten(),
    table.hline(stroke: 0.5pt + navy),
  ),
)

#acronyms((
  ([*AAA*],     [Adaptive Antoulas--Anderson rational approximation. A barycentric rational-approximation algorithm that builds near-optimal poles automatically; the basis of the modern lightning solvers (@ch-coord-transforms, @ch-special-tricks).]),
  ([*BVP*],     [Boundary value problem. The continuous problem $cal(L) u = f$ on a bounded domain with prescribed conditions on $partial Omega$; the central object of @ch-bvp.]),
  ([*CFL*],     [Courant--Friedrichs--Lewy condition. The stability bound on the time step of an explicit method for a hyperbolic PDE; on Fourier grids it scales as $Delta t = cal(O)(N^(-1))$, on Chebyshev grids as $cal(O)(N^(-2))$ for waves and $cal(O)(N^(-4))$ for diffusion (@ch-time-stability).]),
  ([*DCT*],     [Discrete cosine transform. The cosine-only special case of the FFT used to compute Chebyshev coefficients via the standard "extended-vector" trick.]),
  ([*DFT*],     [Discrete Fourier transform. The mathematical operation underlying the FFT; for periodic data on $N$ equispaced points, the DFT is exactly the trapezoidal-rule approximation of the Fourier-series integral (@ch-periodic-trap).]),
  ([*DNS*],     [Direct numerical simulation. The fully resolved time-stepping computation of a turbulent flow with no subgrid model; the *2/3 dealiasing rule* of Étude 11.7 is its workhorse.]),
  ([*ETD*],     [Exponential time differencing. The family of time-stepping schemes that treats the linear part of a stiff PDE exactly via the matrix exponential and the nonlinear part by Runge--Kutta or multistep extrapolation.]),
  ([*ETDRK4*],  [Fourth-order ETD Runge--Kutta. The Cox--Matthews scheme, used for Kuramoto--Sivashinsky in Étude 11.6 and 14.9.]),
  ([*FD*],      [Finite difference. The local-stencil discretisation of derivatives on a grid; the algebraic-convergence baseline against which spectral methods are compared throughout (@ch-differentiation).]),
  ([*FE*],      [Finite element. The Galerkin-on-piecewise-polynomials family; mentioned mainly in contrast with spectral collocation in @ch-bvp and @ch-higher-order.]),
  ([*FFT*],     [Fast Fourier transform. The Cooley--Tukey $cal(O)(N log N)$ algorithm for the DFT; the workhorse of every Fourier and FFT-Chebyshev étude in @ch-fourier-grids -- @ch-fourier-pseudo.]),
  ([*FNO*],     [Fourier neural operator. A family of neural-network architectures that learn maps between function spaces by parameterising convolution in Fourier space.]),
  ([*IF*],      [Integrating factor. The time-stepping family that absorbs the linear part of a stiff ODE/PDE into a multiplicative exponential; used for KdV, NLS, and KS in @ch-fourier-pseudo and @ch-time-stability.]),
  ([*IMEX*],    [Implicit--explicit. A time-stepping family that splits the right-hand side into a stiff linear part (treated implicitly) and a non-stiff nonlinear part (treated explicitly); used for Allen--Cahn in Étude 11.10.]),
  ([*KdV*],     [Korteweg--de Vries equation, $u_t + 6 u u_x + u_(x x x) = 0$. The canonical integrable dispersive wave equation; Études 11.4 and 17.7.]),
  ([*KS*],      [Kuramoto--Sivashinsky equation, $u_t + u u_x + u_(x x) + u_(x x x x) = 0$. A fourth-order spatiotemporal-chaos model; Étude 11.6 and Étude 14.9.]),
  ([*KTE*],     [Kosloff--Tal-Ezer arcsine map. A coordinate transformation that flattens the Chebyshev grid toward equispacing to relax the CFL constraint; (@ch-coord-transforms §19.8).]),
  ([*KTL*],     [Kosloff--Tal-Ezer least-squares (De Marchi--Marchetti--Perracchione 2023). A modern refinement of the KTE map that decouples the polynomial degree from the number of nodes via a least-squares fit, restoring spectral accuracy with the aggressive map.]),
  ([*LARS*],    [Lightning--AAA Rational Stokes algorithm (Xue--Waters--Trefethen 2024). Pole-clustering rational approximation for two-dimensional Stokes flows in polygonal geometries.]),
  ([*MSC*],     [Mathematics Subject Classification. The American Mathematical Society's classification scheme; the four codes most relevant to this volume are 65N35, 65T50, 65L10, 65L60.]),
  ([*NeRF*],    [Neural radiance field. A neural-network parametrisation of view-dependent volumetric scene functions; the aliasing analysis underlying mip-NeRF (@ch-fourier-grids literature review) is a recent application of classical sampling theory.]),
  ([*NLS*],     [Nonlinear Schrödinger equation, $i u_t + u_(x x) + 2 |u|^2 u = 0$. The canonical envelope-soliton equation; Étude 11.5.]),
  ([*NS*],      [Navier--Stokes equations. The continuum equations of viscous fluid flow; the 2D vorticity formulation drives Étude 11.7.]),
  ([*NUFFT*],   [Non-uniform fast Fourier transform. An $cal(O)(N log N)$ algorithm for evaluating Fourier sums at arbitrary off-grid points; relevant whenever the discretisation grid does not coincide with the sample grid (@ch-fourier-grids literature review).]),
  ([*ODE*],     [Ordinary differential equation. Single independent variable.]),
  ([*PDE*],     [Partial differential equation. Multiple independent variables.]),
  ([*PINN*],    [Physics-informed neural network. A neural-network ansatz that includes the residual of a PDE in its loss function; the spectral-bias of standard activation functions is the connection to aliasing (@ch-fourier-grids literature review).]),
  ([*QNM*],     [Quasinormal mode. A complex-eigenvalue mode of a non-self-adjoint operator with outgoing-wave boundary conditions; central to Étude 13.6 and to the gravitational-wave-astronomy literature reviewed in the Afterword.]),
  ([*QR*],      [Orthogonal--triangular factorisation $bold(A) = bold(Q) bold(R)$, and the iterative QR algorithm built upon it. The standard dense-matrix eigensolver (@appendix-linalg).]),
  ([*QZ*],      [Generalised Schur algorithm for matrix pencils $bold(A) bold(u) = lambda bold(B) bold(u)$. The standard generalised-eigenproblem solver; absorbs the subtle handling of rank-deficient $bold(B)$ via the homogeneous $(alpha, beta)$ representation (@appendix-linalg, Étude 18.7).]),
  ([*RK4*],     [Classical fourth-order explicit Runge--Kutta method.]),
  ([*RKC*],     [Runge--Kutta--Chebyshev. A stabilised explicit scheme with a long stability interval along the negative real axis; used for diffusion-dominated problems where implicit solves are expensive.]),
  ([*sFFT*],    [Sparse fast Fourier transform. An $cal(O)(K log N)$ algorithm for the DFT of an $N$-vector with at most $K$ nonzero coefficients in frequency space.]),
  ([*SVD*],     [Singular value decomposition (@appendix-linalg). The factorisation underlying pseudospectra, condition numbers, and rank-revealing computations.]),
  ([*TLS*],     [Total least squares. A regression formulation that allows errors in *both* the design matrix and the right-hand side; the natural fit for spectral parameter estimation when both data and model are noisy.]),
  ([*WKB*],     [Wentzel--Kramers--Brillouin (also Liouville--Green) asymptotic approximation. The semi-classical leading-order behaviour of solutions of $u'' + Q(x) u = 0$; the structural ansatz behind the asymptotic-augmentation strategy of Étude 20.11.]),
  ([*WKS*],     [Whittaker--Shannon--Kotelnikov sampling theorem. The reconstruction-from-samples theorem for band-limited functions; the analytical foundation of the sinc interpolation of @ch-fourier-grids and the rational-Chebyshev sinc method of @ch-unbounded.]),
))

#v(0.6em)
A handful of secondary acronyms appear once or twice and are defined inline at first use: *FEAST* (Polizzi's contour-integral eigensolver), *KTE5* / *KTE-DDD* (specific KTE-map parameter conventions), *MIT* / *AMS* / *SIAM* / *PNAS* (publishers and learned societies), and a few language-specific identifiers (*BLAS*, *LAPACK*, *MKL*, *OpenBLAS*, *FFTW*) that any reader who has installed the dependencies will recognise.
