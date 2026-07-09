// textbook/chapters/polar_coordinates.typ
// Chapter 12: Spectral Methods in Polar Coordinates
// Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
// Email: denys.dutykh@ku.ac.ae
// Homepage: https://www.denys-dutykh.com/
// Last modified: February 2026

#import "../styles/template.typ": dropcap, num, format-table, etude-conclusion, idx, chapter-abstract, exercise, hint-for

// Enable equation numbering for this chapter

= Spectral Methods in Polar Coordinates <ch-polar>

#chapter-abstract(keywords: [Polar coordinates · Coordinate singularity · Doubling trick · Polar Laplacian · Bessel functions · Kronecker products])[
Many problems of scientific interest live not on rectangles but on circular domains: vibrating drumheads, optical fibres, heated wafers, quantum corrals, and pipe flow. This chapter extends the spectral framework to the unit disk by working in polar coordinates and combining Chebyshev discretisation in the radial direction with Fourier discretisation in the angular direction. The principal obstacle is the coordinate singularity at the origin, where factors of $1\/r$ and $1\/r^2$ appear in the polar Laplacian even though the solution itself is perfectly regular there. A naive Chebyshev mapping onto $[0, 1]$ both wastes resolution near the centre and invites severe ill-conditioning. The remedy is the doubling trick: the radius is extended to $[-1, 1]$ and the angular symmetry relating opposite radii is imposed, so that the grid never places a point at the origin at all. The discrete polar Laplacian is then assembled through Kronecker products of block-decomposed radial and angular operators. Computational études verify spectral accuracy by matching computed membrane eigenvalues against the zeros of Bessel functions, and solve the Poisson and heat equations on the disk.
]

#dropcap[The spectral machinery assembled in the preceding chapters --- Chebyshev differentiation for bounded, non-periodic intervals and Fourier methods for periodic domains --- has so far been deployed on rectangular geometries: intervals, squares, and periodic strips. But many problems of scientific interest live on _circular_ domains: vibrating drumheads, optical fibres with circular cross-section, semiconductor wafers under localised heating, quantum particles confined to circular corrals, and incompressible flow in pipes. In this chapter, we extend the spectral framework to the _unit disk_ by introducing polar coordinates#idx("polar coordinates") $(r, theta)$ and combining Chebyshev discretisation in the radial direction with Fourier discretisation in the angular direction. The principal new challenge is the _coordinate singularity#idx("coordinate singularity")_ at $r = 0$, which we overcome with an elegant algebraic device: the doubling trick#idx("doubling trick").]

By the end of this chapter, you should be able to:

1. Understand the coordinate singularity at $r = 0$ and why a naive Chebyshev mapping to $[0, 1]$ wastes resolution near the origin.
2. Apply the _doubling trick_ --- extending $r$ from $[0, 1]$ to $[-1, 1]$ with the symmetry condition $u(r, theta) = u(-r, theta + pi)$ --- to obtain a well-conditioned discretisation that avoids the origin entirely.
3. Assemble the discrete polar Laplacian#idx("polar Laplacian") via Kronecker products of block-decomposed Chebyshev matrices and Fourier angular operators.
4. Compute eigenmodes of the Laplacian on the disk and compare numerical eigenvalues with the zeros of Bessel function#idx("Bessel function")s to verify spectral accuracy.
5. Solve elliptic (Poisson) and parabolic (heat) PDEs on circular domains using the polar spectral framework.

== Why Polar? Motivation and Challenges <sec-polar-motivation>

=== Physical Domains with Circular Geometry

Circular and cylindrical geometries arise across science and engineering. The vibration of a circular drumhead#idx("circular drumhead") --- one of the oldest problems in mathematical physics --- requires solving the wave equation on a disk @MorseIngard1968. In fibre optics, the fundamental modes of a cylindrical waveguide are governed by the Helmholtz equation in a circular cross-section. Heat conduction in semiconductor wafers, which are manufactured as thin circular disks, naturally calls for solving the heat equation in polar coordinates. In quantum mechanics, the so-called _quantum corral_ --- a circular arrangement of atoms on a surface confining electron standing waves --- yields striking eigenmode patterns that can be directly compared with Bessel function predictions @Crommie1993. Incompressible pipe flow, the Couette--Taylor instability between rotating cylinders, vortex dynamics in geophysical flows, and the high-precision modelling of gravitational waves in numerical relativity all involve polar, cylindrical, or spherical domains.

All of these problems share a common feature: the geometry is _intrinsically circular_, and rectangular Chebyshev grids from @ch-spectral-pde are a poor fit. We need a coordinate system that respects the circular symmetry.

=== The Polar Laplacian

Polar coordinates on the unit disk are defined by the transformation
$ x = r cos theta, quad y = r sin theta, quad 0 lt.eq.slant r lt.eq.slant 1, quad 0 lt.eq.slant theta < 2 pi. $ <eq-polar-coords>
In these coordinates, the Laplacian takes the form
$ Delta u = u_(r r) + frac(1, r) u_r + frac(1, r^2) u_(theta theta). $ <eq-polar-laplacian>
The three terms have clear physical interpretations: $u_(r r)$ is the radial curvature, the factor $1\/r$ in the second term is a geometric focusing effect arising from the shrinking of angular arcs as $r arrow 0$, and $u_(theta theta) \/ r^2$ is the angular curvature scaled by the local arc length.

Equation @eq-polar-laplacian is the workhorse of this chapter: every PDE we solve (eigenmodes, Poisson, heat) reduces to discretising and inverting this operator.

=== The Coordinate Singularity at $r = 0$ <sec-pole-problem>

The polar Laplacian @eq-polar-laplacian contains the factors $1\/r$ and $1\/r^2$, which diverge as $r arrow 0$. This is the _pole problem_ (or _coordinate singularity_). It is important to emphasise that the singularity is _not_ physical: the solution $u(x, y)$ can be perfectly smooth at the origin. The singularity is purely an artefact of the polar coordinate representation.

To see this concretely, consider the function $u(x, y) = 1 - x^2 - y^2 = 1 - r^2$, which is smooth everywhere on the disk. Its polar Laplacian is
$ Delta u = -2 + frac(1, r)(-2 r) + frac(1, r^2) dot 0 = -4, $
a well-behaved constant. But computing the individual terms requires evaluating $u_r \/ r = -2 r \/ r = -2$, which at $r = 0$ presents a $0\/0$ indeterminate form. Historically, resolving this singularity required expanding the solution in Fourier series and explicitly enforcing _analytical pole conditions_: physical smoothness at the origin dictates that the $m$-th Fourier coefficient must decay as $O(r^(|m|))$ as $r arrow 0$ @Orszag1974. Any numerical method that places a grid point at $r = 0$ must somehow evaluate these singular-looking terms, and doing so stably is a nontrivial challenge.

=== The Naive Approach: Chebyshev on $[0, 1]$ <sec-naive-grid>

The most natural first attempt is to map the standard Chebyshev grid from $x in [-1, 1]$ to $r in [0, 1]$ via the affine transformation
$ r = frac(x + 1, 2). $ <eq-naive-map>
Combined with an equispaced Fourier grid in $theta$, this produces a valid polar grid on the disk. However, it has two serious drawbacks.

*Drawback 1: Wasted resolution.* The Chebyshev points cluster near both endpoints of $[-1, 1]$, which after the mapping @eq-naive-map means clustering near both $r = 0$ and $r = 1$. The clustering near the boundary $r = 1$ is desirable (it resolves boundary layers), but the clustering near $r = 0$ is wasteful: the small central region near the origin has relatively few degrees of freedom in the angular direction and contributes little to the overall solution accuracy. As noted by Trefethen @Trefethen2000, half the grid points lie inside a circle that encloses only about 25% of the total disk area.

*Drawback 2: Severe CFL restriction.* For time-dependent problems, the minimum spacing between adjacent grid points determines the maximum stable time step through the CFL condition. Near the origin, the Chebyshev clustering produces radial spacings of order $O(N_r^(-2))$, and the combined effect of radial and angular resolution gives a CFL constraint of $Delta t tilde.op O(N_r^(-4))$ for explicit wave-equation solvers --- far more restrictive than the $O(N_r^(-2))$ scaling on a uniform grid.

Both drawbacks motivate the search for a better approach.

=== Computational Étude 12.1: Grid Geometry and the Cost of Clustering <etude-polar-grids>

To make the comparison quantitative, we construct both the naive grid (Chebyshev on $r in [0, 1]$) and the doubled grid (Chebyshev on $r in [-1, 1]$, described in the next section) and compare their properties.

The grid construction is straightforward in all three languages.

#figure(
  image("../figures/ch12/python/polar_grids.pdf", width: 95%),
  caption: [Two spectral grids on the unit disk with $N_r = 25$ radial and $M = 20$ angular points. _Left_: the naive grid maps Chebyshev points to $r in [0, 1]$, resulting in heavy clustering near the origin. The dashed circle (containing half the grid points) encloses only 25% of the disk area. _Right_: the doubled grid uses Chebyshev points on $r in [-1, 1]$ with odd $N_r$, keeping only the positive-$r$ half. No point falls at the origin, and the resolution is more uniformly distributed. After the polar-grid illustration of Trefethen @Trefethen2000[Chap. 11].],
) <fig-polar-grids>

@fig-polar-grids reveals the dramatic difference in grid point distribution. The naive grid devotes an excessive number of points to the small neighbourhood of the origin, while the doubled grid distributes resolution more evenly. @fig-polar-area-cfl quantifies the consequences: the area per radial annulus is far more uniform for the doubled grid, and the maximum stable time step scales as $O(N_r^(-2))$ instead of $O(N_r^(-4))$.

#figure(
  image("../figures/ch12/python/polar_area_cfl.pdf", width: 95%),
  caption: [_Left_: area of each radial annulus for the naive and doubled grids. The naive grid allocates disproportionate area to the thin annuli near $r = 0$. _Right_: CFL time-step scaling. The naive grid forces $Delta t_(max) tilde.op N_r^(-4)$ (blue circles), while the doubled grid gives $Delta t_(max) tilde.op N_r^(-2)$ (red squares), a substantial improvement.],
) <fig-polar-area-cfl>

#etude-conclusion[
  The visual contrast is striking: the naive Chebyshev mapping to $[0, 1]$ concentrates roughly *half* its grid points inside a circle enclosing barely a third of the disk area. This is not aesthetic --- it has direct consequences for time-stepping. The CFL-limited time step scales as $cal(O)(N_r^(-4))$ for the naive grid versus the far milder $cal(O)(N_r^(-2))$ for the doubled grid: doubling the radial resolution on the naive grid forces a *sixteen-fold* reduction in $Delta t$, making explicit time integration prohibitively expensive even at moderate resolution. The doubled grid achieves its advantage by a simple algebraic device --- keeping $r in [-1, 1]$ and discarding the redundant half --- yet the savings are profound. The étude demonstrates a recurring theme: *the choice of coordinate mapping is not a mere technicality but can dominate the practical efficiency of a method*.
]

The code generating @fig-polar-grids and @fig-polar-area-cfl is available in:
- `codes/python/ch12/polar_grid_geometry.py`
- `codes/matlab/ch12/polar_grid_geometry.m`
- `codes/julia/ch12/polar_grid_geometry.jl`

== The Doubling Trick: Extending $r$ to $[-1, 1]$ <sec-doubling-trick>

The idea, due to Fornberg @Fornberg1995 @Fornberg1996, is beguilingly simple: instead of restricting $r$ to $[0, 1]$, we take
$ theta in [0, 2 pi), quad r in [-1, 1]. $ <eq-doubled-domain>
The Chebyshev grid now spans $[-1, 1]$ in $r$ without any affine rescaling, and its natural clustering occurs at $r = plus.minus 1$ (the boundary) rather than at $r = 0$.

=== The 2-to-1 Map and Single-Valuedness <sec-symmetry-condition>

What does it mean for $r$ to be negative? The polar coordinate mapping @eq-polar-coords extends naturally: the point $(r, theta)$ with $r < 0$ maps to the same Cartesian point as $(-r, theta + pi)$. In other words, "going backwards" in the radial direction from angle $theta$ is the same as going forwards from angle $theta + pi$. The map from $(r, theta)$ to $(x, y)$ is therefore _2-to-1_: every point in the disk (except the origin) has two preimages in the $(r, theta)$ rectangle.

For a function $u(r, theta)$ to represent a single-valued function of $(x, y)$, it must satisfy the _symmetry condition_
$ u(r, theta) = u(-r, (theta + pi) mod 2 pi). $ <eq-symmetry-condition>
This is the key constraint that makes the doubling trick work. It tells us that the function values in the $r < 0$ half of the $(r, theta)$ rectangle are completely determined by the values in the $r > 0$ half --- they are redundant copies, shifted by $pi$ in angle. The discrete system can therefore be reduced to half its apparent size.

=== Why $N_r$ Should Be Odd <sec-odd-nr>

To avoid ever evaluating the $1\/r$ and $1\/r^2$ terms at $r = 0$, we choose the number of Chebyshev points $N_r$ to be _odd_. The Chebyshev--Gauss--Lobatto points are $r_j = cos(j pi \/ N_r)$ for $j = 0, 1, dots, N_r$. When $N_r$ is odd, none of these points equals zero, since $cos(j pi \/ N_r) = 0$ requires $j = N_r \/ 2$, which is not an integer when $N_r$ is odd. The coordinate singularity is simply never encountered.

#block(
  stroke: (left: 2pt + rgb("#142D6E").lighten(60%)),
  inset: (left: 12pt, y: 6pt),
)[
  *Remark.* An alternative convention, used by some authors, places the radial nodes at the _interior_ Chebyshev (Gauss) points $r_j = cos((2 j - 1) pi \/ (2 N_r))$, $j = 1, dots, N_r$, which exclude the endpoints $r = plus.minus 1$ and, when $N_r$ is even, skip $r = 0$ as well. Both families work; we follow the odd-$N_r$ Gauss--Lobatto convention of Trefethen @Trefethen2000 and Fornberg @Fornberg1996 throughout this chapter. @ex-polar-even-odd-nr explores the consequences of using even $N_r$ on the Gauss--Lobatto grid, which then places a node exactly at the origin. It must be noted, however, that while the doubling trick bypasses physical pole evaluation, high-resolution simulations of nonlinear PDEs can still suffer from aliasing and stability issues near the origin unless explicit parity restrictions are enforced.
]

=== The Block Decomposition of Chebyshev Matrices <sec-block-decomposition>

The symmetry condition @eq-symmetry-condition turns a $2 N$-dimensional problem into an $N$-dimensional one. To implement this algebraically, we partition the Chebyshev differentiation matrices into blocks.

Consider the Chebyshev first-derivative matrix $D$ and second-derivative matrix $D^((2)) = D^2$ on the grid $r_0 = 1, r_1, dots, r_(N_r) = -1$. After removing the boundary rows and columns (corresponding to $r = 1$ and $r = -1$, where Dirichlet conditions are imposed), we have interior matrices of size $(N_r - 1) times (N_r - 1)$. The interior points split naturally into two groups:

- *Positive-$r$ points*: $r_1, r_2, dots, r_(N_2)$ where $N_2 = (N_r - 1) \/ 2$.
- *Negative-$r$ points*: $r_(N_2 + 1), dots, r_(N_r - 1)$.

The $j$-th positive point $r_j$ pairs with the point $-r_j = r_(N_r - j)$ under the symmetry condition. We partition the second-derivative matrix into four blocks:

$ tilde(D)_r^((2)) = mat(D_1, D_2; D_3, D_4), $ <eq-D2-blocks>

where $D_1$ ($N_2 times N_2$) couples positive-$r$ rows to positive-$r$ columns, $D_2$ couples positive-$r$ rows to negative-$r$ columns (with column indices _reversed_ to match the pairing $r_j arrow.l.r -r_j$), and $D_3$, $D_4$ form the lower blocks which we discard by the symmetry condition. Similarly, the first-derivative matrix is partitioned as

$ tilde(D)_r = mat(E_1, E_2; E_3, E_4). $ <eq-D1-blocks>

By the symmetry condition, the degrees of freedom in the $r < 0$ half are determined by those in the $r > 0$ half. When computing the Laplacian applied to a function satisfying @eq-symmetry-condition, only the first block row (positive-$r$ rows) is needed, but both column blocks contribute: $D_1$ and $E_1$ act on the solution at the _same_ angular positions, while $D_2$ and $E_2$ act on the solution at positions shifted by $pi$ in angle.

== Assembling the Discrete Polar Laplacian <sec-polar-laplacian>

=== Angular Discretisation: Fourier in $theta$ <sec-angular-fourier>

The angular variable $theta in [0, 2 pi)$ is periodic, so we discretise it with an equispaced Fourier grid of $M$ points (with $M$ even):
$ theta_m = frac(2 pi m, M), quad m = 0, 1, dots, M - 1. $ <eq-theta-grid>

For the second angular derivative $partial^2 \/ partial theta^2$, we use the $M times M$ Toeplitz matrix $D_theta^((2))$ whose entries are given by the spectral differentiation formula for periodic functions @Trefethen2000:
$ (D_theta^((2)))_(j j) = -frac(pi^2, 3 h^2) - frac(1, 6), quad (D_theta^((2)))_(j k) = frac((-1)^(j - k + 1), 2 sin^2((j - k) h \/ 2)), quad j eq.not k, $
where $h = 2 pi \/ M$ is the angular spacing. This matrix was encountered in @ch-differentiation in the context of periodic differentiation and can equivalently be computed via FFT multiplication by $-k^2$ in Fourier space, as in @ch-fourier-grids. In advanced Galerkin frameworks, global bases such as _Zernike polynomials_ @Zernike1934 are sometimes used in place of the decoupled Fourier--Chebyshev grids employed here: they are intrinsically orthogonal over the unit disk and inherently satisfy all polar regularity conditions at the origin @BoydYu2011.

=== Radial Operators and the $1\/r$ Weight <sec-radial-operators>

The discrete radial operators are built from the blocks $D_1$, $D_2$, $E_1$, $E_2$ extracted in @sec-block-decomposition, together with the diagonal weight matrix
$ R = op("diag")(1 \/ r_j), quad j = 1, dots, N_2, $ <eq-R-matrix>
which accounts for the $1\/r$ and $1\/r^2$ factors in the polar Laplacian @eq-polar-laplacian.

The effective radial second-derivative operator on the positive-$r$ half is $(D_1 + R E_1)$, which acts on the solution at the same angular positions. The coupling to the $pi$-shifted angular positions comes through $(D_2 + R E_2)$.

=== The Kronecker Product Assembly <sec-kronecker-polar>

Combining radial and angular discretisations via Kronecker products, the discrete polar Laplacian on the positive-$r$ interior grid takes the form
$ L = (D_1 + R E_1) times.o I_M + (D_2 + R E_2) times.o S + R^2 times.o D_theta^((2)), $ <eq-polar-laplacian-discrete>
where $I_M$ is the $M times M$ identity, $S$ is the _block swap matrix#idx("swap matrix")_
$ S = mat(0, I_(M\/2); I_(M\/2), 0) $ <eq-swap-matrix>
that implements the angular $pi$-shift from the symmetry condition, and $R^2 = op("diag")(1 \/ r_j^2)$.

The matrix $L$ has dimensions $(N_2 M) times (N_2 M)$. For typical parameters ($N_r = 25$, $M = 20$, giving $N_2 = 12$), this is a $240 times 240$ system --- modest enough for dense linear algebra, yet large enough to achieve spectral accuracy.

The structure of @eq-polar-laplacian-discrete closely parallels the tensor-product Laplacian on a rectangle from @ch-spectral-pde: the first two Kronecker products handle radial differentiation (with the symmetry folding), and the third handles angular differentiation. The only novelty is the swap matrix $S$, which encodes the 2-to-1 map.

=== Implementation: The `laplacian_polar` Function <sec-laplacian-impl>

We encapsulate the assembly into a reusable function.

Note the critical detail in the block extraction: the negative-$r$ columns are indexed in _reverse order_ (`Nr-1:N2:-1` in 0-based Python, `Nr:-1:N2+2` in 1-based MATLAB/Julia). This reversal ensures that the $j$-th positive-$r$ point is correctly paired with its symmetric counterpart $-r_j$.

The `laplacian_polar` function is available in:
- `codes/python/ch12/laplacian_polar.py`
- `codes/matlab/ch12/laplacian_polar.m`
- `codes/julia/ch12/laplacian_polar.jl`

== Eigenmodes of the Circular Membrane <sec-disk-eigenmodes>

With the polar Laplacian assembled, we turn to the classical eigenproblem: finding the normal modes of vibration of a circular drumhead. This is perhaps the most beautiful application of spectral methods on the disk, and it connects directly to the theory of Bessel functions.

=== Analytical Background: Bessel Functions and Separation of Variables <sec-bessel-background>

The eigenvalue problem on the unit disk with homogeneous Dirichlet boundary conditions reads
$ -Delta u = lambda u, quad u|_(r = 1) = 0. $ <eq-disk-eigenproblem>
Separation of variables $u(r, theta) = R(r) Theta(theta)$ leads to the angular equation $Theta'' + m^2 Theta = 0$ (solved by $cos m theta$ and $sin m theta$ for integer $m gt.eq.slant 0$) and the radial Bessel equation
$ r^2 R'' + r R' + (k^2 r^2 - m^2) R = 0, $
where $lambda = k^2$. The bounded solution is the Bessel function of the first kind, $R(r) = J_m (k r)$. The Dirichlet condition $R(1) = 0$ requires $J_m (k) = 0$, so the eigenvalues are
$ lambda_(m, n) = j_(m, n)^2, $ <eq-bessel-eigenvalues>
where $j_(m, n)$ denotes the $n$-th positive zero of $J_m$. For $m = 0$, the eigenmodes are radially symmetric; for $m gt.eq.slant 1$, each eigenvalue has _multiplicity two_ (corresponding to the $cos m theta$ and $sin m theta$ modes), reflecting the rotational symmetry of the disk.

The first few Bessel zeros#idx("Bessel zeros") and eigenvalues are:
- $j_(0, 1) approx 2.4048$, so $lambda_(0, 1) approx 5.7832$ (fundamental mode, no nodal line#idx("nodal line")s inside the disk).
- $j_(1, 1) approx 3.8317$, so $lambda_(1, 1) approx 14.6820$ (first non-radial mode, one nodal diameter).
- $j_(2, 1) approx 5.1356$, so $lambda_(2, 1) approx 26.3746$ (two nodal diameters).
- $j_(0, 2) approx 5.5201$, so $lambda_(0, 2) approx 30.4713$ (one nodal circle, no nodal diameters).

The spectral method does _not_ use Bessel functions: it discretises the Laplacian as a matrix and computes eigenvalues by dense linear algebra. The fact that the numerical eigenvalues match the Bessel zeros to full floating-point precision is a striking confirmation of spectral accuracy.

=== Computational Étude 12.2: Vibrations of a Circular Drumhead <etude-disk-eigenmodes>

We compute the first 25 eigenpairs of the Laplacian on the unit disk with $N_r = 25$ and $M = 20$, then compare the numerical eigenvalues with the exact values @eq-bessel-eigenvalues.

#figure(
  image("../figures/ch12/python/disk_eigenmodes.pdf", width: 95%),
  caption: [Six eigenmodes of the Laplacian on the unit disk, computed with $N_r = 25$ radial and $M = 20$ angular points. The modes display the characteristic Bessel function patterns: the fundamental mode ($m = 0$, $n = 1$) is a simple dome, while higher modes exhibit nodal diameters (radial lines) and nodal circles. The mode selection follows Program 28 of Trefethen @Trefethen2000.],
) <fig-disk-eigenmodes>

#figure(
  image("../figures/ch12/python/eigenvalue_convergence.pdf", width: 95%),
  caption: [_Left_: the first 25 computed eigenvalues (circles) agree with the exact Bessel zeros squared (crosses) to plotting accuracy. The degenerate pairs ($m gt.eq.slant 1$) appear as overlapping markers. _Right_: the error in the fundamental eigenvalue $lambda_(0,1)$ decreases exponentially with $N_r$, confirming spectral convergence.],
) <fig-eigenvalue-convergence>

The code generating @fig-disk-eigenmodes and @fig-eigenvalue-convergence is available in:
- `codes/python/ch12/disk_eigenmodes.py`
- `codes/matlab/ch12/disk_eigenmodes.m`
- `codes/julia/ch12/disk_eigenmodes.jl`

#etude-conclusion[
  The figure displays the rich modal structure of the circular membrane --- from the smooth axisymmetric dome of the fundamental mode to the intricate nodal patterns of higher modes (radial node lines, concentric nodal circles). These patterns, analytically predicted by Bessel-function theory since the 19th century, are recovered here *without any explicit reference to Bessel functions*: spectral discretisation of the Laplacian plus a call to a standard eigensolver suffice. The convergence is exponential, with agreement to 10+ digits against the squared Bessel zeros $j_(m, n)^2$ --- a powerful cross-validation that both the doubling trick and the Kronecker assembly are implemented correctly. *Watch for degenerate pairs*: for each angular wavenumber $m gt.eq.slant 1$, the eigenvalue $lambda_(m, n)$ appears twice ($cos m theta$ and $sin m theta$ modes), and the *individual* eigenvectors returned by the eigensolver are arbitrary rotations within the 2D eigenspace. This degeneracy is the topic of the next étude.
]

=== Computational Étude 12.3: Eigenvalue Degeneracy and Nodal Rotations <etude-degenerate-modes>

For each $m gt.eq.slant 1$, the eigenvalue $lambda_(m, n) = j_(m, n)^2$ has multiplicity two: the two eigenmodes $J_m (j_(m, n) r) cos m theta$ and $J_m (j_(m, n) r) sin m theta$ span a two-dimensional eigenspace. Any linear combination
$ u_phi = cos phi dot u_1 + sin phi dot u_2 $ <eq-rotation-combination>
is also an eigenmode for the same eigenvalue. As the parameter $phi$ varies from $0$ to $pi$, the nodal lines of $u_phi$ rotate continuously through the disk, tracing out the full family of modes in the eigenspace.

This étude visualises this phenomenon for the first degenerate pair ($m = 1$, $n = 1$, $lambda approx 14.682$). The two numerically computed eigenvectors are combined according to @eq-rotation-combination for six equally spaced values of $phi$, and the zero contours (nodal lines) are plotted on the disk.

#figure(
  image("../figures/ch12/python/nodal_rotation.pdf", width: 95%),
  caption: [Nodal line rotation in the first degenerate eigenspace ($lambda approx 14.68$). As the mixing angle $phi$ varies, the single nodal diameter rotates continuously through the disk. The faint background shading shows the sign of the eigenmode (red positive, blue negative). This continuous family of modes reflects the underlying rotational symmetry of the circular domain.],
) <fig-nodal-rotation>

The code generating @fig-nodal-rotation is available in:
- `codes/python/ch12/disk_nodal_rotation.py`
- `codes/matlab/ch12/disk_nodal_rotation.m`
- `codes/julia/ch12/disk_nodal_rotation.jl`

#etude-conclusion[
  The figure reveals a phenomenon invisible in eigenvalue tables but physically fundamental: *within a degenerate eigenspace, the nodal-line pattern is not fixed but can rotate freely*. As the mixing angle $phi$ sweeps from $0$ to $pi$, the single nodal diameter glides smoothly around the disk, and every orientation is an equally valid eigenmode with the same frequency. This is a direct consequence of the circular symmetry of the domain. From a computational perspective, the étude highlights an important subtlety of numerical eigensolvers: *when an eigenvalue has multiplicity > 1, the returned eigenvectors are not uniquely determined* --- different runs, algorithms, or even floating-point orderings may produce different rotations within the eigenspace. This ambiguity is *not a defect* but a faithful reflection of the underlying physics. Any experiment that breaks circular symmetry (boundary deformation, non-uniform density, localised perturbation) would lift the degeneracy and select a preferred orientation.
]

== The Poisson Equation on the Disk <sec-disk-poisson>

=== Problem Formulation

Having validated the Laplacian through the eigenproblem, we now solve the Poisson equation on the disk:
$ -Delta u = f(r, theta), quad 0 lt.eq.slant r lt.eq.slant 1, quad u|_(r = 1) = 0. $ <eq-poisson-disk>
This models, for instance, the static deflection of a clamped circular membrane under a distributed load $f$, or the electrostatic potential in a grounded cylindrical conductor with a charge distribution.

We choose a localised off-centre Gaussian source that simulates pressing on a drumhead with a finger at position $(x_0, y_0) = (0.4, 0.2)$:
$ f(r, theta) = exp(-beta [(r cos theta - x_0)^2 + (r sin theta - y_0)^2]), $ <eq-gaussian-source>
with $beta = 40$. This forcing is physically meaningful and breaks the radial symmetry, producing an asymmetric deflection pattern.

=== Computational Étude 12.4: An Off-Centre Load on a Membrane <etude-disk-poisson>

With the Laplacian matrix $L$ already assembled, the Poisson equation @eq-poisson-disk reduces to a linear system: evaluate $f$ at the interior grid points, flatten into a vector $bold(f)$, and solve
$ L bold(u) = -bold(f). $
Here $L$ is the discrete Laplacian (so $L bold(u) approx Delta u$), and the minus sign accounts for the sign convention in @eq-poisson-disk.

#figure(
  image("../figures/ch12/python/poisson_disk_surface.pdf", width: 80%),
  caption: [Three-dimensional surface plot of the Poisson solution on the unit disk with the off-centre Gaussian load @eq-gaussian-source. The membrane deflects most strongly near the load position $(0.4, 0.2)$, with the deflection smoothly decaying to zero at the clamped boundary $r = 1$.],
) <fig-poisson-disk-surface>

#figure(
  image("../figures/ch12/python/poisson_disk_contour.pdf", width: 65%),
  caption: [Filled contour plot of the Poisson solution. The red star marks the centre of the Gaussian load. The equipotential lines (level curves) are distorted toward the load location, reflecting the broken radial symmetry.],
) <fig-poisson-disk-contour>

As a final diagnostic, we solve the same problem with a _centred_ Gaussian source ($x_0 = y_0 = 0$) and compare radial profiles $u(r, theta_k)$ at eight equally spaced angles $theta_k$. For the centred source, exact radial symmetry requires all profiles to coincide identically; any numerical deviation would betray artificial anisotropy in the discretisation. For the off-centre source, the radial profiles serve as a quantitative fingerprint of the source location: the angle along which the profile peaks most strongly identifies the source direction, and the radial position of that peak estimates the source's distance from the origin.

#figure(
  image("../figures/ch12/python/radial_symmetry_test.pdf", width: 95%),
  caption: [Radial profiles $u(r, theta_k)$ at eight equally spaced angles $theta_k$, with $r = 1$ (clamped boundary) on the left and $r arrow 0$ (origin) on the right. All profiles vanish at $r = 1$ by the homogeneous Dirichlet condition. _Left_: for the centred source ($x_0 = y_0 = 0$), all eight profiles are indistinguishable to machine precision, collapsing onto a single monotone curve that peaks near the origin. This exact collapse — not merely visual agreement — confirms that the spectral discretisation introduces no artificial anisotropy. _Right_: for the off-centre source at $(0.4, 0.2)$, the profiles fan out. The source lies at radial distance $r_s = sqrt(0.4^2 + 0.2^2) approx 0.45$ from the origin and at angle $theta_s = arctan(0.2 \/ 0.4) approx 0.46$ rad. The orange profile ($theta approx 0.63$ rad, the sampled angle nearest to the source direction) peaks near $r approx 0.45$, consistent with the source distance. Profiles at angles pointing away from the source remain substantially lower throughout, and all profiles reconverge to zero as $r arrow 1$.],
) <fig-radial-symmetry-test>

The code generating @fig-poisson-disk-surface, @fig-poisson-disk-contour, and @fig-radial-symmetry-test is available in:
- `codes/python/ch12/disk_poisson.py`
- `codes/matlab/ch12/disk_poisson.m`
- `codes/julia/ch12/disk_poisson.jl`

#etude-conclusion[
  Poisson with an off-centre Gaussian source is a *stringent test* of the polar spectral framework beyond eigenproblems, because the off-centre source breaks every symmetry the disk possesses. The spectral method resolves the localised peak and the smooth decay to the clamped boundary without spurious oscillations. The radial-symmetry test is an especially valuable *diagnostic*: when the source is centred, all radial profiles at different angles collapse onto a single curve to machine precision --- confirming that the discretisation introduces no artificial anisotropy and that the swap matrix $S$ + Kronecker assembly are correct. A method that mishandles the coordinate singularity would produce angle-dependent artefacts even for a radially symmetric problem. For the off-centre source, the *fan-shaped* family of radial profiles encodes the same geometric information as the contour plot, in a form more amenable to quantitative reading.
]

== Time-Dependent PDEs on the Disk <sec-disk-time-dependent>

The polar Laplacian $L$ can be used as a building block for time-dependent PDEs, just as the rectangular Chebyshev and Fourier Laplacians were used in @ch-spectral-pde and @ch-fourier-pseudo. The method of lines applies directly: discretise space with the polar spectral method, then advance in time with a suitable ODE integrator.

=== Heat Equation: Cooling of a Semiconductor Wafer <sec-disk-heat>

We consider the heat equation on the unit disk with a steady localised source,
$ u_t = alpha Delta u + S(r, theta), quad u|_(r = 1) = 0, quad u(r, theta, 0) = 0, $ <eq-heat-disk>
modelling a circular semiconductor wafer heated by a localised laser spot and cooled through its boundary (held at ambient temperature). Here $alpha = 0.1$ is the thermal diffusivity and
$ S(r, theta) = 2 exp(-30 [(r cos theta - 0.3)^2 + (r sin theta)^2]) $ <eq-heat-source>
is a Gaussian heat source centred at $(0.3, 0)$.

After spatial discretisation, @eq-heat-disk becomes the ODE system
$ dot(bold(u)) = alpha L bold(u) + bold(S), $
which we integrate with the _Crank--Nicolson_ scheme, following the implicit time-stepping framework developed in @ch-spectral-pde:
$ (I - frac(alpha Delta t, 2) L) bold(u)^(n+1) = (I + frac(alpha Delta t, 2) L) bold(u)^n + Delta t bold(S). $ <eq-crank-nicolson-polar>
The left-hand side matrix is constant (independent of $n$), so we factorise it once using LU decomposition and reuse the factors at every time step.

=== Computational Étude 12.5: Heat Dissipation with a Localised Source <etude-disk-heat>

Starting from a cold initial condition $u = 0$, the laser source @eq-heat-source heats the wafer. The temperature rises rapidly near the source, then diffuses outward and eventually reaches a steady state where heat input balances boundary cooling.

#figure(
  image("../figures/ch12/python/heat_disk_snapshots.pdf", width: 95%),
  caption: [Heat dissipation on a circular wafer with a localised laser source @eq-heat-source. Snapshots at six times show the temperature evolving from a cold start ($t = 0$) to steady state ($t = 20$). The hot spot develops near the source location $(0.3, 0)$ and gradually diffuses outward, constrained by the zero-temperature boundary.],
) <fig-heat-disk-snapshots>

#figure(
  image("../figures/ch12/python/heat_disk_energy.pdf", width: 95%),
  caption: [_Left_: maximum temperature as a function of time, showing the approach to steady state. _Right_: comparison of the Crank--Nicolson steady state (solid) with the direct Poisson solve $alpha Delta u = -S$ (dashed) at $theta = 0$. The two curves agree, confirming that the time-dependent solution has converged to the correct equilibrium.],
) <fig-heat-disk-energy>

The code generating @fig-heat-disk-snapshots and @fig-heat-disk-energy is available in:
- `codes/python/ch12/disk_heat.py`
- `codes/matlab/ch12/disk_heat.m`
- `codes/julia/ch12/disk_heat.jl`

#etude-conclusion[
  This étude ties together the spatial machinery of the chapter with the time-stepping techniques of @ch-spectral-pde: the polar Laplacian $L$ slots *seamlessly* into the method-of-lines framework. The snapshots tell a physically compelling story --- a localised laser source rapidly heats a small region near $(0.3, 0)$, thermal diffusion spreads the heat outward, and the zero-temperature boundary condition acts as a heat sink that ultimately balances the input. The most informative validation is the cross-check between the *time-marched* steady state and a *direct* solve of the steady Poisson equation $alpha Delta u = -S$: the two computations exercise very different code paths (LU factorisation + repeated back-substitution vs single linear-system inversion), and their agreement confirms both the temporal integrator and the spatial discretisation. *Crank--Nicolson* is deliberate here: unconditional stability for diffusion problems means $Delta t$ is limited only by accuracy, not by the stiffness of the discrete Laplacian.
]

== Extensions and Generalisations <sec-polar-extensions>

=== Annular Domains

The methods of this chapter adapt readily to the _annulus_ $r_("in") lt.eq.slant r lt.eq.slant r_("out")$, which models, for example, the space between two concentric pipes or the cross-section of a fibre-optic cladding. The Chebyshev grid on $x in [-1, 1]$ is mapped to $[r_("in"), r_("out")]$ by the affine transformation
$ r = frac(r_("out") + r_("in"), 2) + frac(r_("out") - r_("in"), 2) x. $
Since $r_("in") > 0$, the coordinate singularity at $r = 0$ is absent, and the doubling trick is not needed. The Kronecker product structure of the Laplacian carries over unchanged, with $R = op("diag")(1\/r_j)$ evaluated at the mapped grid points. Dirichlet (or Neumann, or Robin) conditions are imposed independently at both boundaries. @ex-polar-annulus-eigs explores this setting.

=== Neumann and Robin Boundary Conditions

Throughout this chapter, we have imposed Dirichlet conditions $u = 0$ at $r = 1$. Other boundary conditions require modifying the radial differentiation matrices. For the Neumann condition $partial u \/ partial r |_(r=1) = 0$, one replaces the boundary rows of the Chebyshev matrix with the entries of the first-derivative matrix $D$, following the _boundary bordering_ technique from @ch-bvp. Robin conditions $alpha u + beta partial u \/ partial r = g$ at $r = 1$ are handled similarly.

=== Cylindrical and Spherical Geometries

Adding a third variable $z in [-H, H]$ to the polar $(r, theta)$ system gives _cylindrical coordinates_ $(r, theta, z)$. If $z$ is bounded and non-periodic, it gets a Chebyshev grid; if periodic, a Fourier grid. The Laplacian becomes a three-way Kronecker product. Recent advancements have optimised this setting: Darrow @Darrow2023 developed a quasi-optimal Chebyshev--Chebyshev--Fourier (CCF) solver that applies the radial doubling trick together with Alternating Direction Implicit (ADI) methods, reducing the computational complexity of 3D cylinder solves from $O(N^(4\/3))$ to $O(N log N)$.

For _spherical coordinates_ $(r, theta, phi)$, the angular part requires _spherical harmonics_ rather than simple Fourier modes. An alternative is to treat the polar angle $phi in [0, pi]$ with a Chebyshev (or Gauss--Legendre) grid and the longitude $theta in [0, 2 pi)$ with Fourier, using parity arguments at the poles analogous to our treatment of $r = 0$.

The Galerkin approach using _Zernike polynomials_ or _Jacobi polynomials_ (specifically, the so-called "one-sided Jacobi" bases) provides another route: basis functions can be chosen to satisfy the pole conditions automatically, eliminating the need for the doubling trick entirely @Shen1995 @Shen1997 @ShenTangWang2011. Furthermore, for complex tensorial systems like the Navier--Stokes equations, _spin-weighted Jacobi polynomials_ act as exact eigenfunctions for covariant derivatives, preventing the catastrophic mixing of vector components at the pole that arises from Christoffel symbols in curvilinear coordinates @Vasil2016. For highly complex domains with piecewise-smooth boundaries, polar spectral-element methods utilising _Log Orthogonal Functions_ are now deployed to capture corner singularities while maintaining spectral accuracy @ChenSunWu2025.

== A non-exhaustive literature overview

The application of spectral methods to polar coordinates and circular domains has a rich, complex, and evolving history, characterised primarily by a continuous mathematical struggle against the coordinate singularity at the origin. The foundational formulation presented in this chapter --- utilising a doubled radial variable $r in [-1, 1]$ paired with a geometric symmetry constraint --- was systematically developed, rigorously analysed, and popularised by Fornberg @Fornberg1995 @Fornberg1996 @Fornberg1996. Trefethen @Trefethen2000 distilled the doubling trick into a highly compact, elegant algebraic implementation (Programs 28 and 29 in his book), which forms the core architecture of the computational études presented here. Historically, the challenge of the pole problem was treated by explicitly imposing analytical regularity conditions on the spectral expansion. As identified in the seminal works of Orszag @Orszag1974 in the context of spherical atmospheric modelling, physical smoothness at the pole dictates that the $m$-th Fourier coefficient must decay as $O(r^(|m|))$. Boyd @Boyd2000 provides a comprehensive survey of these early methodologies, cataloguing dozens of approaches ranging from intricate coordinate transformations to the explicit, and often ill-conditioned, enforcement of parity boundary constraints at the origin. Matsushima and Marcus @Matsushima1995 provided one of the first highly robust spectral methods specifically tailored for vortex dynamics in polar coordinates, demonstrating that proper polynomial selection prevents the severe $O(N_r^(-4))$ CFL restrictions typically associated with radial point clustering.

An important analytical alternative to the collocation approach is the Galerkin method, which utilises global basis functions constructed to satisfy the pole conditions inherently. The most natural basis for the unit disk is the _Zernike polynomial_ family. Originally formulated to quantify optical wavefront aberrations @Zernike1934, Zernike polynomials naturally embed the requisite $O(r^(|m|))$ pole conditions. A thorough evaluation of spectral bases on the disk, including Zernike, Logan--Shepp ridge polynomials, and shifted Chebyshev series, was conducted by Boyd and Yu @BoydYu2011, highlighting the persistent trade-offs between mathematical elegance and computational matrix conditioning. To resolve the issue of dense differentiation matrices, Shen @Shen1997 introduced the use of "one-sided Jacobi" polynomials to construct optimally sparse Galerkin operators, drastically reducing the computational cost of direct solvers in cylindrical domains; the comprehensive framework of Shen, Tang, and Wang @ShenTangWang2011 extended these ideas further, yielding strictly banded operator matrices. These Galerkin methods are preferred in large-scale three-dimensional computations where the $O(N^3)$ cost of dense linear algebra becomes prohibitive.

The modern era of spectral methods (post-2010) has been defined by the pursuit of strictly sparse, well-conditioned operators that scale quasi-optimally in high dimensions. Olver and Townsend @OlverTownsend2013 reshaped non-periodic spectral linear algebra by introducing the _ultraspherical spectral method_: by shifting the polynomial basis after differentiation --- recognising that the derivative of a Chebyshev polynomial is an ultraspherical polynomial --- they generated strictly banded operators, reducing the dense $O(N^3)$ complexity of traditional Chebyshev boundary value problems to a highly efficient $O(N log N)$. This philosophy of banded sparsity was masterfully extended to polar and spherical geometries by Vasil, Burns, Lecoanet, Olver, Brown, and Oishi @Vasil2016, who utilised spin-weighted Jacobi polynomials to perform exact tensor calculus at the pole, avoiding the catastrophic mixing of vector components induced by Christoffel symbols. Their methodology forms the mathematical backbone of the widely used _Dedalus_ simulation framework. Recently, Darrow @Darrow2023 applied these sparse techniques to the three-dimensional closed cylinder, yielding a quasi-optimal Chebyshev--Chebyshev--Fourier solver that manages the radial doubling trick efficiently via Alternating Direction Implicit (ADI) methods.

Finally, the cutting edge of research continues to resolve specific structural and geometric challenges. To ensure stability and unitarity in equations like the linear Schrödinger equation, Gao and Iserles @GaoIserles2025 introduced a novel framework using W-functions and affine space splitting to generate skew-symmetric differentiation matrices in $d$-dimensional unit balls, preventing the artificial numerical dissipation of energy. Simultaneously, for engineering applications involving highly irregular boundaries, Chen, Sun, and Wu @ChenSunWu2025 advanced the spectral-element method in polar coordinates, utilising _Log Orthogonal Functions_ (LOFs) to maintain exponential spectral accuracy even in the presence of corner singularities and piecewise-smooth geometries. The analytical eigenmodes of the disk involve Bessel functions, whose theory was developed in the monumental treatise of Watson @Watson1922 and tabulated by Abramowitz and Stegun @Abramowitz1972. The physical problem of vibrating circular membranes was studied systematically by Morse and Ingard @MorseIngard1968. Applications of spectral methods on disks and in polar coordinates span a wide range of fields: quantum corrals and electron standing waves on metal surfaces @Crommie1993, vortex dynamics in geophysical flows, acoustic scattering from circular obstacles, and optical mode computation in fibre waveguides. Together, these advancements highlight that while the doubling trick remains an elegant pedagogical tool, the broader landscape of polar spectral methods now relies on sophisticated polynomial mapping, sparse linear algebra, and strict structure-preserving operators.

== Summary <sec-polar-summary>

This chapter has developed spectral methods for partial differential equations on the unit disk, using polar coordinates with a Chebyshev--Fourier discretisation. The main themes are:

+ *The coordinate singularity at $r = 0$* is an artefact of the polar representation, not a physical feature of the solution. It manifests as $1\/r$ and $1\/r^2$ factors in the Laplacian that diverge at the origin.

+ *The doubling trick* --- extending $r$ from $[0, 1]$ to $[-1, 1]$ and imposing the symmetry condition $u(r, theta) = u(-r, theta + pi)$ --- sidesteps the singularity entirely. With odd $N_r$, no grid point falls at $r = 0$, and the $1\/r$ division is never evaluated.

+ *Block decomposition* of the Chebyshev matrices implements the symmetry condition algebraically, halving the effective problem size. The reversed column indexing for the negative-$r$ block pairs each positive-$r$ point with its symmetric counterpart.

+ *The discrete polar Laplacian* is assembled via Kronecker products @eq-polar-laplacian-discrete, with structure closely paralleling the rectangular case from @ch-spectral-pde. The only novelty is the block swap matrix $S$ that encodes the angular $pi$-shift.

+ *With this machinery*, eigenmodes (matching Bessel zeros to 10 digits), Poisson problems, and heat equations on the disk are solved in compact code --- each étude requiring fewer than 50 lines per language. For strictly time-dependent PDEs requiring the preservation of unitarity (_e.g._, the linear Schrödinger equation), standard collocation differentiation matrices are non-normal and can introduce artificial dissipation; skew-symmetric differentiation matrices @GaoIserles2025 offer a structure-preserving alternative for such applications.

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
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Method*],
      ),
      table.hline(stroke: 0.5pt + luma(180)),
      [12.1], [Grid geometry on the disk], [Naive vs doubled grids], [Chebyshev mapping],
      [12.2], [$-Delta u = lambda u$ (Helmholtz)], [Polar Laplacian eigenproblem], [Dense `eig`],
      [12.3], [Degenerate eigenspace], [Nodal line rotation], [Linear combination],
      [12.4], [$-Delta u = f$ (Poisson)], [Off-centre Gaussian source], [Direct solve],
      [12.5], [$u_t = alpha Delta u + S$ (heat)], [Localised laser on wafer], [Crank--Nicolson],
    ),
  ),
  caption: [Summary of computational études in this chapter.],
) <tbl-ch12-summary>

== Exercises <sec-polar-exercises>

The exercises below progress from pencil-and-paper properties of the polar Laplacian and the doubling trick, through numerical experiments that reproduce and extend the études of this chapter, to open-ended projects that reach into the current research literature. The computational problems may be carried out in any of the book's three languages; the named scripts under `codes/` give a starting point.

=== Conceptual Exercises

#exercise(title: [The Pole Singularity Is Not Physical])[
  The polar Laplacian @eq-polar-laplacian carries the factors $1\/r$ and $1\/r^2$, yet a function that is smooth across the origin in Cartesian coordinates has a perfectly finite Laplacian there. (a) For $u = 1 - r^2$ evaluate each term of @eq-polar-laplacian separately, exhibit the indeterminate $0\/0$ form in $u_r \/ r$ at $r = 0$, and confirm that the assembled result is the constant $Delta u = -4$. (b) For a general smooth radial profile $u = g(r)$ with $g$ even and $g' (0) = 0$, use a Taylor expansion to show that $u_r \/ r arrow.r g'' (0)$ as $r arrow.r 0$, so the apparent singularity cancels. (c) State the pole regularity condition obeyed by the $m$-th angular Fourier coefficient, $u_m (r) tilde.op r^(|m|)$ as $r arrow.r 0$, and explain why it forces the potentially singular $u_(theta theta) \/ r^2$ term either to vanish (when $m = 0$), to stay bounded (when $|m| gt.eq.slant 2$), or to cancel against the equally singular $u_r \/ r$ term (when $|m| = 1$), so that the assembled Laplacian remains finite.
] <ex-polar-pole-not-physical>

#exercise(title: [Deriving the Polar Laplacian])[
  Derive @eq-polar-laplacian from the Cartesian Laplacian $Delta u = u_(x x) + u_(y y)$ using the change of variables @eq-polar-coords. (a) From $r = sqrt(x^2 + y^2)$ and $theta = arctan(y \/ x)$, compute $r_x, r_y, theta_x, theta_y$ and express $partial_x$ and $partial_y$ in terms of $partial_r$ and $partial_theta$. (b) Apply the first-order operators twice and collect terms to obtain $Delta u = u_(r r) + r^(-1) u_r + r^(-2) u_(theta theta)$. (c) Identify which step produces the $1\/r$ term and explain geometrically why it represents the spreading of radial flux as $r$ increases.
] <ex-polar-laplacian-derivation>

#exercise(title: [Symmetry Forces Angular Parity])[
  Expand a single-valued function on the doubled disk in an angular Fourier series, $u(r, theta) = sum_m u_m (r) thin e^(i m theta)$ with $r in [-1, 1]$. (a) Substitute into the symmetry condition @eq-symmetry-condition and use $e^(i m (theta + pi)) = (-1)^m e^(i m theta)$ to show that each radial coefficient must satisfy the parity relation $u_m (-r) = (-1)^m u_m (r)$. (b) Explain how this even/odd parity in $r$ is exactly the bookkeeping performed by the block swap matrix $S$ of @eq-swap-matrix when it pairs $r_j$ with $-r_j$. (c) Connect the parity $u_m (-r) = (-1)^m u_m (r)$ to the origin behaviour $u_m (r) tilde.op r^(|m|)$ required for smoothness.
] <ex-polar-parity>

#hint-for(<ex-polar-parity>)[Substitute the Fourier series into @eq-symmetry-condition and use $e^(i m (theta + pi)) = (-1)^m e^(i m theta)$; matching the coefficient of $e^(i m theta)$ on each side isolates $u_m (-r) = (-1)^m u_m (r)$.]

#exercise(title: [Why Odd $N_r$ Misses the Origin])[
  The Chebyshev--Gauss--Lobatto nodes are $r_j = cos(j pi \/ N_r)$ for $j = 0, 1, dots, N_r$, a set of $N_r + 1$ points on $[-1, 1]$. (a) Show that $r_j = 0$ is possible only when $j = N_r \/ 2$, and deduce that the origin is a node if and only if $N_r$ is even. (b) Conclude that the odd-$N_r$ convention of this chapter never evaluates the $1\/r$ weight at $r = 0$, whereas even $N_r$ places a node exactly there. (c) Reconcile this with the remark in @sec-odd-nr that an alternative even-node convention also avoids the origin, by identifying which family of points (Gauss--Lobatto with endpoints included versus interior Gauss nodes) each statement assumes.
] <ex-polar-node-origin>

#exercise(title: [Action of the Folded Radial Operator])[
  Let a grid function obey the symmetry condition @eq-symmetry-condition, so that its values on the negative-$r$ nodes are copies of the positive-$r$ values shifted by $pi$ in angle. Starting from the interior second-derivative matrix partitioned as @eq-D2-blocks, show that the second radial derivative at the positive-$r$ nodes is reproduced by $D_1$ acting on the same-angle data together with $D_2$ acting on the $pi$-shifted data. (a) Explain why the negative-$r$ columns must be indexed in reverse, so that column $j$ corresponds to $-r_j = r_(N_r - j)$. (b) Justify discarding the lower block row of @eq-D2-blocks entirely. (c) Repeat the argument for the first-derivative blocks @eq-D1-blocks and the weight $R$ to recover the same-angle operator $(D_1 + R E_1)$ and its $pi$-shifted partner $(D_2 + R E_2)$ of @sec-radial-operators.
] <ex-polar-block-action>

#exercise(title: [Properties of the Block Swap Matrix])[
  The swap matrix $S$ of @eq-swap-matrix encodes the angular $pi$-shift of the doubling trick. (a) Show that $S$ is a symmetric permutation matrix that squares to the identity, $S = S^top$ and $S^2 = I_M$, hence an involution with eigenvalues $plus.minus 1$. (b) Applied to the vector of nodal values on the Fourier grid @eq-theta-grid, show that $S$ carries $theta_m$ to the antipodal angle $theta_(m + M \/ 2)$, and explain why $M$ must be even for this to be an exact permutation of nodes. (c) Deduce that $S$ commutes with the periodic angular second-derivative matrix $D_theta^((2))$, and interpret this as the rotational consistency of the assembled Laplacian @eq-polar-laplacian-discrete.
] <ex-polar-swap-properties>

#exercise(title: [Bessel Eigenvalues and Their Multiplicities])[
  Carry out the separation of variables $u = R(r) thin Theta(theta)$ for the disk eigenproblem @eq-disk-eigenproblem in full. (a) Show that single-valuedness in $theta$ forces the separation constant to equal $m^2$ for a non-negative integer $m$, with $Theta(theta) in {cos m theta, sin m theta}$. (b) Reduce the radial equation to Bessel's equation and impose $R(1) = 0$ to obtain the eigenvalues @eq-bessel-eigenvalues, $lambda_(m, n) = j_(m, n)^2$. (c) Explain why every $m gt.eq.slant 1$ eigenvalue has multiplicity two while the $m = 0$ eigenvalues are simple, and state the number of nodal diameters and nodal circles of the mode $(m, n)$.
] <ex-polar-bessel-degeneracy>

#exercise(title: [Rotation Within a Degenerate Eigenspace])[
  For a degenerate pair with $m gt.eq.slant 1$, let $u_1 = J_m (j_(m, n) r) cos m theta$ and $u_2 = J_m (j_(m, n) r) sin m theta$. Using the mixing @eq-rotation-combination, (a) apply the identity $cos phi cos m theta + sin phi sin m theta = cos(m theta - phi)$ to show that $u_phi = J_m (j_(m, n) r) cos(m theta - phi)$. (b) Deduce that the $m$ nodal diameters rotate rigidly by $phi \/ m$, so that as $phi$ runs from $0$ to $pi$ the nodal pattern is carried back onto itself. (c) Explain why a numerical eigensolver, confronted with this two-dimensional eigenspace, returns an essentially arbitrary value of $phi$, and why this reflects the rotational symmetry of the disk faithfully rather than signalling a defect, as illustrated in @etude-degenerate-modes.
] <ex-polar-eigenspace-rotation>

#hint-for(<ex-polar-eigenspace-rotation>)[Combine the two modes with the addition formula $cos phi cos m theta + sin phi sin m theta = cos(m theta - phi)$; the mixing angle $phi$ then enters only as a rigid phase shift of the angular factor, hence a rotation of the nodal lines.]

=== Computational Exercises

#exercise(title: [The Even-versus-Odd $N_r$ Trap])[
  Modify the eigenmode computation of @etude-disk-eigenmodes (script `disk_eigenmodes`) to use _even_ values of $N_r$, so that a grid point falls exactly at $r = 0$. Plot the condition number of the Laplacian matrix $L$ versus $N_r$ for both even and odd $N_r$. Explain why even $N_r$ leads to a singular or nearly singular system.
] <ex-polar-even-odd-nr>

#exercise(title: [CFL Scaling for the Wave Equation])[
  Implement the wave equation $u_(t t) = c^2 Delta u$ on the disk using leapfrog time stepping and the polar Laplacian from this chapter. For each $N_r in {11, 15, 21, 25, 31}$, experimentally determine the maximum stable time step $Delta t_(max)$ by running the solver until either $t = 1$ or the solution blows up. Plot $Delta t_(max)$ versus $N_r$ on a log-log scale and verify the $O(N_r^(-2))$ scaling anticipated in @etude-polar-grids.
] <ex-polar-cfl-wave>

#exercise(title: [Eigenvalues on an Annulus])[
  Compute the first 10 eigenvalues of $-Delta u = lambda u$ on the annulus#idx("annulus") $0.5 lt.eq.slant r lt.eq.slant 1$ with $u = 0$ on both boundaries. The standard Chebyshev grid on $[-1, 1]$ maps to $[0.5, 1]$ via $r = 0.75 + 0.25 x$. No doubling trick is needed. Compare your results with the analytical eigenvalues obtained via separation of variables, which involve both $J_m$ and $Y_m$ Bessel functions.
] <ex-polar-annulus-eigs>

#exercise(title: [Variable-Coefficient Eigenproblem])[
  Suppose the eigenvalue problem is modified to $-Delta u = lambda (1 + r cos theta \/ 2) u$, modelling a membrane of non-uniform density. Compute the first six eigenvalues and plot the corresponding eigenmodes. How does the spatial variation of density affect the nodal line patterns compared with the uniform case?
] <ex-polar-variable-density>

#exercise(title: [A Manufactured Solution for the Poisson Solver])[
  Construct an exact test for the Poisson solver of @etude-disk-poisson (script `disk_poisson`) by the method of manufactured solutions#idx("manufactured solution"). Choose $u_("exact") (x, y) = (1 - x^2 - y^2) thin v(x, y)$ for a smooth factor $v$ of your choosing, so that $u_("exact")$ vanishes on $r = 1$, and compute the forcing $f = -Delta u_("exact")$ analytically. (a) Take $v equiv 1$ (radially symmetric) and verify that $||u_h - u_("exact")||_infinity$ falls to machine precision and is essentially independent of the angular resolution $M$. (b) Take $v = e^x cos y$, which carries genuine angular structure, and plot the error against $N_r$ at fixed large $M$ and against $M$ at fixed large $N_r$ on semilogarithmic axes. (c) Identify the resolution at which the radial and angular errors balance, and relate the two exponential rates to the smoothness of $u_("exact")$ in $r$ and in $theta$.
] <ex-polar-poisson-mms>

#exercise(title: [Stability and Steady State of the Heat Solver])[
  Using the Crank--Nicolson heat solver of @etude-disk-heat (script `disk_heat`), study the interplay of stability and accuracy for @eq-heat-disk. (a) Confirm the unconditional stability of @eq-crank-nicolson-polar by integrating with a sequence of increasingly large time steps $Delta t$ and checking that the late-time field always settles to the same steady state, namely the direct solution of $alpha Delta u = -S$ from a single linear solve. (b) Replace Crank--Nicolson by explicit forward Euler and locate, by bisection in $Delta t$, the stability threshold; verify that it scales like $O(N_r^(-4))$ as the radial resolution grows. (c) Compare the wall-clock cost of reaching steady state with the two integrators at a tolerance of your choice, and explain why the implicit scheme wins despite its per-step linear solve.
] <ex-polar-heat-stability>

=== Project-Style Exercises

#exercise(title: [Banded Operators on the Disk])[
  The dense Kronecker Laplacian @eq-polar-laplacian-discrete costs $O(N^3)$ to factorise. The ultraspherical spectral method#idx("ultraspherical spectral method") of @OlverTownsend2013, extended to the disk through the sparse Jacobi constructions of Shen, Tang, and Wang @ShenTangWang2011 and the spin-weighted bases of Vasil and collaborators @Vasil2016, replaces it by strictly banded operators. (a) Study one of these formulations and implement a banded radial operator for the disk Poisson problem. (b) Compare its conditioning and factorisation cost against the dense collocation Laplacian of this chapter as $N$ grows. (c) Report the empirical scaling of the solve time and contrast it with the dense $O(N^3)$ baseline.
] <ex-polar-banded>

#hint-for(<ex-polar-banded>)[The enabling fact is that differentiating a Chebyshev polynomial yields an ultraspherical polynomial, turning the dense radial derivative into a banded map between consecutive Gegenbauer bases; compose it with the sparse conversion operators that step between adjacent bases.]

#exercise(title: [A Zernike Galerkin Solver])[
  Zernike polynomials#idx("Zernike polynomials") @Zernike1934 are orthogonal on the unit disk and embed the pole conditions $u_m (r) tilde.op r^(|m|)$ automatically, so a Galerkin discretisation in this basis needs no doubling trick. Following the comparative study of Boyd and Yu @BoydYu2011 and the sparse Galerkin construction of Shen @Shen1997, (a) build the Zernike basis and assemble the disk Laplacian in it. (b) Solve the membrane eigenproblem @eq-disk-eigenproblem and confirm that the computed eigenvalues reproduce the squared Bessel zeros @eq-bessel-eigenvalues to spectral accuracy. (c) Compare the conditioning and sparsity of the Galerkin operator against the collocation Laplacian of @sec-polar-laplacian.
] <ex-polar-zernike>

#hint-for(<ex-polar-zernike>)[Write each Zernike mode as $R_n^m (r) thin e^(i m theta)$ with $R_n^m$ a Jacobi polynomial in $2 r^2 - 1$; the built-in $r^(|m|)$ factor cancels the coordinate singularity, so the radial mass and stiffness integrals can be evaluated exactly by Gauss--Jacobi quadrature.]

#exercise(title: [A Three-Dimensional Cylinder Solver])[
  Extend the polar framework of @sec-polar-extensions to the closed cylinder $(r, theta, z)$ by adjoining a Chebyshev grid in the bounded axial variable $z$, so that the Laplacian becomes a three-way Kronecker product with the radial doubling trick retained. (a) Assemble the three-dimensional operator and solve the heat equation in the cylinder with a localised source. (b) Following Darrow @Darrow2023, replace the dense three-dimensional solve by an Alternating Direction Implicit (ADI) sweep and study how the cost moves toward the quasi-optimal $O(N log N)$ scaling. (c) Validate the steady state against a direct Poisson solve and report the break-even resolution at which the ADI approach overtakes dense factorisation.
] <ex-polar-cylinder-adi>

#hint-for(<ex-polar-cylinder-adi>)[ADI splits each implicit step into successive one-dimensional solves along $r$, $theta$, and $z$, every one of them banded or diagonalisable, so the per-step cost grows almost linearly in the number of unknowns instead of cubically.]

#exercise(title: [Structure-Preserving Operators for the Schrödinger Equation])[
  The linear Schrödinger equation $i u_t = -Delta u$ on the disk requires the discrete evolution to preserve the $L^2$ norm, yet the collocation Laplacian @eq-polar-laplacian-discrete is non-normal and can leak energy. (a) Quantify the non-normality of $L$ (for instance through its departure from normality, or the conditioning of its eigenvector matrix) and track how it grows with $N_r$. (b) Following Gao and Iserles @GaoIserles2025, and in the spirit of the structure-preserving operators noted in @sec-polar-summary, build a skew-symmetric#idx("skew-symmetric matrix") differentiation matrix on the doubled radial grid and integrate the Schrödinger equation with a norm-conserving (Cayley) scheme. (c) Track the $L^2$ norm over many periods for both discretisations and quantify the artificial dissipation introduced by the non-skew operator.
] <ex-polar-schrodinger-skew>

#hint-for(<ex-polar-schrodinger-skew>)[Norm conservation follows when the spatial operator is skew-Hermitian, since the Cayley (Crank--Nicolson) update $(I - tau L \/ 2)^(-1) (I + tau L \/ 2)$ is then unitary; the W-function construction of @GaoIserles2025 produces exactly such a skew-symmetric radial derivative on the doubled grid.]
