// textbook/chapters/polar_coordinates.typ
// Chapter 12: Spectral Methods in Polar Coordinates
// Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
// Email: denys.dutykh@ku.ac.ae
// Homepage: https://www.denys-dutykh.com/
// Last modified: February 2026

#import "../styles/template.typ": dropcap, num, format-table

// Enable equation numbering for this chapter

= Spectral Methods in Polar Coordinates <ch-polar>

#dropcap[The spectral machinery assembled in the preceding chapters --- Chebyshev differentiation for bounded, non-periodic intervals and Fourier methods for periodic domains --- has so far been deployed on rectangular geometries: intervals, squares, and periodic strips. But many problems of scientific interest live on _circular_ domains: vibrating drumheads, optical fibres with circular cross-section, semiconductor wafers under localised heating, quantum particles confined to circular corrals, and incompressible flow in pipes. In this chapter, we extend the spectral framework to the _unit disk_ by introducing polar coordinates $(r, theta)$ and combining Chebyshev discretisation in the radial direction with Fourier discretisation in the angular direction. The principal new challenge is the _coordinate singularity_ at $r = 0$, which we overcome with an elegant algebraic device: the doubling trick.]

By the end of this chapter, you should be able to:

1. Understand the coordinate singularity at $r = 0$ and why a naive Chebyshev mapping to $[0, 1]$ wastes resolution near the origin.
2. Apply the _doubling trick_ --- extending $r$ from $[0, 1]$ to $[-1, 1]$ with the symmetry condition $u(r, theta) = u(-r, theta + pi)$ --- to obtain a well-conditioned discretisation that avoids the origin entirely.
3. Assemble the discrete polar Laplacian via Kronecker products of block-decomposed Chebyshev matrices and Fourier angular operators.
4. Compute eigenmodes of the Laplacian on the disk and compare numerical eigenvalues with the zeros of Bessel functions to verify spectral accuracy.
5. Solve elliptic (Poisson) and parabolic (heat) PDEs on circular domains using the polar spectral framework.

== Why Polar? Motivation and Challenges <sec-polar-motivation>

=== Physical Domains with Circular Geometry

Circular and cylindrical geometries arise across science and engineering. The vibration of a circular drumhead --- one of the oldest problems in mathematical physics --- requires solving the wave equation on a disk @MorseIngard1968. In fibre optics, the fundamental modes of a cylindrical waveguide are governed by the Helmholtz equation in a circular cross-section. Heat conduction in semiconductor wafers, which are manufactured as thin circular disks, naturally calls for solving the heat equation in polar coordinates. In quantum mechanics, the so-called _quantum corral_ --- a circular arrangement of atoms on a surface confining electron standing waves --- yields striking eigenmode patterns that can be directly compared with Bessel function predictions @Crommie1993. Incompressible pipe flow, the Couette--Taylor instability between rotating cylinders, vortex dynamics in geophysical flows, and the high-precision modelling of gravitational waves in numerical relativity all involve polar, cylindrical, or spherical domains.

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

The grid construction is straightforward in all three languages. In Python:

```python
D, x = cheb_matrix(Nr)
r_naive = (x + 1) / 2        # naive: map to [0, 1]
r_doubled = x                # doubled: keep [-1, 1]
N2 = (Nr - 1) // 2
r_pos = r_doubled[1:N2+1]    # positive-r interior points only
```

The MATLAB equivalent:

```matlab
[D, x] = cheb_matrix(Nr);
r_naive = (x + 1) / 2;
r_doubled = x;
N2 = (Nr - 1) / 2;
r_pos = r_doubled(2:N2+1);
```

And in Julia:

```julia
D, x = cheb_matrix(Nr)
r_naive = (x .+ 1) / 2
r_doubled = x
N2 = (Nr - 1) ÷ 2
r_pos = r_doubled[2:N2+1]
```

#figure(
  image("../figures/ch12/python/polar_grids.pdf", width: 95%),
  caption: [Two spectral grids on the unit disk with $N_r = 25$ radial and $M = 20$ angular points. _Left_: the naive grid maps Chebyshev points to $r in [0, 1]$, resulting in heavy clustering near the origin. The dashed circle (containing half the grid points) encloses only 25% of the disk area. _Right_: the doubled grid uses Chebyshev points on $r in [-1, 1]$ with odd $N_r$, keeping only the positive-$r$ half. No point falls at the origin, and the resolution is more uniformly distributed.],
) <fig-polar-grids>

@fig-polar-grids reveals the dramatic difference in grid point distribution. The naive grid devotes an excessive number of points to the small neighbourhood of the origin, while the doubled grid distributes resolution more evenly. @fig-polar-area-cfl quantifies the consequences: the area per radial annulus is far more uniform for the doubled grid, and the maximum stable time step scales as $O(N_r^(-2))$ instead of $O(N_r^(-4))$.

#figure(
  image("../figures/ch12/python/polar_area_cfl.pdf", width: 95%),
  caption: [_Left_: area of each radial annulus for the naive and doubled grids. The naive grid allocates disproportionate area to the thin annuli near $r = 0$. _Right_: CFL time-step scaling. The naive grid forces $Delta t_(max) tilde.op N_r^(-4)$ (blue circles), while the doubled grid gives $Delta t_(max) tilde.op N_r^(-2)$ (red squares), a substantial improvement.],
) <fig-polar-area-cfl>

The code generating @fig-polar-grids and @fig-polar-area-cfl is available in:
- `codes/python/ch12/polar_grid_geometry.py`
- `codes/matlab/ch12/polar_grid_geometry.m`
- `codes/julia/ch12/polar_grid_geometry.jl`

=== Discussion

The visual contrast in @fig-polar-grids is striking and immediately conveys what abstract complexity estimates alone cannot: the naive Chebyshev mapping to $[0, 1]$ concentrates roughly half its grid points inside a circle enclosing barely a third of the disk area. This is not merely an aesthetic deficiency --- it has direct consequences for time-stepping. @fig-polar-area-cfl quantifies the penalty: the CFL-limited time step scales as $O(N_r^(-4))$ for the naive grid versus the far milder $O(N_r^(-2))$ for the doubled grid. In practice, this means that doubling the radial resolution on the naive grid forces a _sixteen-fold_ reduction in the time step, making explicit time integration prohibitively expensive even at moderate resolution.

The doubled grid achieves its advantage by a simple algebraic device --- keeping $r in [-1, 1]$ and discarding the redundant half --- yet the computational savings are profound. This étude demonstrates a recurring theme in scientific computing: the choice of coordinate mapping is not a mere technicality but can dominate the practical efficiency of a method. The quantitative comparison here provides the motivation for the doubling trick developed in the next section.

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
  *Remark.* An alternative convention, used by some authors, takes $N_r$ _even_ and avoids $r = 0$ for the same reason: with even $N_r$, the points $r_j = cos(j pi \/ N_r)$ also skip zero. Both conventions work; we follow the odd-$N_r$ convention of Trefethen @Trefethen2000 and Fornberg @Fornberg1996 throughout this chapter. Exercise 12.1 explores the consequences of using even $N_r$. It must be noted, however, that while the doubling trick bypasses physical pole evaluation, high-resolution simulations of nonlinear PDEs can still suffer from aliasing and stability issues near the origin unless explicit parity restrictions are enforced.
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
where $I_M$ is the $M times M$ identity, $S$ is the _block swap matrix_
$ S = mat(0, I_(M\/2); I_(M\/2), 0) $ <eq-swap-matrix>
that implements the angular $pi$-shift from the symmetry condition, and $R^2 = op("diag")(1 \/ r_j^2)$.

The matrix $L$ has dimensions $(N_2 M) times (N_2 M)$. For typical parameters ($N_r = 25$, $M = 20$, giving $N_2 = 12$), this is a $240 times 240$ system --- modest enough for dense linear algebra, yet large enough to achieve spectral accuracy.

The structure of @eq-polar-laplacian-discrete closely parallels the tensor-product Laplacian on a rectangle from @ch-spectral-pde: the first two Kronecker products handle radial differentiation (with the symmetry folding), and the third handles angular differentiation. The only novelty is the swap matrix $S$, which encodes the 2-to-1 map.

=== Implementation: The `laplacian_polar` Function <sec-laplacian-impl>

We encapsulate the assembly into a reusable function. The Python implementation:

```python
def laplacian_polar(Nr, M):
    """Build the discrete polar Laplacian on the unit disk."""
    N2 = (Nr - 1) // 2
    D, r = cheb_matrix(Nr)
    D2full = D @ D

    # Block decomposition (reversed negative-r columns)
    D1 = D2full[1:N2+1, 1:N2+1]
    D2b = D2full[1:N2+1, Nr-1:N2:-1]
    E1 = D[1:N2+1, 1:N2+1]
    E2 = D[1:N2+1, Nr-1:N2:-1]

    r_pos = r[1:N2+1]
    R = np.diag(1.0 / r_pos)

    # Fourier second derivative (Toeplitz matrix)
    dt = 2 * np.pi / M
    col = np.zeros(M)
    col[0] = -np.pi**2 / (3*dt**2) - 1/6
    for j in range(1, M):
        col[j] = 0.5 * (-1)**(j+1) / np.sin(dt*j/2)**2
    D2t = toeplitz(col)

    # Kronecker assembly
    M2 = M // 2
    S = np.block([[np.zeros((M2,M2)), np.eye(M2)],
                  [np.eye(M2), np.zeros((M2,M2))]])
    L = (np.kron(D1 + R@E1, np.eye(M))
         + np.kron(D2b + R@E2, S)
         + np.kron(R**2, D2t))
    return L, r_pos, theta
```

The MATLAB equivalent:

```matlab
function [L, r_pos, theta] = laplacian_polar(Nr, M)
    N2 = (Nr - 1) / 2;
    [D, r] = cheb_matrix(Nr);
    D2 = D^2;
    D1 = D2(2:N2+1, 2:N2+1);
    D2b = D2(2:N2+1, Nr:-1:N2+2);
    E1 = D(2:N2+1, 2:N2+1);
    E2 = D(2:N2+1, Nr:-1:N2+2);
    r_pos = r(2:N2+1);
    R = diag(1 ./ r_pos);
    dt = 2*pi/M;  M2 = M/2;
    col = [-pi^2/(3*dt^2)-1/6, ...
           0.5*(-1).^(2:M)./sin(dt*(1:M-1)/2).^2]';
    D2t = toeplitz(col);
    Z = zeros(M2); I = eye(M2);
    L = kron(D1+R*E1,eye(M)) + kron(D2b+R*E2,[Z I;I Z]) ...
        + kron(R^2, D2t);
    theta = dt*(0:M-1)';
end
```

And in Julia:

```julia
function laplacian_polar(Nr::Int, M::Int)
    N2 = (Nr - 1) ÷ 2
    D, r = cheb_matrix(Nr)
    D2full = D * D
    D1  = D2full[2:N2+1, 2:N2+1]
    D2b = D2full[2:N2+1, Nr:-1:N2+2]
    E1  = D[2:N2+1, 2:N2+1]
    E2  = D[2:N2+1, Nr:-1:N2+2]
    r_pos = r[2:N2+1]
    R = diagm(1.0 ./ r_pos)
    dt = 2π / M;  M2 = M ÷ 2
    col = zeros(M)
    col[1] = -π^2/(3*dt^2) - 1/6
    for j in 1:M-1
        col[j+1] = 0.5*(-1)^(j+1)/sin(dt*j/2)^2
    end
    D2t = [col[abs(i-j)+1] for i in 1:M, j in 1:M]
    Z = zeros(M2, M2);  Ih = I(M2, M2)
    S = [Z Ih; Ih Z]
    L = kron(D1+R*E1, I(M,M)) + kron(D2b+R*E2, S) + kron(R^2, D2t)
    return L, r_pos, collect(range(0, step=dt, length=M))
end
```

Note the critical detail in the block extraction: the negative-$r$ columns are indexed in _reverse order_ (`Nr-1:N2:-1` in 0-based Python, `Nr:-1:N2+2` in 1-based MATLAB/Julia). This reversal ensures that the $j$-th positive-$r$ point is correctly paired with its symmetric counterpart $-r_j$.

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

The first few Bessel zeros and eigenvalues are:
- $j_(0, 1) approx 2.4048$, so $lambda_(0, 1) approx 5.7832$ (fundamental mode, no nodal lines inside the disk).
- $j_(1, 1) approx 3.8317$, so $lambda_(1, 1) approx 14.6820$ (first non-radial mode, one nodal diameter).
- $j_(2, 1) approx 5.1356$, so $lambda_(2, 1) approx 26.3746$ (two nodal diameters).
- $j_(0, 2) approx 5.5201$, so $lambda_(0, 2) approx 30.4713$ (one nodal circle, no nodal diameters).

The spectral method does _not_ use Bessel functions: it discretises the Laplacian as a matrix and computes eigenvalues by dense linear algebra. The fact that the numerical eigenvalues match the Bessel zeros to full floating-point precision is a striking confirmation of spectral accuracy.

=== Computational Étude 12.2: Vibrations of a Circular Drumhead <etude-disk-eigenmodes>

We compute the first 25 eigenpairs of the Laplacian on the unit disk with $N_r = 25$ and $M = 20$, then compare the numerical eigenvalues with the exact values @eq-bessel-eigenvalues.

The core computation in Python:

```python
Nr, M = 25, 20
L, r_pos, theta = laplacian_polar(Nr, M)
eigenvalues, eigenvectors = np.linalg.eig(-L)
eigenvalues = np.sort(np.real(eigenvalues))
```

In MATLAB:

```matlab
Nr = 25; M = 20;
[L, r_pos, theta] = laplacian_polar(Nr, M);
[V, Lam] = eig(-L);
[lam, idx] = sort(diag(Lam));
V = V(:, idx);
```

In Julia:

```julia
Nr, M = 25, 20
L, r_pos, theta = laplacian_polar(Nr, M)
F = eigen(-L)
idx = sortperm(real.(F.values))
eigenvalues = real.(F.values[idx])
```

#figure(
  image("../figures/ch12/python/disk_eigenmodes.pdf", width: 95%),
  caption: [Six eigenmodes of the Laplacian on the unit disk, computed with $N_r = 25$ radial and $M = 20$ angular points. The modes display the characteristic Bessel function patterns: the fundamental mode ($m = 0$, $n = 1$) is a simple dome, while higher modes exhibit nodal diameters (radial lines) and nodal circles.],
) <fig-disk-eigenmodes>

#figure(
  image("../figures/ch12/python/eigenvalue_convergence.pdf", width: 95%),
  caption: [_Left_: the first 25 computed eigenvalues (circles) agree with the exact Bessel zeros squared (crosses) to plotting accuracy. The degenerate pairs ($m gt.eq.slant 1$) appear as overlapping markers. _Right_: the error in the fundamental eigenvalue $lambda_(0,1)$ decreases exponentially with $N_r$, confirming spectral convergence.],
) <fig-eigenvalue-convergence>

The code generating @fig-disk-eigenmodes and @fig-eigenvalue-convergence is available in:
- `codes/python/ch12/disk_eigenmodes.py`
- `codes/matlab/ch12/disk_eigenmodes.m`
- `codes/julia/ch12/disk_eigenmodes.jl`

=== Discussion

@fig-disk-eigenmodes displays the rich modal structure of the circular membrane: from the smooth axisymmetric dome of the fundamental mode to the intricate nodal patterns of higher modes featuring radial node lines and concentric nodal circles. These patterns, though analytically predicted by Bessel function theory since the 19th century, are here recovered _without_ any explicit reference to Bessel functions --- the spectral discretisation of the Laplacian and a call to a standard eigensolver suffice. The convergence plot in @fig-eigenvalue-convergence is particularly telling: the error in the fundamental eigenvalue drops exponentially with $N_r$, reaching machine precision with remarkably few grid points. This spectral (exponential) convergence rate is characteristic of Chebyshev methods applied to smooth problems, as discussed in @ch-spectral-pde, and it provides a powerful cross-validation: the agreement of numerical eigenvalues with the squared Bessel zeros $j_(m,n)^2$ to 10 or more digits confirms that both the doubling trick and the Kronecker product assembly are implemented correctly.

A subtle but important observation concerns the degenerate eigenvalue pairs visible in @fig-eigenvalue-convergence (left panel): for each angular wave number $m gt.eq.slant 1$, the eigenvalue $lambda_(m,n)$ appears twice, corresponding to the $cos m theta$ and $sin m theta$ modes. The numerical eigensolver returns these as a pair of nearly identical eigenvalues, but the _individual_ eigenvectors are arbitrary rotations within the two-dimensional eigenspace. This degeneracy is explored further in the next étude.

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

=== Discussion

@fig-nodal-rotation reveals a phenomenon that is invisible in eigenvalue tables but physically fundamental: within a degenerate eigenspace, the nodal line pattern is not fixed but can rotate freely. As the mixing angle $phi$ sweeps from $0$ to $pi$ in @eq-rotation-combination, the single nodal diameter glides smoothly around the disk, and every orientation is an equally valid eigenmode with the same frequency. This continuous rotational freedom is a direct consequence of the circular symmetry of the domain --- the Laplacian commutes with rotations, and the two-dimensional eigenspace for each $m gt.eq.slant 1$ carries an irreducible representation of the rotation group.

From a computational perspective, this étude highlights an important subtlety of numerical eigensolvers: when an eigenvalue has multiplicity greater than one, the returned eigenvectors are _not_ uniquely determined. Different runs, different algorithms, or even different floating-point orderings may produce different rotations within the eigenspace. The visualisation in @fig-nodal-rotation demonstrates that this ambiguity is not a defect but a faithful reflection of the underlying physics. Any experiment that breaks the circular symmetry --- a slight deformation of the boundary, a non-uniform density, or a localised perturbation --- would lift the degeneracy and select a preferred orientation, as explored in Exercise 12.4.

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

The core solve in Python:

```python
Nr, M = 31, 40
L, r_pos, theta = laplacian_polar(Nr, M)
RR, TT = np.meshgrid(r_pos, theta, indexing='ij')
XX, YY = RR * np.cos(TT), RR * np.sin(TT)
f = np.exp(-40 * ((XX - 0.4)**2 + (YY - 0.2)**2))
u = np.linalg.solve(L, -f.flatten()).reshape(N2, M)
```

In MATLAB:

```matlab
Nr = 31; M = 40;
[L, r_pos, theta] = laplacian_polar(Nr, M);
[RR, TT] = meshgrid(r_pos, theta);
f = exp(-40*((RR'.*cos(TT') - 0.4).^2 + (RR'.*sin(TT') - 0.2).^2));
u = reshape(L \ (-f(:)), [], M);
```

In Julia:

```julia
Nr, M = 31, 40
L, r_pos, theta = laplacian_polar(Nr, M)
RR = [r for r in r_pos, θ in theta]
TT = [θ for r in r_pos, θ in theta]
f = exp.(-40 .* ((RR.*cos.(TT) .- 0.4).^2 .+ (RR.*sin.(TT) .- 0.2).^2))
u = reshape(L \ (-vec(f)), :, M)
```

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

=== Discussion

The Poisson equation on the disk provides a stringent test of the polar spectral framework beyond eigenproblems, because the off-centre Gaussian source @eq-gaussian-source breaks every symmetry the disk possesses. @fig-poisson-disk-surface and @fig-poisson-disk-contour show that the spectral method resolves the localised peak and the smooth decay to the clamped boundary without spurious oscillations --- a hallmark of spectral accuracy for smooth solutions. The distortion of the equipotential lines toward the source location in @fig-poisson-disk-contour is physically intuitive: a finger pressing off-centre on a clamped drumhead deflects the membrane asymmetrically, with steeper gradients on the side closer to the boundary.

The radial symmetry test in @fig-radial-symmetry-test is an especially valuable diagnostic. When the source is centred, all radial profiles at different angles collapse onto a single curve to machine precision, confirming that the spectral discretisation introduces no artificial anisotropy. This is a non-trivial check: a method that mishandles the coordinate singularity or the angular coupling through the swap matrix $S$ would produce angle-dependent artefacts even for a radially symmetric problem. The clean collapse in @fig-radial-symmetry-test validates both the Kronecker product assembly @eq-polar-laplacian-discrete and the block decomposition from @sec-block-decomposition. For the off-centre source, the fan of radial profiles in @fig-radial-symmetry-test encodes the same geometric information as the contour plot in @fig-poisson-disk-contour, but in a form more amenable to quantitative reading: the angle at which the fan is widest identifies the source direction, and the radial coordinate at which the dominant profile peaks recovers the source's distance $r_s = sqrt(x_0^2 + y_0^2) approx 0.45$ from the origin. This correspondence between the analytical geometry of the source and the shape of the numerical solution confirms that the method captures not only global accuracy but also local directional information correctly.

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

The Crank--Nicolson time loop in Python:

```python
I_mat = np.eye(ndof)
A = I_mat - alpha * dt / 2 * L
B = I_mat + alpha * dt / 2 * L
lu, piv = lu_factor(A)        # factorise once
u = np.zeros(ndof)             # cold start
for step in range(n_steps):
    u = lu_solve((lu, piv), B @ u + dt * S_vec)
```

The MATLAB equivalent:

```matlab
A = eye(ndof) - alpha*dt/2 * L;
B = eye(ndof) + alpha*dt/2 * L;
[Lf, Uf] = lu(A);
u = zeros(ndof, 1);
for step = 1:n_steps
    u = Uf \ (Lf \ (B*u + dt*S_vec));
end
```

And in Julia:

```julia
A = I(ndof) - α * dt / 2 * L
B = I(ndof) + α * dt / 2 * L
F = lu(A)
u = zeros(ndof)
for step in 1:n_steps
    u = F \ (B * u .+ dt .* S_vec)
end
```

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

=== Discussion

This étude ties together the spatial machinery of the chapter with the time-stepping techniques from @ch-spectral-pde, demonstrating that the polar Laplacian $L$ slots seamlessly into the method-of-lines framework. The snapshots in @fig-heat-disk-snapshots tell a physically compelling story: the localised laser source rapidly heats a small region near $(0.3, 0)$, thermal diffusion spreads the heat outward, and the zero-temperature boundary condition acts as a heat sink that ultimately balances the input. The approach to steady state is quantified in @fig-heat-disk-energy (left panel), where the maximum temperature saturates after approximately $t approx 15$.

The most informative validation appears in @fig-heat-disk-energy (right panel): the steady-state radial profile obtained by marching the Crank--Nicolson scheme to large times agrees with the profile obtained by directly solving the steady Poisson equation $alpha Delta u = -S$. This cross-check is significant because the two computations exercise very different parts of the code --- the time-dependent solver requires LU factorisation and repeated forward/back substitution, while the direct solve is a single linear system inversion. Their agreement confirms both the temporal integrator and the spatial discretisation. The choice of Crank--Nicolson is deliberate: its unconditional stability for diffusion problems means the time step is limited only by accuracy, not by the stiffness of the discrete Laplacian, which is a practical advantage of implicit methods emphasised in @ch-spectral-pde.

== Extensions and Generalisations <sec-polar-extensions>

=== Annular Domains

The methods of this chapter adapt readily to the _annulus_ $r_("in") lt.eq.slant r lt.eq.slant r_("out")$, which models, for example, the space between two concentric pipes or the cross-section of a fibre-optic cladding. The Chebyshev grid on $x in [-1, 1]$ is mapped to $[r_("in"), r_("out")]$ by the affine transformation
$ r = frac(r_("out") + r_("in"), 2) + frac(r_("out") - r_("in"), 2) x. $
Since $r_("in") > 0$, the coordinate singularity at $r = 0$ is absent, and the doubling trick is not needed. The Kronecker product structure of the Laplacian carries over unchanged, with $R = op("diag")(1\/r_j)$ evaluated at the mapped grid points. Dirichlet (or Neumann, or Robin) conditions are imposed independently at both boundaries. Exercise 12.3 explores this setting.

=== Neumann and Robin Boundary Conditions

Throughout this chapter, we have imposed Dirichlet conditions $u = 0$ at $r = 1$. Other boundary conditions require modifying the radial differentiation matrices. For the Neumann condition $partial u \/ partial r |_(r=1) = 0$, one replaces the boundary rows of the Chebyshev matrix with the entries of the first-derivative matrix $D$, following the _boundary bordering_ technique from @ch-bvp. Robin conditions $alpha u + beta partial u \/ partial r = g$ at $r = 1$ are handled similarly.

=== Cylindrical and Spherical Geometries

Adding a third variable $z in [-H, H]$ to the polar $(r, theta)$ system gives _cylindrical coordinates_ $(r, theta, z)$. If $z$ is bounded and non-periodic, it gets a Chebyshev grid; if periodic, a Fourier grid. The Laplacian becomes a three-way Kronecker product. Recent advancements have optimised this setting: Darrow @Darrow2023 developed a quasi-optimal Chebyshev--Chebyshev--Fourier (CCF) solver that applies the radial doubling trick together with Alternating Direction Implicit (ADI) methods, reducing the computational complexity of 3D cylinder solves from $O(N^(4\/3))$ to $O(N log N)$.

For _spherical coordinates_ $(r, theta, phi)$, the angular part requires _spherical harmonics_ rather than simple Fourier modes. An alternative is to treat the polar angle $phi in [0, pi]$ with a Chebyshev (or Gauss--Legendre) grid and the longitude $theta in [0, 2 pi)$ with Fourier, using parity arguments at the poles analogous to our treatment of $r = 0$.

The Galerkin approach using _Zernike polynomials_ or _Jacobi polynomials_ (specifically, the so-called "one-sided Jacobi" bases) provides another route: basis functions can be chosen to satisfy the pole conditions automatically, eliminating the need for the doubling trick entirely @Shen1995 @Shen1997 @ShenTangWang2011. Furthermore, for complex tensorial systems like the Navier--Stokes equations, _spin-weighted Jacobi polynomials_ act as exact eigenfunctions for covariant derivatives, preventing the catastrophic mixing of vector components at the pole that arises from Christoffel symbols in curvilinear coordinates @Vasil2016. For highly complex domains with piecewise-smooth boundaries, polar spectral-element methods utilising _Log Orthogonal Functions_ are now deployed to capture corner singularities while maintaining spectral accuracy @ChenSunWu2025.

== A Non-Exhaustive Literature Overview

The application of spectral methods to polar coordinates and circular domains has a rich, complex, and evolving history, characterised primarily by a continuous mathematical struggle against the coordinate singularity at the origin. The foundational formulation presented in this chapter --- utilising a doubled radial variable $r in [-1, 1]$ paired with a geometric symmetry constraint --- was systematically developed, rigorously analysed, and popularised by Fornberg @Fornberg1995 @Fornberg1996 @Fornberg1996. Trefethen @Trefethen2000 distilled the doubling trick into a highly compact, elegant algebraic implementation (Programs 28 and 29 in his book), which forms the core architecture of the computational études presented here. Historically, the challenge of the pole problem was treated by explicitly imposing analytical regularity conditions on the spectral expansion. As identified in the seminal works of Orszag @Orszag1974 in the context of spherical atmospheric modelling, physical smoothness at the pole dictates that the $m$-th Fourier coefficient must decay as $O(r^(|m|))$. Boyd @Boyd2000 provides a comprehensive survey of these early methodologies, cataloguing dozens of approaches ranging from intricate coordinate transformations to the explicit, and often ill-conditioned, enforcement of parity boundary constraints at the origin. Matsushima and Marcus @Matsushima1995 provided one of the first highly robust spectral methods specifically tailored for vortex dynamics in polar coordinates, demonstrating that proper polynomial selection prevents the severe $O(N_r^(-4))$ CFL restrictions typically associated with radial point clustering.

An important analytical alternative to the collocation approach is the Galerkin method, which utilises global basis functions constructed to satisfy the pole conditions inherently. The most natural basis for the unit disk is the _Zernike polynomial_ family. Originally formulated to quantify optical wavefront aberrations @Zernike1934, Zernike polynomials naturally embed the requisite $O(r^(|m|))$ pole conditions. A thorough evaluation of spectral bases on the disk, including Zernike, Logan--Shepp ridge polynomials, and shifted Chebyshev series, was conducted by Boyd and Yu @BoydYu2011, highlighting the persistent trade-offs between mathematical elegance and computational matrix conditioning. To resolve the issue of dense differentiation matrices, Shen @Shen1997 introduced the use of "one-sided Jacobi" polynomials to construct optimally sparse Galerkin operators, drastically reducing the computational cost of direct solvers in cylindrical domains; the comprehensive framework of Shen, Tang, and Wang @ShenTangWang2011 extended these ideas further, yielding strictly banded operator matrices. These Galerkin methods are preferred in large-scale three-dimensional computations where the $O(N^3)$ cost of dense linear algebra becomes prohibitive.

The modern era of spectral methods (post-2010) has been defined by the pursuit of strictly sparse, well-conditioned operators that scale quasi-optimally in high dimensions. Olver and Townsend @OlverTownsend2013 revolutionised the field by introducing the _ultraspherical spectral method_: by shifting the polynomial basis after differentiation --- recognising that the derivative of a Chebyshev polynomial is an ultraspherical polynomial --- they generated strictly banded operators, reducing the dense $O(N^3)$ complexity of traditional Chebyshev boundary value problems to a highly efficient $O(N log N)$. This philosophy of banded sparsity was masterfully extended to polar and spherical geometries by Vasil, Burns, Lecoanet, Olver, Brown, and Oishi @Vasil2016, who utilised spin-weighted Jacobi polynomials to perform exact tensor calculus at the pole, avoiding the catastrophic mixing of vector components induced by Christoffel symbols. Their methodology forms the mathematical backbone of the widely used _Dedalus_ simulation framework. Recently, Darrow @Darrow2023 applied these sparse techniques to the three-dimensional closed cylinder, yielding a quasi-optimal Chebyshev--Chebyshev--Fourier solver that manages the radial doubling trick efficiently via Alternating Direction Implicit (ADI) methods.

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

*Exercise 12.1* (_The even-vs-odd $N_r$ trap_). Modify the eigenmode computation from Étude 12.2 to use _even_ values of $N_r$ (so that a grid point falls at $r = 0$). Plot the condition number of the Laplacian matrix $L$ versus $N_r$ for both even and odd $N_r$. Explain why even $N_r$ leads to a singular or nearly singular system.

*Exercise 12.2* (_CFL scaling for the wave equation_). Implement the wave equation $u_(t t) = c^2 Delta u$ on the disk using leapfrog time stepping and the polar Laplacian from this chapter. For each $N_r in {11, 15, 21, 25, 31}$, experimentally determine the maximum stable time step $Delta t_max$ by running the solver until either $t = 1$ or the solution blows up. Plot $Delta t_max$ versus $N_r$ on a log-log scale and verify the scaling.

*Exercise 12.3* (_Eigenvalues on an annulus_). Compute the first 10 eigenvalues of $-Delta u = lambda u$ on the annulus $0.5 lt.eq.slant r lt.eq.slant 1$ with $u = 0$ on both boundaries. The standard Chebyshev grid on $[-1, 1]$ maps to $[0.5, 1]$ via $r = 0.75 + 0.25 x$. No doubling trick is needed. Compare your results with the analytical eigenvalues obtained via separation of variables, which involve both $J_m$ and $Y_m$ Bessel functions.

*Exercise 12.4* (_Variable-coefficient eigenproblem_). Suppose the eigenvalue problem is modified to $-Delta u = lambda (1 + r cos theta \/ 2) u$, modelling a membrane of non-uniform density. Compute the first six eigenvalues and plot the corresponding eigenmodes. How does the spatial variation of density affect the nodal line patterns compared with the uniform case?
