// textbook/chapters/appendix_b_complex.typ
// Appendix B: Complex Analysis for Convergence Theory
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Last modified: April 2026

#import "../styles/template.typ": dropcap

= Complex Analysis for Convergence Theory <appendix-complex>

#dropcap[The convergence theory of spectral methods is a chapter of complex analysis dressed in numerical clothing.] Almost every quantitative claim about *how fast* a Chebyshev or Fourier expansion converges --- the $rho^(-N)$ Bernstein rate, the strip-analytic exponential decay $e^(-a N)$, the geometric trapezoidal-rule envelope $r^N$ of a $sech$-shaped integrand, the algebraic $1 \/ N^p$ rate near a corner singularity --- ultimately reduces to *the location of the nearest complex-plane singularity* of the underlying analytic continuation. This appendix collects the complex-analytic objects that the body of the volume names without pausing to define.

== Why this appendix <appendix-complex-why>

The rate-determining theorems of @ch-smoothness, @ch-fourier-grids, @ch-quadrature, and @ch-periodic-trap are stated as black boxes and proved (when at all) only in passing. The coordinate transformations of @ch-coord-transforms move complex-plane singularities deliberately --- closer to the real axis to expose them, away from it to cure them. The asymptotic-augmentation strategies of @ch-unbounded match the basis to the singular behaviour of the target function. The complex-plane *detour* around an interior pole, central to Étude 18.11, is itself a complex-analytic technique. And the Hilbert-transform identity of Étude 21.8 is the boundary-value theorem of Sokhotski and Plemelj. A reader who knows where each of these objects sits in the broader complex-analytic landscape will read the body more profitably; this appendix is a brief cartographic reference.

The treatment is informal: theorems are stated, references provided, but proofs are given only when they fit on a page or two. The standard book-length references are @Trefethen2013 (approximation theory and complex analysis as a unit), @Boyd2000 (the encyclopaedic survey), and @TrefethenWeideman2014 for the trapezoidal-rule perspective.

== The Bernstein ellipse <appendix-complex-bernstein>

For the classical Chebyshev approximation problem on $[-1, 1]$, the basic geometric object is the *Bernstein ellipse*. Let $rho > 1$. The Bernstein ellipse $E_rho$ is the image of the circle $|w| = rho$ under the Joukowski map (@appendix-complex-joukowski below):
$ E_rho = {z = (w + 1\/w) \/ 2 : |w| = rho} subset.eq CC. $ <eq-bernstein>
$E_rho$ is an ellipse with foci at $z = plus.minus 1$, semi-major axis $(rho + rho^(-1)) \/ 2$ along the real axis, and semi-minor axis $(rho - rho^(-1)) \/ 2$ along the imaginary axis. The parameter $rho$, sometimes called the *Bernstein parameter* of the ellipse, is the *sum of the semi-axes plus one*.

The relevance to spectral methods is the *Bernstein theorem*:

#block(stroke: 0.8pt + rgb(20, 45, 110), radius: 3pt, inset: 10pt,
  fill: rgb(248, 250, 254))[
  *Bernstein theorem (best polynomial approximation).*  Let $f$ be analytic in the open Bernstein ellipse $E_rho$ and continuous on its closure. Then the Chebyshev best-polynomial-approximation error of degree $N$ satisfies
  $ E_N (f) eq.def inf_(p_N in cal(P)_N) ||f - p_N||_(infinity, [-1, 1]) lt.eq.slant frac(2 thin M(rho), rho^N (rho - 1)), $ <eq-bernstein-theorem>
  where $M(rho) = max_(z in partial E_rho) |f(z)|$.
]

Two consequences of this theorem are quoted dozens of times in the body. First, Chebyshev *interpolation* on the Lobatto grid converges at the same geometric rate $cal(O)(rho^(-N))$, because the ratio between best-polynomial and Lobatto-interpolation error is bounded by the Lebesgue constant $Lambda_N tilde.op log N$ (see @ch-geometry). Second, the Bernstein parameter $rho$ is *largest* when the ellipse is just inside the nearest singularity of $f$: pushing singularities further into the complex plane (e.g.\ via the tanh map of @ch-coord-transforms) directly enlarges $rho$ and hence accelerates convergence.

For the Witch of Agnesi $1 \/ (1 + 4 x^2)$ of @ch-chebyshev (poles at $x = plus.minus i \/ 2$) the Bernstein ellipse passes through $plus.minus i \/ 2$, giving $rho approx 1 + 1 \/ 2 = 1.5$, and the predicted geometric rate $1.5^(-N)$ matches the measured convergence to within a constant. For the Runge function $1 \/ (1 + 25 x^2)$ of @ch-geometry (poles at $plus.minus i \/ 5$) the Bernstein parameter is $rho approx 1.105$, far closer to one and explaining the famously slow Chebyshev convergence on this benchmark.

== Strip analyticity <appendix-complex-strip>

For *periodic* problems on $[0, 2 pi]$ the analogous geometric object is the horizontal *strip of analyticity*. Let $a > 0$. A function $f: RR arrow CC$ is *strip-analytic of width $a$* if it extends to an analytic function in the open horizontal strip
$ S_a = {z in CC : |op("Im") z| < a} subset.eq CC. $ <eq-strip>
The Paley--Wiener theorem (1934) establishes that strip analyticity implies *exponential* decay of Fourier coefficients:

#block(stroke: 0.8pt + rgb(20, 45, 110), radius: 3pt, inset: 10pt,
  fill: rgb(248, 250, 254))[
  *Theorem (Paley--Wiener, periodic version).*  Let $f$ be $2 pi$-periodic and analytic in $S_a$ with $sup_(|op("Im") z| < a) ||f(z)||_2 = M(a) < infinity$. Then the Fourier coefficients of $f$ satisfy
  $ |hat(f)_k| lt.eq.slant M(a) thin e^(-a |k|), quad k in ZZ. $ <eq-paley-wiener>
]

The implications for the trapezoidal rule and Fourier pseudospectral methods are direct: the *interpolation* error and the *trapezoidal-rule integration* error are both bounded by the *tail sum* of $|hat(f)_k|$ for $|k| > N \/ 2$, which @eq-paley-wiener bounds by $cal(O)(e^(-a N \/ 2))$. The *strip width* $a$ is the geometric quantity that determines the rate.

Strip width is set by the imaginary distance to the nearest non-real singularity. The Poisson kernel $1 \/ (2 - cos x)$ of Étude 16.4 has poles at $x = plus.minus i log(2 + sqrt(3)) approx plus.minus 1.317 i$, giving $a = log(2 + sqrt(3)) approx 1.317$ and a per-step error reduction of $e^(-a) = 2 - sqrt(3) approx 0.268$, which @ch-periodic-trap §16.4 measures to all displayed digits. The connection between strip width and Bernstein parameter is the *Joukowski bridge* of the next section.

== The Joukowski map <appendix-complex-joukowski>

The Joukowski map (Joukowski 1910)
$ z = phi(w) = (w + 1 \/ w) \/ 2 $ <eq-joukowski>
is the *single* conformal change of variable that links the periodic Fourier theory of the previous section to the bounded-interval Chebyshev theory of @ch-chebyshev. It maps the exterior of the unit disc $|w| > 1$ in the $w$-plane bijectively onto the exterior of the segment $[-1, 1]$ in the $z$-plane, and the unit circle $|w| = 1$ onto $[-1, 1]$ traversed twice.

Three identities make the Joukowski map central to spectral analysis:

+ *Chebyshev as a cosine in disguise.* On the unit circle $w = e^(i t)$, the map @eq-joukowski reduces to $z = cos t$, and the Chebyshev polynomials become cosines: $T_n (z) = cos(n t)$. This is the identity that powers the FFT-based Chebyshev differentiation of @sec-chebfft.

+ *Bernstein ellipses as Joukowski images of circles.* The image of the circle $|w| = rho > 1$ under @eq-joukowski is precisely the Bernstein ellipse $E_rho$ of @eq-bernstein. The polynomial-approximation theory on $[-1, 1]$ thereby pulls back to a Laurent-series theory on the annulus $1 < |w| < rho$.

+ *Strip width $a$ and Bernstein parameter $rho$.* Under the change of variable $w = e^(i t)$ followed by the Joukowski map, the horizontal strip $|op("Im") t| < a$ in the $t$-plane maps to the annulus $e^(-a) < |w| < e^a$, and hence to the region inside the Bernstein ellipse $E_(e^a)$ minus the segment $[-1, 1]$. The *quantitative bridge* is therefore $rho = e^a$, or equivalently $a = log rho$. A strip-analytic function of width $a$ on $RR$ corresponds to a Bernstein-ellipse-analytic function of parameter $rho = e^a$ on $[-1, 1]$, and the two convergence rates $e^(-a N)$ and $rho^(-N)$ are *the same number*.

This last identity is why the body of the volume can switch fluidly between "the Chebyshev coefficients decay as $rho^(-N)$" and "the Fourier coefficients of the periodic version decay as $e^(-a N)$" without further commentary.

== Cauchy's residue theorem and error formulas <appendix-complex-residue>

Spectral-method error formulas almost always take the form of a contour integral, and the *residue theorem* extracts the explicit decay rate. The cleanest example is Trefethen and Weideman's view of the trapezoidal-rule error @TrefethenWeideman2014.

#block(stroke: 0.8pt + rgb(20, 45, 110), radius: 3pt, inset: 10pt,
  fill: rgb(248, 250, 254))[
  *Theorem (Trefethen--Weideman).*  Let $f$ be $2 pi$-periodic, meromorphic in the strip $|op("Im") z| < a$ with finitely many poles $z_1, dots, z_m$ in the strip. The error of the $N$-point trapezoidal rule on $[0, 2 pi]$ is
  $ E_N (f) eq.def integral_0^(2 pi) f(x) dif x - h sum_(j = 1)^N f(j h) = -2 pi i sum_(k = 1)^m frac(op("Res")_(z = z_k) f(z), 1 - e^(- i N (z_k - z_(k - 1)))) + cal(O)(e^(-a N)). $ <eq-tw-error>
]

The leading-order behaviour is a *sum of residues weighted by geometric series*, and each pole contributes a term of magnitude $cal(O)(|e^(-i N z_k)|) = cal(O)(e^(-N op("Im") z_k))$. The *closest* pole to the real axis dominates, and its imaginary distance is precisely the strip width $a$. The remainder term $cal(O)(e^(-a N))$ is exactly the Paley--Wiener bound of @appendix-complex-strip.

The same residue calculus, transferred to the Joukowski $w$-plane, yields the Hunter--Wang style *complex-plane error portrait* of Étude 15.7. Each contour level in the portrait is a curve along which the residue contribution is constant; the contours are nested ellipses for the Bernstein-parameter formulation and circles for the strip-analyticity formulation, related by @eq-joukowski.

== Sokhotski--Plemelj formulas <appendix-complex-plemelj>

For singular integrals along the real line, the controlling identity is the *Sokhotski--Plemelj formula*. Let $phi.alt$ be a sufficiently smooth function on $RR$ with adequate decay at infinity. The boundary values of the Cauchy integral
$ Phi(z) = frac(1, 2 pi i) integral_(-infinity)^(infinity) frac(phi.alt(t), t - z) dif t $ <eq-cauchy-integral>
are
$ Phi(x plus.minus i 0) = plus.minus frac(1, 2) phi.alt(x) + frac(1, 2 pi i) op("PV") integral_(-infinity)^(infinity) frac(phi.alt(t), t - x) dif t. $ <eq-plemelj>

The Hilbert transform is precisely twice the principal-value piece:
$ H{phi.alt}(x) = frac(1, pi) op("PV") integral_(-infinity)^(infinity) frac(phi.alt(t), x - t) dif t. $ <eq-hilbert-via-plemelj>

The Sokhotski--Plemelj formulas are the analytical foundation of the *FFT-multiplier* algorithm for $H$ used in Étude 21.8. On the periodic interval, the transform reduces to multiplication by $-i op("sgn")(k)$ in Fourier space; the derivation goes via @eq-plemelj applied to the Cauchy integral of the periodised target. The same identity appears in the *radiation-basis* construction of Étude 21.5, where the boundary-value behaviour of an outgoing wave at infinity is encoded as a one-sided Cauchy integral.

== The complex-plane detour <appendix-complex-detour>

Some differential operators have *interior* singularities --- coefficient functions with poles or branch points on the real interval over which the eigenproblem is posed. The conventional discretisation then fails: the numerical eigenvalues do not converge under refinement, because the discrete operator has no well-defined real spectrum. The *complex-plane detour*, central to Étude 18.11, is the analytical-and-numerical remedy.

The idea is short. Replace the real interval $[a, b]$ by a contour $Gamma subset.eq CC$ that *avoids* the interior singularity, with the same endpoints $a, b$. Discretise the differential operator on the contour --- the Chebyshev grid points are placed at $x_j = (a + b) \/ 2 + ((b - a) \/ 2) gamma(t_j)$, where $t_j$ are the standard Lobatto points and $gamma$ is a smooth deformation of the real axis. The resulting eigenvalues are *complex*: the operator is intrinsically non-self-adjoint, because the original singular operator was. The eigenvalues converge under refinement of $N$, and the choice of detour amplitude $Delta$ provides a second tuning parameter (small $Delta$: too close to the pole; large $Delta$: too much imaginary distortion --- a familiar V-shape in the convergence study).

The detour is *not* an approximation to a real spectrum; it is a *direct computation* of the complex spectrum the continuous operator has. The eigenfunctions live on the deformed contour, not on the real axis, and reconstructing them at real $y$ requires analytic continuation. The technique is analysed in @Boyd2000 and @ColbrookTownsend2025; a recent set-valued AAA refinement appears in @LietaertEtAl2022.

The complex-plane detour is also the unifying vantage point for several seemingly separate strategies in the body of the volume: the *radiation basis* of Étude 21.5 (which factors out an asymptotic non-decaying mode), the *singularity subtraction* of Étude 21.3 (which removes a known singular term analytically), and the *asymptotic augmentation* of Étude 20.11 (which adds known oscillatory carriers to the basis) are all examples of the same principle: *teach the basis the complex-analytic structure of the target function before asking the basis to approximate it*.
