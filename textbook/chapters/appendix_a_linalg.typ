// textbook/chapters/appendix_a_linalg.typ
// Appendix A: Linear-Algebra Essentials for Spectral Methods
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Last modified: April 2026

#import "../styles/template.typ": dropcap

= Linear-Algebra Essentials for Spectral Methods <appendix-linalg>

#dropcap[The body of this volume cites a small set of dense linear-algebra workhorses --- QR, QZ, the singular value decomposition, the Schur form, condition numbers, Krylov subspace methods --- as black-box machinery available from any modern numerical library.] This appendix collects the essentials in one place. The treatment is deliberately compact: definitions and key properties are stated, but proofs are sketched only where they illuminate the spectral-methods context, and the reader is referred to standard texts (Trefethen and Bau, Golub and Van Loan, Saad, Watkins) for the full development.

== Why this appendix <appendix-linalg-why>

Several places in the body assume the reader is comfortable with the algorithms listed above. @ch-linear-eigen leans on QR (for regular eigenproblems), QZ (for generalised pencils), the homogeneous $(alpha, beta)$ representation, inverse iteration with shift, and the Arnoldi recurrence. @ch-bvp and @ch-spectral-pde build BVP solvers from dense LU factorisations and refer to the singular value decomposition when discussing the achievable accuracy floor. @ch-time-stability invokes matrix exponentials in integrating-factor schemes. @ch-coord-transforms and @ch-higher-order quote condition-number scalings $kappa(D^p) tilde.op N^(2 p)$ that govern how high-order derivative matrices saturate to the rounding floor.

A reader who is fluent with these tools may skip this appendix; one who has met them in passing but never had cause to nail down the details will find the next eight sections enough to read the body of the volume without leaving the page.

== The QR algorithm <appendix-linalg-qr>

The QR algorithm is the standard route to all the eigenvalues of a dense matrix $bold(A) in CC^(n times n)$. The bare iteration is breathtakingly short:

#block(stroke: 0.8pt + rgb(20, 45, 110), radius: 3pt, inset: 10pt,
  fill: rgb(248, 250, 254))[
  *Algorithm (unshifted QR).*  Set $bold(A)_0 = bold(A)$. For $k = 0, 1, 2, dots$:
  - Factor $bold(A)_k = bold(Q)_k bold(R)_k$ (orthogonal $bold(Q)_k$, upper-triangular $bold(R)_k$).
  - Form $bold(A)_(k + 1) = bold(R)_k bold(Q)_k$.

  Under mild hypotheses, $bold(A)_k$ converges to (real or complex) Schur form: an upper-triangular matrix whose diagonal entries are the eigenvalues of $bold(A)$.
]

Three refinements turn this skeleton into the LAPACK workhorse used by every spectral eigensolver in the body of the volume. _Hessenberg reduction_ first conjugates $bold(A)$ by orthogonal similarities to a matrix that is zero below the first sub-diagonal; this is a one-shot $cal(O)(n^3)$ pre-processing step that shrinks every subsequent QR step to $cal(O)(n^2)$. _Implicit shifts_ replace the explicit factorise-and-multiply update with a single bulge chase along the sub-diagonal, eliminating round-off contamination of the orthogonality of $bold(Q)_k$. _Wilkinson's shift_ chooses the eigenvalue of the trailing $2 times 2$ block closer to $bold(A)_k [n, n]$ as the spectral target, giving locally cubic convergence on real symmetric problems and quadratic convergence in general. The implicit-Q theorem guarantees that the bulge-chased iterate is the same matrix one would obtain from the explicit shifted iteration, up to a sign convention on each column.

The total cost of computing all eigenvalues by Hessenberg-shifted QR is $cal(O)(c thin n^3)$ with a small constant $c approx 10$ in practice @GottliebOrszag1977; eigenvectors cost a comparable amount more. @Boyd2000 emphasises that the constant matters: for the moderate $N approx 30$--$200$ that suffices for spectrally-resolved problems, the dense QR pipeline is faster than any iterative alternative even though its asymptotic complexity is steeper.

== The QZ algorithm for generalised eigenproblems <appendix-linalg-qz>

A pencil $bold(A) bold(u) = lambda bold(B) bold(u)$, with $bold(A), bold(B) in CC^(n times n)$ and $bold(B)$ possibly singular, is the natural object that arises when boundary bordering or basis recombination is applied to a differential eigenproblem. The QZ algorithm of Moler and Stewart (1973) is the QR algorithm's generalisation to pencils.

#block(stroke: 0.8pt + rgb(20, 45, 110), radius: 3pt, inset: 10pt,
  fill: rgb(248, 250, 254))[
  *Algorithm (QZ outline).*  Reduce $bold(A)$ to upper Hessenberg and $bold(B)$ to upper triangular by simultaneous orthogonal transformations:
  $ bold(Q)^* bold(A) bold(Z) = bold(H), quad bold(Q)^* bold(B) bold(Z) = bold(T), $
  with $bold(H)$ upper Hessenberg and $bold(T)$ upper triangular. Then chase a bulge through $bold(H)$ while preserving the triangularity of $bold(T)$, until $bold(H)$ becomes (quasi-)triangular.  The eigenvalues are read off as the homogeneous pairs $(alpha_i, beta_i) = (bold(H)[i, i], bold(T)[i, i])$.
]

The *homogeneous representation* $(alpha_i, beta_i)$ is more than book-keeping. When $bold(B)$ has a non-trivial null space --- which is the *generic* case for tau-bordered or row-replaced spectral pencils, where boundary rows of $bold(B)$ are deliberately set to zero --- the conventional eigenvalue $alpha_i \/ beta_i$ becomes infinite for the indices that lie in $ker bold(B)$. A naive `isfinite(alpha/beta)` filter is fooled by QZ rounding into accepting these as huge but finite numbers of magnitude $cal(O)(||bold(A)|| \/ epsilon_("mach"))$, corrupting any subsequent statistic taken over the spectrum. The principled remedy, which @ch-linear-eigen Étude 18.7 turns into a sharp diagnostic, is to test $|beta_i| < tau_("inf") max(|alpha_i|, |beta_i|)$ for some small $tau_("inf")$ (typically $10^(-10)$) and to discard the matching pairs as algebraic infinities. A summary of pencil pathologies and a survey of nonlinear and rational extensions are given in @Tisseur2001 and @GuttelTisseur2017.

== The singular value decomposition <appendix-linalg-svd>

Every matrix $bold(A) in CC^(m times n)$ admits a factorisation
$ bold(A) = bold(U) bold(Sigma) bold(V)^*, $ <eq-svd>
with $bold(U) in CC^(m times m)$ and $bold(V) in CC^(n times n)$ orthogonal and $bold(Sigma) in RR^(m times n)$ diagonal with non-negative entries $sigma_1 gt.eq.slant sigma_2 gt.eq.slant dots gt.eq.slant sigma_(min(m, n)) gt.eq.slant 0$, called the _singular values_ of $bold(A)$.

Two consequences of the SVD shape spectral practice. The smallest singular value is the *distance to singularity*: $sigma_min(bold(A) - z bold(I))$ measures how close to the spectrum of $bold(A)$ the complex number $z$ sits. The level set
$ Lambda_epsilon (bold(A)) = {z in CC : sigma_min(bold(A) - z bold(I)) lt.eq.slant epsilon} $ <eq-app-pseudospectrum>
is the *$epsilon$-pseudospectrum* of $bold(A)$, the central object of @ch-time-stability §17.5 and @ch-higher-order Étude 14.8 @TrefethenEmbree2005. For non-normal operators the pseudospectral contours can extend far beyond the eigenvalues themselves, explaining transient growth in advection-dominated discretisations and the subcritical transition observed in shear-flow stability problems.

The SVD's algorithmic role is more pedestrian: when a discretisation matrix becomes nearly singular (because a boundary row has been overwritten with a row of $bold(D)$ whose entries scale as $cal(O)(N^2)$, or because a high-order operator has $kappa(bold(D)^p) tilde.op N^(2 p)$), the safe linear-system solve uses a truncated SVD pseudo-inverse rather than the QR or LU factorisation that ignores the rank-deficiency.

== Condition number and conditioning <appendix-linalg-cond>

For a non-singular square matrix $bold(A)$, the *spectral condition number* is
$ kappa(bold(A)) = frac(sigma_max(bold(A)), sigma_min(bold(A))) = ||bold(A)|| dot ||bold(A)^(-1)||. $ <eq-condition>
Its operational meaning is sharp: in floating-point arithmetic the relative error in the computed solution of $bold(A) bold(u) = bold(f)$ is bounded by $kappa(bold(A)) epsilon_("mach")$, where $epsilon_("mach")$ is machine epsilon. A discretisation matrix with $kappa = 10^(12)$ caps the achievable relative accuracy at $10^(-4)$ in IEEE 754 double precision, *no matter how exact the discretisation is*.

For Chebyshev pseudospectral matrices the conditioning is famously aggressive. The Chebyshev first-derivative matrix on the Lobatto grid satisfies $kappa(bold(D)) tilde.op N^2$, and squaring gives $kappa(bold(D)^2) tilde.op N^4$. The general scaling is $kappa(bold(D)^p) tilde.op N^(2 p)$, so the biharmonic operator already has $kappa tilde.op N^8$ at $N approx 30$, which means even modest fourth-order BVPs run out of double-precision digits long before they run out of approximation accuracy. The cures are by now classical: the *Heinrichs basis* $(1 - x^2)^k T_j$ that bakes the boundary conditions into the trial functions (@ch-higher-order, @ch-linear-eigen Étude 18.8), the coupled second-order reformulation that replaces $cal(O)(N^(2 p))$ with $cal(O)(N^4)$ at the cost of doubling the system size, and the *ultraspherical* method of @OlverTownsend2013 that produces almost-banded operators with bounded condition number.

== The Schur decomposition <appendix-linalg-schur>

The Schur form is the algorithmic target of QR. Every matrix $bold(A) in CC^(n times n)$ is unitarily similar to an upper-triangular matrix:
$ bold(A) = bold(Q) bold(T) bold(Q)^*, quad bold(T) "upper triangular," quad bold(Q) "unitary." $ <eq-schur>
The diagonal of $bold(T)$ holds the eigenvalues of $bold(A)$. Over the reals the analogous *real Schur form* is block-upper-triangular with $1 times 1$ and $2 times 2$ diagonal blocks; the $2 times 2$ blocks correspond to complex-conjugate eigenvalue pairs.

Two facts make the Schur form the right choice for spectral applications. First, it is *backward stable*: the computed $hat(bold(T))$ is the exact Schur form of $bold(A) + Delta bold(A)$ with $||Delta bold(A)|| = cal(O)(epsilon_("mach") ||bold(A)||)$, so an eigenvalue of magnitude $|lambda|$ is computed to absolute accuracy $cal(O)(epsilon_("mach") ||bold(A)||)$ regardless of the eigenvector conditioning. Second, it preserves invariant subspaces: the leading $k$ columns of $bold(Q)$ span the same subspace as the eigenvectors associated with the leading $k$ eigenvalues, allowing partial-spectrum calculations without recomputing the full decomposition. The QZ generalisation produces a *generalised Schur form* $(bold(S), bold(T))$ for the pencil $(bold(A), bold(B))$ with the same backward-stable, invariant-subspace-preserving properties.

== Krylov methods at a glance <appendix-linalg-krylov>

When $N$ grows beyond a few hundred or when only a handful of eigenpairs is wanted (for example in a parameter sweep), the dense $cal(O)(n^3)$ cost of QR/QZ becomes prohibitive. *Krylov subspace methods* exploit the fact that the action of $bold(A)$ on a vector is often much cheaper than the assembled matrix --- in spectral methods, the FFT-based $cal(O)(N log N)$ application of $bold(D)$ is the canonical example.

Three Krylov methods recur in the body of the volume:

- *The power method.* Iterates $bold(u)_(k + 1) = bold(A) bold(u)_k \/ ||bold(A) bold(u)_k||$. Converges to the eigenvector of largest modulus eigenvalue, with rate $|lambda_2 \/ lambda_1|$. Used in @ch-linear-eigen Étude 18.9 to find a single dominant mode at low cost.

- *Inverse iteration with shift.* Iterates $bold(u)_(k + 1) = (bold(A) - mu bold(I))^(-1) bold(u)_k$ followed by normalisation. Converges to the eigenvector whose eigenvalue is closest to the shift $mu$, with rate $|lambda_("target") - mu| \/ |lambda_("nearest" eq.not "target") - mu|$. The textbook strategy in @ch-linear-eigen §§18.9--18.10: pick a shift close to the target, do a few inverse-iteration steps, and audit the result against a coarse-grid QR anchor.

- *Arnoldi iteration.* Builds an orthonormal basis $bold(Q)_m$ of the Krylov subspace $cal(K)_m (bold(A), bold(b)) = "span"{bold(b), bold(A) bold(b), bold(A)^2 bold(b), dots, bold(A)^(m - 1) bold(b)}$ by modified Gram--Schmidt. The projected matrix $bold(H)_m = bold(Q)_m^* bold(A) bold(Q)_m$ is upper Hessenberg, and its eigenvalues (the *Ritz values*) approximate the extremal eigenvalues of $bold(A)$. For Hermitian $bold(A)$ the Hessenberg reduces to tridiagonal and the iteration is *Lanczos*. Implicit-restart variants (ARPACK, KrylovKit) are the standard tool for large sparse eigenproblems.

A modern complement to direct Krylov is *contour integration*: the FEAST eigensolver of @Polizzi2009 and the Beyn algorithm of @Beyn2012 extract every eigenvalue inside a user-prescribed contour by approximating the spectral projector via numerical quadrature over the contour. These are the methods of choice when the user knows roughly *where* the wanted eigenvalues lie but not how many; they are also embarrassingly parallel across disjoint contours.

== Practical advice for spectral practitioners <appendix-linalg-practical>

The recurring choices a spectral practitioner faces, condensed:

#block(stroke: 0.8pt + rgb(20, 45, 110), radius: 3pt, inset: 10pt,
  fill: rgb(248, 250, 254))[
  *Decision guide for the linear-algebra layer.*

  + *All eigenvalues, $N lt.eq.slant 1000$:* use dense QR (regular) or QZ (generalised). The constant in $cal(O)(n^3)$ is small and library implementations (LAPACK / Intel MKL / OpenBLAS) are already maximally tuned.

  + *Generalised pencil with rank-deficient $bold(B)$:* always extract eigenvalues via the homogeneous $(alpha, beta)$ pair, never via $alpha \/ beta$; filter algebraic infinities with $|beta_i| < tau_("inf") max(|alpha_i|, |beta_i|)$ for $tau_("inf") tilde.op 10^(-10)$.

  + *A few eigenvalues near a target, $N$ large:* use inverse iteration with shift, anchored on a coarse-grid QR computation that supplies starting guesses and audits the result.

  + *Many eigenvalues over a wide spectral region:* use FEAST or contour-integral Beyn over disjoint sub-contours; both are parallel-friendly and avoid the $cal(O)(n^3)$ cost of full QR.

  + *Linear-system solve with a stiff differential operator:* prefer the *ultraspherical* representation @OlverTownsend2013 if available, or a coupled second-order reformulation (@ch-higher-order §14.5), over the naive $bold(D)^p$ that suffers from $cal(O)(N^(2 p))$ conditioning.

  + *Pseudospectra:* the SVD-based formula @eq-app-pseudospectrum is robust but expensive ($cal(O)(n^3)$ per grid point); for resolvent-norm contour plots use the projection-based algorithm of @TrefethenEmbree2005 chapter 39.
]

The single most important habit is to *never trust eigenvalues that have not been cross-validated under refinement*. The drift-with-$N$ diagnostic of @ch-linear-eigen §18.5 is the operational expression of that habit: a backward-stable QR / QZ result is trustworthy as a property of the discretisation matrix, but the discretisation itself can be wrong, and only resolution-doubling can tell.
