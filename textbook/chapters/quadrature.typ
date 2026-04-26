// textbook/chapters/quadrature.typ
// Chapter 15: Quadrature in Spectral Methods: When Exactness Misleads
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Email: denys.dutykh@ku.ac.ae
// Homepage: https://www.denys-dutykh.com/
// Last modified: April 2026

#import "../styles/template.typ": dropcap, num, format-table, etude-conclusion, idx

// Enable equation numbering for this chapter

= Quadrature in Spectral Methods: When Exactness Misleads <ch-quadrature>

#dropcap[Every numerical analysis textbook follows the same script when presenting numerical integration. First comes the Newton--Cotes family, based on equispaced nodes, with polynomial exactness#idx("polynomial exactness") degree $n - 1$. Then comes Gauss quadrature, based on optimal nodes, with the doubled exactness degree $2n - 1$. The implicit message is that Gauss is "twice as good." This chapter is about why that narrative is incomplete, and sometimes deeply misleading. Drawing on two landmark papers by Trefethen#idx("Trefethen") @TrefethenCC2008 @TrefethenExactness2022, we shall see that the _exactness principle_ (the idea that a quadrature formula should be exact for polynomials of a given degree) is a useful tool for _designing_ formulas, but an unreliable guide to their _actual accuracy_. The story will take us from the spectacular failure of Newton--Cotes to the surprising competitiveness of Clenshaw--Curtis, and finally to the paradox of Gauss--Hermite quadrature on unbounded domains. Along the way, we shall encounter aliasing, rational approximation in the complex plane, and the fundamental lesson of spectral methods: _the choice of approximation space matters more than the degree of polynomial exactness_.]

By the end of this chapter, you should be able to:

1. Explain the principle of interpolatory quadrature and distinguish Newton--Cotes, Gauss--Legendre, and Clenshaw--Curtis rules by their node distributions and degrees of precision.
2. Demonstrate the catastrophic failure of Newton--Cotes quadrature for large $n$ and connect it to the Runge phenomenon studied in @ch-geometry.
3. Construct Gauss--Legendre nodes and weights via the Golub--Welsch eigenvalue algorithm and Clenshaw--Curtis weights via the FFT.
4. Reproduce the empirical observation that Clenshaw--Curtis quadrature nearly matches Gauss--Legendre for all smooth functions, despite having half the degree of precision#idx("degree of precision").
5. Explain the aliasing mechanism in Chebyshev space that accounts for the unexpected competitiveness of Clenshaw--Curtis.
6. Visualise quadrature accuracy from the complex-plane perspective, interpreting the quadrature remainder as a contour integral involving rational approximants.
7. Identify situations where Gauss quadrature is suboptimal (unbounded domains, the Gauss--Hermite paradox) and explain why matching the quadrature to the approximation space matters more than maximising the degree of precision.
8. Connect quadrature principles to the spectral methods developed throughout this textbook, recognising the periodic trapezoidal rule#idx("trapezoidal rule") as the quadrature foundation of Fourier spectral methods.

// ============================================================================
== Interpolatory Quadrature: The Exactness Principle <sec-interpolatory>
// ============================================================================

=== The General Framework

We consider the problem of approximating the definite integral
$ I(f) = integral_(-1)^1 f(x) dif x $ <eq-quad-integral>
by a weighted sum of function values
$ I_n (f) = sum_(j=0)^n w_j f(x_j), $ <eq-quad-formula>
where the _nodes_ ${x_j}$ and _weights_ ${w_j}$ depend on $n$ but not on $f$. The formula @eq-quad-formula is called _interpolatory_ if $I_n (f)$ equals the integral of the unique polynomial interpolant of degree $lt.eq.slant n$ through the data ${(x_j, f(x_j))}$. In this case the weights are determined uniquely by
$ w_j = integral_(-1)^1 L_j (x) dif x, quad L_j (x) = product_(k eq.not j) frac(x - x_k, x_j - x_k), $ <eq-quad-weights>
where $L_j$ is the $j$-th Lagrange basis polynomial. By construction, any interpolatory formula with $n + 1$ points integrates polynomials of degree $lt.eq.slant n$ exactly: if $f in PP_n$, then $f$ coincides with its own interpolant, so $I_n (f) = I(f)$.

The _degree of precision_ (or _exactness degree_) of a quadrature formula is the largest integer $d$ such that $I_n (f) = I(f)$ for all $f in PP_d$. For an interpolatory formula on $n + 1$ points, $d gt.eq.slant n$. The question is whether one can do better by a clever choice of nodes.

=== Newton--Cotes, Gauss, and Clenshaw--Curtis Nodes

Three classical families of interpolatory quadrature formulas differ in their node placement:

- *Newton--Cotes*: equispaced nodes $x_j = -1 + 2j\/n$, $j = 0, 1, dots, n$. Exactness degree $n$ (or $n + 1$ when $n$ is even, thanks to symmetry).

- *Gauss--Legendre*: nodes chosen as roots of the Legendre polynomial $P_(n+1)(x)$, maximising the exactness degree to $2n + 1$. This is the highest degree achievable with $n + 1$ points.

- *Clenshaw--Curtis*: Chebyshev extrema $x_j = cos(j pi \/ n)$, $j = 0, 1, dots, n$. Exactness degree $n$ (the same as Newton--Cotes).

@fig-node-visualization shows these three node sets for $n = 32$. The contrast is striking: Newton--Cotes places nodes uniformly, while both Gauss and Clenshaw--Curtis cluster nodes near the endpoints $plus.minus 1$. This clustering, which we studied in @ch-geometry in the context of polynomial interpolation, is essential for avoiding the Runge phenomenon. Indeed, both Gauss and Clenshaw--Curtis nodes share the same asymptotic density $n pi^(-1)(1 - x^2)^(-1\/2)$ as $n arrow infinity$.

#figure(
  image("../figures/ch15/python/node_visualization.pdf", width: 95%),
  caption: [Newton--Cotes (equispaced), Gauss--Legendre, and Clenshaw--Curtis quadrature nodes on $[-1, 1]$ for $n = 32$. The latter two sets cluster near $plus.minus 1$, sharing the same asymptotic distribution.],
) <fig-node-visualization>

The code generating @fig-node-visualization is available in:
- `codes/python/ch15/quad_node_visualization.py`
- `codes/matlab/ch15/quad_node_visualization.m`
- `codes/julia/ch15/quad_node_visualization.jl`

=== Polynomial Exactness and Degree of Precision

From the standpoint of exactness, Gauss quadrature appears to be the clear winner: it integrates polynomials up to degree $2n + 1$ exactly, compared to degree $n$ for both Newton--Cotes and Clenshaw--Curtis. A standard textbook argument suggests that Gauss should therefore be "twice as efficient" as its competitors.

This chapter will demonstrate that this conclusion is profoundly incomplete. The exactness degree tells us about performance on _polynomials_, but real integrands are not polynomials. As we shall see:

- Newton--Cotes can fail _catastrophically_ despite its degree $n$ exactness.
- Clenshaw--Curtis nearly _matches_ Gauss despite having half its exactness degree.
- Even Gauss can be grossly _suboptimal_ when the approximation space is wrong.

The resolution of these paradoxes will lead us to a deeper understanding of quadrature that connects directly to the spectral methods developed throughout this textbook.

// ============================================================================
== Computational Étude 15.1: Node Distributions and Polynomial Exactness <sec-etude-nodes>
// ============================================================================

Our first étude visualises the three node families and tests their polynomial exactness. For each rule with $n + 1 = 33$ points, we compute the quadrature error $|I_n (x^k) - I(x^k)|$ for $k = 0, 1, dots, 2n$. The results confirm the theoretical degrees of precision: Newton--Cotes and Clenshaw--Curtis integrate up to $x^n$ exactly (within rounding), while Gauss integrates up to $x^(2n+1)$ exactly.

In Python, the core computation is:

```python
def newton_cotes_weights(n):
    """Compute Newton-Cotes weights via Vandermonde system."""
    x = np.linspace(-1, 1, n + 1)
    # Integrate Lagrange basis polynomials
    V = np.vander(x, increasing=True).T
    # Right-hand side: integrals of x^k from -1 to 1
    rhs = np.array([(1 - (-1)**(k+1)) / (k + 1) for k in range(n + 1)])
    w = np.linalg.solve(V, rhs)
    return x, w
```

The equivalent MATLAB implementation:

```matlab
function [x, w] = newton_cotes_weights(n)
    x = linspace(-1, 1, n+1)';
    V = zeros(n+1);
    for k = 0:n
        V(k+1,:) = x.^k;
    end
    rhs = arrayfun(@(k) (1 - (-1)^(k+1))/(k+1), (0:n)');
    w = V \ rhs;
end
```

The Julia implementation:

```julia
function newton_cotes_weights(n)
    x = range(-1, 1, length=n+1) |> collect
    V = [x[j]^k for k in 0:n, j in 1:n+1]
    rhs = [(1 - (-1)^(k+1)) / (k + 1) for k in 0:n]
    w = V \ rhs
    return x, w
end
```

The result of the polynomial exactness test is shown in @fig-polynomial-exactness. Two features stand out at a glance. First, both Newton--Cotes and Clenshaw--Curtis sit at the machine-precision floor for $k lt.eq.slant n = 32$ and only depart from it for $k gt.eq.slant 33$, exactly as the degree-of-precision theory predicts; the slow lift-off of the monomial errors past $k = n$ is the same "monomial illusion" that we shall expose more dramatically in @sec-etude-exactness-table. Second, the Gauss--Legendre curve remains pinned to the floor across the entire range $0 lt.eq.slant k lt.eq.slant 64$, since $2n + 1 = 65$ exceeds every degree we test: this is the doubled exactness in action.

#figure(
  image("../figures/ch15/python/polynomial_exactness.pdf", width: 95%),
  caption: [Monomial quadrature errors $|E_n (x^k)| = |I_n (x^k) - I(x^k)|$ for the three rules with $n + 1 = 33$ points, plotted against the monomial degree $k = 0, 1, dots, 2n$. Newton--Cotes (coral) and Clenshaw--Curtis (teal) are exact for $k lt.eq.slant n = 32$ and lift off afterwards; Gauss--Legendre (navy) stays at machine precision throughout because its degree of precision $2n + 1 = 65$ exceeds the entire test range. Odd $k$ are zero by symmetry and have been clamped to the machine-epsilon floor for the log scale.],
) <fig-polynomial-exactness>

#etude-conclusion[
  The figure validates the *degree-of-precision theory* in three rules at once. Newton--Cotes and Clenshaw--Curtis sit at machine precision for $k lt.eq.slant n = 32$ and lift off afterwards (their degree of precision is $n$); Gauss--Legendre stays pinned to the floor across the entire range $0 lt.eq.slant k lt.eq.slant 64$ because $2 n + 1 = 65$ exceeds every test degree --- *doubled exactness in action*. The smooth lift-off of monomial errors past $k = n$ for Newton--Cotes and Clenshaw--Curtis is the start of the "monomial illusion": tested only on monomials, both rules look gentle when in fact one is catastrophic, as Étude 15.3 will expose.
]

The code generating @fig-polynomial-exactness is available in:
- `codes/python/ch15/quad_polynomial_exactness.py`
- `codes/matlab/ch15/quad_polynomial_exactness.m`
- `codes/julia/ch15/quad_polynomial_exactness.jl`

// ============================================================================
== The Failure of Newton--Cotes <sec-newton-cotes>
// ============================================================================

=== Weight Growth and Divergence

The Newton--Cotes formula has a fatal flaw: for large $n$, some weights become negative and grow exponentially in magnitude. The sum of absolute values $sum |w_j|$ grows as $cal(O)(2^n)$, while the sum of the weights themselves remains $sum w_j = 2$ (the length of the interval). This means that for a typical non-polynomial function, the positive and negative contributions to $I_n (f)$ are enormous but nearly cancelling, leaving the result at the mercy of rounding errors and the wild oscillations of equispaced polynomial interpolation.

Pólya proved in 1933 @Polya1933 that convergence $I_n (f) arrow I(f)$ as $n arrow infinity$ for all $f in C([-1, 1])$ occurs if and only if $sum |w_j|$ remains bounded @TrefethenExactness2022. Since Newton--Cotes weights grow exponentially, the formula diverges for most integrands, even analytic ones.

=== The Runge Function Catastrophe

It is worth noting that the failure of Newton--Cotes does not condemn all equispaced quadrature. The much older method of Gregory --- predating both Newton--Cotes and the first descriptions of calculus by Leibniz and Newton --- applies end corrections to the trapezoidal rule on equispaced grids and can be pushed to well above tenth order without any weights becoming negative @Fornberg2025[Sections 7.3 and F.2]. The contemporary revival of end-corrected trapezoidal rules, discussed further in @ch-periodic-trap, demonstrates that equispaced grids are not inherently unsuitable for high-accuracy quadrature; the pathology lies specifically in the uncorrected Newton--Cotes approach.

The connection to the Runge phenomenon studied in @ch-geometry is direct. The Newton--Cotes formula integrates the equispaced polynomial interpolant, and we know from @ch-geometry that this interpolant oscillates wildly near the endpoints for functions like $f(x) = 1\/(1 + 25x^2)$.

Consider the specific case $n = 29$ (thirty equispaced points). The Newton--Cotes formula integrates $x^k$ exactly for $k lt.eq.slant 29$. Yet when applied to the Runge function, it returns $I_(30)(f) approx -21.8$, while the true integral is $I(f) = (2\/5) arctan(5) approx 0.5495$. The sign is wrong, and the magnitude is off by a factor of 40. This is not a gradual accumulation of error; it is a spectacular, qualitative failure.

// ============================================================================
== Computational Étude 15.2: Newton--Cotes on the Runge Function <sec-etude-runge-quad>
// ============================================================================

This étude reproduces the catastrophic failure of Newton--Cotes quadrature. We integrate $f(x) = 1\/(1 + 25x^2)$ using Newton--Cotes, Gauss--Legendre, and Clenshaw--Curtis for increasing $n$, and plot the absolute error on a semilogarithmic scale.

The core computation in Python:

```python
def runge(x):
    """The Runge function."""
    return 1.0 / (1.0 + 25.0 * x**2)

I_exact = (2.0 / 5.0) * np.arctan(5.0)

# Newton-Cotes errors grow exponentially
for n in n_values:
    x_nc, w_nc = newton_cotes_weights(n)
    err_nc = abs(w_nc @ runge(x_nc) - I_exact)
```

The equivalent MATLAB implementation:

```matlab
runge = @(x) 1 ./ (1 + 25*x.^2);
I_exact = (2/5) * atan(5);

for k = 1:length(n_values)
    n = n_values(k);
    [x_nc, w_nc] = newton_cotes_weights(n);
    err_nc(k) = abs(w_nc' * runge(x_nc) - I_exact);
end
```

The Julia implementation:

```julia
runge(x) = 1 / (1 + 25x^2)
I_exact = (2 / 5) * atan(5)

for (k, n) in enumerate(n_values)
    x_nc, w_nc = newton_cotes_weights(n)
    err_nc[k] = abs(dot(w_nc, runge.(x_nc)) - I_exact)
end
```

@fig-newton-cotes-runge shows the result. Newton--Cotes errors grow exponentially (roughly as $2^n$), while both Gauss and Clenshaw--Curtis converge rapidly. This is a vivid demonstration that polynomial exactness degree alone does not determine practical accuracy.

#figure(
  image("../figures/ch15/python/newton_cotes_runge.pdf", width: 85%),
  caption: [Quadrature errors for $f(x) = 1\/(1 + 25x^2)$ on $[-1, 1]$. Newton--Cotes (red) diverges exponentially while Gauss--Legendre (navy) and Clenshaw--Curtis (teal) converge rapidly.],
) <fig-newton-cotes-runge>

#etude-conclusion[
  Newton--Cotes errors *grow exponentially* (roughly as $2^n$) on the Runge function $1 \/ (1 + 25 x^2)$, while Gauss--Legendre and Clenshaw--Curtis converge rapidly. This is a *vivid demonstration that polynomial exactness alone does not determine practical accuracy*: equispaced Newton--Cotes has the same degree of precision as Clenshaw--Curtis, yet diverges. The culprit is weight growth: Newton--Cotes weights become large and alternating in sign, amplifying round-off and the inherent error in the equispaced polynomial interpolant (the Runge phenomenon of @ch-geometry surfacing here as quadrature divergence). The lesson: *use clustered nodes (Chebyshev-Lobatto for Clenshaw--Curtis, Gauss--Legendre roots for Gauss) regardless of degree-of-precision parity*.
]

The code generating @fig-newton-cotes-runge is available in:
- `codes/python/ch15/quad_newton_cotes_runge.py`
- `codes/matlab/ch15/quad_newton_cotes_runge.m`
- `codes/julia/ch15/quad_newton_cotes_runge.jl`

// ============================================================================
== Exactness in Different Bases <sec-exactness-bases>
// ============================================================================

=== Monomials Flatter, Chebyshev Polynomials Reveal

The failure of Newton--Cotes is well known, but what is less appreciated is that the _standard tests_ for quadrature accuracy can be misleading. Table 2.1 of @TrefethenExactness2022 reveals a striking phenomenon: the Newton--Cotes formula with $n = 30$ points not only integrates $x^k$ exactly for $k lt.eq.slant 29$, but the errors for $k = 30, 32, 34, dots$ beyond the exactness range are _extremely small_ (of order $10^(-4)$ to $10^(-7)$). A naive tester might conclude that Newton--Cotes is performing brilliantly.

The illusion vanishes when we switch from monomials to Chebyshev polynomials. For the same formula and the same degrees, the errors $|E_n (T_k)|$ for $k gt.eq.slant 30$ are _enormous_ (hundreds to tens of thousands). The formula that appeared nearly exact for $x^(34)$ produces an error of about 8923 for $T_(34)(x)$.

=== The Lesson: Basis Matters

The resolution is that $x^k$ has a _numerical degree_ that is much smaller than its algebraic degree. Specifically, for any $epsilon > 0$, $x^k$ can be approximated to accuracy $epsilon$ on $[-1, 1]$ by polynomials of degree $O(k^(1\/2))$ @TrefethenExactness2022. Intuitively, $x^k$ is nearly flat on most of $[-1, 1]$ and only varies significantly near $x = plus.minus 1$, so it "looks like" a low-degree polynomial to a quadrature formula. In contrast, $T_k (x) = cos(k arccos x)$ oscillates genuinely at frequency $k$ across the entire interval, making it a much more demanding test integrand.

This observation has a direct connection to spectral methods: the Chebyshev basis reveals the true resolution requirements of a function, while the monomial basis can hide them. In spectral computations, we always work in orthogonal bases (Chebyshev, Legendre, Fourier) precisely because they provide honest measures of approximation quality.

// ============================================================================
== Computational Étude 15.3: The Monomial--Chebyshev Comparison Table <sec-etude-exactness-table>
// ============================================================================

In this étude we reproduce the key observation from Table 2.1 of @TrefethenExactness2022. For $n = 30$ equispaced points, we compute the Newton--Cotes quadrature errors for both $x^k$ and $T_k (x)$:

```python
def compare_monomial_chebyshev(n):
    """Compare Newton-Cotes errors for x^k vs T_k(x)."""
    x, w = newton_cotes_weights(n)
    results = []
    for k in range(n - 4, n + 12, 2):  # even k near n
        # Monomial integral: int_{-1}^{1} x^k dx
        I_mono = 2.0 / (k + 1)  # zero for odd k
        err_mono = abs(w @ (x**k) - I_mono)
        # Chebyshev integral: int_{-1}^{1} T_k(x) dx
        I_cheb = 2.0 / (1.0 - k**2) if k % 2 == 0 else 0.0
        Tk = np.cos(k * np.arccos(x))
        err_cheb = abs(w @ Tk - I_cheb)
        results.append((k, err_mono, err_cheb))
    return results
```

The equivalent MATLAB implementation:

```matlab
function results = compare_monomial_chebyshev(n)
    [x, w] = newton_cotes_weights(n);
    results = [];
    for k = n-4:2:n+10
        I_mono = 2 / (k + 1);
        err_mono = abs(w' * (x.^k) - I_mono);
        I_cheb = 2 / (1 - k^2);
        Tk = cos(k * acos(x));
        err_cheb = abs(w' * Tk - I_cheb);
        results = [results; k, err_mono, err_cheb];
    end
end
```

The Julia implementation:

```julia
function compare_monomial_chebyshev(n)
    x, w = newton_cotes_weights(n)
    results = []
    for k in (n-4):2:(n+10)
        I_mono = 2 / (k + 1)
        err_mono = abs(dot(w, x .^ k) - I_mono)
        I_cheb = iseven(k) ? 2 / (1 - k^2) : 0.0
        Tk = cos.(k * acos.(x))
        err_cheb = abs(dot(w, Tk) - I_cheb)
        push!(results, (k, err_mono, err_cheb))
    end
    return results
end
```

The results are shown in @tab-exactness-comparison. For $k lt.eq.slant 28$, both bases give zero errors (exact integration). For $k gt.eq.slant 30$, the monomial errors are tiny while the Chebyshev errors are catastrophic. This contrast exposes the fundamental inadequacy of monomial tests as a diagnostic for quadrature quality.

#figure(
  block(
    stroke: (top: 1.5pt + rgb("#142D6E"), bottom: 1.5pt + rgb("#142D6E")),
    inset: 0pt,
    {
      show table: format-table(auto, auto, auto)
      table(
        columns: 3,
        align: (center, center, center),
        inset: (x: 1em, y: 0.6em),
        stroke: none,
        table.hline(stroke: 0.75pt + rgb("#142D6E")),
        table.header(
          table.cell(fill: rgb("#142D6E").lighten(85%))[*$k$*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*$|E_n (x^k)|$*],
          table.cell(fill: rgb("#142D6E").lighten(85%))[*$|E_n (T_k)|$*],
        ),
        table.hline(stroke: 0.5pt + luma(180)),
        [26], num[0], num[0],
        [28], num[0], num[0],
        table.hline(stroke: 0.5pt + luma(180)),
        [30], num[7e-7], num[399.5],
        [32], num[5e-6], num[2711.1],
        [34], num[2e-4], num[8923.1],
        [36], num[4e-5], num[18765.9],
        [38], num[9e-5], num[27812.9],
      )
    },
  ),
  caption: [Newton--Cotes quadrature errors for monomials $x^k$ and Chebyshev polynomials $T_k (x)$ with $n = 30$ points. Below the exactness degree ($k lt.eq.slant 29$), both give zero. Above it, monomials appear well-integrated while Chebyshev polynomials expose the catastrophic ill-conditioning. Adapted from @TrefethenExactness2022.],
) <tab-exactness-comparison>

#etude-conclusion[
  The same Newton--Cotes rule that looked competent on monomials looks *catastrophic* on Chebyshev polynomials: at $k = 38$ the monomial error is $9 times 10^(-5)$ while the Chebyshev-polynomial error is $approx 2.8 times 10^4$. The two columns measure *the same operator* applied to two equivalent bases of the polynomial space, so the gap is not a property of the test functions but of the *basis chosen for diagnosis*. Monomials $x^k$ are a notoriously ill-conditioned basis: their values stay $O(1)$ on $[-1, 1]$ but the linear combinations representing typical smooth functions involve large cancelling coefficients. Chebyshev polynomials $T_k$ are the right basis to diagnose quadrature stability --- as Trefethen @TrefethenExactness2022 emphasises, *test exactness in $T_k$, not $x^k$*.
]

// ============================================================================
== Clenshaw--Curtis Quadrature via the FFT <sec-cc-construction>
// ============================================================================

=== Integrating the Chebyshev Interpolant

The _Clenshaw--Curtis formula_, introduced in 1960 @TrefethenCC2008, computes $I_n (f)$ by integrating the degree $n$ polynomial interpolant of $f$ through the Chebyshev extrema
$ x_j = cos(j pi \/ n), quad j = 0, 1, dots, n. $ <eq-cc-nodes>
These are precisely the Chebyshev points we have used throughout this textbook for spectral differentiation (@ch-chebyshev) and boundary value problems (@ch-bvp). The Clenshaw--Curtis formula is therefore the _natural quadrature companion_ of Chebyshev spectral methods.

If $f$ is expanded in Chebyshev series $f(x) = sum'_(j=0)^infinity a_j T_j (x)$ (where the prime indicates the first term is halved), then truncating to degree $n$ and integrating term by term gives
$ I_n (f) = sum_(j = 0, j "even")^n a_j frac(2, 1 - j^2), $ <eq-cc-chebyshev-moments>
since $integral_(-1)^1 T_j (x) dif x = 2\/(1 - j^2)$ for even $j$ and zero for odd $j$. The Chebyshev coefficients $a_j$ can be computed from the function values $f(x_0), dots, f(x_n)$ via the discrete cosine transform (DCT), which is itself a special case of the FFT. The total cost is therefore $cal(O)(n log n)$, compared to $cal(O)(n^2)$ for computing Gauss--Legendre nodes and weights via the Golub--Welsch algorithm @GolubWelsch1969.

=== The FFT-Based Algorithm

The implementation of Clenshaw--Curtis quadrature follows directly from the seven-line MATLAB code given in @TrefethenCC2008. The algorithm proceeds in four steps:

1. Evaluate $f$ at the Chebyshev points @eq-cc-nodes.
2. Compute the Chebyshev coefficients via the FFT (by extending the data to a $2n$-periodic sequence).
3. Integrate term by term using the Chebyshev moments $2\/(1 - j^2)$.
4. Sum the weighted coefficients.

=== Gauss--Legendre via the Golub--Welsch Algorithm

Gauss--Legendre nodes and weights can be computed elegantly using the Golub--Welsch algorithm @GolubWelsch1969. The key insight is that the Gauss nodes are the eigenvalues of the symmetric tridiagonal Jacobi matrix
$ J_n = mat(delim: "[",
  0, beta_1, , ;
  beta_1, 0, beta_2, ;
  , beta_2, dots.down, ;
  , , , 0
) $
where $beta_j = j \/ sqrt(4j^2 - 1)$ are the recurrence coefficients for the Legendre polynomials. The weights are $w_j = 2 v_(j,1)^2$, where $v_(j,1)$ is the first component of the $j$-th eigenvector. Since the matrix is symmetric tridiagonal, the eigenvalue decomposition costs $cal(O)(n^2)$ operations. Modern asymptotic algorithms by Hale and Townsend @HaleTownsend2013 and Bogaert @Bogaert2014 reduce this cost to $cal(O)(n)$ for very large $n$, but for the moderate values of $n$ that suffice in practice, Golub--Welsch remains the standard.

// ============================================================================
== Computational Étude 15.4: Building Quadrature Rules from Scratch <sec-etude-construction>
// ============================================================================

In this étude we implement both the Golub--Welsch algorithm for Gauss--Legendre quadrature and the FFT-based Clenshaw--Curtis algorithm. The goal is for the reader to see that Clenshaw--Curtis is algorithmically native to spectral methods.

The Gauss--Legendre implementation in Python:

```python
def gauss_legendre(n):
    """Gauss-Legendre nodes and weights via Golub-Welsch."""
    j = np.arange(1, n + 1)
    beta = j / np.sqrt(4.0 * j**2 - 1.0)
    # Symmetric tridiagonal Jacobi matrix
    J = np.diag(beta, 1) + np.diag(beta, -1)
    x, V = np.linalg.eigh(J)
    w = 2.0 * V[0, :]**2
    return x, w
```

The Clenshaw--Curtis implementation in Python:

```python
def clenshaw_curtis_fft(n):
    """Clenshaw-Curtis nodes and weights via DCT-I / FFT."""
    x = np.cos(np.pi * np.arange(n + 1) / n)
    # Chebyshev moments: mu_k = int_{-1}^{1} T_k(x) dx
    c = np.zeros(n + 1)
    c[0::2] = 2.0 / (1.0 - np.arange(0, n + 1, 2)**2)
    # DCT-I via FFT of length 2n (mirror the sequence)
    v = np.concatenate([c, c[-2:0:-1]])
    f = np.real(np.fft.fft(v))
    # Extract weights with endpoint halving
    w = f[:n + 1] / n
    w[0] /= 2; w[-1] /= 2
    return x, w
```

The equivalent MATLAB implementation (adapted from @TrefethenCC2008):

```matlab
function [x, w] = gauss_legendre(n)
    % Golub-Welsch algorithm
    j = (1:n)';
    beta = j ./ sqrt(4*j.^2 - 1);
    J = diag(beta, 1) + diag(beta, -1);
    [V, D] = eig(J);
    x = diag(D);
    [x, idx] = sort(x);
    w = 2 * V(1, idx).^2;
end
```

```matlab
function [x, w] = clenshaw_curtis_weights(n)
    % Clenshaw-Curtis via FFT (Trefethen, 2008)
    x = cos(pi*(0:n)'/n);
    c = zeros(n+1, 1);
    c(1:2:end) = 2 ./ (1 - (0:2:n)'.^2);
    f = real(ifft([c; c(n:-1:2)]));
    w = [f(1)/2; f(2:n); f(n+1)/2];
end
```

The Julia implementation:

```julia
function gauss_legendre(n)
    j = 1:n
    beta = j ./ sqrt.(4j .^ 2 .- 1)
    J = diagm(1 => beta, -1 => beta)
    F = eigen(Symmetric(J))
    idx = sortperm(F.values)
    x = F.values[idx]
    w = 2 .* F.vectors[1, idx] .^ 2
    return x, w
end
```

```julia
function clenshaw_curtis_fft(n)
    # Clenshaw-Curtis nodes and weights via DCT-I / FFT.
    x = [cos(j * π / n) for j in 0:n]
    c = zeros(n + 1)
    for k in 0:2:n
        c[k+1] = 2.0 / (1.0 - k^2)
    end
    v = vcat(c, c[n:(-1):2])  # length 2n, mirror
    f = real.(fft(v))
    w = f[1:n+1] / n
    w[1] /= 2
    w[end] /= 2
    return x, w
end
```

@fig-gauss-cc-construction validates both implementations by comparing their errors against reference values for the integral of $e^x$.

#figure(
  image("../figures/ch15/python/gauss_cc_construction.pdf", width: 85%),
  caption: [Validation of our Gauss--Legendre (Golub--Welsch) and Clenshaw--Curtis (FFT) implementations. The error in integrating $e^x$ on $[-1, 1]$ converges to machine precision for both methods, confirming correctness.],
) <fig-gauss-cc-construction>

#etude-conclusion[
  The étude is a *correctness check*: both rules are constructed from scratch (Gauss--Legendre via the Golub--Welsch eigenvalue decomposition of the Jacobi matrix; Clenshaw--Curtis via the FFT-based Chebyshev moments) and both reach machine precision on the smooth integrand $integral_(-1)^1 e^x dif x$. Beyond validation, the étude makes the *cost asymmetry* explicit: Gauss--Legendre needs $cal(O)(n^2)$ via Golub--Welsch (or $cal(O)(n)$ via Hale--Townsend / Bogaert for very large $n$), while Clenshaw--Curtis costs only $cal(O)(n log n)$ via the FFT. Both are spectrally accurate; the FFT path is preferable whenever the Chebyshev grid coexists with other spectral machinery in the calculation.
]

The code generating @fig-gauss-cc-construction is available in:
- `codes/python/ch15/quad_gauss_cc_construction.py`
- `codes/matlab/ch15/quad_gauss_cc_construction.m`
- `codes/julia/ch15/quad_gauss_cc_construction.jl`

// ============================================================================
== The Empirical Surprise: Clenshaw--Curtis Nearly Ties Gauss <sec-cc-vs-gauss>
// ============================================================================

=== Six Benchmark Functions

With both quadrature rules implemented, we can now test them on a range of integrands. Following @TrefethenCC2008, we consider six functions of decreasing regularity:

1. $f_1 (x) = x^(20)$: a polynomial (band-limited in Chebyshev space).
2. $f_2 (x) = e^x$: entire, with growth $cal(O)(e^(|z|))$ in the complex plane.
3. $f_3 (x) = e^(-x^2)$: entire, with growth $cal(O)(e^(|z|^2))$.
4. $f_4 (x) = 1\/(1 + 16x^2)$: analytic on $[-1, 1]$ with poles at $x = plus.minus i\/4$.
5. $f_5 (x) = e^(-1\/x^2)$: $C^infinity$ but not analytic at $x = 0$.
6. $f_6 (x) = |x|^3$: only $C^2$, with a cusp at $x = 0$.

These functions are ordered by regularity. Standard approximation theory predicts that Gauss should outperform Clenshaw--Curtis by a factor related to the doubled exactness degree. The surprise, documented experimentally in @TrefethenCC2008 and proved theoretically in @TrefethenCC2008 (Theorem 5.1), is that this advantage rarely materialises.

=== The Convergence Race

For $f_1 = x^(20)$, Gauss achieves machine precision at $n = 10$ while Clenshaw--Curtis needs $n = 20$: the factor of 2 is visible. For $f_2 = e^x$, both methods converge exponentially, but Gauss is slightly more efficient (needing roughly $n$ versus $2n$ for a given accuracy); the advantage is less than a factor of 2.

For the remaining four functions, the results are striking. For $f_4 = 1\/(1 + 16x^2)$ and the less smooth functions $f_5$ and $f_6$, the two methods are essentially indistinguishable in accuracy. Clenshaw--Curtis quadrature, with its $n$-point exactness, matches the $2n + 1$-point exactness of Gauss for all practical purposes.

As @TrefethenCC2008 summarises: "Gauss quadrature significantly outperforms Clenshaw--Curtis only for functions analytic in a sizable neighbourhood of $[-1, 1]$. For such functions, the convergence of both methods is so fast that the difference hardly matters."

// ============================================================================
== Computational Étude 15.5: The Six-Function Convergence Race <sec-etude-convergence-race>
// ============================================================================

This étude reproduces the landmark Figure 2 of @TrefethenCC2008. We compute the quadrature error $|I - I_n|$ for both Gauss and Clenshaw--Curtis applied to each of the six benchmark functions, for $n$ ranging from 1 to 32 (or higher).

The core loop in Python:

```python
functions = {
    r'$x^{20}$':     (lambda x: x**20, 2/21),
    r'$e^x$':         (lambda x: np.exp(x), np.exp(1) - np.exp(-1)),
    r'$e^{-x^2}$':    (lambda x: np.exp(-x**2), np.sqrt(np.pi)*erf(1)),
    r'$1/(1+16x^2)$': (lambda x: 1/(1+16*x**2), np.arctan(4)/2),
    r'$e^{-1/x^2}$':  (lambda x: np.where(x==0, 0, np.exp(-1/x**2)),
                        I_exact_5),
    r'$|x|^3$':       (lambda x: np.abs(x)**3, 0.5),
}
for n in range(1, n_max + 1):
    xg, wg = gauss_legendre(n)
    xc, wc = clenshaw_curtis_weights(n)
    for name, (f, I_exact) in functions.items():
        err_gauss[name].append(abs(wg @ f(xg) - I_exact))
        err_cc[name].append(abs(wc @ f(xc) - I_exact))
```

The equivalent MATLAB implementation:

```matlab
for n = 1:n_max
    [xg, wg] = gauss_legendre(n);
    [xc, wc] = clenshaw_curtis_weights(n);
    for k = 1:6
        err_gauss(k, n) = abs(wg * f{k}(xg) - I_exact(k));
        err_cc(k, n) = abs(wc' * f{k}(xc) - I_exact(k));
    end
end
```

The Julia implementation:

```julia
for n in 1:n_max
    xg, wg = gauss_legendre(n)
    xc, wc = clenshaw_curtis_weights(n)
    for (k, (f, Ie)) in enumerate(functions)
        err_gauss[k, n] = abs(dot(wg, f.(xg)) - Ie)
        err_cc[k, n] = abs(dot(wc, f.(xc)) - Ie)
    end
end
```

#figure(
  image("../figures/ch15/python/convergence_race.pdf", width: 95%),
  caption: [Convergence of Gauss (filled circles) and Clenshaw--Curtis (open circles) quadrature for six functions on $[-1, 1]$, for $n$ ranging from 1 to 32. The first function is polynomial, the next two are entire, the fourth is analytic, the fifth is $C^infinity$, and the sixth is $C^2$. For the less smooth functions, the two methods are essentially identical. Adapted from @TrefethenCC2008.],
) <fig-convergence-race>

#etude-conclusion[
  The race against six functions of decreasing smoothness exposes the empirical surprise of the chapter: *Clenshaw--Curtis is essentially as accurate as Gauss in practice*, despite having half the formal degree of precision. For analytic and entire functions the curves are virtually indistinguishable; for the $C^infinity$ and $C^2$ functions they remain glued together because both rules' error envelopes are dominated by the same Chebyshev-coefficient decay of the integrand. Gauss only has a clean two-fold edge for *exact polynomial* test cases, which are not what real applications integrate. The sociological consequence: there is *no good reason* to prefer Gauss--Legendre over Clenshaw--Curtis for a generic smooth integrand on $[-1, 1]$ --- the FFT route is faster to construct and integrates equally well.
]

The code generating @fig-convergence-race is available in:
- `codes/python/ch15/quad_convergence_race.py`
- `codes/matlab/ch15/quad_convergence_race.m`
- `codes/julia/ch15/quad_convergence_race.jl`

// ============================================================================
== Aliasing in Chebyshev Space <sec-aliasing>
// ============================================================================

=== The Chebyshev Coefficient Perspective

Why does Clenshaw--Curtis perform so much better than the exactness principle predicts? The first and most intuitive explanation comes from aliasing in Chebyshev expansions, an idea that goes back to O'Hara and Smith in 1968 @OHaraSmith1968 and was sharpened in @TrefethenCC2008 and @TrefethenExactness2022.

Any Lipschitz continuous $f in C[-1, 1]$ has an absolutely convergent Chebyshev series
$ f(x) = sum_(j=0)^infinity a_j T_j (x). $ <eq-cheb-series>
The Clenshaw--Curtis quadrature error can be decomposed term by term:
$ E_n (f) = I(f) - I_n (f) = sum_(j = n, j "even")^infinity a_j E_n (T_j), $ <eq-cc-error-decomposition>
where $E_n (T_j) = I(T_j) - I_n (T_j)$ is the quadrature error for the individual Chebyshev polynomial $T_j$.

The crucial fact, proved in @TrefethenCC2008 (Theorem 5.2), is that on the Chebyshev grid @eq-cc-nodes, the polynomial $T_(n+p)$ takes the same values as $T_(n-p)$:
$ T_(n+p)(x_j) = T_(n-p)(x_j) quad (0 lt.eq.slant j lt.eq.slant n) $ <eq-aliasing-identity>
for any integer $p$ with $0 lt.eq.slant p lt.eq.slant n$. This is the _aliasing_ phenomenon: on the $n$-point Chebyshev grid, $T_(n+p)$ is indistinguishable from $T_(n-p)$. Consequently, Clenshaw--Curtis quadrature gives the _same_ result for both:
$ I_n (T_(n+p)) = I_n (T_(n-p)) = I(T_(n-p)), $
since $T_(n-p)$ has degree $lt.eq.slant n$ and is therefore integrated exactly.

=== Why Aliasing Rescues Clenshaw--Curtis

The error for $T_(n+p)$ is
$ E_n (T_(n+p)) = frac(8p n, n^4 - 2(p^2 + 1)n^2 + (p^2 - 1)^2) quad (n plus.minus p "even"), $ <eq-cc-Tk-error>
which is of order $cal(O)(n^(-2))$ for moderate $p$. As $p$ ranges from 1 to $n$, the errors remain small, only growing to order 1 as $p$ approaches $n$ (i.e., as $k = n + p$ approaches $2n$).

Now consider the contributions to $E_n (f)$ via @eq-cc-error-decomposition. The "shaded" coefficients $a_j$ with $n lt.eq.slant j lt.eq.slant 2n$ are multiplied by these small quadrature errors. For a function that is _not_ analytic (so that the $a_j$ decrease slowly), the dominant terms come from $j approx 2n$ and beyond, where $E_n (T_j)$ is of order 1 but $a_j$ is already decaying. The net effect is that the Clenshaw--Curtis error matches the Gauss rate. For an _analytic_ function, the $a_j$ decay exponentially, and the shaded coefficients do contribute a visible difference; but by then both methods are already so accurate that the practical distinction is negligible.

This is precisely the content of Theorem 5.1 from @TrefethenCC2008: for functions with $k$ bounded derivatives, Clenshaw--Curtis converges at the same algebraic rate as Gauss:
$ |I(f) - I_n (f)| lt.eq.slant frac(32 V, 15 pi k (2n + 1 - k)^k), $ <eq-cc-theorem-5-1>
which is the _same_ bound that holds for Gauss quadrature (with $n$ replaced by $k\/2$).

// ============================================================================
== Computational Étude 15.6: Aliasing and Chebyshev Coefficients <sec-etude-aliasing>
// ============================================================================

This étude visualises the aliasing mechanism. For $n = 30$, we compute and plot:

1. The Clenshaw--Curtis quadrature error $|E_n (T_k)|$ as a function of even $k$, showing the small errors for $n lt.eq.slant k lt.eq.slant 2n$ and the rise to $cal(O)(1)$ near $k = 2n$.
2. For a specific non-analytic function $f(x) = |x|^3$, we plot the Chebyshev coefficients $|a_j|$ and the products $|a_j E_n (T_j)|$, showing how the small quadrature errors in the shaded region combine with slowly decaying coefficients to produce a small total error.

The core computation in Python:

```python
n = 30
# Clenshaw-Curtis errors for T_k
k_vals = np.arange(0, 200, 2)  # even k only
E_Tk = np.zeros(len(k_vals))
for i, k in enumerate(k_vals):
    I_exact_Tk = 2.0 / (1.0 - k**2) if k > 0 else 2.0
    # Compute I_n(T_k) by evaluating at CC nodes
    xc = np.cos(np.pi * np.arange(n + 1) / n)
    Tk_vals = np.cos(k * np.arccos(xc))
    _, wc = clenshaw_curtis_weights(n)
    E_Tk[i] = abs(wc @ Tk_vals - I_exact_Tk)
```

The equivalent MATLAB implementation:

```matlab
n = 30;
k_vals = 0:2:200;
E_Tk = zeros(size(k_vals));
[xc, wc] = clenshaw_curtis_weights(n);
for i = 1:length(k_vals)
    k = k_vals(i);
    I_exact = 2 / (1 - k^2);
    Tk = cos(k * acos(xc));
    E_Tk(i) = abs(wc' * Tk - I_exact);
end
```

The Julia implementation:

```julia
n = 30
k_vals = 0:2:200
E_Tk = zeros(length(k_vals))
xc, wc = clenshaw_curtis_weights(n)
for (i, k) in enumerate(k_vals)
    I_exact = k == 0 ? 2.0 : 2.0 / (1 - k^2)
    Tk = cos.(k * acos.(xc))
    E_Tk[i] = abs(dot(wc, Tk) - I_exact)
end
```

#figure(
  image("../figures/ch15/python/aliasing_chebyshev.pdf", width: 95%),
  caption: [Left: Clenshaw--Curtis quadrature errors $|E_n (T_k)|$ for even $k$ with $n = 30$. The errors are zero for $k lt.eq.slant 28$, small ($cal(O)(n^(-2))$) in the shaded region $30 lt.eq.slant k lt.eq.slant 60$, and of order 1 beyond. Right: for $f(x) = |x|^3$, the products $|a_k E_n (T_k)|$ show that the shaded coefficients do not dominate the total error. Adapted from @TrefethenExactness2022.],
) <fig-aliasing-chebyshev>

#etude-conclusion[
  This étude *resolves the apparent paradox* of Clenshaw--Curtis quadrature. The left panel shows that beyond exactness ($k > n$), Clenshaw--Curtis errors $|E_n (T_k)|$ remain small in a *shaded transitional band* ($n < k lt.eq.slant 2 n$, error $cal(O)(n^(-2))$) before becoming $cal(O)(1)$. The right panel demonstrates that for a typical low-regularity function ($|x|^3$), the products $|a_k thin E_n (T_k)|$ are *not* dominated by the shaded coefficients --- the small-error zone is large enough to absorb the polynomially decaying Chebyshev coefficients. Aliasing in the Chebyshev expansion (errors of order one only above $k = 2 n$) is therefore *forgiving*: practical Clenshaw--Curtis convergence comes from *both* exact integration of low modes and the modest errors of moderately aliased modes.
]

The code generating @fig-aliasing-chebyshev is available in:
- `codes/python/ch15/quad_aliasing_chebyshev.py`
- `codes/matlab/ch15/quad_aliasing_chebyshev.m`
- `codes/julia/ch15/quad_aliasing_chebyshev.jl`

// ============================================================================
== The Complex-Plane Perspective <sec-complex-plane>
// ============================================================================

=== Quadrature as Rational Approximation

A second, more geometric explanation for the near-equivalence of Gauss and Clenshaw--Curtis comes from the complex plane @TrefethenCC2008. If $f$ is analytic in a neighbourhood of $[-1, 1]$, the Cauchy integral formula gives
$ I = integral_(-1)^1 f(x) dif x = frac(1, 2 pi i) integral.cont_Gamma f(z) phi(z) dif z, $ <eq-cauchy-quad>
where $Gamma$ is a contour enclosing $[-1, 1]$ and
$ phi(z) = integral_(-1)^1 frac(dif x, z - x) = log frac(z + 1, z - 1). $ <eq-phi-function>
Similarly, the quadrature approximation can be written as
$ I_n = sum_(k=0)^n w_k f(x_k) = frac(1, 2 pi i) integral.cont_Gamma f(z) r_n (z) dif z, $ <eq-quad-rational>
where $r_n (z) = sum_(k=0)^n w_k \/ (z - x_k)$ is a rational function of type $(n, n+1)$.

The quadrature error is therefore
$ I - I_n = frac(1, 2 pi i) integral.cont_Gamma f(z) (phi(z) - r_n (z)) dif z. $ <eq-quad-error-contour>
The accuracy of the quadrature formula depends on how well $r_n$ approximates $phi$ in the complex plane near $[-1, 1]$.

=== Padé Approximants and Contour Integrals

Trefethen @TrefethenCC2008 proved that Gauss quadrature corresponds to the _Padé approximant_ of $phi(z)$ at $z = infinity$. This means all the interpolation "power" of Gauss is concentrated at the single point infinity: the function $phi(z) - r_n (z)$ has a zero of order $2n + 3$ at infinity, explaining the exactness degree $2n + 1$.

Clenshaw--Curtis, by contrast, corresponds to a _multipoint Padé approximant_: its $n - 2$ interpolation points are spread along a football-shaped curve in the complex plane near $[-1, 1]$. The approximation at infinity has only order $n + 2$ (matching the exactness degree $n$), but the approximation _near the interval_ is nearly as good as for Gauss.

This is visible in @fig-complex-plane: the contours of $|phi(z) - r_n (z)|$ for Gauss and Clenshaw--Curtis are virtually identical near $[-1, 1]$, diverging only in the far field. For functions with singularities near the interval, the contour integral @eq-quad-error-contour is dominated by the near-field behaviour, and the two methods give nearly identical errors.

// ============================================================================
== Computational Étude 15.7: Complex-Plane Error Portraits <sec-etude-complex>
// ============================================================================

This étude reproduces the striking contour plots of @TrefethenCC2008 (Figure 4). For each quadrature formula, we compute $|phi(z) - r_n (z)|$ on a grid of complex $z$ values and display the contour levels.

```python
def quad_error_complex(z, x, w):
    """Compute |phi(z) - r_n(z)| for a quadrature rule."""
    phi = np.log((z + 1) / (z - 1))  # exact kernel
    r_n = np.sum(w / (z - x))         # rational approximant
    return np.abs(phi - r_n)
```

```matlab
function err = quad_error_complex(z, x, w)
    phi = log((z + 1) ./ (z - 1));
    r_n = sum(w ./ (z - x));
    err = abs(phi - r_n);
end
```

```julia
function quad_error_complex(z, x, w)
    phi = log((z + 1) / (z - 1))
    r_n = sum(w ./ (z .- x))
    return abs(phi - r_n)
end
```

#figure(
  image("../figures/ch15/python/complex_plane.pdf", width: 95%),
  caption: [Contour plots of $|phi(z) - r_n (z)|$ in the complex plane for Newton--Cotes, Gauss, and Clenshaw--Curtis with $n = 16$. The innermost contour is $|phi - r_n| = 1$; successive contours are $10^(-1), 10^(-2), dots$ Gauss and Clenshaw--Curtis are virtually identical near $[-1, 1]$, explaining their similar accuracy for functions with nearby singularities. Adapted from @TrefethenCC2008.],
) <fig-complex-plane>

#etude-conclusion[
  The complex-plane portraits *unify* the chapter's themes. The contour level $|phi - r_n| = c$ predicts the convergence rate for any analytic function with its nearest singularity sitting on that contour: closer singularities $arrow$ slower convergence. Newton--Cotes contours are *visibly squeezed* near the real axis (the equispaced node distribution amplifies error near $[-1, 1]$), explaining its weakness; *Gauss and Clenshaw--Curtis are virtually identical* near $[-1, 1]$, explaining the empirical near-tie of Étude 15.5. The portraits provide a single, geometric language in which to read the convergence rate of *any* quadrature rule on *any* analytic integrand.
]

The code generating @fig-complex-plane is available in:
- `codes/python/ch15/quad_complex_plane.py`
- `codes/matlab/ch15/quad_complex_plane.m`
- `codes/julia/ch15/quad_complex_plane.jl`

// ============================================================================
== Beyond Gauss: When "Optimal" Misleads <sec-beyond-gauss>
// ============================================================================

=== The Gauss--Hermite Paradox

We now turn to the most surprising example in the quadrature literature, drawn from Section 5 of @TrefethenExactness2022. The _Gauss--Hermite formula_ is the standard method for integrals of the form
$ I = integral_(-infinity)^infinity e^(-x^2) f(x) dif x, $ <eq-gauss-hermite-integral>
where $f$ is a bounded function. The formula is defined by choosing nodes ${x_j}$ and weights ${w_j}$ so that $I_n (f) = I(f)$ for all polynomials $f$ of degree $lt.eq.slant 2n - 1$. The nodes are the roots of the Hermite polynomial $H_n (x)$, and the formula achieves the maximum possible exactness degree.

Despite this theoretical optimality, Gauss--Hermite quadrature is _terribly inefficient_ in practice. The root cause is that the Gauss--Hermite nodes spread out to $cal(O)(sqrt(n))$ along the real axis. To integrate polynomials of degree $2n - 1$ multiplied by $e^(-x^2)$, the formula must sample the integrand at these distant points. But the weight function $e^(-x^2)$ has already decayed to negligible values there: for $|x| > 6$, we have $e^(-x^2) < 10^(-16)$.

=== Wasted Weights on Unbounded Domains

The consequence is that for large $n$, the vast majority of Gauss--Hermite quadrature weights fall _below machine precision_. For $n = 100$, about 48 of the 100 weights are below $10^(-16)$. For $n = 1000$, the number is 836 out of 1000. These function evaluations are computational dead weight: they contribute nothing to the integral but consume resources @TrefethenExactness2022.

Worse, the convergence is merely _root-exponential_: errors decrease as $cal(O)(e^(-C sqrt(n)))$ for some $C > 0$. Even with $n = 1000$ points, the error for a simple integrand like $f(x) = cos(x^3)$ is no smaller than $10^(-13)$. This empirical observation has since been promoted to a rigorous theorem: Kazashi, Suzuki and Goda @KazashiSuzukiGoda2023 proved matching lower and upper bounds showing that Gauss--Hermite quadrature decays at exactly half the optimal rate ($n^(-alpha\/2)$ instead of $n^(-alpha)$) in weighted Sobolev spaces of order $alpha$, and the bottleneck is the $1\/sqrt(n)$ node spacing rather than any choice of weights.

// ============================================================================
== Computational Étude 15.8: Gauss--Hermite Wasted Weights <sec-etude-hermite-weights>
// ============================================================================

This étude visualises the inefficiency of Gauss--Hermite quadrature by plotting the weights as a function of node location.

```python
from scipy.special import roots_hermite

def plot_hermite_weights(n):
    """Visualise Gauss-Hermite weight distribution."""
    x, w = roots_hermite(n)
    # Count weights below machine epsilon
    n_wasted = np.sum(w < 1e-16)
    print(f"n={n}: {n_wasted}/{n} weights below machine eps")
    return x, w, n_wasted
```

```matlab
function [x, w, n_wasted] = plot_hermite_weights(n)
    [x, w] = hermpts(n);  % requires Chebfun, or use custom GW
    n_wasted = sum(w < 1e-16);
    fprintf('n=%d: %d/%d weights below machine eps\n', n, n_wasted, n);
end
```

```julia
using FastGaussQuadrature
function plot_hermite_weights(n)
    x, w = gausshermite(n)
    n_wasted = count(w .< 1e-16)
    println("n=$n: $n_wasted/$n weights below machine eps")
    return x, w, n_wasted
end
```

#figure(
  image("../figures/ch15/python/gauss_hermite_weights.pdf", width: 95%),
  caption: [Left: Gauss--Hermite quadrature weights for $n = 100$. The dashed line marks machine precision ($10^(-16)$). About half of the weights lie below this threshold. Right: fraction of weights below machine precision as a function of $n$, approaching 100% as $n arrow infinity$. Adapted from @TrefethenExactness2022.],
) <fig-gauss-hermite-weights>

#etude-conclusion[
  About *half the Gauss--Hermite weights are below machine precision* at $n = 100$, and the fraction approaches $100 %$ as $n arrow infinity$. These weights are not "very small but useful" --- they are *zero in floating-point arithmetic*, and the corresponding nodes are wasted: the rule formally has $n$ nodes but contributes meaningfully only with $cal(O)(sqrt(n))$ of them. This is the core of the warning about "optimal" rules: optimality of the *real-line-with-Gaussian-weight* problem is achieved at the cost of a node distribution whose tails exceed the dynamic range of double precision. For the $integral_(-infinity)^infinity e^(-x^2) f(x) dif x$ integral, *truncated trapezoidal* (Étude 15.9) far outperforms Gauss--Hermite at the same arithmetic budget --- a striking inversion of textbook expectations.
]

The code generating @fig-gauss-hermite-weights is available in:
- `codes/python/ch15/quad_gauss_hermite_weights.py`
- `codes/matlab/ch15/quad_gauss_hermite_weights.m`
- `codes/julia/ch15/quad_gauss_hermite_weights.jl`

// ============================================================================
== Computational Étude 15.9: Gauss--Hermite vs. Truncation <sec-etude-hermite-truncation>
// ============================================================================

The most powerful demonstration of the Gauss--Hermite paradox comes from comparing it with a simple alternative: truncate the real line to a finite interval $[-L, L]$ and apply an ordinary (unweighted) quadrature rule there. We take $f(x) = cos(x^3)$ and compare:

- Gauss--Hermite quadrature on $(-infinity, infinity)$ with $n$ points.
- The periodic trapezoidal rule on $[-6, 6]$ with $n$ points (incorporating $e^(-x^2)$ into the integrand).

```python
from scipy.special import roots_hermite

def gauss_hermite_vs_truncation(n_values, L=6.0):
    """Compare Gauss-Hermite with truncated trapezoidal."""
    f = lambda x: np.cos(x**3)
    I_ref = 1.2254167024651776  # high-precision reference
    err_gh, err_trap = [], []
    for n in n_values:
        # Gauss-Hermite
        x_gh, w_gh = roots_hermite(n)
        err_gh.append(abs(w_gh @ f(x_gh) - I_ref))
        # Trapezoidal on [-L, L] with weight e^{-x^2}
        h = 2 * L / n
        x_tr = np.linspace(-L, L, n + 1)
        g = np.exp(-x_tr**2) * f(x_tr)
        err_trap.append(abs(h * (0.5*g[0] + g[1:-1].sum()
                         + 0.5*g[-1]) - I_ref))
    return err_gh, err_trap
```

```matlab
function [err_gh, err_trap] = gauss_hermite_vs_truncation(n_values, L)
    f = @(x) cos(x.^3);
    I_ref = 1.2254167024651776;
    for k = 1:length(n_values)
        n = n_values(k);
        [x_gh, w_gh] = hermpts(n);
        err_gh(k) = abs(w_gh * f(x_gh) - I_ref);
        h = 2*L / n;
        x_tr = linspace(-L, L, n+1)';
        g = exp(-x_tr.^2) .* f(x_tr);
        err_trap(k) = abs(h * (g(1)/2 + sum(g(2:end-1)) + g(end)/2) - I_ref);
    end
end
```

```julia
function gauss_hermite_vs_truncation(n_values; L=6.0)
    f(x) = cos(x^3)
    I_ref = 1.2254167024651776
    err_gh = Float64[]
    err_trap = Float64[]
    for n in n_values
        x_gh, w_gh = gausshermite(n)
        push!(err_gh, abs(dot(w_gh, f.(x_gh)) - I_ref))
        h = 2L / n
        x_tr = range(-L, L, length=n+1)
        g = exp.(-x_tr .^ 2) .* f.(x_tr)
        push!(err_trap, abs(h * (g[1]/2 + sum(g[2:end-1]) + g[end]/2) - I_ref))
    end
    return err_gh, err_trap
end
```

#figure(
  image("../figures/ch15/python/gauss_hermite_failure.pdf", width: 85%),
  caption: [Quadrature errors for $integral_(-infinity)^infinity e^(-x^2) cos(x^3) dif x$. Gauss--Hermite (navy) converges root-exponentially ($cal(O)(e^(-C sqrt(n)))$, parabola on semilog scale), while the trapezoidal rule on $[-6, 6]$ (coral) converges much faster, reaching machine precision around $n = 80$. Adapted from @TrefethenExactness2022.],
) <fig-gauss-hermite-failure>

#etude-conclusion[
  *Gauss--Hermite converges root-exponentially* ($cal(O)(e^(-C sqrt(n)))$, plotting as a parabola on a semilog scale), while the truncated trapezoidal rule on $[-6, 6]$ achieves *machine precision around $n = 80$*. The gap is the chapter's most surprising experimental fact: a "naive" trapezoidal sum on a wide-enough finite interval beats the formally optimal Gauss--Hermite rule for an integral against $e^(-x^2)$, by orders of magnitude. The reason: $e^(-x^2)$ kills the integrand outside $[-6, 6]$ to *machine precision*, turning the real-line problem into an effectively periodic one where the trapezoidal rule of @ch-periodic-trap converges geometrically. *Optimality is local to a problem class*, and Gauss--Hermite's optimality on polynomial-times-Gaussian integrands does not extend to oscillatory factors like $cos(x^3)$.
]

The code generating @fig-gauss-hermite-failure is available in:
- `codes/python/ch15/quad_gauss_hermite_failure.py`
- `codes/matlab/ch15/quad_gauss_hermite_failure.m`
- `codes/julia/ch15/quad_gauss_hermite_failure.jl`

// ============================================================================
== Approximation Spaces and True Convergence Rates <sec-approx-spaces>
// ============================================================================

=== Theorem 5.1: The Truncation Advantage

The Gauss--Hermite paradox has a precise mathematical explanation, given by Theorem 5.1 of @TrefethenExactness2022. The theorem states:

#block(
  fill: rgb("#142D6E").lighten(92%),
  stroke: (left: 3pt + rgb("#142D6E")),
  inset: (left: 12pt, y: 10pt, right: 10pt),
  width: 100%,
)[
*Theorem* (Trefethen, 2022). Let $f(x)$ be analytic and bounded for $x in (-infinity, infinity)$, and suppose $e^(-x^2)f(x)$ extends to a bounded analytic function in the infinite strip $-a < "Im" x < a$ for some $a > 0$. Let $L > 0$ be fixed, and for each $n gt.eq.slant 1$, let $I_n$ be the estimate of @eq-gauss-hermite-integral obtained by applying Gauss--Legendre, Clenshaw--Curtis, or equispaced trapezoidal quadrature on the truncated interval $[-L n^(1\/3), L n^(1\/3)]$. Then for some $C > 0$,
$ |I - I_n| = cal(O)(e^(-C n^(2\/3))), quad n arrow infinity. $
]

Compare this with the Gauss--Hermite rate of $cal(O)(e^(-C n^(1\/2)))$. The truncated method is _exponentially faster_. The exponent $2\/3$ versus $1\/2$ arises because:

- Gauss--Hermite uses polynomials of degree $n - 1$ on an interval of width $cal(O)(n^(1\/2))$, giving a resolution per unit length of $cal(O)(n^(1\/2))$.
- The truncated method uses $n$ points on an interval of width $cal(O)(n^(1\/3))$, giving a resolution per unit length of $cal(O)(n^(2\/3))$.

The wider interval forced by Gauss--Hermite wastes approximation power on regions where the integrand is negligible.

=== Matching Quadrature to the Function Space

This is the deepest spectral methods lesson in both papers. The issue is not the quadrature formula _per se_ but the _approximation space_ it implicitly uses. Gauss--Hermite implicitly approximates $f$ by weighted polynomials $p(x) e^(-x^2)$, and these polynomials must control their behaviour on the entire interval $[-cal(O)(n^(1\/2)), cal(O)(n^(1\/2))]$. If instead one truncates to $[-cal(O)(n^(1\/3)), cal(O)(n^(1\/3))]$ and uses unweighted polynomials (or trigonometric polynomials), the approximation space becomes much more efficient because it only needs to represent $f$ where the weight function is non-negligible.

The moral: *the right approximation space matters more than the nominal degree of polynomial exactness*. This is precisely the principle that underlies the entire spectral methods enterprise. We choose Chebyshev expansions over monomial expansions, Fourier expansions for periodic problems, and spherical harmonics for problems on the sphere, not because of any exactness principle, but because these bases provide efficient representations of the functions we actually encounter.

// ============================================================================
== Computational Étude 15.10: Experimental Convergence Rates <sec-etude-rates>
// ============================================================================

In this final étude, we verify Theorem 5.1 experimentally. For $f(x) = cos(x^3)$, we plot the errors of Gauss--Hermite and truncated Gauss--Legendre on a log scale against both $n^(1\/2)$ and $n^(2\/3)$. A straight line on the $n^(1\/2)$ plot for Gauss--Hermite and on the $n^(2\/3)$ plot for the truncated method confirms the predicted convergence rates.

```python
n_vals = np.arange(10, 501, 10)
# Compute errors for both methods
err_gh, err_trunc = gauss_hermite_vs_truncation(n_vals, L=6.0)
# Plot against n^{1/2} and n^{2/3}
fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.semilogy(np.sqrt(n_vals), err_gh, 'o')
ax1.set_xlabel(r'$n^{1/2}$')
ax2.semilogy(n_vals**(2/3), err_trunc, 's')
ax2.set_xlabel(r'$n^{2/3}$')
```

```matlab
n_vals = 10:10:500;
[err_gh, err_trunc] = gauss_hermite_vs_truncation(n_vals, 6);
subplot(1,2,1); semilogy(sqrt(n_vals), err_gh, 'o');
xlabel('n^{1/2}');
subplot(1,2,2); semilogy(n_vals.^(2/3), err_trunc, 's');
xlabel('n^{2/3}');
```

```julia
n_vals = 10:10:500
err_gh, err_trunc = gauss_hermite_vs_truncation(n_vals; L=6.0)
p1 = scatter(sqrt.(n_vals), err_gh, yscale=:log10)
xlabel!(p1, L"n^{1/2}")
p2 = scatter(n_vals .^ (2/3), err_trunc, yscale=:log10)
xlabel!(p2, L"n^{2/3}")
```

#figure(
  image("../figures/ch15/python/convergence_rates.pdf", width: 95%),
  caption: [Verification of the asymptotic convergence rates. Left: Gauss--Hermite error plotted against $n^(1\/2)$ shows a straight line on a semilog scale, confirming $cal(O)(e^(-C n^(1\/2)))$. Right: truncated Gauss--Legendre error plotted against $n^(2\/3)$ shows a straight line, confirming the faster rate $cal(O)(e^(-C n^(2\/3)))$ of Theorem 5.1.],
) <fig-quad-convergence-rates>

*How to read these figures.* Both panels are semilogarithmic, with the vertical axis showing $log_10 |I - I_n|$. The horizontal axes, however, are not the raw number of nodes $n$ but the *rescaled* variables $n^(1\/2)$ (left) and $n^(2\/3)$ (right). This rescaling is the whole point of the experiment: an estimate of the form $|I - I_n| approx A thin e^(-C n^alpha)$ becomes
$
  log_10 |I - I_n| approx log_10 A - (C / (ln 10)) thin n^alpha,
$
which is a *straight line* when plotted against $n^alpha$, but a curve against $n$. Seeing a clean line in each panel is therefore the visual signature that we have guessed the correct exponent $alpha$: $alpha = 1\/2$ for Gauss--Hermite and $alpha = 2\/3$ for the truncated Gauss--Legendre rule. Had we plotted Gauss--Hermite against $n^(2\/3)$, the points would bend upward (the rule is slower than that scaling); plotted against $n^(1\/3)$, they would bend downward.

*Interpreting the slopes.* The dashed red lines are least-squares fits to the asymptotic part of each curve, and their slopes give a direct estimate of the constant $C$ in the exponent. With base-10 logarithms,
$
  C = -("slope") times ln 10 approx -("slope") times 2.303 thin .
$
The left panel reports a slope of about $-0.42$, hence $C_("GH") approx 0.97$, so the Gauss--Hermite error decays roughly like $e^(-0.97 thin sqrt(n))$: doubling $n$ improves the error by a factor close to $e^(-0.97 (sqrt(2) - 1) sqrt(n))$, which at $n = 100$ is only about one decade. The right panel has a much shallower slope, near $-0.10$, but the exponent is now $n^(2\/3)$ rather than $n^(1\/2)$, so the constant $C_("GL") approx 0.23$ multiplies a much larger function of $n$. Already at $n approx 60$ the truncated rule has overtaken Gauss--Hermite by several decades, and the gap widens monotonically: the steeper exponent $alpha = 2\/3$ always wins over the larger constant in front of $sqrt(n)$. Finally, the curves flatten and become noisy near $10^(-15)$; this is the machine-precision floor $epsilon_("mach")$, not a failure of the theory, and the fits are deliberately restricted to the region above it.

#etude-conclusion[
  The étude turns the abstract convergence-rate predictions into *measured slopes on rescaled axes*. Plotting $log E$ against $n^(1 \/ 2)$ (Gauss--Hermite) or $n^(2 \/ 3)$ (truncated Gauss--Legendre) produces *straight lines*, confirming the predicted exponents. The fitted slopes give the constants in $E approx A thin e^(-C n^alpha)$ directly: $C_("GH") approx 0.97$, $C_("GL") approx 0.23$. Despite the smaller constant, the steeper exponent $alpha = 2 \/ 3$ wins: by $n approx 60$ the truncated rule has overtaken Gauss--Hermite by several decades. *Steeper exponent always wins eventually over a larger constant in front of $sqrt(n)$* --- a structural lesson that recurs whenever subgeometric rates are compared.
]

The code generating @fig-quad-convergence-rates is available in:
- `codes/python/ch15/quad_convergence_rates.py`
- `codes/matlab/ch15/quad_convergence_rates.m`
- `codes/julia/ch15/quad_convergence_rates.jl`

// ============================================================================
== Bridge to Spectral Methods <sec-quad-bridge>
// ============================================================================

=== The Periodic Trapezoidal Rule and Trigonometric Exactness

We close this chapter by returning to Fourier spectral methods, studied in @ch-fourier-grids and @ch-fourier-pseudo. The periodic trapezoidal rule on $[0, 2pi)$,
$ I_n (f) = frac(2 pi, n) sum_(j=0)^(n-1) f(x_j), quad x_j = frac(2 pi j, n), $ <eq-periodic-trap>
is _exact_ for trigonometric polynomials of degree $lt.eq.slant n\/2 - 1$. From the standpoint of piecewise-linear interpolation (the usual way to think about the trapezoidal rule), the formula has exactness degree 1, which is unimpressive. But the correct exactness space is not piecewise-linear functions; it is _trigonometric_ polynomials.

For a smooth, periodic function $f$, the Fourier coefficients $hat(f)_k$ decay rapidly, and the trapezoidal rule converges _exponentially_: if $f$ is analytic in a strip $|"Im"(x)| < a$, the error is $cal(O)(e^(-a n))$. This exponential convergence, which drives the success of Fourier pseudospectral methods, would be completely invisible if one judged the trapezoidal rule only by its piecewise-linear exactness.

=== Choose the Space First, Then the Quadrature

The examples in this chapter all point to the same conclusion, which is the central lesson of both @TrefethenCC2008 and @TrefethenExactness2022:

#block(
  fill: rgb("#142D6E").lighten(92%),
  stroke: (left: 3pt + rgb("#142D6E")),
  inset: (left: 12pt, y: 10pt, right: 10pt),
  width: 100%,
)[
*The Grand Moral.* Polynomial exactness is a useful _design_ principle for deriving quadrature formulas, but it is an unreliable _predictor_ of practical accuracy. What matters is the _approximation space_: the collection of basis functions in which the integrand is efficiently represented. In spectral methods, the right choice is always dictated by the structure of the problem: Chebyshev for non-periodic problems on bounded intervals, Fourier for periodic problems, and problem-specific bases for unbounded domains.
]

This principle explains every paradox we have encountered:

- *Newton--Cotes fails* because equispaced polynomial interpolation is the wrong approximation space for smooth functions on $[-1, 1]$ (the Runge phenomenon).
- *Clenshaw--Curtis matches Gauss* because both use Chebyshev-type nodes, and the aliasing structure of the Chebyshev basis compensates for the lower exactness degree.
- *Gauss--Hermite is suboptimal* because weighted polynomial approximation on $(-infinity, infinity)$ is the wrong approximation space for functions that are effectively supported on a finite interval.
- *The periodic trapezoidal rule is exponentially accurate* because trigonometric polynomials are the right approximation space for periodic functions.

// ============================================================================
== A Non-Exhaustive Literature Overview <sec-quad-literature>
// ============================================================================

The history of numerical quadrature spans three centuries, but the questions it raises are sharper and more topical today than at any moment since the work of Gauss. The story has a clear arc: from the early triumph of exactness as an algebraic design principle, through its experimental humiliation on simple functions like Runge's, to the modern realisation that exactness is at best a heuristic and at worst an actively misleading metric. The following overview traces this arc from the foundational theorems of the early twentieth century to the rational-approximation and approximation-space programmes that define the field in the mid-2020s.

=== From Newton--Cotes to Pólya: The Convergence Question

Interpolatory quadrature originates in the work of Newton, Cotes, Simpson, and Bode in the seventeenth and eighteenth centuries. These rules were designed to be exact for polynomials of as high a degree as the number of nodes would permit, and for two centuries this exactness was treated as proof of accuracy. The first crack in this picture came in 1933, when Pólya @Polya1933 proved that an interpolatory rule converges as $n arrow infinity$ for every continuous function on $[-1, 1]$ if and only if the sums $sum_(j=0)^n |w_j|$ remain uniformly bounded. Since the standard sum of weights satisfies $sum w_j = 2$, any quadrature rule with strictly positive weights is unconditionally convergent for continuous integrands. Newton--Cotes, however, develops alternating-sign weights for $n gt.eq.slant 8$, and the absolute sum diverges as $cal(O)(2^n \/ (n log n))$. Pólya's theorem therefore predicted, decades before computers could observe it, that high-order Newton--Cotes was doomed to diverge. The history of interpolatory quadrature in the early electronic age, including the practical reaction to Pólya's diagnosis, is recounted in detail in @TrefethenExactness2022.

The same period saw the construction of the Gauss--Legendre rule and its extensions to other weight functions. The classical computational bottleneck in Gauss quadrature--computing the nodes and weights from the Jacobi matrix--was elegantly resolved by Golub and Welsch @GolubWelsch1969 via a tridiagonal eigenvalue computation in $cal(O)(n^2)$ operations. Half a century later, Hale and Townsend @HaleTownsend2013 and Bogaert @Bogaert2014 derived asymptotic expansions of the Legendre polynomial roots that drive the cost down to $cal(O)(n)$ for very large $n$, making Gauss--Legendre rules with millions of points routine in modern spectral codes.

=== The Two Trefethen Papers: Aliasing and the End of Exactness

The empirical refutation of polynomial exactness as a predictor of practical accuracy crystallised in two landmark papers by Trefethen. The first, _Is Gauss Quadrature Better than Clenshaw--Curtis?_ @TrefethenCC2008, showed by direct experiment that the doubled exactness degree of Gauss--Legendre confers no measurable advantage over Clenshaw--Curtis for any function that is not analytic in a sizable Bernstein neighbourhood of $[-1, 1]$. Two complementary explanations were given. The first is the aliasing decomposition originally noted by O'Hara and Smith @OHaraSmith1968: on the $(n+1)$-point Chebyshev grid, the polynomials $T_(n+p)$ and $T_(n-p)$ are indistinguishable, so the Clenshaw--Curtis quadrature error for the "shaded" coefficients $a_j$ with $n lt.eq.slant j lt.eq.slant 2n$ is suppressed by an explicit $cal(O)(n^(-2))$ factor. The second is the complex-plane perspective: Gauss quadrature corresponds to the Padé approximant of $log((z+1)\/(z-1))$ at $z = infinity$, while Clenshaw--Curtis corresponds to a multipoint Padé approximant on a "football-shaped" elliptic contour around $[-1, 1]$, and the two contour integrals are virtually identical in the near-field region that controls the integration error of analytic functions with nearby singularities.

The second paper, _Exactness of Quadrature Formulas_ @TrefethenExactness2022, broadened the indictment. Trefethen showed that monomial tests of the form $|E_n (x^k)|$ are systematically misleading because $x^k$ has _numerical degree_ only $cal(O)(k^(1\/2))$: it is essentially flat on $[-1, 1]$ except in a thin layer near $x = plus.minus 1$, and so a low-order quadrature rule processes it as if it were a low-degree polynomial. Replacing $x^k$ by $T_k$ exposes the ill-conditioning immediately, with errors several orders of magnitude larger for the same rule. This basis-dependent diagnostic gap was the first item of evidence; the second was the Gauss--Hermite paradox on the unbounded domain $(-infinity, infinity)$, which showed that even an "optimal" exactness degree can be globally inefficient when most nodes fall in regions where the weight function is below machine precision. The 2022 paper articulates the moral that has since reshaped the field: spectral practitioners must choose the approximation space first and let the quadrature follow.

This empirical critique of Gauss--Hermite was promoted to a sharp theorem by Kazashi, Suzuki and Goda @KazashiSuzukiGoda2023. They proved matching upper and lower bounds for Gauss--Hermite quadrature in weighted Sobolev spaces of order $alpha$, showing that the error decays at exactly half the optimal rate, $cal(O)(n^(-alpha\/2))$ rather than the best-possible $cal(O)(n^(-alpha))$. Crucially, the obstruction is geometric: it is the $1\/sqrt(n)$ minimum spacing between consecutive Gauss--Hermite nodes, not any choice of weights. By constructing explicit "fooling functions", the authors showed that no modification of the weights, however clever, can rescue the rule from this rate. The same paper established that a suitably truncated trapezoidal rule (with $T tilde sqrt(log n)$) achieves the optimal Sobolev rate.

A parallel line of work has used the same approximation-space lens to relax exactness in hyperinterpolation: An and Wu @AnWu2022 showed that the classical hyperinterpolation framework remains stable and convergent even when the underlying quadrature rule is _not_ exact for the discrete inner products, opening the door to a much wider family of admissible quadrature rules on spheres and other compact manifolds.

=== The Rational Revolution: AAA, the Cauchy Transform, and Bespoke Quadrature

The most dramatic development of the past decade is the maturation of rational approximation as a practical tool, and its application to quadrature. Historically, rational approximation (Padé, Remez) was avoided because of the notorious instability of pole-zero arithmetic and the appearance of spurious "Froissart doublets". The Adaptive Antoulas--Anderson (AAA) algorithm of Nakatsukasa, Sète and Trefethen @NakatsukasaSeteTrefethen2018 cured this defect at a stroke by working with a barycentric representation that bypasses the quotient $p(z)\/q(z)$ entirely. AAA combines a greedy choice of support points with a small linear least-squares solve at each step, producing rational approximants of essentially machine-precision accuracy for arbitrary functions on arbitrary point sets in the complex plane. A continuum version of the algorithm was developed by Driscoll, Nakatsukasa, and Trefethen @DriscollNakatsukasaTrefethen2024, removing the dependence on a discrete sample set.

The link between rational approximation and quadrature was forged in Horning and Trefethen @HorningTrefethen2025. They showed that for any integrable weight function $w(z)$ on a Jordan arc $gamma$, applying AAA to the Cauchy transform $C(s) = (2 pi i)^(-1) integral_gamma w(z)\/(s - z) dif z$ on an enclosing contour $Gamma$ produces a rational approximant whose poles are exactly the desired quadrature nodes and whose residues are the corresponding weights. This single mechanism recovers Gauss--Legendre as a special case (the AAA poles for $w equiv 1$ on $[-1, 1]$ collapse onto the zeros of the Legendre polynomials) and generates highly efficient bespoke rules for non-classical weights, weights with algebraic singularities, weights supported on disjoint intervals, and even integrands analytic in $epsilon$-neighbourhoods of $[-1, 1]$. For the latter class, AAA can be optimised to converge a factor $pi\/2$ faster than Gauss while producing more uniformly distributed nodes that avoid the CFL stiffness of strongly clustered rules. Polynomial exactness plays no role anywhere in the construction.

=== Quadrature by Expansion and the Layer-Potential Frontier

A second domain in which exactness has been conclusively superseded is the evaluation of layer potentials in boundary integral equations (BIEs). The integrands here are weakly singular, singular, or hypersingular, and the smoothness assumptions underpinning Gauss--Legendre exactness bounds simply do not hold near the singular point. Klöckner, Barnett, Greengard and O'Neil @KlocknerEtAl2013 introduced _Quadrature by Expansion_ (QBX), which sidesteps the singularity entirely by evaluating the layer potential at a series of expansion centres displaced a safe distance from the boundary, fitting a local harmonic (Taylor or spherical-harmonic) expansion at each centre, and then evaluating the harmonic expansion analytically back at the singular target. Because the off-boundary potential is smooth, standard high-order rules apply unchanged, and QBX integrates seamlessly with the Fast Multipole Method. The QBX approach is the prototype of an approximation-space substitution: rather than fight the singularity in the original space, replace the space with one in which the integrand is smooth.

=== Sparse Grids and the Curse of Dimensionality

In high-dimensional integration--ubiquitous in uncertainty quantification, stochastic Galerkin methods, and Bayesian inference--full tensor-product extensions of Gauss--Legendre suffer the curse of dimensionality, requiring $n^d$ evaluations in $d$ dimensions. The classical remedy is the Smolyak construction @Smolyak1963, which combines truncated tensor products of one-dimensional rules at carefully chosen levels to preserve the asymptotic convergence rate of the full grid while drastically reducing the node count. The contemporary state of the art is summarised in the Garcke--Griebel monograph @GarckeGriebel2013. Constantine, Eldred and Phipps @ConstantineEldredPhipps2012 and Conrad and Marzouk @ConradMarzouk2013 then exposed a subtle hazard known as _internal aliasing_: when sparse grids are used to compute high-degree polynomial coefficients, the missing nodes contaminate the projection unless the one-dimensional exactness degree is carefully tied to the highest represented polynomial degree. They derived sharp criteria of the form $q(m) = floor(a(m)\/2)$ that eliminate the aliasing error, and they pointed out that nested rules such as Clenshaw--Curtis are far better suited to sparse grids than Gauss--Legendre, because their nodes form perfect subsets across grid levels and previous function evaluations can be reused. This is one of the rare contexts in which the supposedly inferior exactness degree of Clenshaw--Curtis becomes a decisive practical advantage.

=== Breaking Exactness in Modern Spectral PDE Solvers

The frontier of the exactness question is now firmly inside spectral discretisations of nonlinear PDEs. Spectral Galerkin and pseudospectral methods have traditionally insisted on quadrature rules that integrate the relevant inner products _exactly_, in the belief that this exactness was necessary for energy stability, conservation, and the preservation of discrete maximum principles. A series of recent papers has shown that this assumption is false--and unhelpfully restrictive. Wu and Yuan @WuYuan2023 presented a spectral method for the Allen--Cahn equation on the sphere $S^2$ that abandons quadrature exactness entirely, replacing it with a _restricted-isometry_ relation drawn from Marcinkiewicz--Zygmund quadrature theory. The method retains the discrete maximum principle and unconditional energy stability, while permitting time steps that no exactness-based scheme can tolerate.

In a parallel development, Glaubitz, Nordström and Öffner @GlaubitzNordstromOffner2023 generalised the Summation-by-Parts (SBP) framework, the cornerstone of provably stable high-order finite difference and spectral methods, from polynomial spaces to arbitrary function spaces. They showed that the SBP property--and with it, energy stability and discrete conservation--is an artifact of the structure of the boundary projection and integration matrices, not of any underlying polynomial exactness. Cai, Cheng, Wu and Qiu @CaiEtAl2026 introduced _transport polynomial exactness_ (TPE), a mesh-motion-independent generalisation of free-stream preservation for moving-mesh methods, in which the relevant condition is the exact advection of degree-$k$ polynomials by the discrete scheme rather than the static integration of fixed polynomials at each time step.

The unifying lesson of these developments is the same one that closes our chapter: stability, conservation, and optimal convergence are intrinsic properties of the chosen approximation space and its discrete inner products, not fragile artefacts of forcing algebraic zeros via exactness. The spectral practitioner of 2026 should pursue the alignment of the approximation space with the analytic and geometric properties of the integrand, and let exactness emerge--or not--as a by-product.

// ============================================================================
== Summary <sec-quad-summary>
// ============================================================================

This chapter has explored the relationship between polynomial exactness and practical quadrature accuracy, guided by two papers of Trefethen @TrefethenCC2008 @TrefethenExactness2022. The key findings are:

- *Newton--Cotes diverges* for most analytic functions, despite degree $n$ exactness, because its equispaced nodes lead to exponentially growing weights and Runge-type oscillations.
- *Monomials are flattering tests*: $x^k$ has numerical degree only $cal(O)(k^(1\/2))$, so it fails to expose the ill-conditioning of Newton--Cotes. Chebyshev polynomials $T_k$ are the honest diagnostic.
- *Clenshaw--Curtis nearly matches Gauss* for all functions encountered in practice. The aliasing identity $T_(n+p)(x_j) = T_(n-p)(x_j)$ on Chebyshev grids ensures that quadrature errors for the "shaded" Chebyshev coefficients are small, effectively doubling the convergence rate.
- *In the complex plane*, Gauss and Clenshaw--Curtis correspond to rational approximants of $log((z+1)\/(z-1))$ that are virtually identical near $[-1, 1]$, explaining their similar performance for functions with nearby singularities.
- *Gauss--Hermite is paradoxically inefficient* on unbounded domains: truncating to $[-cal(O)(n^(1\/3)), cal(O)(n^(1\/3))]$ and using standard quadrature gives $cal(O)(e^(-C n^(2\/3)))$ convergence, faster than Gauss--Hermite's $cal(O)(e^(-C n^(1\/2)))$.
- *The periodic trapezoidal rule* is exponentially accurate for smooth periodic functions, a fact invisible to piecewise-linear exactness analysis but natural from the trigonometric perspective.

The overarching lesson for spectral methods practitioners: *choose the approximation space first, then the quadrature*. Exactness can mislead; aliasing explains why Clenshaw--Curtis works; approximation spaces explain what really matters.

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
        table.cell(fill: rgb("#142D6E").lighten(85%))[*Highlight*],
      ),
      table.hline(stroke: 0.5pt + luma(180)),
      [15.1], [Node distributions and polynomial exactness], [Vandermonde, FFT, Golub--Welsch], [Visualises the $n$ vs.\u{a0}$2n + 1$ exactness gap],
      [15.2], [Newton--Cotes on the Runge function], [Ill-conditioned weights], [$cal(O)(2^n)$ divergence vs.\u{a0}rapid convergence],
      [15.3], [Monomial vs.\u{a0}Chebyshev exactness table], [$|E_n (x^k)|$ vs.\u{a0}$|E_n (T_k)|$], [Monomials hide $cal(O)(10^4)$ Chebyshev errors],
      [15.4], [Building Gauss--Legendre and Clenshaw--Curtis], [Golub--Welsch + DCT-I via FFT], [Two seven-line algorithms validated],
      [15.5], [Six-function convergence race], [Gauss vs.\u{a0}Clenshaw--Curtis benchmark], [Trefethen 2008, Figure 2 reproduced],
      [15.6], [Aliasing and Chebyshev coefficients], [$T_(n+p) = T_(n-p)$ on the grid], [Why CC matches Gauss for $|x|^3$],
      [15.7], [Complex-plane error portraits], [Padé approximant of $log((z+1)\/(z-1))$], [Near-field contours dominate the error],
      [15.8], [Gauss--Hermite wasted weights], [Counting weights $< 10^(-16)$], [836 of 1000 weights vanish at $n = 1000$],
      [15.9], [Gauss--Hermite vs.\u{a0}truncated trapezoidal], [Domain truncation $[-6, 6]$], [Trapezoidal beats Hermite by orders of magnitude],
      [15.10], [Experimental convergence rates], [Plot vs.\u{a0}$n^(1\/2)$ and $n^(2\/3)$], [Confirms Theorem 5.1 of Trefethen 2022],
    ),
  ),
  caption: [Summary of computational études in this chapter.],
) <tbl-ch15-summary>

// ============================================================================
== Exercises <sec-quadrature-exercises>
// ============================================================================

*Exercise 15.1* (_Newton--Cotes Weight Growth_). For $n = 10, 20, 30, 40, 50$, compute the Newton--Cotes weights $w_0, dots, w_n$ and plot $sum_(j=0)^n |w_j|$ versus $n$. Verify that this sum grows approximately as $2^n \/ (n log n)$. Then compute the quadrature error for $f(x) = e^x$ and show that, despite the exact integral being $e - e^(-1) approx 2.3504$, the Newton--Cotes estimate deteriorates catastrophically for $n gt.eq.slant 30$.

*Exercise 15.2* (_Simpson's Rule Bonus_). Simpson's rule is Newton--Cotes with $n = 2$, which has exactness degree 3 (one higher than the expected degree 2). (a) Prove that any interpolatory quadrature formula on $n + 1$ symmetrically placed nodes has exactness degree $n + 1$ when $n$ is even. (b) Show that this "bonus degree" does not save Newton--Cotes from divergence for large $n$.

*Exercise 15.3* (_Clenshaw--Curtis via the DCT_). Implement Clenshaw--Curtis quadrature using the Type-I discrete cosine transform (DCT-I) rather than the FFT. Compare the accuracy and timing of both implementations for $n = 2^k$, $k = 2, 3, dots, 14$. Which is faster?

*Exercise 15.4* (_The Kink Phenomenon_). For $f(x) = 1\/(1 + 16x^2)$ (poles at $x = plus.minus i\/4$), compute Clenshaw--Curtis quadrature errors for $n = 1, 2, dots, 120$. Show that around $n approx 50$, the convergence rate changes abruptly (a "kink"). Explain this in terms of the interpolation curve of the Clenshaw--Curtis rational approximant crossing the poles. (Hint: see Figure 7 of @TrefethenCC2008.)

*Exercise 15.5* (_Band-Limited Quadrature_). Implement the band-limited quadrature formula described in @TrefethenExactness2022, using the parameter $c = pi n - 12 log n$. For $n = 50$, compute the errors in integrating $e^(i k x)$ for $k = 0, 1, dots, 200$ and compare with Gauss quadrature. Verify that the band-limited formula is more efficient by approximately a factor $pi\/2$ in convergence rate.

*Exercise 15.6* (_Gauss--Laguerre Paradox_). Repeat the analysis of @sec-beyond-gauss for Gauss--Laguerre quadrature, which targets integrals $integral_0^infinity e^(-x) f(x) dif x$. For $f(x) = cos(x^2)$, compare the convergence of Gauss--Laguerre with that of Gauss--Legendre on a truncated interval $[0, L]$ with $L = 36$. Show that Gauss--Laguerre suffers from the same wasted-weight phenomenon as Gauss--Hermite.

*Exercise 15.7* (_Spectral Methods Connection_). Consider the Poisson equation $u''(x) = f(x)$ on $[-1, 1]$ with $u(plus.minus 1) = 0$, using a Chebyshev spectral method with $N$ interior points. The right-hand side integral $integral_(-1)^1 f(x) T_k (x) \/ sqrt(1 - x^2) dif x$ arising in the Galerkin formulation can be evaluated by Gauss--Chebyshev quadrature (exact for degree $2N - 1$) or by Clenshaw--Curtis-type quadrature (exact for degree $N$). Implement both approaches for $f(x) = e^x$ and show that the difference in the solution error is negligible. This demonstrates that the "half exactness" of Clenshaw--Curtis does not degrade spectral method accuracy.

*Exercise 15.8* (_Reflection Essay_). Write a two-page reflection titled: "When should a spectral methods practitioner prefer Chebyshev, Legendre, Hermite, or Fourier quadrature?" Your argument should be framed in terms of approximation spaces, aliasing, analyticity, and node placement, not just exactness degree.
