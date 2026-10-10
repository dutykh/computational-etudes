// textbook/chapters/colophon.typ
// Final printed backmatter: descriptive blurb, keywords + MSC, author bio.
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Last modified: April 2026

#import "../styles/template.typ": dropcap

#heading(level: 1, numbering: none)[About this Book]

#dropcap[*Computational Études: A Spectral Approach* is a graduate-level textbook on spectral methods for differential equations, organised around 132 computational *études* --- short, self-contained studies that pair a focused mathematical concept with a runnable implementation in Python, MATLAB, and Julia.] The book traces a single foundation spine from interpolation and differentiation matrices through Chebyshev and Fourier discretisations to a verification-as-ethic treatment of linear spectral eigenproblems, then branches out to coordinate transformations, unbounded domains, polar geometry, advanced boundary conditions, and the special tricks of the spectral researcher. Every figure, table, and printed numerical statement is reproducible from code that ships alongside the manuscript. The pedagogy follows the *étude* model literally: each chapter assigns a small, complete experiment whose outcome teaches a single principle, and each étude closes with a short *Takeaway* paragraph that names the lesson the experiment proves.

#v(1.5em)

#heading(level: 2, numbering: none)[Keywords]

- spectral methods
- pseudospectral collocation
- Chebyshev polynomials
- Fourier methods
- scientific computing
- numerical analysis

#heading(level: 2, numbering: none)[Subject classification (MSC 2020)]

#let _navy = rgb(20, 45, 110)
#block(
  stroke: (top: 1pt + _navy, bottom: 1pt + _navy),
  inset: 0pt,
  width: 100%,
  table(
    columns: (auto, 1fr),
    align: (right + horizon, left + horizon),
    inset: (x: 1em, y: 0.45em),
    stroke: none,
    table.hline(stroke: 0.5pt + _navy),
    [*65N35*], [Spectral, collocation, and related methods for boundary value problems involving PDEs (*primary*)],
    [*65T50*], [Discrete and fast Fourier transforms],
    [*65L10*], [Numerical solution of boundary value problems for ODEs],
    [*65L60*], [Galerkin and finite element methods for IVPs],
    [*65F15*], [Numerical computation of eigenvalues and eigenvectors of matrices],
    [*30E20*], [Integration in the complex plane (background for Appendix B)],
    table.hline(stroke: 0.5pt + _navy),
  ),
)

#v(2em)

#heading(level: 1, numbering: none)[About the Author]

#dropcap[*Denys Dutykh* is Associate Professor of Mathematics at Khalifa University of Science and Technology in Abu Dhabi, UAE.] He was born in Polohy, Ukraine, and graduated as valedictorian of the master's programme in numerical methods for continuum mechanics at the École Normale Supérieure de Cachan (now ENS Paris-Saclay) in 2005. He earned his PhD in two years at the Centre de Mathématiques et de Leurs Applications (now Centre Borelli) under Frédéric Dias, defending a thesis on the mathematical modelling of tsunami waves in 2007, and a habilitation in environmental mathematical modelling at the University of Savoie in 2010.

In 2008 he joined the *Centre National de la Recherche Scientifique* (CNRS) as a permanent researcher, attached to the Institut National des Sciences Mathématiques et de leurs Interactions and the Laboratory of Mathematics of the University of Savoie Mont Blanc. In 2012--2013 he held a visiting position at University College Dublin's School of Mathematics and Statistics, contributing to the ERC Advanced Grant project *MULTIWAVE*. In 2022 he moved to Khalifa University, where he continues research and teaching in applied mathematics with a focus on spectral methods, water waves, and the mathematical modelling of physical phenomena.

His research interests span the analytical and numerical analysis of nonlinear dispersive waves, hydroelasticity, tsunami modelling, and the spectral methods that are the subject of this volume. He is the author of more than one hundred research articles in applied mathematics, mathematical physics, and computational science.

#v(1.5em)
#align(right)[#text(size: 0.9em, fill: _navy.lighten(20%), style: "italic")[#datetime.today().display("[month repr:long] [year]") --- Abu Dhabi, UAE]]
