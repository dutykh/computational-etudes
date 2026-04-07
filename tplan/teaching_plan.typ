// tplan/teaching_plan.typ
// A Tentative Teaching Plan for MATH 794 - Spring 2026

#set document(
  title: "MATH 794 - Tentative Teaching Plan",
  author: "Dr. Denys Dutykh",
)

#set page(
  paper: "a4",
  margin: (x: 2.5cm, y: 2.5cm),
  header: context {
    if counter(page).get().first() > 1 [
      #set text(size: 9pt, fill: luma(100))
      #h(1fr) MATH 794 -- Spring 2026 #h(1fr)
      #line(length: 100%, stroke: 0.5pt + luma(200))
    ]
  },
  footer: context {
    set text(size: 9pt, fill: luma(100))
    line(length: 100%, stroke: 0.5pt + luma(200))
    v(0.3em)
    [
      #h(1fr)
      #counter(page).display("1 / 1", both: true)
      #h(1fr)
    ]
  },
)

#set text(
  font: "New Computer Modern",
  size: 11pt,
)

#set par(
  justify: true,
  leading: 0.65em,
)

// Colors
#let accent-color = rgb("#1a5f7a")
#let light-accent = rgb("#e8f4f8")
#let completed-color = rgb("#d4edda")
#let current-color = rgb("#fff3cd")
#let exam-color = rgb("#f8d7da")

// Title block
#align(center)[
  #block(
    width: 100%,
    inset: 1.5em,
    radius: 8pt,
    fill: accent-color,
  )[
    #set text(fill: white)
    #text(size: 24pt, weight: "bold")[Computational Études]
    #v(0.3em)
    #text(size: 14pt)[A Spectral Approach]
    #v(1em)
    #line(length: 60%, stroke: 1pt + white.transparentize(50%))
    #v(1em)
    #text(size: 18pt, weight: "semibold")[A Tentative Teaching Plan]
  ]
]

#v(1em)

// Course info box
#grid(
  columns: (1fr, 1fr),
  gutter: 1em,
  block(
    width: 100%,
    inset: 1em,
    radius: 6pt,
    stroke: 1pt + luma(200),
  )[
    #set text(size: 10pt)
    #text(weight: "bold", fill: accent-color)[Course Information]
    #v(0.5em)
    *Course:* MATH 794 -- 01 (CRN 21296) \
    *Semester:* Spring 2026 \
    *Schedule:* Mon & Wed, 16:30 -- 17:45 \
    *Location:* G-G02013
  ],
  block(
    width: 100%,
    inset: 1em,
    radius: 6pt,
    stroke: 1pt + luma(200),
  )[
    #set text(size: 10pt)
    #text(weight: "bold", fill: accent-color)[Instructor]
    #v(0.5em)
    *Name:* Dr. Denys Dutykh \
    *Department:* Mathematics \
    *Institution:* Khalifa University \
    *Location:* Abu Dhabi, UAE
  ],
)

#v(1em)

// Legend
#block(
  width: 100%,
  inset: 0.8em,
  radius: 4pt,
  fill: luma(245),
)[
  #set text(size: 9pt)
  #text(weight: "bold")[Legend:]
  #h(1em)
  #box(width: 1em, height: 1em, fill: completed-color, radius: 2pt)
  #h(0.3em) Completed
  #h(1.5em)
  #box(width: 1em, height: 1em, fill: current-color, radius: 2pt)
  #h(0.3em) Current
  #h(1.5em)
  #box(width: 1em, height: 1em, fill: exam-color, radius: 2pt)
  #h(0.3em) Examination
  #h(1.5em)
  #box(width: 1em, height: 1em, fill: white, stroke: 0.5pt + luma(200), radius: 2pt)
  #h(0.3em) Planned
]

#v(1em)

// Helper function for status styling
#let status-fill(status) = {
  if status == "completed" { completed-color }
  else if status == "current" { current-color }
  else { white }
}

// Teaching plan table
#set table(
  stroke: (x, y) => (
    left: 0.5pt + luma(200),
    right: 0.5pt + luma(200),
    top: if y == 0 { 1.5pt + accent-color } else { 0.5pt + luma(200) },
    bottom: 0.5pt + luma(200),
  ),
)

#show figure: set block(breakable: true)

#figure(
  table(
    columns: (auto, auto, 1fr, 1fr),
    align: (center, center, left, left),
    inset: 0.7em,

    // Header (repeats on page breaks)
    table.header(
      table.cell(fill: accent-color)[#text(fill: white, weight: "bold")[Week]],
      table.cell(fill: accent-color)[#text(fill: white, weight: "bold")[Dates]],
      table.cell(fill: accent-color)[#text(fill: white, weight: "bold")[Monday Lecture]],
      table.cell(fill: accent-color)[#text(fill: white, weight: "bold")[Wednesday Lecture]],
    ),

    // Week 1 (Jan 12-15)
    table.cell(fill: completed-color)[*1*],
    table.cell(fill: completed-color)[Jan 12--14],
    table.cell(fill: completed-color)[
      Course introduction \
      #text(size: 9pt, fill: luma(100))[Syllabus & overview]
    ],
    table.cell(fill: completed-color)[
      *Ch. 2:* Classical PDEs \
      #text(size: 9pt, fill: luma(100))[Heat equation -- separation of variables]
    ],

    // Week 2 (Jan 19-21)
    table.cell(fill: completed-color)[*2*], table.cell(fill: completed-color)[Jan 19--21],
    table.cell(fill: completed-color)[
      *Ch. 3:* Mise en bouche \
      #text(size: 9pt, fill: luma(100))[Method of weighted residuals, collocation example, collocation vs Galerkin -- Project I assigned]
    ],
    table.cell(fill: completed-color)[
      *Ch. 4:* The Geometry of Nodes \
      #text(size: 9pt, fill: luma(100))[Lagrange interpolation, Runge phenomenon, potential theory]
    ],

    // Week 3 (Jan 26-29)
    table.cell(fill:completed-color)[*3*], table.cell(fill:completed-color)[Jan 26--28],
    table.cell(fill:completed-color)[
      *Ch. 5:* Differentiation Matrices \
      #text(size: 9pt, fill: luma(100))[FD stencils, spectral matrices, Fornberg algorithm]
    ],
    table.cell(fill:completed-color)[
      *Ch. 6:* Smoothness and Spectral Accuracy \
      #text(size: 9pt, fill: luma(100))[Smoothness, spectral convergence, exponential accuracy]
    ],

    // Week 4 (Feb 2-5)
    table.cell(fill:completed-color)[*4*], table.cell(fill:completed-color)[Feb 2--4],
    table.cell(fill:completed-color)[
      *Ch. 7:* Chebyshev Differentiation Matrices \
      #text(size: 9pt, fill: luma(100))[Chebyshev nodes, $D_N$ matrix, negative sum trick]
    ],
    table.cell(fill:completed-color)[
      *Ch. 8:* Boundary Value Problems \
      #text(size: 9pt, fill: luma(100))[1D/2D BVPs, matrix stripping, Newton iteration -- Project II assigned]
    ],

    // Week 5 (Feb 9-12)
    table.cell(fill:completed-color)[*5*], table.cell(fill:completed-color)[Feb 9--11],
    table.cell(fill:completed-color)[
      *Project I:* Students Presentations \
      #text(size: 9pt, fill: luma(100))[Oral presentations and discussions of Project I results]
    ],
    table.cell(fill:completed-color)[
      *Ch. 9:* Fourier Space on Grids \
      #text(size: 9pt, fill: luma(100))[Fourier transform, DFT/FFT, aliasing, sinc interpolation]
    ],

    // Week 6 (Feb 16-19)
    table.cell(fill:completed-color)[*6*], table.cell(fill:completed-color)[Feb 16--18],
    table.cell(fill:completed-color)[
      *Ch. 10:* Spectral PDE Solvers \
      #text(size: 9pt, fill: luma(100))[Chebyshev--Fourier recap, Chebyshev differentiation via FFT, method of lines]
    ],
    table.cell(fill:completed-color)[
      *Ch. 11:* Pseudo-spectral PDE Solvers \
      #text(size: 9pt, fill: luma(100))[Fourier Pseudospectral Methods -- Project III assigned]
    ],

    // Week 7 (Feb 23-26)
    table.cell(fill:completed-color)[*7*], table.cell(fill:completed-color)[Feb 23--25],
    table.cell(fill:completed-color)[
      *Ch. 12:* Spectral Methods in Polar Coordinates \
      #text(size: 9pt, fill: luma(100))[Coordinate singularity, Fornberg's doubling trick, Chebyshev--Fourier disk Laplacian, Bessel eigenmodes]
    ],
    table.cell(fill:completed-color)[
      *Ch. 13:* Advanced Boundary Conditions \
      #text(size: 9pt, fill: luma(100))[Lifting (Method I) vs tau (Method II), Robin/Neumann BCs, nonlinear BCs, QEP and quasinormal modes]
    ],

    // Week 8 (Mar 2-5)
    table.cell(fill: exam-color)[*8*],
    table.cell(fill: exam-color)[Mar 2--4],
    table.cell(fill: exam-color)[
      *Project II:* Students Presentations \
      #text(size: 9pt, fill: luma(100))[Oral presentations and discussions of Project II results]
    ],
    table.cell(fill: exam-color)[
      *Case Study* assigned \
      #text(size: 9pt, fill: luma(100))[Replaces Midterm Examination]
    ],

    // Spring Break: Mar 9-22 (includes Eid Al Fitr)
    table.cell(fill: luma(240))[--],
    table.cell(fill: luma(240))[Mar 9--22],
    table.cell(fill: luma(240), colspan: 2)[
      #align(center)[
        #text(style: "italic", fill: luma(100))[Spring Break & Eid Al Fitr -- No classes]
      ]
    ],

    // Week 9 (Mar 23-25)
    table.cell(fill: completed-color)[*9*], table.cell(fill: completed-color)[Mar 23--25],
    table.cell(fill: completed-color)[
      *Ch. 14:* Higher-Order BVPs \
      #text(size: 9pt, fill: luma(100))[Polynomial trick for clamped BCs, $u^((4)) = e^x$ (Étude 14.1), beam eigenmodes (Étude 14.2)]
    ],
    table.cell(fill: completed-color)[
      *Project III:* Students Presentations \
      #text(size: 9pt, fill: luma(100))[Oral presentations and discussions of Project III results]
    ],

    // Week 10 (Mar 30 - Apr 2)
    table.cell(fill: completed-color)[*10*], table.cell(fill: completed-color)[Mar 30--Apr 1],
    table.cell(fill: completed-color)[
      Spectral Method Convergence for QNMs (I) \
      #text(size: 9pt, fill: luma(100))[Overview of the convergence proof for spectral methods applied to quasinormal modes -- Slides: `slides/QNM-SpectralConvergence.pdf`]
    ],
    table.cell(fill: completed-color)[
      Spectral Method Convergence for QNMs (II) \
      #text(size: 9pt, fill: luma(100))[Continuation and discussion of the convergence proof -- Project IV assigned]
    ],

    // Week 11 (Apr 6-9)
    table.cell(fill: current-color)[*11*], table.cell(fill: current-color)[Apr 6--8],
    table.cell(fill: current-color)[
      *Ch. 15:* Quadrature in Spectral Methods \
      #text(size: 9pt, fill: luma(100))[Exactness principle, Newton--Cotes failure, Clenshaw--Curtis vs Gauss, aliasing in Chebyshev space]
    ],
    table.cell(fill: current-color)[
      TBA \
      #text(size: 9pt, fill: luma(100))[To be announced]
    ],

    // Week 12 (Apr 13-16)
    [*12*], [Apr 13--15],
    table.cell(fill: exam-color)[
      *Case Study:* Oral Presentations \
      #text(size: 9pt, fill: luma(100))[Presentations and report submissions -- Mon, Apr 13]
    ],
    [
      TBA \
      #text(size: 9pt, fill: luma(100))[To be announced]
    ],

    // Week 13 (Apr 20-23)
    [*13*], [Apr 20--22],
    [
      TBA \
      #text(size: 9pt, fill: luma(100))[To be announced]
    ],
    [
      TBA \
      #text(size: 9pt, fill: luma(100))[To be announced]
    ],

    // Week 14 (Apr 27-30)
    [*14*], [Apr 27--29],
    table.cell(fill: exam-color)[
      *Project IV:* Oral Presentations \
      #text(size: 9pt, fill: luma(100))[Presentations and report submissions -- Mon, Apr 27]
    ],
    [
      TBA \
      #text(size: 9pt, fill: luma(100))[To be announced]
    ],

    // Finals Week (May 4-14)
    table.cell(fill: exam-color)[--],
    table.cell(fill: exam-color)[May 4--14],
    table.cell(fill: exam-color, colspan: 2)[
      #align(center)[
        #text(weight: "bold")[Final Examinations Period]
      ]
    ],
  ),
  caption: [Teaching schedule for MATH 794 -- Spring 2026],
) <tab-schedule>

#v(1em)

// Notes section
#block(
  width: 100%,
  inset: 1em,
  radius: 6pt,
  fill: light-accent,
)[
  #text(weight: "bold", fill: accent-color)[Notes]
  #v(0.5em)
  #set text(size: 10pt)
  - This schedule is tentative and may be adjusted based on class progress.
  - Chapter numbers refer to _Computational Études: A Spectral Approach_.
  - Office hours are available by appointment.
  - All course materials are available in the course repository.
]

#pagebreak()

// Course Projects section
#text(size: 16pt, weight: "bold", fill: accent-color)[Course Projects]
#v(0.5em)

#block(
  width: 100%,
  inset: 1em,
  radius: 6pt,
  stroke: 1pt + luma(200),
)[
  #text(weight: "bold", fill: accent-color)[Project I: Spectral Methods for Boundary Value Problems]
  #h(1fr)
  #text(size: 9pt, fill: luma(100))[(Based on Chapter 3 -- Assigned: Jan 19, 2026)]
  #v(0.5em)
  #set text(size: 10pt)

  Implement spectral methods (collocation and/or Galerkin) to solve a boundary value problem of your choice. Requirements:

  + Find a linear BVP of at least second order with *non-homogeneous* boundary conditions.
  + Construct an exact solution to this problem (or design the problem around a known exact solution).
  + Choose appropriate basis functions that satisfy the boundary conditions identically.
  + Apply the *Collocation* method, *Galerkin* method, or both to obtain a numerical solution.
    - If using *only* the Collocation method, the BVP must have *non-constant coefficients*.
    - If applying *both* methods, the problem may have constant coefficients (though non-constant coefficients are also welcome).
  + Present results in a table that includes the error at the collocation points.
  + Provide graphical representations showing: the exact solution, the numerical approximation, and the difference (error) between them.

  *Deliverables:* A written report, a short oral presentation (15 min max), and a live demonstration of your computational codes. Grading is based on relative performance within the class.
]

#v(1em)

#block(
  width: 100%,
  inset: 1em,
  radius: 6pt,
  stroke: 1pt + luma(200),
)[
  #text(weight: "bold", fill: accent-color)[Project II: Spectral Collocation for BVPs]
  #h(1fr)
  #text(size: 9pt, fill: luma(100))[(Based on Chapters 6--8 -- Assigned: Feb 2, 2026)]
  #v(0.5em)
  #set text(size: 10pt)

  Implement Chebyshev spectral collocation to solve boundary value problems. Requirements:

  + Implement the Chebyshev differentiation matrix $D_N$ using the negative sum trick.
  + Choose a second-order BVP (linear or nonlinear) with known exact solution.
  + Solve the BVP using spectral collocation with matrix stripping for boundary conditions.
  + Present a convergence study: error vs. $N$ on a semilog plot.
  + For bonus: solve a 2D problem using tensor products.

  *Deliverables:* A written report, a short oral presentation (15 min max), and a live demonstration of your computational codes. Grading is based on relative performance within the class.
]

#v(1em)

#block(
  width: 100%,
  inset: 1em,
  radius: 6pt,
  stroke: 1pt + luma(200),
)[
  #text(weight: "bold", fill: accent-color)[Project III: Fourier Pseudospectral Simulation of Periodic PDEs]
  #h(1fr)
  #text(size: 9pt, fill: luma(100))[(Based on Chapters 9--11 -- Assigned: Feb 18, 2026)]
  #v(0.5em)
  #set text(size: 10pt)

  Implement a Fourier pseudospectral method to simulate the dynamics of a PDE of your choice on a periodic domain. Requirements:

  + Choose a time-dependent PDE in $(x, t)$ variables posed on a periodic spatial domain.
    - The PDE may be linear (e.g. advection, diffusion, Schrödinger) or nonlinear (e.g. Korteweg--de Vries, Burgers, Allen--Cahn, Kuramoto--Sivashinsky).
  + Discretize spatial derivatives using the Fourier pseudospectral method (via the FFT).
  + Advance in time using an appropriate ODE solver (e.g. Runge--Kutta, exponential integrator, or split-step method).
  + Validate your implementation against an exact solution, a conserved quantity, or a known qualitative behaviour of the chosen PDE.
  + Present space-time visualizations of the computed dynamics (e.g. waterfall plots, animations, or snapshots at selected times).
  + For bonus: extend the simulation to two or three spatial dimensions.

  *Deliverables:* A written report, a short oral presentation (15 min max), and a live demonstration of your computational codes. Grading is based on relative performance within the class.
]

#v(1em)

#block(
  width: 100%,
  inset: 1em,
  radius: 6pt,
  stroke: 1pt + luma(200),
)[
  #text(weight: "bold", fill: accent-color)[Project IV: Open Topic in Spectral Methods]
  #h(1fr)
  #text(size: 9pt, fill: luma(100))[(Based on Chapters 12--14 -- Assigned: Apr 1, 2026)]
  #v(0.5em)
  #set text(size: 10pt)

  Choose any topic from the textbook _Computational Étude: A Spectral Approach_, starting from Chapter 12 and until the end of the book. The choice of topic is left to each student according to their personal taste and scientific preferences. Requirements:

  + Select a problem, method, or application from Chapters 12--14 of the textbook that interests you.
    - Examples include (but are not limited to): spectral methods in polar coordinates, eigenmodes on the disk, advanced boundary conditions (Robin, radiation, quasinormal modes), higher-order BVPs, the biharmonic operator, Orr--Sommerfeld stability, pseudospectra, or the Kuramoto--Sivashinsky equation.
  + Implement the chosen method or solve the chosen problem computationally, in any of the three programming languages adopted in our course (Python, MATLAB, or Julia).
  + Validate your implementation against known analytical results, reference solutions, or established benchmarks from the literature.
  + Present your findings through appropriate visualizations (convergence plots, mode shapes, space-time diagrams, _etc._).
  + Discuss the strengths and limitations of the spectral approach for your chosen problem.

  *Deliverables:* A written report, a short oral presentation (15 min max), and a live demonstration of your computational codes. Grading is based on relative performance within the class.
]

#v(2em)

#align(right)[
  #set text(size: 10pt, fill: luma(100))
  _Last updated: #datetime.today().display("[month repr:long] [day], [year]")_
]
