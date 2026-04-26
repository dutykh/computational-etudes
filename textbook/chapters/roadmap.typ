// textbook/chapters/roadmap.typ
// How to Read This Book / Roadmap
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Last modified: April 2026

#import "../styles/template.typ": dropcap

#heading(level: 1, numbering: none)[How to Read This Book]

#dropcap[The book is written for linear reading: each chapter draws on the techniques and intuitions developed in the chapters that precede it. But a text of this size will inevitably be read non-linearly --- by lecturers selecting a one-semester course, by self-learners interested only in a particular family of methods, and by researchers consulting individual études as standalone references.] This roadmap makes the chapter-to-chapter dependencies explicit, proposes three suggested course tracks, and identifies the chapters that may be safely skipped.

#let navy = rgb(20, 45, 110)
#let sky = rgb(120, 150, 210)

#heading(level: 2, numbering: none)[Chapter dependencies]

Most chapters lie on a single _foundation spine_, Ch 1 $arrow$ Ch 2 $arrow$ $dots$ $arrow$ Ch 8; after Chapter 8, the topic chapters branch out. The table records each chapter's *immediate* prerequisites only: implicitly, every chapter assumes the chapters its prerequisites assume. Chapter numbers in the right column are the *minimal* set required to follow the étude-level prose; the full literature reviews at the end of each chapter sometimes assume more.

#block(
  stroke: (top: 1.5pt + navy, bottom: 1.5pt + navy),
  inset: 0pt,
  table(
    columns: (auto, 1fr),
    align: (left, left),
    inset: (x: 1em, y: 0.45em),
    stroke: none,
    table.hline(stroke: 0.75pt + navy),
    table.header(
      table.cell(fill: navy.lighten(85%))[*Chapter*],
      table.cell(fill: navy.lighten(85%))[*Immediate prerequisites*],
    ),
    table.hline(stroke: 0.5pt + luma(180)),
    [1.~Introduction],                                [---],
    [2.~Classical second-order PDEs],                 [1],
    [3.~Mise en bouche (collocation/Galerkin)],       [1, 2],
    [4.~Geometry of nodes],                           [1, 3],
    [5.~Differentiation matrices],                    [4],
    [6.~Smoothness and spectral accuracy],            [4, 5],
    [7.~Chebyshev differentiation],                   [4, 5, 6],
    [8.~Boundary value problems],                     [7],
    [9.~Fourier grids on the line and the circle],    [5, 6],
    [10.~Spectral PDE solvers (Chebyshev + Fourier)], [7, 8, 9],
    [11.~Fourier pseudospectral for nonlinear PDEs],  [9, 10],
    [12.~Polar coordinates],                          [8, 10],
    [13.~Advanced boundary conditions],               [8],
    [14.~Higher-order BVPs],                          [8, 13],
    [15.~Quadrature: when exactness misleads],        [4, 6, 7],
    [16.~Periodic quadrature (trapezoidal rule)],     [6, 9],
    [17.~Time stepping, stability, CFL],              [10, 11],
    [18.~Linear spectral eigenproblems],              [8, 14],
    [19.~Coordinate transformations],                 [7],
    [20.~Spectral methods on unbounded intervals],    [7, 19],
    [21.~Special tricks (bonus)],                     [any of 7--20],
  ),
)

#heading(level: 2, numbering: none)[Three course tracks]

The textbook supports three natural lecture-course formats. Each track is given as a list of chapters in the order in which they should be read; the optional chapters in parentheses can be added if time allows or skipped if not.

#let track-block(title, weeks, chapters, blurb) = block(
  stroke: (left: 2pt + navy),
  inset: (left: 12pt, top: 6pt, bottom: 6pt, right: 0pt),
  spacing: 0.7em,
  width: 100%,
  [
    #text(weight: "semibold", fill: navy)[#title (#weeks).] #h(0.4em) #blurb \
    #v(0.2em) #text(weight: "medium")[Read:] #h(0.4em) #chapters
  ]
)

#track-block(
  "Foundations track",
  "6 weeks, intensive",
  [Ch 1 -- 2 -- 3 -- 4 -- 5 -- 7 -- 8],
  [The minimal Chebyshev-pseudospectral course. Drops Ch 6 if students bring approximation-theory background; otherwise keep it. Output: a working understanding of one-dimensional spectral BVPs.],
)

#track-block(
  "Standard intro track",
  "9 weeks, semester half",
  [Ch 1 -- 2 -- 3 -- 4 -- 5 -- 6 -- 7 -- 8 -- 9 -- 10 -- (11)],
  [Adds the Fourier story (Ch 9 -- 11) and the Chebyshev-FFT PDE solvers of Ch 10. Étude 11.6 (Kuramoto--Sivashinsky) makes a good capstone if Ch 11 is included.],
)

#track-block(
  "Full semester track",
  "13 weeks",
  [Ch 1 -- 10, then a selection from {11, 13, 14, 15 + 16, 17, 18, 20}],
  [Foundations + Fourier + a topic block of the lecturer's choice. Two coherent topic packs: _stability and time integration_ (Ch 11, 17, 21) or _hard problems_ (Ch 13, 14, 18, 20). Ch 12 (polar coordinates) and Ch 19 (coordinate transformations) are optional throughout.],
)

#heading(level: 2, numbering: none)[What you can skip]

The following chapters lie outside the foundation spine and may be visited or omitted independently:

- *Ch 12 (Polar coordinates).* Self-contained once Ch 8 and 10 are read; useful if the audience cares about disk- or annulus-shaped domains. Skip without loss of continuity if the audience works exclusively on intervals or rectangles.
- *Ch 15 + 16 (Quadrature and periodic quadrature).* These two read as a coherent pair, but they are needed elsewhere in the book only as references for spectral-trapezoidal-rule analysis. Skip both if integration is not the audience's primary concern.
- *Ch 19 (Coordinate transformations)* is a prerequisite of Ch 20, but Ch 20 itself is optional. Both belong to a "physics of unbounded domains" topic block that is independent of the rest of the volume.
- *Ch 21 (Special tricks)* is *supplementary throughout*: every étude in Ch 21 stands on a single technique from Ch 7 -- 20 and can be read immediately after the corresponding source chapter. Treat it as a reference rather than a sequential read.
- The level-2 *Literature overview* sections that close most chapters are deliberately optional --- they survey the modern (2020--2026) research landscape and are aimed at the reader who has just finished the chapter and wants to know where the topic is going next.

#heading(level: 2, numbering: none)[A note on the études]

The 132 Computational Études are the *active learning* component of the book. A reader who follows any of the tracks above is encouraged to *run* the scripts of the études in their preferred language, modify a parameter (the resolution $N$, the map parameter $ell$, the boundary condition type), and re-render the figure. The Takeaway block at the end of each étude is the *minimum* one should carry forward; running the code makes the rest of the étude a learnable laboratory rather than an inert paragraph.
