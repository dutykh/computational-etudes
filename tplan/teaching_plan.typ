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

#figure(
  table(
    columns: (auto, auto, 1fr, 1fr),
    align: (center, center, left, left),
    inset: 0.7em,

    // Header
    table.cell(fill: accent-color)[#text(fill: white, weight: "bold")[Week]],
    table.cell(fill: accent-color)[#text(fill: white, weight: "bold")[Dates]],
    table.cell(fill: accent-color)[#text(fill: white, weight: "bold")[Monday Lecture]],
    table.cell(fill: accent-color)[#text(fill: white, weight: "bold")[Wednesday Lecture]],

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
      #text(size: 9pt, fill: luma(100))[Method of weighted residuals, collocation example, collocation vs Galerkin]
    ],
    table.cell(fill: completed-color)[
      *Ch. 4:* The Geometry of Nodes \
      #text(size: 9pt, fill: luma(100))[Lagrange interpolation, Runge phenomenon, potential theory]
    ],

    // Week 3 (Jan 26-29)
    table.cell(fill:current-color)[*3*], table.cell(fill:current-color)[Jan 26--28],
    table.cell(fill:current-color)[], table.cell(fill:current-color)[],
    
    // Week 4 (Feb 2-5)
    [*4*], [Feb 2--4],
    [], [],

    // Week 5 (Feb 9-12)
    [*5*], [Feb 9--11],
    [], [],

    // Week 6 (Feb 16-19)
    [*6*], [Feb 16--18],
    [], [],

    // Week 7 (Feb 23-26)
    [*7*], [Feb 23--25],
    [], [],

    // Week 8 (Mar 2-5)
    [*8*], [Mar 2--4],
    [], [],

    // Week 9 (Mar 9-12)
    [*9*], [Mar 9--11],
    [], [],

    // Spring Break: Mar 16-27 (includes Eid Al Fitr Mar 19-20)
    table.cell(fill: luma(240))[--],
    table.cell(fill: luma(240))[Mar 16--27],
    table.cell(fill: luma(240), colspan: 2)[
      #align(center)[
        #text(style: "italic", fill: luma(100))[Spring Break & Eid Al Fitr -- No classes]
      ]
    ],

    // Week 10 (Mar 30 - Apr 2)
    [*10*], [Mar 30--Apr 1],
    [], [],

    // Week 11 (Apr 6-9)
    [*11*], [Apr 6--8],
    [], [],

    // Week 12 (Apr 13-16)
    [*12*], [Apr 13--15],
    [], [],

    // Week 13 (Apr 20-23)
    [*13*], [Apr 20--22],
    [], [],

    // Week 14 (Apr 27-30)
    [*14*], [Apr 27--29],
    [], [],

    // Finals Week (May 4-14)
    table.cell(fill: luma(240))[--],
    table.cell(fill: luma(240))[May 4--14],
    table.cell(fill: luma(240), colspan: 2)[
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
]

#v(1em)

#grid(
  columns: (1fr, 1fr),
  gutter: 1em,
  block(
    width: 100%,
    inset: 1em,
    radius: 6pt,
    stroke: 1pt + luma(200),
    fill: luma(250),
  )[
    #text(weight: "bold", fill: accent-color)[Project II]
    #v(0.3em)
    #set text(size: 10pt, style: "italic", fill: luma(100))
    To be announced
  ],
  block(
    width: 100%,
    inset: 1em,
    radius: 6pt,
    stroke: 1pt + luma(200),
    fill: luma(250),
  )[
    #text(weight: "bold", fill: accent-color)[Project III]
    #v(0.3em)
    #set text(size: 10pt, style: "italic", fill: luma(100))
    To be announced
  ],
)

#v(1em)

#block(
  width: 50%,
  inset: 1em,
  radius: 6pt,
  stroke: 1pt + luma(200),
  fill: luma(250),
)[
  #text(weight: "bold", fill: accent-color)[Project IV]
  #v(0.3em)
  #set text(size: 10pt, style: "italic", fill: luma(100))
  To be announced
]

#v(2em)

#align(right)[
  #set text(size: 10pt, fill: luma(100))
  _Last updated: #datetime.today().display("[month repr:long] [day], [year]")_
]
