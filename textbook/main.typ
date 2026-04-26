// main.typ
#import "styles/template.typ": project, dropcap, etude-conclusion, idx

#show: project.with(
  title: "Computational Études",
  subtitle: "A Spectral Approach",
  author: "Dr. Denys Dutykh",
  affiliation: [
    Mathematics Department \
    College of Computing and Mathematical Sciences \
    Khalifa University of Science and Technology \
    Abu Dhabi, UAE
  ],
  date: datetime.today().display("[year]"),
)

// --- Book Content ---

#include "chapters/preface.typ"

#include "chapters/acknowledgements.typ"

#include "chapters/roadmap.typ"

#include "chapters/notation.typ"

// Reset counter for main chapters if needed, though the template handles page numbering
#include "chapters/introduction.typ"

#include "chapters/classical_pdes.typ"

#include "chapters/mise_en_bouche.typ"

#include "chapters/geometry_of_nodes.typ"

#include "chapters/differentiation_matrices.typ"

#include "chapters/smoothness_accuracy.typ"

#include "chapters/chebyshev_differentiation.typ"

#include "chapters/boundary_value_problems.typ"

#include "chapters/fourier_grids.typ"

#include "chapters/spectral_pde_solvers.typ"

#include "chapters/fourier_pseudospectral.typ"

#include "chapters/polar_coordinates.typ"

#include "chapters/advanced_boundary_conditions.typ"

#include "chapters/higher_order_bvps.typ"

#include "chapters/quadrature.typ"

#include "chapters/periodic_quadrature.typ"

#include "chapters/time_stability.typ"

#include "chapters/linear_eigenproblems.typ"

#include "chapters/coordinate_transformations.typ"

#include "chapters/unbounded_intervals.typ"

#include "chapters/special_tricks.typ"

#include "chapters/afterword.typ"

// --- Appendices ---
// The Springer-style numbering and "APPENDIX" banner override are scoped
// to the content block below; the bibliography heading that follows
// reverts to the body-chapter rules.
#{
  let _navy = rgb(20, 45, 110)
  let _sky = rgb(120, 150, 210)
  [
    #counter(heading).update(0)
    #counter(figure.where(kind: image)).update(0)
    #counter(figure.where(kind: table)).update(0)
    #counter(math.equation).update(0)

    #set heading(numbering: "A.1")

    #set math.equation(numbering: n => {
      let chapters = counter(heading).get()
      if chapters.len() > 0 and chapters.first() > 0 {
        numbering("(A.1)", chapters.first(), n)
      } else {
        numbering("(1)", n)
      }
    })

    #set figure(numbering: n => {
      let chapters = counter(heading).get()
      if chapters.len() > 0 and chapters.first() > 0 {
        numbering("A.1", chapters.first(), n)
      } else {
        numbering("1", n)
      }
    })

    #show heading.where(level: 1): it => {
      pagebreak(weak: true)
      counter(math.equation).update(0)
      counter(figure.where(kind: image)).update(0)
      counter(figure.where(kind: table)).update(0)
      let number = if it.numbering == none {
        none
      } else {
        counter(heading).display(it.numbering)
      }
      v(3cm)
      grid(
        columns: (6mm, 1fr),
        gutter: 5mm,
        rect(width: 100%, height: 100%, fill: _navy),
        block[
          #if number != none [
            #text(size: 1.2em, weight: "bold", fill: _sky)[APPENDIX #number]
            #v(0.3em)
          ]
          #set text(hyphenate: false)
          #text(size: 2.0em, weight: 700, fill: _navy, it.body)
        ],
      )
      v(1.5cm)
    }

    #include "chapters/appendix_a_linalg.typ"
    #include "chapters/appendix_b_complex.typ"
    #include "chapters/appendix_c_software.typ"
    #include "chapters/appendix_d_acronyms.typ"
  ]
}

// --- Bibliography ---
#bibliography("biblio/library.bib")

// --- Subject Index ---
#include "chapters/subject_index.typ"

// --- Colophon (blurb, keywords, MSC, About the Author) ---
#include "chapters/colophon.typ"
