// main.typ
#import "styles/template.typ": project, dropcap

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

#include "chapters/afterword.typ"

// --- Bibliography ---
#bibliography("biblio/library.bib")