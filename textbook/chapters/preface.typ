// textbook/chapters/preface.typ
// Preface
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Email: denys.dutykh@ku.ac.ae
// Homepage: https://www.denys-dutykh.com/
// Last modified: February 2026
#import "../styles/template.typ": dropcap

#heading(level: 1, numbering: none)[Preface]

#dropcap[The purpose of this book is twofold: to explain to students why spectral methods work and why they are so remarkably efficient, and to teach them how to implement these methods using modern computational tools.]

This text is conceived as a collection of *Computational Études*. In musical education, an étude is a composition designed to perfect a specific technique (be it rapid scales, complex arpeggios, or delicate phrasing) while remaining a pleasing piece of music in its own right. Similarly, each chapter here presents a focused mathematical concept paired with its computational realization: a study that is both instructive and complete.

#heading(level: 3, numbering: none)[Who Is This Book For?]

We have written primarily for graduate students in applied mathematics, physics, and engineering who seek a hands-on understanding of spectral methods. The reader should be comfortable with calculus, linear algebra, and basic programming. Prior exposure to numerical methods is helpful but not essential; we build the necessary foundations as we proceed. In short, this is a book for the classroom and the laptop, not a monograph for the specialist's shelf, a distinction the next section makes explicit.

#heading(level: 3, numbering: none)[What This Book Is, and What It Is Not]

A book is best read with its purpose in mind, so we state ours plainly. _Computational Études_ is a graduate classroom textbook. It is written to be read, and run, from cover to cover by a student who does not yet have a spectral-methods specialist looking over their shoulder.

*This book is.* A hands-on, classroom-tested introduction for interdisciplinary graduate students, in fluid dynamics, solid mechanics, electromagnetics, geophysics, and beyond, who each need a shared, workable foundation before extracting the one technique their own thesis requires. It favors _breadth_ of working technique over exhaustive depth: a panorama of methods, each delivered as a short, runnable étude in Python, Julia, and MATLAB. It carries the philosophy of Trefethen, Boyd, and Fornberg into a modern computational stack, so that every script runs, unchanged, on the reader's own laptop. And it lingers, deliberately, on the small practical things (why a differentiation matrix refuses to converge, how to place collocation nodes, where aliasing bites) that save a newcomer days of debugging.

*This book is not.* A research monograph, nor an advanced reference for specialists who already practise spectral methods. It is not an exhaustive or fully rigorous treatise: it favors clarity and working intuition over complete theorems and convergence proofs, and it stays, by design, in one space dimension, where the essential ideas already live (a choice we explain in the Introduction). For the complete theory, the state of the art, and the multi-dimensional and complex-geometry extensions, the reader is in far better hands with the standard references: Boyd @Boyd2000, Canuto et al. @Canuto2006, Shen, Tang, and Wang @ShenTangWang2011, and Trefethen @Trefethen2000.

In short: if you are a specialist seeking the frontier of the field, this is not the book you need, and the references above will serve you better. If you are a graduate student who must get a spectral solver working this semester, and understand _why_ it works, you are in exactly the right place.

#heading(level: 3, numbering: none)[How Is the Book Organized?]

Each chapter is designed to be largely self-contained. We begin with fundamental concepts (interpolation and differentiation) before advancing to time-stepping schemes and applications. The mathematical exposition is deliberately concise, favoring clarity over exhaustive rigor. Proofs are included when they illuminate; otherwise, we direct the reader to authoritative references.

#heading(level: 3, numbering: none)[Intellectual Debt]

This book owes much to Lloyd N. Trefethen's _Spectral Methods in MATLAB_ @Trefethen2000, which some twenty years ago served as my own gateway into the world of spectral methods. The reader will recognize its genetic imprint throughout: in the emphasis on eigenvalue portraits as diagnostic tools, in the primacy of Chebyshev grids, and in the conviction that a well-chosen fifty-line script can illuminate an idea more vividly than ten pages of formal proof. Where Trefethen's text was a masterclass in MATLAB-driven numerical exploration, the present volume attempts to carry the same spirit forward --- updating the computational toolkit, broadening the range of applications, and adding the layer of mathematical context that turns a recipe into understanding. Anyone who has benefited from Trefethen's book should find themselves on familiar ground here; those who have not yet read it are warmly encouraged to do so.

A second debt, less visible at the outset of the book but very heavy from the eighteenth chapter onwards, is owed to John P. Boyd's encyclopaedic _Chebyshev and Fourier Spectral Methods_ @Boyd2000. Boyd's book is the closest thing the field has to a complete reference: it covers basis choice, mapping strategies, behavioural boundary conditions, the analysis of singular and singularly-perturbed problems, the special tricks that turn a stalled computation into an elegant one, and the symbolic small-$N$ regime in which the usual numerical instincts must be reversed. The chapters of the present text on linear eigenproblems, coordinate transformations, unbounded intervals, and the closing bonus chapter on special tricks all draw heavily on Boyd's exposition; many of the études in those chapters are computational realisations of examples Boyd himself describes. Where Trefethen taught me _how_ to write a spectral solver, Boyd taught me _why_ a particular spectral solver is the right one for a particular problem.

#heading(level: 3, numbering: none)[Reproducible Science]

This book is also an experiment in *reproducible science*. Every figure, every table, every numerical result you see in these pages was generated by code available in the accompanying repository. We provide implementations in three languages --- Python, MATLAB, and Julia --- allowing readers to choose their preferred environment. The Python code emphasizes accessibility and integration with the open-source scientific ecosystem; the MATLAB code leverages its historical significance in numerical computing and, where appropriate, the Advanpix Multiprecision Computing Toolbox for extended precision arithmetic; the Julia code offers a modern alternative that combines the readability of a high-level language with performance approaching that of compiled code.

#heading(level: 3, numbering: none)[How to Use This Book]

We invite you to treat this book not as a static reference, but as a workshop. Clone the repository, run the scripts, modify the parameters, break the code, and fix it. That is the only way to truly master the art of spectral methods.

The complete source code and the Typst manuscript are available at: \
#link("https://github.com/dutykh/computational-etudes/")[https://github.com/dutykh/computational-etudes/]

#v(3em)

#block(width: 100%)[
  *Dr. Denys Dutykh* #h(1fr) Abu Dhabi, UAE
  #align(right)[#datetime.today().display("[month repr:long] [year]")]
]
