# Computational Études: A Spectral Approach

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
[![Built with Typst](https://img.shields.io/badge/Built%20with-Typst-239DAD.svg)](https://typst.app/)

> A pedagogical textbook on spectral methods for differential equations, featuring dual-language implementations in Python and MATLAB.

**Author:** Dr. Denys Dutykh
**Affiliation:** Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE

---

## About This Book

Spectral methods are a powerful class of numerical techniques for solving differential equations. Unlike finite difference or finite element methods, which use local approximations, spectral methods represent solutions as linear combinations of global basis functions—typically Fourier series or Chebyshev polynomials. For smooth problems, this approach yields *spectral accuracy*: errors that decay exponentially fast with the number of degrees of freedom.

This book takes a hands-on, pedagogical approach inspired by musical *études*—short compositions designed to develop specific technical skills while remaining artistically complete. Each chapter focuses on a single mathematical concept and explores it through compact, runnable code implementations.

### Key Features

- **Étude-based pedagogy** — Each chapter is a self-contained study combining theory with implementation
- **Dual-language implementations** — All examples provided in both Python and MATLAB
- **Fully reproducible** — Every figure and result generated from the accompanying code
- **Focus on 1D problems** — Keeps code readable while covering all essential concepts
- **Beautiful typography** — Professionally typeset using Typst with bibliography backreferences

---

## Table of Contents

- **Preface** — Purpose, audience, and how to use this book
- **Acknowledgements** — Thanks to contributing students
1. **Introduction** — The spectral promise, philosophy of études, collocation methods, and modern workflows
2. **Classical Second Order PDEs and Separation of Variables** — Heat, wave, and Laplace equations; separation of variables as the foundation for spectral methods
3. **Mise en Bouche** — A first taste of spectral methods: method of weighted residuals, collocation vs. Galerkin with low-dimensional examples
4. **The Geometry of Nodes** — Runge phenomenon, potential theory, Chebyshev points, Lebesgue constants, and barycentric interpolation

*Additional chapters in development.*

---

## Getting Started

### Prerequisites

**For building the book:**
- [Typst](https://typst.app/) (modern typesetting system)

**For running Python code:**
- Python 3.8+
- NumPy
- SciPy
- Matplotlib

**For running MATLAB code (optional):**
- MATLAB R2020a or later
- [Advanpix Multiprecision Computing Toolbox](https://www.advanpix.com/) (optional, for extended precision)

### Building the Book

```bash
# Clone the repository
git clone https://github.com/dutykh/computational-etudes.git
cd computational-etudes

# Build the PDF textbook
make textbook

# Build the teaching plan
make tplan

# Build everything
make all
```

The compiled PDF will be available at `textbook/build/DD-Computational-Etudes-2026.pdf`.

### Running the Code

**Python:**
```bash
# Chapter 2: Classical PDEs
python codes/python/ch02_classical_pdes/heat_equation_evolution.py
python codes/python/ch02_classical_pdes/heat_equation_waterfall.py
python codes/python/ch02_classical_pdes/wave_equation_evolution.py
python codes/python/ch02_classical_pdes/wave_equation_waterfall.py
python codes/python/ch02_classical_pdes/laplace_equation_2d.py

# Chapter 3: Mise en Bouche
python codes/python/ch03_mise_en_bouche/collocation_example1.py
python codes/python/ch03_mise_en_bouche/collocation_vs_galerkin.py

# Chapter 4: The Geometry of Nodes
python codes/python/ch04_geometry_of_nodes/runge_phenomenon.py
python codes/python/ch04_geometry_of_nodes/chebyshev_success.py
python codes/python/ch04_geometry_of_nodes/chebyshev_points_circle.py
python codes/python/ch04_geometry_of_nodes/equipotential_curves.py
python codes/python/ch04_geometry_of_nodes/lagrange_basis.py
python codes/python/ch04_geometry_of_nodes/lebesgue_functions.py
python codes/python/ch04_geometry_of_nodes/convergence_comparison.py
```

**MATLAB:**
```matlab
% Navigate to the codes directory
cd codes/matlab/ch02_classical_pdes
heat_equation_evolution
heat_equation_waterfall
wave_equation_evolution
wave_equation_waterfall
laplace_equation_2d

cd ../ch03_mise_en_bouche
collocation_example1
collocation_vs_galerkin

cd ../ch04_geometry_of_nodes
runge_phenomenon
chebyshev_success
chebyshev_points_circle
equipotential_curves
lagrange_basis
lebesgue_functions
convergence_comparison
```

---

## Repository Structure

```
computational-etudes/
├── textbook/                    # Typst source for the textbook
│   ├── main.typ                 # Main entry point
│   ├── chapters/                # Chapter content
│   │   ├── preface.typ
│   │   ├── acknowledgements.typ
│   │   ├── introduction.typ
│   │   ├── classical_pdes.typ
│   │   ├── mise_en_bouche.typ
│   │   └── geometry_of_nodes.typ
│   ├── styles/                  # Typography and layout
│   │   └── template.typ
│   ├── biblio/                  # Bibliography
│   │   └── library.bib
│   ├── figures/                 # Generated figures
│   │   ├── ch02/
│   │   │   ├── python/
│   │   │   └── matlab/
│   │   ├── ch03/
│   │   │   ├── python/
│   │   │   └── matlab/
│   │   └── ch04/
│   │       ├── python/
│   │       └── matlab/
│   └── build/                   # Compiled PDF output
├── codes/                       # Code implementations
│   ├── python/
│   │   ├── ch02_classical_pdes/
│   │   ├── ch03_mise_en_bouche/
│   │   └── ch04_geometry_of_nodes/
│   ├── matlab/
│   │   ├── ch02_classical_pdes/
│   │   ├── ch03_mise_en_bouche/
│   │   └── ch04_geometry_of_nodes/
│   └── README.md
├── tplan/                       # Teaching plan (MATH 794)
│   ├── teaching_plan.typ
│   ├── Makefile
│   └── build/
├── Makefile                     # Build automation
├── CLAUDE.md                    # Project conventions
├── LICENSE                      # CC BY-NC-SA 4.0
└── README.md
```

---

## Typst Packages Used

The textbook uses the following Typst packages:
- **[codly](https://typst.app/universe/package/codly)** — Beautiful code blocks with syntax highlighting
- **[retrofit](https://typst.app/universe/package/retrofit)** — Bibliography backreferences showing citation locations

---

## Teaching Materials

This book is used for **MATH 794** at Khalifa University. A tentative teaching plan is available in the `tplan/` directory, which includes:
- Weekly lecture schedule
- Chapter-to-lecture mapping
- Progress tracking

Build the teaching plan:
```bash
make tplan
```

---

## Citation

If you use this book in your research or teaching, please cite it as:

```bibtex
@book{dutykh2026etudes,
  author    = {Dutykh, Denys},
  title     = {Computational Études: A Spectral Approach},
  year      = {2026},
  publisher = {Self-published},
  url       = {https://github.com/dutykh/computational-etudes}
}
```

---

## License

This work is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](https://creativecommons.org/licenses/by-nc-sa/4.0/).

You are free to:
- **Share** — copy and redistribute the material in any medium or format
- **Adapt** — remix, transform, and build upon the material

Under the following terms:
- **Attribution** — You must give appropriate credit to the author
- **NonCommercial** — You may not use the material for commercial purposes
- **ShareAlike** — If you remix or transform the material, you must distribute your contributions under the same license

---

## Contributing

Contributions are welcome! If you find errors, have suggestions, or want to contribute code examples, please:

1. Open an issue to discuss your proposed changes
2. Fork the repository
3. Create a pull request with your improvements

---

## Contact

For questions or feedback, please open an issue on GitHub or contact the author through the university.

---

*This book is a work in progress. Check back for updates as new chapters are added.*
