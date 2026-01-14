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
- **Beautiful typography** — Professionally typeset using Typst

---

## Table of Contents

1. **Introduction** — The spectral promise, philosophy of études, and modern workflows
2. **Classical Second Order PDEs and Separation of Variables** — Heat, wave, and Laplace equations; separation of variables as the foundation for spectral methods

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

# Build the PDF
make textbook
```

The compiled PDF will be available at `textbook/build/DD-Computational-Etudes.pdf`.

### Running the Code

**Python:**
```bash
cd python
python chapter01_introduction.py
```

**MATLAB:**
```matlab
cd matlab
chapter01_introduction
```

---

## Repository Structure

```
computational-etudes/
├── textbook/                # Typst source for the textbook
│   ├── main.typ             # Main entry point
│   ├── chapters/            # Chapter content
│   │   ├── preface.typ
│   │   └── introduction.typ
│   ├── styles/              # Typography and layout
│   │   └── template.typ
│   └── build/               # Compiled PDF output
├── codes/                   # Code implementations
│   ├── python/              # Python implementations
│   └── matlab/              # MATLAB implementations
├── Makefile                 # Build automation
├── LICENSE                  # CC BY-NC-SA 4.0
└── README.md
```

---

## Citation

If you use this book in your research or teaching, please cite it as:

```bibtex
@book{dutykh2025etudes,
  author    = {Dutykh, Denys},
  title     = {Computational Études: A Spectral Approach},
  year      = {2025},
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
