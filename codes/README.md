# Computational Études: Code Repository

This directory contains the Python and MATLAB implementations accompanying the book *Computational Études: A Spectral Approach* by Dr. Denys Dutykh.

## Directory Structure

```
codes/
├── python/
│   ├── ch02_classical_pdes/           # Chapter 2: Classical PDEs
│   │   ├── heat_equation_evolution.py  # Heat equation time evolution
│   │   ├── heat_equation_waterfall.py  # Heat equation 3D waterfall plot
│   │   ├── wave_equation_evolution.py  # Wave equation oscillations
│   │   ├── wave_equation_waterfall.py  # Wave equation 3D waterfall plot
│   │   └── laplace_equation_2d.py      # Laplace equation in 2D strip
│   └── ch03_mise_en_bouche/           # Chapter 3: Mise en Bouche
│       ├── collocation_example1.py     # Three-coefficient collocation example
│       └── collocation_vs_galerkin.py  # Comparison of collocation and Galerkin methods
└── matlab/
    ├── ch02_classical_pdes/
    │   ├── heat_equation_evolution.m
    │   ├── heat_equation_waterfall.m
    │   ├── wave_equation_evolution.m
    │   ├── wave_equation_waterfall.m
    │   └── laplace_equation_2d.m
    └── ch03_mise_en_bouche/
        ├── collocation_example1.m
        └── collocation_vs_galerkin.m
```

## Requirements

### Python

- Python 3.8+
- NumPy
- SciPy
- Matplotlib

Install dependencies:
```bash
pip install numpy scipy matplotlib
```

### MATLAB

- MATLAB R2020a or later (for `exportgraphics` function)
- No additional toolboxes required for basic examples

## Running the Codes

### Python

From the repository root:
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
```

### MATLAB

From MATLAB, navigate to the script directory and run:
```matlab
cd codes/matlab/ch02_classical_pdes
heat_equation_evolution
heat_equation_waterfall
wave_equation_evolution
wave_equation_waterfall
laplace_equation_2d

cd ../ch03_mise_en_bouche
collocation_example1
collocation_vs_galerkin
```

Or add the path and run:
```matlab
addpath('codes/matlab/ch02_classical_pdes')
addpath('codes/matlab/ch03_mise_en_bouche')
heat_equation_evolution
collocation_example1
```

## Output

Figures are saved in `textbook/figures/` organized by chapter and language:

```
textbook/figures/
├── ch02/
│   ├── python/               # Python-generated figures (used in published textbook)
│   │   ├── heat_evolution.pdf
│   │   ├── heat_waterfall.pdf
│   │   ├── wave_evolution.pdf
│   │   ├── wave_waterfall.pdf
│   │   └── laplace_solution.pdf
│   └── matlab/               # MATLAB-generated figures
│       ├── heat_evolution.pdf
│       ├── heat_waterfall.pdf
│       ├── wave_evolution.pdf
│       ├── wave_waterfall.pdf
│       └── laplace_solution.pdf
└── ch03/
    ├── python/
    │   ├── collocation_example1.pdf
    │   └── collocation_vs_galerkin.pdf
    └── matlab/
        ├── collocation_example1.pdf
        └── collocation_vs_galerkin.pdf
```

The **Python figures** are used in the published textbook. MATLAB figures are provided for users who prefer that environment.

## Chapter 2: Classical PDEs

The codes in `ch02_classical_pdes/` visualize the analytical solutions derived in Chapter 2:

| Script | PDE | Description |
|--------|-----|-------------|
| `heat_equation_evolution` | Heat equation | Time evolution of triangle wave initial condition showing smoothing effect |
| `heat_equation_waterfall` | Heat equation | 3D waterfall visualization of temperature evolution |
| `wave_equation_evolution` | Wave equation | Oscillation of plucked string showing standing wave patterns |
| `wave_equation_waterfall` | Wave equation | 3D waterfall visualization of wave dynamics |
| `laplace_equation_2d` | Laplace equation | 2D harmonic function in periodic strip showing mode decay |

## Chapter 3: Mise en Bouche

The codes in `ch03_mise_en_bouche/` implement the low-dimensional spectral methods introduced in Chapter 3:

| Script | Method | Description |
|--------|--------|-------------|
| `collocation_example1` | Collocation | Three-coefficient polynomial approximation for a BVP with u'' - (4x² + 2)u = 0 |
| `collocation_vs_galerkin` | Both | Comparison of collocation and Galerkin methods for a reaction-diffusion problem |

## Reproducibility

Every figure in the book can be regenerated from the codes in this directory. To regenerate all figures:

```bash
# Using Make (from repository root)
make figures-python    # Generate all Python figures
make figures-matlab    # Generate all MATLAB figures
make figures          # Generate both

# Or manually (Python)
find codes/python -name "*.py" -exec python {} \;
```

## Code Style

- **Python**: PEP 8 compliant, NumPy/SciPy ecosystem
- **MATLAB**: Vectorized operations, clear section comments

Both implementations are designed to be readable and educational, prioritizing clarity over brevity.

## Author

Dr. Denys Dutykh
Mathematics Department
Khalifa University of Science and Technology
Abu Dhabi, UAE
