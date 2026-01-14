# Computational Études: Code Repository

This directory contains the Python and MATLAB implementations accompanying the book *Computational Études: A Spectral Approach* by Dr. Denys Dutykh.

## Directory Structure

```
codes/
├── python/
│   └── ch02_classical_pdes/           # Chapter 2: Classical PDEs
│       ├── heat_equation_evolution.py  # Heat equation time evolution
│       ├── wave_equation_evolution.py  # Wave equation oscillations
│       └── laplace_equation_2d.py      # Laplace equation in 2D strip
└── matlab/
    └── ch02_classical_pdes/
        ├── heat_equation_evolution.m
        ├── wave_equation_evolution.m
        └── laplace_equation_2d.m
```

## Requirements

### Python

- Python 3.8+
- NumPy
- Matplotlib

Install dependencies:
```bash
pip install numpy matplotlib
```

### MATLAB

- MATLAB R2020a or later (for `exportgraphics` function)
- No additional toolboxes required for basic examples

## Running the Codes

### Python

From the repository root:
```bash
# Heat equation
python codes/python/ch02_classical_pdes/heat_equation_evolution.py

# Wave equation
python codes/python/ch02_classical_pdes/wave_equation_evolution.py

# Laplace equation
python codes/python/ch02_classical_pdes/laplace_equation_2d.py
```

### MATLAB

From MATLAB, navigate to the script directory and run:
```matlab
cd codes/matlab/ch02_classical_pdes
heat_equation_evolution
wave_equation_evolution
laplace_equation_2d
```

Or add the path and run:
```matlab
addpath('codes/matlab/ch02_classical_pdes')
heat_equation_evolution
wave_equation_evolution
laplace_equation_2d
```

## Output

Figures are saved in `textbook/figures/` organized by chapter and language:

```
textbook/figures/
└── ch02/
    ├── python/          # Python-generated figures (used in published textbook)
    │   ├── heat_evolution.pdf
    │   ├── wave_evolution.pdf
    │   └── laplace_solution.pdf
    └── matlab/          # MATLAB-generated figures
        ├── heat_evolution.pdf
        ├── wave_evolution.pdf
        └── laplace_solution.pdf
```

The **Python figures** are used in the published textbook. MATLAB figures are provided for users who prefer that environment.

## Chapter 2: Classical PDEs

The codes in `ch02_classical_pdes/` visualize the analytical solutions derived in Chapter 2:

| Script | PDE | Description |
|--------|-----|-------------|
| `heat_equation_evolution` | Heat equation | Time evolution of triangle wave initial condition showing smoothing effect |
| `wave_equation_evolution` | Wave equation | Oscillation of plucked string showing standing wave patterns |
| `laplace_equation_2d` | Laplace equation | 2D harmonic function in periodic strip showing mode decay |

## Reproducibility

Every figure in the book can be regenerated from the codes in this directory. To regenerate all figures:

```bash
# Python (from repository root)
find codes/python -name "*.py" -exec python {} \;

# MATLAB (from MATLAB command window)
% Run each .m file in the codes/matlab directory
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
