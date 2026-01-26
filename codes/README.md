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
│   ├── ch03_mise_en_bouche/           # Chapter 3: Mise en Bouche
│   │   ├── collocation_example1.py     # Three-coefficient collocation example
│   │   └── collocation_vs_galerkin.py  # Comparison of collocation and Galerkin methods
│   ├── ch04_geometry_of_nodes/        # Chapter 4: The Geometry of Nodes
│   │   ├── runge_phenomenon.py         # Runge phenomenon visualization
│   │   ├── chebyshev_success.py        # Chebyshev interpolation success
│   │   ├── chebyshev_points_circle.py  # Geometric construction of Chebyshev points
│   │   ├── equipotential_curves.py     # Potential theory equipotential curves
│   │   ├── lagrange_basis.py           # Lagrange basis functions comparison
│   │   ├── lebesgue_functions.py       # Lebesgue functions and constants
│   │   └── convergence_comparison.py   # Convergence rate comparison
│   ├── ch05_differentiation_matrices/  # Chapter 5: Differentiation Matrices
│   │   ├── fdweights.py                # Fornberg's algorithm for FD weights
│   │   ├── spectral_matrix_periodic.py # Periodic spectral differentiation matrix
│   │   ├── fd_matrix_bandwidth.py      # FD matrix sparsity visualization
│   │   ├── spectral_matrix_structure.py # Spectral matrix structure visualization
│   │   ├── stencil_pyramid.py          # Fornberg recursion pyramid diagram
│   │   ├── periodic_cardinal_functions.py # Periodic cardinal function visualization
│   │   ├── spectral_derivatives_demo.py # Spectral differentiation demonstration
│   │   ├── higher_order_derivatives.py # Higher-order derivatives (up to 4th)
│   │   └── convergence_comparison.py   # FD vs spectral convergence comparison
│   ├── ch06_chebyshev_differentiation/ # Chapter 6: Chebyshev Differentiation
│   │   ├── cheb_matrix.py              # Chebyshev differentiation matrix construction
│   │   ├── cheb_diff_demo.py           # Chebyshev differentiation demonstration
│   │   ├── cheb_convergence.py         # Convergence analysis (waterfall plot)
│   │   ├── cheb_grid_comparison.py     # Periodic vs Chebyshev grid comparison
│   │   ├── cheb_matrix_structure.py    # Chebyshev matrix structure visualization
│   │   └── cheb_cardinal.py            # Chebyshev cardinal function visualization
│   └── ch07_boundary_value_problems/   # Chapter 7: Boundary Value Problems
│       ├── bvp_linear.py               # 1D Poisson equation with spectral collocation
│       ├── bvp_variable_coeff.py       # Variable coefficient BVP
│       ├── bvp_eigenvalue.py           # Eigenvalue problems (quantum harmonic oscillator)
│       ├── bvp_2d_poisson.py           # 2D Poisson equation with tensor product grids
│       ├── bvp_helmholtz.py            # Helmholtz equation
│       └── bvp_nonlinear.py            # Nonlinear BVP (Bratu equation)
└── matlab/
    ├── ch02_classical_pdes/
    │   ├── heat_equation_evolution.m
    │   ├── heat_equation_waterfall.m
    │   ├── wave_equation_evolution.m
    │   ├── wave_equation_waterfall.m
    │   └── laplace_equation_2d.m
    ├── ch03_mise_en_bouche/
    │   ├── collocation_example1.m
    │   └── collocation_vs_galerkin.m
    ├── ch04_geometry_of_nodes/
    │   ├── runge_phenomenon.m
    │   ├── chebyshev_success.m
    │   ├── chebyshev_points_circle.m
    │   ├── equipotential_curves.m
    │   ├── lagrange_basis.m
    │   ├── lebesgue_functions.m
    │   └── convergence_comparison.m
    ├── ch05_differentiation_matrices/
    │   ├── fdweights.m                 # Fornberg's algorithm for FD weights
    │   ├── fd_matrix_periodic.m        # Periodic FD matrix construction
    │   ├── spectral_matrix_periodic.m  # Periodic spectral differentiation matrix
    │   ├── fd_matrix_bandwidth.m       # FD matrix sparsity visualization
    │   ├── spectral_matrix_structure.m # Spectral matrix structure visualization
    │   ├── stencil_pyramid.m           # Fornberg recursion pyramid diagram
    │   ├── periodic_cardinal_functions.m # Periodic cardinal function visualization
    │   ├── spectral_derivatives_demo.m # Spectral differentiation demonstration
    │   ├── higher_order_derivatives.m  # Higher-order derivatives (up to 4th)
    │   └── convergence_comparison.m    # FD vs spectral convergence comparison
    ├── ch06_chebyshev_differentiation/
    │   ├── cheb_matrix.m               # Chebyshev differentiation matrix construction
    │   ├── verify_cheb_matrix.m        # Verification of Chebyshev matrix accuracy
    │   ├── cheb_grid_comparison.m      # Equispaced vs Chebyshev grid comparison
    │   ├── cheb_matrix_structure.m     # Matrix heatmap and row profiles
    │   ├── cheb_cardinal.m             # Cardinal functions visualization
    │   ├── cheb_diff_demo.m            # Differentiation demo with Witch of Agnesi
    │   └── cheb_convergence.m          # Convergence rates for test functions
    └── ch07_boundary_value_problems/
        ├── bvp_linear.m                # 1D Poisson equation solver
        ├── bvp_variable_coeff.m        # Variable coefficient BVP
        ├── bvp_eigenvalue.m            # Eigenvalue problem solver
        ├── bvp_2d_poisson.m            # 2D Poisson with tensor product grids
        ├── bvp_helmholtz.m             # Helmholtz equation near resonance
        └── bvp_nonlinear.m             # Bratu equation with Newton iteration
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

# Chapter 4: The Geometry of Nodes
python codes/python/ch04_geometry_of_nodes/runge_phenomenon.py
python codes/python/ch04_geometry_of_nodes/chebyshev_success.py
python codes/python/ch04_geometry_of_nodes/chebyshev_points_circle.py
python codes/python/ch04_geometry_of_nodes/equipotential_curves.py
python codes/python/ch04_geometry_of_nodes/lagrange_basis.py
python codes/python/ch04_geometry_of_nodes/lebesgue_functions.py
python codes/python/ch04_geometry_of_nodes/convergence_comparison.py

# Chapter 5: Differentiation Matrices
python codes/python/ch05_differentiation_matrices/fd_matrix_bandwidth.py
python codes/python/ch05_differentiation_matrices/spectral_matrix_structure.py
python codes/python/ch05_differentiation_matrices/stencil_pyramid.py
python codes/python/ch05_differentiation_matrices/periodic_cardinal_functions.py
python codes/python/ch05_differentiation_matrices/spectral_derivatives_demo.py
python codes/python/ch05_differentiation_matrices/higher_order_derivatives.py
python codes/python/ch05_differentiation_matrices/convergence_comparison.py

# Chapter 6: Chebyshev Differentiation Matrices
python codes/python/ch06_chebyshev_differentiation/cheb_matrix.py
python codes/python/ch06_chebyshev_differentiation/cheb_diff_demo.py
python codes/python/ch06_chebyshev_differentiation/cheb_convergence.py
python codes/python/ch06_chebyshev_differentiation/cheb_grid_comparison.py
python codes/python/ch06_chebyshev_differentiation/cheb_matrix_structure.py
python codes/python/ch06_chebyshev_differentiation/cheb_cardinal.py

# Chapter 7: Boundary Value Problems
python codes/python/ch07_boundary_value_problems/bvp_linear.py
python codes/python/ch07_boundary_value_problems/bvp_variable_coeff.py
python codes/python/ch07_boundary_value_problems/bvp_eigenvalue.py
python codes/python/ch07_boundary_value_problems/bvp_2d_poisson.py
python codes/python/ch07_boundary_value_problems/bvp_helmholtz.py
python codes/python/ch07_boundary_value_problems/bvp_nonlinear.py
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

cd ../ch04_geometry_of_nodes
runge_phenomenon
chebyshev_success
chebyshev_points_circle
equipotential_curves
lagrange_basis
lebesgue_functions
convergence_comparison

cd ../ch05_differentiation_matrices
fd_matrix_bandwidth
spectral_matrix_structure
stencil_pyramid
periodic_cardinal_functions
spectral_derivatives_demo
higher_order_derivatives
convergence_comparison

cd ../ch06_chebyshev_differentiation
cheb_matrix
verify_cheb_matrix
cheb_grid_comparison
cheb_matrix_structure
cheb_cardinal
cheb_diff_demo
cheb_convergence

cd ../ch07_boundary_value_problems
bvp_linear
bvp_variable_coeff
bvp_eigenvalue
bvp_2d_poisson
bvp_helmholtz
bvp_nonlinear
```

Or add the path and run:
```matlab
addpath('codes/matlab/ch02_classical_pdes')
addpath('codes/matlab/ch03_mise_en_bouche')
addpath('codes/matlab/ch04_geometry_of_nodes')
addpath('codes/matlab/ch05_differentiation_matrices')
addpath('codes/matlab/ch06_chebyshev_differentiation')
addpath('codes/matlab/ch07_boundary_value_problems')
heat_equation_evolution
collocation_example1
runge_phenomenon
fd_matrix_bandwidth
cheb_matrix
bvp_linear
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
├── ch03/
│   ├── python/
│   │   ├── collocation_example1.pdf
│   │   └── collocation_vs_galerkin.pdf
│   └── matlab/
│       ├── collocation_example1.pdf
│       └── collocation_vs_galerkin.pdf
├── ch04/
│   ├── python/
│   │   ├── runge_phenomenon.pdf
│   │   ├── chebyshev_success.pdf
│   │   ├── chebyshev_points_circle.pdf
│   │   ├── equipotential_curves.pdf
│   │   ├── lagrange_basis.pdf
│   │   ├── lebesgue_functions.pdf
│   │   └── convergence_comparison.pdf
│   └── matlab/
│       ├── runge_phenomenon.pdf
│       ├── chebyshev_success.pdf
│       ├── chebyshev_points_circle.pdf
│       ├── equipotential_curves.pdf
│       ├── lagrange_basis.pdf
│       ├── lebesgue_functions.pdf
│       └── convergence_comparison.pdf
├── ch05/
│   ├── python/
│   │   ├── fd_stencil_schematic.pdf
│   │   ├── fd_matrix_bandwidth.pdf
│   │   ├── spectral_matrix_structure.pdf
│   │   ├── stencil_pyramid.pdf
│   │   ├── periodic_cardinal_functions.pdf
│   │   ├── spectral_derivatives_demo.pdf
│   │   ├── higher_order_derivatives.pdf
│   │   ├── d2_comparison.pdf
│   │   └── convergence_comparison.pdf
│   └── matlab/
│       ├── fd_stencil_schematic.pdf
│       ├── fd_matrix_bandwidth.pdf
│       ├── spectral_matrix_structure.pdf
│       ├── stencil_pyramid.pdf
│       ├── periodic_cardinal_functions.pdf
│       ├── spectral_derivatives_demo.pdf
│       ├── higher_order_derivatives.pdf
│       ├── d2_comparison.pdf
│       └── convergence_comparison.pdf
├── ch06/
│   ├── python/
│   │   ├── cheb_diff_demo.pdf
│   │   ├── convergence_waterfall.pdf
│   │   ├── grid_comparison.pdf
│   │   ├── cheb_matrix_structure.pdf
│   │   └── cheb_cardinal.pdf
│   └── matlab/
│       ├── cheb_diff_demo.pdf
│       ├── convergence_waterfall.pdf
│       ├── grid_comparison.pdf
│       ├── cheb_matrix_structure.pdf
│       └── cheb_cardinal.pdf
└── ch07/
    ├── python/
    │   ├── poisson_1d.pdf
    │   ├── variable_coeff.pdf
    │   ├── eigenvalue_problem.pdf
    │   ├── tensor_grid.pdf
    │   ├── poisson_2d.pdf
    │   ├── laplacian_sparsity.pdf
    │   ├── helmholtz.pdf
    │   └── bratu.pdf
    └── matlab/
        ├── poisson_1d.pdf
        ├── variable_coeff.pdf
        ├── eigenvalue_problem.pdf
        ├── tensor_grid.pdf
        ├── poisson_2d.pdf
        ├── laplacian_sparsity.pdf
        ├── helmholtz.pdf
        └── bratu.pdf
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

## Chapter 4: The Geometry of Nodes

The codes in `ch04_geometry_of_nodes/` explore polynomial interpolation, the Runge phenomenon, and the superiority of Chebyshev nodes:

| Script | Topic | Description |
|--------|-------|-------------|
| `runge_phenomenon` | Runge phenomenon | Demonstrates divergence of equispaced polynomial interpolation |
| `chebyshev_success` | Chebyshev interpolation | Shows successful Chebyshev interpolation of the Runge function |
| `chebyshev_points_circle` | Geometric construction | Visualizes Chebyshev points as projections from the unit circle |
| `equipotential_curves` | Potential theory | Compares equipotential curves for uniform and Chebyshev distributions |
| `lagrange_basis` | Basis functions | Compares Lagrange basis functions for equispaced vs Chebyshev nodes |
| `lebesgue_functions` | Lebesgue constants | Visualizes Lebesgue functions and their asymptotic growth rates |
| `convergence_comparison` | Convergence rates | Compares convergence of equispaced vs Chebyshev interpolation |

## Chapter 5: Differentiation Matrices

The codes in `ch05_differentiation_matrices/` implement finite difference and spectral differentiation matrices, demonstrating Fornberg's algorithm and comparing convergence rates:

| Script | Topic | Description |
|--------|-------|-------------|
| `fdweights` | Fornberg's algorithm | Computes FD weights for arbitrary node distributions and derivative orders |
| `fd_matrix_periodic` | Periodic FD matrix | Constructs periodic finite difference matrices of various orders (MATLAB only) |
| `spectral_matrix_periodic` | Spectral matrix | Constructs the periodic spectral differentiation matrix using cotangent formula |
| `fd_matrix_bandwidth` | Sparsity patterns | Visualizes bandwidth growth from 2nd order FD to spectral (dense) |
| `spectral_matrix_structure` | Matrix structure | Shows Toeplitz structure and skew-symmetry of the spectral matrix |
| `stencil_pyramid` | Fornberg recursion | Illustrates the recursive structure of Fornberg's weight computation |
| `periodic_cardinal_functions` | Cardinal functions | Visualizes periodic cardinal functions (discrete Dirichlet kernel) |
| `spectral_derivatives_demo` | Derivatives demo | Demonstrates spectral accuracy for first and second derivatives |
| `higher_order_derivatives` | Higher derivatives | Higher-order derivatives (up to 4th) with convergence analysis |
| `convergence_comparison` | Convergence rates | Compares algebraic (FD) vs exponential (spectral) convergence |

## Chapter 6: Chebyshev Differentiation Matrices

The codes in `ch06_chebyshev_differentiation/` implement spectral differentiation on bounded (non-periodic) domains using Chebyshev-Gauss-Lobatto points:

| Script | Topic | Description |
|--------|-------|-------------|
| `cheb_matrix` | Core algorithm | Chebyshev differentiation matrix construction with negative sum trick |
| `cheb_diff_demo` | Differentiation demo | Demonstrates spectral accuracy for derivatives on bounded domains |
| `cheb_convergence` | Convergence analysis | Waterfall plot showing exponential convergence vs polynomial degree |
| `cheb_grid_comparison` | Grid comparison | Compares periodic (equispaced) vs non-periodic (Chebyshev) grids |
| `cheb_matrix_structure` | Matrix structure | Visualizes the dense structure and symmetry of the Chebyshev matrix |
| `cheb_cardinal` | Cardinal functions | Chebyshev cardinal function visualization (Lagrange interpolants) |

## Chapter 7: Boundary Value Problems

The codes in `ch07_boundary_value_problems/` solve various BVPs using Chebyshev spectral collocation:

| Script | Problem | Description |
|--------|---------|-------------|
| `bvp_linear` | 1D Poisson | Solves u_xx = f(x) with Dirichlet boundary conditions |
| `bvp_variable_coeff` | Variable coefficients | BVP with spatially varying coefficients |
| `bvp_eigenvalue` | Eigenvalue problem | Quantum harmonic oscillator eigenvalues via spectral collocation |
| `bvp_2d_poisson` | 2D Poisson | Tensor product grid solution of the 2D Poisson equation |
| `bvp_helmholtz` | Helmholtz equation | u_xx + k²u = f with applications to wave scattering |
| `bvp_nonlinear` | Nonlinear BVP | Bratu equation solved with Newton iteration |

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
