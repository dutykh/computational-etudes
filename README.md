# Computational Études: A Spectral Approach

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
[![Built with Typst](https://img.shields.io/badge/Built%20with-Typst-239DAD.svg)](https://typst.app/)

![Computational Études: A Spectral Approach](assets/Illustration1_compressed.png)

> A pedagogical textbook on spectral methods for differential equations, featuring triple-language implementations in Python, MATLAB, and Julia.

**Author:** Dr. Denys Dutykh
**Affiliation:** Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE

---

## About This Book

Spectral methods are a powerful class of numerical techniques for solving differential equations. Unlike finite difference or finite element methods, which use local approximations, spectral methods represent solutions as linear combinations of global basis functions—typically Fourier series or Chebyshev polynomials. For smooth problems, this approach yields *spectral accuracy*: errors that decay exponentially fast with the number of degrees of freedom.

This book takes a hands-on, pedagogical approach inspired by musical *études*—short compositions designed to develop specific technical skills while remaining artistically complete. Each chapter focuses on a single mathematical concept and explores it through compact, runnable code implementations.

### Key Features

- **Étude-based pedagogy** — Each chapter is a self-contained study combining theory with implementation
- **Triple-language implementations** — All examples provided in Python, MATLAB, and Julia
- **Fully reproducible** — Every figure and result generated from the accompanying code
- **Focus on 1D problems** — Keeps code readable while covering all essential concepts
- **Beautiful typography** — Professionally typeset using Typst with bibliography backreferences

---

## What You Will Learn

By working through this book, you will acquire the following skills:

- **Spectral interpolation and approximation** — Understand why Chebyshev and Fourier bases yield exponential convergence for smooth functions, and how the choice of collocation points determines accuracy.
- **Differentiation matrices** — Build and manipulate the discrete operators that turn differential equations into linear algebra problems.
- **Boundary value problems** — Solve second-order ODEs and PDEs with Dirichlet, Neumann, Robin, and nonlinear boundary conditions using spectral collocation.
- **Time-dependent PDEs** — Combine spectral spatial discretization with method-of-lines time stepping (explicit, implicit, and IMEX schemes) to simulate wave propagation, diffusion, and reaction-diffusion dynamics.
- **Fourier pseudospectral methods** — Exploit the FFT for periodic problems, including dealiasing strategies and conservation properties.
- **Multi-dimensional problems** — Extend spectral methods to 2D via Kronecker products (Cartesian grids) and Chebyshev--Fourier decompositions (polar coordinates).
- **Advanced boundary treatment** — Master both the lifting technique (Method I) and the tau/row-replacement approach (Method II), including frequency-dependent boundary conditions that lead to quadratic eigenvalue problems.
- **Eigenvalue problems** — Compute spectra of differential operators, from classical Sturm--Liouville problems to the quasinormal modes of black holes.

---

## How to Use This Book

This textbook is designed so that *reading and coding go hand in hand*. Here is a recommended workflow:

1. **Read the theory.** Each chapter opens with the mathematical background you need — no more, no less. Definitions and key formulas are highlighted so you can find them later.

2. **Run the code.** Every computational *etude* comes with a self-contained script that produces one or more figures. Run the script, inspect the output, and compare it with the figure in the text. This is the single most effective way to build intuition.

3. **Modify and experiment.** Change a parameter (grid size *N*, time step, physical constant) and observe what happens. Does convergence improve? Does the solution blow up? The études are deliberately short so that experimentation is painless.

4. **Attempt the exercises.** Each chapter ends with exercises that extend the études in new directions. Some ask you to apply the same technique to a different equation; others push you toward open-ended exploration.

5. **Follow the dependency graph.** The chapters are ordered so that each one builds on the previous. Chapters 1--4 lay the mathematical and computational foundations. Chapters 5--8 develop the core spectral machinery on bounded domains. Chapters 9--11 cover Fourier methods for periodic problems. Chapters 12--13 treat special geometries and advanced boundary conditions. You may skip ahead if you already have the prerequisites, but the intended path is sequential.

### The Three-Language Approach

Every étude in this book is implemented in three programming languages: **Python**, **MATLAB**, and **Julia**. The implementations are kept deliberately parallel — the same variable names, the same algorithmic structure, the same numerical parameters — so that switching between languages is effortless.

Why three languages? Because different communities have different preferences, and we want every reader to feel at home:

- **Python** (with NumPy, SciPy, and Matplotlib) is the *lingua franca* of scientific computing today and the natural choice for readers coming from machine learning or data science.
- **MATLAB** remains the workhorse of applied mathematics and engineering. Its matrix-first syntax makes spectral methods particularly natural to express, and the author's personal preference leans toward MATLAB for day-to-day numerical exploration.
- **Julia** combines the readability of MATLAB with the performance of compiled languages. It is a compelling choice for readers who want to scale up from the didactic examples in this book to production-grade simulations.

Pick the language you are most comfortable with and use the other two as references. If you have never written a line of code before, start with Python or MATLAB — the scripts in this book are short enough that you will pick up the syntax as you go.

---

## Table of Contents

- **Preface** — Purpose, audience, and how to use this book
- **Acknowledgements** — Thanks to contributing students
1. **Introduction** — The spectral promise, philosophy of études, collocation methods, and modern workflows
2. **Classical Second Order PDEs and Separation of Variables** — Heat, wave, and Laplace equations; separation of variables as the foundation for spectral methods
3. **Mise en Bouche** — A first taste of spectral methods: method of weighted residuals, collocation vs. Galerkin with low-dimensional examples
4. **The Geometry of Nodes** — Runge phenomenon, potential theory, Chebyshev points, Lebesgue constants, and barycentric interpolation
5. **Differentiation Matrices** — Finite differences as sparse approximations, periodic spectral matrices, Fornberg's algorithm, and spectral convergence
6. **Smoothness and Spectral Accuracy** — Fourier coefficient decay, aliasing, convergence theorems, and the quantum harmonic oscillator
7. **Chebyshev Differentiation Matrices** — Non-periodic spectral methods on bounded domains, Chebyshev-Gauss-Lobatto points, explicit matrix formulas
8. **Boundary Value Problems** — Spectral collocation for BVPs, matrix surgery for boundary conditions, eigenvalue problems, 2D Poisson equation, nonlinear problems
9. **Physical and Fourier Space on Grids** — Continuous, semidiscrete, and discrete Fourier transforms; aliasing and the Nyquist interval; band-limited interpolation; FFT spectra and smoothness
10. **Spectral PDE Solvers with Chebyshev and Fourier Grids** — Chebyshev differentiation via FFT; method of lines; 1D and 2D wave equations; 1D and 2D heat equations; Poisson and Helmholtz equations; variable coefficient transport
11. **Fourier Pseudospectral Methods for Periodic PDEs** — Linear advection, Burgers equation with dealiasing, KdV solitons and Zabusky--Kruskal recurrence, Allen--Cahn phase separation, nonlinear Schrödinger recurrence, Kuramoto--Sivashinsky chaos, 2D Navier--Stokes vorticity dynamics
12. **Spectral Methods in Polar Coordinates** — Coordinate singularity at the origin, Fornberg's doubling trick, Chebyshev--Fourier discretisation of the disk Laplacian, Bessel function eigenmodes, Poisson equation on the disk, heat equation with Crank--Nicolson
13. **Advanced Boundary Conditions for Spectral Methods** — Method I (lifting) and Method II (tau/row replacement), inhomogeneous Dirichlet, Neumann and Robin conditions, Allen--Cahn metastability, nonlinear radiation BCs with Newton iteration, 2D Kronecker product problems, quadratic eigenvalue problems and quasinormal modes of black holes
14. **Higher-Order Boundary Value Problems** — Polynomial trick for clamped BCs, clamped beam under exponential load, beam eigenmodes, coupled second-order systems, 2D biharmonic operator via Kronecker products, clamped plate eigenmodes, quarter-plate symmetry reduction, Orr--Sommerfeld equation and hydrodynamic stability, pseudospectra and non-normality, Kuramoto--Sivashinsky equation with ETDRK4
15. **Quadrature in Spectral Methods: When Exactness Misleads** — Newton--Cotes failure and the Runge function, Clenshaw--Curtis via FFT, Golub--Welsch algorithm for Gauss--Legendre, the six-function convergence race, aliasing in Chebyshev space, complex-plane error portraits, Gauss--Hermite paradox on unbounded domains, approximation spaces and convergence rates, periodic trapezoidal rule and trigonometric exactness
16. **Integration of Periodic Functions: Why the Trapezoidal Rule Becomes Spectral** — Poisson's 1820s ellipse paradox, trigonometric exactness via Fourier series and aliasing, the five-class taxonomy of convergence rates (band-limited, algebraic, geometric, supergeometric, subgeometric), strip-analyticity theorem of Trefethen--Weideman, the doubled-rate observation, real-line trapezoidal rule on the Gaussian, FFT computation of Fourier coefficients
17. **Time Stepping, Stability, and the CFL Constraint** — The catastrophe gallery (deliberate blow-ups), stability regions of classical integrators, Fourier CFL scaling (constant and variable coefficients), Chebyshev stiffness and outlier eigenvalues, pseudospectra for nonnormal operators, fair comparison of six time-stepping cures (FE, BE, CN, IF-RK4, Dufort--Frankel, RKC), KdV integrating-factor RK4

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

**For running Julia code:**
- Julia 1.10+
- CairoMakie.jl (plotting)
- FFTW.jl (Fast Fourier Transform)

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

# Generate Python figures
make figures-python

# Generate Julia figures
make figures-julia

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
python codes/python/ch02/heat_equation_evolution.py
python codes/python/ch02/heat_equation_waterfall.py
python codes/python/ch02/wave_equation_evolution.py
python codes/python/ch02/wave_equation_waterfall.py
python codes/python/ch02/laplace_equation_2d.py

# Chapter 3: Mise en Bouche
python codes/python/ch03/collocation_example1.py
python codes/python/ch03/collocation_vs_galerkin.py

# Chapter 4: The Geometry of Nodes
python codes/python/ch04/runge_phenomenon.py
python codes/python/ch04/chebyshev_success.py
python codes/python/ch04/chebyshev_points_circle.py
python codes/python/ch04/equipotential_curves.py
python codes/python/ch04/lagrange_basis.py
python codes/python/ch04/lebesgue_functions.py
python codes/python/ch04/convergence_comparison.py
python codes/python/ch04/convergence_zoom.py
python codes/python/ch04/lebesgue_constants_zoom.py
python codes/python/ch04/lebesgue_random_nodes.py
python codes/python/ch04/lebesgue_random_chebyshev.py

# Chapter 5: Differentiation Matrices
python codes/python/ch05/fd_matrix_bandwidth.py
python codes/python/ch05/spectral_matrix_structure.py
python codes/python/ch05/stencil_pyramid.py
python codes/python/ch05/periodic_cardinal_functions.py
python codes/python/ch05/spectral_derivatives_demo.py
python codes/python/ch05/higher_order_derivatives.py
python codes/python/ch05/convergence_comparison.py

# Chapter 6: Smoothness and Spectral Accuracy
python codes/python/ch06/fourier_decay.py
python codes/python/ch06/aliasing_demo.py
python codes/python/ch06/convergence_rates.py
python codes/python/ch06/harmonic_oscillator.py

# Chapter 7: Chebyshev Differentiation Matrices
python codes/python/ch07/cheb_matrix.py
python codes/python/ch07/cheb_diff_demo.py
python codes/python/ch07/cheb_convergence.py
python codes/python/ch07/cheb_grid_comparison.py
python codes/python/ch07/cheb_matrix_structure.py
python codes/python/ch07/cheb_cardinal.py

# Chapter 8: Boundary Value Problems
python codes/python/ch08/bvp_linear.py
python codes/python/ch08/bvp_variable_coeff.py
python codes/python/ch08/bvp_eigenvalue.py
python codes/python/ch08/bvp_2d_poisson.py
python codes/python/ch08/bvp_helmholtz.py
python codes/python/ch08/bvp_nonlinear.py

# Chapter 9: Physical and Fourier Space on Grids
python codes/python/ch09/two_views_function.py
python codes/python/ch09/aliasing_demo.py
python codes/python/ch09/sinc_interpolation.py
python codes/python/ch09/fft_aliasing.py
python codes/python/ch09/smoothness_spectra.py
python codes/python/ch09/zero_padding_interpolation.py

# Chapter 10: Spectral PDE Solvers
python codes/python/ch10/chebfft.py
python codes/python/ch10/cheb_fourier_geometry.py
python codes/python/ch10/chebfft_accuracy.py
python codes/python/ch10/wave1d_cheb.py
python codes/python/ch10/wave2d_cheb.py
python codes/python/ch10/heat1d_cheb.py
python codes/python/ch10/heat2d_cheb.py
python codes/python/ch10/poisson2d_cheb.py
python codes/python/ch10/transport_variable.py

# Chapter 11: Fourier Pseudospectral Methods
python codes/python/ch11/fourier_diff_comparison.py
python codes/python/ch11/advection_fourier.py
python codes/python/ch11/burgers_fourier.py
python codes/python/ch11/kdv_soliton.py
python codes/python/ch11/kdv_zabusky_kruskal.py
python codes/python/ch11/allen_cahn.py
python codes/python/ch11/schrodinger.py
python codes/python/ch11/nls_recurrence.py
python codes/python/ch11/kuramoto_sivashinsky.py
python codes/python/ch11/navier_stokes_2d.py
python codes/python/ch11/transport_variable.py

# Chapter 12: Spectral Methods in Polar Coordinates
python codes/python/ch12/polar_grid_geometry.py
python codes/python/ch12/disk_eigenmodes.py
python codes/python/ch12/disk_nodal_rotation.py
python codes/python/ch12/disk_poisson.py
python codes/python/ch12/disk_heat.py

# Chapter 13: Advanced Boundary Conditions
python codes/python/ch13/bc_inhom_lifting.py
python codes/python/ch13/bc_helmholtz_robin.py
python codes/python/ch13/bc_allen_cahn.py
python codes/python/ch13/bc_radiative.py
python codes/python/ch13/bc_laplace_2d.py
python codes/python/ch13/bc_qnm_poschl_teller.py
python codes/python/ch13/bc_vibrating_string.py

# Chapter 14: Higher-Order Boundary Value Problems
python codes/python/ch14/ho_clamped_beam.py
python codes/python/ch14/ho_beam_eigenmodes.py
python codes/python/ch14/ho_coupled_comparison.py
python codes/python/ch14/ho_plate_eigenmodes.py
python codes/python/ch14/ho_quarter_plate.py
python codes/python/ch14/ho_orr_sommerfeld.py
python codes/python/ch14/ho_pseudospectra.py
python codes/python/ch14/ho_kuramoto_sivashinsky.py

# Chapter 15: Quadrature in Spectral Methods
python codes/python/ch15/quad_node_visualization.py
python codes/python/ch15/quad_polynomial_exactness.py
python codes/python/ch15/quad_newton_cotes_runge.py
python codes/python/ch15/quad_exactness_table.py
python codes/python/ch15/quad_gauss_cc_construction.py
python codes/python/ch15/quad_convergence_race.py
python codes/python/ch15/quad_aliasing_chebyshev.py
python codes/python/ch15/quad_complex_plane.py
python codes/python/ch15/quad_gauss_hermite_weights.py
python codes/python/ch15/quad_gauss_hermite_failure.py
python codes/python/ch15/quad_convergence_rates.py

# Chapter 16: Integration of Periodic Functions
python codes/python/ch16/trap_poisson_ellipse.py
python codes/python/ch16/trap_band_limited.py
python codes/python/ch16/trap_algebraic_decay.py
python codes/python/ch16/trap_poisson_kernel.py
python codes/python/ch16/trap_supergeometric.py
python codes/python/ch16/trap_subgeometric.py
python codes/python/ch16/trap_real_line.py
python codes/python/ch16/trap_fft_coefficients.py

# Chapter 17: Time Stepping, Stability, and the CFL Constraint
python codes/python/ch17/catastrophe_gallery.py
python codes/python/ch17/stability_regions.py
python codes/python/ch17/fourier_cfl.py
python codes/python/ch17/fourier_cfl_variable.py
python codes/python/ch17/chebyshev_stiffness.py
python codes/python/ch17/pseudospectra_demo.py
python codes/python/ch17/fair_comparison.py
python codes/python/ch17/kdv_rk_comparison.py
```

**Julia:**
```bash
# Chapter 2: Classical PDEs
julia codes/julia/ch02/heat_equation_evolution.jl
julia codes/julia/ch02/heat_equation_waterfall.jl
julia codes/julia/ch02/wave_equation_evolution.jl
julia codes/julia/ch02/wave_equation_waterfall.jl
julia codes/julia/ch02/laplace_equation_2d.jl

# Chapter 3: Mise en Bouche
julia codes/julia/ch03/collocation_example1.jl
julia codes/julia/ch03/collocation_vs_galerkin.jl

# Chapter 4: The Geometry of Nodes
julia codes/julia/ch04/runge_phenomenon.jl
julia codes/julia/ch04/chebyshev_success.jl
julia codes/julia/ch04/chebyshev_points_circle.jl
julia codes/julia/ch04/equipotential_curves.jl
julia codes/julia/ch04/lagrange_basis.jl
julia codes/julia/ch04/lebesgue_functions.jl
julia codes/julia/ch04/convergence_comparison.jl
julia codes/julia/ch04/convergence_zoom.jl
julia codes/julia/ch04/lebesgue_constants_zoom.jl

# Chapter 5: Differentiation Matrices
julia codes/julia/ch05/fd_matrix_bandwidth.jl
julia codes/julia/ch05/spectral_matrix_structure.jl
julia codes/julia/ch05/stencil_pyramid.jl
julia codes/julia/ch05/spectral_derivatives_demo.jl
julia codes/julia/ch05/convergence_comparison.jl

# Chapter 6: Smoothness and Spectral Accuracy
julia codes/julia/ch06/fourier_decay.jl
julia codes/julia/ch06/aliasing_demo.jl
julia codes/julia/ch06/convergence_rates.jl

# Chapter 7: Chebyshev Differentiation Matrices
julia codes/julia/ch07/cheb_matrix.jl
julia codes/julia/ch07/cheb_grid_comparison.jl
julia codes/julia/ch07/cheb_matrix_structure.jl
julia codes/julia/ch07/cheb_cardinal.jl
julia codes/julia/ch07/cheb_diff_demo.jl
julia codes/julia/ch07/cheb_convergence.jl

# Chapter 8: Boundary Value Problems
julia codes/julia/ch08/bvp_linear.jl
julia codes/julia/ch08/bvp_variable_coeff.jl
julia codes/julia/ch08/bvp_nonlinear.jl
julia codes/julia/ch08/bvp_eigenvalue.jl
julia codes/julia/ch08/bvp_2d_poisson.jl
julia codes/julia/ch08/bvp_helmholtz.jl
julia codes/julia/ch08/harmonic_oscillator.jl

# Chapter 9: Physical and Fourier Space on Grids
julia codes/julia/ch09/two_views_function.jl
julia codes/julia/ch09/aliasing_demo.jl
julia codes/julia/ch09/sinc_interpolation.jl
julia codes/julia/ch09/fft_aliasing.jl
julia codes/julia/ch09/smoothness_spectra.jl
julia codes/julia/ch09/zero_padding_interpolation.jl

# Chapter 10: Spectral PDE Solvers
julia codes/julia/ch10/cheb_fourier_geometry.jl
julia codes/julia/ch10/chebfft_accuracy.jl
julia codes/julia/ch10/wave1d_cheb.jl
julia codes/julia/ch10/wave2d_cheb.jl
julia codes/julia/ch10/heat1d_cheb.jl
julia codes/julia/ch10/heat2d_cheb.jl
julia codes/julia/ch10/poisson2d_cheb.jl
julia codes/julia/ch10/transport_variable.jl

# Chapter 11: Fourier Pseudospectral Methods
julia codes/julia/ch11/fourier_diff_comparison.jl
julia codes/julia/ch11/advection_fourier.jl
julia codes/julia/ch11/burgers_fourier.jl
julia codes/julia/ch11/kdv_soliton.jl
julia codes/julia/ch11/kdv_zabusky_kruskal.jl
julia codes/julia/ch11/allen_cahn.jl
julia codes/julia/ch11/schrodinger.jl
julia codes/julia/ch11/nls_recurrence.jl
julia codes/julia/ch11/kuramoto_sivashinsky.jl
julia codes/julia/ch11/navier_stokes_2d.jl
julia codes/julia/ch11/transport_variable.jl

# Chapter 12: Spectral Methods in Polar Coordinates
julia codes/julia/ch12/polar_grid_geometry.jl
julia codes/julia/ch12/disk_eigenmodes.jl
julia codes/julia/ch12/disk_nodal_rotation.jl
julia codes/julia/ch12/disk_poisson.jl
julia codes/julia/ch12/disk_heat.jl

# Chapter 13: Advanced Boundary Conditions
julia codes/julia/ch13/bc_inhom_lifting.jl
julia codes/julia/ch13/bc_helmholtz_robin.jl
julia codes/julia/ch13/bc_allen_cahn.jl
julia codes/julia/ch13/bc_radiative.jl
julia codes/julia/ch13/bc_laplace_2d.jl
julia codes/julia/ch13/bc_qnm_poschl_teller.jl
julia codes/julia/ch13/bc_vibrating_string.jl

# Chapter 14: Higher-Order Boundary Value Problems
julia codes/julia/ch14/ho_clamped_beam.jl
julia codes/julia/ch14/ho_beam_eigenmodes.jl
julia codes/julia/ch14/ho_coupled_comparison.jl
julia codes/julia/ch14/ho_plate_eigenmodes.jl
julia codes/julia/ch14/ho_quarter_plate.jl
julia codes/julia/ch14/ho_orr_sommerfeld.jl
julia codes/julia/ch14/ho_pseudospectra.jl
julia codes/julia/ch14/ho_kuramoto_sivashinsky.jl

# Chapter 15: Quadrature in Spectral Methods
julia codes/julia/ch15/quad_node_visualization.jl
julia codes/julia/ch15/quad_polynomial_exactness.jl
julia codes/julia/ch15/quad_newton_cotes_runge.jl
julia codes/julia/ch15/quad_exactness_table.jl
julia codes/julia/ch15/quad_gauss_cc_construction.jl
julia codes/julia/ch15/quad_convergence_race.jl
julia codes/julia/ch15/quad_aliasing_chebyshev.jl
julia codes/julia/ch15/quad_complex_plane.jl
julia codes/julia/ch15/quad_gauss_hermite_weights.jl
julia codes/julia/ch15/quad_gauss_hermite_failure.jl
julia codes/julia/ch15/quad_convergence_rates.jl

# Chapter 16: Integration of Periodic Functions
julia codes/julia/ch16/trap_poisson_ellipse.jl
julia codes/julia/ch16/trap_band_limited.jl
julia codes/julia/ch16/trap_algebraic_decay.jl
julia codes/julia/ch16/trap_poisson_kernel.jl
julia codes/julia/ch16/trap_supergeometric.jl
julia codes/julia/ch16/trap_subgeometric.jl
julia codes/julia/ch16/trap_real_line.jl
julia codes/julia/ch16/trap_fft_coefficients.jl

# Chapter 17: Time Stepping, Stability, and the CFL Constraint
julia codes/julia/ch17/catastrophe_gallery.jl
julia codes/julia/ch17/stability_regions.jl
julia codes/julia/ch17/fourier_cfl.jl
julia codes/julia/ch17/fourier_cfl_variable.jl
julia codes/julia/ch17/chebyshev_stiffness.jl
julia codes/julia/ch17/pseudospectra_demo.jl
julia codes/julia/ch17/fair_comparison.jl
julia codes/julia/ch17/kdv_rk_comparison.jl
```

**MATLAB:**
```matlab
% Navigate to the codes directory
cd codes/matlab/ch02
heat_equation_evolution
heat_equation_waterfall
wave_equation_evolution
wave_equation_waterfall
laplace_equation_2d

cd ../ch03
collocation_example1
collocation_vs_galerkin

cd ../ch04
runge_phenomenon
chebyshev_success
chebyshev_points_circle
equipotential_curves
lagrange_basis
lebesgue_functions
convergence_comparison
convergence_zoom
lebesgue_constants_zoom

cd ../ch05
fd_matrix_bandwidth
spectral_matrix_structure
stencil_pyramid
periodic_cardinal_functions
spectral_derivatives_demo
higher_order_derivatives
convergence_comparison

cd ../ch06
fourier_decay
aliasing_demo
convergence_rates
harmonic_oscillator

cd ../ch07
cheb_matrix
verify_cheb_matrix
cheb_grid_comparison
cheb_matrix_structure
cheb_cardinal
cheb_diff_demo
cheb_convergence

cd ../ch08
bvp_linear
bvp_variable_coeff
bvp_eigenvalue
bvp_2d_poisson
bvp_helmholtz
bvp_nonlinear

cd ../ch09
two_views_function
aliasing_demo
sinc_interpolation
fft_aliasing
smoothness_spectra
zero_padding_interpolation

cd ../ch10
cheb_fourier_geometry
chebfft_accuracy
wave1d_cheb
wave2d_cheb
heat1d_cheb
heat2d_cheb
poisson2d_cheb
transport_variable

cd ../ch11
fourier_diff_comparison
advection_fourier
burgers_fourier
kdv_soliton
kdv_zabusky_kruskal
allen_cahn
schrodinger
nls_recurrence
kuramoto_sivashinsky
navier_stokes_2d
transport_variable

cd ../ch12
polar_grid_geometry
disk_eigenmodes
disk_nodal_rotation
disk_poisson
disk_heat

cd ../ch13
bc_inhom_lifting
bc_helmholtz_robin
bc_allen_cahn
bc_radiative
bc_laplace_2d
bc_qnm_poschl_teller
bc_vibrating_string

cd ../ch14
ho_clamped_beam
ho_beam_eigenmodes
ho_coupled_comparison
ho_plate_eigenmodes
ho_quarter_plate
ho_orr_sommerfeld
ho_pseudospectra
ho_kuramoto_sivashinsky

cd ../ch15
quad_node_visualization
quad_polynomial_exactness
quad_newton_cotes_runge
quad_exactness_table
quad_gauss_cc_construction
quad_convergence_race
quad_aliasing_chebyshev
quad_complex_plane
quad_gauss_hermite_weights
quad_gauss_hermite_failure
quad_convergence_rates

cd ../ch16
trap_poisson_ellipse
trap_band_limited
trap_algebraic_decay
trap_poisson_kernel
trap_supergeometric
trap_subgeometric
trap_real_line
trap_fft_coefficients

cd ../ch17
catastrophe_gallery
stability_regions
fourier_cfl
fourier_cfl_variable
chebyshev_stiffness
pseudospectra_demo
fair_comparison
kdv_rk_comparison
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
│   │   ├── geometry_of_nodes.typ
│   │   ├── differentiation_matrices.typ
│   │   ├── smoothness_accuracy.typ
│   │   ├── chebyshev_differentiation.typ
│   │   ├── boundary_value_problems.typ
│   │   ├── fourier_grids.typ
│   │   ├── spectral_pde_solvers.typ
│   │   ├── fourier_pseudospectral.typ
│   │   ├── polar_coordinates.typ
│   │   └── advanced_boundary_conditions.typ
│   ├── styles/                  # Typography and layout
│   │   └── template.typ
│   ├── biblio/                  # Bibliography
│   │   └── library.bib
│   ├── figures/                 # Generated figures
│   │   ├── ch02/
│   │   │   ├── python/
│   │   │   ├── matlab/
│   │   │   └── julia/
│   │   ├── ch03/
│   │   │   ├── python/
│   │   │   ├── matlab/
│   │   │   └── julia/
│   │   ├── ch04/
│   │   │   ├── python/
│   │   │   ├── matlab/
│   │   │   └── julia/
│   │   ├── ch05/
│   │   │   ├── python/
│   │   │   ├── matlab/
│   │   │   └── julia/
│   │   ├── ch06/
│   │   │   ├── python/
│   │   │   ├── matlab/
│   │   │   └── julia/
│   │   ├── ch07/
│   │   │   ├── python/
│   │   │   ├── matlab/
│   │   │   └── julia/
│   │   ├── ch08/
│   │   │   ├── python/
│   │   │   ├── matlab/
│   │   │   └── julia/
│   │   ├── ch09/
│   │   │   ├── python/
│   │   │   ├── matlab/
│   │   │   └── julia/
│   │   ├── ch10/
│   │   │   ├── python/
│   │   │   ├── matlab/
│   │   │   └── julia/
│   │   ├── ch11/
│   │   │   ├── python/
│   │   │   ├── matlab/
│   │   │   └── julia/
│   │   ├── ch12/
│   │   │   ├── python/
│   │   │   ├── matlab/
│   │   │   └── julia/
│   │   ├── ch13/
│   │   │   ├── python/
│   │   │   ├── matlab/
│   │   │   └── julia/
│   │   ├── ch14/
│   │   │   ├── python/
│   │   │   ├── matlab/
│   │   │   └── julia/
│   │   ├── ch15/
│   │   │   ├── python/
│   │   │   ├── matlab/
│   │   │   └── julia/
│   │   ├── ch16/
│   │   │   ├── python/
│   │   │   ├── matlab/
│   │   │   └── julia/
│   │   └── ch17/
│   │       ├── python/
│   │       ├── matlab/
│   │       └── julia/
│   └── build/                   # Compiled PDF output
├── codes/                       # Code implementations
│   ├── python/
│   │   ├── ch02/                # Classical PDEs
│   │   ├── ch03/                # Mise en Bouche
│   │   ├── ch04/                # Geometry of Nodes
│   │   ├── ch05/                # Differentiation Matrices
│   │   ├── ch06/                # Smoothness and Spectral Accuracy
│   │   ├── ch07/                # Chebyshev Differentiation
│   │   ├── ch08/                # Boundary Value Problems
│   │   ├── ch09/                # Physical and Fourier Space on Grids
│   │   ├── ch10/                # Spectral PDE Solvers
│   │   ├── ch11/                # Fourier Pseudospectral Methods
│   │   ├── ch12/                # Spectral Methods in Polar Coordinates
│   │   ├── ch13/                # Advanced Boundary Conditions
│   │   ├── ch14/                # Higher-Order Boundary Value Problems
│   │   ├── ch15/                # Quadrature in Spectral Methods
│   │   ├── ch16/                # Integration of Periodic Functions
│   │   └── ch17/                # Time Stepping and CFL Constraint
│   ├── matlab/
│   │   ├── ch02/                # Classical PDEs
│   │   ├── ch03/                # Mise en Bouche
│   │   ├── ch04/                # Geometry of Nodes
│   │   ├── ch05/                # Differentiation Matrices
│   │   ├── ch06/                # Smoothness and Spectral Accuracy
│   │   ├── ch07/                # Chebyshev Differentiation
│   │   ├── ch08/                # Boundary Value Problems
│   │   ├── ch09/                # Physical and Fourier Space on Grids
│   │   ├── ch10/                # Spectral PDE Solvers
│   │   ├── ch11/                # Fourier Pseudospectral Methods
│   │   ├── ch12/                # Spectral Methods in Polar Coordinates
│   │   ├── ch13/                # Advanced Boundary Conditions
│   │   ├── ch14/                # Higher-Order Boundary Value Problems
│   │   ├── ch15/                # Quadrature in Spectral Methods
│   │   ├── ch16/                # Integration of Periodic Functions
│   │   └── ch17/                # Time Stepping and CFL Constraint
│   ├── julia/
│   │   ├── ch02/                # Classical PDEs
│   │   ├── ch03/                # Mise en Bouche
│   │   ├── ch04/                # Geometry of Nodes
│   │   ├── ch05/                # Differentiation Matrices
│   │   ├── ch06/                # Smoothness and Spectral Accuracy
│   │   ├── ch07/                # Chebyshev Differentiation
│   │   ├── ch08/                # Boundary Value Problems
│   │   ├── ch09/                # Physical and Fourier Space on Grids
│   │   ├── ch10/                # Spectral PDE Solvers
│   │   ├── ch11/                # Fourier Pseudospectral Methods
│   │   ├── ch12/                # Spectral Methods in Polar Coordinates
│   │   ├── ch13/                # Advanced Boundary Conditions
│   │   ├── ch14/                # Higher-Order Boundary Value Problems
│   │   ├── ch15/                # Quadrature in Spectral Methods
│   │   ├── ch16/                # Integration of Periodic Functions
│   │   └── ch17/                # Time Stepping and CFL Constraint
│   └── README.md
├── slides/                      # Lecture slides
│   └── QNM-SpectralConvergence.pdf  # Spectral method convergence proof for QNMs
├── tplan/                       # Teaching plan (MATH 794)
│   ├── teaching_plan.typ
│   ├── Makefile
│   └── build/
├── Makefile                     # Build automation
├── LICENSE                      # CC BY-NC-SA 4.0
└── README.md
```

---

## Typst Packages Used

The textbook uses the following Typst packages:
- **[codly](https://typst.app/universe/package/codly)** — Beautiful code blocks with syntax highlighting
- **[retrofit](https://typst.app/universe/package/retrofit)** — Bibliography backreferences showing citation locations

---

## Related Resources

For readers who want to quickly experience the power of spectral methods without writing code from scratch, we highly recommend:

- **[Chebfun](https://www.chebfun.org/)** — An open-source MATLAB package that extends MATLAB's capabilities to continuous functions. Chebfun automatically determines the number of Chebyshev points needed to represent a function to machine precision, making spectral methods as easy to use as built-in MATLAB commands. Operations like differentiation, integration, and solving differential equations become one-liners. Chebfun is an excellent tool for exploring the concepts in this book interactively.

---

## Teaching Materials

This book is used for **MATH 794** at Khalifa University. A tentative teaching plan is available in the `tplan/` directory, which includes:
- Weekly lecture schedule
- Chapter-to-lecture mapping
- Progress tracking

Lecture slides are available in the `slides/` directory:
- `QNM-SpectralConvergence.pdf` — Overview of the convergence proof for spectral methods applied to quasinormal modes

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
