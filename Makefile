.PHONY: textbook clean figures figures-python figures-matlab figures-julia all tplan

# Tools
TYPST ?= typst
PYTHON ?= python3
MATLAB ?= matlab
JULIA ?= $(HOME)/.juliaup/bin/julia

# Textbook compilation
SRC = textbook/main.typ
OUT_DIR = textbook/build
OUT = $(OUT_DIR)/DD-Computational-Etudes-2026.pdf

# Teaching plan compilation
TPLAN_SRC = tplan/teaching_plan.typ
TPLAN_OUT_DIR = tplan/build
TPLAN_OUT = $(TPLAN_OUT_DIR)/teaching_plan.pdf

# Python scripts - Chapter 2
PY_CH02 = codes/python/ch02
PY_SCRIPTS_CH02 = $(PY_CH02)/heat_equation_evolution.py \
                  $(PY_CH02)/heat_equation_waterfall.py \
                  $(PY_CH02)/wave_equation_evolution.py \
                  $(PY_CH02)/wave_equation_waterfall.py \
                  $(PY_CH02)/laplace_equation_2d.py

# Python scripts - Chapter 3
PY_CH03 = codes/python/ch03
PY_SCRIPTS_CH03 = $(PY_CH03)/collocation_example1.py \
                  $(PY_CH03)/collocation_vs_galerkin.py

# Python scripts - Chapter 4
PY_CH04 = codes/python/ch04
PY_SCRIPTS_CH04 = $(PY_CH04)/runge_phenomenon.py \
                  $(PY_CH04)/chebyshev_success.py \
                  $(PY_CH04)/chebyshev_points_circle.py \
                  $(PY_CH04)/equipotential_curves.py \
                  $(PY_CH04)/lagrange_basis.py \
                  $(PY_CH04)/lebesgue_functions.py \
                  $(PY_CH04)/lebesgue_constants_zoom.py \
                  $(PY_CH04)/lebesgue_random_nodes.py \
                  $(PY_CH04)/convergence_comparison.py \
                  $(PY_CH04)/convergence_zoom.py

# Python scripts - Chapter 5
PY_CH05 = codes/python/ch05
PY_SCRIPTS_CH05 = $(PY_CH05)/fd_matrix_bandwidth.py \
                  $(PY_CH05)/spectral_matrix_structure.py \
                  $(PY_CH05)/fd_stencil_schematic.py \
                  $(PY_CH05)/stencil_pyramid.py \
                  $(PY_CH05)/convergence_comparison.py \
                  $(PY_CH05)/spectral_derivatives_demo.py

# Python scripts - Chapter 6 (Smoothness and Spectral Accuracy)
PY_CH06 = codes/python/ch06
PY_SCRIPTS_CH06 = $(PY_CH06)/fourier_decay.py \
                  $(PY_CH06)/aliasing_demo.py \
                  $(PY_CH06)/convergence_rates.py

# Python scripts - Chapter 7 (Chebyshev Differentiation)
PY_CH07 = codes/python/ch07
PY_SCRIPTS_CH07 = $(PY_CH07)/cheb_matrix.py \
                  $(PY_CH07)/cheb_grid_comparison.py \
                  $(PY_CH07)/cheb_matrix_structure.py \
                  $(PY_CH07)/cheb_cardinal.py \
                  $(PY_CH07)/cheb_diff_demo.py \
                  $(PY_CH07)/cheb_convergence.py

# Python scripts - Chapter 8 (Boundary Value Problems)
PY_CH08 = codes/python/ch08
PY_SCRIPTS_CH08 = $(PY_CH08)/bvp_linear.py \
                  $(PY_CH08)/bvp_variable_coeff.py \
                  $(PY_CH08)/bvp_nonlinear.py \
                  $(PY_CH08)/bvp_eigenvalue.py \
                  $(PY_CH08)/bvp_2d_poisson.py \
                  $(PY_CH08)/bvp_helmholtz.py \
                  $(PY_CH08)/harmonic_oscillator.py

# Python scripts - Chapter 9 (Fourier Grids)
PY_CH09 = codes/python/ch09
PY_SCRIPTS_CH09 = $(PY_CH09)/two_views_function.py \
                  $(PY_CH09)/aliasing_demo.py \
                  $(PY_CH09)/sinc_interpolation.py \
                  $(PY_CH09)/fft_aliasing.py \
                  $(PY_CH09)/smoothness_spectra.py \
                  $(PY_CH09)/zero_padding_interpolation.py

# Python scripts - Chapter 10 (Spectral PDE Solvers)
PY_CH10 = codes/python/ch10
PY_SCRIPTS_CH10 = $(PY_CH10)/chebfft.py \
                  $(PY_CH10)/cheb_fourier_geometry.py \
                  $(PY_CH10)/chebfft_accuracy.py \
                  $(PY_CH10)/wave1d_cheb.py \
                  $(PY_CH10)/wave2d_cheb.py \
                  $(PY_CH10)/heat1d_cheb.py \
                  $(PY_CH10)/heat2d_cheb.py \
                  $(PY_CH10)/poisson2d_cheb.py \
                  $(PY_CH10)/transport_variable.py

# MATLAB scripts - Chapter 2
M_CH02 = codes/matlab/ch02
M_SCRIPTS_CH02 = $(M_CH02)/heat_equation_evolution.m \
                 $(M_CH02)/heat_equation_waterfall.m \
                 $(M_CH02)/wave_equation_evolution.m \
                 $(M_CH02)/wave_equation_waterfall.m \
                 $(M_CH02)/laplace_equation_2d.m

# MATLAB scripts - Chapter 3
M_CH03 = codes/matlab/ch03
M_SCRIPTS_CH03 = $(M_CH03)/collocation_example1.m \
                 $(M_CH03)/collocation_vs_galerkin.m

# MATLAB scripts - Chapter 4
M_CH04 = codes/matlab/ch04
M_SCRIPTS_CH04 = $(M_CH04)/runge_phenomenon.m \
                 $(M_CH04)/chebyshev_success.m \
                 $(M_CH04)/chebyshev_points_circle.m \
                 $(M_CH04)/equipotential_curves.m \
                 $(M_CH04)/lagrange_basis.m \
                 $(M_CH04)/lebesgue_functions.m \
                 $(M_CH04)/lebesgue_random_nodes.m \
                 $(M_CH04)/convergence_comparison.m

# MATLAB scripts - Chapter 5
M_CH05 = codes/matlab/ch05
M_SCRIPTS_CH05 = $(M_CH05)/fd_matrix_bandwidth.m \
                 $(M_CH05)/spectral_matrix_structure.m \
                 $(M_CH05)/fd_stencil_schematic.m \
                 $(M_CH05)/stencil_pyramid.m \
                 $(M_CH05)/convergence_comparison.m \
                 $(M_CH05)/spectral_derivatives_demo.m

# MATLAB scripts - Chapter 6 (Smoothness and Spectral Accuracy)
M_CH06 = codes/matlab/ch06
M_SCRIPTS_CH06 = $(M_CH06)/fourier_decay.m \
                 $(M_CH06)/aliasing_demo.m \
                 $(M_CH06)/convergence_rates.m

# MATLAB scripts - Chapter 7 (Chebyshev Differentiation)
M_CH07 = codes/matlab/ch07
M_SCRIPTS_CH07 = $(M_CH07)/cheb_matrix.m \
                 $(M_CH07)/cheb_grid_comparison.m \
                 $(M_CH07)/cheb_matrix_structure.m \
                 $(M_CH07)/cheb_cardinal.m \
                 $(M_CH07)/cheb_diff_demo.m \
                 $(M_CH07)/cheb_convergence.m

# MATLAB scripts - Chapter 8 (Boundary Value Problems)
M_CH08 = codes/matlab/ch08
M_SCRIPTS_CH08 = $(M_CH08)/bvp_linear.m \
                 $(M_CH08)/bvp_variable_coeff.m \
                 $(M_CH08)/bvp_nonlinear.m \
                 $(M_CH08)/bvp_eigenvalue.m \
                 $(M_CH08)/bvp_2d_poisson.m \
                 $(M_CH08)/bvp_helmholtz.m \
                 $(M_CH08)/harmonic_oscillator.m

# MATLAB scripts - Chapter 9 (Fourier Grids)
M_CH09 = codes/matlab/ch09
M_SCRIPTS_CH09 = $(M_CH09)/two_views_function.m \
                 $(M_CH09)/aliasing_demo.m \
                 $(M_CH09)/sinc_interpolation.m \
                 $(M_CH09)/fft_aliasing.m \
                 $(M_CH09)/smoothness_spectra.m \
                 $(M_CH09)/zero_padding_interpolation.m

# MATLAB scripts - Chapter 10 (Spectral PDE Solvers)
M_CH10 = codes/matlab/ch10
M_SCRIPTS_CH10 = $(M_CH10)/chebfft.m \
                 $(M_CH10)/cheb_matrix.m \
                 $(M_CH10)/cheb_fourier_geometry.m \
                 $(M_CH10)/chebfft_accuracy.m \
                 $(M_CH10)/wave1d_cheb.m \
                 $(M_CH10)/wave2d_cheb.m \
                 $(M_CH10)/heat1d_cheb.m \
                 $(M_CH10)/heat2d_cheb.m \
                 $(M_CH10)/poisson2d_cheb.m \
                 $(M_CH10)/transport_variable.m

# Julia scripts - Chapter 2
JL_CH02 = codes/julia/ch02
JL_SCRIPTS_CH02 = $(JL_CH02)/heat_equation_evolution.jl \
                  $(JL_CH02)/heat_equation_waterfall.jl \
                  $(JL_CH02)/wave_equation_evolution.jl \
                  $(JL_CH02)/wave_equation_waterfall.jl \
                  $(JL_CH02)/laplace_equation_2d.jl

# Julia scripts - Chapter 3
JL_CH03 = codes/julia/ch03
JL_SCRIPTS_CH03 = $(JL_CH03)/collocation_example1.jl \
                  $(JL_CH03)/collocation_vs_galerkin.jl

# Julia scripts - Chapter 4
JL_CH04 = codes/julia/ch04
JL_SCRIPTS_CH04 = $(JL_CH04)/runge_phenomenon.jl \
                  $(JL_CH04)/chebyshev_success.jl \
                  $(JL_CH04)/chebyshev_points_circle.jl \
                  $(JL_CH04)/equipotential_curves.jl \
                  $(JL_CH04)/lagrange_basis.jl \
                  $(JL_CH04)/lebesgue_functions.jl \
                  $(JL_CH04)/lebesgue_constants_zoom.jl \
                  $(JL_CH04)/lebesgue_random_nodes.jl \
                  $(JL_CH04)/convergence_comparison.jl \
                  $(JL_CH04)/convergence_zoom.jl

# Julia scripts - Chapter 5
JL_CH05 = codes/julia/ch05
JL_SCRIPTS_CH05 = $(JL_CH05)/fd_matrix_bandwidth.jl \
                  $(JL_CH05)/spectral_matrix_structure.jl \
                  $(JL_CH05)/fd_stencil_schematic.jl \
                  $(JL_CH05)/stencil_pyramid.jl \
                  $(JL_CH05)/convergence_comparison.jl \
                  $(JL_CH05)/spectral_derivatives_demo.jl

# Julia scripts - Chapter 6 (Smoothness and Spectral Accuracy)
JL_CH06 = codes/julia/ch06
JL_SCRIPTS_CH06 = $(JL_CH06)/fourier_decay.jl \
                  $(JL_CH06)/aliasing_demo.jl \
                  $(JL_CH06)/convergence_rates.jl

# Julia scripts - Chapter 7 (Chebyshev Differentiation)
JL_CH07 = codes/julia/ch07
JL_SCRIPTS_CH07 = $(JL_CH07)/cheb_matrix.jl \
                  $(JL_CH07)/cheb_grid_comparison.jl \
                  $(JL_CH07)/cheb_matrix_structure.jl \
                  $(JL_CH07)/cheb_cardinal.jl \
                  $(JL_CH07)/cheb_diff_demo.jl \
                  $(JL_CH07)/cheb_convergence.jl

# Julia scripts - Chapter 8 (Boundary Value Problems)
JL_CH08 = codes/julia/ch08
JL_SCRIPTS_CH08 = $(JL_CH08)/bvp_linear.jl \
                  $(JL_CH08)/bvp_variable_coeff.jl \
                  $(JL_CH08)/bvp_nonlinear.jl \
                  $(JL_CH08)/bvp_eigenvalue.jl \
                  $(JL_CH08)/bvp_2d_poisson.jl \
                  $(JL_CH08)/bvp_helmholtz.jl \
                  $(JL_CH08)/harmonic_oscillator.jl

# Julia scripts - Chapter 9 (Fourier Grids)
JL_CH09 = codes/julia/ch09
JL_SCRIPTS_CH09 = $(JL_CH09)/two_views_function.jl \
                  $(JL_CH09)/aliasing_demo.jl \
                  $(JL_CH09)/sinc_interpolation.jl \
                  $(JL_CH09)/fft_aliasing.jl \
                  $(JL_CH09)/smoothness_spectra.jl \
                  $(JL_CH09)/zero_padding_interpolation.jl

# Julia scripts - Chapter 10 (Spectral PDE Solvers)
JL_CH10 = codes/julia/ch10
JL_SCRIPTS_CH10 = $(JL_CH10)/chebfft.jl \
                  $(JL_CH10)/cheb_fourier_geometry.jl \
                  $(JL_CH10)/chebfft_accuracy.jl \
                  $(JL_CH10)/wave1d_cheb.jl \
                  $(JL_CH10)/wave2d_cheb.jl \
                  $(JL_CH10)/heat1d_cheb.jl \
                  $(JL_CH10)/heat2d_cheb.jl \
                  $(JL_CH10)/poisson2d_cheb.jl \
                  $(JL_CH10)/transport_variable.jl

# Figure outputs - Chapter 2
FIG_DIR_CH02 = textbook/figures/ch02
PY_FIGS_CH02 = $(FIG_DIR_CH02)/python/heat_evolution.pdf \
               $(FIG_DIR_CH02)/python/heat_waterfall.pdf \
               $(FIG_DIR_CH02)/python/wave_evolution.pdf \
               $(FIG_DIR_CH02)/python/wave_waterfall.pdf \
               $(FIG_DIR_CH02)/python/laplace_solution.pdf
M_FIGS_CH02 = $(FIG_DIR_CH02)/matlab/heat_evolution.pdf \
              $(FIG_DIR_CH02)/matlab/heat_waterfall.pdf \
              $(FIG_DIR_CH02)/matlab/wave_evolution.pdf \
              $(FIG_DIR_CH02)/matlab/wave_waterfall.pdf \
              $(FIG_DIR_CH02)/matlab/laplace_solution.pdf

# Figure outputs - Chapter 3
FIG_DIR_CH03 = textbook/figures/ch03
PY_FIGS_CH03 = $(FIG_DIR_CH03)/python/collocation_example1.pdf \
               $(FIG_DIR_CH03)/python/collocation_vs_galerkin.pdf
M_FIGS_CH03 = $(FIG_DIR_CH03)/matlab/collocation_example1.pdf \
              $(FIG_DIR_CH03)/matlab/collocation_vs_galerkin.pdf

# Figure outputs - Chapter 4
FIG_DIR_CH04 = textbook/figures/ch04
PY_FIGS_CH04 = $(FIG_DIR_CH04)/python/runge_phenomenon.pdf \
               $(FIG_DIR_CH04)/python/chebyshev_success.pdf \
               $(FIG_DIR_CH04)/python/chebyshev_points_circle.pdf \
               $(FIG_DIR_CH04)/python/equipotential_curves.pdf \
               $(FIG_DIR_CH04)/python/lagrange_basis.pdf \
               $(FIG_DIR_CH04)/python/lebesgue_functions.pdf \
               $(FIG_DIR_CH04)/python/lebesgue_constants_zoom.pdf \
               $(FIG_DIR_CH04)/python/lebesgue_random_nodes.pdf \
               $(FIG_DIR_CH04)/python/convergence_comparison.pdf \
               $(FIG_DIR_CH04)/python/convergence_zoom.pdf
M_FIGS_CH04 = $(FIG_DIR_CH04)/matlab/runge_phenomenon.pdf \
              $(FIG_DIR_CH04)/matlab/chebyshev_success.pdf \
              $(FIG_DIR_CH04)/matlab/chebyshev_points_circle.pdf \
              $(FIG_DIR_CH04)/matlab/equipotential_curves.pdf \
              $(FIG_DIR_CH04)/matlab/lagrange_basis.pdf \
              $(FIG_DIR_CH04)/matlab/lebesgue_functions.pdf \
              $(FIG_DIR_CH04)/matlab/lebesgue_random_nodes.pdf \
              $(FIG_DIR_CH04)/matlab/convergence_comparison.pdf

# Figure outputs - Chapter 5
FIG_DIR_CH05 = textbook/figures/ch05
PY_FIGS_CH05 = $(FIG_DIR_CH05)/python/fd_matrix_bandwidth.pdf \
               $(FIG_DIR_CH05)/python/spectral_matrix_structure.pdf \
               $(FIG_DIR_CH05)/python/fd_stencil_schematic.pdf \
               $(FIG_DIR_CH05)/python/stencil_pyramid.pdf \
               $(FIG_DIR_CH05)/python/convergence_comparison.pdf \
               $(FIG_DIR_CH05)/python/spectral_derivatives_demo.pdf
M_FIGS_CH05 = $(FIG_DIR_CH05)/matlab/fd_matrix_bandwidth.pdf \
              $(FIG_DIR_CH05)/matlab/spectral_matrix_structure.pdf \
              $(FIG_DIR_CH05)/matlab/fd_stencil_schematic.pdf \
              $(FIG_DIR_CH05)/matlab/stencil_pyramid.pdf \
              $(FIG_DIR_CH05)/matlab/convergence_comparison.pdf \
              $(FIG_DIR_CH05)/matlab/spectral_derivatives_demo.pdf

# Figure outputs - Chapter 6 (Smoothness and Spectral Accuracy)
FIG_DIR_CH06 = textbook/figures/ch06
PY_FIGS_CH06 = $(FIG_DIR_CH06)/python/decay_hierarchy.pdf \
               $(FIG_DIR_CH06)/python/aliasing_visualization.pdf \
               $(FIG_DIR_CH06)/python/convergence_rates.pdf

# Figure outputs - Chapter 7 (Chebyshev Differentiation)
FIG_DIR_CH07 = textbook/figures/ch07
PY_FIGS_CH07 = $(FIG_DIR_CH07)/python/grid_comparison.pdf \
               $(FIG_DIR_CH07)/python/cheb_matrix_structure.pdf \
               $(FIG_DIR_CH07)/python/cheb_cardinal.pdf \
               $(FIG_DIR_CH07)/python/cheb_diff_demo.pdf \
               $(FIG_DIR_CH07)/python/convergence_waterfall.pdf

# Figure outputs - Chapter 8 (Boundary Value Problems)
FIG_DIR_CH08 = textbook/figures/ch08
PY_FIGS_CH08 = $(FIG_DIR_CH08)/python/poisson_1d.pdf \
               $(FIG_DIR_CH08)/python/variable_coeff.pdf \
               $(FIG_DIR_CH08)/python/bratu.pdf \
               $(FIG_DIR_CH08)/python/eigenvalue_problem.pdf \
               $(FIG_DIR_CH08)/python/tensor_grid.pdf \
               $(FIG_DIR_CH08)/python/poisson_2d.pdf \
               $(FIG_DIR_CH08)/python/laplacian_sparsity.pdf \
               $(FIG_DIR_CH08)/python/helmholtz.pdf \
               $(FIG_DIR_CH08)/python/harmonic_oscillator.pdf

# MATLAB figure outputs - Chapter 6 (Smoothness and Spectral Accuracy)
M_FIGS_CH06 = $(FIG_DIR_CH06)/matlab/decay_hierarchy.pdf \
              $(FIG_DIR_CH06)/matlab/aliasing_visualization.pdf \
              $(FIG_DIR_CH06)/matlab/convergence_rates.pdf

# MATLAB figure outputs - Chapter 7 (Chebyshev Differentiation)
M_FIGS_CH07 = $(FIG_DIR_CH07)/matlab/grid_comparison.pdf \
              $(FIG_DIR_CH07)/matlab/cheb_matrix_structure.pdf \
              $(FIG_DIR_CH07)/matlab/cheb_cardinal.pdf \
              $(FIG_DIR_CH07)/matlab/cheb_diff_demo.pdf \
              $(FIG_DIR_CH07)/matlab/convergence_waterfall.pdf

# MATLAB figure outputs - Chapter 8 (Boundary Value Problems)
M_FIGS_CH08 = $(FIG_DIR_CH08)/matlab/poisson_1d.pdf \
              $(FIG_DIR_CH08)/matlab/variable_coeff.pdf \
              $(FIG_DIR_CH08)/matlab/bratu.pdf \
              $(FIG_DIR_CH08)/matlab/eigenvalue_problem.pdf \
              $(FIG_DIR_CH08)/matlab/tensor_grid.pdf \
              $(FIG_DIR_CH08)/matlab/poisson_2d.pdf \
              $(FIG_DIR_CH08)/matlab/laplacian_sparsity.pdf \
              $(FIG_DIR_CH08)/matlab/helmholtz.pdf \
              $(FIG_DIR_CH08)/matlab/harmonic_oscillator.pdf

# Figure outputs - Chapter 9 (Fourier Grids)
FIG_DIR_CH09 = textbook/figures/ch09
PY_FIGS_CH09 = $(FIG_DIR_CH09)/python/two_views_function.pdf \
               $(FIG_DIR_CH09)/python/aliasing_demo.pdf \
               $(FIG_DIR_CH09)/python/sinc_interpolation.pdf \
               $(FIG_DIR_CH09)/python/fft_aliasing.pdf \
               $(FIG_DIR_CH09)/python/smoothness_spectra.pdf \
               $(FIG_DIR_CH09)/python/zero_padding_interpolation.pdf

# MATLAB figure outputs - Chapter 9 (Fourier Grids)
M_FIGS_CH09 = $(FIG_DIR_CH09)/matlab/two_views_function.pdf \
              $(FIG_DIR_CH09)/matlab/aliasing_demo.pdf \
              $(FIG_DIR_CH09)/matlab/sinc_interpolation.pdf \
              $(FIG_DIR_CH09)/matlab/fft_aliasing.pdf \
              $(FIG_DIR_CH09)/matlab/smoothness_spectra.pdf \
              $(FIG_DIR_CH09)/matlab/zero_padding_interpolation.pdf

# Figure outputs - Chapter 10 (Spectral PDE Solvers)
FIG_DIR_CH10 = textbook/figures/ch10
PY_FIGS_CH10 = $(FIG_DIR_CH10)/python/cheb_fourier_geometry.pdf \
               $(FIG_DIR_CH10)/python/chebfft_accuracy.pdf \
               $(FIG_DIR_CH10)/python/wave1d_waterfall.pdf \
               $(FIG_DIR_CH10)/python/wave2d_snapshots.pdf \
               $(FIG_DIR_CH10)/python/heat1d_evolution.pdf \
               $(FIG_DIR_CH10)/python/heat2d_snapshots.pdf \
               $(FIG_DIR_CH10)/python/poisson2d_solution.pdf \
               $(FIG_DIR_CH10)/python/transport_variable.pdf

# MATLAB figure outputs - Chapter 10 (Spectral PDE Solvers)
M_FIGS_CH10 = $(FIG_DIR_CH10)/matlab/cheb_fourier_geometry.pdf \
              $(FIG_DIR_CH10)/matlab/chebfft_accuracy.pdf \
              $(FIG_DIR_CH10)/matlab/wave1d_waterfall.pdf \
              $(FIG_DIR_CH10)/matlab/wave2d_snapshots.pdf \
              $(FIG_DIR_CH10)/matlab/heat1d_evolution.pdf \
              $(FIG_DIR_CH10)/matlab/heat2d_snapshots.pdf \
              $(FIG_DIR_CH10)/matlab/poisson2d_solution.pdf \
              $(FIG_DIR_CH10)/matlab/transport_variable.pdf

# Julia figure outputs - Chapter 2
JL_FIGS_CH02 = $(FIG_DIR_CH02)/julia/heat_evolution.pdf \
               $(FIG_DIR_CH02)/julia/heat_waterfall.pdf \
               $(FIG_DIR_CH02)/julia/wave_evolution.pdf \
               $(FIG_DIR_CH02)/julia/wave_waterfall.pdf \
               $(FIG_DIR_CH02)/julia/laplace_solution.pdf

# Julia figure outputs - Chapter 3
JL_FIGS_CH03 = $(FIG_DIR_CH03)/julia/collocation_example1.pdf \
               $(FIG_DIR_CH03)/julia/collocation_vs_galerkin.pdf

# Julia figure outputs - Chapter 4
JL_FIGS_CH04 = $(FIG_DIR_CH04)/julia/runge_phenomenon.pdf \
               $(FIG_DIR_CH04)/julia/chebyshev_success.pdf \
               $(FIG_DIR_CH04)/julia/chebyshev_points_circle.pdf \
               $(FIG_DIR_CH04)/julia/equipotential_curves.pdf \
               $(FIG_DIR_CH04)/julia/lagrange_basis.pdf \
               $(FIG_DIR_CH04)/julia/lebesgue_functions.pdf \
               $(FIG_DIR_CH04)/julia/lebesgue_constants_zoom.pdf \
               $(FIG_DIR_CH04)/julia/lebesgue_random_nodes.pdf \
               $(FIG_DIR_CH04)/julia/convergence_comparison.pdf \
               $(FIG_DIR_CH04)/julia/convergence_zoom.pdf

# Julia figure outputs - Chapter 5
JL_FIGS_CH05 = $(FIG_DIR_CH05)/julia/fd_matrix_bandwidth.pdf \
               $(FIG_DIR_CH05)/julia/spectral_matrix_structure.pdf \
               $(FIG_DIR_CH05)/julia/fd_stencil_schematic.pdf \
               $(FIG_DIR_CH05)/julia/stencil_pyramid.pdf \
               $(FIG_DIR_CH05)/julia/convergence_comparison.pdf \
               $(FIG_DIR_CH05)/julia/spectral_derivatives_demo.pdf

# Julia figure outputs - Chapter 6 (Smoothness and Spectral Accuracy)
JL_FIGS_CH06 = $(FIG_DIR_CH06)/julia/decay_hierarchy.pdf \
               $(FIG_DIR_CH06)/julia/aliasing_visualization.pdf \
               $(FIG_DIR_CH06)/julia/convergence_rates.pdf

# Julia figure outputs - Chapter 7 (Chebyshev Differentiation)
JL_FIGS_CH07 = $(FIG_DIR_CH07)/julia/grid_comparison.pdf \
               $(FIG_DIR_CH07)/julia/cheb_matrix_structure.pdf \
               $(FIG_DIR_CH07)/julia/cheb_cardinal.pdf \
               $(FIG_DIR_CH07)/julia/cheb_diff_demo.pdf \
               $(FIG_DIR_CH07)/julia/convergence_waterfall.pdf

# Julia figure outputs - Chapter 8 (Boundary Value Problems)
JL_FIGS_CH08 = $(FIG_DIR_CH08)/julia/poisson_1d.pdf \
               $(FIG_DIR_CH08)/julia/variable_coeff.pdf \
               $(FIG_DIR_CH08)/julia/bratu.pdf \
               $(FIG_DIR_CH08)/julia/eigenvalue_problem.pdf \
               $(FIG_DIR_CH08)/julia/tensor_grid.pdf \
               $(FIG_DIR_CH08)/julia/poisson_2d.pdf \
               $(FIG_DIR_CH08)/julia/laplacian_sparsity.pdf \
               $(FIG_DIR_CH08)/julia/helmholtz.pdf \
               $(FIG_DIR_CH08)/julia/harmonic_oscillator.pdf

# Julia figure outputs - Chapter 9 (Fourier Grids)
JL_FIGS_CH09 = $(FIG_DIR_CH09)/julia/two_views_function.pdf \
               $(FIG_DIR_CH09)/julia/aliasing_demo.pdf \
               $(FIG_DIR_CH09)/julia/sinc_interpolation.pdf \
               $(FIG_DIR_CH09)/julia/fft_aliasing.pdf \
               $(FIG_DIR_CH09)/julia/smoothness_spectra.pdf \
               $(FIG_DIR_CH09)/julia/zero_padding_interpolation.pdf

# Julia figure outputs - Chapter 10 (Spectral PDE Solvers)
JL_FIGS_CH10 = $(FIG_DIR_CH10)/julia/cheb_fourier_geometry.pdf \
               $(FIG_DIR_CH10)/julia/chebfft_accuracy.pdf \
               $(FIG_DIR_CH10)/julia/wave1d_waterfall.pdf \
               $(FIG_DIR_CH10)/julia/wave2d_snapshots.pdf \
               $(FIG_DIR_CH10)/julia/heat1d_evolution.pdf \
               $(FIG_DIR_CH10)/julia/heat2d_snapshots.pdf \
               $(FIG_DIR_CH10)/julia/poisson2d_solution.pdf \
               $(FIG_DIR_CH10)/julia/transport_variable.pdf

# Combined figure variables
PY_FIGS = $(PY_FIGS_CH02) $(PY_FIGS_CH03) $(PY_FIGS_CH04) $(PY_FIGS_CH05) $(PY_FIGS_CH06) $(PY_FIGS_CH07) $(PY_FIGS_CH08) $(PY_FIGS_CH09) $(PY_FIGS_CH10)
M_FIGS = $(M_FIGS_CH02) $(M_FIGS_CH03) $(M_FIGS_CH04) $(M_FIGS_CH05) $(M_FIGS_CH06) $(M_FIGS_CH07) $(M_FIGS_CH08) $(M_FIGS_CH09) $(M_FIGS_CH10)
JL_FIGS = $(JL_FIGS_CH02) $(JL_FIGS_CH03) $(JL_FIGS_CH04) $(JL_FIGS_CH05) $(JL_FIGS_CH06) $(JL_FIGS_CH07) $(JL_FIGS_CH08) $(JL_FIGS_CH09) $(JL_FIGS_CH10)

# Default target: build everything
all: figures textbook tplan

# Build textbook (depends on figures)
textbook: $(OUT)

$(OUT): $(SRC) textbook/chapters/preface.typ textbook/chapters/introduction.typ textbook/chapters/classical_pdes.typ textbook/chapters/mise_en_bouche.typ textbook/chapters/geometry_of_nodes.typ textbook/chapters/differentiation_matrices.typ textbook/chapters/smoothness_accuracy.typ textbook/chapters/chebyshev_differentiation.typ textbook/chapters/boundary_value_problems.typ textbook/chapters/fourier_grids.typ textbook/chapters/spectral_pde_solvers.typ textbook/styles/template.typ $(PY_FIGS)
	mkdir -p $(OUT_DIR)
	$(TYPST) compile $(SRC) $(OUT)

# Generate all figures (Python is primary, MATLAB is optional)
figures: figures-python

figures-python: $(PY_FIGS)

figures-matlab: $(M_FIGS)

figures-julia: $(JL_FIGS)

# Python figure generation rules - Chapter 2
$(FIG_DIR_CH02)/python/heat_evolution.pdf: $(PY_CH02)/heat_equation_evolution.py
	@mkdir -p $(FIG_DIR_CH02)/python
	$(PYTHON) $<

$(FIG_DIR_CH02)/python/wave_evolution.pdf: $(PY_CH02)/wave_equation_evolution.py
	@mkdir -p $(FIG_DIR_CH02)/python
	$(PYTHON) $<

$(FIG_DIR_CH02)/python/laplace_solution.pdf: $(PY_CH02)/laplace_equation_2d.py
	@mkdir -p $(FIG_DIR_CH02)/python
	$(PYTHON) $<

$(FIG_DIR_CH02)/python/heat_waterfall.pdf: $(PY_CH02)/heat_equation_waterfall.py
	@mkdir -p $(FIG_DIR_CH02)/python
	$(PYTHON) $<

$(FIG_DIR_CH02)/python/wave_waterfall.pdf: $(PY_CH02)/wave_equation_waterfall.py
	@mkdir -p $(FIG_DIR_CH02)/python
	$(PYTHON) $<

# Python figure generation rules - Chapter 3
$(FIG_DIR_CH03)/python/collocation_example1.pdf: $(PY_CH03)/collocation_example1.py
	@mkdir -p $(FIG_DIR_CH03)/python
	$(PYTHON) $<

$(FIG_DIR_CH03)/python/collocation_vs_galerkin.pdf: $(PY_CH03)/collocation_vs_galerkin.py
	@mkdir -p $(FIG_DIR_CH03)/python
	$(PYTHON) $<

# Python figure generation rules - Chapter 4
$(FIG_DIR_CH04)/python/runge_phenomenon.pdf: $(PY_CH04)/runge_phenomenon.py
	@mkdir -p $(FIG_DIR_CH04)/python
	$(PYTHON) $<

$(FIG_DIR_CH04)/python/chebyshev_success.pdf: $(PY_CH04)/chebyshev_success.py
	@mkdir -p $(FIG_DIR_CH04)/python
	$(PYTHON) $<

$(FIG_DIR_CH04)/python/chebyshev_points_circle.pdf: $(PY_CH04)/chebyshev_points_circle.py
	@mkdir -p $(FIG_DIR_CH04)/python
	$(PYTHON) $<

$(FIG_DIR_CH04)/python/equipotential_curves.pdf: $(PY_CH04)/equipotential_curves.py
	@mkdir -p $(FIG_DIR_CH04)/python
	$(PYTHON) $<

$(FIG_DIR_CH04)/python/lagrange_basis.pdf: $(PY_CH04)/lagrange_basis.py
	@mkdir -p $(FIG_DIR_CH04)/python
	$(PYTHON) $<

$(FIG_DIR_CH04)/python/lebesgue_functions.pdf: $(PY_CH04)/lebesgue_functions.py
	@mkdir -p $(FIG_DIR_CH04)/python
	$(PYTHON) $<

$(FIG_DIR_CH04)/python/lebesgue_constants_zoom.pdf: $(PY_CH04)/lebesgue_constants_zoom.py
	@mkdir -p $(FIG_DIR_CH04)/python
	$(PYTHON) $<

$(FIG_DIR_CH04)/python/convergence_comparison.pdf: $(PY_CH04)/convergence_comparison.py
	@mkdir -p $(FIG_DIR_CH04)/python
	$(PYTHON) $<

$(FIG_DIR_CH04)/python/convergence_zoom.pdf: $(PY_CH04)/convergence_zoom.py
	@mkdir -p $(FIG_DIR_CH04)/python
	$(PYTHON) $<

$(FIG_DIR_CH04)/python/lebesgue_random_nodes.pdf: $(PY_CH04)/lebesgue_random_nodes.py
	@mkdir -p $(FIG_DIR_CH04)/python
	$(PYTHON) $<

# Python figure generation rules - Chapter 5
$(FIG_DIR_CH05)/python/fd_matrix_bandwidth.pdf: $(PY_CH05)/fd_matrix_bandwidth.py
	@mkdir -p $(FIG_DIR_CH05)/python
	$(PYTHON) $<

$(FIG_DIR_CH05)/python/spectral_matrix_structure.pdf: $(PY_CH05)/spectral_matrix_structure.py
	@mkdir -p $(FIG_DIR_CH05)/python
	$(PYTHON) $<

$(FIG_DIR_CH05)/python/fd_stencil_schematic.pdf: $(PY_CH05)/fd_stencil_schematic.py
	@mkdir -p $(FIG_DIR_CH05)/python
	$(PYTHON) $<

$(FIG_DIR_CH05)/python/stencil_pyramid.pdf: $(PY_CH05)/stencil_pyramid.py
	@mkdir -p $(FIG_DIR_CH05)/python
	$(PYTHON) $<

$(FIG_DIR_CH05)/python/convergence_comparison.pdf: $(PY_CH05)/convergence_comparison.py
	@mkdir -p $(FIG_DIR_CH05)/python
	$(PYTHON) $<

$(FIG_DIR_CH05)/python/spectral_derivatives_demo.pdf: $(PY_CH05)/spectral_derivatives_demo.py
	@mkdir -p $(FIG_DIR_CH05)/python
	$(PYTHON) $<

# Python figure generation rules - Chapter 6 (Smoothness and Spectral Accuracy)
$(FIG_DIR_CH06)/python/decay_hierarchy.pdf: $(PY_CH06)/fourier_decay.py
	@mkdir -p $(FIG_DIR_CH06)/python
	$(PYTHON) $<

$(FIG_DIR_CH06)/python/aliasing_visualization.pdf: $(PY_CH06)/aliasing_demo.py
	@mkdir -p $(FIG_DIR_CH06)/python
	$(PYTHON) $<

$(FIG_DIR_CH06)/python/convergence_rates.pdf: $(PY_CH06)/convergence_rates.py
	@mkdir -p $(FIG_DIR_CH06)/python
	$(PYTHON) $<

# Python figure generation rules - Chapter 7 (Chebyshev Differentiation)
$(FIG_DIR_CH07)/python/grid_comparison.pdf: $(PY_CH07)/cheb_grid_comparison.py $(PY_CH07)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH07)/python
	$(PYTHON) $<

$(FIG_DIR_CH07)/python/cheb_matrix_structure.pdf: $(PY_CH07)/cheb_matrix_structure.py $(PY_CH07)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH07)/python
	$(PYTHON) $<

$(FIG_DIR_CH07)/python/cheb_cardinal.pdf: $(PY_CH07)/cheb_cardinal.py $(PY_CH07)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH07)/python
	$(PYTHON) $<

$(FIG_DIR_CH07)/python/cheb_diff_demo.pdf: $(PY_CH07)/cheb_diff_demo.py $(PY_CH07)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH07)/python
	$(PYTHON) $<

$(FIG_DIR_CH07)/python/convergence_waterfall.pdf: $(PY_CH07)/cheb_convergence.py $(PY_CH07)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH07)/python
	$(PYTHON) $<

# Python figure generation rules - Chapter 8 (Boundary Value Problems)
$(FIG_DIR_CH08)/python/poisson_1d.pdf: $(PY_CH08)/bvp_linear.py $(PY_CH07)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH08)/python
	$(PYTHON) $<

$(FIG_DIR_CH08)/python/variable_coeff.pdf: $(PY_CH08)/bvp_variable_coeff.py $(PY_CH07)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH08)/python
	$(PYTHON) $<

$(FIG_DIR_CH08)/python/bratu.pdf: $(PY_CH08)/bvp_nonlinear.py $(PY_CH07)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH08)/python
	$(PYTHON) $<

$(FIG_DIR_CH08)/python/eigenvalue_problem.pdf: $(PY_CH08)/bvp_eigenvalue.py $(PY_CH07)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH08)/python
	$(PYTHON) $<

$(FIG_DIR_CH08)/python/tensor_grid.pdf: $(PY_CH08)/bvp_2d_poisson.py $(PY_CH07)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH08)/python
	$(PYTHON) $<

$(FIG_DIR_CH08)/python/poisson_2d.pdf: $(PY_CH08)/bvp_2d_poisson.py $(PY_CH07)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH08)/python
	$(PYTHON) $<

$(FIG_DIR_CH08)/python/laplacian_sparsity.pdf: $(PY_CH08)/bvp_2d_poisson.py $(PY_CH07)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH08)/python
	$(PYTHON) $<

$(FIG_DIR_CH08)/python/helmholtz.pdf: $(PY_CH08)/bvp_helmholtz.py $(PY_CH07)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH08)/python
	$(PYTHON) $<

$(FIG_DIR_CH08)/python/harmonic_oscillator.pdf: $(PY_CH08)/harmonic_oscillator.py
	@mkdir -p $(FIG_DIR_CH08)/python
	$(PYTHON) $<

# Python figure generation rules - Chapter 9 (Fourier Grids)
$(FIG_DIR_CH09)/python/two_views_function.pdf: $(PY_CH09)/two_views_function.py
	@mkdir -p $(FIG_DIR_CH09)/python
	$(PYTHON) $<

$(FIG_DIR_CH09)/python/aliasing_demo.pdf: $(PY_CH09)/aliasing_demo.py
	@mkdir -p $(FIG_DIR_CH09)/python
	$(PYTHON) $<

$(FIG_DIR_CH09)/python/sinc_interpolation.pdf: $(PY_CH09)/sinc_interpolation.py
	@mkdir -p $(FIG_DIR_CH09)/python
	$(PYTHON) $<

$(FIG_DIR_CH09)/python/fft_aliasing.pdf: $(PY_CH09)/fft_aliasing.py
	@mkdir -p $(FIG_DIR_CH09)/python
	$(PYTHON) $<

$(FIG_DIR_CH09)/python/smoothness_spectra.pdf: $(PY_CH09)/smoothness_spectra.py
	@mkdir -p $(FIG_DIR_CH09)/python
	$(PYTHON) $<

$(FIG_DIR_CH09)/python/zero_padding_interpolation.pdf: $(PY_CH09)/zero_padding_interpolation.py
	@mkdir -p $(FIG_DIR_CH09)/python
	$(PYTHON) $<

# Python figure generation rules - Chapter 10 (Spectral PDE Solvers)
$(FIG_DIR_CH10)/python/cheb_fourier_geometry.pdf: $(PY_CH10)/cheb_fourier_geometry.py $(PY_CH10)/chebfft.py
	@mkdir -p $(FIG_DIR_CH10)/python
	$(PYTHON) $<

$(FIG_DIR_CH10)/python/chebfft_accuracy.pdf: $(PY_CH10)/chebfft_accuracy.py $(PY_CH10)/chebfft.py
	@mkdir -p $(FIG_DIR_CH10)/python
	$(PYTHON) $<

$(FIG_DIR_CH10)/python/wave1d_waterfall.pdf: $(PY_CH10)/wave1d_cheb.py $(PY_CH10)/chebfft.py
	@mkdir -p $(FIG_DIR_CH10)/python
	$(PYTHON) $<

$(FIG_DIR_CH10)/python/wave2d_snapshots.pdf: $(PY_CH10)/wave2d_cheb.py $(PY_CH10)/chebfft.py
	@mkdir -p $(FIG_DIR_CH10)/python
	$(PYTHON) $<

$(FIG_DIR_CH10)/python/heat1d_evolution.pdf: $(PY_CH10)/heat1d_cheb.py $(PY_CH10)/chebfft.py
	@mkdir -p $(FIG_DIR_CH10)/python
	$(PYTHON) $<

$(FIG_DIR_CH10)/python/heat2d_snapshots.pdf: $(PY_CH10)/heat2d_cheb.py $(PY_CH10)/chebfft.py
	@mkdir -p $(FIG_DIR_CH10)/python
	$(PYTHON) $<

$(FIG_DIR_CH10)/python/poisson2d_solution.pdf: $(PY_CH10)/poisson2d_cheb.py $(PY_CH10)/chebfft.py
	@mkdir -p $(FIG_DIR_CH10)/python
	$(PYTHON) $<

$(FIG_DIR_CH10)/python/transport_variable.pdf: $(PY_CH10)/transport_variable.py
	@mkdir -p $(FIG_DIR_CH10)/python
	$(PYTHON) $<

# MATLAB figure generation rules - Chapter 2
$(FIG_DIR_CH02)/matlab/heat_evolution.pdf: $(M_CH02)/heat_equation_evolution.m
	@mkdir -p $(FIG_DIR_CH02)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH02)/matlab/heat_waterfall.pdf: $(M_CH02)/heat_equation_waterfall.m
	@mkdir -p $(FIG_DIR_CH02)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH02)/matlab/wave_evolution.pdf: $(M_CH02)/wave_equation_evolution.m
	@mkdir -p $(FIG_DIR_CH02)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH02)/matlab/wave_waterfall.pdf: $(M_CH02)/wave_equation_waterfall.m
	@mkdir -p $(FIG_DIR_CH02)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH02)/matlab/laplace_solution.pdf: $(M_CH02)/laplace_equation_2d.m
	@mkdir -p $(FIG_DIR_CH02)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

# MATLAB figure generation rules - Chapter 3
$(FIG_DIR_CH03)/matlab/collocation_example1.pdf: $(M_CH03)/collocation_example1.m
	@mkdir -p $(FIG_DIR_CH03)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH03)/matlab/collocation_vs_galerkin.pdf: $(M_CH03)/collocation_vs_galerkin.m
	@mkdir -p $(FIG_DIR_CH03)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

# MATLAB figure generation rules - Chapter 4
$(FIG_DIR_CH04)/matlab/runge_phenomenon.pdf: $(M_CH04)/runge_phenomenon.m
	@mkdir -p $(FIG_DIR_CH04)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH04)/matlab/chebyshev_success.pdf: $(M_CH04)/chebyshev_success.m
	@mkdir -p $(FIG_DIR_CH04)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH04)/matlab/chebyshev_points_circle.pdf: $(M_CH04)/chebyshev_points_circle.m
	@mkdir -p $(FIG_DIR_CH04)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH04)/matlab/equipotential_curves.pdf: $(M_CH04)/equipotential_curves.m
	@mkdir -p $(FIG_DIR_CH04)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH04)/matlab/lagrange_basis.pdf: $(M_CH04)/lagrange_basis.m
	@mkdir -p $(FIG_DIR_CH04)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH04)/matlab/lebesgue_functions.pdf: $(M_CH04)/lebesgue_functions.m
	@mkdir -p $(FIG_DIR_CH04)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH04)/matlab/convergence_comparison.pdf: $(M_CH04)/convergence_comparison.m
	@mkdir -p $(FIG_DIR_CH04)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH04)/matlab/lebesgue_random_nodes.pdf: $(M_CH04)/lebesgue_random_nodes.m
	@mkdir -p $(FIG_DIR_CH04)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

# MATLAB figure generation rules - Chapter 5
$(FIG_DIR_CH05)/matlab/fd_matrix_bandwidth.pdf: $(M_CH05)/fd_matrix_bandwidth.m
	@mkdir -p $(FIG_DIR_CH05)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH05)/matlab/spectral_matrix_structure.pdf: $(M_CH05)/spectral_matrix_structure.m
	@mkdir -p $(FIG_DIR_CH05)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH05)/matlab/fd_stencil_schematic.pdf: $(M_CH05)/fd_stencil_schematic.m
	@mkdir -p $(FIG_DIR_CH05)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH05)/matlab/stencil_pyramid.pdf: $(M_CH05)/stencil_pyramid.m
	@mkdir -p $(FIG_DIR_CH05)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH05)/matlab/convergence_comparison.pdf: $(M_CH05)/convergence_comparison.m
	@mkdir -p $(FIG_DIR_CH05)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH05)/matlab/spectral_derivatives_demo.pdf: $(M_CH05)/spectral_derivatives_demo.m
	@mkdir -p $(FIG_DIR_CH05)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

# MATLAB figure generation rules - Chapter 6 (Smoothness and Spectral Accuracy)
$(FIG_DIR_CH06)/matlab/decay_hierarchy.pdf: $(M_CH06)/fourier_decay.m
	@mkdir -p $(FIG_DIR_CH06)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH06)/matlab/aliasing_visualization.pdf: $(M_CH06)/aliasing_demo.m
	@mkdir -p $(FIG_DIR_CH06)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH06)/matlab/convergence_rates.pdf: $(M_CH06)/convergence_rates.m
	@mkdir -p $(FIG_DIR_CH06)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

# MATLAB figure generation rules - Chapter 7 (Chebyshev Differentiation)
$(FIG_DIR_CH07)/matlab/grid_comparison.pdf: $(M_CH07)/cheb_grid_comparison.m $(M_CH07)/cheb_matrix.m
	@mkdir -p $(FIG_DIR_CH07)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH07)/matlab/cheb_matrix_structure.pdf: $(M_CH07)/cheb_matrix_structure.m $(M_CH07)/cheb_matrix.m
	@mkdir -p $(FIG_DIR_CH07)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH07)/matlab/cheb_cardinal.pdf: $(M_CH07)/cheb_cardinal.m $(M_CH07)/cheb_matrix.m
	@mkdir -p $(FIG_DIR_CH07)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH07)/matlab/cheb_diff_demo.pdf: $(M_CH07)/cheb_diff_demo.m $(M_CH07)/cheb_matrix.m
	@mkdir -p $(FIG_DIR_CH07)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH07)/matlab/convergence_waterfall.pdf: $(M_CH07)/cheb_convergence.m $(M_CH07)/cheb_matrix.m
	@mkdir -p $(FIG_DIR_CH07)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

# MATLAB figure generation rules - Chapter 8 (Boundary Value Problems)
$(FIG_DIR_CH08)/matlab/poisson_1d.pdf: $(M_CH08)/bvp_linear.m $(M_CH07)/cheb_matrix.m
	@mkdir -p $(FIG_DIR_CH08)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH08)/matlab/variable_coeff.pdf: $(M_CH08)/bvp_variable_coeff.m $(M_CH07)/cheb_matrix.m
	@mkdir -p $(FIG_DIR_CH08)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH08)/matlab/bratu.pdf: $(M_CH08)/bvp_nonlinear.m $(M_CH07)/cheb_matrix.m
	@mkdir -p $(FIG_DIR_CH08)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH08)/matlab/eigenvalue_problem.pdf: $(M_CH08)/bvp_eigenvalue.m $(M_CH07)/cheb_matrix.m
	@mkdir -p $(FIG_DIR_CH08)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH08)/matlab/tensor_grid.pdf: $(M_CH08)/bvp_2d_poisson.m $(M_CH07)/cheb_matrix.m
	@mkdir -p $(FIG_DIR_CH08)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH08)/matlab/poisson_2d.pdf: $(M_CH08)/bvp_2d_poisson.m $(M_CH07)/cheb_matrix.m
	@mkdir -p $(FIG_DIR_CH08)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH08)/matlab/laplacian_sparsity.pdf: $(M_CH08)/bvp_2d_poisson.m $(M_CH07)/cheb_matrix.m
	@mkdir -p $(FIG_DIR_CH08)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH08)/matlab/helmholtz.pdf: $(M_CH08)/bvp_helmholtz.m $(M_CH07)/cheb_matrix.m
	@mkdir -p $(FIG_DIR_CH08)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH08)/matlab/harmonic_oscillator.pdf: $(M_CH08)/harmonic_oscillator.m
	@mkdir -p $(FIG_DIR_CH08)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

# MATLAB figure generation rules - Chapter 9 (Fourier Grids)
$(FIG_DIR_CH09)/matlab/two_views_function.pdf: $(M_CH09)/two_views_function.m
	@mkdir -p $(FIG_DIR_CH09)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH09)/matlab/aliasing_demo.pdf: $(M_CH09)/aliasing_demo.m
	@mkdir -p $(FIG_DIR_CH09)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH09)/matlab/sinc_interpolation.pdf: $(M_CH09)/sinc_interpolation.m
	@mkdir -p $(FIG_DIR_CH09)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH09)/matlab/fft_aliasing.pdf: $(M_CH09)/fft_aliasing.m
	@mkdir -p $(FIG_DIR_CH09)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH09)/matlab/smoothness_spectra.pdf: $(M_CH09)/smoothness_spectra.m
	@mkdir -p $(FIG_DIR_CH09)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH09)/matlab/zero_padding_interpolation.pdf: $(M_CH09)/zero_padding_interpolation.m
	@mkdir -p $(FIG_DIR_CH09)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

# MATLAB figure generation rules - Chapter 10 (Spectral PDE Solvers)
$(FIG_DIR_CH10)/matlab/cheb_fourier_geometry.pdf: $(M_CH10)/cheb_fourier_geometry.m $(M_CH10)/chebfft.m $(M_CH10)/cheb_matrix.m
	@mkdir -p $(FIG_DIR_CH10)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH10)/matlab/chebfft_accuracy.pdf: $(M_CH10)/chebfft_accuracy.m $(M_CH10)/chebfft.m $(M_CH10)/cheb_matrix.m
	@mkdir -p $(FIG_DIR_CH10)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH10)/matlab/wave1d_waterfall.pdf: $(M_CH10)/wave1d_cheb.m $(M_CH10)/chebfft.m $(M_CH10)/cheb_matrix.m
	@mkdir -p $(FIG_DIR_CH10)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH10)/matlab/wave2d_snapshots.pdf: $(M_CH10)/wave2d_cheb.m $(M_CH10)/chebfft.m $(M_CH10)/cheb_matrix.m
	@mkdir -p $(FIG_DIR_CH10)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH10)/matlab/heat1d_evolution.pdf: $(M_CH10)/heat1d_cheb.m $(M_CH10)/chebfft.m $(M_CH10)/cheb_matrix.m
	@mkdir -p $(FIG_DIR_CH10)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH10)/matlab/heat2d_snapshots.pdf: $(M_CH10)/heat2d_cheb.m $(M_CH10)/chebfft.m $(M_CH10)/cheb_matrix.m
	@mkdir -p $(FIG_DIR_CH10)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH10)/matlab/poisson2d_solution.pdf: $(M_CH10)/poisson2d_cheb.m $(M_CH10)/chebfft.m $(M_CH10)/cheb_matrix.m
	@mkdir -p $(FIG_DIR_CH10)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR_CH10)/matlab/transport_variable.pdf: $(M_CH10)/transport_variable.m
	@mkdir -p $(FIG_DIR_CH10)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

# Julia figure generation rules - Chapter 2
$(FIG_DIR_CH02)/julia/heat_evolution.pdf: $(JL_CH02)/heat_equation_evolution.jl
	@mkdir -p $(FIG_DIR_CH02)/julia
	$(JULIA) $<

$(FIG_DIR_CH02)/julia/heat_waterfall.pdf: $(JL_CH02)/heat_equation_waterfall.jl
	@mkdir -p $(FIG_DIR_CH02)/julia
	$(JULIA) $<

$(FIG_DIR_CH02)/julia/wave_evolution.pdf: $(JL_CH02)/wave_equation_evolution.jl
	@mkdir -p $(FIG_DIR_CH02)/julia
	$(JULIA) $<

$(FIG_DIR_CH02)/julia/wave_waterfall.pdf: $(JL_CH02)/wave_equation_waterfall.jl
	@mkdir -p $(FIG_DIR_CH02)/julia
	$(JULIA) $<

$(FIG_DIR_CH02)/julia/laplace_solution.pdf: $(JL_CH02)/laplace_equation_2d.jl
	@mkdir -p $(FIG_DIR_CH02)/julia
	$(JULIA) $<

# Julia figure generation rules - Chapter 3
$(FIG_DIR_CH03)/julia/collocation_example1.pdf: $(JL_CH03)/collocation_example1.jl
	@mkdir -p $(FIG_DIR_CH03)/julia
	$(JULIA) $<

$(FIG_DIR_CH03)/julia/collocation_vs_galerkin.pdf: $(JL_CH03)/collocation_vs_galerkin.jl
	@mkdir -p $(FIG_DIR_CH03)/julia
	$(JULIA) $<

# Julia figure generation rules - Chapter 4
$(FIG_DIR_CH04)/julia/runge_phenomenon.pdf: $(JL_CH04)/runge_phenomenon.jl
	@mkdir -p $(FIG_DIR_CH04)/julia
	$(JULIA) $<

$(FIG_DIR_CH04)/julia/chebyshev_success.pdf: $(JL_CH04)/chebyshev_success.jl
	@mkdir -p $(FIG_DIR_CH04)/julia
	$(JULIA) $<

$(FIG_DIR_CH04)/julia/chebyshev_points_circle.pdf: $(JL_CH04)/chebyshev_points_circle.jl
	@mkdir -p $(FIG_DIR_CH04)/julia
	$(JULIA) $<

$(FIG_DIR_CH04)/julia/equipotential_curves.pdf: $(JL_CH04)/equipotential_curves.jl
	@mkdir -p $(FIG_DIR_CH04)/julia
	$(JULIA) $<

$(FIG_DIR_CH04)/julia/lagrange_basis.pdf: $(JL_CH04)/lagrange_basis.jl
	@mkdir -p $(FIG_DIR_CH04)/julia
	$(JULIA) $<

$(FIG_DIR_CH04)/julia/lebesgue_functions.pdf: $(JL_CH04)/lebesgue_functions.jl
	@mkdir -p $(FIG_DIR_CH04)/julia
	$(JULIA) $<

$(FIG_DIR_CH04)/julia/lebesgue_constants_zoom.pdf: $(JL_CH04)/lebesgue_constants_zoom.jl
	@mkdir -p $(FIG_DIR_CH04)/julia
	$(JULIA) $<

$(FIG_DIR_CH04)/julia/convergence_comparison.pdf: $(JL_CH04)/convergence_comparison.jl
	@mkdir -p $(FIG_DIR_CH04)/julia
	$(JULIA) $<

$(FIG_DIR_CH04)/julia/convergence_zoom.pdf: $(JL_CH04)/convergence_zoom.jl
	@mkdir -p $(FIG_DIR_CH04)/julia
	$(JULIA) $<

$(FIG_DIR_CH04)/julia/lebesgue_random_nodes.pdf: $(JL_CH04)/lebesgue_random_nodes.jl
	@mkdir -p $(FIG_DIR_CH04)/julia
	$(JULIA) $<

# Julia figure generation rules - Chapter 5
$(FIG_DIR_CH05)/julia/fd_matrix_bandwidth.pdf: $(JL_CH05)/fd_matrix_bandwidth.jl
	@mkdir -p $(FIG_DIR_CH05)/julia
	$(JULIA) $<

$(FIG_DIR_CH05)/julia/spectral_matrix_structure.pdf: $(JL_CH05)/spectral_matrix_structure.jl
	@mkdir -p $(FIG_DIR_CH05)/julia
	$(JULIA) $<

$(FIG_DIR_CH05)/julia/fd_stencil_schematic.pdf: $(JL_CH05)/fd_stencil_schematic.jl
	@mkdir -p $(FIG_DIR_CH05)/julia
	$(JULIA) $<

$(FIG_DIR_CH05)/julia/stencil_pyramid.pdf: $(JL_CH05)/stencil_pyramid.jl
	@mkdir -p $(FIG_DIR_CH05)/julia
	$(JULIA) $<

$(FIG_DIR_CH05)/julia/convergence_comparison.pdf: $(JL_CH05)/convergence_comparison.jl $(JL_CH05)/fdweights.jl
	@mkdir -p $(FIG_DIR_CH05)/julia
	$(JULIA) $<

$(FIG_DIR_CH05)/julia/spectral_derivatives_demo.pdf: $(JL_CH05)/spectral_derivatives_demo.jl $(JL_CH05)/fdweights.jl
	@mkdir -p $(FIG_DIR_CH05)/julia
	$(JULIA) $<

# Julia figure generation rules - Chapter 6 (Smoothness and Spectral Accuracy)
$(FIG_DIR_CH06)/julia/decay_hierarchy.pdf: $(JL_CH06)/fourier_decay.jl
	@mkdir -p $(FIG_DIR_CH06)/julia
	$(JULIA) $<

$(FIG_DIR_CH06)/julia/aliasing_visualization.pdf: $(JL_CH06)/aliasing_demo.jl
	@mkdir -p $(FIG_DIR_CH06)/julia
	$(JULIA) $<

$(FIG_DIR_CH06)/julia/convergence_rates.pdf: $(JL_CH06)/convergence_rates.jl
	@mkdir -p $(FIG_DIR_CH06)/julia
	$(JULIA) $<

# Julia figure generation rules - Chapter 7 (Chebyshev Differentiation)
$(FIG_DIR_CH07)/julia/grid_comparison.pdf: $(JL_CH07)/cheb_grid_comparison.jl $(JL_CH07)/cheb_matrix.jl
	@mkdir -p $(FIG_DIR_CH07)/julia
	$(JULIA) $<

$(FIG_DIR_CH07)/julia/cheb_matrix_structure.pdf: $(JL_CH07)/cheb_matrix_structure.jl $(JL_CH07)/cheb_matrix.jl
	@mkdir -p $(FIG_DIR_CH07)/julia
	$(JULIA) $<

$(FIG_DIR_CH07)/julia/cheb_cardinal.pdf: $(JL_CH07)/cheb_cardinal.jl $(JL_CH07)/cheb_matrix.jl
	@mkdir -p $(FIG_DIR_CH07)/julia
	$(JULIA) $<

$(FIG_DIR_CH07)/julia/cheb_diff_demo.pdf: $(JL_CH07)/cheb_diff_demo.jl $(JL_CH07)/cheb_matrix.jl
	@mkdir -p $(FIG_DIR_CH07)/julia
	$(JULIA) $<

$(FIG_DIR_CH07)/julia/convergence_waterfall.pdf: $(JL_CH07)/cheb_convergence.jl $(JL_CH07)/cheb_matrix.jl
	@mkdir -p $(FIG_DIR_CH07)/julia
	$(JULIA) $<

# Julia figure generation rules - Chapter 8 (Boundary Value Problems)
$(FIG_DIR_CH08)/julia/poisson_1d.pdf: $(JL_CH08)/bvp_linear.jl $(JL_CH07)/cheb_matrix.jl
	@mkdir -p $(FIG_DIR_CH08)/julia
	$(JULIA) $<

$(FIG_DIR_CH08)/julia/variable_coeff.pdf: $(JL_CH08)/bvp_variable_coeff.jl $(JL_CH07)/cheb_matrix.jl
	@mkdir -p $(FIG_DIR_CH08)/julia
	$(JULIA) $<

$(FIG_DIR_CH08)/julia/bratu.pdf: $(JL_CH08)/bvp_nonlinear.jl $(JL_CH07)/cheb_matrix.jl
	@mkdir -p $(FIG_DIR_CH08)/julia
	$(JULIA) $<

$(FIG_DIR_CH08)/julia/eigenvalue_problem.pdf: $(JL_CH08)/bvp_eigenvalue.jl $(JL_CH07)/cheb_matrix.jl
	@mkdir -p $(FIG_DIR_CH08)/julia
	$(JULIA) $<

$(FIG_DIR_CH08)/julia/tensor_grid.pdf: $(JL_CH08)/bvp_2d_poisson.jl $(JL_CH07)/cheb_matrix.jl
	@mkdir -p $(FIG_DIR_CH08)/julia
	$(JULIA) $<

$(FIG_DIR_CH08)/julia/poisson_2d.pdf: $(JL_CH08)/bvp_2d_poisson.jl $(JL_CH07)/cheb_matrix.jl
	@mkdir -p $(FIG_DIR_CH08)/julia
	$(JULIA) $<

$(FIG_DIR_CH08)/julia/laplacian_sparsity.pdf: $(JL_CH08)/bvp_2d_poisson.jl $(JL_CH07)/cheb_matrix.jl
	@mkdir -p $(FIG_DIR_CH08)/julia
	$(JULIA) $<

$(FIG_DIR_CH08)/julia/helmholtz.pdf: $(JL_CH08)/bvp_helmholtz.jl $(JL_CH07)/cheb_matrix.jl
	@mkdir -p $(FIG_DIR_CH08)/julia
	$(JULIA) $<

$(FIG_DIR_CH08)/julia/harmonic_oscillator.pdf: $(JL_CH08)/harmonic_oscillator.jl
	@mkdir -p $(FIG_DIR_CH08)/julia
	$(JULIA) $<

# Julia figure generation rules - Chapter 9 (Fourier Grids)
$(FIG_DIR_CH09)/julia/two_views_function.pdf: $(JL_CH09)/two_views_function.jl
	@mkdir -p $(FIG_DIR_CH09)/julia
	$(JULIA) $<

$(FIG_DIR_CH09)/julia/aliasing_demo.pdf: $(JL_CH09)/aliasing_demo.jl
	@mkdir -p $(FIG_DIR_CH09)/julia
	$(JULIA) $<

$(FIG_DIR_CH09)/julia/sinc_interpolation.pdf: $(JL_CH09)/sinc_interpolation.jl
	@mkdir -p $(FIG_DIR_CH09)/julia
	$(JULIA) $<

$(FIG_DIR_CH09)/julia/fft_aliasing.pdf: $(JL_CH09)/fft_aliasing.jl
	@mkdir -p $(FIG_DIR_CH09)/julia
	$(JULIA) $<

$(FIG_DIR_CH09)/julia/smoothness_spectra.pdf: $(JL_CH09)/smoothness_spectra.jl
	@mkdir -p $(FIG_DIR_CH09)/julia
	$(JULIA) $<

$(FIG_DIR_CH09)/julia/zero_padding_interpolation.pdf: $(JL_CH09)/zero_padding_interpolation.jl
	@mkdir -p $(FIG_DIR_CH09)/julia
	$(JULIA) $<

# Julia figure generation rules - Chapter 10 (Spectral PDE Solvers)
$(FIG_DIR_CH10)/julia/cheb_fourier_geometry.pdf: $(JL_CH10)/cheb_fourier_geometry.jl $(JL_CH10)/chebfft.jl
	@mkdir -p $(FIG_DIR_CH10)/julia
	$(JULIA) $<

$(FIG_DIR_CH10)/julia/chebfft_accuracy.pdf: $(JL_CH10)/chebfft_accuracy.jl $(JL_CH10)/chebfft.jl
	@mkdir -p $(FIG_DIR_CH10)/julia
	$(JULIA) $<

$(FIG_DIR_CH10)/julia/wave1d_waterfall.pdf: $(JL_CH10)/wave1d_cheb.jl $(JL_CH10)/chebfft.jl
	@mkdir -p $(FIG_DIR_CH10)/julia
	$(JULIA) $<

$(FIG_DIR_CH10)/julia/wave2d_snapshots.pdf: $(JL_CH10)/wave2d_cheb.jl $(JL_CH10)/chebfft.jl
	@mkdir -p $(FIG_DIR_CH10)/julia
	$(JULIA) $<

$(FIG_DIR_CH10)/julia/heat1d_evolution.pdf: $(JL_CH10)/heat1d_cheb.jl $(JL_CH10)/chebfft.jl
	@mkdir -p $(FIG_DIR_CH10)/julia
	$(JULIA) $<

$(FIG_DIR_CH10)/julia/heat2d_snapshots.pdf: $(JL_CH10)/heat2d_cheb.jl $(JL_CH10)/chebfft.jl
	@mkdir -p $(FIG_DIR_CH10)/julia
	$(JULIA) $<

$(FIG_DIR_CH10)/julia/poisson2d_solution.pdf: $(JL_CH10)/poisson2d_cheb.jl $(JL_CH10)/chebfft.jl
	@mkdir -p $(FIG_DIR_CH10)/julia
	$(JULIA) $<

$(FIG_DIR_CH10)/julia/transport_variable.pdf: $(JL_CH10)/transport_variable.jl
	@mkdir -p $(FIG_DIR_CH10)/julia
	$(JULIA) $<

# Teaching plan compilation
tplan: $(TPLAN_OUT)

$(TPLAN_OUT): $(TPLAN_SRC)
	mkdir -p $(TPLAN_OUT_DIR)
	$(TYPST) compile $(TPLAN_SRC) $(TPLAN_OUT)

# Clean targets
clean:
	rm -f $(OUT)

clean-tplan:
	rm -f $(TPLAN_OUT)

clean-figures:
	rm -f $(PY_FIGS) $(M_FIGS) $(JL_FIGS)
	rm -f $(FIG_DIR_CH02)/python/*.png $(FIG_DIR_CH02)/matlab/*.png $(FIG_DIR_CH02)/julia/*.png
	rm -f $(FIG_DIR_CH03)/python/*.png $(FIG_DIR_CH03)/matlab/*.png $(FIG_DIR_CH03)/julia/*.png
	rm -f $(FIG_DIR_CH04)/python/*.png $(FIG_DIR_CH04)/matlab/*.png $(FIG_DIR_CH04)/julia/*.png
	rm -f $(FIG_DIR_CH05)/python/*.png $(FIG_DIR_CH05)/matlab/*.png $(FIG_DIR_CH05)/julia/*.png
	rm -f $(FIG_DIR_CH06)/python/*.png $(FIG_DIR_CH06)/matlab/*.png $(FIG_DIR_CH06)/julia/*.png
	rm -f $(FIG_DIR_CH07)/python/*.png $(FIG_DIR_CH07)/matlab/*.png $(FIG_DIR_CH07)/julia/*.png
	rm -f $(FIG_DIR_CH08)/python/*.png $(FIG_DIR_CH08)/matlab/*.png $(FIG_DIR_CH08)/julia/*.png
	rm -f $(FIG_DIR_CH09)/python/*.png $(FIG_DIR_CH09)/matlab/*.png $(FIG_DIR_CH09)/julia/*.png
	rm -f $(FIG_DIR_CH10)/python/*.png $(FIG_DIR_CH10)/matlab/*.png $(FIG_DIR_CH10)/julia/*.png

clean-all: clean clean-tplan clean-figures
