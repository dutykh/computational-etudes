.PHONY: textbook clean figures figures-python figures-matlab all tplan

# Tools
TYPST ?= typst
PYTHON ?= python3
MATLAB ?= matlab

# Textbook compilation
SRC = textbook/main.typ
OUT_DIR = textbook/build
OUT = $(OUT_DIR)/DD-Computational-Etudes-2026.pdf

# Teaching plan compilation
TPLAN_SRC = tplan/teaching_plan.typ
TPLAN_OUT_DIR = tplan/build
TPLAN_OUT = $(TPLAN_OUT_DIR)/teaching_plan.pdf

# Python scripts - Chapter 2
PY_CH02 = codes/python/ch02_classical_pdes
PY_SCRIPTS_CH02 = $(PY_CH02)/heat_equation_evolution.py \
                  $(PY_CH02)/heat_equation_waterfall.py \
                  $(PY_CH02)/wave_equation_evolution.py \
                  $(PY_CH02)/wave_equation_waterfall.py \
                  $(PY_CH02)/laplace_equation_2d.py

# Python scripts - Chapter 3
PY_CH03 = codes/python/ch03_mise_en_bouche
PY_SCRIPTS_CH03 = $(PY_CH03)/collocation_example1.py \
                  $(PY_CH03)/collocation_vs_galerkin.py

# Python scripts - Chapter 4
PY_CH04 = codes/python/ch04_geometry_of_nodes
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
PY_CH05 = codes/python/ch05_differentiation_matrices
PY_SCRIPTS_CH05 = $(PY_CH05)/fd_matrix_bandwidth.py \
                  $(PY_CH05)/spectral_matrix_structure.py \
                  $(PY_CH05)/fd_stencil_schematic.py \
                  $(PY_CH05)/stencil_pyramid.py \
                  $(PY_CH05)/convergence_comparison.py \
                  $(PY_CH05)/spectral_derivatives_demo.py

# Python scripts - Chapter 6
PY_CH06 = codes/python/ch06_chebyshev_differentiation
PY_SCRIPTS_CH06 = $(PY_CH06)/cheb_matrix.py \
                  $(PY_CH06)/cheb_grid_comparison.py \
                  $(PY_CH06)/cheb_matrix_structure.py \
                  $(PY_CH06)/cheb_cardinal.py \
                  $(PY_CH06)/cheb_diff_demo.py \
                  $(PY_CH06)/cheb_convergence.py

# Python scripts - Chapter 7
PY_CH07 = codes/python/ch07_boundary_value_problems
PY_SCRIPTS_CH07 = $(PY_CH07)/bvp_linear.py \
                  $(PY_CH07)/bvp_variable_coeff.py \
                  $(PY_CH07)/bvp_nonlinear.py \
                  $(PY_CH07)/bvp_eigenvalue.py \
                  $(PY_CH07)/bvp_2d_poisson.py \
                  $(PY_CH07)/bvp_helmholtz.py

# MATLAB scripts - Chapter 2
M_CH02 = codes/matlab/ch02_classical_pdes
M_SCRIPTS_CH02 = $(M_CH02)/heat_equation_evolution.m \
                 $(M_CH02)/heat_equation_waterfall.m \
                 $(M_CH02)/wave_equation_evolution.m \
                 $(M_CH02)/wave_equation_waterfall.m \
                 $(M_CH02)/laplace_equation_2d.m

# MATLAB scripts - Chapter 3
M_CH03 = codes/matlab/ch03_mise_en_bouche
M_SCRIPTS_CH03 = $(M_CH03)/collocation_example1.m \
                 $(M_CH03)/collocation_vs_galerkin.m

# MATLAB scripts - Chapter 4
M_CH04 = codes/matlab/ch04_geometry_of_nodes
M_SCRIPTS_CH04 = $(M_CH04)/runge_phenomenon.m \
                 $(M_CH04)/chebyshev_success.m \
                 $(M_CH04)/chebyshev_points_circle.m \
                 $(M_CH04)/equipotential_curves.m \
                 $(M_CH04)/lagrange_basis.m \
                 $(M_CH04)/lebesgue_functions.m \
                 $(M_CH04)/lebesgue_random_nodes.m \
                 $(M_CH04)/convergence_comparison.m

# MATLAB scripts - Chapter 5
M_CH05 = codes/matlab/ch05_differentiation_matrices
M_SCRIPTS_CH05 = $(M_CH05)/fd_matrix_bandwidth.m \
                 $(M_CH05)/spectral_matrix_structure.m \
                 $(M_CH05)/fd_stencil_schematic.m \
                 $(M_CH05)/stencil_pyramid.m \
                 $(M_CH05)/convergence_comparison.m \
                 $(M_CH05)/spectral_derivatives_demo.m

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

# Figure outputs - Chapter 6
FIG_DIR_CH06 = textbook/figures/ch06
PY_FIGS_CH06 = $(FIG_DIR_CH06)/python/grid_comparison.pdf \
               $(FIG_DIR_CH06)/python/cheb_matrix_structure.pdf \
               $(FIG_DIR_CH06)/python/cheb_cardinal.pdf \
               $(FIG_DIR_CH06)/python/cheb_diff_demo.pdf \
               $(FIG_DIR_CH06)/python/convergence_waterfall.pdf

# Figure outputs - Chapter 7
FIG_DIR_CH07 = textbook/figures/ch07
PY_FIGS_CH07 = $(FIG_DIR_CH07)/python/poisson_1d.pdf \
               $(FIG_DIR_CH07)/python/variable_coeff.pdf \
               $(FIG_DIR_CH07)/python/bratu.pdf \
               $(FIG_DIR_CH07)/python/eigenvalue_problem.pdf \
               $(FIG_DIR_CH07)/python/tensor_grid.pdf \
               $(FIG_DIR_CH07)/python/poisson_2d.pdf \
               $(FIG_DIR_CH07)/python/laplacian_sparsity.pdf \
               $(FIG_DIR_CH07)/python/helmholtz.pdf

# Combined figure variables
PY_FIGS = $(PY_FIGS_CH02) $(PY_FIGS_CH03) $(PY_FIGS_CH04) $(PY_FIGS_CH05) $(PY_FIGS_CH06) $(PY_FIGS_CH07)
M_FIGS = $(M_FIGS_CH02) $(M_FIGS_CH03) $(M_FIGS_CH04) $(M_FIGS_CH05)

# Default target: build everything
all: figures textbook tplan

# Build textbook (depends on figures)
textbook: $(OUT)

$(OUT): $(SRC) textbook/chapters/preface.typ textbook/chapters/introduction.typ textbook/chapters/classical_pdes.typ textbook/chapters/mise_en_bouche.typ textbook/chapters/geometry_of_nodes.typ textbook/chapters/differentiation_matrices.typ textbook/chapters/chebyshev_differentiation.typ textbook/chapters/boundary_value_problems.typ textbook/styles/template.typ $(PY_FIGS)
	mkdir -p $(OUT_DIR)
	$(TYPST) compile $(SRC) $(OUT)

# Generate all figures (Python is primary, MATLAB is optional)
figures: figures-python

figures-python: $(PY_FIGS)

figures-matlab: $(M_FIGS)

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

# Python figure generation rules - Chapter 6
$(FIG_DIR_CH06)/python/grid_comparison.pdf: $(PY_CH06)/cheb_grid_comparison.py $(PY_CH06)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH06)/python
	$(PYTHON) $<

$(FIG_DIR_CH06)/python/cheb_matrix_structure.pdf: $(PY_CH06)/cheb_matrix_structure.py $(PY_CH06)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH06)/python
	$(PYTHON) $<

$(FIG_DIR_CH06)/python/cheb_cardinal.pdf: $(PY_CH06)/cheb_cardinal.py $(PY_CH06)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH06)/python
	$(PYTHON) $<

$(FIG_DIR_CH06)/python/cheb_diff_demo.pdf: $(PY_CH06)/cheb_diff_demo.py $(PY_CH06)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH06)/python
	$(PYTHON) $<

$(FIG_DIR_CH06)/python/convergence_waterfall.pdf: $(PY_CH06)/cheb_convergence.py $(PY_CH06)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH06)/python
	$(PYTHON) $<

# Python figure generation rules - Chapter 7
$(FIG_DIR_CH07)/python/poisson_1d.pdf: $(PY_CH07)/bvp_linear.py $(PY_CH06)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH07)/python
	$(PYTHON) $<

$(FIG_DIR_CH07)/python/variable_coeff.pdf: $(PY_CH07)/bvp_variable_coeff.py $(PY_CH06)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH07)/python
	$(PYTHON) $<

$(FIG_DIR_CH07)/python/bratu.pdf: $(PY_CH07)/bvp_nonlinear.py $(PY_CH06)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH07)/python
	$(PYTHON) $<

$(FIG_DIR_CH07)/python/eigenvalue_problem.pdf: $(PY_CH07)/bvp_eigenvalue.py $(PY_CH06)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH07)/python
	$(PYTHON) $<

$(FIG_DIR_CH07)/python/tensor_grid.pdf: $(PY_CH07)/bvp_2d_poisson.py $(PY_CH06)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH07)/python
	$(PYTHON) $<

$(FIG_DIR_CH07)/python/poisson_2d.pdf: $(PY_CH07)/bvp_2d_poisson.py $(PY_CH06)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH07)/python
	$(PYTHON) $<

$(FIG_DIR_CH07)/python/laplacian_sparsity.pdf: $(PY_CH07)/bvp_2d_poisson.py $(PY_CH06)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH07)/python
	$(PYTHON) $<

$(FIG_DIR_CH07)/python/helmholtz.pdf: $(PY_CH07)/bvp_helmholtz.py $(PY_CH06)/cheb_matrix.py
	@mkdir -p $(FIG_DIR_CH07)/python
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
	rm -f $(PY_FIGS) $(M_FIGS)
	rm -f $(FIG_DIR_CH02)/python/*.png $(FIG_DIR_CH02)/matlab/*.png
	rm -f $(FIG_DIR_CH03)/python/*.png $(FIG_DIR_CH03)/matlab/*.png
	rm -f $(FIG_DIR_CH04)/python/*.png $(FIG_DIR_CH04)/matlab/*.png
	rm -f $(FIG_DIR_CH05)/python/*.png $(FIG_DIR_CH05)/matlab/*.png
	rm -f $(FIG_DIR_CH06)/python/*.png
	rm -f $(FIG_DIR_CH07)/python/*.png

clean-all: clean clean-tplan clean-figures
