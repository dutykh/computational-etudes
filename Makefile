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

# Combined figure variables
PY_FIGS = $(PY_FIGS_CH02) $(PY_FIGS_CH03)
M_FIGS = $(M_FIGS_CH02) $(M_FIGS_CH03)

# Default target: build everything
all: figures textbook tplan

# Build textbook (depends on figures)
textbook: $(OUT)

$(OUT): $(SRC) textbook/chapters/preface.typ textbook/chapters/introduction.typ textbook/chapters/classical_pdes.typ textbook/chapters/mise_en_bouche.typ textbook/styles/template.typ $(PY_FIGS)
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

clean-all: clean clean-tplan clean-figures
