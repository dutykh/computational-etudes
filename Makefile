.PHONY: book clean figures figures-python figures-matlab all

# Tools
TYPST ?= typst
PYTHON ?= python3
MATLAB ?= matlab

# Book compilation
SRC = book/main.typ
OUT_DIR = book/build
OUT = $(OUT_DIR)/DD-Computational-Etudes.pdf

# Python scripts
PY_CH02 = codes/python/ch02_classical_pdes
PY_SCRIPTS = $(PY_CH02)/heat_equation_evolution.py \
             $(PY_CH02)/heat_equation_waterfall.py \
             $(PY_CH02)/wave_equation_evolution.py \
             $(PY_CH02)/wave_equation_waterfall.py \
             $(PY_CH02)/laplace_equation_2d.py

# MATLAB scripts
M_CH02 = codes/matlab/ch02_classical_pdes
M_SCRIPTS = $(M_CH02)/heat_equation_evolution.m \
            $(M_CH02)/heat_equation_waterfall.m \
            $(M_CH02)/wave_equation_evolution.m \
            $(M_CH02)/wave_equation_waterfall.m \
            $(M_CH02)/laplace_equation_2d.m

# Figure outputs
FIG_DIR = book/figures/ch02
PY_FIGS = $(FIG_DIR)/python/heat_evolution.pdf \
          $(FIG_DIR)/python/heat_waterfall.pdf \
          $(FIG_DIR)/python/wave_evolution.pdf \
          $(FIG_DIR)/python/wave_waterfall.pdf \
          $(FIG_DIR)/python/laplace_solution.pdf
M_FIGS = $(FIG_DIR)/matlab/heat_evolution.pdf \
         $(FIG_DIR)/matlab/heat_waterfall.pdf \
         $(FIG_DIR)/matlab/wave_evolution.pdf \
         $(FIG_DIR)/matlab/wave_waterfall.pdf \
         $(FIG_DIR)/matlab/laplace_solution.pdf

# Default target: build everything
all: figures book

# Build book (depends on figures)
book: $(OUT)

$(OUT): $(SRC) book/chapters/preface.typ book/chapters/introduction.typ book/chapters/classical_pdes.typ book/styles/template.typ $(PY_FIGS)
	mkdir -p $(OUT_DIR)
	$(TYPST) compile $(SRC) $(OUT)

# Generate all figures (Python is primary, MATLAB is optional)
figures: figures-python

figures-python: $(PY_FIGS)

figures-matlab: $(M_FIGS)

# Python figure generation rules
$(FIG_DIR)/python/heat_evolution.pdf: $(PY_CH02)/heat_equation_evolution.py
	@mkdir -p $(FIG_DIR)/python
	$(PYTHON) $<

$(FIG_DIR)/python/wave_evolution.pdf: $(PY_CH02)/wave_equation_evolution.py
	@mkdir -p $(FIG_DIR)/python
	$(PYTHON) $<

$(FIG_DIR)/python/laplace_solution.pdf: $(PY_CH02)/laplace_equation_2d.py
	@mkdir -p $(FIG_DIR)/python
	$(PYTHON) $<

# Python waterfall figure rules
$(FIG_DIR)/python/heat_waterfall.pdf: $(PY_CH02)/heat_equation_waterfall.py
	@mkdir -p $(FIG_DIR)/python
	$(PYTHON) $<

$(FIG_DIR)/python/wave_waterfall.pdf: $(PY_CH02)/wave_equation_waterfall.py
	@mkdir -p $(FIG_DIR)/python
	$(PYTHON) $<

# MATLAB figure generation rules
$(FIG_DIR)/matlab/heat_evolution.pdf: $(M_CH02)/heat_equation_evolution.m
	@mkdir -p $(FIG_DIR)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR)/matlab/heat_waterfall.pdf: $(M_CH02)/heat_equation_waterfall.m
	@mkdir -p $(FIG_DIR)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR)/matlab/wave_evolution.pdf: $(M_CH02)/wave_equation_evolution.m
	@mkdir -p $(FIG_DIR)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR)/matlab/wave_waterfall.pdf: $(M_CH02)/wave_equation_waterfall.m
	@mkdir -p $(FIG_DIR)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

$(FIG_DIR)/matlab/laplace_solution.pdf: $(M_CH02)/laplace_equation_2d.m
	@mkdir -p $(FIG_DIR)/matlab
	$(MATLAB) -nodisplay -nosplash -batch "run('$<')"

# Clean targets
clean:
	rm -f $(OUT)

clean-figures:
	rm -f $(PY_FIGS) $(M_FIGS)
	rm -f $(FIG_DIR)/python/*.png $(FIG_DIR)/matlab/*.png

clean-all: clean clean-figures
