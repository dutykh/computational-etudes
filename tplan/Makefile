# tplan/Makefile
# Makefile for Teaching Plan compilation
#
# Author: Dr. Denys Dutykh
#         Mathematics Department
#         Khalifa University of Science and Technology
#         Abu Dhabi, UAE

.PHONY: all clean

# Tools
TYPST ?= typst

# Source and output
SRC = teaching_plan.typ
OUT_DIR = build
OUT = $(OUT_DIR)/teaching_plan.pdf

# Default target
all: $(OUT)

# Compile teaching plan
$(OUT): $(SRC)
	@mkdir -p $(OUT_DIR)
	$(TYPST) compile $(SRC) $(OUT)
	@echo "Teaching plan compiled: $(OUT)"

# Clean target
clean:
	rm -f $(OUT)
	@echo "Cleaned teaching plan build files"
