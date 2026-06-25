# ==========================================
# Makefile for the Response to the Reviewers' Reports
# Manuscript: "Computational Études: A Spectral Approach"
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology,
#         Abu Dhabi, UAE)
# ==========================================

TYPST := typst
DOC   := ResponseToReviewers
SRC   := $(DOC).typ
PDF   := DD-Computational-Etudes-Response-to-Reviewers.pdf

ASSETS := assets/KUlogo.png assets/CCMSLogo.png assets/SignDD.png

.PHONY: all clean

all: $(PDF)

$(PDF): $(SRC) $(ASSETS)
	$(TYPST) compile $(SRC) $(PDF)

clean:
	rm -f $(PDF)
