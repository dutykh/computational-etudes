.PHONY: book clean

TYPST ?= typst
SRC = book/main.typ
OUT_DIR = book/build
OUT = $(OUT_DIR)/DD-Computational-Etudes.pdf

book: $(OUT)

$(OUT): $(SRC) book/chapters/preface.typ book/chapters/introduction.typ book/styles/template.typ
	mkdir -p $(OUT_DIR)
	$(TYPST) compile $(SRC) $(OUT)

clean:
	rm -f $(OUT)
