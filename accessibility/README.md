# Alt Text for Springer submission

This folder contains the **alternative text (Alt Text)** for every figure and
table in *Computational Études: A Spectral Approach*, prepared for Springer
Nature's accessibility (EU Accessibility Act / WCAG) requirement.

## The deliverable

**`alt_text.xlsx`** — submit this file to Springer together with the manuscript.
`alt_text.csv` is an identical UTF‑8 copy (opens natively in Excel) as a backup.

One row per figure/table — **222 rows** total:

| Count | Kind |
|------:|------|
| 180 | image figures (`figN_M`) |
| 3 | drawn vector illustrations — Fig. 3.3, 3.4, 9.2 (no separate image file) |
| 38 | numbered tables (`tblN_M`) |
| 1 | afterword comparison table |

Figures and tables use **separate counters**, matching the book: image/drawn
figures are "Fig. N.M", tables are "Table N.M".

### Columns
`Chapter` · `Type` (Figure/Table) · `Number` (e.g. *Fig. 8.3*) ·
`Springer ID` (e.g. *fig8_3* / *tbl8_2*) · `Source file` (manuscript-relative
path, or a note for drawn/table entries) · `Caption` · `Alt Text`.

This is a superset — hide/remove any columns Springer's template does not need.
If Springer provided an Excel template, the columns can be re-mapped to it.

## How it was produced

1. `extract_figures.py` — parses every `#figure(...)` block in
   `textbook/chapters/*.typ` with a Typst-aware bracket matcher, classifies kind
   (image / table / drawn), assigns book numbers + Springer ids, and writes
   `figures.json`. Counts were cross-checked against `typst query` (183 image +
   39 table = 222).
2. `make_payloads.py` — splits `figures.json` into per-chapter files under
   `payloads/`.
3. Alt Text was written **per chapter** from the captions + surrounding prose
   (caption-grounded; images were not rendered), saved under `alt_out/chN.json`.
   Style follows Springer's "How to Write Good Alt Text" / WCAG: describes the
   visual content and takeaway (not a copy of the caption), no "Image of…"
   openers, math spelled out in words, plain text.
4. `build_xlsx.py` — merges everything into `alt_text.xlsx` / `alt_text.csv`,
   and verifies 100 % coverage (every figure has Alt Text, no orphan entries).

## Regenerating after figure changes

```bash
python3 accessibility/extract_figures.py     # re-scan chapters
python3 accessibility/make_payloads.py        # refresh per-chapter payloads
# (re-write Alt Text only for changed chapters in accessibility/alt_out/)
python3 accessibility/build_xlsx.py           # rebuild the spreadsheet
```

The build prints a warning if any figure lacks Alt Text or any Alt Text key has
no matching figure, so drift is caught automatically.
