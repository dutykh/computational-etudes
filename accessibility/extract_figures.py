#!/usr/bin/env python3
"""
Extract every #figure(...) block from the Computational Etudes Typst sources and
assign book figure/table numbers + Springer-style ids.

Produces accessibility/figures.json — one record per figure block, in document
order, with: chapter, kind (image|table|drawn), book number ("Fig. 3.1" /
"Table 8.2"), springer id (fig3_1 / tbl8_2), source image path(s), label, and a
lightly-cleaned plain-text caption.

The parser is mode-aware: it tracks Typst code vs. markup contexts, strings,
math ($...$), and comments, because:
  * inside a caption [...] (markup) the chars ( ) are *literal text*, only [ ] { }
    are structural;
  * inside image("..") (code) " starts a string;
  * inside $...$ math, brackets/parens are literal symbols.
Naive bracket counting would mis-match on captions like "$x in [-1, 1)$" or
"(left panel)".
"""
import json
import os
import re

BASE = "/home/user/workspace/computational-etudes/textbook"
CHAPTERS = os.path.join(BASE, "chapters")
OUT = os.path.join(os.path.dirname(__file__), "figures.json")

# (filename, chapter number) in book order. Introduction = Ch 1.
NUMBERED = [
    ("introduction.typ", 1),
    ("classical_pdes.typ", 2),
    ("mise_en_bouche.typ", 3),
    ("geometry_of_nodes.typ", 4),
    ("differentiation_matrices.typ", 5),
    ("smoothness_accuracy.typ", 6),
    ("chebyshev_differentiation.typ", 7),
    ("boundary_value_problems.typ", 8),
    ("fourier_grids.typ", 9),
    ("spectral_pde_solvers.typ", 10),
    ("fourier_pseudospectral.typ", 11),
    ("polar_coordinates.typ", 12),
    ("advanced_boundary_conditions.typ", 13),
    ("higher_order_bvps.typ", 14),
    ("quadrature.typ", 15),
    ("periodic_quadrature.typ", 16),
    ("time_stability.typ", 17),
    ("linear_eigenproblems.typ", 18),
    ("coordinate_transformations.typ", 19),
    ("unbounded_intervals.typ", 20),
    ("special_tricks.typ", 21),
]

OPEN = {"(": ")", "[": "]", "{": "}"}


def find_match(text, i):
    """text[i] is an opening bracket; return index of its matching close char.

    Interior context is 'markup' for '[' and 'code' for '(' / '{'.
    """
    open_ch = text[i]
    close_ch = OPEN[open_ch]
    interior = "markup" if open_ch == "[" else "code"
    n = len(text)
    j = i + 1
    while j < n:
        c = text[j]
        # comments (both contexts)
        if c == "/" and j + 1 < n and text[j + 1] == "/":
            k = text.find("\n", j)
            j = n if k == -1 else k + 1
            continue
        if c == "/" and j + 1 < n and text[j + 1] == "*":
            k = text.find("*/", j + 2)
            j = n if k == -1 else k + 2
            continue
        # math mode $...$ (opaque)
        if c == "$":
            j += 1
            while j < n:
                if text[j] == "\\":
                    j += 2
                    continue
                if text[j] == "$":
                    j += 1
                    break
                j += 1
            continue
        if interior == "code":
            if c == '"':
                j += 1
                while j < n:
                    if text[j] == "\\":
                        j += 2
                        continue
                    if text[j] == '"':
                        j += 1
                        break
                    j += 1
                continue
            if c == close_ch:
                return j
            if c in OPEN:
                m = find_match(text, j)
                if m == -1:
                    return -1
                j = m + 1
                continue
            j += 1
            continue
        else:  # markup: only [ ] { } structural; ( ) are literal text
            if c == "\\":
                j += 2
                continue
            if c == close_ch:
                return j
            if c == "[" or c == "{":
                m = find_match(text, j)
                if m == -1:
                    return -1
                j = m + 1
                continue
            j += 1
            continue
    return -1


def clean_caption(raw):
    """Light cleanup of Typst caption markup into readable plain text."""
    s = raw
    # strip index macros: #idx("term") -> term
    s = re.sub(r'#idx\(\s*"([^"]*)"\s*\)', r"\1", s)
    # unescape common backslash escapes
    for esc, rep in [("\\/", "/"), ("\\#", "#"), ("\\_", "_"),
                     ("\\&", "&"), ("\\@", "@"), ("\\$", "$"),
                     ("\\*", "*"), ("\\~", "~")]:
        s = s.replace(esc, rep)
    # collapse whitespace/newlines
    s = re.sub(r"\s+", " ", s).strip()
    return s


def extract_caption(body):
    m = re.search(r"caption:\s*", body)
    if not m:
        return ""
    # find first '[' after caption:
    k = body.find("[", m.end())
    if k == -1:
        return ""
    end = find_match(body, k)
    if end == -1:
        return ""
    return clean_caption(body[k + 1:end])


def extract_figures_from_file(path, chapter):
    with open(path, encoding="utf-8") as f:
        text = f.read()
    figs = []
    for m in re.finditer(r"#figure\(", text):
        paren = m.end() - 1  # index of '('
        close = find_match(text, paren)
        if close == -1:
            raise RuntimeError(f"unmatched #figure( in {path} at offset {paren}")
        body = text[paren + 1:close]
        line = text[:m.start()].count("\n") + 1
        # label after the closing paren
        tail = text[close + 1:close + 80]
        lm = re.match(r"\s*<([^>]+)>", tail)
        label = lm.group(1) if lm else None
        # image paths
        images = re.findall(r'image\(\s*"([^"]+)"', body)
        # kind
        if images:
            kind = "image"
        elif re.search(r"\btable\(", body):
            kind = "table"
        else:
            kind = "drawn"
        figs.append({
            "chapter": chapter,
            "file": os.path.basename(path),
            "line": line,
            "kind": kind,
            "label": label,
            "images": images,
            "caption": extract_caption(body),
        })
    return figs


def main():
    records = []
    for fname, ch in NUMBERED:
        path = os.path.join(CHAPTERS, fname)
        figs = extract_figures_from_file(path, ch)
        fig_n = 0
        tbl_n = 0
        for fr in figs:
            if fr["kind"] == "table":
                tbl_n += 1
                fr["number"] = f"Table {ch}.{tbl_n}"
                fr["springer_id"] = f"tbl{ch}_{tbl_n}"
            else:  # image or drawn -> "Figure" counter
                fig_n += 1
                fr["number"] = f"Fig. {ch}.{fig_n}"
                fr["springer_id"] = f"fig{ch}_{fig_n}"
            records.append(fr)

    # Afterword (back matter) — has 1 figure; number it separately.
    aw = os.path.join(CHAPTERS, "afterword.typ")
    aw_figs = extract_figures_from_file(aw, "Afterword")
    for idx, fr in enumerate(aw_figs, 1):
        if fr["kind"] == "table":
            fr["number"] = f"Afterword Table {idx}"
            fr["springer_id"] = f"figAW_tbl{idx}"
        else:
            fr["number"] = f"Afterword Fig. {idx}"
            fr["springer_id"] = f"figAW_{idx}"
        records.append(fr)

    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(records, f, ensure_ascii=False, indent=2)

    # summary
    from collections import Counter
    by_kind = Counter(r["kind"] for r in records)
    print(f"wrote {len(records)} records -> {OUT}")
    print("by kind:", dict(by_kind))
    per_ch = Counter(r["chapter"] for r in records)
    for ch in list(range(1, 22)) + ["Afterword"]:
        if per_ch.get(ch):
            print(f"  ch {ch}: {per_ch[ch]}")


if __name__ == "__main__":
    main()
