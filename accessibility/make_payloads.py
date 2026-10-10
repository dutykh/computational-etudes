#!/usr/bin/env python3
"""Split figures.json into per-chapter payload files for the alt-text agents."""
import json
import os
from collections import defaultdict

HERE = os.path.dirname(__file__)
recs = json.load(open(os.path.join(HERE, "figures.json"), encoding="utf-8"))

os.makedirs(os.path.join(HERE, "payloads"), exist_ok=True)
os.makedirs(os.path.join(HERE, "alt_out"), exist_ok=True)

by_ch = defaultdict(list)
for r in recs:
    by_ch[r["chapter"]].append(r)


def slug(ch):
    return "AW" if ch == "Afterword" else str(ch)


for ch, figs in by_ch.items():
    payload = {
        "chapter": ch,
        "file": figs[0]["file"],
        "figures": [
            {
                "springer_id": f["springer_id"],
                "number": f["number"],
                "kind": f["kind"],
                "label": f["label"],
                "primary_image": (f["images"][0] if f["images"] else None),
                "n_images": len(f["images"]),
                "caption": f["caption"],
            }
            for f in figs
        ],
    }
    p = os.path.join(HERE, "payloads", f"ch{slug(ch)}.json")
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)
    print(f"ch{slug(ch)}: {len(figs)} figures -> {p}")
