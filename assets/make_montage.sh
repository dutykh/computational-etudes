#!/usr/bin/env bash
# Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
#
# The strip of three figures of the book shown at the top of README.md, made from the figures that the
# Python codes write into textbook/figures (chapters 4 and 14). Needs ImageMagick. Run from anywhere:
#   bash assets/make_montage.sh
set -euo pipefail
cd "$(dirname "$0")/.."
F=textbook/figures
montage "$F/ch04/python/chebyshev_points_circle.png" \
        "$F/ch14/python/kuramoto_sivashinsky.png" \
        "$F/ch14/python/orr_sommerfeld.png" \
        -tile 3x1 -geometry x900+24+12 -background white miff:- |
    convert - -resize 1800x -strip -quality 88 assets/etudes-figures.jpg
