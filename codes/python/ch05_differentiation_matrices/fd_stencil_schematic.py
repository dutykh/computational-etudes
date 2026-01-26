#!/usr/bin/env python3
"""
fd_stencil_schematic.py

Creates a schematic illustration of finite difference stencils in 1D,
showing how the stencil width increases with the order of accuracy,
culminating in the spectral method that uses all nodes.

Illustrates:
    - 3-point stencil (2nd order): uses x_{i-1}, x_i, x_{i+1}
    - 5-point stencil (4th order): uses x_{i-2}, ..., x_{i+2}
    - 7-point stencil (6th order): uses x_{i-3}, ..., x_{i+3}
    - Spectral stencil: uses ALL nodes

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Études: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Rectangle
from pathlib import Path

# -----------------------------------------------------------------------------
# Publication-quality matplotlib configuration
# -----------------------------------------------------------------------------
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman', 'CMU Serif', 'DejaVu Serif']
rcParams['mathtext.fontset'] = 'cm'
rcParams['axes.labelsize'] = 11
rcParams['axes.titlesize'] = 12
rcParams['xtick.labelsize'] = 10
rcParams['ytick.labelsize'] = 10
rcParams['legend.fontsize'] = 10
rcParams['figure.dpi'] = 150
rcParams['savefig.dpi'] = 300
rcParams['text.usetex'] = False
rcParams['axes.linewidth'] = 0.8

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
# Book color scheme
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'
LIGHT_GRAY = '#BBBBBB'
VERY_LIGHT = '#DDDDDD'

# Output path
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch05' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'fd_stencil_schematic.pdf'


# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
def main():
    # Grid setup
    n_points = 11  # Total points to show (odd number for symmetry)
    center_idx = 5  # Index of central point (0-indexed)

    # Stencil configurations: (half_width, label, is_spectral)
    # half_width = -1 means spectral (all points)
    stencils = [
        (1, '2nd order FD (3 points)', False),
        (2, '4th order FD (5 points)', False),
        (3, '6th order FD (7 points)', False),
        (-1, 'Spectral (all points)', True),
    ]

    # Create figure
    fig, axes = plt.subplots(4, 1, figsize=(11, 7.5))

    # Point positions
    x_left = 0.08
    x_right = 0.72
    x_positions = np.linspace(x_left, x_right, n_points)
    h = x_positions[1] - x_positions[0]

    for ax_idx, (ax, (half_width, label, is_spectral)) in enumerate(zip(axes, stencils)):
        ax.set_xlim(0, 1)
        ax.set_ylim(-1, 1)
        ax.axis('off')

        y_base = 0.0

        # Draw horizontal grid line
        ax.axhline(y=y_base, xmin=0.05, xmax=0.77, color=VERY_LIGHT,
                   linewidth=1.2, linestyle='--', zorder=0)

        # For spectral: all points are in the stencil
        if is_spectral:
            active_indices = set(range(n_points))
        else:
            active_indices = set(range(center_idx - half_width,
                                       center_idx + half_width + 1))

        # Draw shaded region for stencil
        if is_spectral:
            rect_left = x_positions[0] - 0.012
            rect_right = x_positions[-1] + 0.012
        else:
            rect_left = x_positions[center_idx - half_width] - 0.012
            rect_right = x_positions[center_idx + half_width] + 0.012

        rect = Rectangle((rect_left, y_base - 0.45),
                          rect_right - rect_left, 0.9,
                          facecolor=SKY, alpha=0.15, edgecolor=SKY,
                          linewidth=1.2, zorder=1)
        ax.add_patch(rect)

        # Draw all grid points
        for i, x in enumerate(x_positions):
            if i == center_idx:
                # Central point
                color = CORAL
                size = 220
                zorder = 10
            elif i in active_indices:
                # Active stencil point
                color = NAVY
                size = 160
                zorder = 5
            else:
                # Inactive point
                color = LIGHT_GRAY
                size = 80
                zorder = 3

            ax.scatter(x, y_base, s=size, c=color, zorder=zorder,
                       edgecolors='white', linewidths=2)

        # Add continuation dots
        ax.text(x_positions[0] - 0.025, y_base, r'$\cdots$',
                ha='right', va='center', fontsize=16, color=LIGHT_GRAY)
        ax.text(x_positions[-1] + 0.025, y_base, r'$\cdots$',
                ha='left', va='center', fontsize=16, color=LIGHT_GRAY)

        # Add labels BELOW the points to avoid overlap
        label_y = y_base - 0.7

        if is_spectral:
            # For spectral: label first, center, and last
            labels_spec = [
                (0, r'$x_0$', NAVY),
                (center_idx, r'$x_i$', CORAL),
                (n_points - 1, r'$x_N$', NAVY),
            ]
            for idx, lbl, col in labels_spec:
                ax.text(x_positions[idx], label_y, lbl,
                        ha='center', va='top', fontsize=12, color=col)
        else:
            # For FD: label boundary and center of stencil
            left_idx = center_idx - half_width
            right_idx = center_idx + half_width

            if half_width == 1:
                left_lbl, right_lbl = r'$x_{i-1}$', r'$x_{i+1}$'
            elif half_width == 2:
                left_lbl, right_lbl = r'$x_{i-2}$', r'$x_{i+2}$'
            else:
                left_lbl, right_lbl = r'$x_{i-3}$', r'$x_{i+3}$'

            ax.text(x_positions[left_idx], label_y, left_lbl,
                    ha='center', va='top', fontsize=12, color=NAVY)
            ax.text(x_positions[center_idx], label_y, r'$x_i$',
                    ha='center', va='top', fontsize=12, color=CORAL)
            ax.text(x_positions[right_idx], label_y, right_lbl,
                    ha='center', va='top', fontsize=12, color=NAVY)

        # Add method label on the right
        ax.text(0.88, y_base, label, ha='left', va='center',
                fontsize=12, color=NAVY, fontweight='bold')

        # Draw bracket below stencil
        bracket_y = y_base + 0.55
        tick_len = 0.12
        ax.plot([rect_left, rect_left], [bracket_y - tick_len, bracket_y],
                color=NAVY, linewidth=1.5, solid_capstyle='butt')
        ax.plot([rect_right, rect_right], [bracket_y - tick_len, bracket_y],
                color=NAVY, linewidth=1.5, solid_capstyle='butt')
        ax.plot([rect_left, rect_right], [bracket_y, bracket_y],
                color=NAVY, linewidth=1.5, solid_capstyle='butt')

        # Grid spacing annotation (first row only)
        if ax_idx == 0:
            hx1 = x_positions[center_idx]
            hx2 = x_positions[center_idx + 1]
            hy = y_base + 0.85

            ax.annotate('', xy=(hx2, hy), xytext=(hx1, hy),
                        arrowprops=dict(arrowstyle='<->', color=TEAL,
                                        lw=2, shrinkA=0, shrinkB=0,
                                        mutation_scale=15))
            ax.text((hx1 + hx2) / 2, hy + 0.15, r'$h$',
                    ha='center', va='bottom', fontsize=14,
                    color=TEAL, fontweight='bold')

    # Main title
    fig.suptitle('Finite Difference Stencils: From Local to Global',
                 fontsize=15, fontweight='bold', color=NAVY, y=0.97)

    # Bottom annotation
    fig.text(0.5, 0.02,
             r"As stencil width increases, accuracy improves. "
             r"The spectral method is the limit: all nodes contribute to $u'(x_i)$.",
             ha='center', fontsize=11, style='italic', color='gray')

    plt.tight_layout(rect=[0, 0.05, 1, 0.94])
    plt.subplots_adjust(hspace=0.35)

    # Save figure
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_FILE, bbox_inches='tight', pad_inches=0.15)
    fig.savefig(OUTPUT_FILE.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.15, dpi=300)

    print(f'Figure saved to: {OUTPUT_FILE.resolve()}')
    plt.close(fig)


if __name__ == '__main__':
    main()
