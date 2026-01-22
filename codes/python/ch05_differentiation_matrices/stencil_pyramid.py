#!/usr/bin/env python3
"""
stencil_pyramid.py

Visualizes Fornberg's recursive algorithm for computing finite difference
weights as a "stencil pyramid" showing how weights are computed progressively
from 1-point to N-point stencils.

The visualization shows:
    - The recursive structure of the algorithm
    - How weights for smaller stencils inform larger stencils
    - The data dependencies in the computation

Key Insight:
    The algorithm builds up weights level by level:
    - Level 0: Single node (weight = 1 for interpolation)
    - Level 1: Two nodes (weights from level 0)
    - Level 2: Three nodes (weights from level 1)
    - ...
    - Level n: n+1 nodes (weights from level n-1)

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
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
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
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'
LIGHT_GRAY = '#E8E8E8'
MEDIUM_GRAY = '#CCCCCC'

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch05' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'stencil_pyramid.pdf'


# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
def main():
    fig, ax = plt.subplots(figsize=(9, 6))

    # Parameters
    n_levels = 5  # Number of levels (m = 0, 1, 2, 3, 4)
    box_width = 0.8
    box_height = 0.5
    h_spacing = 1.0
    v_spacing = 1.2

    # Draw boxes for each level
    boxes = {}  # Store box positions for drawing arrows

    for m in range(n_levels):  # m = stencil size - 1 (number of intervals)
        n_boxes = m + 1  # Number of nodes at this level
        y = (n_levels - 1 - m) * v_spacing

        # Center the row
        total_width = n_boxes * box_width + (n_boxes - 1) * (h_spacing - box_width)
        x_start = -total_width / 2 + box_width / 2

        for k in range(n_boxes):  # k = node index
            x = x_start + k * h_spacing

            # Choose color based on position
            if m == n_levels - 1:
                # Final level (result) - highlight
                color = NAVY
                text_color = 'white'
                alpha = 1.0
            elif m == n_levels - 2:
                # Previous level - medium highlight
                color = SKY
                text_color = 'white'
                alpha = 0.9
            else:
                # Earlier levels - faded
                color = LIGHT_GRAY
                text_color = NAVY
                alpha = 0.8

            # Draw box
            box = FancyBboxPatch(
                (x - box_width/2, y - box_height/2),
                box_width, box_height,
                boxstyle="round,pad=0.02,rounding_size=0.1",
                facecolor=color,
                edgecolor=NAVY,
                linewidth=1.5,
                alpha=alpha
            )
            ax.add_patch(box)

            # Add text label
            if m <= 2:
                label = f'$\\delta_{{{m}}}^{{({k},{m})}}$'
            else:
                label = f'$\\delta^{{({k},{m})}}$'

            ax.text(x, y, label, ha='center', va='center',
                    fontsize=9, color=text_color, fontweight='bold')

            boxes[(m, k)] = (x, y)

    # Draw arrows showing dependencies
    arrow_style = "Simple, tail_width=0.3, head_width=2, head_length=3"

    for m in range(1, n_levels):
        n_boxes = m + 1
        for k in range(n_boxes):
            x_to, y_to = boxes[(m, k)]

            # Arrows from previous level
            if k < m:
                # This node existed in previous stencil
                x_from, y_from = boxes[(m-1, k)]
                # Draw arrow
                arrow = FancyArrowPatch(
                    (x_from, y_from - box_height/2 - 0.05),
                    (x_to, y_to + box_height/2 + 0.05),
                    arrowstyle='->', mutation_scale=10,
                    color=TEAL, alpha=0.6, linewidth=1.2
                )
                ax.add_patch(arrow)

            if k == m and m > 0:
                # This is the "new" node - depends on the rightmost node of previous level
                x_from, y_from = boxes[(m-1, m-1)]
                arrow = FancyArrowPatch(
                    (x_from + box_width/4, y_from - box_height/2 - 0.05),
                    (x_to - box_width/4, y_to + box_height/2 + 0.05),
                    arrowstyle='->', mutation_scale=10,
                    color=CORAL, alpha=0.7, linewidth=1.5,
                    connectionstyle="arc3,rad=0.2"
                )
                ax.add_patch(arrow)

    # Add level labels on the left
    for m in range(n_levels):
        y = (n_levels - 1 - m) * v_spacing
        x = -4.5
        ax.text(x, y, f'$m = {m}$', ha='left', va='center',
                fontsize=11, color=NAVY)
        ax.text(x + 1.2, y, f'({m+1} node{"s" if m > 0 else ""})',
                ha='left', va='center', fontsize=9, color='gray')

    # Add title and annotations
    ax.set_title("Fornberg's Recursive Algorithm: The Stencil Pyramid",
                 fontsize=13, fontweight='bold', pad=20)

    # Add explanation text at bottom
    explanation = (
        "Each box represents a weight $\\delta_j^{(k,m)}$ for node $k$ "
        "in a stencil of $m+1$ nodes.\n"
        "Arrows show data dependencies: weights at level $m$ "
        "depend on weights from level $m-1$."
    )
    ax.text(0, -0.8, explanation, ha='center', va='top',
            fontsize=10, color='gray', style='italic',
            transform=ax.transData)

    # Add legend for arrow colors
    ax.plot([], [], '-', color=TEAL, linewidth=2, alpha=0.6,
            label='Update existing node')
    ax.plot([], [], '-', color=CORAL, linewidth=2, alpha=0.7,
            label='Compute new node (with $\\beta$ factor)')
    ax.legend(loc='upper right', fontsize=9, framealpha=0.9)

    # Set axis limits and remove axes
    ax.set_xlim(-5.5, 5.5)
    ax.set_ylim(-1.5, (n_levels - 1) * v_spacing + 1)
    ax.set_aspect('equal')
    ax.axis('off')

    plt.tight_layout()

    # Save figure
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_FILE, bbox_inches='tight', pad_inches=0.1)
    fig.savefig(OUTPUT_FILE.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.1, dpi=300)

    print(f'Figure saved to: {OUTPUT_FILE.resolve()}')

    plt.close(fig)


if __name__ == '__main__':
    main()
