#!/usr/bin/env python3
"""
chebyshev_points_circle.py

Visualizes the geometric construction of Chebyshev points from the unit circle.

Chebyshev-Gauss-Lobatto points are obtained by:
1. Placing N+1 equally spaced points on the upper half of the unit circle
2. Projecting these points vertically onto the x-axis

This geometric interpretation reveals why Chebyshev points cluster near ±1:
equal spacing on the circle maps to denser spacing near the boundaries when
projected horizontally.

The point density on [-1,1] follows:
    ρ(x) ≈ N / (π√(1-x²))

which diverges at the boundaries, exactly compensating for the growth of
interpolation error there.

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
from matplotlib.patches import Arc, FancyArrowPatch
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
rcParams['text.usetex'] = False  # Set to True if LaTeX is available
rcParams['axes.linewidth'] = 0.8
rcParams['xtick.major.width'] = 0.8
rcParams['ytick.major.width'] = 0.8

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
N = 8  # Number of intervals (N+1 points)

# Book color scheme
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'
LIGHT_GRAY = '#CCCCCC'

# Output path
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch04' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'chebyshev_points_circle.pdf'


# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
def main():
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 5))

    # Draw unit circle (upper half)
    theta_circle = np.linspace(0, np.pi, 200)
    x_circle = np.cos(theta_circle)
    y_circle = np.sin(theta_circle)
    ax.plot(x_circle, y_circle, color=NAVY, linewidth=1.5, label='Unit circle')

    # Draw x-axis
    ax.axhline(y=0, color='#666666', linewidth=0.8, linestyle='-')
    ax.plot([-1.3, 1.3], [0, 0], color='#666666', linewidth=0.8)

    # Chebyshev points: angles equally spaced on [0, π]
    j = np.arange(N + 1)
    theta_j = j * np.pi / N

    # Points on the circle
    x_circle_pts = np.cos(theta_j)
    y_circle_pts = np.sin(theta_j)

    # Points on the x-axis (Chebyshev nodes)
    x_cheb = np.cos(theta_j)
    y_cheb = np.zeros_like(x_cheb)

    # Draw projection lines
    for i in range(N + 1):
        ax.plot([x_circle_pts[i], x_cheb[i]], [y_circle_pts[i], y_cheb[i]],
                color=LIGHT_GRAY, linewidth=0.8, linestyle='--', zorder=1)

    # Draw points on circle
    ax.scatter(x_circle_pts, y_circle_pts, color=SKY, s=60, zorder=3,
               edgecolors=NAVY, linewidths=1, label='Points on circle')

    # Draw Chebyshev points on axis
    ax.scatter(x_cheb, y_cheb, color=CORAL, s=80, zorder=3,
               edgecolors='#B03A2E', linewidths=1.5, marker='s',
               label='Chebyshev points')

    # Draw a single clear angle arc for j=4 (middle point)
    arc_idx = N // 2  # j=4 for N=8
    arc_radius = 0.4
    arc = Arc((0, 0), 2*arc_radius, 2*arc_radius, angle=0, theta1=0,
              theta2=np.degrees(theta_j[arc_idx]), color=TEAL, linewidth=1.2)
    ax.add_patch(arc)

    # Annotate the angle with better positioning
    arc_label_angle = theta_j[arc_idx] / 2  # midpoint of the arc
    arc_label_x = (arc_radius + 0.15) * np.cos(arc_label_angle)
    arc_label_y = (arc_radius + 0.15) * np.sin(arc_label_angle)
    ax.annotate(r'$\theta_j = \frac{j\pi}{N}$', xy=(arc_label_x, arc_label_y),
                fontsize=10, color=TEAL, ha='left', va='center')

    # Label selected points - avoid crowding by choosing well-spaced indices
    # For N=8: label j=0, 2, 4, 6, 8 but position them carefully
    label_indices = [0, 2, N//2, N-2, N]
    for i in label_indices:
        if i < len(x_cheb):
            # Position labels below the x-axis for all points
            if i == N // 2:  # Center point (j=4)
                offset_y = -0.15
                offset_x = 0.08
            elif i == 0:  # j=0 (rightmost)
                offset_y = -0.18
                offset_x = 0
            elif i == N:  # j=N (leftmost)
                offset_y = -0.18
                offset_x = 0
            else:
                offset_y = 0.18
                offset_x = 0
            ax.annotate(f'$j = {i}$', xy=(x_cheb[i], y_cheb[i]),
                        xytext=(x_cheb[i] + offset_x, offset_y),
                        fontsize=9, ha='center', color='#444444')

    # Add formula box
    textstr = (r'$x_j = \cos\left(\frac{j\pi}{N}\right)$' + '\n' +
               r'$j = 0, 1, \ldots, N$')
    props = dict(boxstyle='round', facecolor='white', alpha=0.9,
                 edgecolor=NAVY, linewidth=0.5)
    ax.text(0.98, 0.95, textstr, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', horizontalalignment='right', bbox=props)

    # Add explanation annotations with arrows for clarity
    # "Equal spacing on circle" - point to an arc segment on the circle
    mid_arc_idx = 3  # Point to arc between j=2 and j=3
    mid_arc_x = np.cos(theta_j[mid_arc_idx] - np.pi/(2*N))
    mid_arc_y = np.sin(theta_j[mid_arc_idx] - np.pi/(2*N))
    ax.annotate('Equal spacing\non circle', xy=(mid_arc_x, mid_arc_y),
                xytext=(-0.35, 0.95), fontsize=9, color=SKY,
                ha='center', style='italic',
                arrowprops=dict(arrowstyle='->', color=SKY, lw=0.8,
                                connectionstyle='arc3,rad=0.2'))

    # "Clustering near boundaries" - point to the dense region near x=1
    ax.annotate('Clustering\nnear boundaries', xy=(0.85, 0),
                xytext=(0.65, -0.28), fontsize=9, color=CORAL,
                ha='center', style='italic',
                arrowprops=dict(arrowstyle='->', color=CORAL, lw=0.8,
                                connectionstyle='arc3,rad=-0.2'))

    # Styling
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-0.45, 1.15)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$x$')
    ax.set_title(f'Geometric Construction of Chebyshev Points ($N = {N}$)',
                 fontsize=11)

    # Remove unnecessary spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # Custom axis
    ax.set_yticks([0, 0.5, 1])
    ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax.tick_params(colors='#444444')

    ax.legend(loc='upper left', frameon=True, fancybox=False,
              edgecolor='none', facecolor='white', framealpha=0.9)

    plt.tight_layout()

    # Save figure
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_FILE, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(OUTPUT_FILE.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {OUTPUT_FILE.resolve()}')

    # Print node positions
    print(f"\nChebyshev-Gauss-Lobatto points for N = {N}:")
    print("-" * 40)
    print(f"{'j':>4} {'θ_j':>12} {'x_j':>12}")
    print("-" * 40)
    for i in range(N + 1):
        theta_str = f"{i}π/{N}" if i > 0 and i < N else ("0" if i == 0 else "π")
        print(f"{i:>4} {theta_str:>12} {x_cheb[i]:>12.6f}")
    print("-" * 40)

    # Show clustering effect
    print("\nSpacing between consecutive points:")
    print("-" * 40)
    spacings = np.diff(np.sort(x_cheb))
    for i, s in enumerate(spacings):
        print(f"Interval {i}: Δx = {s:.6f}")
    print("-" * 40)
    print(f"Ratio (boundary/center): {spacings[0]/spacings[N//2]:.2f}")

    plt.close(fig)


if __name__ == '__main__':
    main()
