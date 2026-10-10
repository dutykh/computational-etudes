#!/usr/bin/env python3
"""
eye_aliasing.py
Chapter 9: Physical and Fourier Space on Grids

Demonstrates visual aliasing: sin(n) for integer n appears to oscillate
slowly when plotted, because the eye aliases the high frequency to a
low-frequency pattern.

Author: Dr. Denys Dutykh
Date: February 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def main():
    # Integer arguments: sin(n) for n = 1, 2, ..., 3000
    n = np.arange(1, 3001)
    y = np.sin(n)

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 3.5))

    ax.plot(n, y, '.', color='#142D6E', markersize=0.8)

    ax.set_xlabel(r'$n$', fontsize=12)
    ax.set_ylabel(r'$\sin(n)$', fontsize=12)
    ax.set_title(r'Visual aliasing: $\sin(n)$ for integer $n$ appears to oscillate slowly',
                 fontsize=12)
    ax.set_xlim(0, 3000)
    ax.set_ylim(-1.3, 1.3)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save figure
    output_dir = Path(__file__).parent.parent.parent.parent / 'textbook' / 'figures' / 'ch09' / 'python'
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / 'eye_aliasing.pdf', bbox_inches='tight')
    print(f"Figure saved to {output_dir / 'eye_aliasing.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
