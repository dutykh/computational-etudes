#!/usr/bin/env python3
"""
laplacian_polar.py - Discrete polar Laplacian on the unit disk

Constructs the spectral Laplacian operator using Chebyshev discretisation
in the radial direction (with the doubling trick r ∈ [-1,1]) and Fourier
discretisation in the angular direction.

Chapter 12: Spectral Methods in Polar Coordinates
Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
Last modified: February 2026
"""

import numpy as np
import sys
import os

# Add ch07 to path for cheb_matrix
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'ch07'))
from cheb_matrix import cheb_matrix


def laplacian_polar(Nr, M):
    """
    Build the discrete polar Laplacian on the unit disk.

    Uses Chebyshev discretisation in r ∈ [-1,1] (doubling trick, odd Nr)
    and Fourier discretisation in θ ∈ [0, 2π). Imposes u = 0 at r = 1
    (Dirichlet boundary condition).

    Parameters
    ----------
    Nr : int
        Number of Chebyshev radial points (must be odd, ≥ 3).
    M : int
        Number of Fourier angular points (must be even, ≥ 4).

    Returns
    -------
    L : ndarray, shape (N2*M, N2*M)
        Discrete Laplacian matrix on the interior grid.
    r_pos : ndarray, shape (N2,)
        Positive radial grid points (r > 0), sorted in decreasing order.
    theta : ndarray, shape (M,)
        Angular grid points θ_m = 2π m / M, m = 0, ..., M-1.
    """
    if Nr % 2 == 0:
        raise ValueError("Nr must be odd to avoid a grid point at r = 0.")
    if M % 2 != 0:
        raise ValueError("M must be even for the Fourier angular grid.")

    N2 = (Nr - 1) // 2  # number of positive-r interior points

    # --- Radial direction: Chebyshev on [-1, 1] ---
    D, r = cheb_matrix(Nr)
    D2full = D @ D  # second derivative matrix, size (Nr+1) x (Nr+1)

    # Block decomposition directly from the full (Nr+1)×(Nr+1) matrices.
    # Chebyshev points: r[0]=1 (boundary), r[1]...r[N2] (positive interior),
    # r[N2+1]...r[Nr-1] (negative interior), r[Nr]=-1 (boundary).
    # The j-th positive-r point r[j] pairs with -r[j] = r[Nr-j] under
    # the symmetry condition, so negative-r columns must be reversed.
    #
    # Positive-r rows:    indices 1, 2, ..., N2
    # Positive-r columns: indices 1, 2, ..., N2
    # Negative-r columns (reversed): indices Nr-1, Nr-2, ..., N2+1

    D1 = D2full[1:N2+1, 1:N2+1]               # r>0 rows, r>0 cols
    D2_blk = D2full[1:N2+1, Nr-1:N2:-1]       # r>0 rows, r<0 cols (reversed)
    E1 = D[1:N2+1, 1:N2+1]                    # first deriv: r>0 rows, r>0 cols
    E2 = D[1:N2+1, Nr-1:N2:-1]                # first deriv: r>0 rows, r<0 cols (reversed)

    # Positive radial points (interior only, excluding r=1 boundary)
    r_pos = r[1:N2+1]  # r > 0 points in decreasing order

    # Diagonal matrix R = diag(1/r_j) for positive r
    R = np.diag(1.0 / r_pos)

    # --- Angular direction: Fourier on [0, 2π) ---
    theta = 2 * np.pi * np.arange(M) / M
    dt = theta[1] - theta[0]  # angular spacing

    # Fourier second derivative matrix (Toeplitz)
    M2 = M // 2
    D2t = _fourier_second_deriv_matrix(M, dt)

    # --- Kronecker product assembly ---
    # L = (D1 + R·E1) ⊗ I_M + (D2_blk + R·E2) ⊗ S + R² ⊗ D_θ²
    #
    # where S is the block swap matrix that maps θ to θ+π:
    # S = [0_{M/2}  I_{M/2}]
    #     [I_{M/2}  0_{M/2}]
    # This implements the angular π-shift from the symmetry condition.

    I_M = np.eye(M)
    Z = np.zeros((M2, M2))
    I_M2 = np.eye(M2)
    S = np.block([[Z, I_M2],
                  [I_M2, Z]])

    L = (np.kron(D1 + R @ E1, I_M)
         + np.kron(D2_blk + R @ E2, S)
         + np.kron(R ** 2, D2t))

    return L, r_pos, theta


def _fourier_second_deriv_matrix(M, dt):
    """
    Build the M × M Fourier second-derivative Toeplitz matrix.

    Parameters
    ----------
    M : int
        Number of equispaced angular points (must be even).
    dt : float
        Angular spacing 2π/M.

    Returns
    -------
    D2t : ndarray, shape (M, M)
        Second-derivative matrix for periodic functions on [0, 2π).
    """
    from scipy.linalg import toeplitz

    # First column of the Toeplitz matrix
    col = np.zeros(M)
    col[0] = -np.pi**2 / (3 * dt**2) - 1.0 / 6.0
    for j in range(1, M):
        col[j] = 0.5 * (-1)**(j + 1) / np.sin(dt * j / 2)**2

    D2t = toeplitz(col)
    return D2t


if __name__ == "__main__":
    # Quick test: build Laplacian and check eigenvalues
    Nr, M = 25, 20
    L, r_pos, theta = laplacian_polar(Nr, M)
    N2 = (Nr - 1) // 2
    print(f"Nr = {Nr}, M = {M}")
    print(f"N2 = {N2} (positive-r interior points)")
    print(f"Laplacian size: {L.shape}")
    print(f"r_pos range: [{r_pos[-1]:.6f}, {r_pos[0]:.6f}]")

    # Eigenvalue check against Bessel function zeros
    eigenvalues = np.sort(np.real(np.linalg.eigvals(-L)))
    print(f"\nFirst 10 eigenvalues of -L:")
    for i in range(10):
        print(f"  λ_{i+1:2d} = {eigenvalues[i]:.10f}")

    # Bessel zeros squared for comparison:
    # j_{0,1}^2 ≈ 5.7832, j_{1,1}^2 ≈ 14.6820 (×2, degenerate),
    # j_{2,1}^2 ≈ 26.3746 (×2), j_{0,2}^2 ≈ 30.4713
    print(f"\nExpected (Bessel zeros squared):")
    print(f"  j_{{0,1}}^2 ≈  5.7832  (multiplicity 1)")
    print(f"  j_{{1,1}}^2 ≈ 14.6820  (multiplicity 2)")
    print(f"  j_{{2,1}}^2 ≈ 26.3746  (multiplicity 2)")
    print(f"  j_{{0,2}}^2 ≈ 30.4713  (multiplicity 1)")
