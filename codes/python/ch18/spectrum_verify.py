#!/usr/bin/env python3
"""
spectrum_verify.py

Chapter 18: Linear Spectral Eigenproblems --- shared utility.

Drift-with-N diagnostic of Boyd (2000, Ch. 7.5). Given two spectra
lam1 (coarse) and lam2 (fine), classify each coarse mode as trusted or
suspect by measuring how much it moves between the two truncations,
scaled by the local intermodal separation.

Implements:
    sigma_j              (Boyd Eq. 7.19) intermodal separation
    delta_j,ordinal      (Boyd Eq. 7.20) scaled ordinal drift
    delta_j,nearest      (Boyd Eq. 7.21) scaled nearest drift

A mode is flagged trusted iff max(delta_ord, delta_nst) < tol.

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np


def intermodal_separation(lam: Sequence[float]) -> np.ndarray:
    """Return sigma_j per Boyd Eq. 7.19 for a sorted real spectrum."""
    lam = np.asarray(lam, dtype=float)
    n = lam.size
    sigma = np.zeros(n)
    if n == 1:
        sigma[0] = abs(lam[0])
        return sigma
    sigma[0] = abs(lam[0] - lam[1])
    if n > 2:
        diffs = np.abs(np.diff(lam))          # length n-1
        sigma[1:-1] = 0.5 * (diffs[:-1] + diffs[1:])
    sigma[-1] = abs(lam[-1] - lam[-2])
    tiny = sigma < 1e-14
    sigma[tiny] = np.maximum(np.abs(lam[tiny]), 1e-14)
    return sigma


def ordinal_drift(lam1: Sequence[float], lam2: Sequence[float]) -> np.ndarray:
    """Scaled ordinal drift delta_j,ordinal per Boyd Eq. 7.20."""
    lam1 = np.asarray(lam1, dtype=float)
    lam2 = np.asarray(lam2, dtype=float)
    sigma = intermodal_separation(lam1)
    m = min(lam1.size, lam2.size)
    delta = np.full(lam1.size, np.inf)
    delta[:m] = np.abs(lam1[:m] - lam2[:m]) / sigma[:m]
    return delta


def nearest_drift(lam1: Sequence[float], lam2: Sequence[float]) -> np.ndarray:
    """Scaled nearest drift delta_j,nearest per Boyd Eq. 7.21.

    For each lam1[j], find the lam2[k] closest to it and divide by sigma_j.
    Vectorised via broadcasting; O(n1 * n2) memory but fine for spectra
    of the size this chapter computes (n <= a few hundred).
    """
    lam1 = np.asarray(lam1, dtype=float)
    lam2 = np.asarray(lam2, dtype=float)
    sigma = intermodal_separation(lam1)
    # |lam1[j] - lam2[k]|  shape (n1, n2)
    diffs = np.abs(lam1[:, None] - lam2[None, :])
    mn = diffs.min(axis=1) if lam2.size else np.full(lam1.size, np.inf)
    return mn / sigma


@dataclass
class SpectrumReport:
    """Output of verify_spectrum. Attributes mirror the Julia struct."""
    lam1: np.ndarray
    lam2: np.ndarray
    sigma: np.ndarray
    delta_ordinal: np.ndarray
    delta_nearest: np.ndarray
    trusted: np.ndarray
    tol: float
    n_trusted: int = field(init=False)

    def __post_init__(self):
        self.n_trusted = int(self.trusted.sum())

    def to_dict(self) -> dict:
        """Serialise to plain dict for JSON dump in the validation harness."""
        return {
            "lam1": self.lam1.tolist(),
            "lam2": self.lam2.tolist(),
            "sigma": self.sigma.tolist(),
            "delta_ordinal": self.delta_ordinal.tolist(),
            "delta_nearest": self.delta_nearest.tolist(),
            "trusted": self.trusted.astype(int).tolist(),
            "tol": self.tol,
            "n_trusted": self.n_trusted,
        }


def verify_spectrum(lam1: Sequence[float],
                    lam2: Sequence[float],
                    tol: float = 1e-3,
                    sort_input: bool = True) -> SpectrumReport:
    """Drift diagnostic. Returns a SpectrumReport.

    A mode is trusted iff both ordinal and nearest drifts are < tol.
    Demanding both is conservative: ordinal alone is fooled by mode
    reordering; nearest alone can be fooled when two bad eigenvalues
    coincidentally sit near each other.

    Input is assumed real and sorted. For complex spectra, sort by real
    part or magnitude upstream and pass the real projection.
    """
    lam1 = np.asarray(lam1, dtype=float)
    lam2 = np.asarray(lam2, dtype=float)
    if sort_input:
        lam1 = np.sort(lam1)
        lam2 = np.sort(lam2)
    sigma = intermodal_separation(lam1)
    d_ord = ordinal_drift(lam1, lam2)
    d_nst = nearest_drift(lam1, lam2)
    trusted = (d_ord < tol) & (d_nst < tol) & np.isfinite(d_ord) & np.isfinite(d_nst)
    return SpectrumReport(
        lam1=lam1, lam2=lam2, sigma=sigma,
        delta_ordinal=d_ord, delta_nearest=d_nst,
        trusted=trusted, tol=float(tol),
    )


if __name__ == "__main__":
    # quick sanity check: Laplacian analytic spectrum vs perturbed spectrum
    j = np.arange(1, 21)
    lam_exact = (j * np.pi / 2) ** 2
    lam_perturbed = lam_exact.copy()
    # leave first 12 modes alone, corrupt the last 8 (emulating under-resolution)
    lam_perturbed[12:] *= 1.5
    report = verify_spectrum(lam_exact, lam_perturbed, tol=1e-2)
    print(f"Trusted modes: {report.n_trusted} of {lam_exact.size}")
    print("Ordinal drift:", np.round(report.delta_ordinal, 4))
