"""Closed Gaussian ring — ``main.tex`` § Closed-Loop Positional Model."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._gaussian import (
    as_array,
    effective_concentration_molar,
    encounter_probability,
    shell_volume_nm3,
)


@dataclass(frozen=True)
class GaussianRing:
    r"""Gaussian ring for BD ↔ 5′ initiation encounter on a closed mRNA.

    Parameters
    ----------
    n_total :
        Total Kuhn segments on the closed ring \(N\).
    kuhn_length_nm :
        Kuhn length \(b\) (nm).
    r_min, r_max :
        BD–initiation capture shell (nm).
    """

    n_total: float
    kuhn_length_nm: float = 2.0
    r_min: float = 1.0
    r_max: float = 8.0

    def __post_init__(self) -> None:
        if self.n_total <= 1.0:
            raise ValueError("n_total must be > 1")
        if self.kuhn_length_nm <= 0:
            raise ValueError("kuhn_length_nm must be positive")
        if self.r_max <= self.r_min:
            raise ValueError("r_max must be > r_min")

    def _clip_n(self, n: ArrayLike) -> NDArray[np.floating]:
        return np.clip(as_array(n), 1e-6, self.n_total - 1e-6)

    def mean_square_separation(self, n: ArrayLike) -> NDArray[np.floating]:
        r"""⟨R²⟩ = b² n (N−n) / N for arc length \(n\) Kuhn segments."""
        n = self._clip_n(n)
        return (self.kuhn_length_nm**2) * n * (self.n_total - n) / self.n_total

    def density_at(self, r: ArrayLike, n: ArrayLike) -> NDArray[np.floating]:
        r"""Three-dimensional ring density \(p_{\\mathrm{ring}}(\\mathbf R\\mid n)\)."""
        r = as_array(r)
        n = self._clip_n(n)
        # Broadcast: allow scalar n with vector r or vice versa.
        s2 = self.mean_square_separation(n)
        pref = (3.0 / (2.0 * np.pi * s2)) ** 1.5
        return pref * np.exp(-3.0 * r**2 / (2.0 * s2))

    def radial_density(self, r: ArrayLike, n: ArrayLike) -> NDArray[np.floating]:
        r"""Radial density \(p_{\\mathrm{ring}}(R\\mid n)\)."""
        r = as_array(r)
        return 4.0 * np.pi * r**2 * self.density_at(r, n)

    def encounter_probability(self, n: ArrayLike) -> NDArray[np.floating]:
        r"""\(P_{\\mathrm{enc}}^{\\mathrm{ring}}(n)\)."""
        return encounter_probability(self.mean_square_separation(n), self.r_min, self.r_max)

    def capture_volume_nm3(self) -> float:
        r"""\(V_{\\mathrm{BD}}\) (nm³)."""
        return shell_volume_nm3(self.r_min, self.r_max)

    def effective_concentration(self, n: ArrayLike) -> NDArray[np.floating]:
        r"""\(C_{\\mathrm{eff}}^{\\mathrm{ring}}(n)\) in mol/L (M)."""
        return effective_concentration_molar(
            self.encounter_probability(n), self.r_min, self.r_max
        )
