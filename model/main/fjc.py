"""Freely jointed chain (open) — ``main.tex`` § Polymer Conformational Model."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._gaussian import (
    effective_concentration_molar,
    encounter_probability,
    shell_volume_nm3,
)


@dataclass(frozen=True)
class FreelyJointedChain:
    r"""Open FJC with Gaussian end-to-end statistics.

    Parameters
    ----------
    n_segments :
        Number of Kuhn segments \(N = L/b\).
    kuhn_length_nm :
        Kuhn length \(b\) (nm).
    r_min, r_max :
        Protein-compatible capture shell (nm).
    """

    n_segments: float
    kuhn_length_nm: float = 2.0
    r_min: float = 1.0
    r_max: float = 8.0

    def __post_init__(self) -> None:
        if self.n_segments <= 0:
            raise ValueError("n_segments must be positive")
        if self.kuhn_length_nm <= 0:
            raise ValueError("kuhn_length_nm must be positive")
        if self.r_max <= self.r_min:
            raise ValueError("r_max must be > r_min")

    @property
    def end_to_end_msq(self) -> float:
        """⟨ρ²⟩ = N b² (nm²)."""
        return float(self.n_segments * self.kuhn_length_nm**2)

    def density_at(self, rho: ArrayLike) -> NDArray[np.floating]:
        r"""Three-dimensional Gaussian density \(p(\\boldsymbol{\\rho})\)."""
        rho = np.asarray(rho, dtype=float)
        nb2 = self.n_segments * self.kuhn_length_nm**2
        pref = (3.0 / (2.0 * np.pi * nb2)) ** 1.5
        return pref * np.exp(-3.0 * rho**2 / (2.0 * nb2))

    def radial_density(self, rho: ArrayLike) -> NDArray[np.floating]:
        r"""Radial density \(p_\\rho(\\rho) = 4\\pi\\rho^2 p(\\boldsymbol{\\rho})\)."""
        rho = np.asarray(rho, dtype=float)
        return 4.0 * np.pi * rho**2 * self.density_at(rho)

    def encounter_probability(self) -> float:
        r"""\(P_{\\mathrm{enc}}^{\\mathrm{loop}} = \\int_{\\rho_{\\min}}^{\\rho_{\\max}} p_\\rho(\\rho)\\,d\\rho\)."""
        return float(encounter_probability(self.end_to_end_msq, self.r_min, self.r_max))

    def capture_volume_nm3(self) -> float:
        r"""\(V_{\\mathrm{loop}}\) (nm³)."""
        return shell_volume_nm3(self.r_min, self.r_max)

    def mean_spatial_density(self) -> float:
        r"""\(\\bar p_{\\mathrm{loop}} = P_{\\mathrm{enc}} / V_{\\mathrm{loop}}\) (nm⁻³)."""
        return self.encounter_probability() / self.capture_volume_nm3()

    def effective_concentration(self) -> float:
        r"""\(C_{\\mathrm{eff}}^{\\mathrm{loop}}\) in mol/L (M)."""
        return float(
            effective_concentration_molar(
                self.encounter_probability(), self.r_min, self.r_max
            )
        )
