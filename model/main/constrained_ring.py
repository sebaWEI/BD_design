"""Constrained closed-loop ring — ``main.tex`` § Constrained Closed-Loop Model."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._gaussian import as_array
from .gaussian_ring import GaussianRing

# kcal/mol at 37 °C
RT_DEFAULT = 0.616


@dataclass(frozen=True)
class ConstrainedRing:
    r"""Gaussian ring × 5′/3′ arm-length match constraint.

    \[
    C_{\\mathrm{eff}}^{\\mathrm{constrained}}(n)
    =
    C_{\\mathrm{eff}}^{\\mathrm{ring}}(n)\\,
    P_{\\mathrm{match}}(n)
    \]

    with
    \[
    P_{\\mathrm{match}}
    =
    \\exp\\!\\left[
    -\\frac{k_{\\mathrm{match}}}{2RT}
    (n_{5A} - n_{B3})^2
    \\right].
    \]

    Parameters
    ----------
    ring :
        Underlying ``GaussianRing`` (defines \(N\), \(b\), capture shell).
    n_5a :
        Effective Kuhn length of the 5′ arm (\(n_{5A}\)).
    k_match :
        Harmonic constraint strength \(k_{\\mathrm{match}}\) (energy units of ``rt``).
    rt :
        Thermal energy \(RT\) (same units as \(k_{\\mathrm{match}}\)).
    """

    ring: GaussianRing
    n_5a: float
    k_match: float
    rt: float = RT_DEFAULT

    def __post_init__(self) -> None:
        if self.n_5a <= 0:
            raise ValueError("n_5a must be positive")
        if self.k_match < 0:
            raise ValueError("k_match must be non-negative")
        if self.rt <= 0:
            raise ValueError("rt must be positive")

    @classmethod
    def from_params(
        cls,
        n_total: float,
        n_5a: float,
        k_match: float,
        *,
        kuhn_length_nm: float = 2.0,
        r_min: float = 1.0,
        r_max: float = 8.0,
        rt: float = RT_DEFAULT,
    ) -> ConstrainedRing:
        return cls(
            ring=GaussianRing(
                n_total=n_total,
                kuhn_length_nm=kuhn_length_nm,
                r_min=r_min,
                r_max=r_max,
            ),
            n_5a=n_5a,
            k_match=k_match,
            rt=rt,
        )

    @property
    def kappa(self) -> float:
        r"""Convenience: \(k_{\\mathrm{match}} / (2RT)\)."""
        return self.k_match / (2.0 * self.rt)

    def match_potential(self, n_b3: ArrayLike) -> NDArray[np.floating]:
        r"""\(U_{\\mathrm{match}}(n) = \\frac{k}{2}(n_{5A} - n_{B3})^2\)."""
        n_b3 = as_array(n_b3)
        return 0.5 * self.k_match * (self.n_5a - n_b3) ** 2

    def match_probability(self, n_b3: ArrayLike) -> NDArray[np.floating]:
        r"""\(P_{\\mathrm{match}}(n)\)."""
        return np.exp(-self.match_potential(n_b3) / self.rt)

    def ring_arc_segments(self, n_b3: ArrayLike) -> NDArray[np.floating]:
        r"""Default ring arc: BD → poly(A) → cap → AUG ≈ \(n_{B3} + n_{5A}\)."""
        return as_array(n_b3) + self.n_5a

    def effective_concentration(
        self,
        n_b3: ArrayLike,
        *,
        n_ring: ArrayLike | None = None,
    ) -> NDArray[np.floating]:
        r"""\(C_{\\mathrm{eff}}^{\\mathrm{constrained}}\) in mol/L (M).

        Parameters
        ----------
        n_b3 :
            Effective 3′-arm Kuhn length used in \(P_{\\mathrm{match}}\).
        n_ring :
            Arc length (Kuhn segments) for the Gaussian-ring \(C_{\\mathrm{eff}}^{\\mathrm{ring}}\).
            Defaults to \(n_{B3} + n_{5A}\).
        """
        if n_ring is None:
            n_ring = self.ring_arc_segments(n_b3)
        c_ring = self.ring.effective_concentration(n_ring)
        return c_ring * self.match_probability(n_b3)

    def site_score(
        self,
        n_b3: ArrayLike,
        delta_g_bd: ArrayLike,
        *,
        n_ring: ArrayLike | None = None,
    ) -> NDArray[np.floating]:
        r"""Trans site score ∝ \(C_{\\mathrm{eff}}^{\\mathrm{constrained}}\\exp(-\\Delta G_{BD}/RT)\)."""
        dg = as_array(delta_g_bd)
        return self.effective_concentration(n_b3, n_ring=n_ring) * np.exp(-dg / self.rt)
