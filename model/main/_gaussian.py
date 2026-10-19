"""Shared Gaussian-chain encounter integrals (main.tex §§1–3)."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.stats import chi2

NA = 6.02214076e23  # Avogadro / mol
# 1 nm³ = 1e-24 L → molarity when dividing by NA * V_L
NM3_TO_L = 1e-24


def as_array(x: ArrayLike) -> NDArray[np.floating]:
    return np.asarray(x, dtype=float)


def shell_volume_nm3(r_min: float, r_max: float) -> float:
    r"""\(V = \\frac{4\\pi}{3}(R_{\\max}^3 - R_{\\min}^3)\) in nm³."""
    if r_max <= r_min:
        raise ValueError("r_max must be > r_min")
    return (4.0 / 3.0) * np.pi * (r_max**3 - r_min**3)


def encounter_probability(s2_nm2: ArrayLike, r_min: float, r_max: float) -> NDArray[np.floating]:
    """P(r_min ≤ R ≤ r_max) for isotropic Gaussian with ⟨R²⟩ = s2.

    Each Cartesian component is N(0, s2/3), so 3R²/s2 ~ χ²_3.
    """
    s2 = np.clip(as_array(s2_nm2), 1e-30, None)
    hi = chi2.cdf(3.0 * r_max**2 / s2, df=3)
    lo = chi2.cdf(3.0 * r_min**2 / s2, df=3)
    return np.clip(hi - lo, 0.0, 1.0)


def effective_concentration_molar(
    p_enc: ArrayLike,
    r_min: float,
    r_max: float,
) -> NDArray[np.floating]:
    r"""Shell-averaged effective concentration in mol/L (M).

    \(C_{\\mathrm{eff}} = P_{\\mathrm{enc}} / (N_A V)\) with \(V\) in litres.
    """
    v_l = shell_volume_nm3(r_min, r_max) * NM3_TO_L
    return as_array(p_enc) / (NA * v_l)
