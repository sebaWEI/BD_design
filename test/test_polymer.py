"""Unit tests for main.tex polymer classes."""

from __future__ import annotations

import numpy as np

from model.main import ConstrainedRing, FreelyJointedChain, GaussianRing


def test_fjc_encounter_and_ceff() -> None:
    chain = FreelyJointedChain(n_segments=100.0, kuhn_length_nm=2.0, r_min=1.0, r_max=8.0)
    p = chain.encounter_probability()
    assert 0.0 < p <= 1.0
    c = chain.effective_concentration()
    assert c > 0.0
    # denser / shorter chains → higher encounter in a fixed shell (rough check)
    short = FreelyJointedChain(n_segments=20.0, kuhn_length_nm=2.0)
    long = FreelyJointedChain(n_segments=500.0, kuhn_length_nm=2.0)
    assert short.encounter_probability() > long.encounter_probability()


def test_fjc_radial_normalizes_roughly() -> None:
    chain = FreelyJointedChain(n_segments=80.0)
    rho = np.linspace(0.0, 80.0, 4000)
    # trapezoid integral of p_ρ should be ≈ 1
    integ = np.trapezoid(chain.radial_density(rho), rho)
    assert abs(integ - 1.0) < 0.02


def test_gaussian_ring_peak_near_half() -> None:
    ring = GaussianRing(n_total=120.0, kuhn_length_nm=2.0)
    n = np.linspace(5.0, 115.0, 45)
    c = ring.effective_concentration(n)
    assert c.min() > 0.0
    # ⟨R²⟩ is maximized near n = N/2 → encounter (fixed shell) is typically lower
    # than near the ends; C_eff should not be flat.
    assert c.max() / c.min() > 1.05


def test_constrained_ring_match_peaks_at_n5a() -> None:
    model = ConstrainedRing.from_params(
        n_total=120.0, n_5a=20.0, k_match=0.024, rt=1.0,
    )
    n = np.array([10.0, 20.0, 30.0])
    p = model.match_probability(n)
    assert int(np.argmax(p)) == 1
    assert float(p[1]) == float(np.max(p))
    c = model.effective_concentration(n)
    assert int(np.argmax(c)) == 1


def test_constrained_reduces_to_ring_when_k_zero() -> None:
    ring = GaussianRing(n_total=100.0)
    constrained = ConstrainedRing(ring=ring, n_5a=25.0, k_match=0.0)
    n_b3 = np.array([10.0, 25.0, 40.0])
    n_ring = n_b3 + 25.0
    np.testing.assert_allclose(
        constrained.effective_concentration(n_b3),
        ring.effective_concentration(n_ring),
    )


def test_site_score_prefers_favorable_dg() -> None:
    model = ConstrainedRing.from_params(n_total=120.0, n_5a=20.0, k_match=0.02, rt=0.616)
    n = np.array([20.0, 20.0])
    dg = np.array([-10.0, -5.0])
    score = model.site_score(n, dg)
    assert score[0] > score[1]
