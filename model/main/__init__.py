r"""Polymer conformational models from ``model/main.tex``.

Classes (increasing structure):

1. ``FreelyJointedChain`` — open FJC → \(C_{\\mathrm{eff}}^{\\mathrm{loop}}\)
2. ``GaussianRing`` — closed Gaussian ring → \(C_{\\mathrm{eff}}^{\\mathrm{ring}}(n)\)
3. ``ConstrainedRing`` — ring × arm-match → \(C_{\\mathrm{eff}}^{\\mathrm{constrained}}(n)\)
"""

from .constrained_ring import ConstrainedRing
from .fjc import FreelyJointedChain
from .gaussian_ring import GaussianRing

__all__ = ["FreelyJointedChain", "GaussianRing", "ConstrainedRing"]
