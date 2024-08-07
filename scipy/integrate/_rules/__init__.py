"""Numerical cubature algorithms"""

from ._base import (
    Rule, FixedRule,
    NestedFixedRule, NestedRule,
    ProductFixed, ProductNestedFixed
)
from ._genz_malik import GenzMalikCubature
from ._newton_cotes import NewtonCotesQuadrature
from ._gauss_kronrod import GaussKronrodQuadrature
from ._gauss_legendre import GaussLegendreQuadrature

__all__ = [s for s in dir() if not s.startswith('_')]
