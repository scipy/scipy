"""Numerical cubature algorithms"""

from ._base import (
    Rule, FixedRule,
    NestedFixedRule,
    ProductNestedFixed,
)
from ._genz_malik import GenzMalikCubature
from ._gauss_kronrod import GaussKronrodQuadrature
from ._gauss_legendre import GaussLegendreQuadrature

__all__ = [
    'Rule', 'FixedRule', 'NestedFixedRule', 'ProductNestedFixed',
    'GenzMalikCubature', 'GaussKronrodQuadrature', 'GaussLegendreQuadrature',
]
