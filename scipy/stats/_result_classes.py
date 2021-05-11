# This module exists only to allow Sphinx to generate docs
# for the result objects returned by some functions in stats.

"""
Result classes
--------------

.. currentmodule:: scipy.stats._result_classes

.. autosummary::
   :toctree: generated/

   BinomTestResult
   OddsRatioResult
   RelativeRiskResult

"""

__all__ = ['BinomTestResult', 'OddsRatioResult', 'RelativeRiskResult']


from ._binomtest import BinomTestResult
from ._odds_ratio import OddsRatioResult
from ._relative_risk import RelativeRiskResult
