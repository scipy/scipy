# This module exists only to allow Sphinx to generate docs
# for the result objects returned by some functions in stats
# _without_ adding them to the main stats documentation page.

"""
Result classes used in :mod:`scipy.stats`
-----------------------------------------

.. currentmodule:: scipy.stats._result_classes

.. autosummary::
   :toctree: generated/

   RelativeRiskResult
   BinomTestResult
   TukeyHSDResult
   PearsonRResult
   FitResult

Warnings / Errors used in :mod:`scipy.stats`
--------------------------------------------

.. currentmodule:: scipy.stats._result_classes

.. autosummary::
   :toctree: generated/

   FitError

Warnings / Errors used in :mod:`scipy.stats`
--------------------------------------------

.. currentmodule:: scipy.stats._result_classes

.. autosummary::
   :toctree: generated/

   BootstrapDegenerateDistributionWarning
   F_onewayConstantInputWarning
   F_onewayBadInputSizesWarning
   PearsonRConstantInputWarning
   PearsonRNearConstantInputWarning
   SpearmanRConstantInputWarning

"""
# Import classes below to allow Sphinx to render the documentation above
#
# Add classes to __all__ to make it easy to import them from the private
# `scipy.stats_result_classes` namespace.
#
# Classes included above but not added to __all__ here have already been added
# to the `stats` namespace. This should be done sparingly, and as with any
# addition to the public API, must be approved by the scipy-dev mailing list.

# Results
__all__ = ['BinomTestResult', 'RelativeRiskResult', 'TukeyHSDResult',
           'PearsonRResult', 'FitResult']
# Warnings
__all__ += ['BootstrapDegenerateDistributionWarning',]

# Results
from ._binomtest import BinomTestResult
from ._relative_risk import RelativeRiskResult
from ._hypotests import TukeyHSDResult
from ._stats_py import PearsonRResult
from ._fit import FitResult

# Warnings
from ._resampling import BootstrapDegenerateDistributionWarning
