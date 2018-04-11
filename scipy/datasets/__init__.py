"""
================================
Datasets (:mod:`scipy.datasets`)
================================

.. currentmodule:: scipy.datasets

.. autosummary::
   :toctree: generated/

   ascent - Get example image for processing
   face - Get example image for processing
   electrocardiogram - Load an example of a one-dimensional signal.

"""

from __future__ import (division as _division,
                        print_function as _print_function,
                        absolute_import as _absolute_import)

from ._funcs import face, ascent, electrocardiogram


__all__ = ['face', 'ascent', 'electrocardiogram']


from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
