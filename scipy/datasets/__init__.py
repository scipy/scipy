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

from ._fetchers import face, ascent, electrocardiogram


__all__ = ['ascent', 'electrocardiogram', 'face']


from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
