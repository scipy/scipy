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
import warnings

try:
    # https://github.com/scipy/scipy/pull/15607#issuecomment-1176457275
    # TODO: Remove warning filter after next certifi release
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=DeprecationWarning)
        import pooch
except ImportError:
    pooch = None

msg = (
    "Missing optional dependency 'pooch' required for scipy.datasets module. "
    "Please use pip or conda to install 'pooch'."
)
if pooch is None:
    raise ImportError(msg)


from ._fetchers import face, ascent, electrocardiogram  # noqa: E402
__all__ = ['ascent', 'electrocardiogram', 'face']


from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
