"""
SciPy: A scientific computing package for Python
================================================

Documentation is available in the docstrings and
online at http://docs.scipy.org.

Contents
--------
SciPy imports all the functions from the NumPy namespace, and in
addition provides:

Subpackages
-----------
::

 odr                          --- Orthogonal Distance Regression [*]
 misc                         --- Various utilities that don't have
                                  another home.
 cluster                      --- Vector Quantization / Kmeans [*]
 fftpack                      --- Discrete Fourier Transform algorithms
                                  [*]
 io                           --- Data input and output [*]
 sparse.linalg.eigen.lobpcg   --- Locally Optimal Block Preconditioned
                                  Conjugate Gradient Method (LOBPCG) [*]
 special                      --- Airy Functions [*]
 lib.blas                     --- Wrappers to BLAS library [*]
 sparse.linalg.eigen          --- Sparse Eigenvalue Solvers [*]
 stats                        --- Statistical Functions [*]
 lib                          --- Python wrappers to external libraries
                                  [*]
 lib.lapack                   --- Wrappers to LAPACK library [*]
 maxentropy                   --- Routines for fitting maximum entropy
                                  models [*]
 integrate                    --- Integration routines [*]
 ndimage                      --- n-dimensional image package [*]
 linalg                       --- Linear algebra routines [*]
 spatial                      --- Spatial data structures and algorithms
                                  [*]
 interpolate                  --- Interpolation Tools [*]
 sparse.linalg                --- Sparse Linear Algebra [*]
 sparse.linalg.dsolve.umfpack --- :Interface to the UMFPACK library: [*]
 sparse.linalg.dsolve         --- Linear Solvers [*]
 optimize                     --- Optimization Tools [*]
 sparse.linalg.eigen.arpack   --- Eigenvalue solver using iterative
                                  methods. [*]
 signal                       --- Signal Processing Tools [*]
 sparse                       --- Sparse Matrices [*]

 [*] - using a package requires explicit import

Global symbols from subpackages
-------------------------------
::

 misc                  --> info, factorial, factorial2, factorialk,
                           comb, who, lena, central_diff_weights,
                           derivative, pade, source
 fftpack               --> fft, fftn, fft2, ifft, ifft2, ifftn,
                           fftshift, ifftshift, fftfreq
 stats                 --> find_repeats
 linalg.dsolve.umfpack --> UmfpackContext

Utility tools
-------------
::

 test              --- Run scipy unittests
 show_config       --- Show scipy build configuration
 show_numpy_config --- Show numpy build configuration
 __version__       --- Scipy version string
 __numpy_version__ --- Numpy version string

"""

__all__ = ['pkgload','test']

from numpy import show_config as show_numpy_config
if show_numpy_config is None:
    raise ImportError,"Cannot import scipy when running from numpy source directory."
from numpy import __version__ as __numpy_version__

# Import numpy symbols to scipy name space
import numpy as _num
from numpy import oldnumeric
from numpy import *
from numpy.random import rand, randn
from numpy.fft import fft, ifft
from numpy.lib.scimath import *

# Emit a warning if numpy is too old
majver, minver = [float(i) for i in _num.version.version.split('.')[:2]]
if majver < 1 or (majver == 1 and minver < 2):
    import warnings
    warnings.warn("Numpy 1.2.0 or above is recommended for this version of " \
                  "scipy (detected version %s)" % _num.version.version,
                  UserWarning)

__all__ += ['oldnumeric']+_num.__all__

__all__ += ['randn', 'rand', 'fft', 'ifft']

del _num
# Remove the linalg imported from numpy so that the scipy.linalg package can be
# imported.
del linalg
__all__.remove('linalg')

try:
    from scipy.__config__ import show as show_config
except ImportError:
    msg = """Error importing scipy: you cannot import scipy while
    being in scipy source directory; please exit the scipy source
    tree first, and relaunch your python intepreter."""
    raise ImportError(msg)
from scipy.version import version as __version__

# Load scipy packages and their global_symbols
from numpy._import_tools import PackageLoader
import os as _os
SCIPY_IMPORT_VERBOSE = int(_os.environ.get('SCIPY_IMPORT_VERBOSE','-1'))
del _os
pkgload = PackageLoader()
pkgload(verbose=SCIPY_IMPORT_VERBOSE,postpone=True)

from numpy.testing import Tester
test = Tester().test
bench = Tester().bench
