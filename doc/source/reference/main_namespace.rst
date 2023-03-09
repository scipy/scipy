========================
The main SciPy namespace
========================

.. currentmodule:: scipy

The main ``scipy`` namespace has very few objects in it by design. Only show
generical functionality related to testing, build info and versioning, and one
class (`LowLevelCallable`) that didn't fit in one of the
submodules, are present:

.. autosummary::
   :toctree: generated/

   LowLevelCallable
   show_config
   test

The one public attribute is:

================== ===============================================
``__version__``    SciPy version string
================== ===============================================


Submodules
----------

=============  ==================================================
`cluster`      Clustering functionality
`constants`    Physical and mathematical constants and units
`datasets`     Load SciPy datasets
`fft`          Discrete Fourier and related transforms
`fftpack`      Discrete Fourier transforms (legacy)
`integrate`    Numerical integration and ODEs
`interpolate`  Interpolation
`io`           Scientific data format reading and writing
`linalg`       Linear algebra functionality
`misc`         Utility routines (deprecated)
`ndimage`      N-dimensional image processing and interpolation
`odr`          Orthogonal distance regression
`optimize`     Numerical optimization
`signal`       Signal processing
`sparse`       Sparse arrays, linear algebra and graph algorithms
`spatial`      Spatial data structures and algorithms
`special`      Special functions
`stats`        Statistical functions
=============  ==================================================
