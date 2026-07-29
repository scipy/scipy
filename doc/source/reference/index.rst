.. _scipy-api:

*************
API Reference
*************

.. A `toctree` should always be placed into the file's top-level section.
   Otherwise, Sphinx will become confused with section hierarchies.

.. toctree::
   :caption: API Reference
   :hidden:

   scipy.cluster <cluster>
   scipy.constants <constants>
   scipy.datasets <datasets>
   scipy.differentiate <differentiate>
   scipy.fft <fft>
   scipy.fftpack <fftpack>
   scipy.integrate <integrate>
   scipy.interpolate <interpolate>
   scipy.io <io>
   scipy.linalg <linalg>
   scipy.ndimage <ndimage>
   scipy.optimize <optimize>
   scipy.signal <signal>
   scipy.sparse <sparse>
   scipy.spatial <spatial>
   scipy.special <special>
   scipy.stats <stats>

SciPy's functionality is organized into :ref:`submodules <submodule_list>`, whereas
the :ref:`main namespace <main_namespace>` (``scipy``) only has a few utility functions.
In SciPy, most functions and classes are self-contained and are straightforward to use,
e.g.:

>>> from scipy.constants import speed_of_light
>>> from scipy.signal.windows import hann
...
>>> print(f"{speed_of_light} m/s is quite fast.")
299792458.0 m/s is quite fast.
>>> hann(7, sym=True)  # 7 sample symmetric Hann window
array([0.  , 0.25, 0.75, 1.  , 0.75, 0.25, 0.  ])

The remainder of this page is contains the list of :ref:`submodules <submodule_list>`,
followed by the content of the :ref:`main namespace <main_namespace>`. The :ref:`design
conventions for SciPy submodules <design_conventions_modules>` are part of the
:ref:`scipy-development`.


.. _submodule_list:

Submodules
==========
The public submodules have the following structure:

+---+---+-------------------------------------------------------------------+
| :ref:`scipy <main_namespace>` Main namespace                              |
+---+---+-------------------------------------------------------------------+
|   | `scipy.cluster` Clustering algorithms                                 |
+---+---+-------------------------------------------------------------------+
|   |   | `scipy.cluster.hierarchy` Hierarchical clustering                 |
+---+---+-------------------------------------------------------------------+
|   |   | `scipy.cluster.vq` K-means clustering and vector quantization     |
+---+---+-------------------------------------------------------------------+
|   | `scipy.constants` Physical and mathematical constants                 |
+---+---+-------------------------------------------------------------------+
|   | `scipy.datasets` Datasets                                             |
+---+---+-------------------------------------------------------------------+
|   | `scipy.differentiate` Finite Difference Differentiation               |
+---+---+-------------------------------------------------------------------+
|   | `scipy.fft` Discrete Fourier transforms                               |
+---+---+-------------------------------------------------------------------+
|   | `scipy.fftpack` Legacy discrete Fourier transforms                    |
+---+---+-------------------------------------------------------------------+
|   | `scipy.integrate` Integration and ODEs                                |
+---+---+-------------------------------------------------------------------+
|   | `scipy.interpolate` Interpolation                                     |
+---+---+-------------------------------------------------------------------+
|   | `scipy.io` Input and output                                           |
+---+---+-------------------------------------------------------------------+
|   |   | `scipy.io.arff` ARFF files                                        |
+---+---+-------------------------------------------------------------------+
|   |   | `scipy.io.matlab` MATLAB® files                                   |
+---+---+-------------------------------------------------------------------+
|   |   | `scipy.io.wavfile` WAV sound files                                |
+---+---+-------------------------------------------------------------------+
|   | `scipy.linalg` Linear algebra                                         |
+---+---+-------------------------------------------------------------------+
|   |   | `scipy.linalg.blas` Low-level BLAS functions                      |
+---+---+-------------------------------------------------------------------+
|   |   | `scipy.linalg.cython_blas` BLAS Functions for Cython              |
+---+---+-------------------------------------------------------------------+
|   |   | `scipy.linalg.interpolative` Interpolative matrix decomposition   |
+---+---+-------------------------------------------------------------------+
|   |   | `scipy.linalg.lapack` Low-level LAPACK functions                  |
+---+---+-------------------------------------------------------------------+
|   |   | `scipy.linalg.cython_lapack` LAPACK functions for Cython          |
+---+---+-------------------------------------------------------------------+
|   | `scipy.ndimage` Multidimensional image processing                     |
+---+---+-------------------------------------------------------------------+
|   | `scipy.optimize` Optimization and root finding                        |
+---+---+-------------------------------------------------------------------+
|   |   | `scipy.optimize.cython_optimize` Cython optimize root finding API |
+---+---+-------------------------------------------------------------------+
|   |   | `scipy.optimize.elementwise` Elementwise Scalar Optimization      |
+---+---+-------------------------------------------------------------------+
|   | `scipy.signal` Signal processing                                      |
+---+---+-------------------------------------------------------------------+
|   |   | `scipy.signal.windows` Window functions                           |
+---+---+-------------------------------------------------------------------+
|   | `scipy.sparse` Sparse linear algebra                                  |
+---+---+-------------------------------------------------------------------+
|   |   | `scipy.sparse.csgraph` Compressed sparse graph routines           |
+---+---+-------------------------------------------------------------------+
|   |   | `scipy.sparse.linalg` Sparse linear algebra                       |
+---+---+-------------------------------------------------------------------+
|   | `scipy.spatial` Spatial algorithms and data structures                |
+---+---+-------------------------------------------------------------------+
|   |   | `scipy.spatial.distance` Distance computations                    |
+---+---+-------------------------------------------------------------------+
|   |   | `scipy.spatial.transform` Spatial Transformations                 |
+---+---+-------------------------------------------------------------------+
|   | `scipy.special` Special functions                                     |
+---+---+-------------------------------------------------------------------+
|   | `scipy.stats` Statistical functions                                   |
+---+---+-------------------------------------------------------------------+
|   |   | `scipy.stats.contingency` Contingency table functions             |
+---+---+-------------------------------------------------------------------+
|   |   | `scipy.stats.mstats` Statistical functions for masked arrays      |
+---+---+-------------------------------------------------------------------+
|   |   | `scipy.stats.qmc` Quasi-Monte Carlo submodule                     |
+---+---+-------------------------------------------------------------------+
|   |   | `scipy.stats.sampling` Random Number Generators                   |
+---+---+-------------------------------------------------------------------+

The ``misc`` submodule is deprecated and does not contain functions anymore.


.. _main_namespace:

Main namespace (``scipy``)
=============================

.. automodule:: scipy
    :exclude-members: LowLevelCallable, show_config
