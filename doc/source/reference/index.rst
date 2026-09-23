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

The remainder of this page contains the list of :ref:`submodules <submodule_list>`,
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


.. _main_namespace:

Main namespace (``scipy``)
=============================

.. automodule:: scipy
    :exclude-members: LowLevelCallable, show_config


.. _array-api-coverage:

Array API coverage
==================

Many SciPy functions supports multiple array types via the 
`Python array API standard <https://data-apis.org/array-api/latest/index.html>`_.

This standard allows users to use any array API compatible array library
with parts of SciPy out of the box, with with the main
principle being *"array type in equals array type out"*.
Currently, SciPy supports NumPy, PyTorch, JAX, Dask and CuPy arrays.

This clustering example shows usage with PyTorch tensors as inputs and return
values:

.. code:: python

    >>> import torch
    >>> from scipy.cluster.vq import vq
    >>> code_book = torch.tensor([[1., 1., 1.],
    ...                           [2., 2., 2.]])
    >>> features  = torch.tensor([[1.9, 2.3, 1.7],
    ...                           [1.5, 2.5, 2.2],
    ...                           [0.8, 0.6, 1.7]])
    >>> code, dist = vq(features, code_book)
    >>> code
    tensor([1, 1, 0], dtype=torch.int32)
    >>> dist
    tensor([0.4359, 0.7348, 0.8307])

The below tables show the current state of alternative backend support across
SciPy's modules. Public functions, function-like callables, and classes are
included in the tables. Parts of the public API which are deemed out-of-scope
are excluded from consideration when calculating coverage percentages. If a
module or submodule contains no in-scope functions, it is excluded from the
tables. For example, `scipy.datasets` is excluded because its contents are
considered out-of-scope.

.. toctree::
   :hidden:

   array_api_modules_tables/cluster_vq
   array_api_modules_tables/cluster_hierarchy
   array_api_modules_tables/constants
   array_api_modules_tables/differentiate
   array_api_modules_tables/fft
   array_api_modules_tables/integrate
   array_api_modules_tables/interpolate
   array_api_modules_tables/linalg
   array_api_modules_tables/linalg_interpolative
   array_api_modules_tables/ndimage
   array_api_modules_tables/optimize
   array_api_modules_tables/optimize_elementwise
   array_api_modules_tables/signal
   array_api_modules_tables/signal_windows
   array_api_modules_tables/sparse
   array_api_modules_tables/sparse_linalg
   array_api_modules_tables/sparse_csgraph
   array_api_modules_tables/spatial
   array_api_modules_tables/spatial_distance
   array_api_modules_tables/spatial_transform
   array_api_modules_tables/special
   array_api_modules_tables/stats
   array_api_modules_tables/stats_contingency
   array_api_modules_tables/stats_qmc

Support on CPU
--------------

.. array-api-support-per-module::
   :backend_type: cpu
   :cluster.vq: array_api_support_cluster_vq_cpu
   :cluster.hierarchy: array_api_support_cluster_hierarchy_cpu
   :constants: array_api_support_constants_cpu
   :differentiate: array_api_support_differentiate_cpu
   :fft: array_api_support_fft_cpu
   :integrate: array_api_support_integrate_cpu
   :interpolate: array_api_support_interpolate_cpu
   :linalg: array_api_support_linalg_cpu
   :linalg.interpolative: array_api_support_linalg_interpolative_cpu
   :ndimage: array_api_support_ndimage_cpu
   :optimize: array_api_support_optimize_cpu
   :optimize.elementwise: array_api_support_optimize_elementwise_cpu
   :signal: array_api_support_signal_cpu
   :signal.windows: array_api_support_signal_windows_cpu
   :sparse: array_api_support_sparse_cpu
   :sparse.linalg: array_api_support_sparse_linalg_cpu
   :sparse.csgraph: array_api_support_sparse_csgraph_cpu
   :spatial: array_api_support_spatial_cpu
   :spatial.distance: array_api_support_spatial_distance_cpu
   :spatial.transform: array_api_support_spatial_transform_cpu
   :special: array_api_support_special_cpu
   :stats: array_api_support_stats_cpu
   :stats.contingency: array_api_support_stats_contingency_cpu
   :stats.qmc: array_api_support_stats_qmc_cpu

Support on GPU
--------------

.. array-api-support-per-module::
   :backend_type: gpu
   :cluster.vq: array_api_support_cluster_vq_gpu
   :cluster.hierarchy: array_api_support_cluster_hierarchy_gpu
   :constants: array_api_support_constants_gpu
   :differentiate: array_api_support_differentiate_gpu
   :fft: array_api_support_fft_gpu
   :integrate: array_api_support_integrate_gpu
   :interpolate: array_api_support_interpolate_gpu
   :linalg: array_api_support_linalg_gpu
   :linalg.interpolative: array_api_support_linalg_interpolative_gpu
   :ndimage: array_api_support_ndimage_gpu
   :optimize: array_api_support_optimize_gpu
   :optimize.elementwise: array_api_support_optimize_elementwise_gpu
   :signal: array_api_support_signal_gpu
   :signal.windows: array_api_support_signal_windows_gpu
   :sparse: array_api_support_sparse_gpu
   :sparse.linalg: array_api_support_sparse_linalg_gpu
   :sparse.csgraph: array_api_support_sparse_csgraph_gpu
   :spatial: array_api_support_spatial_gpu
   :spatial.distance: array_api_support_spatial_distance_gpu
   :spatial.transform: array_api_support_spatial_transform_gpu
   :special: array_api_support_special_gpu
   :stats: array_api_support_stats_gpu
   :stats.contingency: array_api_support_stats_contingency_gpu
   :stats.qmc: array_api_support_stats_qmc_gpu

Support with JIT
----------------

.. array-api-support-per-module::
   :backend_type: jit
   :cluster.vq: array_api_support_cluster_vq_jit
   :cluster.hierarchy: array_api_support_cluster_hierarchy_jit
   :constants: array_api_support_constants_jit
   :differentiate: array_api_support_differentiate_jit
   :fft: array_api_support_fft_jit
   :integrate: array_api_support_integrate_jit
   :interpolate: array_api_support_interpolate_jit
   :linalg: array_api_support_linalg_jit
   :linalg.interpolative: array_api_support_linalg_interpolative_jit
   :ndimage: array_api_support_ndimage_jit
   :optimize: array_api_support_optimize_jit
   :optimize.elementwise: array_api_support_optimize_elementwise_jit
   :signal: array_api_support_signal_jit
   :signal.windows: array_api_support_signal_windows_jit
   :sparse: array_api_support_sparse_jit
   :sparse.linalg: array_api_support_sparse_linalg_jit
   :sparse.csgraph: array_api_support_sparse_csgraph_jit
   :spatial: array_api_support_spatial_jit
   :spatial.distance: array_api_support_spatial_distance_jit
   :spatial.transform: array_api_support_spatial_transform_jit
   :special: array_api_support_special_jit
   :stats: array_api_support_stats_jit
   :stats.contingency: array_api_support_stats_contingency_jit
   :stats.qmc: array_api_support_stats_qmc_jit
