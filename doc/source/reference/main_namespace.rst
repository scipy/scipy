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
`optimize`     Numerical optimization
`signal`       Signal processing
`sparse`       Sparse arrays, linear algebra and graph algorithms
`spatial`      Spatial data structures and algorithms
`special`      Special functions
`stats`        Statistical functions
=============  ==================================================



.. _arrayapi_coverage:

Array API Coverage
------------------


Many SciPy functions supports multiple array types via the 
`Python array API standard <https://data-apis.org/array-api/latest/index.html>`_.

This standard allows users to use any array API compatible array library
with parts of SciPy out of the box, with with the main
principle being *"array type in equals array type out"*.
Currently, SciPy supports NumPy, PyTorch, JAX, Dask and CuPy arrays.

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
   array_api_modules_tables/io
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
``````````````

.. array-api-support-per-module::
   :backend_type: cpu
   :cluster.vq: array_api_support_cluster_vq_cpu
   :cluster.hierarchy: array_api_support_cluster_hierarchy_cpu
   :constants: array_api_support_constants_cpu
   :differentiate: array_api_support_differentiate_cpu
   :fft: array_api_support_fft_cpu
   :integrate: array_api_support_integrate_cpu
   :interpolate: array_api_support_interpolate_cpu
   :io: array_api_support_io_cpu
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
``````````````

.. array-api-support-per-module::
   :backend_type: gpu
   :cluster.vq: array_api_support_cluster_vq_gpu
   :cluster.hierarchy: array_api_support_cluster_hierarchy_gpu
   :constants: array_api_support_constants_gpu
   :differentiate: array_api_support_differentiate_gpu
   :fft: array_api_support_fft_gpu
   :integrate: array_api_support_integrate_gpu
   :interpolate: array_api_support_interpolate_gpu
   :io: array_api_support_io_gpu
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
````````````````

.. array-api-support-per-module::
   :backend_type: jit
   :cluster.vq: array_api_support_cluster_vq_jit
   :cluster.hierarchy: array_api_support_cluster_hierarchy_jit
   :constants: array_api_support_constants_jit
   :differentiate: array_api_support_differentiate_jit
   :fft: array_api_support_fft_jit
   :integrate: array_api_support_integrate_jit
   :interpolate: array_api_support_interpolate_jit
   :io: array_api_support_io_jit
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

