.. _scipy_parallel_execution:

Parallel execution support in SciPy
===================================

SciPy aims to provide functionality that is performant, i.e. has good execution
speed. On modern computing hardware, CPUs often have many CPU cores - and hence
users may benefit from parallel execution. This page aims to give a brief
overview of the options available to employ parallel execution.

Some key points related to parallelism:

- SciPy itself defaults to single-threaded execution.

- The exception to that single-threaded default is code that calls into a BLAS
  or LAPACK library for linear algebra functionality (either direct or via
  NumPy). BLAS/LAPACK libraries almost always default to multi-threaded execution,
  typically using all available CPU cores.

  - Users can control the threading behavior of the BLAS/LAPACK library that
    SciPy and NumPy are linked with through
    `threadpoolctl <https://github.com/joblib/threadpoolctl>`__.

- SciPy functionality may provide parallel execution in an opt-in manner. This
  is exposed through a ``workers=`` keyword in individual APIs, which takes an
  integer for the number of threads or processes to use, and in some cases also
  a map-like callable (e.g., ``multiprocessing.Pool``). See `scipy.fft.fft` and
  `scipy.optimize.differential_evolution` for examples.

  - SciPy-internal threading is done with OS-level thread pools. OpenMP is not
    used within SciPy.

- SciPy works well with `multiprocessing` and with `threading`. The former has
  higher overhead than the latter, but is widely used and robust. The latter may
  offer performance benefits for some usage scenarios - however, please read
  :ref:`scipy_thread_safety`.

- SciPy has *experimental* support for free-threaded CPython, starting with
  SciPy 1.15.0 (and Python 3.13.0, NumPy 2.1.0).

- SciPy has *experimental* support in a growing number of submodules and
  functions for array libraries other than NumPy, such as PyTorch, CuPy and
  JAX. Those libraries default to parallel execution and may offer significant
  performance benefits (and GPU execution). See :ref:`dev-arrayapi` for more
  details.
