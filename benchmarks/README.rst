..  -*- rst -*-

================
SciPy benchmarks
================

Benchmarking Scipy with Airspeed Velocity.


Usage
-----

Run a benchmark against currently checked out Scipy version (don't record the
result)::

    python runtests.py --bench sparse.Arithmetic

Compare change in benchmark results to another branch::

    python runtests.py --bench-compare master sparse.Arithmetic

Run ASV commands::

    ./run.py run --skip-existing-commits --steps 10 ALL
    ./run.py publish
    ./run.py preview

The ``run.py`` script sets up some environment variables and does other minor
maintenance jobs for you, including cloning the results repository of Scipy
benchmarks.

More on how to use ``asv`` can be found in `ASV documentation`_

.. _ASV documentation: https://spacetelescope.github.io/asv/


Writing benchmarks
------------------

See `ASV documentation`_ for basics on how to write benchmarks.

Some things to consider:

- When importing things from Scipy on the top of the test files, do it as::

      try:
          from scipy.sparse.linalg import onenormest
      except ImportError:
          pass

  The benchmark files need to be importable also when benchmarking old versions
  of Scipy. The benchmarks themselves don't need any guarding against missing
  features --- only the top-level imports.

- Try to keep the runtime of the benchmark reasonable.

- Use ASV's ``time_`` methods for benchmarking times rather than cooking up
  time measurements via ``time.clock``, even if it requires some juggling when
  writing the benchmark.

- Preparing arrays etc. should generally be put in the ``setup`` method rather
  than the ``time_`` methods, to avoid counting preparation time together with
  the time of the benchmarked operation.

- Use ``run_monitored`` from ``common.py`` if you need to measure memory usage.

