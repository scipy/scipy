..  -*- rst -*-

================
SciPy benchmarks
================

Benchmarking SciPy with Airspeed Velocity.


Usage
-----

Airspeed Velocity manages building and Python virtualenvs by itself,
unless told otherwise. Some of the benchmarking features in
``runtests.py`` also tell ASV to use the SciPy compiled by
``runtests.py``. To run the benchmarks, you do not need to install a
development version of SciPy to your current Python environment.

Run a benchmark against currently checked out SciPy version (don't record the
result)::

    python runtests.py --bench sparse.Arithmetic

Compare change in benchmark results to another branch::

    python runtests.py --bench-compare master sparse.Arithmetic

Run benchmarks against the system-installed SciPy rather than rebuilding::

    python runtests.py -n --bench sparse.Arithmetic

Run ASV commands::

    cd benchmarks
    ./run.py run --skip-existing-commits --steps 10 ALL
    ./run.py publish
    ./run.py preview

The ``run.py`` script sets up some environment variables and does other minor
maintenance jobs for you. The benchmark suite is runnable directly using the
``asv`` command.

More on how to use ``asv`` can be found in `ASV documentation`_
Command-line help is available as usual via ``asv --help`` and
``asv run --help``.

.. _ASV documentation: https://asv.readthedocs.io/


Writing benchmarks
------------------

See `ASV documentation`_ for basics on how to write benchmarks.

Some things to consider:

- When importing things from SciPy on the top of the test files, do it as::

      try:
          from scipy.sparse.linalg import onenormest
      except ImportError:
          pass

  The benchmark files need to be importable also when benchmarking old versions
  of SciPy. The benchmarks themselves don't need any guarding against missing
  features --- only the top-level imports.

- Try to keep the runtime of the benchmark reasonable.

- Use ASV's ``time_`` methods for benchmarking times rather than cooking up
  time measurements via ``time.clock``, even if it requires some juggling when
  writing the benchmark.

- Preparing arrays etc. should generally be put in the ``setup`` method rather
  than the ``time_`` methods, to avoid counting preparation time together with
  the time of the benchmarked operation.

- Use ``run_monitored`` from ``common.py`` if you need to measure memory usage.

- Benchmark versioning: by default ``asv`` invalidates old results
  when there is any code change in the benchmark routine or in
  setup/setup_cache.

  This can be controlled manually by setting a fixed benchmark version
  number, using the ``version`` attribute. See `ASV documentation`_
  for details.

  If set manually, the value needs to be changed manually when old
  results should be invalidated. In case you want to preserve previous
  benchmark results when the benchmark did not previously have a
  manual ``version`` attribute, the automatically computed default
  values can be found in ``results/benchmark.json``.
