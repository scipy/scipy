:orphan:

.. highlight:: console

.. _benchmarking-with-asv:

Benchmarking SciPy with Airspeed Velocity
=========================================

*This document introduces benchmarking, including reviewing SciPy
benchmark test results online, writing a benchmark test, and running it
locally. For a video run-through of writing a test and running it
locally, see* \ `Benchmarking SciPy`_\ *.*

As written in the `airspeed velocity (asv) documentation`_:

 airspeed velocity (asv) is a tool for benchmarking Python packages over their
 lifetime. Runtime, memory consumption and even custom-computed values
 may be tracked. The results are displayed in an interactive web frontend
 that requires only a basic static webserver to host.

To see what this means, take a look at `airspeed velocity of an unladen
scipy`_. Each plot summarizes the execution time of a particular test
over the commit history of the project; that is, as each commit is
merged, the benchmark test is run, its execution time is measured, and
the elapsed time is plotted. In addition to tracking the performance of
the code a commit is *intended* to affect, running *all* benchmarks for
each commit is helpful for identifying unintentional regressions:
significant increases in the execution time of one or more benchmark
tests. As SciPy is a web of interconnected code, the repercussions of a
small change may not be immediately obvious to a contributor, so this
benchmark suite makes it easier to detect regressions and identify the
commit that caused them. When you contribute a substantial new feature -
or notice a feature that doesn’t already have a benchmark test - please
consider writing benchmarks.

Writing benchmarks
------------------

*The* \ :ref:`Writing benchmarks <asv:writing-benchmarks>` \ *section of the
airspeed velocity documentation is the definitive guide to writing benchmarks.
Please see also the* \ `SciPy benchmarks readme`_\ *.*

To see how benchmarks are written, take a look at
|optimize-linprog-py|_. Each subclass of
``Benchmark`` defines a benchmark test. For example, the ``KleeMinty``
class defines a benchmark test based on the `Klee-Minty hypercube
problem`_, a diabolical test of the simplex algorithm for linear
programming. The class has four parts:

-  ``setup`` prepares the benchmark to run. The execution time of this
   function is *not* counted in the benchmark results, so this is a good
   place to set up all variables that define the problem. In the ``KleeMinty``
   example, this involves generating arrays ``c``, ``A_ub``, and ``b_ub``
   corresponding with a Klee-Minty hypercube in ``dims`` dimensions and
   storing them as instance variables.
-  ``time_klee_minty`` actually runs the benchmark test. This function
   executes after a ``KleeMinty`` object has been instantiated and
   ``setup`` has run, so it gets the arrays defining the problem from
   ``self``. Note that the prefix ``time`` in the function name
   indicates to ``asv`` that the execution time of this function *is* to
   be counted in the benchmark results.
-  ``params`` is a list of lists defining parameters of the test.
   Benchmarks are run for all possible combinations of these parameters.
   For example, the first time the benchmark is run, the first element
   of ``methods`` (``simplex``) is passed into ``setup`` and
   ``time_klee_minty`` as the first argument, ``meth``, and the first
   element of ``[3, 6, 9]`` (``3``) is passed into ``setup`` and
   ``time_klee_minty`` as the second argument, ``dims``. The next time
   the benchmark is run, ``setup`` and ``time_klee_minty`` are passed
   ``revised simplex`` and ``6`` as arguments, and so this continues
   until all combinations of parameters have been used.
-  ``param_names`` is a list of human-readable names for each element of
   the ``params`` list. These are used for presenting results.

Results of this benchmark over the past few years are available by
clicking on the `KleeMinty.time_klee_minty`_ link at `airspeed velocity
of an unladen scipy`_. Note that each trace of the plot corresponds with
a combination of benchmark parameters and environment settings
(e.g. Cython version), and that the visibility of the traces can be
toggled using the control panel on the left.

Running benchmarks locally
--------------------------

*Before beginning, ensure that* \ `airspeed velocity`_ \ *is
installed.*

After contributing new benchmarks, you should test them locally before
submitting a pull request.

To run all benchmarks, navigate to the root SciPy directory at the
command line and execute::

   python runtests.py --bench

where ``--bench`` activates the benchmark suite instead of the test
suite. This builds SciPy and runs the benchmarks. (*Note: this could
take a while. Benchmarks often take longer to run than unit tests, and
each benchmark is run multiple times to measure the distribution in
execution times.*)

To run benchmarks from a particular benchmark module, such as
``optimize_linprog.py``, simply append the filename without the
extension::

   python runtests.py --bench optimize_linprog

To run a benchmark defined in a class, such as ``KleeMinty`` from
``optimize_linprog.py``::

   python runtests.py --bench optimize_linprog.KleeMinty

To compare benchmark results between the active branch and another, such
as ``master``::

   python runtests.py --bench-compare master optimize_linprog.KleeMinty

All of the commands above display the results in plain text in the
console, and the results are not saved for comparison with future
commits. For greater control, a graphical view, and to have results
saved for future comparison, you can use |run-py|_
to run the benchmarks.

``run.py`` is a "wrapper" for the ``asv`` terminal command, the use of
which is described in `Using airspeed velocity`_. ``run.py`` simply sets
some environment variables for you, then sends all remaining command
line arguments to ``asv``.

To use it, navigate to ``scipy/benchmarks`` in the console and then
execute::

   python run.py run

Just like ``asv run`` from the asv documentation, this command runs the
whole benchmark suite and saves the results for comparison against
future commits.

To run only a single benchmark, such as ``KleeMinty`` from
``optimize_linprog.py``::

   python run.py run --bench optimize_linprog.KleeMinty

(Note that the ``--bench`` option here is sent to ``asv``; its meaning
is different than the ``--bench`` option for ``runtests.py``.)

One great feature of ``asv`` is that it can automatically run a
benchmark not just for the current commit, but for every commit in a
range. ``linprog`` ``method='interior-point'`` was merged into SciPy
with commit |7fa17f2369e0e5ad055b23cc1a5ee079f9e8ca32|_, so let’s
run the ``KleeMinty`` benchmark for 10 commits between then and now to
track its performance over time::

   python run.py run --bench optimize_linprog.KleeMinty --steps 10 7fa17f..

.. note::

   This will take a while, because SciPy has to be rebuilt for each
   commit! For more information about specifying ranges of commits, see
   the `git revisions documentation`_.

To "publish" the results (prepare them to be viewed) and "preview" them
in an interactive console::

   python run.py publish
   python run.py preview

ASV will report that it is running a server. Using any browser, you can
review the results by navigating to http://127.0.0.1:8080 (local
machine, port 8080).

For much more information about the ``asv`` commands accessible via
``run.py``, see the airspeed velocity `Commands`_ documentation. (Tip:
check out the ``asv find`` command and the ``--quick``,
``--skip-existing-commits``, and ``--profile`` options for ``asv run``.)

.. _git revisions documentation: https://git-scm.com/docs/gitrevisions#_specifying_ranges
.. _Commands: https://asv.readthedocs.io/en/stable/commands.html#commands
.. _airspeed velocity: https://github.com/airspeed-velocity/asv
.. _Using airspeed velocity: https://asv.readthedocs.io/en/stable/using.html#running-benchmarks
.. _Benchmarking SciPy: https://youtu.be/edLQ8KRpupQ
.. _airspeed velocity (asv) documentation: https://asv.readthedocs.io/en/stable/
.. _airspeed velocity of an unladen scipy: https://pv.github.io/scipy-bench/
.. _SciPy benchmarks readme: https://github.com/scipy/scipy/blob/master/benchmarks/README.rst
.. _Klee-Minty hypercube problem: https://en.wikipedia.org/wiki/Klee%E2%80%93Minty_cube
.. _KleeMinty.time_klee_minty: https://pv.github.io/scipy-bench/#optimize_linprog.KleeMinty.time_klee_minty

.. |optimize-linprog-py| replace:: ``scipy/benchmarks/benchmarks/optimize_linprog.py``
.. _optimize-linprog-py: https://github.com/scipy/scipy/blob/master/benchmarks/benchmarks/optimize_linprog.py

.. |run-py| replace:: ``scipy/benchmarks/run.py``
.. _run-py: https://github.com/scipy/scipy/blob/master/benchmarks/run.py

.. |7fa17f2369e0e5ad055b23cc1a5ee079f9e8ca32| replace:: ``7fa17f2369e0e5ad055b23cc1a5ee079f9e8ca32``
.. _7fa17f2369e0e5ad055b23cc1a5ee079f9e8ca32: https://github.com/scipy/scipy/commit/7fa17f2369e0e5ad055b23cc1a5ee079f9e8ca32
