.. _continuous-integration:

======================
Continuous Integration
======================

Continuous integration (CI) is part of our development process and ensure that
every piece of code or documentation which is contributed to SciPy is working
and does not have unforeseen effects.

.. note:: Before submitting or updating your PR, please ensure that you tested
          your changes locally. See :ref:`pr-checklist` and :ref:`devpy-test`.

Workflows
=========

We run more than 20 different workflows with different versions of the
dependencies, different architectures, etc. A PR must pass all these checks
before it can be merged as to ensure a sustainable state of the project.

Apart from the unit tests, the documentation and examples in the docstrings are
also checked. These are common failing workflows as Sphinx and doctests have
very strict rules. These aspects are very important as documentation and
examples are user facing elements. Ensures that these elements are properly
rendered.

The logs can be long, but you will always find out why your build/test did not
pass a check. Simply click on ``Details`` to access the logs.

Following is a list of all the different workflows in use. They are grouped
by CI resources providers.

GitHub Actions
--------------
* ``Lint``: PEP8 and code style
* ``Windows Tests``: test suite runs for Windows
* ``Linux Tests``: test suite runs for Linux
* ``macOS Tests``: test suite runs for macOS (``x86_64``)
* ``Wheels builder``: builds wheels for SciPy releases as well as *nightly* builds.
* ``Check the rendered docs here!``: live preview of the documentation
* ``prerelease_deps_coverage_64bit_blas``: use pre-released version of the
  dependencies and check coverage
* ``gcc-9``: build with minimal supported version of GCC, install the wheel,
  then run the test suite with `python -OO`
* ``Array API``: test Array API support

The test suite runs on GitHub Actions and other platforms cover a range of
test/environment conditions: Python and NumPy versions
(lowest-supported to nightly builds), 32-bit vs. 64-bit, different compilers,
and more - for details, see the ``.yml`` configuration files.

CircleCI
--------
* ``build_docs``: build the documentation
* ``build_scipy``
* ``run_benchmarks``: verify how the changes impact performance
* ``refguide_check``: doctests from examples and benchmarks

.. _skip-ci:

Skipping
========

Being an open-source project, we have access to a quota of CI resources.
Ultimately, resources are limited and we should use them with care. This is
why we ask you to verify your changes locally before pushing them.

Depending on the proposed change, you might want to skip part of the checks.
It will be at the discretion of a maintainer to re-run some tests before
integration.

Skipping CI can be achieved by adding a special text in the commit message:

* ``[skip actions]``: will skip GitHub Actions
* ``[skip circle]``: will skip CircleCI
* ``[docs only]``: will skip *all but* the CircleCI checks and the linter
* ``[lint only]``: will skip *all but* the linter
* ``[skip ci]``: will skip *all* CI

Of course, you can combine these to skip multiple workflows.

This skip information should be placed on a new line. In this example, we
just updated a ``.rst`` file in the documentation and ask to skip all but the
relevant docs checks (skip GitHub Actions' workflows)::

    DOC: improve QMCEngine examples.

    [docs only]

Failures due to test duration
=============================

Some CI jobs install |pytest-fail-slow|_ and report failures when the test
execution time exceeds a threshold duration.

- By default, all tests are subject to a 5 second limit; i.e., the option
  ``--fail-slow=5.0`` is used in a "full" test job.
- All tests not marked ``slow`` (``@pytest.mark.slow``) are subject to a
  1 second limit; i.e. the option ``--fail-slow=1.0`` is used in a "fast"
  test job.
- Exceptions are made using the ``pytest.mark.fail_slow`` decorator; e.g.
  a test can be marked ``@pytest.mark.fail_slow(10)`` to give it a ten
  second limit regardless of whether it is part of the "fast" or "full"
  test suite.

If a test fails by exceeding the time limit at any point during the
development of a PR, please adjust the test to ensure that it does
not fail in the future. Even if new tests do not fail, please check
the details of workflows that include "fail slow" in their name
before PRs merge. These include lists of tests that are approaching
(or have exceeded) their time limit. Due to variation in execution
times, tests with execution times near the threshold should be adjusted
to avoid failure even if their execution time were to increase by 50%;
typical tests should have much greater margin (at least 400%).
Adjustment options include:

- Making the test faster.
- Marking the test as ``slow``, if it is acceptable to run the test
  on a reduced set of platforms.
- Marking the test as ``xslow``, if it is acceptable to run the test
  only occasionally.
- Breaking up the test or parameterizing it, and possible marking
  parts of it as slow. Note that this does not reduce the total
  test duration, so other options are preferred.
- For truly critical tests that are unavoidably slow, add an exception
  using ``pytest.mark.fail_slow``.

See :ref:`devpy-test` for more information about working with slow tests
locally.

Wheel builds
============

Wheels for SciPy releases and
`*nightly* <https://anaconda.org/scientific-python-nightly-wheels/scipy>`_ builds are built
using cibuildwheel in a
`Github Action <https://github.com/scipy/scipy/blob/main/.github/workflows/wheels.yml>`_.
The Action runs:

* when the commit message contains the text ``[wheel build]``
* on a scheduled basis once a week
* when it is started manually.
* when there is a push to the repository with a GitHub reference starting with ``refs/tags/v`` (and not ending with ``dev0``)

The action does not run on forks of the main SciPy repository. The wheels that
are created are available as artifacts associated with a successful run of the
Action. When the Action runs on a schedule, or is manually started, the wheels
are uploaded to the
`*scientific-python-nightly-wheels* <https://anaconda.org/scientific-python-nightly-wheels/scipy>`_
repository.

It is not advised to use cibuildwheel to build scipy wheels on your own system
as it will automatically install gfortran compilers and various other
dependencies. Instead, one could use an isolated Docker container to build
Linux wheels.

.. |pytest-fail-slow| replace:: ``pytest-fail-slow``
.. _pytest-fail-slow: https://github.com/jwodder/pytest-fail-slow