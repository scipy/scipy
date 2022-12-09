:orphan:

.. _devpy-runtest:

===========================
Running SciPy Tests Locally
===========================

Basic test writing and execution from within the Python interpreter is
documented in the `NumPy/SciPy Testing Guidelines`_. This page includes
information about running tests from the command line using SciPy's
``dev.py`` command line tool. *Note: Before beginning, ensure that* |pytest|_
*is installed.*

To run all tests, navigate to the root SciPy directory at the command
line and execute

::

   python dev.py test -v

where ``-v`` activates the ``--verbose`` option. This builds SciPy (or
updates an existing build) and runs the tests.

To run tests on a particular submodule, such as ``optimize``, use the
``--submodule`` option:

::

   python dev.py test -v -s optimize

To run a particular test module, use the ``--test`` option:

::

   python dev.py test -v -t scipy.<module>.tests.<test_file>

Example for |test-linprog|_ file tests, run:

::

   python dev.py test -v -t scipy.optimize.tests.test_linprog

To run a test class:

::

   python dev.py test -v -t scipy.<module>.tests.<test_file>::<TestClass>

Example for ``TestLinprogRSCommon`` class from ``test_linprog.py``:

::

   python dev.py test -v -t scipy.optimize.tests.test_linprog::TestLinprogRSCommon

To run a particular test:

::

   python dev.py test -v -t scipy.<module>.tests.<test_file>::<test_name>

Example for ``test_unknown_solvers_and_options`` from ``test_linprog.py``:

::

   python dev.py test -v -t scipy.optimize.tests.test_linprog::test_unknown_solvers_and_options

For tests within a class, you need to specify the class name and the test name:

::

   python dev.py test -v -t scipy.<module>.tests.<test_file>::<TestClass>::<test_name>

Example:

::

   python dev.py test -v -t scipy.optimize.tests.test_linprog::TestLinprogRSCommon::test_nontrivial_problem_with_guess


Other useful options include:

-  ``--coverage`` to generate a test coverage report in
   ``scipy/build/coverage/index.html``. *Note:* |pytest-cov|_ *must be
   installed.*
-  ``--doc`` to build the docs in ``scipy/doc/build``. By default,
   docs are built only in the ``html`` format, but you can
   change this by appending the name of the desired format.
-  ``python dev.py refguide-check`` to check whether the objects in a SciPy
   submodule's ``__all__`` dict correspond to the objects included in the
   reference guide. It also checks the validity of code samples in docstrings.
   Run ``python dev.py refguide-check --help`` for more information.
-  ``--bench`` to run all benchmarks. See :ref:`benchmarking-with-asv`.
-  ``--pep8`` to perform pep8 check.
-  ``--mypy`` to run *mypy* on the codebase.
-  ``-n`` or ``--no-build`` to prevent SciPy from updating the build
   before testing
-  ``-j`` or ``--parallel`` *n* to engage *n* cores when building SciPy;
   e.g. \ ``python dev.py test -j 4`` engages four cores. As of `#10172`_
   this also runs the tests on four cores if |pytest-xdist|_ is installed.
-  ``-m`` or ``--mode`` ``full`` to run the full test suite, including slow
   tests. For example, ``python dev.py test -m full``.
-  ``--`` to send remaining command line arguments to ``pytest`` instead of
   ``dev.py test``. For instance, while ``-n`` sent to ``pytest.py`` activates
   the ``--no-build`` option, ``-n`` sent to ``pytest`` runs the tests on
   multiple cores; e.g. \ ``python dev.py test -- -n 4`` runs tests using
   four cores. *Note:* |pytest-xdist|_ *must be installed for testing on
   multiple cores.*

Other options not documented here are listed in the ``main`` function of
the source code for |runtests-py|_. For much more information about
``pytest``, see the ``pytest``
`documentation <https://docs.pytest.org/en/latest/usage.html>`_.

Tips:
-----

If you built SciPy from source but are having trouble running tests
after a change to the codebase, try deleting the ``scipy/build``
directory. This forces ``dev.py`` to completely rebuild SciPy before
performing tests.

There is an additional level of very slow tests (several minutes),
which are disabled even when calling ``python dev.py test -m full``.
They can be enabled by setting the environment variable ``SCIPY_XSLOW=1``
before running the test suite.

.. |runtests-py| replace:: ``runtests.py``
.. _runtests-py: https://github.com/scipy/scipy/blob/main/runtests.py

.. |pytest-cov| replace:: ``pytest-cov``
.. _pytest-cov: https://pypi.org/project/pytest-cov/

.. _#10172: https://github.com/scipy/scipy/pull/10172

.. |pytest-xdist| replace:: ``pytest-xdist``
.. _pytest-xdist: https://pypi.org/project/pytest-xdist/

.. _NumPy/SciPy Testing Guidelines: https://github.com/numpy/numpy/blob/main/doc/TESTS.rst

.. |pytest| replace:: ``pytest``
.. _pytest: https://docs.pytest.org/en/latest/

.. |test-linprog| replace:: ``scipy/optimize/tests/test_linprog.py``
.. _test-linprog: https://github.com/scipy/scipy/blob/main/scipy/optimize/tests/test_linprog.py
