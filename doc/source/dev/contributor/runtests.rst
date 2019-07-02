:orphan:

.. _runtests:

===========================
Running SciPy Tests Locally
===========================

Basic test writing and execution from within the Python interpreter is
documented in the `NumPy/SciPy Testing Guidelines`_. This page includes
information about running tests from the command line using SciPy’s
``runtests.py``, which permits greater control. *Note: Before beginning,
ensure that* |pytest|_ *is installed.*

To run all tests, navigate to the root SciPy directory at the command
line and execute

::

   python runtests.py -v

where ``-v`` activates the ``--verbose`` option. This builds SciPy (or
updates an existing build) and runs the tests.

To run tests on a particular submodule, such as ``optimize``, use the
``--submodule`` option:

::

   python runtests.py -v -s optimize

To run a particular test module, such as
|test-linprog|_, use the ``--test`` option:

::

   python runtests.py -v -t scipy/optimize/tests/test_linprog.py

To run a test class, such as ``TestLinprogIPDense`` from
``test_linprog.py``:

::

   python runtests.py -v -t scipy/optimize/tests/test_linprog.py::TestLinprogIPDense

To run a particular test, such as ``test_unknown_solver`` from
``test_linprog.py``:

::

   python runtests.py -v -t scipy/optimize/tests/test_linprog.py::test_unknown_solver

For tests within a class, you need to specify the class name and the
test name:

::

   python runtests.py -v -t scipy/optimize/tests/test_linprog.py::TestLinprogIPDense::test_nontrivial_problem

Other useful options include:

-  ``-c`` or ``--coverage`` to generate a test coverage report in
   ``scipy/build/coverage/index.html``. *Note:* |pytest-cov|_ *must be
   installed.*
-  ``-n`` or ``--no-build`` to prevent SciPy from updating the build
   before testing
-  ``-j`` or ``--parallel`` *n* to engage *n* cores when building SciPy;
   e.g. \ ``python runtests.py -j 4`` engages four cores. As of `#10172`_
   this also runs the tests on four cores if |pytest-xdist|_ is installed.
-  ``-m`` or ``--mode`` ``full`` to run the full test suite, including slow
   tests. For example, ``python runtests.py -m full``.
-  ``--`` to send remaining command line arguments to ``pytest`` instead of
   ``runtest.py``. For instance, while ``-n`` sent to ``pytest.py`` activates
   the ``--no-build`` option, ``-n`` sent to ``pytest`` runs the tests on
   multiple cores; e.g. \ ``python runtests.py -- -n 4`` runs tests using
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
directory. This forces ``runtest.py`` to completely rebuild SciPy before
performing tests.

There is an additional level of very slow tests (several minutes),
which are disabled even when calling ``python runtests.py -m full``.
They can be enabled by setting the environment variable ``SCIPY_XSLOW=1``
before running the test suite.

.. |runtests-py| replace:: ``runtests.py``
.. _runtests-py: https://github.com/scipy/scipy/blob/master/runtests.py

.. |pytest-cov| replace:: ``pytest-cov``
.. _pytest-cov: https://pypi.org/project/pytest-cov/

.. _#10172: https://github.com/scipy/scipy/pull/10172

.. |pytest-xdist| replace:: ``pytest-xdist``
.. _pytest-xdist: https://pypi.org/project/pytest-xdist/

.. _NumPy/SciPy Testing Guidelines: https://github.com/numpy/numpy/blob/master/doc/TESTS.rst.txt

.. |pytest| replace:: ``pytest``
.. _pytest: https://docs.pytest.org/en/latest/

.. |test-linprog| replace:: ``scipy/optimize/tests/test_linprog.py``
.. _test-linprog: https://github.com/scipy/scipy/blob/master/scipy/optimize/tests/test_linprog.py
