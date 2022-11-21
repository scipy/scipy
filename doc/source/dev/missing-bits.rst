:orphan:

.. _missing-bits:

Code and Documentation Style Guide - The Missing Bits
=====================================================

This is a collection of coding and documentation guidelines for SciPy that
are not explicitly stated in the existing guidelines and standards, including

* `PEP-8 <https://www.python.org/dev/peps/pep-0008>`_ Style Guide for Python Code
* `PEP-257 <https://www.python.org/dev/peps/pep-0257>`_ Docstring Conventions
* `NumPy docstring standard
  <https://numpydoc.readthedocs.io/en/latest/format.html>`_
* :doc:`NumPy Testing Guidelines <numpy:reference/testing>`

Some of these are trivial, and might not seem worth discussing, but in many
cases, the issue has come up in a pull request review in either the SciPy
or NumPy repositories.  If a style issue is important enough that a reviewer
will require a change before merging, then it is important enough to be
documented--at least for cases where the issue can be resolved with a simple
rule.


Coding Style and Guidelines
---------------------------

Required keyword names
~~~~~~~~~~~~~~~~~~~~~~
For new functions or methods with more than a few arguments, all parameters
after the first few "obvious" ones should *require* the use of the keyword
when given.  This is implemented by including ``*`` at the appropriate point
in the signature.

For example, a function ``foo`` that operates on a single array but that has
several optional parameters (say ``method``, ``flag``, ``rtol`` and ``atol``)
would be defined as::

    def foo(x, *, method='basic', flag=False, rtol=1.5e-8, atol=1-12):
        ...

To call ``foo``, all parameters other than ``x`` must be given with an
explicit keyword, e.g. ``foo(arr, rtol=1e-12, method='better')``.

This forces callers to give explicit keyword parameters (which most users
would probably do anyway even without the use of ``*``), *and* it means
additional parameters can be added to the function anywhere after the
``*``; new parameters do not have to be added after the existing parameters.


Test functions from `numpy.testing`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In new code, don't use `assert_almost_equal`, `assert_approx_equal` or
`assert_array_almost_equal`. This is from the docstrings of these
functions::

    It is recommended to use one of `assert_allclose`,
    `assert_array_almost_equal_nulp` or `assert_array_max_ulp`
    instead of this function for more consistent floating point
    comparisons.

For more information about writing unit tests, see the
:doc:`NumPy Testing Guidelines <numpy:reference/testing>`.


Testing that expected exceptions are raised
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When writing a new test that a function call raises a certain exception,
the preferred style is to use ``pytest.raises`` as a context manager, with
the code that is supposed to raise the exception in the code block defined
by the context manager.  The `match` keyword argument of ``pytest.raises``
is given with enough of the expected error message attached to the exception
to ensure that the expected exception is raised.

For example, the function `scipy.stats.zmap` is supposed to raise a
``ValueError`` if the input contains ``nan`` and ``nan_policy`` is ``"raise"``.
A test for this is::

    scores = np.array([1, 2, 3])
    compare = np.array([-8, -3, 2, 7, 12, np.nan])
    with pytest.raises(ValueError, match='input contains nan'):
        stats.zmap(scores, compare, nan_policy='raise')

The ``match`` argument ensures that the test doesn't pass by raising
a ``ValueError`` that is not related to the input containing ``nan``.
