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


Return Objects
~~~~~~~~~~~~~~
For new functions or methods that return two or more conceptually distinct
elements, return the elements in an object type that is not iterable. In
particular, do not return a ``tuple``, ``namedtuple``, or a "bunch" produced
by ``scipy._lib._bunch.make_tuple_bunch``, the latter being reserved for adding
new attributes to iterables returned by existing functions. Instead, use an
existing return class (e.g. `~scipy.optimize.OptimizeResult`), a new, custom
return class.

This practice of returning non-iterable objects forces callers to be more
explicit about the element of the returned object that they wish to access,
and it makes it easier to extend the function or method in a backward
compatible way.

If the return class is simple and not public (i.e. importable from a public
module), it may be documented like::

    Returns
    -------
    res : MyResultObject
        An object with attributes:

        attribute1 : ndarray
            Customized description of attribute 1.
        attribute2 : ndarray
            Customized description of attribute 2.

Here "MyResultObject" above does not link to external documentation because it
is simple enough to fully document all attributes immediately below its name.

Some return classes are sufficiently complex to deserve their own rendered
documentation. This is fairly standard if the return class is public, but
return classes should only be public if 1) they are intended to be imported by
end-users and 2) if they have been approved by the forum. For complex,
private return classes, please see  how `~scipy.stats.binomtest` summarizes
`~scipy.stats._result_classes.BinomTestResult` and links to its documentation,
and note that ``BinomTestResult`` cannot be imported from `~scipy.stats`.

Depending on the complexity of "MyResultObject", a normal class or a dataclass
can be used. When using dataclasses, do not use ``dataclasses.make_dataclass``,
instead use a proper declaration. This allows autocompletion to list all
the attributes of the result object and improves static analysis.
Finally, hide private attributes if any::

    @dataclass
    class MyResultObject:
        statistic: np.ndarray
        pvalue: np.ndarray
        confidence_interval: ConfidenceInterval
        _rho: np.ndarray = field(repr=False)


Test functions from `numpy.testing`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In new code, don't use `assert_almost_equal`, `assert_approx_equal` or
`assert_array_almost_equal`. This is from the docstrings of these functions::

    It is recommended to use one of `assert_allclose`,
    `assert_array_almost_equal_nulp` or `assert_array_max_ulp`
    instead of this function for more consistent floating point
    comparisons.

For more information about writing unit tests, see the
:doc:`NumPy Testing Guidelines <numpy:reference/testing>`.


Testing expected exceptions/ warnings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When writing a new test that a function call raises an exception or emits a
warning, the preferred style is to use ``pytest.raises``/``pytest.warns`` as
a context manager, with the code that is supposed to raise the exception in
the code block defined by the context manager. The ``match`` keyword argument
is given with enough of the expected message attached to the exception/warning
to distinguish it from other exceptions/warnings of the same class. Do not use
``np.testing.assert_raises`` or ``np.testing.assert_warns``, as they do not
support a ``match`` parameter.

For example, the function `scipy.stats.zmap` is supposed to raise a
``ValueError`` if the input contains ``nan`` and ``nan_policy`` is ``"raise"``.
A test for this is::

    scores = np.array([1, 2, 3])
    compare = np.array([-8, -3, 2, 7, 12, np.nan])
    with pytest.raises(ValueError, match='input contains nan'):
        stats.zmap(scores, compare, nan_policy='raise')

The ``match`` argument ensures that the test doesn't pass by raising
a ``ValueError`` that is not related to the input containing ``nan``.
