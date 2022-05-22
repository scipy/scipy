
Code and Documentation Style Guide - The Missing Bits
=====================================================

This is a collection of coding and documentation guidelines for SciPy that
are not explicitly stated in the existing guidelines and standards, including

* `PEP-8 <https://www.python.org/dev/peps/pep-0008>`_ Style Guide for Python Code
* `PEP-257 <https://www.python.org/dev/peps/pep-0257>`_ Docstring Conventions
* `NumPy docstring standard
  <https://numpydoc.readthedocs.io/en/latest/format.html>`_
* NumPy `Testing Guidelines
  <https://numpy.org/doc/stable/reference/testing.html>`_

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

For more information about writing unit tests, see `Testing Guidelines
<https://numpy.org/doc/stable/reference/testing.html>`_.

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


Documentation Guidelines
------------------------

Use "must", not "should"
~~~~~~~~~~~~~~~~~~~~~~~~
When specifying a required condition on the input parameters, the
word "must" is preferable to "should".  For many English speakers,
"must" implies a stronger constraint than "should",  e.g. "I must
have oxygen to live" versus "I should exercise more".

    Yes::

            Parameters
            ----------
            x : float
                x must be nonnegative.

    No::

            Parameters
            ----------
            x : float
                x should be nonnegative.


Use of the 'versionadded' markup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* For a new function, the 'versionadded' markup goes in the "Notes" section,
  *not* in the description at the beginning of the docstring.
* For a new argument added to an existing function,  two locations have been
  used for the the 'versionadded' markup, [TBD: which is preferred?]:

  * At the end of the description of the argument in the "Parameters" section
  * In the "Notes" section.  In this case, the `versionadded` markup
    wouldn't be used.  Instead, the new addition is noted with a plain
    text comment.


Citing wikipedia articles in the "References" section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
It is acceptable to use wikipedia articles as references.
When creating the citation for the reference, include the article title,
the name "Wikipedia" (similar to how one gives a journal title), and the
URL.

    Yes::

        .. [1] "Zeta Distribution", Wikipedia,
               https://en.wikipedia.org/wiki/Zeta_distribution

    No::

        .. [1] https://en.wikipedia.org/wiki/Zeta_distribution    


Use of ``np`` in the "Examples" section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Do not include ``import numpy as np`` in the code given in the "Examples"
section.  The NumPy Docstring standard says "The examples may assume that
``import numpy as np`` is executed before the example code in numpy."  That
statement makes the import *optional*; this guideline says explicitly
that the import statement must not be included.


Bulleted lists
~~~~~~~~~~~~~~
This is not so much a guideline as it is a reminder of the Sphinx markup
for bulleted lists.  The incorrect use of indentation is common enough
that it is worthwhile mentioning it here.

When creating a bulleted list:

* Don't end the preceding line with `::`.
* Don't indent the bullets.
* Include a blank line before and after the list.

Some examples:

    Yes::

        Some text that precedes this interesting list:

        * The first item in the list.
        * The second item in the list.
        * You get the idea.

        Some text that follows the list.

    No::

        Some text that precedes this interesting list:

          * The first item in the list.
          * The second item in the list.
          * You get the idea.

        Some text that follows the list.

    No::

        Some text that precedes this interesting list:
        * The first item in the list.
        * The second item in the list.
        * You get the idea.
        Some text that follows the list.

