
Code and Documentation Style Guide - The Missing Bits
=====================================================

This is a collection of coding and documentation guidelines for SciPy that
are not explicitly stated in the existing guidelines and standards, including

* `PEP-8 <https://www.python.org/dev/peps/pep-0008>`_ Style Guide for Python Code
* `PEP-257 <https://www.python.org/dev/peps/pep-0257>`_ Docstring Conventions
* `NumPy docstring standard <https://numpydoc.readthedocs.io/en/latest/format.html>`_

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


A parameter called lambda
~~~~~~~~~~~~~~~~~~~~~~~~~
There are functions where convention dictates that (ideally) a parameter
should be called 'lambda'.  In Python, 'lambda' is a reserved word that
cannot be used as a variable name.  If a new function or class has a
parameter where a name that is similar to 'lambda' is really the best
choice, the preferred name to use is [TBD: ``lam``, ``lmb``, ``lamb``,
``lmbda``, ``lambda_``?].

There has been a lack of consistency in the past:

============================== =============
Function or class              Variable name
============================== =============
`scipy.stats.boxcox`           ``lmbda``
`scipy.stats.boxcox_llf`       ``lmb``
`scipy.stats.yeojohnson`       ``lmbda``
`scipy.stats.yeojohnson_llf`   ``lmb``
`scipy.stats.chi2_contingency` ``lambda_``
`scipy.stats.boltzmann`        ``lambda_``
`scipy.stats.planck`           ``lambda_``
`scipy.stats.tukeylambda`      ``lam``
`scipy.signal.cspline1d`       ``lamb``
`scipy.signal.qspline1d`       ``lamb``
`scipy.signal.cspline2d`       ``lambda`` [*]
`scipy.signal.qspline2d`       ``lambda`` [*]
`scipy.signal.spline_filter`   ``lmbda``
`scipy.special.tklmbda`        ``lmbda``
`scipy.special.boxcox`         ``lmbda``
`scipy.special.boxcox1p`       ``lmbda``
`scipy.special.inv_boxcox`     ``lmbda``
`scipy.special.inv_boxcox1p`   ``lmbda``
============================== =============

.. [*] This function accepts positional arguments only, but the docstring
       documents the argument as `lambda`.


Implicit string continuation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When an implicit continuation is used to create a long literal string, the
string is generally split at a space, and the space is included in the
second part of the string.

  Yes::

      long_string = ("This is a very long string that demonstrates"
                     " this guideline.")

  No::

      long_string = ("This is a very long string that demonstrates "
                     "this guideline.")


Documentation Guidelines
------------------------

Use of LaTeX
~~~~~~~~~~~~
[TBD.  In the past, we had an unwritten rule that LaTex markup should be
restricted to the Notes section.  We don't follow that rule these days.
Should we?]


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


Last line of docstring is a blank line
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Include a blank line at the end of multiline docstrings.  (This refers
to a final blank line *within* the docstring, not a blank line in the
code that follows the docstring.)

This is not part of [PEP-257](https://www.python.org/dev/peps/pep-0257/),
and the examples in that PEP do not include the blank line.  It is
also not explicitly stated in the 
[NumPy docstring standard](https://numpydoc.readthedocs.io/en/latest/format.html),
but the few complete examples shown there all include a blank line at
the end.
