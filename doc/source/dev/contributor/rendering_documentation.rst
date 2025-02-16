.. _contributing-docs:

=======================================
Contributing to the SciPy documentation
=======================================

We're eager to hear about and fix doc defects. But to attack the biggest
problems we end up having to defer or overlook some bug reports. Here are the
best defects to go after.

Top priority goes to **technical inaccuracies** – a docstring missing a
parameter, a faulty description of a function/parameter/method, and so on. Other
“structural” defects like broken links also get priority. All these fixes are
easy to confirm and put in place. You can
:ref:`submit a pull request (PR)<editing-workflow>` with the fix, if you know
how to do that; otherwise please
`open an issue <https://github.com/scipy/scipy/issues/new/choose>`__.

**Typos and misspellings** fall on a lower rung; we welcome hearing about them
but may not be able to fix them promptly. These too can be handled as pull
requests or issues.

Obvious **wording** mistakes (like leaving out a “not”) fall into the typo
category, but other rewordings – even for grammar – require a judgment call,
which raises the bar. One can imagine cases where it is not clear that a "fix"
should be made, e.g.:

* Attempting to "fix" all the uses (or lack thereof) of the Oxford comma.
* Cases where the correctness of common usage is evolving, e.g. "comprised of"

The easiest fixes to accept are those where the original version is clearly and
unambiguously wrong; changes that require subtle editorial judgement are
probably best avoided. (But note that this is not about updating documentation
to fix confusing statements or otherwise deal with documentation problems
reported by users.)

.. note::

   As a general guideline, try to accumulate small documentation changes (such
   as typos) instead of sending them one by one. Whenever possible, also make
   sure to use the correct commands to :ref:`skip CI checks <skip-ci>` on
   documentation changes.

Some functions/objects defined in C or Fortran extension modules have their
docstrings defined separately from the actual code. Make sure to do a search for
the function docstring you are looking for using either `grep` or other similar
tools.

.. _rendering-documentation:

Rendering documentation locally with Sphinx
-------------------------------------------

SciPy docstrings are rendered to HTML using `Sphinx`_ and the
`PyData Sphinx theme`_. Writing
docstrings is covered in the :ref:`numpy:howto-document`; this document
explains how to check that docstrings render properly.

*For a video walkthrough, please see* \ `Rendering SciPy Documentation
with Sphinx`_ \ *.*

To render the documentation on your own machine:

0. Ensure that you have a working SciPy build (see :ref:`building-from-source`).
#. Then run ``python dev.py doc`` to build the documentation.
   This can take a while the first time, but subsequent documentation builds
   are typically much faster.
#. View the documentation in ``doc/build/html/``. You can start
   with ``index.html`` and browse, or you can jump straight to the file you’re
   interested in.

.. note::

   - Changes to certain documents do not take effect when Sphinx documentation
     is rebuilt. In this case, you can build from scratch by deleting the
     directories ``scipy/doc/build`` and ``source/reference/generated``, or by
     running ``python dev.py doc clean`` then building the docs again.

   - In case the SciPy version found by the above command is different from
     that of the latest commit in the repo, you will see a message like::

         installed scipy 5fd20ec1aa != current repo git version '35fd20ec1a'

     This indicates that you're likely picking up the wrong SciPy install,
     check with ``python -c "import scipy; print(scipy.__file__)"``.

   - The SciPy documentation contains interactive examples rendered with the
     help of the ``jupyterlite-sphinx`` extension. For more details, see
     :ref:`adding-notebooks` and :ref:`interactive-docs`.

.. _rendering-documentation-cloud:

Checking Documentation on the Cloud
-----------------------------------

Once a PR is opened, you can check that documentation renders correctly
on the cloud.

#. Log in to `GitHub`_.
#. Log in `CircleCI`_ using your GitHub account.
#. Back in GitHub, at the bottom of the PR, select “Show all Checks”.
#. Next to “Check the rendered docs here!”, select “Details”.

.. _docs-guidelines:

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
                `x` must be nonnegative.

    No::

            Parameters
            ----------
            x : float
                `x` should be nonnegative.


Use of the 'versionadded' markup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* For a new function, the 'versionadded' markup goes in the "Notes" section,
  *not* in the description at the beginning of the docstring.
* For a new argument added to an existing function,  the 'versionadded' markup
  is placed at the end of the description of the argument in the "Parameters"
  section.


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


DOIs in references
~~~~~~~~~~~~~~~~~~
The use of DOIs in references is strongly recommended.
There is special Sphinx syntax for DOIs: ``:doi:``. For example::

    .. [2] D. Fishkind, S. Adali, H. Patsolic, L. Meng, D. Singh, V. Lyzinski,
           C. Priebe, "Seeded graph matching", Pattern Recognit. 87 (2019):
           203-215, :doi:`10.1016/j.patcog.2018.09.014`

(arXiv articles also have special markup available: ``:arxiv:``.)


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


Self-contained examples
~~~~~~~~~~~~~~~~~~~~~~~
Each "Example" section (both in docstrings and general documentation)
must be self-contained. This means that all imports
must be explicit, the data used must be defined, and the code should "just
work" when copy-pasted into a fresh Python interpreter.

    Yes::

        >>> import numpy as np
        >>> rng = np.random.default_rng()

    No::

        >>> rng = np.random.default_rng()

What is possible (and recommended) is to intersperse blocks of code with
explanations. Blank lines must separate each code block from the explanatory
text.

    Yes::

        Some initial text

        >>> import numpy as np
        >>> rng = np.random.default_rng()

        This is some explanation

        >>> rng.random(10)


Examples and randomness
~~~~~~~~~~~~~~~~~~~~~~~
In the continuous integration (CI) suite, examples are executed and the output
is compared against the provided reference. The main goal is to ensure that
the *example* is correct; a failure warns us that the example may need to be
adjusted (e.g. because the API has changed since it was written).
Doctests are not meant to be used as unit tests of underlying implementation.

In case a random number generator is needed, `np.random.Generator` must be
used. The canonical way to create a NumPy ``Generator`` is to use
`np.random.default_rng`.

    Yes::

        >>> import numpy as np
        >>> rng = np.random.default_rng()
        >>> sample = rng.random(10)

    Yes::

        >>> import numpy as np
        >>> rng = np.random.default_rng(102524723947864966825913730119128190984)
        >>> sample = rng.random(10)

    No::

        >>> import numpy as np
        >>> sample = np.random.random(10)

Seeding the generator object is optional. If a seed is used, avoid common numbers and
instead generate a seed with ``np.random.SeedSequence().entropy``.
If no seed is provided, the default value
``1638083107694713882823079058616272161``
is used when doctests are executed. In either case, the rendered
documentation will not show the seed. The intent is to discourage users from
copy/pasting seeds in their code and instead make an explicit decision about
the use of a seed in their program. The consequence is that users cannot
reproduce the results of the example exactly, so examples using random data
should not refer to precise numerical values based on random data or rely on
them to make their point.

Legacy directive
~~~~~~~~~~~~~~~~

If a function, module or API is in *legacy* mode, meaning that it is kept around
for backwards compatibility reasons, but is not recommended to use in new code,
you can use the ``.. legacy::`` directive.

By default, if used with no arguments, the legacy directive will generate the
following output:

.. legacy::


We strongly recommend that you also add a custom message, such as a new API to
replace the old one. This message will be appended to the default message::

   .. legacy::

      New code should use :mod:`scipy.fft`.

will create the following output:

.. legacy::

   New code should use :mod:`scipy.fft`.

Finally, if you want to mention a function, method (or any custom object)
instead of a *submodule*, you can use an optional argument::

    .. legacy:: function

This will create the following output:

.. legacy:: function

---

.. _GitHub: https://github.com/
.. _CircleCI: https://circleci.com/vcs-authorize/
.. _Sphinx: https://www.sphinx-doc.org/en/master/
.. _PyData Sphinx theme: https://pydata-sphinx-theme.readthedocs.io/en/latest/
.. _Sphinx-Design: https://sphinx-design.readthedocs.io
.. _numpydoc: https://numpydoc.readthedocs.io
.. _matplotlib: https://www.matplotlib.org/
.. _Rendering SciPy Documentation with Sphinx: https://youtu.be/kGSYU39EhJQ
.. _git submodules: https://git-scm.com/book/en/v2/Git-Tools-Submodules
.. _Make build automation tool: https://en.wikipedia.org/wiki/Make_(software)
