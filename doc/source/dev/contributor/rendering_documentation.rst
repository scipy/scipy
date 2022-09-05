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

0. Ensure that you have a working SciPy :ref:`dev-env` active.
   You need to be able to ``import scipy`` regardless of Python's working
   directory; the ``python setup.py develop`` and ``conda develop`` commands
   from the :ref:`quickstart <dev-env>` guides make this possible.
#. Install `Sphinx`_, `PyData Sphinx theme`_, `Sphinx-Design`_
   and `matplotlib`_. For example, if you're using the Anaconda distribution of
   Python, enter in a terminal window::
      
      conda install sphinx pydata-sphinx-theme sphinx-design matplotlib --channel conda-forge
      
   The list of requirements is in ``scipy/doc_requirements.txt``.
#. In a terminal window, browse to the ``scipy/doc`` directory. Note the
   presence of the file ``Makefile``.
#. Execute ``git submodule update --init``.
   Some of the documentation theme files are not distributed
   with the main ``scipy`` repository; this keeps them up to date using
   `git submodules`_.
#. Enter ``make html``. If you have multiple version of Python on
   your path, you can choose which version to use by appending
   ``PYTHON=python3.9`` to this command, where ``python3.9`` is to be
   replaced with the name of the Python you use for SciPy development.
   This uses the `Make build automation tool`_
   to execute the documentation build instructions from the ``Makefile``.
   This can take a while the first time, but subsequent documentation builds
   are typically much faster.
#. View the documentation in ``scipy/doc/build/html``. You can start
   with ``index.html`` and browse, or you can jump straight to the file you’re
   interested in.

.. note::

   Changes to certain documents do not take effect when Sphinx documentation
   is rebuilt. In this case, you can build from scratch by deleting the
   ``scipy/doc/build`` directory, then building again.

.. _rendering-documentation-cloud:

Checking Documentation on the Cloud
-----------------------------------

Once a PR is opened, you can check that documentation renders correctly
on the cloud.

#. Log in to `GitHub`_.
#. Log in `CircleCI`_ using your GitHub account.
#. Back in GitHub, at the bottom of the PR, select “Show all Checks”.
#. Next to “ci/circleci: build_docs artifact”, select “Details”.

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


.. _GitHub: https://github.com/
.. _CircleCI: https://circleci.com/vcs-authorize/
.. _Sphinx: https://www.sphinx-doc.org/en/master/
.. _PyData Sphinx theme: https://pydata-sphinx-theme.readthedocs.io/en/latest/
.. _Sphinx-Design: https://sphinx-design.readthedocs.io
.. _matplotlib: https://www.matplotlib.org/
.. _Rendering SciPy Documentation with Sphinx: https://youtu.be/kGSYU39EhJQ
.. _git submodules: https://git-scm.com/book/en/v2/Git-Tools-Submodules
.. _Make build automation tool: https://en.wikipedia.org/wiki/Make_(software)
