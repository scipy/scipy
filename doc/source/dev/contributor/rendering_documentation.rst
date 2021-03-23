:orphan:

.. _rendering-documentation:

===================================
Rendering Documentation with Sphinx
===================================

SciPy docstrings are rendered to HTML using `Sphinx`_ and the
`PyData Sphinx theme`_. Writing
docstrings is covered in the :ref:`numpy:howto-document`; this document
explains how to check that docstrings render properly.

*For a video walkthrough, please see* \ `Rendering SciPy Documentation
with Sphinx`_ \ *.*

.. _rendering-documentation-locally:

Rendering Documentation Locally
-------------------------------

To render the documentation on your own machine:

0. Ensure that you have a working SciPy :ref:`dev-env` active.
   You need to be able to ``import scipy`` regardless of Python's working
   directory; the ``python setup.py develop`` and ``conda develop`` commands
   from the :ref:`quickstart <dev-env>` guides make this possible.
#. Install `Sphinx`_, `PyData Sphinx theme`_ and `matplotlib`_. For
   example, if you're using the Anaconda distribution of Python, enter in a
   terminal window ``conda install sphinx pydata-sphinx-theme matplotlib --channel conda-forge``.
   The list of requirements is in ``scipy/doc_requirements.txt``.
#. In a terminal window, browse to the ``scipy/doc`` directory. Note the
   presence of the file ``Makefile``.
#. Execute ``git submodule update --init``.
   Some of the documentation theme files are not distributed
   with the main ``scipy`` repository; this keeps them up to date using
   `git submodules`_.
#. Enter ``make html-scipyorg``. If you have multiple version of Python on
   your path, you can choose which version to use by appending
   ``PYTHON=python3.7`` to this command, where ``python3.7`` is to be
   replaced with the name of the Python you use for SciPy development.
   This uses the `Make build automation tool`_
   to execute the documentation build instructions from the ``Makefile``.
   This can take a while the first time, but subsequent documentation builds
   are typically much faster.
#. View the documentation in ``scipy/doc/build/html-scipyorg``. You can start
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

.. _GitHub: https://github.com/
.. _CircleCI: https://circleci.com/vcs-authorize/
.. _Sphinx: https://www.sphinx-doc.org/en/master/
.. _PyData Sphinx theme: https://pydata-sphinx-theme.readthedocs.io/en/latest/
.. _matplotlib: https://www.matplotlib.org/
.. _Rendering SciPy Documentation with Sphinx: https://youtu.be/kGSYU39EhJQ
.. _git submodules: https://git-scm.com/book/en/v2/Git-Tools-Submodules
.. _Make build automation tool: https://en.wikipedia.org/wiki/Make_(software)
