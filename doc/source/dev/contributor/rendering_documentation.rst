:orphan:

.. _rendering-documentation:

===================================
Rendering Documentation with Sphinx
===================================

SciPy docstrings are rendered to html using `Sphinx`_. Writing
docstrings is covered in the :ref:`numpy:howto-document`; this document
explains how to check that docstrings render properly.

*For a video walkthrough, please see* \ `Rendering SciPy Documentation
with Sphinx`_ \ *.*

.. _rendering-documentation-locally:

Rendering Documentation Locally
-------------------------------

To render the documentation on your own machine (macOS or Linux):

#. Install `Sphinx`_ or ensure that your installation is up to date. For
   example, if you're using the Anaconda distribution of Python, enter in a
   terminal window ``conda install sphinx`` or ``conda update sphinx``.
#. In a terminal window, browse to the ``scipy/doc`` directory. Note the
   presence of the file ``Makefile``.
#. If this is your first time building the docs, execute ``git submodule
   init``. After you’ve initialized for the first time, enter ``git submodule
   update`` instead. Some of the documentation theme files are not distributed
   with the main ``scipy`` repository; this keeps them up to date with
   `git submodules`_.
#. Enter ``make html-scipyorg``. This uses the `Make build automation tool`_
   to execute the documentation build instructions from the ``Makefile``.
   This can take a while the first time, but subsequent documentation builds
   are typically much faster. *Note: If you use a virtual environment for
   development, activate it first.*
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
#. Next to “ci/circleci: build_docs”, select “Details”.
#. Select the “Artifacts” tab and find ``index.html``.

.. _GitHub: https://github.com/
.. _CircleCI: https://circleci.com/vcs-authorize/
.. _Sphinx: http://www.sphinx-doc.org/en/master/
.. _Rendering SciPy Documentation with Sphinx: https://youtu.be/kGSYU39EhJQ
.. _git submodules: https://git-scm.com/book/en/v2/Git-Tools-Submodules
.. _Make build automation tool: https://en.wikipedia.org/wiki/Make_(software)
