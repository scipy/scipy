:orphan:

.. highlight:: console

.. _recommended-development-setup:

=============================
Recommended development setup
=============================

*This document does not include detailed explanations. For more step-by-step
procedures, see* :ref:`quickstart-mac` *or* :ref:`building`.

Since SciPy contains parts written in C, C++, and Fortran that need to be
compiled before use, make sure you have the necessary compilers and Python
development headers installed.  Having compiled code also means that importing
SciPy from the development sources needs some additional steps, which are
explained below.

First fork a copy of the main SciPy repository in Github onto your own
account and then create your local repository via::

    $ git clone git@github.com:YOURUSERNAME/scipy.git scipy
    $ cd scipy
    $ git submodule update --init
    $ git remote add upstream git://github.com/scipy/scipy.git

Second to code review pull requests it is helpful to have a local copy of the
code changes in the pull request. The preferred method to bring a PR from the
github repository to your local repo in a new branch::

    $ git fetch upstream pull/PULL_REQUEST_ID/head:NEW_BRANCH_NAME

The value of ``PULL_REQUEST_ID`` will be the PR number and the
``NEW_BRANCH_NAME`` will be the name of the branch in your local repository
where the diffs will reside.

Now you have a branch in your local development area to code review in Python.

To build the development version of SciPy and run tests, spawn
interactive shells with the Python import paths properly set up etc.,
do one of::

    $ python runtests.py -v
    $ python runtests.py -v -s optimize
    $ python runtests.py -v -t scipy.special.tests.test_basic::test_xlogy
    $ python runtests.py --ipython
    $ python runtests.py --python somescript.py
    $ python runtests.py --bench

This builds SciPy first, so the first time it may take some time.  If
you specify ``-n``, the tests are run against the version of SciPy (if
any) found on current PYTHONPATH.  *Note: if you run into a build issue,
more detailed build documentation can be found in* :ref:`building`.

Some of the tests in SciPy are very slow and need to be separately
enabled. See :ref:`runtests` for details.