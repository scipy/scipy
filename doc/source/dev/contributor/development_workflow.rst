:orphan:

.. _development-workflow:

====================
Development workflow
====================

*Note: consider watching* `SciPy Development Workflow`_ *before or after
reading to see an example of fixing a bug and submitting a pull request.*

This guide assumes that you have created your own fork (copy) of the SciPy
repository, cloned the repository on your own machine, and built SciPy from this
source code. If you haven't, check the :ref:`dev-env` pages appropriate to your
system. Before getting started here, there are two other things you need to do
just once before you start modifying SciPy.

#. In a terminal, introduce yourself to Git::

      git config --global user.email you@yourdomain.com
      git config --global user.name "Your Name"

   This information credits you for your work, but note that it will become
   publicly available if you "push" your work to GitHub. See
   `Setting your commit email address in Git`_ for more information.

#. Navigate to the root directory of your local SciPy repository and enter::

      git remote add upstream https://github.com/scipy/scipy.git

   This associates the name ``upstream`` with the official SciPy repository
   located at `https://github.com/scipy/scipy.git <https://github.com/scipy/scipy.git>`_.
   Note that when you cloned your fork of the SciPy repository, Git already
   associated the name ``origin`` with your fork. The reason you need both of
   these `"remotes"`_ is that you will typically start with the latest version
   of SciPy from the official repository ``upstream``, make changes, "push"
   your changes to your fork of the repository ``origin``, and then submit
   a "pull request" asking SciPy to "pull" your changes from your fork into
   the official repository.

#. Initialize git submodules::

      git submodule update --init

   This fetches and updates any submodules that SciPy needs (such as `Boost`).

Basic workflow
##############

In short:

1. Start a new *feature branch* for each set of edits that you do.
   See :ref:`below <making-a-new-feature-branch>`.

2. Hack away! See :ref:`below <editing-workflow>`.

3. When finished:

   - *Contributors*: push your feature branch to your own Github repo, and
     :ref:`create a pull request <asking-for-merging>`.

   - *Core developers* If you want to push changes without
     further review, see the notes :ref:`below <pushing-to-main>`.

This way of working helps to keep work well organized and the history
as clear as possible.

.. seealso::

   There are many online tutorials to help you `learn git`_. For discussions
   of specific git workflows, see these discussions on `linux git workflow`_,
   and `ipython git workflow`_.

.. _making-a-new-feature-branch:

Making a new feature branch
===========================

First, navigate to the SciPy root directory in your terminal and fetch new
commits from the ``upstream`` repository::

   git fetch upstream

Then, create a new branch based on the main branch of the upstream
repository::

   git checkout -b my-new-feature upstream/main

Equivalently, you might want to keep the main branch of your own repository
up to date and create a new branch based on that::

   git checkout main
   git rebase upstream/main
   git checkout -b my-new-feature

In order, these commands

#. ensure that the ``main`` branch of your local repository is checked out,

#. apply all the latest changes from the ``upstream/main`` (main SciPy
   repository main branch) to your local ``main`` branch, and

#. create and check out a new branch (``-b``) based on your local ``main``
   branch.

In any case, it's important that your feature branch include the latest
changes from the upstream main to help avoid
`merge conflicts <https://help.github.com/en/articles/resolving-a-merge-conflict-using-the-command-line>`_
when it's time to submit a pull request.

It's also a good idea to build this branch and run tests before continuing.
Assuming you've followed one of the :ref:`dev-env` pages to set up your
development environment, you'll need to activate your development virtual
environment, perform an in-place build, and run tests::

   conda activate name-of-your-virtual-environment
   python setup.py build_ext --inplace
   python runtests.py -v

Otherwise, see :ref:`building`, :ref:`runtests` for more information.

.. _editing-workflow:

The editing workflow
====================

Overview
--------

::

   # hack hack
   git status # Optional
   git diff # Optional
   git add modified_file
   git commit
   # push the branch to your own Github repo
   git push origin my-new-feature

In more detail
--------------

#. Make some changes. When you feel that you've made a complete, working set
   of related changes, move on to the next steps.

#. Optional: Check which files have changed with ``git status`` (see `git
   status`_). You'll see a listing like this one::

     # On branch my-new-feature
     # Changed but not updated:
     #   (use "git add <file>..." to update what will be committed)
     #   (use "git checkout -- <file>..." to discard changes in working directory)
     #
     #	modified:   README
     #
     # Untracked files:
     #   (use "git add <file>..." to include in what will be committed)
     #
     #	INSTALL
     no changes added to commit (use "git add" and/or "git commit -a")

#. Optional: Compare the changes with the previous version using with ``git
   diff`` (`git diff`_). This brings up a simple text browser interface that
   highlights the difference between your files and the previous version.

#. Add any relevant modified or new files using  ``git add modified_file``
   (see `git add`_). This puts the files into a staging area, which is a queue
   of files that will be added to your next commit. Only add files that have
   related, complete changes. Leave files with unfinished changes for later
   commits.

#. To commit the staged files into the local copy of your repo, do ``git
   commit``. At this point, a text editor will open up to allow you to write a
   commit message. Read the :ref:`commit message
   section<writing-the-commit-message>` to be sure that you are writing a
   properly formatted and sufficiently detailed commit message. After saving
   your message and closing the editor, your commit will be saved. For trivial
   commits, a short commit message can be passed in through the command line
   using the ``-m`` flag. For example, ``git commit -am "ENH: Some message"``.

   In some cases, you will see this form of the commit command: ``git commit
   -a``. The extra ``-a`` flag automatically commits all modified files and
   removes all deleted files. This can save you some typing of numerous ``git
   add`` commands; however, it can add unwanted changes to a commit if you're
   not careful. For more information, see `why the -a flag?`_ - and the
   helpful use-case description in the `tangled working copy problem`_.

#. Push the changes to your forked repo on github_::

      git push origin my-new-feature

   For more information, see `git push`_.

.. note::

   Assuming you have followed the instructions in these pages, git will create
   a default link to your github_ repo called ``origin``. In git >= 1.7, you
   can ensure that the link to origin is permanently set by using the
   ``--set-upstream`` option::

      git push --set-upstream origin my-new-feature

   From now on, git_ will know that ``my-new-feature`` is related to the
   ``my-new-feature`` branch in your own github_ repo. Subsequent push calls
   are then simplified to the following::

      git push

   You have to use ``--set-upstream`` for each new branch that you create.


It may be the case that while you were working on your edits, new commits have
been added to ``upstream`` that affect your work. In this case, follow the
:ref:`rebasing-on-main` instructions to apply those changes to your branch.

.. _writing-the-commit-message:

Writing the commit message
--------------------------

Commit messages should be clear and follow a few basic rules.  Example::

   ENH: add functionality X to SciPy.<submodule>.

   The first line of the commit message starts with a capitalized acronym
   (options listed below) indicating what type of commit this is. Then a blank
   line, then more text if needed.  Lines shouldn't be longer than 72
   characters.  If the commit is related to a ticket, indicate that with
   "See #3456", "See ticket 3456", "Closes #3456", or similar.

Describing the motivation for a change, the nature of a bug for bug fixes or
some details on what an enhancement does are also good to include in a commit
message. Messages should be understandable without looking at the code
changes. A commit message like ``MAINT: fixed another one`` is an example of
what not to do; the reader has to go look for context elsewhere.

Standard acronyms to start the commit message with are::

   API: an (incompatible) API change
   BENCH: changes to the benchmark suite
   BLD: change related to building SciPy
   BUG: bug fix
   DEP: deprecate something, or remove a deprecated object
   DEV: development tool or utility
   DOC: documentation
   ENH: enhancement
   MAINT: maintenance commit (refactoring, typos, etc.)
   REV: revert an earlier commit
   STY: style fix (whitespace, PEP8)
   TST: addition or modification of tests
   REL: related to releasing SciPy

.. note:: You can add some markers to skip part of the continuous integration.
          See :ref:`continuous-integration`.

.. _asking-for-merging:

Asking for your changes to be merged with the main repo
-------------------------------------------------------

When you feel your work is finished, you can create a pull request (PR). Github
has a nice help page that outlines the process for `filing pull requests`_.

If your changes involve modifications to the API or addition/modification of a
function, you should initiate a code review. This involves sending an email to
the `SciPy mailing list`_ with a link to your PR along with a description of
and a motivation for your changes.

.. _pr-checklist:

Checklist before submitting a PR
--------------------------------

-  Did you check that the code can be distributed under a BSD license? See
   :ref:`license-considerations`.
-  Are there unit tests with good code coverage? See
   `NumPy/SciPy Testing Guidelines`_.
-  Do all unit tests pass locally? See :ref:`runtests`.
-  Do all public function have docstrings including examples? See the
   `numpydoc docstring guide`_.
-  Does the documentation render correctly? See :ref:`rendering-documentation`.
-  Is the code style correct? See :ref:`pep8-scipy`.
-  Are there benchmarks? See :ref:`benchmarking-with-asv`.
-  Is the commit message :ref:`formatted correctly <numpy:writing-the-commit-message>`?
-  Is the docstring of the new functionality tagged with
   ``.. versionadded:: X.Y.Z`` (where ``X.Y.Z`` is the version number of the
   next release? See the ``updating``, ``workers``, and ``constraints``
   documentation of |differential_evolution|_, for example. You can get the
   next version number from the most recent release notes on `the wiki`_ or
   from the ``MAJOR`` and ``MINOR`` version number variables in |setup.py|_.
-  In case of larger additions, is there a tutorial or more extensive
   module-level description? Tutorial files are in ``doc/source/tutorial``.
-  If compiled code is added, is it integrated correctly via ``setup.py``?
   See :ref:`compiled-code` for more information.

.. include:: ../gitwash/git_links.inc

.. _Scipy Development Workflow: https://youtu.be/HgU01gJbzMY

.. _Setting your commit email address in Git: https://help.github.com/en/articles/setting-your-commit-email-address-in-git

.. _"remotes": https://help.github.com/en/categories/managing-remotes

.. _NumPy/SciPy Testing Guidelines: https://docs.scipy.org/doc/numpy/reference/testing.html

.. _numpydoc docstring guide: https://numpydoc.readthedocs.io/en/latest/format.html

.. _the wiki: https://github.com/scipy/scipy/wiki

.. |differential_evolution| replace:: ``differential_evolution``
.. _differential_evolution: https://github.com/scipy/scipy/blob/main/scipy/optimize/_differentialevolution.py

.. |setup.py| replace:: ``setup.py``
.. _setup.py: https://github.com/scipy/scipy/blob/main/setup.py
