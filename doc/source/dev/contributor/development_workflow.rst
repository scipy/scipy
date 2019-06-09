.. _development-workflow:

====================
Development workflow
====================

You already have your own forked copy of the `SciPy repository`_, by
following :ref:`forking`, :ref:`set-up-fork`, you have configured git_
by following :ref:`configure-git`, and have linked the upstream
repository as explained in :ref:`linking-to-upstream`.

What is described below is a recommended workflow with Git.

Basic workflow
##############

In short:

1. Start a new *feature branch* for each set of edits that you do.
   See :ref:`below <making-a-new-feature-branch>`.

2. Hack away! See :ref:`below <editing-workflow>`

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

First, fetch new commits from the ``upstream`` repository:

::

   git fetch upstream

Then, create a new branch based on the master branch of the upstream
repository::

   git checkout -b my-new-feature upstream/master


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
   status`_).  You'll see a listing like this one::

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
   a default link to your github_ repo called ``origin``.  In git >= 1.7 you
   can ensure that the link to origin is permanently set by using the
   ``--set-upstream`` option::

      git push --set-upstream origin my-new-feature

   From now on git_ will know that ``my-new-feature`` is related to the
   ``my-new-feature`` branch in your own github_ repo. Subsequent push calls
   are then simplified to the following::

      git push

   You have to use ``--set-upstream`` for each new branch that you create.


It may be the case that while you were working on your edits, new commits have
been added to ``upstream`` that affect your work. In this case, follow the
:ref:`rebasing-on-master` instructions to apply those changes to your branch.

.. _writing-the-commit-message:

Writing the commit message
--------------------------

Commit messages should be clear and follow a few basic rules.  Example::

   ENH: add functionality X to SciPy.<submodule>.

   The first line of the commit message starts with a capitalized acronym
   (options listed below) indicating what type of commit this is.  Then a blank
   line, then more text if needed.  Lines shouldn't be longer than 72
   characters.  If the commit is related to a ticket, indicate that with
   "See #3456", "See ticket 3456", "Closes #3456" or similar.

Describing the motivation for a change, the nature of a bug for bug fixes or
some details on what an enhancement does are also good to include in a commit
message.  Messages should be understandable without looking at the code
changes.  A commit message like ``MAINT: fixed another one`` is an example of
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

-  Are there unit tests with good code coverage?
-  Do all public function have docstrings including examples?
-  Is the code style correct (PEP8, pyflakes)
-  Is the commit message `formatted correctly`_?
-  Is the new functionality tagged with ``.. versionadded:: X.Y.Z`` (with
   X.Y.Z the version number of the next release - can be found in setup.py)?
-  Is the new functionality mentioned in the release notes of the next
   release?
-  Is the new functionality added to the reference guide?
-  In case of larger additions, is there a tutorial or more extensive
   module-level description?
-  In case compiled code is added, is it integrated correctly via setup.py
-  If you are a first-time contributor, did you add yourself to THANKS.txt?
   Please note that this is perfectly normal and desirable - the aim is to
   give every single contributor credit, and if you don't add yourself it's
   simply extra work for the reviewer (or worse, the reviewer may forget).
-  Did you check that the code can be distributed under a BSD license?

.. include:: ../gitwash/git_links.inc

.. _formatted correctly: https://docs.scipy.org/doc/numpy/dev/gitwash/development_workflow.html#writing-the-commit-message
