:orphan:

.. _reviewing-prs:

=======================
Reviewing Pull Requests
=======================

.. _pull-request-workflow-features:

Using workflow features
-----------------------

When reviewing pull requests, please use workflow tracking features on
Github as appropriate:

1. After you have finished reviewing, and want to ask for the submitter
   to make the changes:

   - Change your review status to "Changes requested".

     This can be done on Github, PR page, ``Files changed`` tab,
     ``Review changes`` (button on top right).

   - Alternatively: add the ``needs-work`` label.

     This can be done on the PR page, ``Labels`` menu on the right.

2. When you re-review the same pull request again, and want to request
   more changes:

   - Do the "Changes requested" thing again, even if the previous status
     was also 'Changes requested'.

   - Alternatively:
     Remove the existing ``needs-work`` label, and then re-add the label
     back again. (Github will add a notice on the page that you did so.)

3. If you're happy about the current status:

   - Mark the pull request as Approved (same way as Changes requested).

   - Alternatively: remove the ``needs-work`` label.

   - Alternatively (for core developers): merge the pull request, if
     you think it is ready to be merged.

This allows automatically tracking which PRs are in need of attention.

The review status is listed at: https://pav.iki.fi/scipy-needs-work/
The page can also be generated using https://github.com/pv/github-needs-work

Some of the information is also visible on Github directly, although
(as of Aug 2019) Github does not show which pull requests have been
updated since the last review.


Code from pull request
----------------------

When you review a pull request created by someone else, it's helpful to have a
copy of their code on your own machine so that you can play with it locally.

One way to do this is to navigate to the SciPy root directory in the terminal
and enter::

   git fetch upstream pull/PULL_REQUEST_ID/head:NEW_BRANCH_NAME

where ``PULL_REQUEST_ID`` is the five digit number corresponding with the
pull request (e.g. ``10286`` for `PR #10286`_) and ``NEW_BRANCH_NAME`` is
whatever name you'd like to use to refer to the author's code (e.g.
``review_10286``).

Now you can check out the branch::

   git checkout NEW_BRANCH_NAME

which converts the code in your local repository to match the author's modified
version of SciPy.

Assuming you set up your development environment according to
:ref:`quickstart-mac` or :ref:`quickstart-ubuntu`, you you can now activate your development environment::

   conda activate scipydev

build the code and test it::

   python setup.py build_ext --inplace
   python runtests.py -v

and if you ``import`` SciPy from Python, you'll be importing the
author's modified version of SciPy.

If you want to collaborate with the author on their PR, you might instead
want to set up a new remote to the author's fork of SciPy::

   git remote add REMOTE_NAME https://github.com/AUTHOR/scipy.git

where ``AUTHOR`` is the author's GitHub user name and ``REMOTE_NAME`` is
whatever name you want to use to refer to this author's repository.

From there, you can view the author's branches::

   git remote show REMOTE_NAME

and create your own branch based on one of them::

   git checkout --track REMOTE_NAME/BRANCH_NAME

where ``BRANCH_NAME`` is the name of the branch you want to start from. This
creates a copy of this branch (with the same name) in your local repository.
If make changes to this branch and push to your GitHub repository
(``origin``), you can then create a pull request to merge your changes with the
author's repository.

.. _PR #10286: https://github.com/scipy/scipy/pull/10286
