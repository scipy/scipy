Development on GitHub 
=====================
SciPy development largely takes place on GitHub; this section describes the
expected way of working for issues, pull requests and managing the main
``scipy`` repository.

.. _labels-and-milestones:

Labels and Milestones
---------------------
Each issue and pull request normally gets at least two labels: one for the
topic or component (``scipy.stats``, ``Documentation``, etc.), and one for the
nature of the issue or pull request (``enhancement``, ``maintenance``,
``defect``, etc.).  Other labels that may be added depending on the situation:

- ``good-first-issue``: for issues suitable to be tackled by new contributors.
- ``needs-work``: for pull requests that have review comments that haven't been
  addressed.
- ``needs-decision``: for issues or pull requests that need a decision.
- ``needs-champion``: for pull requests that were not finished by the original
  author, but are worth resurrecting.
- ``backport-candidate``: bugfixes that should be considered for backporting by
  the release manager.

A milestone is created for each version number for which a release is planned.
Issues that need to be addressed and pull requests that need to be merged for a
particular release should be set to the corresponding milestone.  After a pull
request is merged, its milestone (and that of the issue it closes) should be
set to the next upcoming release - this makes it easy to get an overview of
changes and to add a complete list of those to the release notes.


Pull request review workflow
----------------------------

When reviewing pull requests, please make use of pull request workflow
features, see :ref:`pull-request-workflow-features`.


Dealing with pull requests
--------------------------

- When merging contributions, a committer is responsible for ensuring that
  those meet the requirements outlined in :ref:`Contributing to SciPy <hacking>`.
  Also check that new features and backwards compatibility breaks were discussed
  on the scipy-dev mailing list.
- New code goes in via a pull request (PR).
- Merge new code with the green button.  In case of merge conflicts, ask the PR
  submitter to rebase (this may require providing some ``git`` instructions).
- Backports and trivial additions to finish a PR (really trivial, like a typo
  or PEP8 fix) can be pushed directly.
- For PRs that add new features or are in some way complex, wait at least a day
  or two before merging it.  That way, others get a chance to comment before
  the code goes in.
- Squashing commits or cleaning up commit messages of a PR that you consider
  too messy is OK. Make sure though to retain the original author name when
  doing this. Squashing is highly recommended whenever commit messages do not
  (roughly) follow the guidelines in :ref:`writing-the-commit-message`.
- Make sure that the labels and milestone on a merged PR are set correctly.
- When you want to reject a PR: if it's very obvious, you can just close it and
  explain why. If it's not obvious, then it's a good idea to first explain why you
  think the PR is not suitable for inclusion in SciPy, then let a second
  committer comment or close.


.. _backporting:

Backporting
-----------
All pull requests (whether they contain enhancements, bug fixes or something else),
should be made against main.  Only bug fixes are candidates for backporting
to a maintenance branch.  The backport strategy for SciPy is to (a) only backport
fixes that are important, and (b) to only backport when it's reasonably sure
that a new bugfix release on the relevant maintenance branch will be made.
Typically, the developer who merges an important bugfix adds the
``backport-candidate`` label and pings the release manager, who decides on
whether and when the backport is done.  After the backport is completed, the
``backport-candidate`` label has to be removed again.

A good strategy for a backport pull request is to combine several main
branch pull requests, to reduce the burden on continuous integration tests
and to reduce the merge commit cluttering of maintenance branch history. It
is generally best to have a single commit for each of the main branch pull
requests represented in the backport pull request. This way, history is clear
and can be reverted in a straightforward manner if needed.

Release notes
-------------
When a PR gets merged, consider if the changes need to be mentioned in the
release notes.  What needs mentioning: new features, backwards incompatible
changes, deprecations, and "other changes" (anything else noteworthy enough,
see older release notes for the kinds of things worth mentioning).

Release note entries are maintained on
`the wiki <https://github.com/scipy/scipy/wiki>`_.
The release manager will gather content from there and integrate it into the html
docs.  We use this mechanism to avoid merge conflicts that would happen if
every PR touched the same file under ``doc/release/`` directly.

Changes can be monitored (`Atom feed <https://github.com/scipy/scipy/wiki.atom>`_)
and pulled (the wiki is a git repo: ``https://github.com/scipy/scipy.wiki.git``).


Other
-----

*Cross-referencing:* Cross-referencing issues and pull requests on GitHub is
often useful.  GitHub allows doing that by using ``gh-xxxx`` or ``#xxxx`` with
``xxxx`` the issue/PR number.  The ``gh-xxxx`` format is strongly preferred,
because it's clear that that is a GitHub link.  Older issues contain ``#xxxx``
which is about Trac (what we used pre-GitHub) tickets.

*PR naming convention:* Pull requests, issues and commit messages usually start
with a three-letter abbreviation like ``ENH:`` or ``BUG:``.  This is useful to
quickly see what the nature of the commit/PR/issue is.  For the full list of
abbreviations, see :ref:`writing the commit message
<numpy:writing-the-commit-message>`.
