Decision making process
=======================
This section documents the way in which decisions about various aspects of the
SciPy project are made.  Note that the below is only documenting the current
way of working; a more formal governance model is expected to be adopted in the
near future.

Code
----
Any significant decisions on adding (or not adding) new features, breaking
backwards compatibility or making other significant changes to the codebase
should be made on the scipy-dev mailing list after a discussion (preferably
with full consensus).

Any non-trivial change (where trivial means a typo, or a one-liner maintenance
commit) has to go in through a pull request (PR).  It has to be reviewed by
another developer.  In case review doesn't happen quickly enough and it is
important that the PR is merged quickly, the submitter of the PR should send a
message to mailing list saying he/she intends to merge that PR without review
at time X for reason Y unless someone reviews it before then.

Changes and new additions should be tested. Untested code is broken code.

Commit rights
-------------
Who gets commit rights is decided by the core development team; changes in
commit rights will then be announced on the scipy-dev mailing list.

Who the core development team is comprised of is a little fuzzy - there are
quite a few people who do have commit rights and would like to keep them but
are no longer active (so they're not in the core team).
To get an idea, look at the output of::

    $ git shortlog --grep="Merge pull request" -a -c -s <current_release_minus_2>..upstream/master|sort -n

and apply some common sense to it (and don't forget people who are still active
but never merge PRs).

Other project aspects
---------------------
All decisions are taken by the core development team on the scipy-dev mailing
list.
