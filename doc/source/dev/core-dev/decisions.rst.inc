Decision making process
=======================

SciPy has a formal governance model, documented in :ref:`governance`.  The
section below documents in an informal way what happens in practice for
decision making about code and commit rights. The formal governance model is
leading, the below is only provided for context.

Code
----

Any significant decisions on adding (or not adding) new features, breaking
backwards compatibility or making other significant changes to the codebase
should be made on the scipy-dev forum after a discussion (preferably
with full consensus).

Any non-trivial change (where trivial means a typo, or a one-liner maintenance
commit) has to go in through a pull request (PR).  It has to be reviewed by
another developer.  In case review doesn't happen quickly enough and it is
important that the PR is merged quickly, the submitter of the PR should send a
message to the forum saying they intend to merge that PR without review
at time X for reason Y unless someone reviews it before then.

Changes and new additions should be tested. Untested code is broken code.

Commit rights
-------------

Who gets commit rights is decided by the SciPy Steering Council; changes in
commit rights will then be announced on the scipy-dev forum.
