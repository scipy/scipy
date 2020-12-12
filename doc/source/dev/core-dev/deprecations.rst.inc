.. _deprecations:

Deprecations
============

There are various reasons for wanting to remove existing functionality: it's
buggy, the API isn't understandable, it's superseded by functionality with
better performance, it needs to be moved to another SciPy submodule, etc.

In general it's not a good idea to remove something without warning users about
that removal first.  Therefore this is what should be done before removing
something from the public API:

#. Propose to deprecate the functionality on the scipy-dev mailing list and get
   agreement that that's OK.
#. Add a ``DeprecationWarning`` for it, which states that the functionality was
   deprecated, and in which release. For Cython APIs, see
   :ref:`deprecating-public-cython-api` for the practical steps.
#. Mention the deprecation in the release notes for that release.
#. Wait till at least 6 months after the release date of the release that
   introduced the ``DeprecationWarning`` before removing the functionality.
#. Mention the removal of the functionality in the release notes.

The 6 months waiting period in practice usually means waiting two releases.
When introducing the warning, also ensure that those warnings are filtered out
when running the test suite so they don't pollute the output.

It's possible that there is reason to want to ignore this deprecation policy
for a particular deprecation; this can always be discussed on the scipy-dev
mailing list.

