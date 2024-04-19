.. _deprecations:

Deprecations
============

There are various reasons for wanting to remove existing functionality: it's
buggy, the API isn't understandable, it's superseded by functionality with
better performance, it needs to be moved to another SciPy submodule, etc.

In general, it's not a good idea to remove something or make other backwards
compatibility breaking changes without warning users about
the change first. Therefore this is what should be done before deprecating
something from the public API:

#. Propose to deprecate the functionality on the scipy-dev mailing list and get
   agreement that that's OK.
#. When changing a function signature (e.g. deprecating a keyword argument), this
   opportunity should be used to convert that function to keyword-only usage,
   see `GH-14714<https://github.com/scipy/scipy/issues/14714>`_.
#. Add a ``DeprecationWarning`` that is triggered by using the deprecated API.
   The ``stacklevel`` must be set correctly so that the warning refers
   to the function’s caller, rather than the ``warnings.warn()`` call itself.
   For Cython APIs, see
   :ref:`deprecating-public-cython-api` for the practical steps.
#. For behavior changes (rather than outright removals), it must be possible for
   the user to silence the warning without doing manual filtering or warnings.
   In other words, a replacement API needs to be ready, and the warning should
   contain clear instructions what needs to be done to avoid the warning.
#. Use the ``.. deprecated::`` directive to document the deprecation. The rest
   of the documentation should be left until the function is removed and the
   function should be left in the toctree. Deprecation warnings causing the doc
   build to fail can be resolved by adding an entry to the
   `conf.py <https://github.com/scipy/scipy/blob/main/doc/source/conf.py>`_
   file.
#. Ensure that both these warning messages contain the following information:

   #. What has been deprecated?
   #. What version the deprecation will be executed?
   #. What the user should use instead?

#. Add a test using ``pytest.deprecated_call()`` to ensure the deprecation warning
   is triggered as expected. Where applicable, also test that the new and old
   API give the same result.
#. Mention the deprecation in the release notes for that release and add it to the
   `deprecation tracker <https://github.com/scipy/scipy/issues/15765>`_.
#. Wait at least two releases from the release that
   introduced the ``DeprecationWarning`` before removing the functionality.
#. Mention the removal of the functionality in the release notes.

There may be reason to want to ignore this deprecation policy
for a particular deprecation; this can always be discussed on the scipy-dev
mailing list.

