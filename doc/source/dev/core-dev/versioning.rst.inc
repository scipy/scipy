.. _version-numbering:

Version numbering
=================
SciPy version numbering complies to `PEP 440`_.  Released final versions, which
are the only versions appearing on `PyPI`_, are numbered ``MAJOR.MINOR.MICRO``
where:

- ``MAJOR`` is an integer indicating the major version.  It changes very
  rarely; a change in ``MAJOR`` indicates large (possibly backwards-incompatible)
  changes.
- ``MINOR`` is an integer indicating the minor version.  Minor versions are
  typically released twice a year and can contain new features, deprecations and
  bug-fixes.
- ``MICRO`` is an integer indicating a bug-fix version.  Bug-fix versions are
  released when needed, typically one or two per minor version.  They cannot
  contain new features or deprecations.

Released alpha, beta and rc (release candidate) versions are numbered
like final versions but with postfixes ``a#``, ``b#`` and ``rc#`` respectively,
with ``#`` an integer.  Development versions are postfixed with ``.dev0+<git-commit-hash>``.

Examples of valid SciPy version strings are::

    0.16.0
    0.15.1
    0.14.0a1
    0.14.0b2
    0.14.0rc1
    0.17.0.dev0+ac53f09

An installed SciPy version contains these version identifiers::

    scipy.__version__            # complete version string, including git commit hash for dev versions
    scipy.version.short_version  # string, only major.minor.micro
    scipy.version.version        # string, same as scipy.__version__
    scipy.version.full_version   # string, same as scipy.__version__
    scipy.version.release        # bool, development or (alpha/beta/rc/final) released version
    scipy.version.git_revision   # string, git commit hash from which scipy was built


.. _PEP 440: https://www.python.org/dev/peps/pep-0440
.. _PyPI: https://pypi.python.org/
