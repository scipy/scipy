.. _continuous-integration:

======================
Continuous Integration
======================

Continuous integration is part of our development process and ensure that
every piece of code or documentation which is contributed to SciPy is working
and does not have unforeseen effects.

.. note:: Before submitting or updating your PR, please ensure that you tested
          locally your changes. See :ref:`pr-checklist`.

Workflows
=========

We run more than 20 different workflows with different version of the
dependencies, different architectures, etc. A PR must pass all these checks
before it can be merged as to ensure a sustainable state of the project.

Apart from the unit tests, the documentation and examples in the docstrings are
also checked. These are common failing workflows as Sphinx and doctests have
very strict rules. These aspects are very important as documentation and
examples are user facing elements. Ensure that these elements are properly
rendered.

The logs can be long, but you will always find out why your did not pass a
check. Simply click on ``Details`` to access the logs.

Following is a list of all the different workflows in use. They are gouped
by CI resources providers.

GitHub Actions
--------------
* ``Linux Tests``: Python 3.7 and nightly Cpython 3.9
* ``macOS Tests``: test matrix from Python 3.7 to 3.9

Azure
-----
* ``Lint``: PEP8 and code style
* ``Windows Python``: test matrix from Python 3.7 to 3.9
* ``Linux_Python_37_32bit_full``
* ``wheel_optimized_gcc48``
* ``source_distribution``
* ``refguide_asv_check``: doctests from examples

CircleCI
--------
* ``build_docs``: build the documentation
* ``build_docs artifact``: live preview of the documentation
* ``build_scipy``
* ``run_benchmarks``: verify how the changes impact performance

Codecov
-------
* ``patch``: the impact on code coverage due to your changes
* ``project``: the coverage of the whole project

Skipping
========

Being an open-source project, we have access to a quota of CI ressources.
Ultimately, resources are limited and we should use them with care. This is
why we ask you to verify your changes locally before pushing them.

Depending on the proposed change, you might want to skip part of the checks.
It will be to the discretion of a maintainer to re-run some tests before
integration.

Skipping CI can be achieve by adding a special text in the commit message:

* ``[skip azp]``: will skip Azure
* ``[skip actions]``: will skip GitHub Actions
* ``[skip ci]``: will skip CicleCI

Of course, you can combine these to skip multiple workflows.

This skip information should be placed on a new line. In this example, we
just updated a ``.rst`` file in the documentation and ask to skip Azure and
GitHub Actions' workflows::

    DOC: improve QMCEngine examples
    [skip azp] [skip actions]
