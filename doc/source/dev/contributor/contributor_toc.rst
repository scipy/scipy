.. _contributor-toc:

=======================
SciPy contributor guide
=======================

This guide is designed to help you quickly find the information you need about
SciPy development after you've reviewed the introductory material in
:ref:`hacking` or :ref:`dev-quickstart`.

You can also watch `SciPy Development Workflow`_, a five-minute video example of
fixing a bug and submitting a pull request (*note: this video is from 2018, so
the build steps are different by now - the overall workflow is still the same
though*).

- :ref:`building-from-source` - how to set up a development environment,
  including installing compilers and SciPy dependencies, cloning the SciPy
  repository on GitHub and updating git submodules, and using the ``dev.py``
  interface for building and running tests.
- :ref:`editing-scipy` - how to edit SciPy Python code, with tips on finding
  which module contains SciPy functionality to be edited, adding new modules to
  SciPy, and complying with PEP8 style standards
- :ref:`unit-tests` - how to write and run unit tests for SciPy with the pytest
  framework
- :ref:`docs` - how to write reStructuredText documentation that complies with
  docstring standards, build documentation locally with Sphinx, and view
  documentation built during continuous integration checks
- :ref:`toc-benchmarking` - how to benchmark code with airspeed velocity
- :ref:`compiled-code` - how to add fast, compiled code to SciPy
- :ref:`continuous-integration` - how does our continuous integration system
  works and how to debug your PR

.. _editing-scipy:

Editing SciPy
-------------
- :ref:`development-workflow` lays out what to do after your development environment is set up
- :ref:`pep8-scipy` gives some tips for ensuring that your code is PEP8 compliant
- :ref:`git-development` is a guide to using ``git``, the distributed version-control system used to manage the changes made to SciPy code from around the world
- :ref:`scipy-api` contains some important notes about how SciPy code is organized and documents the structure of the SciPy API; if you are going to import other SciPy code, read this first
- :ref:`reviewing-prs` explains how to review another author's SciPy code locally
- :ref:`triaging` explains how to curate issues and PRs, as well as how GitHub team permissions work for SciPy
- :ref:`adding-new` has information on how to add new methods, functions and classes
- :ref:`core-dev-guide` has background information including how decisions are made and how a release is prepared; it's geared toward :ref:`Core Developers <governance>`, but contains useful information for all contributors
- :ref:`missing-bits` - code and documentation style guide


.. _unit-tests:

Testing
-------
- :doc:`numpy:reference/testing` is the definitive guide to writing unit tests
  of NumPy or SciPy code (part of the NumPy documentation)
- :ref:`devpy-test` documents ``dev.py test``, the command to build SciPy and
  run tests locally
- :ref:`debugging-linalg-issues` 

.. _docs:

Documentation
-------------
- :ref:`numpy:howto-document` contains everything you need to know about writing docstrings, which are rendered to produce HTML documentation using `Sphinx`_ (part of the NumPy documentation)
- :ref:`contributing-docs` contains information on how to contribute to the SciPy documentation
- :ref:`rendering-documentation` it's important to check how changes to the documentation render before merging a PR; this document explains how you can do that

.. _toc-benchmarking:

Benchmarks
----------
- :ref:`benchmarking-with-asv` explains how to add benchmarks to SciPy using `airspeed velocity`_


.. _compiled-code:

Compiled code
-------------
- :ref:`adding-cython` extending and compiling Python code with `Cython`_ can significantly improve its performance; this document helps you get started
- :ref:`other-languages` discusses the use of C, C++, and Fortran code in SciPy
- :ref:`public-cython-api` on guidelines on exposing public Cython APIs

.. _Scipy Development Workflow: https://youtu.be/HgU01gJbzMY

.. _Sphinx: http://www.sphinx-doc.org/en/master/

.. _Airspeed Velocity: https://asv.readthedocs.io/en/stable/

.. _Cython: https://cython.org/

.. |*| replace:: \ :sup:`*` \

.. toctree::
    :hidden:

    development_workflow
    pep8
    ../gitwash/gitwash
    ../../reference/index
    reviewing_prs
    ../triage
    adding_new
    ../core-dev/index
    ../missing-bits
    NumPy testing guidelines <https://numpy.org/devdocs/reference/testing.html>
    devpy_test
    debugging_linalg_issues
    How to contribute documentation <https://numpy.org/devdocs/dev/howto-docs.html>
    rendering_documentation
    benchmarking
    cython
    compiled_code
    public_cython_api
