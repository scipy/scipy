.. _contributor-toc:

=======================
SciPy contributor guide
=======================

This guide is designed to help you quickly find the information you need about
SciPy development after you've reviewed the introductory material in
:ref:`hacking` or :ref:`dev-quickstart`.

You can also watch `SciPy Development Workflow`_, a five-minute video example of
fixing a bug and submitting a pull request.

- :ref:`dev-env` - how to set up and maintain a development environment, including installing compilers and SciPy dependencies, creating a personal fork of the SciPy repository on GitHub, using git to manage a local repository with development branches, performing an in-place build of SciPy, and creating a virtual environment that adds this development version of SciPy to the Python path
- :ref:`editing-scipy` - how to edit SciPy Python code, with tips on finding which module contains SciPy functionality to be edited, adding new modules to SciPy, and complying with PEP8 style standards
- :ref:`unit-tests` - how to write and run unit tests for SciPy with the pytest framework
- :ref:`docs` - how to write reStructuredText documentation that complies with docstring standards, build documentation locally with Sphinx, and view documentation built during continuous integration checks
- :ref:`toc-benchmarking` - how to benchmark code with airspeed velocity
- :ref:`toc-cython` - how to add fast, compiled code to SciPy
- :ref:`continuous-integration` - how does our continuous integration system works and how to debug your PR

.. _dev-env:

Development environment
-----------------------

- :ref:`system-level` shows how to install system-level dependencies for Linux, Mac or Windows (needed if you're not using conda).
- :ref:`conda-guide` presents a step-by-step process for setting up a convenient SciPy development environment with conda *(recommended)*.
- :ref:`ubuntu-guide` presents a step-by-step process for setting up a convenient SciPy development environment in Ubuntu Linux.
- :ref:`quickstart-docker` presents a step-by-step process for building SciPy using Docker; if you have trouble with the instructions above, this may be your best option
- :ref:`quickstart-gitpod` presents a step-by-step process for using Gitpod for SciPy development; this process requires minimal setup and is newcomer friendly

.. _editing-scipy:

Editing SciPy
-------------
- :ref:`development-workflow` lays out what to do after your development environment is set up
- :ref:`building` has details on building from sources on Linux, Mac and Windows
- :ref:`meson` for how to use the Meson build system
- :ref:`pep8-scipy` gives some tips for ensuring that your code is PEP8 compliant
- :ref:`git-development` is a guide to using ``git``, the distributed version-control system used to manage the changes made to SciPy code from around the world
- :ref:`scipy-api` contains some important notes about how SciPy code is organized and documents the structure of the SciPy API; if you are going to import other SciPy code, read this first
- :ref:`reviewing-prs` explains how to review another author's SciPy code locally
- :ref:`adding-new` has information on how to add new methods, functions and classes
- :ref:`core-dev-guide` has background information including how decisions are made and how a release is prepared; it's geared toward :ref:`Core Developers <governance>`, but contains useful information for all contributors
- :ref:`missing-bits` - code and documentation style guide


.. _unit-tests:

Unit tests
----------
- :doc:`numpy:reference/testing` is the definitive guide to writing unit tests of NumPy or SciPy code (part of the NumPy documentation)
- :ref:`runtests` documents ``runtests.py``, a convenient script for building SciPy and running tests locally

.. _docs:

Documentation
-------------
- :ref:`numpy:howto-document` contains everything you need to know about writing docstrings, which are rendered to produce HTML documentation using `Sphinx`_ (part of the NumPy documentation)
- :ref:`rendering-documentation` it's important to check how changes to the documentation render before merging a PR; this document explains how you can do that

.. _toc-benchmarking:

Benchmarks
----------
- :ref:`benchmarking-with-asv` explains how to add benchmarks to SciPy using `airspeed velocity`_


.. _toc-cython:

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
