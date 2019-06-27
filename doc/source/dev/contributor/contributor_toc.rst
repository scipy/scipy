.. _contributor-toc:

=======================
SciPy Contributor Guide
=======================

This guide is designed to help you quickly find the information you need about SciPy development after you've reviewed the introductory material in :ref:`hacking`. If you're new to this and want to start coding ASAP, you've found the right place.

- :ref:`dev-env` - how to set up and maintain a development environment, including installing compilers and SciPy dependencies, creating a personal fork of the SciPy repository on GitHub, using git to manage a local repository with development branches, performing an in-place build of SciPy, and creating a virtual environment that adds this development version of SciPy to the Python path
- :ref:`editing-scipy` - how to edit SciPy Python code, with tips on finding which module contains SciPy functionality to be edited, adding new modules to SciPy, and complying with PEP8 style standards
- :ref:`unit-tests` - how to write and run unit tests for SciPy with the pytest framework
- :ref:`docs` - how to write reStructuredText documentation that complies with docstring standards, build documentation locally with Sphinx, and view documentation built during continuous integration checks
- :ref:`toc-benchmarking` - how to benchmark code with AirSpeed Velocity
- :ref:`toc-cython` - how to add fast, compiled code to SciPy

.. _dev-env:

Development Environment
-----------------------
- :ref:`quickstart-mac` presents a step-by-step process for setting up a convenient SciPy development environment in macOS.
- :ref:`quickstart-ubuntu` presents a step-by-step process for setting up a convenient SciPy development environment in Ubuntu.
- :ref:`building` - If you don't have macOS or Ubuntu, try these instructions to help you build SciPy on your operating system.
- :ref:`recommended-development-setup` includes additional notes about the development setup. All of this information is contained elsewhere, but it is retained as a legacy document.

.. _editing-scipy:

Editing SciPy
-------------
- :ref:`development-workflow` lays out what to do after your development environment is set up.
- `SciPy Development Workflow`_ is a five-minute video example of fixing a bug and submitting a pull request.
- :ref:`pep8-scipy` gives some tips for ensuring that your code is PEP8 compliant.
- :ref:`git-development` is a guide to using ``git``, the distributed version-control system used to manage the changes made to SciPy code from around the world.
- :ref:`scipy-api` contains some important notes about how SciPy code is organized and documents the structure of the SciPy API. If you are going to import other SciPy code, read this first.
- :ref:`reviewing-prs` explains how to review another author's SciPy code locally.
- :doc:`numpy:reference/distutils_guide` - Check this out before adding any new files to SciPy.
- :ref:`core-dev-guide` has background information including how decisions are made and how a release is prepared. It's geared toward :ref:`Core Developers<governance>`, but contains useful information for all contributors.

.. _unit-tests:

Unit Tests
----------
- :doc:`numpy:reference/testing` is the definitive guide to writing unit tests of SciPy code.
- :ref:`runtests` documents ``runtests.py``, a convenient script for building SciPy and running tests locally.

.. _docs:

Documentation
-------------
- :ref:`numpy:howto-document` contains everything you need to know about writing docstrings, which are rendered to produce HTML documentation using `Sphinx`_.
- :ref:`rendering-documentation` - It's important to check how changes to the documentation render before merging a PR; this document explains how you can do that.

.. _toc-benchmarking:

Benchmarks
----------
- :ref:`benchmarking-with-asv` explains how to add benchmarks to SciPy using `Airspeed Velocity`_.


.. _toc-cython:

.. _compiled-code:

Compiled Code
-------------
- :ref:`adding-cython` - Extending and compiling Python code with `Cython`_ can significantly improve its performance. This document helps you get started.
- :ref:`other-languages` discusses the use of C, C++, and Fortran code in SciPy.

.. _Scipy Development Workflow: https://youtu.be/HgU01gJbzMY

.. _Sphinx: http://www.sphinx-doc.org/en/master/

.. _Airspeed Velocity: https://asv.readthedocs.io/en/stable/

.. _Cython: https://cython.org/

.. |*| replace:: \ :sup:`*` \
