.. _contributor-toc:

=============================
Contributor Table of Contents
=============================

This table of contents is designed to help you quickly find the information you need about SciPy development after you've reviewed the introductory material in :ref:`hacking`. If you're new to this and want to start coding ASAP, you've found the right place. Start with the bold links (in order), then expand your knowledge with the rest.

- :ref:`dev-env` - how to set up and maintain a development environment, including installing compilers and SciPy dependencies, creating a personal fork of the SciPy repository on GitHub, using git to manage a local repository with development branches, performing an in-place build of SciPy, and creating a virtual environment that adds this development version of SciPy to the Python path
- :ref:`editing-scipy` - how to edit SciPy Python code, with tips on finding which module contains SciPy functionality to be edited, adding new modules to SciPy, and complying with PEP8 style standards
- :ref:`unit-tests` - how to write and run unit tests for SciPy with the pytest framework
- :ref:`docs` - how to write reStructuredText documentation that complies with docstring standards, build documentation locally, and view documentation built during continuous integration checks
- :ref:`benchmarking` - how to benchmark code with AirSpeed Velocity
- :ref:`cython` - how to add fast, compiled code to SciPy

.. _dev-env:

Development Environment
-----------------------
- :ref:`quickstart-mac` presents a step-by-step process for setting up a convenient SciPy development environment on macOS.
- :ref:`building` - If you don't have macOS, try these instructions to help you build SciPy on your operating system.
- :ref:`recommended-development-setup` includes additional notes about the development setup. All of this information is contained elsewhere, but it is retained as a legacy document.

.. _editing-scipy:

Editing SciPy
-------------
- :ref:`development-workflow` lays out what to do after your development environment is set up.
- `SciPy Development Workflow`_ is also the name of a five-minute video example of fixing a bug and submitting a pull request.
- `PEP8 and SciPy`_ gives some tips for ensuring that your code is PEP8 compliant.
- :ref:`git-development` is a guide to using ``git``, the distributed version-control system used to manage the changes made to SciPy code from around the world. It is written in the context of NumPy, but working with SciPy is essentially identical - just substitute the text SciPy/``scipy`` wherever NumPy/``numpy`` appears.
 
.. _unit-tests:

Unit Tests
----------

.. _docs:

Documentation
-------------

.. _benchmarks:

Benchmarks
----------

.. _cython:

Cython
------

.. _Scipy Development Workflow: https://youtu.be/HgU01gJbzMY

.. _PEP8 and SciPy: https://github.com/scipy/scipy/wiki/PEP8-and-SciPy
