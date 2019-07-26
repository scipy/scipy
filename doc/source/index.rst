SciPy
=====

:Release: |release|
:Date: |today|

SciPy (pronounced "Sigh Pie") is open-source software for mathematics,
science, and engineering.

.. toctree::
   :maxdepth: 1

   install_upgrade
   api
   release

Tutorial
--------

Tutorials with worked examples and background information for most SciPy
submodules.

.. toctree::
   :maxdepth: 2

   tutorial/index.rst

Developer Documentation
-----------------------

If you're interested in contributing to SciPy, start here:

.. toctree::
   :maxdepth: 1

   dev/conduct/code_of_conduct
   hacking
   dev/contributor/contributor_toc

To get an overview of where help or new features are desired or planned, see
the roadmap:

.. toctree::
   :maxdepth: 1

   roadmap
   roadmap-detailed
   toolchain

For a more detailed look at how the SciPy project works:

.. toctree::
   :maxdepth: 1

   dev/core-dev/index
   dev/governance/governance

.. This toctree defines previous/next for contributor guide documents
.. toctree::
   :hidden:

   dev/contributor/quickstart_mac
   dev/contributor/quickstart_ubuntu
   dev/contributor/development_workflow
   dev/contributor/pep8
   dev/contributor/rendering_documentation
   dev/contributor/runtests
   dev/contributor/benchmarking
   dev/contributor/cython

.. These files are not intended to be in any toctree. because they have not
   been maintained.They should only be reached via the contributor guide if
   they are specifically sought, not via next/previous links.
..   building/index
..   dev/gitwash/gitwash
..   dev/contributor/recommended_development_setup
..   dev/contributor/compiled_code


.. module:: scipy

API Reference
-------------

The exact API of all functions and classes, as given by the docstrings. The API
documents expected types and allowed features for all functions, and all
parameters available for the algorithms.

.. toctree::
   :maxdepth: 1

   cluster
   constants
   fft
   fftpack
   integrate
   interpolate
   io
   linalg
   misc
   ndimage
   odr
   optimize
   signal
   sparse
   sparse.linalg
   sparse.csgraph
   spatial
   special
   stats
   stats.mstats
   ccallback
