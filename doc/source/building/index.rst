:orphan:

.. _building:

Building from sources
=====================

.. note::

   If you are only trying to install SciPy, see
   :doc:`../install_upgrade`.

Build instructions for different operating systems and an FAQ:

.. toctree::
   :maxdepth: 2

   linux
   windows
   macosx
   faq


Reference for build options
===========================

Scipy has several tunable build-time options, which can be set.

- ``site.cfg``: build-time library configuration file, see
  ``site.cfg.example`` for details.

- Environment variables ``NPY_LAPACK_ORDER``, ``NPY_BLAS_ORDER``, ``OPENBLAS``,
  ``ATLAS``, etc., also controlling library configuration.
  See `Numpy documentation <numpy-blasdoc>`_ for more details.

- Environment variable ``NPY_USE_BLAS_ILP64=1``: build using 64-bit
  integer size (ILP64) BLAS+LAPACK libraries.

  Note that even when this is set, Scipy requires *also* 32-bit
  integer size (LP64) BLAS+LAPACK libraries to be available and
  configured. This is because only some components in Scipy make use
  of the 64-bit capabilities.

.. _numpy-blasdoc: https://numpy.org/devdocs/user/building.html#accelerated-blas-lapack-libraries
