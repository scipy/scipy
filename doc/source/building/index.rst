:orphan:

.. _building:

Building from sources
=====================

.. note::

   If you are only trying to install SciPy, see
   `Installation <https://scipy.org/install>`__.

Build instructions for different operating systems and an FAQ.

.. _system-level:

System-level dependencies
-------------------------

SciPy uses compiled code for speed, which means you might need extra
dependencies to build it on your system.

.. tabs::

  .. tab:: Linux
  
    If you want to use the system Python and ``pip``, you will need:

    * C, C++, and Fortran compilers (typically ``gcc``, ``g++``, and ``gfortran``).

    * Python header files (typically a package named ``python3-dev`` or
      ``python3-devel``)

    * BLAS and LAPACK libraries. `OpenBLAS <https://github.com/xianyi/OpenBLAS/>`__
      is the SciPy default; other variants include
      `ATLAS <http://math-atlas.sourceforge.net/>`__ and
      `MKL <https://software.intel.com/en-us/intel-mkl>`__. 

    Typically, you will want to install all of the above from packages supplied by
    your Linux distribution, as building them yourself is complicated. If you need
    to use specific BLAS/LAPACK libraries, you can do

    ::

        export BLAS=/path/to/libblas.so
        export LAPACK=/path/to/liblapack.so
        export ATLAS=/path/to/libatlas.so
        python setup.py ............

    * The `Cython <https://cython.org/>`__ and
      `Pythran <https://pythran.readthedocs.io>`__ ahead-of-time compilers are also
      necessary, as is ``pybind11``. It is recommended to install these packages
      with ``pip`` or ``conda``, because it is possible (even likely) that you need
      newer versions of these packages than the ones that are available in your
      Linux distribution.

    If you are using conda, these dependencies can be installed in the conda
    environment itself. See :ref:`conda-guide` for more details.

    See :ref:`build-linux` for specific details on installing these dependencies
    for popular Linux distros.

  .. tab:: OSX

    Install Apple Developer Tools. An easy way to do this is to
    `open a terminal window <https://blog.teamtreehouse.com/introduction-to-the-mac-os-x-command-line>`_,
    enter the command

    ::

        xcode-select --install

    and follow the prompts. Apple Developer Tools includes
    `git <https://git-scm.com/>`_, the software we need to download and manage the
    SciPy source code.

    We recommend :ref:`using conda <conda-guide>` to set up your environment on a
    Mac, as this will also allow you to install all necessary dependencies in a
    conda environment, instead of manually. If you prefer not using conda, see
    :ref:`build-osx`.

  .. tab:: Windows

    We recommend :ref:`using conda <conda-guide>` to set up your environment on
    Windows, as this will also allow you to install all necessary dependencies
    in a conda environment, instead of manually. If you prefer not using conda
    or need more details, see :ref:`build-windows`.

Detailed instructions
---------------------

.. toctree::
   :maxdepth: 1

   linux
   osx
   windows
   meson
   meson_advanced
   faq


Reference for build options
===========================

SciPy has several tunable build-time options, which can be set.

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
