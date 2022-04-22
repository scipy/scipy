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

.. note::

    You can skip these steps if you are using ``conda``, as these dependencies
    will be installed automatically.

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

    .. tabs::

      .. tab:: Debian/Ubuntu Linux

        .. note::

            These instructions have been tested on Ubuntu Linux 16.04, 18.04, and
            20.04.

        Python should be available in your system via the ``python3`` command. To
        install the remaining system-level dependencies, run::

          sudo apt install -y gcc g++ gfortran libopenblas-dev liblapack-dev pkg-config
          sudo apt install -y python3-pip python3-dev

        Alternatively, you can do::

          sudo apt build-dep scipy

        This command installs whatever is needed to build SciPy, with the
        advantage that new dependencies or updates to required versions are
        handled by the package managers.

        See also :ref:`ubuntu-guide`.


      .. tab:: Fedora

        To install SciPy build requirements, you can do::

          sudo dnf install gcc-gfortran python3-devel openblas-devel lapack-devel

        Alternatively, you can do::

          sudo dnf builddep scipy

        This command installs whatever is needed to build SciPy, with the
        advantage that new dependencies or updates to required versions are
        handled by the package managers.

      .. tab:: CentOS/RHEL

        To install SciPy build requirements, you can do::

          sudo yum install gcc-gfortran python3-devel openblas-devel lapack-devel

        Alternatively, you can do::

          sudo yum-builddep scipy

        This command installs whatever is needed to build SciPy, with the
        advantage that new dependencies or updates to required versions are
        handled by the package managers.

      .. tab:: Arch

        To install SciPy build requirements, you can do::

          sudo pacman -S gcc-fortran openblas

    All further work should proceed in a virtual environment. Popular options
    include the standard library ``venv`` module or a separate ``virtualenv``
    package. 

    * The `Cython <https://cython.org/>`__ and
      `Pythran <https://pythran.readthedocs.io>`__ ahead-of-time compilers are also
      necessary, as is ``pybind11``. It is recommended to install these packages
      with ``pip``, because it is possible (even likely) that you need newer
      versions of these packages than the ones that are available in your Linux
      distribution.

    If you are using conda, these dependencies can be installed in the conda
    environment itself. See :ref:`conda-guide` for more details.

  .. tab:: OSX

    Install Apple Developer Tools. An easy way to do this is to
    `open a terminal window <https://blog.teamtreehouse.com/introduction-to-the-mac-os-x-command-line>`_,
    enter the command

    ::

        xcode-select --install

    and follow the prompts. Apple Developer Tools includes
    `git <https://git-scm.com/>`_, the software we need to download and manage the
    SciPy source code. See also :ref:`build-osx`.

  .. tab:: Windows

    See :ref:`build-windows`.

Detailed instructions
---------------------

.. toctree::
   :maxdepth: 1

   osx
   windows
   meson
   meson_advanced
   building_faq


Reference for build options
===========================

SciPy has several tunable build-time options, which can be set.

.. warning::

    This content is for the old `numpy.distutils`-based build and doesn't apply
    to the Meson build (i.e., when building with ``python dev.py``).

- ``site.cfg``: build-time library configuration file, see
  ``site.cfg.example`` for details.

- Environment variables ``NPY_LAPACK_ORDER``, ``NPY_BLAS_ORDER``, ``OPENBLAS``,
  ``ATLAS``, etc., also controlling library configuration.
  See `Numpy documentation <numpy-blasdoc>`_ for more details.

- Environment variable ``NPY_USE_BLAS_ILP64=1``: build using 64-bit
  integer size (ILP64) BLAS+LAPACK libraries.

  Note that even when this is set, SciPy requires *also* 32-bit
  integer size (LP64) BLAS+LAPACK libraries to be available and
  configured. This is because only some components in SciPy make use
  of the 64-bit capabilities.

.. _numpy-blasdoc: https://numpy.org/devdocs/user/building.html#accelerated-blas-lapack-libraries
