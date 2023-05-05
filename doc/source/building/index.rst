.. _building_from_source_ref:

Building from source
====================

.. note::

   If you are only trying to install SciPy, we recommend using binaries - see
   `Installation <https://scipy.org/install>`__ for details on that.

Building SciPy from source requires setting up system-level dependencies
(compilers, BLAS/LAPACK libraries, etc.) first, and then invoking a build. The
build may be done in order to install SciPy for local usage, develop SciPy
itself, or build redistributable binary packages. And it may be desired to
customize aspects of how the build is done. This guide will cover all these
aspects. In addition, it provides background information on how the SciPy build
works, and links to up-to-date guides for generic Python build & packaging
documentation that is relevant.


.. _system-level:

System-level dependencies
-------------------------

SciPy uses compiled code for speed, which means you might need extra
dependencies to build it on your system.

.. note::

    You can skip these steps if you are using ``conda``, as these dependencies
    will be installed automatically.

.. tab-set::

  .. tab-item:: Linux

    If you want to use the system Python and ``pip``, you will need:

    * C, C++, and Fortran compilers (typically ``gcc``, ``g++``, and ``gfortran``).

    * Python header files (typically a package named ``python3-dev`` or
      ``python3-devel``)

    * BLAS and LAPACK libraries. `OpenBLAS <https://github.com/xianyi/OpenBLAS/>`__
      is the SciPy default; other variants include
      `ATLAS <http://math-atlas.sourceforge.net/>`__ and
      `MKL <https://software.intel.com/en-us/intel-mkl>`__.

    .. tab-set::

      .. tab-item:: Debian/Ubuntu Linux

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


      .. tab-item:: Fedora

        To install SciPy build requirements, you can do::

          sudo dnf install gcc-gfortran python3-devel openblas-devel lapack-devel pkgconfig

        Alternatively, you can do::

          sudo dnf builddep scipy

        This command installs whatever is needed to build SciPy, with the
        advantage that new dependencies or updates to required versions are
        handled by the package managers.

      .. tab-item:: CentOS/RHEL

        To install SciPy build requirements, you can do::

          sudo yum install gcc-gfortran python3-devel openblas-devel lapack-devel pkgconfig

        Alternatively, you can do::

          sudo yum-builddep scipy

        This command installs whatever is needed to build SciPy, with the
        advantage that new dependencies or updates to required versions are
        handled by the package managers.

      .. tab-item:: Arch

        To install SciPy build requirements, you can do::

          sudo pacman -S gcc-fortran openblas pkgconf

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

  .. tab-item:: macOS

    Install Apple Developer Tools. An easy way to do this is to
    `open a terminal window <https://blog.teamtreehouse.com/introduction-to-the-mac-os-x-command-line>`_,
    enter the command::

        xcode-select --install

    and follow the prompts. Apple Developer Tools includes Git, the Clang C/C++
    compilers, and other development utilities that may be required.

    Do *not* use the macOS system Python. Instead, install Python
    with `the python.org installer <https://www.python.org/downloads/>`__ or
    with a package manager like Homebrew, MacPorts or Fink.

    The other system dependencies you need are a Fortran compiler, BLAS and
    LAPACK libraries, and pkg-config. They're easiest to install with
    `Homebrew <https://brew.sh/>`__::

        brew install gfortran openblas pkg-config

    .. note::

        As of SciPy >=1.2.0, we do not support compiling against the system
        Accelerate library for BLAS and LAPACK. It does not support a sufficiently
        recent LAPACK interface. This is planned to change in 2023, because macOS
        13.3 introduced a major upgrade to Accelerate which resolved all known
        issues.

  .. tab-item:: Windows

    A compatible set of C, C++ and Fortran compilers is needed to build SciPy.
    This is tricker on Windows than on other platforms, because MSVC does not
    support Fortran, and gfortran and MSVC can't be used together. You will
    need one of these sets of compilers:

    1. MSVC + Intel Fortran (``ifort``)
    2. Intel compilers (``icc``, ``ifort``)
    3. Mingw-w64 compilers (``gcc``, ``g++``, ``gfortran``)

    Compared to macOS and Linux, building SciPy on Windows is more difficult,
    largely due to the difficulty of setting up these compilers. It is not
    possible to just call a one-liner on the command prompt as you would on
    other platforms.

    First, install Microsoft Visual Studio - the 2019 Community Edition or any
    newer version will work (see the
    `Visual Studio download site <https://visualstudio.microsoft.com/downloads/>`__). 
    This is needed even if you use the MinGW-w64 or Intel compilers, in order
    to ensure you have the Windows Universal C Runtime.

    .. tab-set::

      .. tab-item:: MSVC

        TODO: MSVC-specific guidance

      .. tab-item:: Intel

        TODO: Intel-specific guidance

      .. tab-item:: MinGW-w64

        TODO: MinGW-w64 specific guidance

        It makes sense to use Windows ``cmd`` or Powershell for the the build as it is
        a more native tool. This requires placing the MinGW compilers on the path.
        Hence, make sure that the following folder (or the folder you have installed
        MSYS to) is on the system path variable sufficiently close to the top.

    .. note::

        Compilers should be on the system path (i.e., the ``PATH`` environment
        variable) in order to be found, with the exception of MSVC which will
        be found automatically if and only if there are no other compilers on
        the ``PATH``. You can use any shell (e.g., Powershell, ``cmd`` or Git
        Bash) to invoke a build. To check that this is the case, try invoking a
        Fortran compiler in the shell you use (e.g., ``gfortran --version`` or
        ``ifort --version``).


Building SciPy from source
--------------------------

If you want to only install SciPy from source once and not do any development
work, then the recommended way to build and install is to use ``pip``::

    # For the latest stable release:
    pip install scipy --no-binary scipy

    # For the latest development version, directly from GitHub:
    pip install https://github.com/scipy/scipy/archive/refs/heads/main.zip

    # If you have a local clone of the SciPy git repository:
    pip install .

If you want to build from source in order to work on SciPy itself, then use
our ``dev.py`` developer interface with::

    python dev.py build

For more details on developing with ``dev.py``, see :ref:`meson`.

Detailed instructions
---------------------

.. toctree::
   :maxdepth: 1

   blas_lapack

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
  See :ref:`NumPy documentation <numpy:accelerated-blas-lapack-libraries>`
  for more details.

- Environment variable ``NPY_USE_BLAS_ILP64=1``: build using 64-bit
  integer size (ILP64) BLAS+LAPACK libraries.

  Note that even when this is set, SciPy requires *also* 32-bit
  integer size (LP64) BLAS+LAPACK libraries to be available and
  configured. This is because only some components in SciPy make use
  of the 64-bit capabilities.
