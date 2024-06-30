.. _building-from-source:

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

SciPy uses compiled code for speed, which means you need compilers and some
other system-level (i.e, non-Python / non-PyPI) dependencies to build it on
your system.

.. note::

    If you are using Conda, you can skip the steps in this section - with the
    exception of installing compilers for Windows or the Apple Developer Tools
    for macOS. All other dependencies will be installed automatically by the
    ``mamba env create -f environment.yml`` command.

.. tab-set::

  .. tab-item:: Linux
    :sync: linux

    If you want to use the system Python and ``pip``, you will need:

    * C, C++, and Fortran compilers (typically ``gcc``, ``g++``, and ``gfortran``).

    * Python header files (typically a package named ``python3-dev`` or
      ``python3-devel``)

    * BLAS and LAPACK libraries. `OpenBLAS <https://github.com/xianyi/OpenBLAS/>`__
      is the SciPy default; other variants include
      `ATLAS <http://math-atlas.sourceforge.net/>`__ and
      `MKL <https://software.intel.com/en-us/intel-mkl>`__.

    * ``pkg-config`` for dependency detection.

    .. tab-set::

      .. tab-item:: Debian/Ubuntu Linux

        To install SciPy build requirements, you can do::

          sudo apt install -y gcc g++ gfortran libopenblas-dev liblapack-dev pkg-config python3-pip python3-dev

        Alternatively, you can do::

          sudo apt build-dep scipy

        This command installs whatever is needed to build SciPy, with the
        advantage that new dependencies or updates to required versions are
        handled by the package managers.

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

  .. tab-item:: macOS
    :sync: macos

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

    To allow the build tools to find OpenBLAS, you must run::

        brew info openblas | grep PKG_CONFIG_PATH

    This will give you a command starting with ``export PKG_CONFIG_PATH=``, which
    you must run.

    .. note::

        As of SciPy 1.14.0, we have added support for the Accelerate library
        for BLAS and LAPACK. It requires macOS 13.3 or greater. To build with
        Accelerate instead of OpenBLAS, see :ref:`blas-lapack-selection`.

  .. tab-item:: Windows
    :sync: windows

    A compatible set of C, C++ and Fortran compilers is needed to build SciPy.
    This is trickier on Windows than on other platforms, because MSVC does not
    support Fortran, and gfortran and MSVC can't be used together. You will
    need one of these sets of compilers:

    1. Mingw-w64 compilers (``gcc``, ``g++``, ``gfortran``) - *recommended,
       because it's easiest to install and is what we use for SciPy's own CI
       and binaries*
    2. MSVC + Intel Fortran (``ifort``)
    3. Intel compilers (``icc``, ``ifort``)

    Compared to macOS and Linux, building SciPy on Windows is a little more
    difficult, due to the need to set up these compilers. It is not possible to
    just call a one-liner on the command prompt as you would on other
    platforms.

    First, install Microsoft Visual Studio - the 2019 Community Edition or any
    newer version will work (see the
    `Visual Studio download site <https://visualstudio.microsoft.com/downloads/>`__).
    This is needed even if you use the MinGW-w64 or Intel compilers, in order
    to ensure you have the Windows Universal C Runtime (the other components of
    Visual Studio are not needed when using Mingw-w64, and can be deselected if
    desired, to save disk space).

    .. tab-set::

      .. tab-item:: MinGW-w64

        There are several sources of binaries for MinGW-w64. We recommend the
        RTools versions, which can be installed with Chocolatey (see
        Chocolatey install instructions `here <https://chocolatey.org/install>`_)::

            choco install rtools -y --no-progress --force --version=4.0.0.20220206

        In case of issues, we recommend using the exact same version as used
        in the `SciPy GitHub Actions CI jobs for Windows
        <https://github.com/scipy/scipy/blob/main/.github/workflows/windows.yml>`__.

      .. tab-item:: MSVC

        The MSVC installer does not put the compilers on the system path, and
        the install location may change. To query the install location, MSVC
        comes with a ``vswhere.exe`` command-line utility. And to make the
        C/C++ compilers available inside the shell you are using, you need to
        run a ``.bat`` file for the correct bitness and architecture (e.g., for
        64-bit Intel CPUs, use ``vcvars64.bat``).

        For detailed guidance, see `Use the Microsoft C++ toolset from the command line
        <https://learn.microsoft.com/en-us/cpp/build/building-on-the-command-line?view=msvc-170>`__.

      .. tab-item:: Intel

        Similar to MSVC, the Intel compilers are designed to be used with an
        activation script (``Intel\oneAPI\setvars.bat``) that you run in the
        shell you are using. This makes the compilers available on the path.
        For detailed guidance, see
        `Get Started with the Intel® oneAPI HPC Toolkit for Windows
        <https://www.intel.com/content/www/us/en/docs/oneapi-hpc-toolkit/get-started-guide-windows/2023-1/overview.html>`__.

    .. note::

        Compilers should be on the system path (i.e., the ``PATH`` environment
        variable should contain the directory in which the compiler executables
        can be found) in order to be found, with the exception of MSVC which
        will be found automatically if and only if there are no other compilers
        on the ``PATH``. You can use any shell (e.g., Powershell, ``cmd`` or
        Git Bash) to invoke a build. To check that this is the case, try
        invoking a Fortran compiler in the shell you use (e.g., ``gfortran
        --version`` or ``ifort --version``).

    .. warning::

        When using a conda environment it is possible that the environment
        creation will not work due to an outdated Fortran compiler. If that
        happens, remove the ``compilers`` entry from ``environment.yml`` and
        try again. The Fortran compiler should be installed as described in
        this section.


Building SciPy from source
--------------------------

If you want to only install SciPy from source once and not do any development
work, then the recommended way to build and install is to use ``pip``.
Otherwise, conda is recommended.

.. note::

    If you don't have a conda installation yet, we recommend using
    Mambaforge_; any conda flavor will work though.

Building from source to use SciPy
`````````````````````````````````

.. tab-set::

  .. tab-item:: Conda env
    :sync: conda

    If you are using a conda environment, ``pip`` is still the tool you use to
    invoke a from-source build of SciPy. It is important to always use the
    ``--no-build-isolation`` flag to the ``pip install`` command, to avoid
    building against a ``numpy`` wheel from PyPI. In order for that to work you
    must first install the remaining build dependencies into the conda
    environment::

      # Either install all SciPy dev dependencies into a fresh conda environment
      mamba env create -f environment.yml

      # Or, install only the required build dependencies
      mamba install python numpy cython pythran pybind11 compilers openblas meson-python pkg-config

      # To build the latest stable release:
      pip install scipy --no-build-isolation --no-binary scipy

      # To build a development version, you need a local clone of the SciPy git repository:
      git clone https://github.com/scipy/scipy.git
      cd scipy
      git submodule update --init
      pip install . --no-build-isolation

  .. tab-item:: Virtual env or system Python
    :sync: pip

    ::

      # To build the latest stable release:
      pip install scipy --no-binary scipy

      # To build a development version, you need a local clone of the SciPy git repository:
      git clone https://github.com/scipy/scipy.git
      cd scipy
      git submodule update --init
      pip install .



.. _the-dev-py-interface:

Building from source for SciPy development
``````````````````````````````````````````

If you want to build from source in order to work on SciPy itself, first clone
the SciPy repository::

      git clone https://github.com/scipy/scipy.git
      cd scipy
      git submodule update --init

Then you want to do the following:

1. Create a dedicated development environment (virtual environment or conda
   environment),
2. Install all needed dependencies (*build*, and also *test*, *doc* and
   *optional* dependencies), 
3. Build SciPy with our ``dev.py`` developer interface.

Step (3) is always the same, steps (1) and (2) are different between conda and
virtual environments:

.. tab-set::

  .. tab-item:: Conda env
    :sync: conda

    To create a ``scipy-dev`` development environment with every required and
    optional dependency installed, run::

        mamba env create -f environment.yml
        mamba activate scipy-dev

  .. tab-item:: Virtual env or system Python
    :sync: pip

    .. note::

       There are many tools to manage virtual environments, like ``venv``,
       ``virtualenv``/``virtualenvwrapper``, ``pyenv``/``pyenv-virtualenv``,
       Poetry, PDM, Hatch, and more. Here we use the basic ``venv`` tool that
       is part of the Python stdlib. You can use any other tool; all we need is
       an activated Python environment.

    Create and activate a virtual environment in a new directory named ``venv`` (
    note that the exact activation command may be different based on your OS and shell
    - see `"How venvs work" <https://docs.python.org/3/library/venv.html#how-venvs-work>`__
    in the ``venv`` docs).

    .. tab-set::

      .. tab-item:: Linux
        :sync: linux

        ::

          python -m venv venv
          source venv/bin/activate

      .. tab-item:: macOS
        :sync: macos

        ::

          python -m venv venv
          source venv/bin/activate

      .. tab-item:: Windows
        :sync: windows

        ::

          python -m venv venv
          .\venv\Scripts\activate

    Then install the Python-level dependencies (see ``pyproject.toml``) from
    PyPI with::

       # All dependencies
       python -m pip install -r requirements/all.txt

       # Alternatively, you can install just the dependencies for certain
       # development tasks:

       # Build and dev dependencies (for `python dev.py {build, lint, mypy}`)
       python -m pip install -r requirements/build.txt -r requirements/dev.txt

       # Doc dependencies (for `python dev.py {doc, refguide-check}`)
       python -m pip install -r requirements/doc.txt

       # Test dependencies (for `python dev.py {test, bench, refguide-check}`)
       python -m pip install -r requirements/test.txt

To build SciPy in an activated development environment, run::

    python dev.py build

This will install SciPy inside the repository (by default in a
``build-install`` directory). You can then run tests (``python dev.py test``),
drop into IPython (``python dev.py ipython``), or take other development steps
like build the html documentation or running benchmarks. The ``dev.py``
interface is self-documenting, so please see ``python dev.py --help`` and
``python dev.py <subcommand> --help`` for detailed guidance.


.. admonition:: IDE support & editable installs

    While the ``dev.py`` interface is our recommended way of working on SciPy,
    it has one limitation: because of the custom install location, SciPy
    installed using ``dev.py`` will not be recognized automatically within an
    IDE (e.g., for running a script via a "run" button, or setting breakpoints
    visually). This will work better with an *in-place build* (or "editable
    install").

    Editable installs are supported. It is important to understand that **you
    may use either an editable install or dev.py in a given repository clone,
    but not both**. If you use editable installs, you have to use ``pytest``
    and other development tools directly instead of using ``dev.py``.

    To use an editable install, ensure you start from a clean repository (run
    ``git clean -xdf`` if you've built with ``dev.py`` before) and have all
    dependencies set up correctly as described higher up on this page. Then
    do::

        # Note: the --no-build-isolation is important! meson-python will
        # auto-rebuild each time SciPy is imported by the Python interpreter.
        pip install -e . --no-build-isolation

        # To run the tests for, e.g., the `scipy.linalg` module:
        pytest scipy/linalg

    When making changes to SciPy code, including to compiled code, there is no
    need to manually rebuild or reinstall. When you run ``git clean -xdf``,
    which removes the built extension modules, remember to also uninstall SciPy
    with ``pip uninstall scipy``.

    See the meson-python_ documentation on editable installs for more details
    on how things work under the hood.


Customizing builds
------------------

.. toctree::
   :maxdepth: 1

   compilers_and_options
   blas_lapack
   cross_compilation
   redistributable_binaries


Background information
----------------------

.. toctree::
   :maxdepth: 1

   understanding_meson
   introspecting_a_build
   distutils_equivalents


.. _Mambaforge: https://github.com/conda-forge/miniforge#mambaforge
.. _meson-python: https://mesonbuild.com/meson-python/
