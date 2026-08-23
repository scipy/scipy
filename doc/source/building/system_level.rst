System-level dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^

SciPy uses compiled code for speed, which means you need compilers and some
other system-level (i.e, non-Python / non-PyPI) dependencies to build it on
your system.

.. note::

    If you are using conda, you can skip the steps in this section - with the
    exception of installing compilers for Windows or the Apple Developer Tools
    for macOS. All other dependencies will be installed automatically by the
    ``conda env create -f environment.yml`` command.

.. tab-set::

  .. tab-item:: Linux
    :sync: linux

    If you want to use the system Python and ``pip``, you will need:

    * C and C++ compilers (typically ``gcc``, ``g++``).

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

          sudo apt install -y gcc g++ libopenblas-dev liblapack-dev pkg-config python3-pip python3-dev

        Alternatively, you can do::

          sudo apt build-dep scipy

        This command installs whatever is needed to build SciPy, with the
        advantage that new dependencies or updates to required versions are
        handled by the package managers.

      .. tab-item:: Fedora / RHEL & CentOS 8+

        To install SciPy build requirements, you can do::

          sudo dnf install gcc python3-devel openblas-devel lapack-devel pkgconfig

        Alternatively, you can do::

          sudo dnf builddep scipy

        This command installs whatever is needed to build SciPy, with the
        advantage that new dependencies or updates to required versions are
        handled by the package managers.

      .. tab-item:: CentOS & RHEL <=7

        To install SciPy build requirements, you can do::

          sudo yum install gcc python3-devel openblas-devel lapack-devel pkgconfig

        Alternatively, you can do::

          sudo yum-builddep scipy

        This command installs whatever is needed to build SciPy, with the
        advantage that new dependencies or updates to required versions are
        handled by the package managers.

      .. tab-item:: Arch

        To install SciPy build requirements, you can do::

          sudo pacman -S gcc openblas pkgconf

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

    The other system dependencies you need are BLAS and
    LAPACK libraries, and pkg-config. They're easiest to install with
    `Homebrew <https://brew.sh/>`__::

        brew install openblas pkg-config

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

    A compatible set of C and C++ compilers is needed to build SciPy.
    You will need one of these sets of compilers:

    1. Mingw-w64 compilers (``gcc``, ``g++``) - *recommended,
       because it's easiest to install and is what we use for SciPy's own CI
       and binaries*
    2. Clang-cl
    3. Intel compilers (``icc``)

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
        Git Bash) to invoke a build.
