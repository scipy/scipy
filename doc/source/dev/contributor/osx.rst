.. _build-osx:

Building from source on Mac
===========================

.. note::

    These instructions have been tested on macOS 11.1.

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
conda environment, instead of manually.

Installing without using conda
------------------------------

Apple ships its own version of Python with OS X. However, we *strongly*
recommend installing the `official Python distribution
<https://www.python.org/downloads/>`__.

Alternatively, use Python from one of the OS X package managers (Homebrew,
MacPorts, Fink).

Compilers (C/C++/FORTRAN/Cython/Pythran)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Though virtually any C/C++ compiler may be used with SciPy, on macOS we
recommend using the default Clang C/C++ compiler. It can be installed as
an Xcode command line tool (this is what the ``xcode-select`` command given
higher up has already done).

The only thing missing is a Fortran compiler.  We recommend you
use ``gfortran``, and install it with `Homebrew <https://brew.sh/>`__.

BLAS/LAPACK Installation
~~~~~~~~~~~~~~~~~~~~~~~~

You will also need to install a library providing the BLAS and LAPACK
interfaces. OpenBLAS, ATLAS and MKL all work. OpenBLAS can be installed
via `Homebrew <https://brew.sh/>`, and is the library we recommend using.

.. note::

    As of SciPy >=1.2.0, we do not support compiling against the system
    Accelerate library for BLAS and LAPACK. It does not support a sufficiently
    recent LAPACK interface.
