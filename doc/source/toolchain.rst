Toolchain Roadmap
=================

The use of the SciPy library requires (or optionally depends upon) several
other libraries in order to operate, the main dependencies being Python
and NumPy. It requires a larger collection of libraries and tools in order
to build the library or to build the documentation.

Of course, the tooling and libraries are themselves not static.
This document aims to provide a guide as to how SciPy's use of
these dynamic dependencies will proceed over time.

SciPy aims to be compatible with a number of releases of its dependent
libraries and tools. Forcing the user base to other components for upgrade
for every release would greatly diminish the value of SciPy. However,
maintaining backwards compatibility with very old tooling/libraries
imposes limitations on which newer functionalities and capabilities
can be incorporated.
SciPy takes a somewhat conservative approach, maintaining compatibility with
several major releases of Python and NumPy on the major platforms.
(That may in and of itself impose further restrictions. See the C Compilers
section for an example.)


- First and foremost, SciPy is a Python project, hence it requires a Python environment.
- BLAS and LAPACK numerical libraries need to be installed.
- Compilers for C, C++, Cython, and Fortran code are needed.
- The Python environment needs the ``NumPy`` package to be installed.
- Testing requires the ``pytest`` Python package.
- Building the documentation requires the ``matplotlib``, Sphinx packages, as well as a LaTeX installation.

The tooling used to build CPython has some implications for the tooling used
in building SciPy.
It also has implications for the examples used in the
documentation (e.g., docstrings for functions),
as these examples can only use functionality present in all supported configurations.


Building SciPy
--------------

Python Versions
^^^^^^^^^^^^^^^

SciPy is compatible with several versions of Python.  When dropping support for
older Python versions, SciPy takes guidance from NEP 29 [10]_.  Python 2.7
support was dropped for SciPy releases numbered 1.3 and above but is still
available in release 1.2.x, which is a long-term support release [1]_, [2]_.

================  =======================================================================
 Date             Pythons supported
================  =======================================================================
 2018              Py2.7, Py3.4+ (SciPy 1.2.x is the last release to support Python 2.7)
 2019              Py3.5+ (but Py2.7-specific code not removed)
 2020              Py3.6+ (removal of Py2.7-specific code permitted)
================  =======================================================================

NumPy
^^^^^

SciPy depends on NumPy but releases of SciPy are not tied to releases of NumPy.
SciPy attempts to be compatible with at least the 4 previous releases of NumPy.
In particular, SciPy cannot rely on features of just the latest NumPy, but
needs to be written using what is common in all of those 4 releases. [1]_, [3]_.

The table shows the NumPy versions suitable for each major Python version
(for SciPy 1.3.x unless otherwise stated).

=================  ========================    ===========================
 Python             Minimum NumPy version       Maximum NumPy version
=================  ========================    ===========================
2.7 (SciPy 1.2)      1.8.2                      1.16.x
3.5 (SciPy 1.4)      1.13.3                     1.18.x
3.6 (SciPy 1.5)      1.14.5                     1.19.x
3.7                  1.16.5                     >= 1.20.x
3.8                  1.17.3                     >= 1.20.x
3.9                  1.19.3                     >= 1.20.x
=================  ========================    ===========================


C Compilers
^^^^^^^^^^^

SciPy is compatible with most modern C compilers (in particular ``clang``).
However, CPython on Windows is built with specific versions of the Microsoft
Visual C++ compiler [7]_, [8]_, [9]_, as is the corresponding build of SciPy.
This has implications for the C language standards that can be supported [6]_.
Starting from MS Visual Studio 16.8, C11/C17 is supported [11]_ (without the
C11 optional features [12]_ like atomics, threading, VLAs & complex types),
though Windows builds of CPython have not yet upgraded this far.

===================   ==============   ===================
CPython               MS Visual C++    C Standard
===================   ==============   ===================
2.7, 3.0, 3.1, 3.2       9.0           C90
3.3, 3.4                10.0           C90 & some of C99
3.5, 3.6                14.0           C90 & most of C99
3.7, 3.8, 3.9           15.7           C90 & most of C99
===================   ==============   ===================



C and C++ Language Standards
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

C and C++ language standards for SciPy are generally guidelines
rather than official decisions. This is particularly true of
attempting to predict adoption timelines for newer standards.

================  =======================================================================
 Date              C Standard
================  =======================================================================
 <= 2018           C90
 2019              C90 for old code, may consider C99 for new
 2020              C99
 2020              C++11
 2021              C++14, C++17
 ?                 C11, C17
 ?                 C++20
================  =======================================================================

The use of MS Visual Studio 9.0 (which doesn't have support C99)
to build Python2.7 has meant that C code in SciPy has had to conform
to the earlier C90 standard for the language and standard library.
With the dropping of Python 2.7 for SciPy 1.3.x, the C90 restriction is no
longer imposed by compilers. For GCC version < 5, an explicit ``-std=c99``
may have to be added by the user if C99 features are used in SciPy code.
*Note: even though C99 has been a standard for 20 years, experience has shown
that not all features are supported equally well across all platforms.*

C17 (occasionally called C18) is a bug fix for C11, so C11 may be skipped entirely.
Microsoft has taken very long to achieve conformance to C99/C11/C17, but as soon as CPython
is built with Visual Studio 16.8 or newer (see above), it will be possible to use C17
(though optional C11 features like atomics & threading are so far not supported in MSVC).


In practice, the C++ feature set that can be used is limited by the
availability in the MS VisualStudio versions that SciPy needs to support.
Since dropping support for Python 2.7, C++11 can be used universally, and
since dropping support for Python 3.6, the same is true also for C++14 & C++17.
This is because the oldest still required version of MS Visual Studio
(Visual Studio 15.7 <-> MSVC 19.14, see [8]_) has effectively full support, see [4]_.
Compiler support for C++20 is still under heavy development.

.. note::

    Developer Note: Some C99 features would be useful for scientific
    programming, in particular better support of IEEE 754 [5]_.
    SciPy has a small include file ``scipy/_lib/_c99compat.h`` which
    provides access to a few functions. Use in conjunction
    with ``<numpy/npy_math.h>``.

    ========================================= ========================================================
     Feature                                  Workaround
    ========================================= ========================================================
    ``isnan()``, ``isinf()``, ``isfinite()``  Use ``sc_isnan()``, ``sc_isinf()``, ``sc_isfinite()``
    ``NAN``                                   Use ``NPY_NAN`` (it is *almost* equivalent)
    inline functions                          Make static functions and place in an include .h file
    mid-block variable declarations           Declare variables at the top of the block
    ========================================= ========================================================


Fortran Compilers
^^^^^^^^^^^^^^^^^

Generally, any well-maintained compiler is likely suitable and can be
used to build SciPy.

======== ==================
 Tool     Version
======== ==================
gfortran   >= 4.8.0
ifort     A recent version
flang     A recent version
======== ==================


Cython Compiler
^^^^^^^^^^^^^^^

SciPy always requires a recent Cython compiler.

======== ============ ===============
 Tool    Tool Version  SciPy version
======== ============ ===============
Cython     >= 0.29.13  1.4.1
Cython     >= 0.29.18  1.5.0
======== ============ ===============



Other Libraries
^^^^^^^^^^^^^^^

Any library conforming to the BLAS/LAPACK interface may be used.
OpenBLAS, ATLAS, MKL, BLIS, and reference Netlib libraries are known to work.

=============== =====================================================
 Library           Minimum version
=============== =====================================================
LAPACK           3.4.1
BLAS             A recent version of OpenBLAS, MKL or ATLAS.
                 The Accelerate BLAS is no longer supported.
=============== =====================================================


There are some additional optional dependencies.

=============== ======== ==========================================
 Library        Version   URL
=============== ======== ==========================================
mpmath          Recent    http://mpmath.org/
scikit-umfpack  Recent    https://pypi.org/project/scikit-umfpack/
=============== ======== ==========================================


Moreover, Scipy supports interaction with other libraries. The test suite
has additional compatibility tests that are run when these are installed:

=========================  ========  ====================================
 Tool                      Version    URL
=========================  ========  ====================================
pydata/sparse              Recent     https://github.com/pydata/sparse/
=========================  ========  ====================================


Testing and Benchmarking
--------------------------

Testing and benchmarking require recent versions of:

=========================  ========  ====================================
 Tool                      Version    URL
=========================  ========  ====================================
pytest                     Recent     https://docs.pytest.org/en/latest/
asv (airspeed velocity)    Recent     https://asv.readthedocs.io/
=========================  ========  ====================================


Building the Documentation
--------------------------

==========   =================================================
 Tool        Version
==========   =================================================
Sphinx       Whatever recent versions work. >= 2.0.
numpydoc     Whatever recent versions work. >= 0.8.0.
matplotlib   Generally suggest >= 2.0.
LaTeX        A recent distribution, such as ``TeX Live 2016``.
==========   =================================================

[The ``numpydoc`` package is also used, but that is currently
packaged in ``doc/sphinxext``.]


.. note::

    Developer Note: The versions of ``numpy`` and ``matplotlib`` required have
    implications for the examples in Python docstrings.
    Examples must be able to be executed both in the environment used to
    build the documentation,
    as well as with any supported versions of ``numpy/matplotlib`` that
    a user may use with this release of SciPy.


Packaging
---------

A Recent version of:

=============  ========  =============================================
 Tool          Version    URL
=============  ========  =============================================
setuptools     Recent     https://pypi.org/project/setuptools/
wheel          Recent     https://pythonwheels.com
multibuild     Recent     https://github.com/matthew-brett/multibuild
=============  ========  =============================================

:ref:`making-a-release` and :ref:`distributing-a-release` contain information on
making and distributing a SciPy release.

References
----------

.. [1] https://docs.scipy.org/doc/scipy/reference/release.1.2.0.html
.. [2] https://python3statement.org
.. [3] https://docs.scipy.org/doc/numpy/release.html
.. [4] https://en.cppreference.com/w/cpp/compiler_support
.. [5] https://en.wikipedia.org/wiki/IEEE_754-1985
.. [6] https://blogs.msdn.microsoft.com/vcblog/2013/07/19/c99-library-support-in-visual-studio-2013/
.. [7] https://pythondev.readthedocs.io/windows.html#python-and-visual-studio-version-matrix
.. [8] https://en.wikipedia.org/wiki/Microsoft_Visual_C%2B%2B#Internal_version_numbering
.. [9] https://wiki.python.org/moin/WindowsCompilers
.. [10] https://numpy.org/neps/nep-0029-deprecation_policy.html
.. [11] https://devblogs.microsoft.com/cppblog/c11-and-c17-standard-support-arriving-in-msvc/
.. [12] https://en.wikipedia.org/wiki/C11_%28C_standard_revision%29#Optional_features
