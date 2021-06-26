.. _toolchain-roadmap:

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
- Compilers for C, C++, Fortran code are needed, as well as for Cython & Pythran (the latter is opt-out currently)
- The Python environment needs the ``NumPy`` package to be installed.
- Testing requires the ``pytest`` Python package.
- Building the documentation requires the ``matplotlib``, Sphinx packages along with PyData theme,
  as well as a LaTeX installation.

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
 2021              Py3.7+
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


Compilers
^^^^^^^^^

Building SciPy requires compilers for C, C++, Fortran, as well as the
python transpilers Cython and Pythran (the latter is an opt-out dependency
as of version 1.7.0).

To maintain compatibility with a large number of platforms & setups, especially
where using the official wheels (or other distribution channels like Anaconda
or conda-forge) is not possible, SciPy keeps compatibility with old compilers.

Official Builds
~~~~~~~~~~~~~~~

Currently, SciPy wheels are being built as follows:

================  ========================  ===========================  ==============================
 Platform          Azure Base Image [14]_    Compilers                    Comment
================  ========================  ===========================  ==============================
Linux (nightly)    ``ubuntu-18.04``          GCC 4.8                      See ``azure-pipelines.yml``
Linux (release)    ``ubuntu-18.04``          GCC 7.5                      Built in separate repo [15]_
OSX                ``macOS-10.14``           LLVM 11.0                    Built in separate repo [15]_
Windows            ``VS2017-Win2016``        Visual Studio 2017 (15.9)    See ``azure-pipelines.yml``
================  ========================  ===========================  ==============================

Note that the OSX wheels additionally vendor gfortran 4.8, see [15]_.


C Compilers
~~~~~~~~~~~

SciPy is compatible with most modern C compilers (in particular ``clang``).
In addition to concerns about compatibility with non-standard platforms,
there was a long-standing restriction that Windows builds of SciPy had to use
the same version of the Microsoft Visual C++ compiler as were used for CPython
itself, for reasons of ABI-compatibility [6]_, [7]_, [8]_, [9]_.

With the introduction of the "Universal C Runtime" [16]_ since the release of
Visual Studio 2015, this restriction has been lifted. For more context, see the
explanations by Steve Dower (member of the CPython-on-Windows core developers)
on this topic [17]_.

The use of MS Visual Studio 9.0 (which doesn't have support C99)
to build Python 2.7 has meant that C code in SciPy has had to conform
to the earlier C90 standard for the language and standard library.
With the dropping of Python 2.7 for SciPy 1.3.x, the C90 restriction is no
longer imposed by compilers. For GCC version < 5, an explicit ``-std=c99``
may have to be added by the user if C99 features are used in SciPy code.

In terms of C language standards, it's relevant to note that C11 has optional
features [12]_ (e.g. atomics, threading), some of which (VLAs & complex types)
were mandatory in the C99 standard. C17 (occasionally called C18) can be
considered a bug fix for C11, so generally, C11 may be skipped entirely.

SciPy has been restricted in the use of more advanced language features by the
available compiler support, and Microsoft in particular has taken very long to
achieve conformance to C99/C11/C17, however starting from MS Visual Studio 16.8,
C11/C17 is supported [11]_ (though without the C11 optional features).

Therefore, using C features beyond C90 is contingent upon updating the windows
toolchain for SciPy, as well as checking compiler support for the desired feature
across all minimally supported compiler versions. In short:

===================   ==============   =============================================
CPython               MS Visual C++    C Standard
===================   ==============   =============================================
2.7, 3.0, 3.1, 3.2       9.0           C90
3.3, 3.4                10.0           C90 & some of C99
3.5, 3.6                14.0           C90 & most of C99
3.7, 3.8, 3.9           15.7           Dependent on MSVC version used to build SciPy
===================   ==============   =============================================


C and C++ Language Standards
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C and C++ language standards for SciPy are generally guidelines
rather than official decisions. This is particularly true of
attempting to predict adoption timelines for newer standards.

================  =======================================================================
 Date              C/C++ Standard
================  =======================================================================
 <= 2018           C90
 2019              C90 for old code, may consider C99 for new
 2020              C99
 2020              C++11
 2021              C++14
 ?                 C11, C17, C++17, C++20
================  =======================================================================

For C, C11/C17 support will be available as soon as the ``vmImage`` for
building SciPy is upgraded to ``windows-2019`` (which is compatible with
currently supported CPython versions and "just" needs to be executed). This is
because GCC & LLVM support all relevant C11 features with the oldest currently
used versions, and C17 is just a bugfix for C11, as mentioned above.

On the C++ side, since dropping support for Python 2.7, C++11 can be used
universally. For C++14, Windows is not a restriction anymore since Visual
Studio 15.9 (<-> _MSC_VER 19.16, see [8]_), has full support (same for C++17),
see [4]_. However, using C++14 still requires bumping the GCC minimal
requirement to 5.x and C++17 will require GCC >= 7 [4]_.
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
~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~

SciPy always requires a recent Cython compiler.

======== ============ ===============
 Tool    Tool Version  SciPy version
======== ============ ===============
Cython     >= 0.29.13  1.4.1
Cython     >= 0.29.18  1.5.0
======== ============ ===============


OpenMP support
^^^^^^^^^^^^^^

For various reasons [13]_, SciPy cannot be distributed with built-in OpenMP support.
When using the optional Pythran support, OpenMP-enabled parallel code can be
generated when building from source.

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

====================  =================================================
 Tool                 Version
====================  =================================================
Sphinx                Whatever recent versions work. >= 2.0.
PyData Sphinx theme   Whatever recent versions work. >= 0.6.1.
Sphinx-Panels         Whatever recent versions work. >= 0.5.2.
numpydoc              Whatever recent versions work. >= 0.8.0.
matplotlib            Generally suggest >= 2.0.
LaTeX                 A recent distribution, such as ``TeX Live 2016``.
====================  =================================================

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
.. [13] https://github.com/scipy/scipy/issues/10239
.. [14] https://docs.microsoft.com/en-us/azure/devops/pipelines/agents/hosted
.. [15] https://github.com/MacPython/scipy-wheels
.. [16] https://docs.microsoft.com/en-gb/cpp/windows/universal-crt-deployment
.. [17] https://discuss.python.org/t/toolchain-upgrade-on-windows/6377/4
