Toolchain and Roadmap
=====================
SciPy requires several tools in order to be built.

- First and foremost, SciPy is a Python project hence it requires a Python environment.
- The LAPACK and OpenBLAS numerical libraries need to be installed.
- Compilers for C, C++, Cython and Fortran code are needed.
- The Python environment needs the ``NumPy`` package to be installed.
- Testing requires the ``pytest`` Python package.
- Building the documentation requires Sphinx and ``numpydoc``.

The tooling used to build CPython has some implications for the tooling used in building SciPy.

Python Versions
---------------

SciPy is compatible with several versions of Python, and some
specific decisions are still under consideration, especially
with regard to future changes.

================  =================================================
 Date             Pythons supported
================  =================================================
 <= 2018           Py2.7, Py3.5 - Py3.7
 2019              Py3.5 - (but Py2.7-specific code not removed)
 2020              Py3.?- (removal of Py2.7-specific code permitted)
================  =================================================


C Compilers
-----------
SciPy is compatible with most modern C compilers.  However CPython on Windows is
built with specific versions of the Microsoft Visual C++ compiler [1]_, as is the
corresponding build of SciPy.  This has implications for the C language standards
that can be supported.

===================   ==============   ===================
CPython               MS Visual C++    C Standard
===================   ==============   ===================
3.5, 3.6                14.0           C90 & most of C99
3.3, 3.4                10.0           C90 & some of C99
2.7, 3.0, 3.1, 3.2       9.0           C90
===================   ==============   ===================


C Language Standards
--------------------
C and C++ language standards for SciPy are generally guidelines
rather than official decisions. This is particularly true of
attempting to predict adoption timelines for newer standards.

================  =========================================
 Date              C Standard
================  =========================================
 <= 2018           C90
 2019              C90 for old code, may allow C99 for new
 2020              C99
 No decision       C11
 No decision       C17, C18
================  =========================================

Some C99 features would be useful for scientific programming, in particular better support of
IEEE 754 [2]_.  Experience has shown that not all features are supported equally well across
all platforms. SciPy has a small include file ``scipy/_lib/_c99compat.h`` which provides
access to a few functions.  Use in conjunction with ``<numpy/npy_math.h>``.

================================  ========================================================
 Feature                           Workaround
================================  ========================================================
 isnan(), isinf(), isfinite()      Use sc_isnan(), sc_isinf(), sc_isfinite()
 NAN                               Use NPY_NAN (*almost* equivalent)
 inline functions                  Make static functions and place in an include .h file
 mid-block variable declarations   Declare variables at the top of block
================================  ========================================================

For C++, we currently suggest that anything beyond C++11 is going to be impossible
for a very long time because of ecosystem support restrictions.

Fortran Compilers
-----------------

Generally, any well-maintained compiler is likely suitable. gfortran >= 4.8.0,
ifort, or flang can be used to build SciPy.

Cython Compiler
---------------

SciPy always requires a recent Cython compiler. As of v1.2.0, the minumum version is 0.29.0

NumPy
-----
SciPy depends on NumPy but releases of SciPy are not tied to releases of NumPy.
SciPy attempts to be compatible with at least the 4 previous releases of NumPy.
In particular, SciPy can not rely on features of just the latest NumPy, but needs to be
written using what is common in all of those 4 releases.

========  ========================    ===========================
 Python     minimum Numpy version     maximum Numpy version
========  ========================    ===========================
3.7         1.13.1                     >= 1.16.x
3.6         1.12.1                     >= 1.16.x
3.5         1.9.3                      >= 1.16.x
2.7         1.8.2                      >= 1.16.x
========  ========================    ===========================



Other Libraries
---------------

- LAPACK: >= 3.4.1
- OpenBLAS: A recent version


Building the Documentation
--------------------------

- Sphinx: whatever recent versions work. >= 1.6.6.
- numpydoc: whatever recent versions work. >=  0.8.0.
- matplotlib: generally suggest >= 2.0

Testing and Benchmarking
------------------------
Recent version of:

- pytest https://docs.pytest.org/en/latest/
- asv (airspeed velocity)  https://asv.readthedocs.io/

Packaging
---------
Recent version of:

- wheel
- multibuild


References
----------

.. [1] https://blogs.msdn.microsoft.com/vcblog/2013/07/19/c99-library-support-in-visual-studio-2013/
.. [2] https://en.wikipedia.org/wiki/IEEE_754-1985
