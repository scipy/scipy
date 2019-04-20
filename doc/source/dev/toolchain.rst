Toolchain and Roadmap
=====================

The use of the SciPy library requires (or optionally depends upon) several
other libraries in order to operate, the main dependency being Python and NumPy.
It requires a larger collection of libraries and tools in order to build the library,
or to build the documentation.

Of course, the tooling and libraries are themselves not static.  This document aims to provide a
guide as to how SciPy's use of these dynamic dependencies will proceed over time.

SciPy aims to be compatible with a number of releases of its dependent libraries and tools.
Forcing the user base to other components for upgrade for every release would greatly
diminish the value of SciPy.  However, maintaining backwards compatibility
with very old tooling/libraries imposes limitations on which newer functionalities and capabilities
can be incorporated.  SciPy takes a somewhat conservative approach, maintaining compatibility with
several major releases of Python and NumPy on the major platforms. (That may in of itself impose
further restrictions.  See the C Compilers section for an example.)


- First and foremost, SciPy is a Python project hence it requires a Python environment.
- The LAPACK and OpenBLAS numerical libraries need to be installed.
- Compilers for C, C++, Cython and Fortran code are needed.
- The Python environment needs the ``NumPy`` package to be installed.
- Testing requires the ``pytest`` Python package.
- Building the documentation requires the ``matplotlib``, Sphinx and ``numpydoc`` packages, as well as a LaTeX installation.

The tooling used to build CPython has some implications for the tooling used in building SciPy.
It also has implications for the documentation in docstrings.


Building SciPy
--------------

Python Versions
^^^^^^^^^^^^^^^

SciPy is compatible with several versions of Python, and some
specific decisions are still under consideration, especially
with regard to future changes.
Python 2.7 support was dropped for SciPy
releases numbered 1.3 and above but is still available in Release 1.2.x,
which is a long-term support release. [1]_, [2]_.

================  =======================================================================
 Date             Pythons supported
================  =======================================================================
 2018              Py2.7, Py3.4+ (SciPy 1.2.x is the last release to support Python 2.7)
 2019              Py3.5+ (but Py2.7-specific code not removed)
 2020              Py3.5+ (removal of Py2.7-specific code permitted)
================  =======================================================================



C Compilers
^^^^^^^^^^^

SciPy is compatible with most modern C compilers.  However CPython on Windows is
built with specific versions of the Microsoft Visual C++ compiler [3]_, as is the
corresponding build of SciPy.  This has implications for the C language standards
that can be supported.

===================   ==============   ===================
CPython               MS Visual C++    C Standard
===================   ==============   ===================
2.7, 3.0, 3.1, 3.2       9.0           C90
3.3, 3.4                10.0           C90 & some of C99
3.5, 3.6                14.0           C90 & most of C99
===================   ==============   ===================



C and C++ Language Standards
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

C and C++ language standards for SciPy are generally guidelines
rather than official decisions. This is particularly true of
attempting to predict adoption timelines for newer standards.

================  ===========================================
 Date              C Standard
================  ===========================================
 <= 2018           C90
 2019              C90 for old code, may consider C99 for new
 2020              C99
 ?                 C11
 ?                 C17, C18
================  ===========================================

C99 has been a standard for 20 years, but experience has shown that
not all features are supported equally well across all platforms.
C18 is a bug fix for C11, so C11 may be skipped entirely.

For C++, we currently suggest that anything beyond C++11 is going to be impossible
for a very long time because of ecosystem support restrictions. See [4]_.

.. note::

    Developer Note: Some C99 features would be useful for scientific programming, in particular better support of
    IEEE 754 [5]_.  SciPy has a small include file ``scipy/_lib/_c99compat.h`` which provides
    access to a few functions.  Use in conjunction with ``<numpy/npy_math.h>``.

    ===================================   ========================================================
     Feature                               Workaround
    ===================================   ========================================================
     `isnan()`, `isinf()`, `isfinite()`   Use `sc_isnan()`, `sc_isinf()`, `sc_isfinite()`
     `NAN`                                Use `NPY_NAN` (*almost* equivalent)
     inline functions                     Make static functions and place in an include .h file
     mid-block variable declarations      Declare variables at the top of the block
    ===================================   ========================================================


Fortran Compilers
^^^^^^^^^^^^^^^^^

Generally, any well-maintained compiler is likely suitable and can be used to build SciPy.

- gfortran >= 4.8.0
- ifort
- flang

Cython Compiler
^^^^^^^^^^^^^^^

SciPy always requires a recent Cython compiler. As of v1.2.0, the minumum version is 0.29.0.

NumPy
^^^^^

SciPy depends on NumPy but releases of SciPy are not tied to releases of NumPy.
SciPy attempts to be compatible with at least the 4 previous releases of NumPy.
In particular, SciPy can not rely on features of just the latest NumPy, but needs to be
written using what is common in all of those 4 releases. [1]_, [6]_.

========  ========================    ===========================
 Python    Minimum NumPy version       Maximum NumPy version
========  ========================    ===========================
3.7         1.13.1                     >= 1.16.x
3.6         1.12.1                     >= 1.16.x
3.5         1.9.3                      >= 1.16.x
2.7         1.8.2                      1.16.x
========  ========================    ===========================


Other Libraries
^^^^^^^^^^^^^^^

- LAPACK: >= 3.4.1
- OpenBLAS: A recent version


Testing and Benchmarking
--------------------------

A Recent version of:

- pytest https://docs.pytest.org/en/latest/
- asv (airspeed velocity)  https://asv.readthedocs.io/
- mpmath http://mpmath.org


Building the Documentation
--------------------------

- Sphinx: whatever recent versions work. >= 2.0.
- numpydoc: whatever recent versions work. >=  0.8.0.
- matplotlib: generally suggest >= 2.0
- LaTeX: A recent distibution.


.. note::

    Developer Note: The version of ``matplotlib`` required has
    implications for the examples in Python docstrings.
    Examples must be able to be executed both in the environment used to build the documentation,
    as well as any supported version of ``matplotlib`` that a user may use with this release of SciPy.


Packaging
---------

A Recent version of:

- setuptools
- wheel  https://pythonwheels.com
- multibuild  https://github.com/matthew-brett/multibuild

:ref:`making-a-release` and :ref:`distributing-a-release` contain information on
making and distributing a SciPy release.

References
----------

.. [1] https://docs.scipy.org/doc/scipy/reference/release.1.2.0.html
.. [2] https://python3statement.org
.. [3] https://blogs.msdn.microsoft.com/vcblog/2013/07/19/c99-library-support-in-visual-studio-2013/
.. [4] https://en.cppreference.com/w/cpp/compiler_support
.. [5] https://en.wikipedia.org/wiki/IEEE_754-1985
.. [6] https://docs.scipy.org/doc/numpy/release.html
