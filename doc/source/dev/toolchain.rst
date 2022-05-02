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
- Building the documentation requires the ``matplotlib``, Sphinx packages along with PyData theme.

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
older Python versions, SciPy takes guidance from NEP 29 [1]_.  Python 2.7
support was dropped starting from SciPy 1.3.

================  =======================================================================
 Date             Pythons supported
================  =======================================================================
 2018              Py2.7, Py3.4+ (SciPy 1.2.x is the last release to support Python 2.7)
 2019              Py3.5+ (but Py2.7-specific code not removed)
 2020              Py3.6+ (removal of Py2.7-specific code permitted)
 2021              Py3.7+
 2022              Py3.8+
================  =======================================================================

NumPy
^^^^^

SciPy depends on NumPy but releases of SciPy are not tied to releases of NumPy.
SciPy attempts to be compatible with at least the 4 previous releases of NumPy.
In particular, SciPy cannot rely on features of just the latest NumPy, but
needs to be written using what is common in all of those 4 releases [2]_.

The table shows the NumPy versions suitable for each major Python version.

=================  ========================    =======================
 SciPy version      Python versions             NumPy versions
=================  ========================    =======================
 1.2                2.7, >=3.4, <=3.7           >=1.8.2, <= 1.16.x
 1.4                >=3.5, <=3.8                >=1.13.3, <= 1.17.3
 1.5                >=3.6, <=3.9                >=1.14.5, <= 1.19.3
 1.6                >=3.7, <=3.9                >=1.16.5, <= 1.20.x
 1.7.0/1            >=3.7, <3.10                >=1.16.5, <1.23.0
 1.7.2-x            >=3.7, <3.11                >=1.16.5, <1.24.0
 1.8                >=3.8, <3.11                >=1.17.3, <1.24.0
=================  ========================    =======================

In specific cases, such as a particular architecture, these requirements
could vary. Please check the release notes [3]_ and the meta-package
``oldest-supported-numpy`` for more info [4]_.


Compilers
^^^^^^^^^

Building SciPy requires compilers for C, C++, Fortran, as well as the
python transpilers Cython and Pythran (the latter is an opt-out dependency
starting from version 1.7.0).

To maintain compatibility with a large number of platforms & setups, especially
where using the official wheels (or other distribution channels like Anaconda
or conda-forge) is not possible, SciPy keeps compatibility with old compilers.

Official Builds
~~~~~~~~~~~~~~~

Currently, SciPy wheels are being built as follows:

================  ========================  ===========================  ==============================
 Platform          Azure Base Image [5]_     Compilers                    Comment
================  ========================  ===========================  ==============================
Linux (nightly)    ``ubuntu-18.04``          GCC 6.5                      See ``azure-pipelines.yml``
Linux (release)    ``ubuntu-18.04``          GCC 7.5                      Built in separate repo [6]_
OSX                ``macOS-10.15``           LLVM 12.0.0                  Built in separate repo [6]_
Windows            ``windows-2019``          Visual Studio 2019 (16.11)   See ``azure-pipelines.yml``
================  ========================  ===========================  ==============================

Note that the OSX wheels additionally vendor gfortran 4.9,
see submodule ``gfortran-install`` in [6]_.


C Compilers
~~~~~~~~~~~

SciPy is compatible with most modern C compilers (in particular ``clang``).
In addition to concerns about compatibility with non-standard platforms,
there was a long-standing restriction that Windows builds of SciPy had to use
the same version of the Microsoft Visual C++ compiler as were used for CPython
itself, for reasons of ABI-compatibility [7]_, [8]_.

With the introduction of the "Universal C Runtime" [9]_ since the release of
Visual Studio 2015, this restriction has been lifted. For more context, see the
explanations by Steve Dower (member of the CPython-on-Windows core developers)
on this topic [10]_.

The use of MS Visual Studio 9.0 (which doesn't have support for C99)
to build Python 2.7 has meant that C code in SciPy has had to conform
to the earlier C90 standard for the language and standard library.
With the dropping of Python 2.7 for SciPy 1.3.x, the C90 restriction is no
longer imposed by compilers.

In terms of C language standards, it's relevant to note that C11 has optional
features [11]_ (e.g. atomics, threading), some of which (VLAs & complex types)
were mandatory in the C99 standard. C17 (occasionally called C18) can be
considered a bug fix for C11, so generally, C11 may be skipped entirely.

SciPy has been restricted in the use of more advanced language features by the
available compiler support, and Microsoft in particular has taken very long to
achieve conformance to C99/C11/C17, however starting from MS Visual Studio 16.8,
C11/C17 is supported [12]_ (though without the C11 optional features).
C99 ``<complex.h>`` would be particularly interesting for SciPy;
MSVC conformance for this is being tracked here [13]_.

Therefore, using C features beyond C90 was only possible insofar there was support on
windows; however, as of as of the end of 2021, a sufficiently recent compiler is used.
This is because GCC & LLVM support all relevant C11 features with the oldest currently
used versions, and C17 is just a bugfix for C11, as mentioned above. In short:

================  =======================================================================
 Date              C Standard
================  =======================================================================
 <= 2018           C90
 2019              C90 for old code, may consider C99 for new
 2020              C99 (no ``<complex.h>``, ``<stdatomic.h>``, ``<threads.h>`` & VLAs)
 2021              C17 (no ``<complex.h>``, ``<stdatomic.h>``, ``<threads.h>`` & VLAs)
 ?                 C23, ``<complex.h>``, ``<stdatomic.h>``, ...
================  =======================================================================


C++ Language Standards
~~~~~~~~~~~~~~~~~~~~~~

C++ language standards for SciPy are generally guidelines
rather than official decisions. This is particularly true of
attempting to predict adoption timelines for newer standards.

================  =======================================================================
 Date              C++ Standard
================  =======================================================================
 <= 2019           C++03
 2020              C++11
 2021              C++14
 ?                 C++17, C++20, C++23
================  =======================================================================

Since dropping support for Python 2.7, C++11 can be used
universally, and since dropping Python 3.6, the Visual Studio version
(that had previously been stuck with 14.0 due to ABI compatibility with
CPython) has been recent enough to support even C++17.

Since the official builds (see above) use a pretty recent version of LLVM,
the bottleneck for C++ support is therefore the oldest supported GCC version,
where SciPy has been constrained mainly by the version in the oldest supported
manylinux versions & images [14]_.

At the end of 2021 (with the final removal of ``manylinux1`` wheels), SciPy
now has a minimum GCC requirement of GCC 6.3, which has full C++14 support
[15]_. This corresponds to the lowest present GCC version in relevant manylinux
versions - somewhat surprisingly, it is not the oldest remaining
``manylinux2010`` that is the most restrictive (due to the ABI-compatible
"RHEL Dev Toolset" backports, it has GCC 8.3), but actually ``manylinux_2_24``
that only comes with GCC 6.3 [16]_.

C++17 _language_ support will require GCC >= 7 (released May 2017). As of the
end of 2021, support for the entirety of the C++17 standard library has not yet
been completed across all compilers; similarly, support for C++20 and C++23
is still under heavy development. [15]_

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


Cython & Pythran
~~~~~~~~~~~~~~~~

SciPy always requires a recent Cython compiler. Since 1.7, Pythran
is a build dependency (currently with the possibility to opt out).


OpenMP support
^^^^^^^^^^^^^^

For various reasons [17]_, SciPy cannot be distributed with built-in OpenMP support.
When using the optional Pythran support, OpenMP-enabled parallel code can be
generated when building from source.

Other Libraries
^^^^^^^^^^^^^^^

Any library conforming to the BLAS/LAPACK interface may be used.
OpenBLAS, ATLAS, MKL, BLIS, and reference Netlib libraries are known to work.

=============== =====================================================
 Library           Minimum version
=============== =====================================================
LAPACK           3.7.1
BLAS             A recent version of OpenBLAS, MKL or ATLAS.
                 The Accelerate BLAS library is no longer supported.
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
PyData Sphinx theme   Whatever recent versions work. >= 0.8.1.
Sphinx-Panels         Whatever recent versions work. >= 0.5.2.
Sphinx-Tabs           Whatever recent versions work. >= 3.2.0.
numpydoc              Whatever recent versions work. >= 0.8.0.
matplotlib            Generally suggest >= 2.0.
====================  =================================================

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

.. [1] https://numpy.org/neps/nep-0029-deprecation_policy.html
.. [2] https://numpy.org/doc/stable/release.html
.. [3] https://scipy.github.io/devdocs/release.html
.. [4] https://github.com/scipy/oldest-supported-numpy
.. [5] https://docs.microsoft.com/en-us/azure/devops/pipelines/agents/hosted
.. [6] https://github.com/MacPython/scipy-wheels
.. [7] https://pythondev.readthedocs.io/windows.html#python-and-visual-studio-version-matrix
.. [8] https://en.wikipedia.org/wiki/Microsoft_Visual_C%2B%2B#Internal_version_numbering
.. [9] https://docs.microsoft.com/en-gb/cpp/windows/universal-crt-deployment
.. [10] https://discuss.python.org/t/toolchain-upgrade-on-windows/6377/4
.. [11] https://en.wikipedia.org/wiki/C11_%28C_standard_revision%29#Optional_features
.. [12] https://devblogs.microsoft.com/cppblog/c11-and-c17-standard-support-arriving-in-msvc/
.. [13] https://developercommunity.visualstudio.com/t/Support-for-C99-Complex-numbers/1409049?space=8&q=complex
.. [14] https://github.com/mayeut/pep600_compliance
.. [15] https://en.cppreference.com/w/cpp/compiler_support
.. [16] https://github.com/pypa/manylinux/issues/1012
.. [17] https://github.com/scipy/scipy/issues/10239
