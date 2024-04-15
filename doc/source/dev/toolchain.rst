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
libraries and tools. Forcing the user base to upgrade other components
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
- The Python environment needs the ``numpy`` package to be installed.
- Testing requires the ``pytest`` and ``hypothesis`` Python packages.
- Building the documentation requires the ``matplotlib``, Sphinx and MyST-NB_ packages along with PyData theme.

.. _MyST-NB: https://myst-nb.readthedocs.io/

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
older Python versions, SciPy takes guidance from [NEP29]_. Generally, support for
the oldest Python version is dropped 42 months after the original release. Following
the acceptance of PEP 602, this mostly happens in April, and gets picked up by the
mid-year release of SciPy.

.. dropdown:: Python version support over time

    Python 2.7 support was dropped starting from SciPy 1.3.

    ================  =======================================================================
     Date             Pythons supported
    ================  =======================================================================
     2024              Py3.10+
     2023              Py3.9+
     2022              Py3.8+
     2021              Py3.7+
     2020              Py3.6+
     2019              Py3.5+
     2018              Py2.7, Py3.4+
    ================  =======================================================================

NumPy
^^^^^

SciPy depends on NumPy but releases of SciPy are not tied to releases of NumPy.
SciPy attempts to be compatible with at least the 4 previous releases of NumPy.
In particular, SciPy cannot rely on features of just the latest NumPy, but
needs to be written using what is common in all of those 4 `NumPy releases`_.

.. _NumPy releases: https://numpy.org/doc/stable/release.html

.. dropdown:: Python and NumPy version support per SciPy version

    The table shows the NumPy versions suitable for each major Python version.
    This table does not distinguish SciPy patch versions (e.g. when a new Python
    version is released, SciPy will generally issue a compatible patch version).

    =================  ========================    =======================
     SciPy version      Python versions             NumPy versions
    =================  ========================    =======================
     1.13               >=3.9, <3.13                >=1.22.4, <2.3.0
     1.12               >=3.9, <3.13                >=1.22.4, <2.0.0
     1.11               >=3.9, <3.13                >=1.21.6, <1.27.0
     1.10               >=3.8, <3.12                >=1.19.5, <1.26.0
     1.9                >=3.8, <3.12                >=1.18.5, <1.26.0
     1.8                >=3.8, <3.11                >=1.17.3, <1.24.0
     1.7                >=3.7, <3.11                >=1.16.5, <1.23.0
     1.6                >=3.7, <3.10                >=1.16.5, <1.21.0
     1.5                >=3.6, <3.10                >=1.14.5, <1.20.0
     1.4                >=3.5, <3.9                 >=1.13.3, <1.18.0
     1.2                2.7, >=3.4, <3.8            >=1.8.2, <1.17.0
    =================  ========================    =======================

In specific cases, such as a particular architecture, these requirements
could vary. Please check the `release notes`_ and the meta-package
``oldest-supported-numpy`` for more info [OSN]_.

.. _release notes: https://scipy.github.io/devdocs/release.html

Compilers
^^^^^^^^^

Building SciPy requires compilers for C, C++, Fortran, as well as the
python transpilers Cython and Pythran (the latter is an opt-out dependency
starting from version 1.7.0).

To maintain compatibility with a large number of platforms & setups, especially
where using the official wheels (or other distribution channels like Anaconda
or conda-forge) is not possible, SciPy tries to keep compatibility with older
compilers, on platforms that have not yet reached their official end-of-life.

As explained in more detail below, the current minimal compiler versions are:

==========  ===========================  ===============================  ============================
 Compiler    Default Platform (tested)    Secondary Platform (untested)    Minimal Version
==========  ===========================  ===============================  ============================
 GCC         Linux                        AIX, Alpine Linux, OSX           GCC 9.x
 LLVM        OSX                          Linux, FreeBSD, Windows          LLVM 14.x
 MSVC        Windows                      -                                Visual Studio 2019 (vc142)
==========  ===========================  ===============================  ============================

Note that the lower bound for LLVM is not enforced. Older versions should
work - as long as they support core (non-stdlib) C++17 -, but no version
below LLVM 14 is tested regularly during development. Please file an issue
if you encounter a problem during compilation.

Official Builds
~~~~~~~~~~~~~~~

Currently, SciPy wheels are being built as follows:

================    ==============================   ==============================   =============================
 Platform            `CI`_ `Base`_ `Images`_          Compilers                        Comment
================    ==============================   ==============================   =============================
Linux x86            ``ubuntu-22.04``                 GCC 10.2.1                       ``cibuildwheel``
Linux arm            ``docker-builder-arm64``         GCC 11.3.0                       ``cibuildwheel``
OSX x86_64           ``macos-latest``                 clang-14/gfortran 13.0           ``cibuildwheel``
OSX arm64            ``macos-14``                     clang-14/gfortran 13.0           ``cibuildwheel``
Windows              ``windows-2019``                 GCC 10.3 (`rtools`_)             ``cibuildwheel``
================    ==============================   ==============================   =============================

.. _CI: https://github.com/actions/runner-images
.. _Base: https://cirrus-ci.org/guide/docker-builder-vm/#under-the-hood
.. _Images: https://github.com/orgs/cirruslabs/packages?tab=packages&q=macos
.. _rtools: https://community.chocolatey.org/packages/rtools#versionhistory

Note that the OSX wheels additionally vendor gfortran 11.3.0 for x86_64,
and gfortran 12.1.0 for arm64. See ``tools/wheels/cibw_before_build_macos.sh``.


C Compilers
~~~~~~~~~~~

SciPy is compatible with most modern C compilers (in particular ``clang``).
Nowadays, there is reasonable support for recent C language standards across
all relevant compilers, though this is very different from how things used to
be. The following paragraphs primarily discuss the evolution of these
constraints; readers who do not care about historical context can skip ahead
to the table at the end.

.. dropdown:: Historical context around ABI vs. compiler support vs. C standards

    In the past, the most restrictive compiler on relevant platforms in terms
    of C support was the Microsoft Visual C++ compiler & toolset (together known
    as MSVC; it has a complicated `version scheme`_) [MSVC]_.
    Up until Visual Studio 2013, each MSVC version came with
    an updated C Runtime (CRT) library that was incompatible with the previous
    ones.

    This lack of compatibility of the Application Binary Interface (ABI) meant
    that all projects wanting to communicate across this interface (e.g. calling a
    function from a shared library) needed to be (re)compiled with the same MSVC
    version. The long support of CPython 2.7 meant that python itself was stuck
    for a long time with VS 2008 (in order not to break the ABI in patch
    releases), and thus SciPy was stuck on that version as well.

    The use of VS 2008 (which doesn't have support for C99) to compile builds for
    CPython 2.7 meant for a long time that C code in SciPy has had to conform
    to the earlier C90 standard for the language and standard library. After
    dropping support for CPython 2.7 in SciPy 1.3.x, that restriction was finally
    lifted (though only gradually at first).

    With the introduction of the "Universal C Runtime" [UCRT]_ since the
    release of Visual Studio 2015, the ABI of C Runtime has been stable, which
    means that the restriction of having to use the same compiler version for
    SciPy as for the underlying CPython version is no longer applicable. This
    stability is not indefinite though: Microsoft has been planning an
    ABI-breaking release - across the compiler resp. C/C++ standard libraries -
    (tentatively called "`vNext`_") for quite a while, but so far it is unclear
    when this will arrive. Once that happens, SciPy will again be restricted to
    at most the last ABI-compatible Visual Studio release (currently VS 2022)
    until all CPython versions supported according to NEP29 have been built
    upstream with vNext-compatible compilers.

    More specifically, there is a distinction between the Microsoft Visual
    Studio version and the version of the targeted "`toolset`_", which is defined
    as "The Microsoft C++ compiler, linker, standard libraries, and related
    utilities". Each version of Visual Studio comes with a default version of the
    MSVC toolset (for example VS2017 with vc141, VS2019 with vc142), but it is
    possible to target older toolsets even in newer versions of Visual Studio.
    Due to the nature of compilers (i.e. split into frontend and backend), it
    depends whether the limiting factor for supporting a given feature (e.g. in C)
    is due to the version of Visual Studio or the toolset, but in general the
    latter is a harder barrier and thus the effective lower bound.

    This is due to the fact that while the ABI stays compatible between toolset
    versions (until vNext), all linking operations must use a toolset at least
    as new as the one used to build any of the involved artefacts, meaning that
    toolset version bumps tend to be "infectious", as in: requiring all consuming
    libraries to also bump their toolset (and probably compiler) version. This is
    more of an issue for NumPy than SciPy, as the latter has only a small C API
    and is compiled against by far fewer projects than NumPy. Additionally, using
    a newer toolset means that users of libraries that compile C++ code (as SciPy
    does) might also need a newer Microsoft Visual C++ `Redistributable`_, which
    might have to be distributed to them.

    Summing up, the minimal requirement for the MSVC compiler resp. toolset per
    SciPy version was determined predominantly by the oldest supported CPython
    version at the time. The first SciPy version to raise the minimal requirement
    beyond that was SciPy 1.9, due to the inclusion of the HiGHS submodule, which
    does not compile with vc141 (and the aggressive removal of VS2017 in public CI
    making it infeasible to keep ensuring that everything everywhere works with
    non-default toolset versions).

.. _version scheme: https://en.wikipedia.org/wiki/Microsoft_Visual_C%2B%2B#Internal_version_numbering
.. _vNext: https://github.com/microsoft/STL/issues/169
.. _toolset: https://docs.microsoft.com/en-us/cpp/build/projects-and-build-systems-cpp#the-msvc-toolset
.. _Redistributable: https://docs.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist

==============  =================  =================  =================
SciPy version    CPython support    MS Visual C++      Toolset version
==============  =================  =================  =================
 Until 1.2       2.7 & 3.4+         VS 2008 (9.0)      vc90
 1.3, 1.4        3.5+               VS 2010 (10.0)     vc100
 1.5             3.6+               VS 2015 (14.0)     vc140
 1.6, 1.7        3.7+               VS 2017 (14.1)     vc141
 1.8             3.8+               VS 2017 (14.1)     vc141
 1.9             3.8+               VS 2019 (14.20)    vc142
==============  =================  =================  =================

In terms of C language standards, it's relevant to note that C11 has `optional features`_
(e.g. atomics, threading), some of which (VLAs & complex types)
were mandatory in the C99 standard. C17 (occasionally called C18) can be
considered a bug fix for C11, so generally, C11 may be skipped entirely.

.. _optional features: https://en.wikipedia.org/wiki/C11_%28C_standard_revision%29#Optional_features

SciPy has been restricted in the use of more advanced language features by the
available compiler support, and Microsoft in particular has taken very long to
achieve conformance to C99/C11/C17, however starting from `Visual Studio 16.8`_,
C11/C17 is supported (though without the C11 optional features).
C99 ``<complex.h>`` `support <https://developercommunity.visualstudio.com/t/714008>`_
would be particularly interesting for SciPy.
However, it's still possible to use complex types on windows, provided that
`windows-specific types`_ are used.

.. _Visual Studio 16.8: https://docs.microsoft.com/en-us/cpp/overview/visual-cpp-language-conformance#c-standard-library-features-1
.. _windows-specific types: https://docs.microsoft.com/en-us/cpp/c-runtime-library/complex-math-support

Therefore, using C features beyond C90 was only possible insofar as there was support on
Windows; however, as of as of the end of 2021, a sufficiently recent compiler is used.
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
 2022              C++17 (core language + universally available stdlib features)
 ?                 C++17 (with full stdlib), C++20, C++23, C++26
================  =======================================================================

.. dropdown:: Historical context for compiler constraints due to manylinux

    Since dropping support for Python 2.7, C++11 can be used
    universally, and since dropping Python 3.6, the Visual Studio version
    (that had previously been stuck with 14.0 due to ABI compatibility with
    CPython) has been recent enough to support even C++17.

    Since the official builds (see above) use a pretty recent version of LLVM,
    the bottleneck for C++ support is therefore the oldest supported GCC version,
    where SciPy has been constrained mainly by the version in the oldest supported
    manylinux versions & images [MANY]_.

    At the end of 2021 (with the final removal of ``manylinux1`` wheels), the
    minimal requirement of GCC moved to 6.3, which has full C++14 support [CPP]_.
    This corresponded to the lowest-present GCC version in relevant manylinux
    versions, though this was still considering the Debian-based "outlier"
    ``manylinux_2_24``, which - in contrast to previous manylinux images based on
    RHEL-derivative CentOS that could benefit from the ABI-compatible GCC backports
    in the "RHEL Dev Toolset" - was stuck with GCC 6.3. That image failed to take
    off not least due to those `outdated compilers`_ and reached its EOL in
    mid-2022. For different reasons, ``manylinux2010`` also reached its EOL
    around the `same time <https://github.com/pypa/manylinux/issues/1281>`_.

    The remaining images ``manylinux2014`` and ``manylinux_2_28`` currently support
    GCC 10 and 12, respectively. The latter will continue to receive updates as new
    GCC versions become available as backports, but the former will likely not
    change since the CentOS project is not responsive anymore about publishing
    aarch64 `backports <https://github.com/pypa/manylinux/issues/1266>`_ of GCC 11.

This leaves all the main platforms and their compilers with comparatively
recent versions. However, SciPy has historically also endeavored to support
less common platforms as well - if not with binary artefacts (i.e. wheels),
then at least by remaining compilable from source - which includes for example
AIX, Alpine Linux and FreeBSD.

.. dropdown:: Platform support and other constraints on compiler

    For AIX 7.2 & 7.3 the default compiler is GCC 10 (AIX 7.1 had its EOL in 2023),
    but GCC 11/12 is installable `side-by-side`_, and similarly, there is the
    LLVM 17-based `Open XL`_ for AIX.

    The oldest currently-supported `Alpine Linux`_ release is 3.16, and already
    `comes with <https://distrowatch.com/table.php?distribution=alpine>`_ GCC 11.
    For `FreeBSD`_, the oldest currently-supported `13.x release`_ comes with
    LLVM 14 (and GCC 13 is available as a `freebsd-port`_).

    Finally there is the question of which machines are widely used by people
    needing to compile SciPy from source for other reasons (e.g. SciPy developers,
    or people wanting to compile for themselves for performance reasons).
    The oldest relevant distributions (without RHEL-style backports) are Ubuntu
    20.04 LTS (which has GCC 9 but also has a backport of GCC 10; Ubuntu 22.04 LTS
    has GCC 11) and Debian Bullseye (with GCC 10; Bookworm has GCC 12).
    This is the weakest restriction for determining the lower bounds of compiler
    support (power users and developers can be expected to keep their systems at
    least somewhat up-to-date, or use backports where available), and gradually
    becomes less important as usage numbers of old distributions dwindle.

.. _outdated compilers: https://github.com/pypa/manylinux/issues/1012
.. _side-by-side: https://www.ibm.com/support/pages/aix-toolbox-open-source-software-downloads-alpha#G
.. _Open XL: https://www.ibm.com/docs/en/openxl-c-and-cpp-aix/17.1.2?topic=new-enhanced-llvm-clang-support
.. _Alpine Linux: https://alpinelinux.org/releases/
.. _FreeBSD: https://www.freebsd.org/releases/
.. _13.x release: https://www.freebsd.org/releases/13.2R/relnotes/
.. _freebsd-port: https://ports.freebsd.org/cgi/ports.cgi?query=gcc

All the currently lowest-supported compiler versions (GCC 9, LLVM 14,
VS2019 with vc142) have full support for the C++17 *core language*,
which can therefore be used unconditionally.
However, as of mid-2024, support for the entirety of the C++17 standard library
has not yet been completed across all compilers [CPP]_, particularly LLVM.
It is therefore necessary to check if a given stdlib-feature is supported by
all compilers before it can be used in SciPy.

C++20 support is stabilizing very slowly, even aside from modules, coroutines
and several not-yet-universally-supported stdlib features. Given how big of a
release the C++20 standard was, it is expected that it will take a `while yet`_
before we can start considering moving our baseline.
Compiler support for C++23 and C++26 is still under heavy development [CPP]_.

.. _while yet: https://discourse.llvm.org/t/rfc-clang-17-0-6-would-be-minimum-version-to-build-llvm-in-c-20/75345/8

Fortran Compilers
~~~~~~~~~~~~~~~~~

Generally, any well-maintained compiler is likely suitable and can be
used to build SciPy. That said, we do not test with old ``gfortran`` versions,
which is why we are matching the lower bound with the one for GCC above.

============= =====================================
 Tool          Version
============= =====================================
gfortran       >= 9.x
ifort/ifx      A recent version (not tested in CI)
flang (LLVM)   >= 17.x
============= =====================================


Cython & Pythran
~~~~~~~~~~~~~~~~

SciPy always requires a recent Cython compiler. Since 1.7, Pythran
is a build dependency (currently with the possibility to opt out).


OpenMP support
^^^^^^^^^^^^^^

For `various reasons`_, SciPy cannot be distributed with built-in OpenMP support.
When using the optional Pythran support, OpenMP-enabled parallel code can be
generated when building from source.

.. _various reasons: https://github.com/scipy/scipy/issues/10239

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
pooch           Recent    https://pypi.org/project/pooch/
=============== ======== ==========================================


Moreover, SciPy supports interaction with other libraries. The test suite
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
Hypothesis                 Recent     https://hypothesis.readthedocs.io/
asv (airspeed velocity)    Recent     https://asv.readthedocs.io/
=========================  ========  ====================================


Building the Documentation
--------------------------

====================  =================================================
 Tool                 Version
====================  =================================================
Sphinx                Whatever recent versions work. >= 5.0.
PyData Sphinx theme   Whatever recent versions work. >= 0.15.2.
Sphinx-Design         Whatever recent versions work. >= 0.4.0.
numpydoc              Whatever recent versions work. >= 1.5.0.
matplotlib            Generally suggest >= 3.5.
MyST-NB               Whatever recent versions work. >= 0.17.1
jupyterlite-sphinx    Whatever recent versions work. >= 0.12.0
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

.. [CPP] https://en.cppreference.com/w/cpp/compiler_support
.. [MANY] https://github.com/mayeut/pep600_compliance
.. [MSVC] https://docs.microsoft.com/en-us/cpp/overview/visual-cpp-in-visual-studio
.. [NEP29] https://numpy.org/neps/nep-0029-deprecation_policy.html
.. [OSN] https://github.com/scipy/oldest-supported-numpy
.. [UCRT] https://docs.microsoft.com/en-gb/cpp/windows/universal-crt-deployment
