Distributing
============

Distributing Python packages is nontrivial - especially for a package with
complex build requirements like Scipy - and subject to change.  For an up-to-date
overview of recommended tools and techniques, see the `Python Packaging User
Guide`_.  This document discusses some of the main issues and considerations for
Scipy.

Dependencies
------------
Dependencies are things that a user has to install in order to use (or
build/test) a package.  They usually cause trouble, especially if they're not
optional.  Scipy tries to keep its dependencies to a minimum; currently they
are:

*Unconditional run-time dependencies:*

- Numpy_

*Conditional run-time dependencies:*

- nose_ (to run the test suite)
- asv_ (to run the benchmarks)
- matplotlib_ (for some functions that can produce plots)
- Pillow_ (for image loading/saving)
- scikits.umfpack_ (optionally used in ``sparse.linalg``)
- mpmath_ (for more extended tests in ``special``)

*Unconditional build-time dependencies:*

- Numpy_
- A BLAS and LAPACK implementation (reference BLAS/LAPACK, ATLAS, OpenBLAS,
  MKL, Accelerate are all known to work)
- (for development versions) Cython_

*Conditional build-time dependencies:*

- setuptools_
- wheel_ (``python setup.py bdist_wheel``)
- Sphinx_ (docs)
- matplotlib_ (docs)
- LaTeX (pdf docs)
- Pillow_ (docs)

Furthermore of course one needs C, C++ and Fortran compilers to build Scipy,
but those we don't consider to be dependencies and are therefore not discussed
here.  For details, see http://scipy.org/scipylib/building/index.html.

When a package provides useful functionality and it's proposed as a new
dependency, consider also if it makes sense to vendor (i.e. ship a copy of it with
scipy) the package instead.  For example, six_ and decorator_ are vendored in
``scipy._lib``.

The only dependency that is reported to pip_  is Numpy_, see
``install_requires`` in Scipy's main ``setup.py``.  The other dependencies
aren't needed for Scipy to function correctly, and the one unconditional build
dependency that pip_ knows how to install (Cython_) we prefer to treat like a
compiler rather than a Python package that pip_ is allowed to upgrade.

Issues with dependency handling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
There are some serious issues with how Python packaging tools handle
dependencies reported by projects.  Because Scipy gets regular bug reports
about this, we go in a bit of detail here.

Scipy only reports its dependency on Numpy via ``install_requires`` if Numpy
isn't installed at all on a system.  This will only change when there are
either 32-bit and 64-bit Windows wheels for Numpy on PyPI or when
``pip upgrade`` becomes available (with sane behavior, unlike ``pip install
-U``, see `this PR
<https://github.com/pypa/pip/pull/3194>`_).  For more details, see
`this summary <http://article.gmane.org/gmane.comp.python.distutils.devel/24218>`_.

The situation with ``setup_requires`` is even worse; pip_ doesn't handle that
keyword at all, while ``setuptools`` has issues (here's a `current one
<https://bitbucket.org/pypa/setuptools/issues/391>`_) and invokes
``easy_install`` which comes with its own set of problems (note that Scipy doesn't
support ``easy_install`` at all anymore; issues specific to it will be closed
as "wontfix").


Supported Python and Numpy versions
-----------------------------------
The Python_ versions that Scipy supports are listed in the list of PyPI
classifiers in ``setup.py``, and mentioned in the release notes for each
release.  All newly released Python versions will be supported as soon as
possible.  The general policy on dropping support for a Python version is that
(a) usage of that version has to be quite low (say <5% of users) and (b) the
version isn't included in an active long-term support release of one of the
main Linux distributions anymore.  Scipy typically follows Numpy, which has a
similar policy.  The final decision on dropping support is always taken on the
scipy-dev mailing list.

The lowest supported Numpy_ version for a Scipy version is mentioned in the
release notes and is encoded in ``scipy/__init__.py`` and the
``install_requires`` field of ``setup.py``.  Typically the latest Scipy release
supports 3 or 4 minor versions of Numpy.  That may become more if the frequency
of Numpy releases increases (it's about 1x/year at the time of writing).
Support for a particular Numpy version is typically dropped if (a) that Numpy
version is several years old, and (b) the maintenance cost of keeping support
is starting to outweigh the benefits.  The final decision on dropping support
is always taken on the scipy-dev mailing list.

Supported versions of optional dependencies and compilers is less clearly
documented, and also isn't tested well or at all by Scipy's Continuous
Integration setup.  Issues regarding this are dealt with as they come up in the
issue tracker or mailing list.


Building binary installers
--------------------------
.. note::

   This section is only about building Scipy binary installers to *distribute*.
   For info on building Scipy on the same machine as where it will be used, see
   `here <http://scipy.org/scipylib/building/index.html>`_.

There are a number of things to take into consideration when building binaries
and distributing them on `PyPI`_ or elsewhere.

**General**

- A binary is specific for a single Python version (because different Python
  versions aren't ABI-compatible, at least up to Python 3.4).
- Build against the lowest Numpy version that you need to support, then it will
  work for all Numpy versions with the same major version number (Numpy does
  maintain backwards ABI compatibility).

**Windows**

- For 64-bit Windows installers built with a free toolchain, use the method
  documented at https://github.com/numpy/numpy/wiki/Mingw-static-toolchain.
  That method will likely be used for Scipy itself once it's clear that the
  maintenance of that toolchain is sustainable long-term.  See the MingwPy_
  project and `this thread
  <http://article.gmane.org/gmane.comp.python.numeric.general/61727>`_ for
  details.
- The other way to produce 64-bit Windows installers is with ``icc``, ``ifort``
  plus ``MKL`` (or ``MSVC`` instead of ``icc``).  For Intel toolchain
  instructions see
  `here <https://software.intel.com/en-us/articles/numpyscipy-with-intel-mkl>`_
  and for (partial) MSVC instructions see
  `here <https://github.com/numpy/numpy/wiki/Building-with-MSVC>`_.
- Older Scipy releases contained a .exe "superpack" installer.  Those contain
  3 complete builds (no SSE, SSE2, SSE3), and were built with
  https://github.com/numpy/numpy-vendor.  That build setup is known to not work
  well anymore and is no longer supported.  It used g77 instead of gfortran,
  due to complex DLL distribution issues (see `gh-2829
  <https://github.com/scipy/scipy/issues/2829>`_).  Because the toolchain is no
  longer supported, g77 support isn't needed anymore and Scipy can now include
  Fortran 90/95 code.

**OS X**

- To produce OS X wheels that work with various Python versions (from
  python.org, Homebrew, MacPython), use the build method provided by
  https://github.com/MacPython/scipy-wheels.
- DMG installers for the Python from python.org on OS X can still be produced
  by ``tools/scipy-macosx-installer/``.  Scipy doesn't distribute those
  installers anymore though, now that there are binary wheels on PyPi.

**Linux**

Besides PyPi not allowing Linux wheels (which is about to change with `PEP 513
<https://www.python.org/dev/peps/pep-0513>`_), there are no specific issues with
building binaries.  To build a set of wheels for a Linux distribution and
providing them in a Wheelhouse_, look at the wheel_ and Wheelhouse_ docs.  A
Wheelhouse for wheels compatible with TravisCI is http://wheels.scipy.org.



.. _Numpy: http://numpy.org
.. _Python: http://python.org
.. _nose: http://nose.readthedocs.org
.. _asv: http://asv.readthedocs.org
.. _matplotlib: http://matplotlib.org
.. _Pillow: http://pillow.readthedocs.org
.. _scikits.umfpack: https://pypi.python.org/pypi/scikit-umfpack
.. _mpmath: http://mpmath.org
.. _Cython: http://cython.org
.. _setuptools: https://bitbucket.org/pypa/setuptools
.. _wheel: wheel.readthedocs.org
.. _pip: http://pip-installer.org
.. _PyPI: http://pypi.python.org/pypi/scipy
.. _Python Packaging User Guide: https://packaging.python.org
.. _Wheelhouse: https://pypi.python.org/pypi/Wheelhouse
.. _MingwPy: https://mingwpy.github.io
