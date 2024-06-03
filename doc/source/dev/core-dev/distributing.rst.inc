.. _distributing-a-release:

Distributing
============

Distributing Python packages is nontrivial - especially for a package with
complex build requirements like SciPy - and subject to change.  For an up-to-date
overview of recommended tools and techniques, see the `Python Packaging User
Guide`_.  This document discusses some of the main issues and considerations for
SciPy.

Dependencies
------------
Dependencies are things that a user has to install in order to use (or
build/test) a package.  They usually cause trouble, especially if they're not
optional.  SciPy tries to keep its dependencies to a minimum; the current required
and optional build time dependencies can be seen in `SciPy's configuration file`_, 
``pyproject.toml``. The only non-optional runtime dependency is NumPy_.

Furthermore, of course one needs C, C++ and Fortran compilers to build SciPy,
but we don't consider those to be dependencies, and therefore they are not discussed
here.  For details, see :ref:`building-from-source`.

When a package provides useful functionality and it's proposed as a new
dependency, consider also if it makes sense to vendor (i.e. ship a copy of it with
SciPy) the package instead.  For example, decorator_ is vendored in
``scipy._lib``.

Issues with dependency handling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
There are some issues with how Python packaging tools handle
dependencies reported by projects.  Because SciPy gets regular bug reports
about this, we go in a bit of detail here.

SciPy only reports its dependency on NumPy via ``install_requires`` if NumPy
isn't installed at all on a system, or when building wheels with
``bdist_wheel``. SciPy no longer uses ``setup_requires`` (which in the past
invoked ``easy_install``); build dependencies are now handled only via
``pyproject.toml``. ``pyproject.toml`` relies on PEP 517; ``pip`` has
``--no-use-pep517`` and ``--no-build-isolation`` flags that may ignore
``pyproject.toml`` or treat it differently - if users use those flags, they
are responsible for installing the correct build dependencies themselves.


.. _numpy-version-ranges:

Version ranges for NumPy and other dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For dependencies it's important to set lower and upper bounds on their
versions. For *build-time* dependencies, they are specified in
``pyproject.toml`` and the versions will *only* apply to the SciPy build
itself. It's fine to specify either a range or a specific version for a
dependency like ``wheel`` or ``setuptools``. For NumPy we have to worry
about ABI compatibility too, hence we specify the version with ``==``
to the lowest supported version (because NumPy's ABI is backward but not
forward compatible).

For *run-time dependencies* (currently only ``numpy``), we specify the range
of versions in ``pyproject.toml`` and in ``install_requires`` in ``setup.py``.
Getting the upper bound right is slightly tricky.  If we don't set any bound, a
too-new version will be pulled in a few years down the line, and NumPy may have
deprecated and removed some API that SciPy depended on by then. On the other
hand if we set the upper bound to the newest already-released version, then as
soon as a new NumPy version is released there will be no matching SciPy version
that works with it. Given that NumPy and SciPy both release in a 6-monthly
cadence and that features that get deprecated in NumPy should stay around for
another two releases, we specify the upper bound as ``<1.xx+3.0`` (where ``xx``
is the minor version of the latest already-released NumPy.


.. _supported-py-numpy-versions:

Supported Python and NumPy versions
-----------------------------------
The Python_ versions that SciPy supports are listed in the list of PyPI
classifiers in ``setup.py``, and mentioned in the release notes for each
release.  All newly released Python versions will be supported as soon as
possible.  For the general policy on dropping support for a Python or NumPy
version, see :ref:`NEP 29 <NEP29>`.  The final decision on dropping support is
always taken on the scipy-dev forum.

The lowest supported Numpy_ version for a SciPy version is mentioned in the
release notes and is encoded in ``pyproject.toml``, ``scipy/__init__.py`` and the
``install_requires`` field of ``setup.py``.  Typically the latest SciPy release
supports ~5-7 minor versions of NumPy: up to 2.5 years' old NumPy versions,
(given that the frequency of NumPy releases is about 2x/year at the time of
writing) plus two versions into the future.

Supported versions of optional dependencies and compilers is documented in
:ref:`toolchain-roadmap`. Note that not all versions of optional dependencies
that are supported are tested well or at all by SciPy's Continuous
Integration setup.  Issues regarding this are dealt with as they come up in the
issue tracker or forum.


Building binary installers
--------------------------
.. note::

   This section is only about building SciPy binary installers to *distribute*.
   For info on building SciPy on the same machine as where it will be used, see
   `this scipy.org page <https://scipy.github.io/devdocs/dev/contributor/building.html>`_.

There are a number of things to take into consideration when building binaries
and distributing them on PyPI or elsewhere.

**General**

- A binary is specific for a single Python version (because different Python
  versions aren't ABI-compatible, at least up to Python 3.4).
- Build against the lowest NumPy version that you need to support, then it will
  work for all NumPy versions with the same major version number (NumPy does
  maintain backwards ABI compatibility).
- The easiest available toolchain for building portable SciPy binaries
  is our ``cibuildwheel`` infrastructure for common platforms, with
  details available in our CI infrastructure code and available via the
  ``cibuildwheel`` command on Windows, Linux, and MacOS, albeit with some extra
  external dependencies required in some cases

**Windows**

- For 64-bit Windows installers built with a free toolchain, use the method
  documented at https://github.com/numpy/numpy/wiki/Mingw-static-toolchain.
  That method will likely be used for SciPy itself once it's clear that the
  maintenance of that toolchain is sustainable long-term.  See the MingwPy_
  project and `this thread
  <https://mail.python.org/pipermail/numpy-discussion/2015-October/074056.html>`_ for
  details.
- The other way to produce 64-bit Windows installers is with ``icc``, ``ifort``
  plus ``MKL`` (or ``MSVC`` instead of ``icc``).  For Intel toolchain
  instructions see
  `this article <https://software.intel.com/en-us/articles/numpyscipy-with-intel-mkl>`_
  and for (partial) MSVC instructions see
  `this wiki page <https://github.com/numpy/numpy/wiki/Building-with-MSVC>`_.
- Older SciPy releases contained a .exe "superpack" installer.  Those contain
  3 complete builds (no SSE, SSE2, SSE3), and were built with
  https://github.com/numpy/numpy-vendor.  That build setup is known to not work
  well anymore and is no longer supported.  It used g77 instead of gfortran,
  due to complex DLL distribution issues (see `gh-2829
  <https://github.com/scipy/scipy/issues/2829>`_).  Because the toolchain is no
  longer supported, g77 support isn't needed anymore and SciPy can now include
  Fortran 90/95 code.

**Linux**

- PyPI-compatible Linux wheels can be produced via the manylinux_ project,
  which is used under the hood by our ``cibuildwheel`` infrastructure.

Other Linux build-setups result in PyPI incompatible wheels, which
would need to be distributed via custom channels, e.g. in a
Wheelhouse_, see at the wheel_ and Wheelhouse_ docs.


.. _`SciPy's configuration file`: https://github.com/scipy/scipy/blob/main/pyproject.toml
.. _Numpy: https://numpy.org
.. _Python: https://www.python.org
.. _nose: https://nose.readthedocs.io
.. _asv: https://asv.readthedocs.org
.. _matplotlib: https://matplotlib.org
.. _Pillow: https://pillow.readthedocs.org
.. _scikits.umfpack: https://pypi.org/project/scikit-umfpack
.. _mpmath: http://mpmath.org
.. _pooch: https://www.fatiando.org/pooch/latest/
.. _Cython: https://cython.org
.. _pybind11: https://github.com/pybind/pybind11
.. _setuptools: https://github.com/pypa/setuptools
.. _wheel: https://wheel.readthedocs.io/
.. _pip: https://pip.pypa.io/en/stable/
.. _Python Packaging User Guide: https://packaging.python.org
.. _Wheelhouse: https://pypi.org/project/Wheelhouse
.. _MingwPy: https://mingwpy.github.io
.. _Sphinx: http://www.sphinx-doc.org/
.. _PyData Sphinx theme: https://pydata-sphinx-theme.readthedocs.io/en/latest/
.. _Sphinx-Design: https://sphinx-design.readthedocs.io
.. _six: https://pypi.org/project/six
.. _decorator: https://github.com/micheles/decorator
.. _manylinux: https://github.com/pypa/manylinux/
.. _threadpoolctl: https://github.com/joblib/threadpoolctl
.. _MyST-NB: https://myst-nb.readthedocs.io/
.. _Hypothesis: https://hypothesis.readthedocs.io/en/latest/index.html
