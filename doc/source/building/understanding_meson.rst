Understanding Meson
===================

Building SciPy relies on the following tools, which can be considered part of
the build system:

- ``meson``: the Meson build system, installable as a pure Python package from
  PyPI or conda-forge
- ``ninja``: the build tool invoked by Meson to do the actual building (e.g.
  invoking compilers). Installable also from PyPI (on all common platforms) or
  conda-forge.
- ``pkg-config``: the tool used for discovering dependencies (in particular
  BLAS/LAPACK). Available on conda-forge (and Homebrew, Chocolatey, and Linux
  package managers), but not packaged on PyPI.
- ``meson-python``: the Python build backend (i.e., the thing that gets invoked
  via a hook in ``pyproject.toml`` by a build frontend like ``pip`` or
  ``pypa/build``). This is a thin layer on top of Meson, with as main roles (a)
  interface with build frontends, and (b) produce sdists and wheels with valid
  file names and metadata.

Building with Meson happens in stages:

- A configure stage (``meson setup``) to detect compilers, dependencies and
  build options, and create the build directory and ``build.ninja`` file,
- A compile stage (``meson compile`` or ``ninja``), where the extension modules
  that are part of a built SciPy package get compiled,
- An install stage (``meson install``) to install the installable files from
  the source and build directories to the target install directory,

Meson has a good build dependency tracking system, so invoking a build for a
second time will rebuild only targets for which any sources or dependencies
have changed.


To learn more about Meson
-------------------------

Meson has `very good documentation <https://mesonbuild.com/>`__;
it pays off to read it, and is often the best source of answers for "how to do
X". Furthermore, an extensive pdf book on Meson can be obtained for free at
https://nibblestew.blogspot.com/2021/12/this-year-receive-gift-of-free-meson.html

To learn more about the design principles Meson uses, the recent talks linked
from `mesonbuild.com/Videos <https://mesonbuild.com/Videos.html>`__ are also a
good resource.


Explanation of build stages
---------------------------

*This is for teaching purposes only; there should be no need to execute these
stages separately. The dev.py scripts in the root of the repo also contains
these steps and may be studied for insights.*

Assume we're starting from a clean repo and a fully set up conda environment::

  git clone git@github.com:scipy/scipy.git
  git submodule update --init
  mamba env create -f environment.yml
  mamba activate scipy-dev

To now run the configure stage of the build and instruct Meson to put the build
artifacts in ``build/`` and a local install under ``build-install/`` relative
to the root of the repo, do::

  meson setup build --prefix=$PWD/build-install

To then run the compile stage of the build, do::

  ninja -C build

In the command above, ``-C`` is followed by the name of the build directory.
You can have multiple build directories at the same time. Meson is fully
out-of-place, so those builds will not interfere with each other. You can for
example have a GCC build, a Clang build and a debug build in different
directories.

To then install SciPy into the prefix (``build-install/`` here, but note that
that's just an arbitrary name we picked here)::

  meson install -C build

It will then install to ``build-install/lib/python3.11/site-packages/scipy``,
which is not on your Python path, so to add it do (*again, this is for learning
purposes, using ``PYTHONPATH`` explicitly is typically not the best idea*)::

  export PYTHONPATH=$PWD/build-install/lib/python3.11/site-packages/

Now we should be able to import ``scipy`` and run the tests. Remembering that
we need to move out of the root of the repo to ensure we pick up the package
and not the local ``scipy/`` source directory::

  cd doc
  python -c "from scipy import constants as s; s.test()"

The above runs the tests for a single module, ``constants``. Other ways of
running the tests should also work, for example::

  pytest --pyargs scipy

The full test suite should pass, without any build warnings on Linux (with the
GCC version for which ``-Werror`` is enforced in CI at least) and with at most
a moderate amount of warnings on other platforms.
