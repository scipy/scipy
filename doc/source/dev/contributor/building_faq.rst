.. _build-install-faq:

=================
Build/Install FAQ
=================

How do I set up multiple versions of SciPy on my machine?
=========================================================

You may want to set up a development version of SciPy in parallel to a released
version that you use to do your job/research.

If you use the ``conda`` package manager, this is covered in the
:ref:`conda-guide`.

Another simple way to achieve this is to install the released version in
site-packages, by using a binary installer or pip, for example, and set
up the development version in a virtualenv. First, install
`virtualenv`_ (optionally, use `virtualenvwrapper`_), then create your
virtualenv (named scipy-dev here) with::

    $ virtualenv scipy-dev

Now, whenever you want to switch to the virtual environment, you can use the
command ``source scipy-dev/bin/activate``, and ``deactivate`` to exit from the
virtual environment and back to your previous shell. With scipy-dev
activated, first install Scipy's dependencies::

    $ pip install numpy pytest cython pybind11

After that, you can install a development version of Scipy, for example, via::

    $ python setup.py install

The installation goes to the virtual environment.

How do I set up an in-place build for development?
==================================================

For development, you can set up an in-place build so that changes made to
``.py`` files have effect without rebuild. First, run::

    $ python setup.py build_ext -i

Then you need to point your PYTHONPATH environment variable to this directory.
Some IDEs (`Spyder`_, for example) have utilities to manage PYTHONPATH. On Linux
and OSX, you can run the command::

    $ export PYTHONPATH=$PWD

and on Windows::

    $ set PYTHONPATH=/path/to/scipy

Now, editing a Python source file in SciPy allows you to immediately
test and use your changes (in ``.py`` files), by simply restarting the
interpreter.

.. _virtualenv: https://virtualenv.pypa.io/

.. _virtualenvwrapper: https://bitbucket.org/dhellmann/virtualenvwrapper/

.. _Spyder: https://www.spyder-ide.org/

How do I checkout a pull request from GitHub locally?
=====================================================

Second to code review pull requests it is helpful to have a local copy of the
code changes in the pull request. The preferred method to bring a PR from the
github repository to your local repo in a new branch::

    $ git fetch upstream pull/PULL_REQUEST_ID/head:NEW_BRANCH_NAME

The value of ``PULL_REQUEST_ID`` will be the PR number and the
``NEW_BRANCH_NAME`` will be the name of the branch in your local repository
where the diffs will reside.

Now you have a branch in your local development area to code review in Python.

How do I deal with Fortran ABI mismatch?
========================================

Some linear algebra libraries are built with g77 ABI and others with
GFortran ABI, and these two ABIs are incompatible. Therefore, if you
build SciPy with ``gfortran`` and link to a linear algebra library, like
MKL, which is built with g77 ABI, then there'll be an exception or a
segfault. SciPy fixes this by using the CBLAS API for the few
functions in the BLAS API that suffers from this issue.

Note that SciPy needs to know at build time, what needs to be done and
the build system will automatically check whether linear algebra
library is MKL and if so, use the CBLAS API instead of the BLAS API.
If autodetection fails or if the user wants to override this
autodetection mechanism, use the following:

*For ``meson`` based builds (new in 1.9.0):*

Use the ``-Duse-g77-abi=true`` build option. E.g.,::

    $ meson setup build -Duse-g77-abi=true

A more complete example, also configuring the BLAS/LAPACK libraries and picking
a better Python install behavior (this is what conda-forge could be using for
example)::

    $ meson setup builddir -Duse-g77-abi=true -Dblas=blas -Dlapack=lapack -Dpython.install_env=auto
    $ meson install -C builddir

*For ``distutils`` based builds:*

Set the environment variable ``SCIPY_USE_G77_ABI_WRAPPER`` to 0 or 1 to disable
or enable using CBLAS API.

How do I use a custom BLAS distribution on Linux?
=================================================

To customize which BLAS is used, you can set up a ``site.cfg`` file. See the
``site.cfg.example`` file in the numpy source for the options you can set.

Note that Debian and Ubuntu package optimized BLAS libraries in an exchangeable
way. You can install libraries, such as ATLAS or OpenBLAS and change the default
one used via the alternatives mechanism::

    $ sudo apt-get install libopenblas-base libatlas3-base
    $ update-alternatives --list libblas.so.3
    /usr/lib/atlas-base/atlas/libblas.so.3
    /usr/lib/libblas/libblas.so.3
    /usr/lib/openblas-base/libopenblas.so.0

    $ sudo update-alternatives --set libblas.so.3 /usr/lib/openblas-base/libopenblas.so.0

See ``/usr/share/doc/libatlas3-base/README.Debian`` for instructions on how to
build optimized ATLAS packages for your specific CPU. The packaged OpenBLAS
chooses the optimal code at runtime so it does not need recompiling unless the
packaged version does not yet support the used CPU.

You can also use a library you built yourself by preloading it. This does not
require administrator rights::

    LD_PRELOAD=/path/to/libatlas.so.3 ./my-application


Version-specific notes
======================

If you have any problems installing SciPy on your Mac
based on these instructions, please check the `scipy-dev mailing list archives
<https://www.scipy.org/mailing-lists>`__
for possible solutions. If you
are still stuck, feel free to join scipy-dev for further
assistance. Please have the following information ready:

* Your OS version

* The versions of gcc and gfortran and where you obtained gfortran

  * ``$ gcc --version``

  * ``$ gfortran --version``

* The versions of NumPy and SciPy that you are trying to install

* The full output of ``$ python setup.py build``
