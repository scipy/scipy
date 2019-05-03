#############################
Building From Source on Linux
#############################

====================
Generic instructions
====================

To build NumPy/SciPy from source, get the `source package
<https://github.com/scipy/scipy>`__, unpack it, and:

::

   python setup.py install --user   # installs to your home directory

or

::

   python setup.py build
   python setup.py install --prefix=$HOME/local

Before building, you will also need to install packages that NumPy and
SciPy depend on

* BLAS and LAPACK libraries (optional but strongly recommended for
  NumPy, required for SciPy): typically `ATLAS
  <http://math-atlas.sourceforge.net/>`__ + `OpenBLAS
  <https://github.com/xianyi/OpenBLAS/>`__, or `MKL
  <https://software.intel.com/en-us/intel-mkl>`__.

* C and Fortran compilers (typically ``gcc`` and ``gfortran``).

* Python header files (typically a package named ``python-dev`` or ``python-devel``)

* Unless you are building from released source packages, the `Cython
  <http://cython.org/>`__ compiler is necessary (typically in a
  package named ``cython``). For building recent SciPy, it is possible
  that you need Cython in a newer version than is available in your
  distribution.

Typically, you will want to install all of the above from packages
supplied by your Linux distribution, as building them yourself is
complicated. If you need to use specific BLAS/LAPACK libraries, you
can do

::

   export BLAS=/path/to/libblas.so
   export LAPACK=/path/to/liblapack.so
   export ATLAS=/path/to/libatlas.so
   python setup.py ............

If you don't want to any LAPACK, just do "``export LAPACK=``".

You will find below additional installation instructions and advice
for many major Linux distributions.


=====================
Specific instructions
=====================

.. contents::
   :local:


Debian / Ubuntu
===============

To build from source the following packages are needed::

   sudo apt-get install gcc gfortran python-dev libopenblas-dev liblapack-dev cython

To customize which BLAS is used, you can setup a `site.cfg` file.  See
the `site.cfg.example` file in the numpy source for the options you
can set.

Note that Debian and Ubuntu package optimized BLAS libraries in a
exchangeable way.  You can install libraries such as ATLAS or OpenBLAS
and change the default one used via the alternatives mechanism:

::

    $ sudo apt-get install libopenblas-base libatlas3-base
    $ update-alternatives --list libblas.so.3
    /usr/lib/atlas-base/atlas/libblas.so.3
    /usr/lib/libblas/libblas.so.3
    /usr/lib/openblas-base/libopenblas.so.0

    $ sudo update-alternatives --set libblas.so.3 /usr/lib/openblas-base/libopenblas.so.0

See /usr/share/doc/libatlas3-base/README.Debian for instructions on
how to build optimized ATLAS packages for your specific CPU.  The
packaged OpenBLAS chooses the optimal code at runtime so it does not
need recompiling unless the packaged version does not yet support the
used CPU.

You can also use a library you built yourself by preloading it. This does not
require administrator rights.

::

    LD_PRELOAD=/path/to/libatlas.so.3 ./my-application


Fedora 26
=========

To install scipy build requirements, you can do::

    sudo dnf install gcc-gfortran python3-devel python2-devel openblas-devel lapack-devel Cython


Intel C compiler and MKL
========================

Intel MKL 11.0 (updated Dec 2012)
---------------------------------

Add the following lines to site.cfg in your top level NumPy directory
to use Intel® MKL for Intel® 64 (or earlier known as em64t)
architecture, considering the default installation path of Intel® MKL
which is bundled with Intel® Composer XE SP1 version on Linux:

::

   [mkl]
   library_dirs = /opt/intel/composer_xe_2013/mkl/lib/intel64
   include_dirs = /opt/intel/composer_xe_2013/mkl/include
   mkl_libs = mkl_intel_lp64,mkl_intel_thread,mkl_core

If you are building NumPy for 32 bit, please add as the following

::

   [mkl]
   library_dirs = /opt/intel/composer_xe_2013/mkl/lib/ia32
   include_dirs = /opt/intel/composer_xe_2013/mkl/include
   mkl_libs = mkl_intel,mkl_intel_thread,mkl_core

Instead of the layered linking approach for the Intel® MKL as shown
above, you may also use the dynamic interface lib mkl_rt.lib. So, for
both the ia32 and intel64 architecture make the change as below

::

   mkl_libs = mkl_rt

Modify cc_exe in numpy/numpy/distutils/intelccompiler.py to be
something like:

::

   cc_exe = 'icc -O2 -g -openmp -avx'

Here we use, default optimizations (-O2), OpenMP threading (-openmp)
and Intel® AVX optimizations for Intel® Xeon E5 or E3 Series which are
based on Intel® SandyBridge Architecture (-avx).  Run icc --help for
more information on processor-specific options.

Compile and install NumPy with the Intel compiler (on 64-bit platforms replace "intel" with "intelem"):

::

   python setup.py config --compiler=intel build_clib --compiler=intel build_ext --compiler=intel install

Compile and install SciPy with the Intel compilers (on 64-bit
platforms replace "intel" with "intelem"):

::

   python setup.py config --compiler=intel --fcompiler=intel build_clib --compiler=intel --fcompiler=intel build_ext --compiler=intel --fcompiler=intel install

You'll have to set LD_LIBRARY_PATH to Intel® MKL libraries (exact
values will depend on your architecture, compiler and library
versions) and OpenMP library for NumPy to work.  If you build NumPy
for Intel® 64 platforms:

::

   $export LD_LIBRARY_PATH=/opt/intel/composer_xe_2013/mkl/lib/intel64: /opt/intel/composer_xe_2013/compiler/lib/intel64:$LD_LIBRARY_PATH

If you build NumPy for ia32 bit platforms:

::

   $export LD_LIBRARY_PATH=/opt/intel/composer_xe_2013/mkl/lib/ia32: /opt/intel/composer_xe_2013/compiler/lib/ia32:$LD_LIBRARY_PATH
