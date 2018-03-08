===============================
Building From Source on Mac OSX
===============================

These instructions describe how to build NumPy and SciPy libraries from
source.

If you just want to use NumPy or SciPy, install pre-built binaries as described
in :doc:`/install`.

Python
------

Apple ships its own version of Python with OS X. However, we
*strongly* recommend installing the `official Python distribution
<https://www.python.org/downloads/>`__.

Alternatively, use Python from one of the OS X package managers 
(Homebrew, MacPorts, Fink).  

Apple's Developer Tools
-----------------------

Apple's Developer Tools provide a number of key libraries,
particularly the `vecLib Framework
<https://developer.apple.com/library/mac/documentation/Performance/Conceptual/vecLib/>`__
, which includes the BLAS and LAPACK libraries for optimizing matrix and vector
operations. The most recent version should be included on your OS X installation CD.
Ensure that all components are installed by choosing **customize** when available
during the install process and selecting all optional packages - at least the X11
development tools and (on OS X 10.6 or lower) the 10.4 SDK.

Compilers (C/C++/FORTRAN/Cython)
--------------------------------

Though virtually any commercial C/C++ compiler may be used with SciPy, OS X
comes with GNU C compilers pre-installed. The only thing missing is the GNU
FORTRAN compiler.

We recommend gfortran; this is a free, open source, F95 compiler. We suggest you
use the following binaries:

* gfortran installed via `Homebrew <https://brew.sh/>`__, or,
* http://r.research.att.com/tools/gcc-42-5666.3-darwin11.pkg (for Xcode
  4.2 or higher)

See `this site <http://r.research.att.com/tools/>`__ for the most recent links.

Unless you are building from released source packages, the `Cython
<http://cython.org/>`__ compiler is also needed.

Version-specific notes
----------------------

This section notes only things specific to one version of OS X or Python. 
The build instructions in :ref:`Obtaining and Building NumPy and SciPy
<osx-obtaining-and-building>` apply to all versions.

.. _osx-obtaining-and-building:

Obtaining and Building NumPy and SciPy
--------------------------------------

You may install NumPy and SciPy either by checking out the source
files from the Git repositories, or unpacking them from a source
archive file from :doc:`/scipylib/download`. If you choose the latter,
simply expand the archive (generally a gzipped tar file), otherwise
check out the following branches from the repository:

::

       $ git clone https://github.com/numpy/numpy.git
       $ git clone https://github.com/scipy/scipy.git

Both NumPy and SciPy are built as follows:

::

       $ python setup.py build
       $ python setup.py install

The above applies to the `official Python distribution
<https://www.python.org/downloads/>`__, which is 32-bit
only for 2.6 while 32/64-bit bundles are available for 2.7 and
3.x. For alternative 64-bit Pythons (either from Apple or home-built)
on Snow Leopard, you may need to extend your build flags to specify
the architecture by setting LDFLAGS and FFLAGS.

Note that with distutils (setup.py) given build flags like LDFLAGS
**do not extend but override the defaults**, so you have to specify
all necessary flags. Only try this if you know what you're doing!

After a successful build, you may try running the built-in unit tests
for SciPy:

::

       $ python
       >>> import numpy as np
       >>> np.test('full')
       >>> import scipy
       >>> scipy.test()

Be sure not to import numpy or scipy while you're in the numpy/scipy
source tree. Change directory first.

If you have any problems installing SciPy on your Mac
based on these instructions, please check the :doc:`scipy-users and
scipy-dev mailing list archives
</scipylib/mailing-lists>` for possible solutions. If you
are still stuck, feel free to join scipy-users for further
assistance. Please have the following information ready:

* Your OS version

* The versions of gcc and gfortran and where you obtained gfortran

  * ``$ gcc --version``

  * ``$ gfortran --version``

* The versions of numpy and scipy that you are trying to install

* The full output of ``$ python setup.py build``
