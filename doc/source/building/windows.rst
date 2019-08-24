===============================
Building From Source on Windows
===============================

.. contents::
   :local:

Overview
--------

Compared to OSX and Linux, building NumPy and SciPy on Windows is more difficult,
largely due to the lack of compatible, open-source libraries like LAPACK_ or
ATLAS_ that are necessary to build both libraries and have them perform
relatively well. It is not possible to just call a one-liner on the command
prompt as you would on other platforms via ``sudo apt-get install`` machinery.

This document describes one option to build OpenBLAS and SciPy from source that 
was validated for scipy 1.0.0. However, in light of all the work currently being 
done, **do not expect** these instructions to be accurate in the long-run and be
sure to check up on any of the open source projects mentioned for the most up-to-date
information. For more information on all of these projects, the Mingwpy_ website
is an excellent source of in-depth information than this document will provide.

.. _Mingwpy: https://mingwpy.github.io/
.. _ATLAS: http://math-atlas.sourceforge.net/
.. _OpenBLAS: https://github.com/xianyi/OpenBLAS
.. _LAPACK: http://www.netlib.org/lapack/


Building the Released SciPy
---------------------------

This section provides the step-by-step process to build the released scipy. If you want
to build completely from source, you should estimate at least three hours to build all
libraries and compile SciPy. Feel free to stop and inspect any step at any time, but
for this section, we'll just mention the steps without providing an in-depth explanation
for the reasons behind them. If you have further questions about what we're doing, more
in-depth documentation is provided in the sections below. Also, please make sure to read
this section before proceeding, as the presence or absence of error messages in general
is not a good indication of whether you've completed a step correctly. Each step creates
particular files, and what ultimately matters is whether you have built the required files
rather than whether error messages appeared in your terminal.

Building OpenBLAS
=================

First, we need to install the software required to build OpenBLAS_, which is the BLAS_
library that we're going to use. Because the software to build OpenBLAS is different than
that required to build SciPy and because OpenBLAS takes a long time to build, we're going
to start building OpenBLAS first and then explain what to do next while the OpenBLAS build
is running. **Alternatively, if you'd rather download a pre-built OpenBLAS, download the
one of the** `pre-built zip files`_ **and skip to the Installing OpenBLAS section below.**

Otherwise, install MSYS2 using `these instructions`_ including the `pacman` update
instructions. Occasionally during the updates the terminal might ask you to close the
terminal but then might refuse to be closed and hang. If this happens you can kill it via
Task Manager and continue with the instructions. Make sure to install the correct
architecture for the SciPy that you want to build (eg. 32 or 64 bit). Now, you have three
options for opening a terminal which are MSYS2, MINGW (32 or 64 bit). After updating all
the packages, now we are ready to install some more package bundles that we will need. 
Open a MSYS2 terminal and type the following depending on the architecture of your 
choice; run the following for a 32-bit build

.. code:: shell

    pacman -S --needed base-devel mingw-w64-i686-toolchain mingw-w64-i686-cmake git

and for 64-bit

.. code:: shell

    pacman -S --needed base-devel mingw-w64-x86_64-toolchain mingw-w64-x86_64-cmake git

It will prompt to whether install everything in these packages and you can simply accept
all via hitting enter key at each step. 

We should be aware of the fact that these tools also install Python2, very similar to 
a virtual environment, which is only usable within an MSYS2 terminal and we are **not**
going to use it at any point. After updating, now we are going to use the build
toolchain that we have installed in the previous step. Depending on 32/64bit choice,
we will switch to another shell that MSYS2 created. In your start menu you should
see three MSYS2 terminal shortcuts. Select the one with either 64 or 32bit indicator.
The reason why we do this is that the toolchain and compilers are available to these
shells and not to the standard MSYS2 terminal.

If you already have a GitHub repository folder where you keep your own repos, it is better 
to use that location to keep things nice and tidy since we are going to clone yet another 
repository to obtain the source code, hence

.. code:: shell

    cd /c/<wherever the GitHub repo folder is>/GitHub

You don't necessarily need to build in that particular location, but it should be somewhere
convenient. To make sure that we're ready to build, type the following in the terminal:

.. code:: shell

   make
   gfortran
   gcc

These commands should give errors as we have not provided any arguments to them.
However an error also implies that they are accessible on the path. Now clone
the repository required to build OpenBLAS:

.. code:: shell

   git clone https://github.com/matthew-brett/build-openblas.git
   cd build-openblas
   git submodule update --init --recursive

If any of these commands fail, you're not ready to build. Go back and make sure that MSYS2
is installed correctly and has the required packages enabled. Now, let's set some
environment variables. In the MSYS2 terminal, type the following.

.. code:: shell

    export OPENBLAS_COMMIT=5f998ef
    export OPENBLAS_ROOT="C:\\opt"
    export BUILD_BITS=64

Please check these variables' purpose for a moment. More specifically, make sure that
you have read/write access to the path that :code:`OPENBLAS_ROOT` points to. The output of the
OpenBLAS build will  be collected in this folder. Make sure that the :code:`OPENBLAS_COMMIT`
points to the correct OpenBLAS commit that you want to build in the cloned repo. In the
future, :code:`build_openblas` repository might get updated and you might want to get those
updates by changing the commit. Make sure that the architecture is correctly set to either
32 or 64 bit. And after you've made sure of that, start the OpenBLAS build with:

.. code:: shell

    ./build_openblas.sh

Building OpenBLAS is challenging. The build may fail with an error after a few
hours but may also fail silently and produce an incorrect binary. Please, if you
have any issues, `report them`_ so that we can save the next person's time.

While you're waiting on OpenBLAS to finish building, go ahead and install `build tools`_
from Microsoft, since these take a while to install and you'll need them later.

After the :code:`build_openblas.sh` script has completed (probably with an error), there
should be an :code:`openblas.a` file somewhere on your system. If :code:`OPENBLAS_ROOT` was
set to :code:`C:\\opt`, then you might see a line like this in the MSYS2 terminal:

.. code:: shell

   Copying the static library to /c/opt/64/lib

Installing OpenBLAS
===================

If you see that line, then you might have OpenBLAS correctly built, even if other failures
might have occurred. Look in that folder for :code:`openblas.a`. If you find a file called
something like :code:`libopenblas_5f998ef_gcc7_2_0.a`, just rename it to :code:`openblas.a`
and continue. If the file isn't there, then poke around and try to find the file elsewhere
in :code:`OPENBLAS_ROOT`. If you don't have that file, you'll probably need to find out what
happened and then build OpenBLAS again. But if you have that file, we'll assume that you've
completed this step correctly. Proceeding on that assumption, let's build SciPy.

**Before continuing, make sure that you don't have other copies of either**
:code:`openblas.lib` **or** :code:`libopenblas.lib` **on your computer elsewhere.
Multiple copies could result in later build errors that will be difficult to debug.
You may verifiy that the openblas library was correctly picked up by looking for
the following in your build log:**

.. code:: shell

   FOUND:
      libraries = ['openblas']
      library_dirs = ['C:\opt\64\lib']
      language = c
      define_macros = [('HAVE_CBLAS', None)]

Building SciPy
==============

Once you have built OpenBLAS, it's time to build SciPy. Before continuing make sure to
install the following software for building on the latest Python version. For building
on other Python versions, see the WindowsCompilers_ page.

1) Install Microsoft Visual Studio 2015 or 2017 Community Edition (use the `build tools`_
   from Microsoft)
2) Finally, install Python from https://python.org/ (make sure to check the box to install
   pip)

After you've installed the required software, open an MSYS2 terminal, change to a good
location to build, and clone SciPy.

.. code:: shell

   cd C:\Users\MyUser\Downloads
   git clone https://github.com/scipy/scipy.git
   cd scipy
   
Now we need to copy the :code:`openblas.a` file that we've built earlier to the correct
location. If your Python is installed somewhere like the following:

.. code:: shell

   C:\Users\<user name>\AppData\Local\Programs\Python\Python36\python.exe


Then you'll need to put the :code:`openblas.a` file somewhere like the following:

.. code:: shell

   C:\Users\<user name>\AppData\Local\Programs\Python\Python36\Lib

Adjust the location accordingly based on where :code:`python.exe` is located. Now for a
sanity check. Type  the following and press enter.

.. code:: shell

    gfortran

If you see an error with the above command, :code:`gfortran` is not correctly installed.
Go back to the "Building OpenBLAS" section and make sure that you have installed the correct
tools.

Now install the dependencies that we need to build and test SciPy. **It's important that you
specify the full path to the native Python interpreter so that the built-in MSYS2 Python will
not be used. Attempting to build with the MSYS2 Python will not work correctly.**

.. code:: shell

    /c/Users/<user name>/AppData/Local/Programs/Python/Python36/python.exe \
         -m pip install numpy>=1.14.0 cython pytest pytest-xdist

Please note that this is a simpler procedure than what is used for the official binaries.
**Your binaries will only work with the latest NumPy (v1.14.0dev and higher)**. For
building against older NumPy versions, see `Building Against an Older NumPy Version`_.
Make sure that you are in the same directory where  ``setup.py`` is (you should be if you 
have not changed directories):

.. code:: shell

    ls setup.py

Assuming that you have set up everything correctly, you should be ready to build. Run
the following commands:

.. code:: shell

    /c/Users/<user name>/AppData/Local/Programs/Python/Python36/python.exe \
         -m pip wheel -v -v -v .
    /c/Users/<user name>/AppData/Local/Programs/Python/Python36/python.exe \
         runtests.py --mode full

Congratulatations, you've built SciPy!

.. _BLAS: https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms
.. _OpenBLAS: https://github.com/xianyi/OpenBLAS
.. _`these instructions`: https://github.com/msys2/msys2/wiki/MSYS2-installation
.. _`build tools`: https://www.visualstudio.com/downloads/#build-tools-for-visual-studio-2017
.. _`report them`: https://github.com/scipy/scipy/issues/new
.. _`pre-built zip files`: https://3f23b170c54c2533c070-1c8a9b3114517dc5fe17b7c3f8c63a43.ssl.cf2.rackcdn.com/
.. _WindowsCompilers: https://wiki.python.org/moin/WindowsCompilers

Building Against an Older NumPy Version
---------------------------------------

If you want to build SciPy to work with an older numpy version, then you will need 
to replace the NumPy "distutils" folder with the folder from the latest numpy.
The following powershell snippet can upgrade NumPy distutils while retaining an older
NumPy ABI_.

.. code:: shell

      $NumpyDir = $((python -c 'import os; import numpy; print(os.path.dirname(numpy.__file__))') | Out-String).Trim()
      rm -r -Force "$NumpyDir\distutils"
      $tmpdir = New-TemporaryFile | %{ rm $_; mkdir $_ }
      git clone -q --depth=1 -b master https://github.com/numpy/numpy.git $tmpdir
      mv $tmpdir\numpy\distutils $NumpyDir

.. _ABI: https://en.wikipedia.org/wiki/Application_binary_interface

Additional Resources
--------------------

As discussed in the overview, this document is not meant to provide extremely detailed explanations on how to build
NumPy and SciPy on Windows. This is largely because currently, there is no single superior way to do so
and because the process for building these libraries on Windows is under development. It is likely that any
information will go out of date relatively soon. If you wish to receive more assistance, please reach out to the NumPy
and SciPy mailing lists, which can be found `here <https://www.scipy.org/scipylib/mailing-lists.html>`__.  There are many
developers out there, working on this issue right now, and they would certainly be happy to help you out!  Google is also
a good resource, as there are many people out there who use NumPy and SciPy on Windows, so it would not be surprising if
your question or problem has already been addressed.
