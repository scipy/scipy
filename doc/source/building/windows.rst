===============================
Building From Source on Windows
===============================

.. contents::
   :local:

Overview
--------

Compared to OSX and Linux, building NumPy and SciPy on Windows is more 
difficult, largely due to the lack of compatible, open-source libraries like 
BLAS/LAPACK_ and open-source compilers that are necessary to build both
libraries and have them perform relatively well. It is not possible to just call
a one-liner on the command prompt as you would on other platforms via ``sudo apt
install`` machinery.

This document describes one option to build OpenBLAS and SciPy from source 
that was validated for scipy 1.4.0. However, in light of all the work 
currently being done, please **do not expect** these instructions to be accurate
in the long-run and be sure to check up on any of the open source projects
mentioned for the most up-to-date information. For more information on all of
these projects, the Mingwpy_ website is an excellent source of in-depth 
information than this document will provide.

.. _Mingwpy: https://mingwpy.github.io/
.. _OpenBLAS: https://github.com/xianyi/OpenBLAS
.. _LAPACK: http://www.netlib.org/lapack/


Building the Released SciPy
---------------------------

This section provides the step-by-step process to build the released scipy. 
If you want to build completely from source, you should estimate at least 
three hours to build all libraries and compile SciPy. Feel free to stop and 
inspect any step at any time, but for this section, we'll just mention the 
steps without providing an in-depth explanation for the reasons behind them. 
If you have further questions about what we're doing, more in-depth
documentation is provided in the sections below. Also, please make sure to
read this section before proceeding, as the presence or absence of error
messages in general is not a good indication of whether you've completed a
step correctly. Each step creates particular files, and what ultimately
matters is whether you have built the required files rather than whether
error messages appeared in your terminal.

Building OpenBLAS
=================

First, we need to install the software required to build OpenBLAS_, which is 
the BLAS_ library that we're going to use. Because the software to build 
OpenBLAS is different than that required to build SciPy and because OpenBLAS 
takes a long time to build, we're going to start building OpenBLAS first and 
then explain what to do next while the OpenBLAS build is running. 
**Alternatively, if you'd rather download a pre-built OpenBLAS, download the 
one of the** `pre-built zip files`_ **and skip to the Installing OpenBLAS 
section below.**. However it is also likely that your version of Windows and the
compiler you wish to use won't be compatible with what these prebuilt binaries
produced with. This is still one of the main pain points of building for
Windows. That's why we will attempt to build our own OpenBLAS.

First install MSYS2 using `these instructions`_ including the `pacman` 
update instructions. Occasionally during the updates the terminal might ask 
you to close the terminal but then might refuse to be closed and hang. If 
this happens you can kill it via Task Manager and continue with the 
instructions. Make sure to install the correct architecture for the SciPy
that you want to build (eg. 32 or 64 bit). Now, in your start menu, you have
three options for opening a terminal which are MSYS2, MSYS2 MinGW (32/64 bit).
Choose the one that applies to you and keep running the same command below until
you get all components are up-to-date.

.. code:: shell

    pacman -Syuu

Now, next step is to install some more package bundles that we will need. Open
a MSYS2 terminal and type the following depending on the architecture of your
choice; run the following for the common 64-bit build

.. code:: shell

    pacman -S --needed base-devel mingw-w64-x86_64-toolchain mingw-w64-x86_64-cmake git

and for 32-bit run instead

.. code:: shell

    pacman -S --needed base-devel mingw-w64-i686-toolchain mingw-w64-i686-cmake git

If you are not sure which one you want, choose 64-bit option in every step. 

It will prompt to whether install everything in these packages and you can 
simply accept all via hitting enter key at each step which also takes some time.

Once you install everything, just to have a fresh start close and reopen the
same terminal.

If you already have a GitHub repository folder where you keep your own repos, 
it is better to use that location to keep things nice and tidy since we are 
going to clone yet another repository to obtain the source code, hence

.. code:: shell

    cd /c/<wherever you want to clone OpenBLAS into>/

It should be somewhere convenient and with write permissions. To make sure that
we're ready to build, type the following in the terminal one-by-one:

.. code:: shell

   make
   gfortran
   gcc
   git

Each of these commands should give errors as we have not provided any arguments
to them. However an error also implies that they are accessible on the path
which is what we wanted to test. If any of these commands fail, you're not ready
to build. Go back and make sure that MSYS2 is installed correctly and has the
required packages enabled. Now time to clone the OpenBLAS repository somewhere
convenient.

.. code:: shell

   cd /c/Users/<user name>/Documents/GitHub
   git clone https://github.com/xianyi/OpenBLAS.git
   cd OpenBLAS
   git submodule update --init --recursive

Now change the directory one level up via :code:`cd..` to get out of the
directory and create a file named `build_openblas.sh`. The easiest way is to
type

.. code:: shell

    touch build_openblas.sh

Then open this file in a text editor like Notepad++ and paste the following
content in this empty file: 

.. code:: shell

    # Adjust the following to your liking and your MSYS installation status
    OPENBLAS_ROOT=/c/opt
    BUILD_BITS=64

    # Print some gcc info that MSYS2 discovered in the path
    which gcc
    gcc --version

    # Get into the repository that we cloned
    cd OpenBLAS

    # Change the following to a specific branch/tag/release you wish
    # Consult the git manual to learn more about your options
    git checkout tags/v0.3.7 -b v0.3.7

    # The following two lines clean up in case we make a mistake and need
    # to run the script again
    git clean -fxd
    rm -rf $OPENBLAS_ROOT/$BUILD_BITS

    # Set architecture flags
    march="x86-64"
    extra="-fno-asynchronous-unwind-tables"
    vc_arch="X64"
    cflags="-O2 -march=$march -mtune=generic $extra"
    fflags="$cflags -frecursive -ffpe-summary=invalid,zero"

    # Build name for output library from gcc version and OpenBLAS commit.
    GCC_TAG="gcc_$(gcc -dumpversion | tr .- _)"
    OPENBLAS_VERSION=$(git describe --tags)
    # Build OpenBLAS
    # Variable used in creating output libraries
    export LIBNAMESUFFIX=${OPENBLAS_VERSION}-${GCC_TAG}
    make BINARY=$BUILD_BITS DYNAMIC_ARCH=1 USE_THREAD=1 USE_OPENMP=0 NO_WARMUP=1 BUILD_LAPACK_DEPRECATED=1 COMMON_OPT="$cflags" FCOMMON_OPT="$fflags"
    make install PREFIX=$OPENBLAS_ROOT/$BUILD_BITS

This is the automation script that will make sure the right variables are used
in the right place. Linux users are very familiar to such scripts but for
Windows users it might be a bit awkward. You can think of these as ``.bat``
files. You can change the variables to your situation. After you've created this
file and you are one directory up the OpenBLAS repo of that, start the OpenBLAS
build with:

.. code:: shell

    ./build_openblas.sh

Building OpenBLAS is challenging. The build may fail with an error after a 
few hours but may also fail silently and produce an incorrect binary. Please, 
if you have any issues, `report them`_ so that we can save the next person's 
time.

While you're waiting on OpenBLAS to finish building, go ahead and install
`build tools`_ from Microsoft, since these take a while to install and you'll 
need them later.

After the :code:`build_openblas.sh` script has completed, there should be an
:code:`libopenblas.....a` as a resulting artifact. If :code:`OPENBLAS_ROOT` was
set to :code:`C:\\opt`, then you might see a line like this in the MSYS2
terminal:

.. code:: shell

   Copying the static library to /c/opt/64/lib
   
If you, by any chance, receive the following error

.. code::shell

    <command-line>:0:4: error: expected identifier or '(' before numeric constant
    
that means you have some header file definition clash and you have to downgrade
certain items. See this
`OpenBLASwiki <https://github.com/xianyi/OpenBLAS/wiki/How-to-use-OpenBLAS-in-Microsoft-Visual-Studio#build-openblas-on-windows-os>`__
page to read on which ones and how to do it. This should be sufficient for you
to build OpenBLAS.

Installing OpenBLAS
===================

If you see the last line mentioning the static library copy, then you might have
OpenBLAS correctly built, even if other failures might have occurred. Look in
the folder you used as a parameter to :code:`OPENBLAS_ROOT\64\lib`.
If you find a file called something like
:code:`libopenblas_v0.3.7-gcc_9_2_0p-r0.2.20.a`, just make a copy and rename it
to :code:`openblas.a`.

If you don't have that file, you'll probably need to find
out what happened and then build OpenBLAS again. We know this is **very** annoying
however unfortunately we have no other alternatives. But if you have that file,
we'll assume that you've completed this step correctly. Proceeding on that
assumption, let's build SciPy.

**Before continuing, make sure that you don't have other copies of either**
:code:`openblas.lib` **or** :code:`libopenblas.lib` **on your computer elsewhere.
Multiple copies could result in later build errors that will be difficult to debug.
You may verify that the openblas library was correctly picked up by looking for
the following in your build log:**

.. code:: shell

   FOUND:
      libraries = ['openblas']
      library_dirs = ['C:\...........\lib']
      language = c
      define_macros = [('HAVE_CBLAS', None)]

Building SciPy
==============

Once you have built OpenBLAS, it's time to build SciPy. Before continuing make
sure to install the following software for building on the latest Python
version. For building on other Python versions, see the WindowsCompilers_ page.

1) Install Microsoft Visual Studio 2015 or 2017 Community Edition (use the
   `build tools`_ from Microsoft)
2) Finally, install Python from https://www.python.org/ (make sure to check
   the box to add Python to path)

Just like before pick a convenient place to clone SciPy. Next to OpenBLAS is
often a convenient option (note: not inside OpenBLAS folder but next to).

.. code:: shell

   cd C:\Users\Ilhan\Documents\GitHub
   git clone https://github.com/scipy/scipy.git
   cd scipy
   
Now we need to copy the :code:`openblas.a` file that we've built earlier to the
correct location. If your Python is installed somewhere like the following:

.. code:: shell

   C:\Users\<user name>\AppData\Local\Programs\Python\Python37\python.exe


Then you'll need to put the :code:`openblas.a` file that we previously copied and
renamed somewhere like the following:

.. code:: shell

   C:\Users\<user name>\AppData\Local\Programs\Python\Python38\Lib

Adjust the location accordingly based on where :code:`python.exe` is located.

At this stage, we are done with the OpenBLAS part and hopefully we will not need
to build OpenBLAS anytime soon. But we tend to build SciPy more often as it is
on a quicker release cycle. Hence it makes sense to use Windows ``cmd`` or 
Powershell for the the build as it is a more native tool. This requires placing
the MinGW compilers on the path. Hence for a sanity check type  the following on
the command prompt or Powershell

.. code:: shell

    gfortran

If you see a missing command error with the above, :code:`gfortran` is not
correctly installed or not on the path. Hence make sure that the following 
folder (or the folder you have installed MSYS to) is on the system path
variables sufficiently high.

.. code:: shell

    C:\MSYS64\MINGW64\BIN

Now install the dependencies that we need to build and test SciPy. 

.. code:: shell

    python -m pip install wheel setuptools numpy>=1.14.0 Cython>=0.29.13 pybind11>=2.2.4 pytest pytest-xdist

The last two are for using SciPy's test suite which is handy if you want to test
some new change locally.

Please note that this is a simpler procedure than what is used for the official binaries.
**Your binaries will only work with the latest NumPy (v1.14.0dev and higher)**. For
building against older NumPy versions, see `Building Against an Older NumPy Version`_.

Assuming that you are in the top of the SciPy repository directory where
``setup.py`` is and assuming that you have set up everything correctly, you
are ready to build. Run the following commands:

.. code:: shell

    python setup.py build

When everything finishes without an error, congratulatations, you've built
SciPy!

You can further install the build SciPy via 

.. code:: shell

    python setup.py install

Just make sure that you uninstalled the existing installation of other SciPy if
there were any.


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
