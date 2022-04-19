.. _build-windows:

===============================
Building from source on Windows
===============================

.. contents::
   :local:

Overview
--------

Compared to OSX and Linux, building NumPy and SciPy on Windows is more
difficult, largely due to the lack of compatible, open-source libraries like
BLAS/LAPACK_ and open-source compilers that are necessary to build both
libraries and have them perform relatively well. It is not possible to just
call a one-liner on the command prompt as you would on other platforms via
``sudo apt install`` machinery.

This document describes one option to build OpenBLAS and SciPy from source
that was validated for SciPy 1.6.0. However, in light of all the work
currently being done, please **do not expect** these instructions to be
accurate in the long-run and be sure to check up on any of the open source
projects mentioned for the most up-to-date information. For more information
on all of these projects, the Mingwpy_ website is an excellent source of
in-depth information than this document will provide.

.. _Mingwpy: https://mingwpy.github.io/
.. _OpenBLAS: https://github.com/xianyi/OpenBLAS
.. _LAPACK: http://www.netlib.org/lapack/


Building the Released SciPy
---------------------------

This section provides the step-by-step process to build the released SciPy.
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
section below.**. However it is also likely that your version of Windows and
the compiler you wish to use won't be compatible with what these prebuilt
binaries produced. This is still one of the main pain points of building
for Windows. That's why we will attempt to build our own OpenBLAS.

We start by installing the MSYS2 platform, on which the OpenBLAS build will take
place. First, download the MSYS2 installer from `msysintaller`_ via choosing
32 or 64 bit. Make sure to install the correct architecture for the SciPy
that you want to build (e.g., 32 or 64 bit). If you are not sure which one to use,
proceed with 64bit. Please follow the installation instructions carefully,
especially step 6 and 7 to update all components.

.. note::

    Occasionally,
    during the updates, the terminal might ask you to close the terminal but then
    might refuse to be closed and hang. If this happens, you can kill it via Task
    Manager and continue with the instructions.

Now, the next step is to install some more package bundles that we will need. Open
a MSYS2 **MinGW** (64 or 32 bit) terminal and type the following depending on the
architecture of your choice; run the following for the common 64-bit build

.. code:: shell

    pacman -S --needed base-devel mingw-w64-x86_64-toolchain mingw-w64-x86_64-cmake git

and for 32-bit run instead

.. code:: shell

    pacman -S --needed base-devel mingw-w64-i686-toolchain mingw-w64-i686-cmake git

Again, if you are not sure which one you want, choose 64-bit option in every
step.

It will prompt to whether install everything in these packages and you can
simply accept all via hitting enter key at each step which also takes some time
to complete. Once you install everything, close and
reopen the MSYS2 MinGW terminal.

If you already have a GitHub repository folder where you keep your own repos,
it is better to use that location to keep things nice and tidy since we are
going to clone yet another repository to obtain the source code. It should be
somewhere convenient and with write permissions. If this is your first time then
you can pick "Documents\GitHub" as a viable option. We will assume that you
picked this folder in the rest of this document. You can create a folder in "My
Documents" using Windows Explorer. To make sure that we're ready to build, type
the following in the terminal one-by-one:

.. code:: shell

   make
   gfortran
   gcc
   git

Each of these commands should fail as we have not provided any arguments
to them. However, an explicit failure from the program rather than from
the command prompt implies that the program is accessible on the path,
which is what we wanted to test. In turn, if an error about the command being
not found is returned, then installation of the packages didn't complete
successfully. If any of these are missing, you're not ready to build. Go back
and make sure that MSYS2 is installed correctly and has the required packages
enabled.

Now it's time to clone the OpenBLAS repository somewhere convenient. Run the
following line-by-line separately, modifying the path to your GitHub repo
folder as appropriate.

.. code:: shell

   cd C:\Users\<user name>\Documents\GitHub
   git clone https://github.com/xianyi/OpenBLAS.git
   cd OpenBLAS
   git submodule update --init --recursive
   git fetch --all --tags --prune

Now we are going to switch to a release of our choice. At the time of writing,
the newest OpenBLAS release version is 0.3.7, hence we will use that.

.. code:: shell

   git checkout tags/v0.3.7 -b v0.3.7

You can see all available options via

.. code:: shell

   git tag

Now change the directory one level up via :code:`cd ..` to get out of the
directory and create a file named `build_openblas.sh`. The easiest way is to
type

.. code:: shell

    touch build_openblas.sh

Of course, you can still also use Windows Explorer to create a new txt file at
that location and then rename it. So the resulting structure would be

.. code:: shell

    my repo folder
        ├─── build_openblas.sh
        ├─── OpenBLAS
                ├─── ...

Then open this file in any text editor, like Notepad++, and paste the following
content in this empty file:

.. code:: shell

    # You may adjust to your preferred output directory
    OPENBLAS_ROOT=/c/opt

    # Adjust to match the MSYS2 version you installed
    BUILD_BITS=64

    # Print some gcc info that MSYS2 discovered in the path
    which gcc
    gcc --version

    # Get into the repository that we cloned
    cd OpenBLAS

    # The following two lines clean up in case we make a mistake and need
    # to run the script again
    git clean -fxd
    git reset --hard
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
    make BINARY=$BUILD_BITS DYNAMIC_ARCH=1 USE_THREAD=1 USE_OPENMP=0 \
        NO_WARMUP=1 BUILD_LAPACK_DEPRECATED=1 \
        COMMON_OPT="$cflags" FCOMMON_OPT="$fflags"
    make install PREFIX=$OPENBLAS_ROOT/$BUILD_BITS

This is the automation script that will make sure the right variables are used
in the right place. Linux users are very familiar to such scripts, but for
Windows users it might be a bit awkward. You can think of these as ``.bat``
files. The script should work as-in for MSYS2 64-bit, but you can change the
variables to your situation as needed. After you've created
this file and you are one directory up the OpenBLAS repo of that, start the
OpenBLAS build with:

.. code:: shell

    ./build_openblas.sh

Building OpenBLAS is challenging and time-consuming. The build may fail with an
error after a few hours but may also fail silently and produce an incorrect
binary. Please, if you have any issues, `report them`_ so that we can save the
next person's time.

One of the known issues is the following; if you, by any chance, receive the
following error

.. code:: shell

    <command-line>:0:4: error: expected identifier or '(' before numeric constant

that means you have some header file definition clash and you have to downgrade
certain items. This is not related to SciPy but let's attempt to provide a
solution. See this
`OpenBLASwiki <https://github.com/xianyi/OpenBLAS/wiki/How-to-use-OpenBLAS-in-Microsoft-Visual-Studio#build-openblas-on-windows-os>`__
page to read on which packages to downgrade and how to do it.
Basically, it involves downloading three files. Then in the MSYS terminal
change the directory to the place where you downloaded the files and run the
commands given in the wiki link. Then come back to the script directory where
`./build_openblas.sh` lives and try again. This should be sufficient for you to
build OpenBLAS.

While you're waiting on OpenBLAS to finish building, go ahead and install
`build tools`_ from Microsoft, since these take a while to install and you'll
need them later.

After the :code:`build_openblas.sh` script has completed, there should be an
:code:`libopenblas.....a` as a resulting artifact. If :code:`OPENBLAS_ROOT` was
set to :code:`C:\\opt`, then you might see a line like this in the MSYS2
terminal:

.. code:: shell

   Copying the static library to /c/opt/64/lib

This is very good news: you have successfully built OpenBLAS!


Installing OpenBLAS
===================

Look for the ``lib`` folder in the folder you used as a parameter to
:code:`OPENBLAS_ROOT`. (It's ``/c/opt/64/lib`` if you didn't change anything in
the script.) You will find three ``.a`` files such as (the names can differ):

.. code:: shell

    libopenblas_v0.2.20-2-g5f998efd-gcc_9_2_0.a
    libopenblas_v0.2.20-2-g5f998efd-gcc_9_2_0.dll.a
    libopenblas_v0.2.20-2-g5f998efd-gcc_9_2_0.p-r0.2.20.a

From these three we are interested only in the first one. Just make a copy and
rename it to :code:`openblas.a`.

If you don't have that file, you'll probably need to find
out what happened and then build OpenBLAS again. We know this is **very**
annoying, but unfortunately we have no other alternatives. The first place
to look for is inside the OpenBLAS directory. If the build succeeds but (for
some reason) auto-moving files to :code:`OPENBLAS_ROOT` fails, the artifacts
will stay inside the OpenBLAS repo
folder. But if you have that file, that's great and we'll assume that you've
completed this step correctly. Proceeding on that assumption, let's build
SciPy.

Before continuing, make sure that you don't have other copies of either
:code:`openblas.a` or :code:`libopenblas.a` from previous attempts or via
previous downloads. Multiple copies could result in later build errors that
will be difficult to debug. If this is the first attempt, you don't need to
worry about this step.

Building SciPy
==============

Once you have built OpenBLAS, it's time to build SciPy. Before continuing, make
sure to install the following software for building on the latest Python
version. For building on other Python versions, see the WindowsCompilers_ page.
We are also assuming that your Python is on the system path. That is to say,
when you type `python` in the Windows command prompt the correct Python is
executed.

Install Microsoft Visual Studio 2017 or 2019 Community Edition (use the
`build tools`_ from Microsoft). If you feel that it is too bloated to install
everything in that bundle (which we do feel a bit so) then here are a subset
which are tested during the build of SciPy 1.6.0 and VS 2019. You can switch
to the individual items view at the top and select only the following

.. code:: shell

    C++ Core Features
    Windows Universal C Runtime
    MSVC v142 - VS 2019 C++ x64/x86 build tools (...)
    Windows 10 SDK (10.0.18362.0)
    C++ 2019 Redistributable Update
    C++ Clang-cl for 142 build tools (x64/x86)
    C++ Clang Compiler for Windows (8.0.1)

Just like before, pick a convenient place to
clone SciPy. Next to OpenBLAS is often a convenient option (note: not inside
OpenBLAS folder but next to). Continuing the example from above

.. code:: shell

    my repo folder
        ├─── build_openblas.sh
        ├─── OpenBLAS
        ├─── SciPy
                ├─── ...

Again using the same generic example folder from above

.. code:: shell

   cd C:\Users\<user name>\Documents\GitHub
   git clone https://github.com/scipy/scipy.git
   cd scipy
   git submodule update --init

Now we need to copy the :code:`openblas.a` file that we've built earlier to the
correct location. If your Python is installed somewhere like the following:

.. code:: shell

   C:\Users\<user name>\AppData\Local\Programs\Python\Python38\python.exe

then you'll need to put the :code:`openblas.a` file that we previously copied
and renamed somewhere like the following:

.. code:: shell

   C:\Users\<user name>\AppData\Local\Programs\Python\Python38\Lib

Adjust the location accordingly based on where :code:`python.exe` is located.

At this stage, we are done with the OpenBLAS part and hopefully we will not need
to build OpenBLAS anytime soon. But we tend to build SciPy more often as it is
on a quicker release cycle. Hence it makes sense to use Windows ``cmd`` or
Powershell for the the build as it is a more native tool. This requires placing
the MinGW compilers on the path.  Hence, make sure that the following
folder (or the folder you have installed MSYS to) is on the system path
variable sufficiently close to the top.

.. code:: shell

    C:\MSYS64\MINGW64\BIN

For a sanity check, restart ``cmd`` or Powershell and type:

.. code:: shell

    gfortran

If you see a missing command error with the above, :code:`gfortran` is not
correctly installed or is still not on the path. However, we assume that it is now
on the path and accessible.

Now install the dependencies that we need to build and test SciPy.

.. code:: shell

    python -m pip install wheel setuptools numpy>=1.16.5 Cython>=0.29.18 pybind11>=2.4.3 pythran>=0.9.12 pytest pytest-xdist

.. note::

    These instructions use ``pip`` as the package manager. You can also use ``conda``
    according to the instructions in the :ref:`conda-guide` with minimal modifications
    (e.g. you don't need to install ``gfortran`` and ``git`` because you already have them).

The last two are for using SciPy's test suite, which is handy if you want to test
some new change locally.

Please note that this is a simpler procedure than what is used for the official
binaries. **Your binaries will only work with the latest NumPy**.
For building against older NumPy versions, see
`Building Against an Older NumPy Version`_.

Assuming that you are in the top of the SciPy repository directory where
``setup.py`` is and assuming that you have set up everything correctly, you
are ready to build. Run the following commands:

.. code:: shell

    python setup.py build

You may verify that the OpenBLAS library was correctly picked up by looking for
the following in your build log:

.. code:: shell

   FOUND:
      libraries = ['openblas']
      library_dirs = ['C:\...........\lib']
      language = f77
      define_macros = [('HAVE_CBLAS', None)]

Notice that there will be multiple lines similar to these. You only need to
track the OpenBLAS one.

When everything finishes without an error, congratulations! You've built
SciPy!

You can further install the built SciPy via

.. code:: shell

    python setup.py install

Just make sure that you uninstalled the existing installation of other SciPy if
there were any (by the regular ``pip uninstall scipy`` machinery).


.. _BLAS: https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms
.. _OpenBLAS: https://github.com/xianyi/OpenBLAS
.. _`msysintaller`: https://www.msys2.org/
.. _`build tools`: https://www.visualstudio.com/downloads/#build-tools-for-visual-studio-2017
.. _`report them`: https://github.com/scipy/scipy/issues/new
.. _`pre-built zip files`: https://3f23b170c54c2533c070-1c8a9b3114517dc5fe17b7c3f8c63a43.ssl.cf2.rackcdn.com/
.. _WindowsCompilers: https://wiki.python.org/moin/WindowsCompilers

Building Against an Older NumPy Version
---------------------------------------

If you want to build SciPy to work with an older NumPy version, then you will need
to replace the NumPy "distutils" folder with the folder from the latest numpy.
The following Powershell snippet can upgrade NumPy distutils while retaining an older
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
and SciPy mailing lists, which can be found `here <https://www.scipy.org/mailing-lists>`__.  There are many
developers out there working on this issue right now, and they would certainly be happy to help you out!  Google is also
a good resource, as there are many people out there who use NumPy and SciPy on Windows, so it would not be surprising if
your question or problem has already been addressed.
