# arpack-ng [![arpack-ng CI/CD](https://github.com/opencollab/arpack-ng/actions/workflows/jobs.yml/badge.svg)](https://github.com/opencollab/arpack-ng/actions/workflows/jobs.yml)

ARPACK-NG is a collection of Fortran77 subroutines designed to solve large scale eigenvalue problems.
| mandatory dependencies | optional dependencies     | category      |
|------------------------|---------------------------|---------------|
| BLAS, LAPACK           | MPI, Eigen3, Boost.Python | LinearAlgebra |

## About the project

This project started as a joint project between Debian, Octave and Scilab in order to provide a common and maintained version of arpack.
This is now a community project maintained by a few volunteers.
Indeed, no single release has been published by Rice university for the last few years and since many software (Octave, Scilab, R, Matlab...)
forked it and implemented their own modifications, arpack-ng aims to tackle this by providing a common repository, maintained versions with a testsuite.
`arpack-ng` is replacing arpack almost everywhere.

## Important Features

- Reverse Communication Interface (RCI).
- Single and Double Precision Real Arithmetic Versions for Symmetric, Non-symmetric, Standard or Generalized Problems.
- Single and Double Precision Complex Arithmetic Versions for Standard or Generalized Problems.
- Routines for Banded Matrices - Standard or Generalized Problems.
- Routines for The Singular Value Decomposition.
- Example driver routines that may be used as templates to implement numerous
- Shift-Invert strategies for all problem types, data types and precision.
- `arpackmm`: utility to test arpack with matrix market files. Note: to run this utility, you need the eigen library (to handle RCI).

## Documentation

Within DOCUMENTS directory there are three files for templates on how to invoke the computational modes of ARPACK.

- ex-sym.doc
- ex-nonsym.doc and
- ex-complex.doc

Also look in the README.MD file for explanations concerning the
other documents.

## ILP64 support

About ILP64 support:

- Sequential arpack supports [ILP64](https://www.intel.com/content/www/us/en/develop/documentation/onemkl-linux-developer-guide/top/linking-your-application-with-onemkl/linking-in-detail/linking-with-interface-libraries/using-the-ilp64-interface-vs-lp64-interface.html), but, parallel arpack doesn't.
- Reminder: you can NOT mix `ILP64` with `LP64`. If you compile `arpack-ng` with `ILP64` (resp. `LP64`) support, you MUST insure your BLAS/LAPACK is compliant with `ILP64` (resp. `LP64`).
- Set `INTERFACE64` at configure time.

Note for F77/F90 developers:

- All files which needs `ILP64` support must include `"arpackicb.h"`.
- When coding, use `i_int` (defined in `arpackicb.h`) instead of `c_int`. `i_int` stands for ISO_C_BINDING integer: it's `#defined` to `c_int` or `c_int64_t` according to the architecture.

Note for C/C++ developers:

- All files which needs `ILP64` support must include `"arpackdef.h"`.
- When coding, use `a_int` (defined in `arpackdef.h`) instead of `int`. Here, `a_int` stands for "architecture int": it's `#defined` to `int` or `int64_t` according to the architecture.

**Example**: to test arpack with sequential `ILP64` MKL assuming you use gnu compilers

```bash
$ ./bootstrap
$ export FFLAGS='-DMKL_ILP64 -I/usr/include/mkl'
$ export FCFLAGS='-DMKL_ILP64 -I/usr/include/mkl'
$ export LIBS='-Wl,--no-as-needed -L/usr/lib/x86_64-linux-gnu -lmkl_sequential -lmkl_core -lpthread -lm -ldl'
$ export INTERFACE64=1
$ ./configure --with-blas=mkl_gf_ilp64 --with-lapack=mkl_gf_ilp64
$ make all check
```

## ISO_C_BINDING support

About ISO_C_BINDING support:

- The install will now provide `arpack.h/hpp`, `parpack.h/hpp` and friends.
- Examples of use can be found in `./TESTS` and` ./PARPACK/TESTS/MPI`.

ISO_C_BINDING is a feature of modern Fortran meant to handle safely interoperability between Fortran and C (in practice, no more need to use ugly tricks to link F77 functions to C code using "underscored" symbols). Basically, ISO_C_BINDING make sure all fortran variables are typed (which may not always be the case when using `implicit` keyword in fortran): this way, C compilers can link properly. For more informations on ISO_C_BINDING, you can checkout the following links:

- <http://fortranwiki.org/fortran/show/ISO_C_BINDING>
- <http://fortranwiki.org/fortran/show/Generating+C+Interfaces>

Using ICB is seamless:

- Compile `arpack-ng` with ISO_C_BINDING: you'll get both old-fashion fortran symbols and new ISO_C_BINDING symbols available for linking.
- Add `#include "arpack.h"` in your C code.
- Replace all [sdcz][ae]upd calls by [sdcz][ae]upd_c: functions suffixed with _c are ISO_C_BINDING compliant (exposing same arguments than original fortran functions).

**Example**: to test arpack with ISO_C_BINDING

```bash
$ ./configure --enable-icb
$ cmake -D ICB=ON
```

## Eigen support

`arpack-ng` provides C++ eigensolver based on both ISO_C_BINDING and `eigen`.

Check out `./EXAMPLES/MATRIX_MARKET/README` for more details.

**Example**: to test arpack with `eigen`

```bash
$ mkdir build
$ cd build
$ cmake -D EXAMPLES=ON -D ICB=ON -D EIGEN=ON ..
$ make all check
```

## Python support

`pyarpack`: python support based on `Boost.Python.Numpy` exposing C++ API.
`pyarpack` exposes in python the `arpack-ng` C++ eigensolver (based on `eigen`).

Check out `./EXAMPLES/PYARPACK/README` for more details.

**Example**: to test arpack with python3

```bash
$ mkdir build
$ cd build
$ cmake -D EXAMPLES=ON -D ICB=ON -D EIGEN=ON -D PYTHON3=ON ..
$ make all check
```

## 📁 Directory structure

- You have successfully unbundled ARPACK-NG` and are now in the ARPACK-NG directory that was created for you.

- The directory SRC contains the top level routines including the highest level **reverse communication interface** routines

  - `ssaupd`, `dsaupd`: symmetric single and double precision
  - `snaupd`, `dnaupd`: non-symmetric single and double precision
  - `cnaupd`, `znaupd`: complex non-symmetric single and double precision
  - The headers of these routines contain full documentation of calling sequence and usage.
  - Additional information is given in the `/DOCUMENTS` directory.

- The directory `PARPACK` contains the Parallel ARPACK routines.

- Example driver programs that illustrate all the computational modes, data types and precisions may be found in the EXAMPLES directory. Upon executing the `ls EXAMPLES` command you should see the following directories

  ```bash
  ├── BAND
  ├── COMPLEX
  ├── Makefile.am
  ├── MATRIX_MARKET
  ├── NONSYM
  ├── PYARPACK
  ├── README
  ├── SIMPLE
  ├── SVD
  └── SYM
  ```

  - Example programs for banded, complex, nonsymmetric, symmetric, and singular value decomposition may be found in the directories BAND, COMPLEX, NONSYM, SYM, SVD respectively.
  - Look at the README file for further information.
  - To get started, get into the SIMPLE directory to see example programs that illustrate the use of ARPACK in the simplest modes of operation for the most commonly posed standard eigenvalue problems.

> Example programs for Parallel ARPACK may be found in the directory `PARPACK/EXAMPLES`. Look at the README file for further information.

## Install 🚀

### Getting arpack-ng

Unlike ARPACK, ARPACK-NG is providing autotools and cmake based build system. In addition, `ARPACK-NG` also provides
ISO_C_BINDING support, which enables to call fortran subroutines natively from C or C++.

First, obtain the source code 📥 from github:

```bash
$ git clone https://github.com/opencollab/arpack-ng.git
$ cd ./arpack-ng
```

If you prefer the ssh to obtain the source code, then use:

```bash
$ git clone git@github.com:opencollab/arpack-ng.git
$ cd ./arpack-ng
```

> Note, It is recommended to install `arpack` at standard location on your system by using your root privilege.

### Using autotools

In the source directory, use the following commands to configure, build and install `arpack-ng`.

```bash
$ sh bootstrap
$ ./configure --enable-mpi
$ make
$ make check
$ sudo make install
```

Congratulations 🎉, you have installed `arpack` lib using autotools (caution: you need `sudo` to install in your system).

The above-mentioned process will build everything including the examples and parallel support using MPI.

### Using cmake

You can install `ARPACK-NG` by using cmake. If you do not have cmake, then please download the binary from `pip` using:

```bash
$ python3 -m pip install cmake
$ which cmake && cmake --version
```

After installing cmake, follow the instruction given below.

Caution: Make sure you are in source directory of ARPACK-NG.

```bash
$ mkdir build
$ cd build
$ cmake -D EXAMPLES=ON -D MPI=ON -D BUILD_SHARED_LIBS=ON ..
$ make
$ sudo make install
```

✨ Congratulations, you have installed `arpack` lib using cmake (caution: you need `sudo` to install in your system).

The above-mentioned process will build everything including the examples and parallel support using MPI.

### Customize build / install

You can also customize the installation of `arpack` using the autotools.

To customize the install directories:

```bash
$ LIBSUFFIX="64" ./configure
$ make all install
```

To enable ILP64 support:

```bash
$ INTERFACE64="1" ITF64SUFFIX="ILP64" ./configure
$ make all install
```

To enable ISO_C_BINDING support:

```bash
$ ./configure --enable-icb
```

You can customize the build by declaring the cmake options during configuration.

To customize the install directories:

```bash
$ cmake -D LIBSUFFIX="64" ..
$ make all install
```

To enable ILP64 support:

```bash
$ cmake -D INTERFACE64=ON -D ITF64SUFFIX="ILP64" ..
$ make all install
```

To enable ISO_C_BINDING support:

```bash
$ cmake -D ICB=ON
```

## Supported Operating Systems:

### Linux support

`arpack-ng` runs on debian-based distros.

### Mac OS support

On mac OS, with GNU compilers, you may need to customize options:

```bash
$ LIBS="-framework Accelerate" FFLAGS="-ff2c -fno-second-underscore" FCFLAGS="-ff2c -fno-second-underscore" ./configure
```

### Windows support

`arpack-ng` can be installed on Windows as a MinGW-w64 package via various distribution, for example through [MSYS2](https://packages.msys2.org/package/mingw-w64-x86_64-arpack) with `pacman -S mingw-w64-x86_64-arpack`. It can also be built and installed through [vcpkg](https://github.com/microsoft/vcpkg) with `vcpkg install arpack-ng`.

## Using arpack-ng from your own codebase

The `*.pc` and `*.cmake` files provided by `arpack-ng` are only pointing to arpack libraries.
If you need other libraries (like MPI), you must add them alongside arpack (see CMake example below).

Typically, if you need

- ARPACK: at compile/link time, you'll need to provide BLAS and LAPACK.

- ARPACK with eigen support (arpackSolver): at compile/link time, you'll need to provide BLAS, LAPACK and Eigen.

- PARPACK: at compile/link time, you'll need to provide BLAS, LAPACK and MPI.

Examples are provided in `tstCMakeInstall.sh` and `tstAutotoolsInstall.sh` generated after running cmake/configure.

### With autotools

First, set `PKG_CONFIG_PATH` to the location in the installation directory where `arpack.pc` lies.

Then, insert the following lines in your `configure.ac`:
```
PKG_CHECK_MODULES([ARPACK], [arpack])
AC_SUBST([ARPACK_CFLAGS])
AC_SUBST([ARPACK_LIBS])
```

Note: make sure you have installed `pkg-config`.

### With CMake

You can use arpack in your CMake builds by using `ARPACK::ARPACK` target. For example,

```cmake
FIND_PACKAGE(arpackng)
ADD_EXECUTABLE(main main.f)
TARGET_INCLUDE_DIRECTORIES(main PUBLIC ARPACK::ARPACK)
TARGET_LINK_LIBRARIES(main ARPACK::ARPACK)
```

To use PARPACK in your Cmake builds, use `PARPACK::PARPACK` target:

```cmake
FIND_PACKAGE(arpackng)
FIND_PACKAGE(MPI REQUIRED COMPONENTS Fortran)
ADD_EXECUTABLE(main main.f)
TARGET_INCLUDE_DIRECTORIES(main PUBLIC PARPACK::PARPACK)
TARGET_LINK_LIBRARIES(main PARPACK::PARPACK)
TARGET_INCLUDE_DIRECTORIES(main PUBLIC MPI::MPI_Fortran)
TARGET_LINK_LIBRARIES(main MPI::MPI_Fortran)
```

Note: Make sure to update `CMAKE_MODULE_PATH` env variable (otheriwse, `find_package` won't find arpack-ng cmake file).

### FAQ

- Where can I find ARPACK user's guide?

  http://li.mit.edu/Archive/Activities/Archive/CourseWork/Ju_Li/MITCourses/18.335/Doc/ARPACK/Lehoucq97.pdf

- Calling arpack's aupd methods returns `info = -9 - Starting vector is zero.`: why?

  Residuals are null. Try to set `resid` to small values (like epsilon machine magnitude) but *not exactly* zero.
  Residuals `resid = A*v - lamdba*v` target *exactly* the zero vector.
  When `resid` is close enough to zero, the iterative procedure stops.

- Say I have an estimate of an eigen value, how to give this information to arpack?

  You need to shift of an amount of about this estimate of `lambda`. Grep `backTransform` in `arpackSolver.hpp` to see an example.
  For more informations, checkout "NUMERICAL METHODS FOR LARGE EIGENVALUE PROBLEMS" by Yousef Saad: https://www-users.cse.umn.edu/~saad/eig_book_2ndEd.pdf (paragraph 4.1.2. and section 4.1.).

- Say I have an estimate of an eigen vector, how to give this information to arpack?

  You need to copy this eigen vector estimate in `v` (not `resid`) and set `info` to 1 before calling aupd methods.
  The `v` vector targets a non-null vector such that `resid = 0`, that is, such that `A*v = lambda*v`.

- Using PARPACK, I get incorrect eigen values.

  Make sure each MPI processor handles a subpart of the eigen system (matrices) only.
  ARPACK handles and solves the whole eigen problem (matrices) at once.
  PARPACK doesn't: each MPI processor must handle and solve a subpart of the eigen system (matrices) only (independently from the other processors).
  See examples for Fortran in folder `PARPACK/EXAMPLES/MPI`, and for C/C++ examples in `PARPACK/TESTS/MPI/icb_parpack_c.c` and `PARPACK/TESTS/MPI/icb_parpack_cpp.cpp`

## Using MKL instead of BLAS / LAPACK

How to use arpack-ng with Intel MKL:

- Let autotools/cmake find MKL for you based on pkg-config files (setting `PKG_CONFIG_PATH`) or cmake options (`BLA_VENDOR=Intel10_64lp` for lp64, `BLA_VENDOR=Intel10_64ilp` for ilp64).
- Refers to the Intel Link Advisor: <https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html>.

## Good luck and enjoy 🎊
