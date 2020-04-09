# Cython Additions

See `PYREADME`.

# HiGHS - Linear optimization software

[![Build Status](https://travis-ci.org/ERGO-Code/HiGHS.svg?branch=master)](https://travis-ci.org/ERGO-Code/HiGHS)

HiGHS is a high performance serial and parallel solver for large scale sparse
linear programming (LP) problems of the form

    Minimize c^Tx subject to L <= Ax <= U; l <= x <= u

It is written in C++ with OpenMP directives, and has been developed and tested on various linux and Windows installations using both the GNU (g++) and Intel (icc) C++ compilers. Note that HiGHS requires (at least) version 4.9 of the GNU compiler. It has no third-party dependencies.

HiGHS is based on the dual revised simplex method implemented in HSOL, which was originally written by Qi Huangfu. Features such as presolve, crash and advanced basis start have been added by Julian Hall, Ivet Galabova. Other features, and interfaces to C, C#, FORTRAN, Julia and Python, have been written by Michael Feldmeier.

Although HiGHS is freely available under the MIT license, we would be pleased to learn about users' experience and give advice via email sent to highsopt@gmail.com.

Reference
---------

Parallelizing the dual revised simplex method
Q. Huangfu and J. A. J. Hall
Mathematical Programming Computation, 10 (1), 119-142, 2018.
DOI: 10.1007/s12532-017-0130-5

http://www.maths.ed.ac.uk/hall/HuHa13/

Performance
-----------

The performance of HiGHS relative to some commercial and open-source simplex solvers may be assessed via the Mittelmann benchmarks on http://plato.asu.edu/ftp/lpsimp.html

Documentation
-------------

The rest of this file gives brief documentation for HiGHS. Comprehensive documentation is available from https://www.highs.dev.

Compilation
-----------

HiGHS uses CMake as build system. First setup
a build folder and call CMake as follows

    mkdir build
    cd build
    cmake ..

Then compile the code using

    make

This installs the executable `bin/highs`.

Testing
-------

To perform a quick test whether the compilation was successful, run

    ctest

Run-time options
----------------

In the following discussion, the name of the executable file generated is
assumed to be `highs`.

HiGHS can read plain text MPS files and LP files and the following command
solves the model in `ml.mps`

    highs ml.mps

HiGHS options
-------------
Usage:
    highs [OPTION...] [file]

      --model_file arg       File of model to solve.
      --presolve arg         Presolve: "choose" by default - "on"/"off" are alternatives.
      --solver arg           Solver: "choose" by default - "simplex"/"ipm" are alternatives.
      --parallel arg         Parallel solve: "choose" by default - "on"/"off" are alternatives.
      --time_limit arg       Run time limit (double).
      --options_file arg     File containing HiGHS options.

  -h, --help                 Print help.

Language interfaces and further documentation
---------------------------------------------

There are HiGHS interfaces for C, C#, FORTRAN, Julia and Python in HiGHS/src/interfaces, with example driver files in HiGHS/examples. Documentation beyond what is in this file is "work in progress", but we expect to have some available before the end of 2019. However, we are happy to give a reasonable level of support via email sent to highsopt@gmail.com.

Parallel code
-------------

At the moment the parallel option is temporarily unavailable due to a large
refactoring in progress. This document will be updated once we have completed
the interface currently being developed.

In order to use OpenMP if available, set`-DOPENMP=ON` during the configuration
step (`cmake ..`).

When compiled with the parallel option on, the number of threads used at run
time is the value of the environment variable `OMP_NUM_THREADS`. For example,
to use HiGHS with eight threads to solve `ml.mps` execute

    export OMP_NUM_THREADS=8
    highs --parallel ml.mps

If `OMP_NUM_THREADS` is not set, either because it has not been set or due to
executing the command

    unset OMP_NUM_THREADS

then all available threads will be used.

If run with `OMP_NUM_THREADS=1`, HiGHS is serial. The `--parallel` run-time
option will cause HiGHS to use serial minor iterations and, although this
could lead to better performance on some problems, performance will typically be
diminished.

When compiled with the parallel option and `OMP_NUM_THREADS>1` or unset, HiGHS
will use multiple threads. If `OMP_NUM_THREADS` is unset, HiGHS will try to use
all available threads so performance may be very slow. Although the best value
will be problem and architecture dependent, `OMP_NUM_THREADS=8` is typically a
good choice. Although HiGHS is slower when run in parallel than in serial for
some problems, it is typically faster in parallel.

HiGHS Library
-------------

HiGHS is compiled in a shared library. Running

`make install`

from the build folder installs the library in `lib/`, as well as all header files in `include/`. For a custom
installation in `install_folder` run

`cmake -DCMAKE_INSTALL_PREFIX=install_folder ..`

and then

`make install`

To use the library from a CMake project use

`find_package(HiGHS)`

and add the correct path to HIGHS_DIR.

Compiling and linking without CMake
-----------------------------------

An executable defined in the file `use_highs.cpp` is linked with the HiGHS library as follows. After running the code above, compile and run with

`g++ -o use_highs use_highs.cpp -I install_folder/include/ -L install_folder/lib/ -lhighs`

`LD_LIBRARY_PATH=install_folder/lib/ ./use_highs`

Interfaces
----------

GAMS
----

Set custom options with `-D<option>=<value>` during the configuration step (`cmake ..`):

- `GAMS_ROOT`:
    path to GAMS system: enables building of GAMS interface

If build with GAMS interface, then HiGHS can be made available as solver
in GAMS by adding an entry for HiGHS to the file gmscmpun.txt in the GAMS
system folder (gmscmpnt.txt on Windows):
```
HIGHS 11 5 0001020304 1 0 2 LP RMIP
gmsgenus.run
gmsgenux.out
/path/to/libhighs.so his 1 1
```

OSI
---
- `OSI_ROOT`:
    path to COIN-OR/Osi build/install folder (OSI_ROOT/lib/pkg-config/osi.pc should exist)
