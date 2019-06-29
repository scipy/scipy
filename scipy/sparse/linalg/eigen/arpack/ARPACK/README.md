ARPACK-NG is a collection of Fortran77 subroutines designed to solve large scale
eigenvalue problems.

<a href="https://travis-ci.org/opencollab/arpack-ng"><img src="https://travis-ci.org/opencollab/arpack-ng.svg"/></a><br/>
[![Coverage Status](https://coveralls.io/repos/github/opencollab/arpack-ng/badge.svg?branch=master)](https://coveralls.io/github/opencollab/arpack-ng?branch=master)

Important Features:

* Reverse Communication Interface.
* Single and Double Precision Real Arithmetic Versions for Symmetric,
  Non-symmetric, Standard or Generalized Problems.
* Single and Double Precision Complex Arithmetic Versions for Standard or
  Generalized Problems.
* Routines for Banded Matrices - Standard or Generalized Problems.
* Routines for The Singular Value Decomposition.
* Example driver routines that may be used as templates to implement numerous
  Shift-Invert strategies for all problem types, data types and precision.
* arpackmm: utility to test arpack with matrix market files.
  Note: to run this utility, you need the eigen library (to handle RCI).
* ILP64 support:
  * users: set INTERFACE64 at configure time.
  * developers:
    * all files which needs ILP64 support must include "arpackdef.h".
    * when coding, use a_int (defined in arpackdef.h) instead of int.
      a_int stands for "architecture int": it's #defined to int or int64_t according
      to the architecture.
  * example: to test arpack with sequential ILP64 MKL assuming you use gnu compilers
    ```$ ./bootstrap
    $ export FFLAGS='-I/usr/include/mkl'
    $ export FCFLAGS='-I/usr/include/mkl'
    $ export LIBS='-Wl,--no-as-needed -lmkl_sequential -lmkl_core -lpthread -lm -ldl'
    $ export INTERFACE64=1
    $ ./configure --with-blas=mkl_gf_ilp64 --with-lapack=mkl_gf_ilp64
    $ make all check```

This project started as a joint project between Debian, Octave and Scilab in order to
provide a common and maintained version of arpack.
This is now a community project maintained by a few volunteers.

Indeed, no single release has been published by Rice university for the last
few years and since many software (Octave, Scilab, R, Matlab...) forked it and
implemented their own modifications, arpack-ng aims to tackle this by providing
a common repository, maintained versions with a testsuite.

arpack-ng is replacing arpack almost everywhere.

1. You have successfully unbundled ARPACK-NG and are now in the ARPACK-NG
   directory that was created for you.

2. The directory SRC contains the top level routines including
   the highest level reverse communication interface routines

* ssaupd, dsaupd - symmetric single and double precision
* snaupd, dnaupd - non-symmetric single and double precision
* cnaupd, znaupd - complex non-symmetric single and double precision

   The headers of these routines contain full documentation of calling
   sequence and usage.  Additional information is in the DOCUMENTS directory.

   The directory PARPACK contains the Parallel ARPACK routines.


3. Example driver programs that illustrate all the computational modes,
   data types and precisions may be found in the EXAMPLES directory.
   Upon executing the 'ls EXAMPLES' command you should see

* BAND
* COMPLEX
* NONSYM
* README
* SIMPLE
* SVD
* SYM

   Example programs for banded, complex, nonsymmetric, symmetric,
   and singular value decomposition may be found in the directories
   BAND, COMPLEX, NONSYM, SYM, SVD respectively.  Look at the README
   file for further information.  To get started, get into the SIMPLE
   directory to see example programs that illustrate the use of ARPACK in
   the simplest modes of operation for the most commonly posed
   standard eigenvalue problems.


   Example programs for Parallel ARPACK may be found in the directory
   PARPACK/EXAMPLES. Look at the README file for further information.

   The following instructions explain how to make the ARPACK library.

4. Unlike ARPACK, ARPACK-NG is providing autotools and cmake based build
   system and iso_c_binding support (which enables to call fortran 
   subroutines natively from C or C++).

Therefore, the classical commands should work as expected:
   
   
    $ sh bootstrap
    $ ./configure
    $ make
    $ make check
    $ make install

Furthermore, ARPACK-NG now provides CMake functionality:
   
    $ mkdir build
    $ cd build
    $ cmake -D EXAMPLES=ON -D MPI=ON -D BUILD_SHARED_LIBS=ON ..
    $ make
    $ make install
   builds everything including examples and parallel support (with MPI).

To build with code coverage:
   
    $ mkdir build
    $ cd build
    $ cmake -DCOVERALLS=ON -DCMAKE_BUILD_TYPE=Debug ..
    $ make all check test coveralls

To get iso_c_binding support:
   
    $ ./configure --enable-icb
    $ cmake -D ICB=ON
   The install will now provide arpack.h/hpp, parpack.h/hpp and friends.
   Examples of use can be found in ./TESTS and ./PARPACK/TESTS/MPI.
   A few related links can be found here:

   * http://fortranwiki.org/fortran/show/ISO_C_BINDING
   * http://fortranwiki.org/fortran/show/Generating+C+Interfaces
   * https://www.roguewave.com/sites/rw/files/attachments/StandardizedMixedLanguageProgrammingforCandFortran.pdf

5. Within DOCUMENTS directory there are three files

* ex-sym.doc
* ex-nonsym.doc and
* ex-complex.doc

   for templates on how to invoke the computational modes of ARPACK.
   Also look in the README.MD file for explanations concerning the
   other documents.
 
 Good luck and enjoy.
