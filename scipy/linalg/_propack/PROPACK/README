PROPACK Version 2.1,                              Stanford, April 2005

OVERVIEW
  This directory contains a Fortran version of the PROPACK software,
which is designed to efficiently compute the singular values and
singular vectors of a large, sparse and/or structured matrix. The
basic Krylov-subspace algorithm used is Lanczos bidiagonalization,
implemented with partial reorthogonalization. The use of partial
reorthogonalization often improves performance significantly compared
to the classic Lanczos algorithm with full reorthogonalization; the
exact amount of improvement depends on the distribution of the
singular values. Two sets of SVD routines are available, on with and
one without implicit restarting. Implicit restarting allows the
computation of a given number of singular values and corresponding
vectors to be done in a fixed amount of memory. The amount of memory
used by the ordinary version is proportional to the number of
iterations required for the singular values to converge, and this is
generally not known in advance, but since the total number of
matrix-vector multiplications needed is usually lower for the
non-restarted version it still can be the method of choice in many
cases.

  The main driver routines DLANSVD and DLANSVD_IRL are found in
"dlansvd.F" and "dlansvd_irl.F", which also contain descriptions of
the input parameters. A set of example programs for computing the SVD
of sparse matrices in several simple formats, including the commonly
used Harwell-Boeing format, are included in the Examples directory.


INSTALLATION

  To install the software follow the steps below:

1. Uncompress and untar the files using 

% gunzip PROPACK77.tar.gz
% tar xvf PROPACK77.tar

2. Edit the make option file make.<plat> in the PROPACK/Make directory
where <plat> corresponds to your platform. Currently make option files
for <plat> = { linux_gcc_ia32 | linux_icc_ia32 | linux_gcc_ia64 |
linux_icc_ia64 | irix | sunos | ibm } are available. In particular you
need to set the variables LINKFLAGS, LINKPATH and BLAS such that the
BLAS library installed on your system is linked correctly (see
below). You can also set various flags passed to the compiler and
linker. After you have done this type

% ./configure

in the PROPACK directory. The configure script determines the platform
you are running on and generates 'make.inc' with all the platform
dependent flags based on the appropriate make.<plat> from the Make
directory. On Intel based platforms (ia32 and ia64) the configure
script takes the optional argument "-icc", which will select the make
configuration in make.linux_icc_ia32 and make.linux_icc_ia64, which
uses the Intel icc and ifc/ifort compilers. If available, the Intel
compilers usually generate significantly faster code than gcc, in
particular for the ia64 platform. On AIX, ia32, ia64, and IRIX
platforms the option "-openmp", passed to the configure script, will
cause a multi-threaded (parallel) version of PROPACK to be built. The
parallelization is done using the OpenMP shared memory programming
model (see http://www.openmp.org/), and the number of threads
(processors) used can be selected by setting the environment variable
OMP_NUM_THREADS to the desired number before running a
program. Warning: The parallelization is very fine grained and thus
mostly suited for large matrices (m,n > 100,000, say) or possibly
smaller matrices when running on (non-distributed) shared memory
computers with low memory latency. The parallel performance on
machines with distributed memory leaves something to be desired
(is very far from linear speedup).

3. Build the libraries by typing

% make

This will build the libraries lib<precision>propack_<PLAT>.a, which
contains the PROPACK routines proper, and
lib<precision>lapack_util_<PLAT>.a, which contains various LAPACK 3.0
routines called by PROPACK. Here <PLAT> refers to the platform name
specified in make.inc, and <precision> is "s", "d", "c", and "z",
corresponding to single (real*4), double precision (real*8), complex
(complex*8) and double complex (complex*16). To use the PROPACK
routines, link your program with lib<precision>propack_<PLAT>.a,
lib<precision>lapack_util_<PLAT>.a and the BLAS library on your
system. The libraries corresponding to the four different precisions
are located in the directories single, double, complex8, and
complex16.


EXAMPLE PROGRAMS

  Two example programs "example.F" and "example_irl.F" are provided
for each of the four precisions in the subdirectory Examples.
"example.F" illustrated how to compute part of the SVD using the
non-restarted algorithm, while "example_irl.F" illustrates the use of
the implicitly restarted version.Build and run them by typing

% cd <precision>/Examples
% make 
% example.<PLAT>.x < example.in
% example_irl.<PLAT>.x < example_irl.in

The example programs read a matrix stored in Harwell-Boeing format
from a file and compute a number of singular values as specified in
the input file. A test matrix from the Harwell-Boeing collection is
provided in the file Examples/illc1850.rra (for single and double) and
Examples/mhd1280b.cua (for complex8 and complex16). For more test
matrices see, e.g., the Matrix Market website:

  http://math.nist.gov/MatrixMarket.  

The example programs can also read matrices stored in diagonal,
coordinate or dense formats (binary or ASCII), which is useful for
testing the algorithms with known test matrices without having to
write new code. See Examples/example.F and Examples/matvec.F for
details. Examples of real matrices stored in coordinate and diagonal
ASCII format are provided in Examples/illc1850.coord and
Examples/illc1850.diag.

WARNING: Matrices stored in binary format are often incompatible 
between machines with different wordsize, e.g. 32-bit vesus 64-bit,
or endianess, i.e. little-endian (x86/Itanium/Alpha) versus 
big-endian (PPC/Power/MIPS).


TESTING THE INSTALLATION

The output produced by the example programs, compiled with the GCC
3.2.2 compiler on a Linux workstation with a Pentium 4 processor, is
provided in the files

  Sigma_200_illc1850.ascii, U_200_illc1850.ascii, 
  V_200_illc1850.ascii

and

  Sigma_IRL_200_illc1850.ascii, U_IRL_200_illc1850.ascii, 
  V_IRL_200_illc1850.ascii,

which are located in the directories <precision>/Examples/Output. 
Typing

% make; make test; make verify

in the top-level PROPACK directory will build the example programs for
all precisions, run them with the provided test matrices and verify
that the results are consistent with those in the files listed above
using the program in Examples/compare.F. The comparison is mainly
meant to catch serious bugs or errors in the installation, so the
error bounds used in the test are quite generous. Small-ish
differences caused by different round-off errors or sloppy floating
point arithmetic on some platforms should not generate any warnings.
For the test examples in double and complex*16 precision the maximal
relative error in the singular values should be of the order
1e-15. For the test examples in single and complex*8 precision the
maximal relative error in the singular values should be of the order
1e-6.


OBTAINING THE BLAS LIBRARY
  If your system does not already have this library installed, we
recommend using the freely available and very fast version by
Kazushige Goto (UT-Austin and the Japan Patent Office), which can be
downloaded here: http://www.cs.utexas.edu/users/flame/goto. Another
set of fast BLAS routines optimized for various platforms is available
from the ATLAS project at the Netlib software repository, see
http://www.netlib.org/atlas. More information about the BLAS as well
as generic (un-optimized) Fortran source code is available at
http://www.netlib.org/blas.


CONTACT INFORMATION
  Questions and comments about PROPACK are welcome and should be
directed to:

Rasmus Munk Larsen      
W.W. Hansen Experimental Physics Laboratory (HEPL), Annex A210
Stanford University,  Stanford, CA 94305-4085
E-mail: rmunk@quake.stanford.edu 

(C) Rasmus Munk Larsen, Stanford University, March 2004.
