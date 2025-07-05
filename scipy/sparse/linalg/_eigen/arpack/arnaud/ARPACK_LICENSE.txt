This is the ARPACK package from
http://www.caam.rice.edu/software/ARPACK/

Specifically the files are from
http://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz
with the patch
http://www.caam.rice.edu/software/ARPACK/SRC/patch.tar.gz

The ARPACK README is at
http://www.caam.rice.edu/software/ARPACK/SRC/readme.arpack

---

ARPACK is a collection of Fortran77 subroutines designed to solve large 
scale eigenvalue problems. 

The package is designed to compute a few eigenvalues and corresponding 
eigenvectors of a general n by n matrix A. It is most appropriate for large 
sparse or structured matrices A where structured means that a matrix-vector
product w <- Av requires order n rather than the usual order n**2 floating 
point operations. This software is based upon an algorithmic variant of the
Arnoldi process called the Implicitly Restarted Arnoldi Method (IRAM). When
the matrix A is symmetric it reduces to a variant of the Lanczos process 
called the Implicitly Restarted Lanczos Method (IRLM). These variants may be
viewed as a synthesis of the Arnoldi/Lanczos process with the Implicitly 
Shifted QR technique that is suitable for large scale problems. For many 
standard problems, a matrix factorization is not required. Only the action
of the matrix on a vector is needed.  ARPACK software is capable of solving
large scale symmetric, nonsymmetric, and generalized eigenproblems from 
significant application areas. The software is designed to compute a few (k)
eigenvalues with user specified features such as those of largest real part 
or largest magnitude.  Storage requirements are on the order of n*k locations.
No auxiliary storage is required. A set of Schur basis vectors for the desired
k-dimensional eigen-space is computed which is numerically orthogonal to working
precision. Numerically accurate eigenvectors are available on request. 

Important Features: 

    o  Reverse Communication Interface. 
    o  Single and Double Precision Real Arithmetic Versions for Symmetric,
       Non-symmetric, Standard or Generalized Problems.
    o  Single and Double Precision Complex Arithmetic Versions for Standard
       or Generalized Problems. 
    o  Routines for Banded Matrices - Standard or Generalized Problems. 
    o  Routines for The Singular Value Decomposition. 
    o  Example driver routines that may be used as templates to implement
       numerous Shift-Invert strategies for all problem types, data types
       and precision. 

---

The ARPACK license is the BSD 3-clause license ("New BSD License")
http://www.caam.rice.edu/software/ARPACK/RiceBSD.txt

---

BSD Software License

Pertains to ARPACK and P_ARPACK

Copyright (c) 1996-2008 Rice University.
Developed by D.C. Sorensen, R.B. Lehoucq, C. Yang, and K. Maschhoff.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.

- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
