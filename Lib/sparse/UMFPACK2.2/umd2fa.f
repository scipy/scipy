        SUBROUTINE UMD2FA (N, NE, JOB, TRANSA, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO)
        INTEGER N, NE, JOB, LVALUE, LINDEX, INDEX (LINDEX), KEEP (20),
     $          ICNTL (20), INFO (40)
        DOUBLE PRECISION
     $          VALUE (LVALUE)
        DOUBLE PRECISION
     $          CNTL (10), RINFO (20)
        LOGICAL TRANSA
        
C=== UMD2FA ============================================================
C
C  Unsymmetric-pattern MultiFrontal Package (UMFPACK). Version 2.2d
C  Copyright (C) 1997, Timothy A. Davis, University of Florida, USA.
C  ALL RIGHTS RESERVED.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  July 7, 1997. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C***********************************************************************
C* NOTICE:  "The UMFPACK Package may be used SOLELY for educational,   *
C* research, and benchmarking purposes by non-profit organizations and *
C* the U.S. government.  Commercial and other organizations may make   *
C* use of UMFPACK SOLELY for benchmarking purposes only.  UMFPACK may  *
C* be modified by or on behalf of the User for such use but at no time *
C* shall UMFPACK or any such modified version of UMFPACK become the    *
C* property of the User.  UMFPACK is provided without warranty of any  *
C* kind, either expressed or implied.  Neither the Authors nor their   *
C* employers shall be liable for any direct or consequential loss or   *
C* damage whatsoever arising out of the use or misuse of UMFPACK by    *
C* the User.  UMFPACK must not be sold.  You may make copies of        *
C* UMFPACK, but this NOTICE and the Copyright notice must appear in    *
C* all copies.  Any other use of UMFPACK requires written permission.  *
C* Your use of UMFPACK is an implicit agreement to these conditions."  *
C*                                                                     *
C* The MA38 Package in Release 12 of the Harwell Subroutine Library    *
C* (HSL) has equivalent functionality (and identical calling interface)*
C* as UMFPACK (the HSL has single and double precision versions only,  *
C* however).  It is available for commercial use.   Technical reports, *
C* information on HSL, and matrices are available via the World Wide   *
C* Web at http://www.cis.rl.ac.uk/struct/ARCD/NUM.html, or by          *
C* anonymous ftp at seamus.cc.rl.ac.uk/pub.  Also contact Dr. Scott    *
C* Roberts, Harwell Subroutine Library, B 552, AEA Technology,         *
C* Harwell, Didcot, Oxon OX11 0RA, England.                            *
C* telephone (44) 1235 434988, fax (44) 1235 434136                    *
C* email Scott.Roberts@aeat.co.uk, who will provide details of price   *
C* and conditions of use.                                              *
C***********************************************************************

C=======================================================================
C  HSL Compatibility:  this routine has the same arguments as MA38A/AD. 

C=======================================================================
C  USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Given a sparse matrix A, find a sparsity-preserving and numerically-
C  acceptable pivot order and compute the LU factors, PAQ = LU.  The
C  matrix is optionally preordered into a block upper triangular form
C  (BTF).  Pivoting is performed within each diagonal block to maintain
C  sparsity and numerical stability.  The method used to factorize the
C  matrix is an unsymmetric-pattern variant of the multifrontal method.
C  Most of the floating-point work is done in the Level-3 BLAS (dense
C  matrix multiply).  In addition, approximate degrees are used in the
C  Markowitz-style pivot search to reduce the symbolic overhead.  For
C  best performance, be sure to use an optimized BLAS library.
C
C  This routine is normally preceded by a call to UMD21I to
C  initialize the default control parameters.  UMD21I need only be
C  called once.  A call to UMD2FA can be followed by any number of
C  calls to UMD2SO, which solves a linear system using the LU factors
C  computed by this routine.  A call to UMD2FA can also be followed by
C  any number of calls to UMD2RF, which factorizes another matrix with
C  the same nonzero pattern as the matrix factorized by UMD2FA (but with
C  different numerical values).
C 
C  For more information, see T. A. Davis and I. S. Duff, "An 
C  unsymmetric-pattern multifrontal method for sparse LU factorization",
C  SIAM J. Matrix Analysis and Applications (to appear), also
C  technical report TR-94-038, CISE Dept., Univ. of Florida,
C  P.O. Box 116120, Gainesville, FL 32611-6120, USA.  The method used
C  here is a modification of that method, described in T. A. Davis,
C  "A combined unifrontal/multifrontal method for unsymmetric sparse
C  matrices," TR-94-005.  (Technical reports are available via WWW at
C  http://www.cis.ufl.edu/).  The appoximate degree update algorithm
C  used here has been incorporated into an approximate minimum degree
C  ordering algorithm, desribed in P. Amestoy, T. A. Davis, and I. S.
C  Duff, "An approximate minimum degree ordering algorithm", SIAM J.
C  Matrix Analysis and Applications (to appear, also TR-94-039).  The
C  approximate minimum degree ordering algorithm is implemented as MC47
C  in the Harwell Subroutine Library (MC47 is not called by
C  UMFPACK).

C=======================================================================
C  INSTALLATION NOTE:
C=======================================================================
C
C  Requires the BLAS (Basic Linear Algebra Subprograms) and two routines
C  from the Harwell Subroutine Library.  Ideally, you should use
C  vendor-optimized BLAS for your computer.  If you do not have them,
C  you may obtain the Fortran BLAS from 1.  Send email to 
C  netlib@ornl.gov with the two-line message:
C               send index from blas
C               send blas.shar from blas
C
C  To obtain the two Harwell Subroutine Library (HSL) routines, send
C  email to netlib@ornl.gov with the message:
C               send mc21b.f mc13e.f from harwell
C  These two routines HSL contain additional licensing restrictions.
C  If you want to run UMFPACK without them, see the "INSTALLATION
C  NOTE:" comment in UMD2FB.
C
C  To permamently disable any diagnostic and/or error printing, see
C  the "INSTALLATION NOTE:" comments in UMD2P1 and UMD2P2.
C
C  To change the default control parameters, see the
C  "INSTALLATION NOTE:" comments in UMD21I

C=======================================================================
C  ARGUMENTS:
C=======================================================================

C           ------------------------------------------------------------
C  n:       An integer variable.
C           Must be set by caller on input (not modified).
C           Order of the matrix.  Restriction:  n >= 1.

C           ------------------------------------------------------------
C  ne:      An integer variable.
C           Must be set by caller on input (not modified).
C           Number of entries in input matrix.  Restriction:  ne => 1.

C           ------------------------------------------------------------
C  job:     An integer variable.
C           Must be set by caller on input (not modified).
C           If job=1, then a column-oriented form of the input matrix
C           is preserved, otherwise, the input matrix is overwritten
C           with its LU factors.  If iterative refinement is to done
C           in UMD2SO, (Icntl (8) > 0), then job must be set to 1.

C           ------------------------------------------------------------
C  transa:  A logical variable.
C           Must be set by caller on input (not modified).
C           If false then A is factorized: PAQ = LU.  Otherwise, A
C           transpose is factorized:  PA'Q = LU.

C           ------------------------------------------------------------
C  lvalue:  An integer variable.
C           Must be set by caller on input (not modified).
C           Size of the Value array.  Restriction:  lvalue >= 2*ne
C           is required to convert the input form of the matrix into
C           the internal represenation.  lvalue >= ne + axcopy is
C           required to start the factorization, where axcopy = ne if
C           job = 1, or axcopy = 0 otherwise.  During factorization,
C           additional memory is required to hold the frontal matrices.
C           The internal representation of the matrix is overwritten
C           with the LU factors, of size (Keep (2) - Keep (1) + 1
C           + axcopy), on output.

C           ------------------------------------------------------------
C  lindex:  An integer variable.
C           Must be set by caller on input (not modified).
C           Size of the Index array.  Restriction: lindex >= 3*ne+2*n+1,
C           is required to convert the input form of the matrix into
C           its internal representation.  lindex >= wlen + alen + acopy
C           is required to start the factorization, where
C           wlen <= 11*n + 3*dn + 8 is the size of the workspaces,
C           dn <= n is the number of columns with more than d
C           entries (d = max (64, sqrt (n)) is the default),
C           alen <= 2*ne + 11*n + 11*dn + dne is the size of the
C           internal representation of the matrix, dne <= ne is the
C           number of entries in such columns with more than d entries,
C           and acopy = ne+n+1 if job = 1, or acopy = 0 otherwize.
C           During factorization, the internal representation of size
C           alen is overwritten with the LU factors, of size
C           luilen = (Keep (5) - Keep (3) + 1 - acopy) on output.
C           Additional memory is also required to hold the unsymmetric
C           quotient graph, but this also overwrites the input matrix.
C           Usually about 7*n additional space is adequate for this
C           purpose.  Just prior to the end of factorization,
C           lindex >= wlen + luilen + acopy is required.

C           ------------------------------------------------------------
C  Value:   A double precision array of size lvalue.
C           Must be set by caller on input.  Modified on output.  On
C           input, Value (1..ne) holds the original matrix in triplet
C           form.  On output, Value holds the LU factors, and
C           (optionally) a column-oriented form of the original matrix
C           - otherwise the input matrix is overwritten with the LU
C           factors.

C           ------------------------------------------------------------
C  Index:   An integer array of size lindex.
C           Must be set by caller on input.  Modified on output.  On
C           input, Index (1..2*ne) holds the original matrix in triplet
C           form.  On output, Index holds the LU factors, and
C           (optionally) a column-oriented form of the original matrix
C           - otherwise the input matrix is overwritten with the LU
C           factors.
C
C           On input the kth triplet (for k = 1...ne) is stored as:
C                       A (row,col) = Value (k)
C                       row         = Index (k)
C                       col         = Index (k+ne)
C           If there is more than one entry for a particular position,
C           the values are accumulated, and the number of such duplicate
C           entries is returned in Info (2), and a warning flag is
C           set.  However, applications such as finite element methods
C           naturally generate duplicate entries which are then
C           assembled (added) together.  If this is the case, then
C           ignore the warning message.
C
C           On output, the LU factors and the column-oriented form
C           of A (if preserved) are stored in:
C               Value (Keep (1)...Keep (2))
C               Index (Keep (3)...Keep (5))
C           where Keep (2) = lvalue, and Keep (5) = lindex.

C           ------------------------------------------------------------
C  Keep:    An integer array of size 20.
C
C           Keep (1 ... 5):  Need not be set by caller on input.
C               Modified on output.
C               Keep (1): LU factors start here in Value
C               Keep (2) = lvalue: LU factors end here in Value
C               Keep (3): LU factors start here in Index
C               Keep (4): LU factors needed for UMD2RF start here
C                             in Index
C               Keep (5) = lindex: LU factors end here in Index
C
C           Keep (6 ... 8):  Must be set by caller on input (not
C               modified).
C               integer control arguments not normally modified by the
C               user.  See UMD21I for details, which sets the defaults.
C               Keep (6) is the largest representable positive
C               integer.  Keep (7) and Keep (8) determine the
C               size of d, where columns with more than d original
C               entries are treated as a priori frontal matrices.
C
C           Keep (9 ... 20): Unused.  Reserved for future releases.

C           ------------------------------------------------------------
C  Cntl:    A double precision array of size 10.
C           Must be set by caller on input (not modified).
C           real control arguments, see UMD21I for a description,
C           which sets the defaults. UMD2FA uses Cntl (1) and Cntl (2).

C           ------------------------------------------------------------
C  Icntl:   An integer array of size 20.
C           Must be set by caller on input (not modified).
C           Integer control arguments, see UMD21I for a description,
C           which sets the defaults.  UMD2FA uses Icntl (1..7).

C           ------------------------------------------------------------
C  Info:    An integer array of size 40.
C           Need not be set by caller on input.  Modified on output.
C           It contains information about the execution of UMD2FA.
C
C           Info (1): zero if no error occurred, negative if
C               an error occurred and the factorization was not
C               completed, positive if a warning occurred (the
C               factorization was completed). 
C
C               These errors cause the factorization to terminate:
C
C               Error   Description
C               -1      n < 1
C               -2      ne < 1
C               -3      lindex too small
C               -4      lvalue too small
C               -5      both lindex and lvalue are too small
C
C               With these warnings the factorization was able to 
C               complete:
C
C               Error   Description
C               1       invalid entries
C               2       duplicate entries
C               3       invalid and duplicate entries
C               4       singular matrix
C               5       invalid entries, singular matrix
C               6       duplicate entries, singular matrix
C               7       invalid and duplicate entries, singular matrix
C
C               Subsequent calls to UMD2RF and UMD2SO can only be made
C               if Info (1) is zero or positive.  If Info (1)
C               is negative, then some or all of the remaining
C               Info and Rinfo arrays may not be valid.
C
C           Info (2): duplicate entries in A.  A warning is set
C               if Info (2) > 0.  However, the duplicate entries
C               are summed and the factorization continues.  Duplicate
C               entries are sometimes intentional - for finite element
C               codes, for example.
C
C           Info (3): invalid entries in A, indices not in 1..n.
C               These entries are ignored and a warning is set
C               in Info (1).
C
C           Info (4): zero.  Used by UMD2RF only.
C
C           Info (5): entries in A after adding duplicates and
C               removing invalid entries.
C
C           Info (6): entries in diagonal blocks of A.
C
C           Info (7): entries in off-diagonal blocks of A.  Zero
C               if Info (9) = 1.
C
C           Info (8): 1-by-1 diagonal blocks.
C
C           Info (9): blocks in block-triangular form.
C
C           Info (10): entries below diagonal in L.
C
C           Info (11): entries below diagonal in U.
C
C           Info (12): entries in L+U+offdiagonal part.
C
C           Info (13): frontal matrices.
C
C           Info (14): garbage collections performed on Index, when
C               memory is exhausted.  Garbage collections are performed
C               to remove external fragmentation.  If Info (14) is
C               excessively high, performance can be degraded.  Try
C               increasing lindex if that occurs.  Note that external
C               fragmentation in *both* Index and Value is removed when
C               either is exhausted.
C
C           Info (15): garbage collections performed on Value.
C
C           Info (16): diagonal pivots chosen.
C
C           Info (17): numerically acceptable pivots found in A.
C               If less than n, then A is singular (or nearly so).
C               The factorization still proceeds, and UMD2SO can still
C               be called.  The zero-rank active submatrix of order
C               n - Info (17) is replaced with the identity matrix
C               (assuming BTF is not in use).  If BTF is in use, then
C               one or more of the diagonal blocks are singular. 
C
C           Info (18): memory used in Index.
C
C           Info (19): minimum memory needed in Index
C               (or minimum recommended).  If lindex is set to
C               Info (19) on a subsequent call, then a moderate
C               number of garbage collections (Info (14)) will
C               occur.
C
C           Info (20): memory used in Value.
C
C           Info (21): minimum memory needed in Value
C               (or minimum recommended).  If lvalue is set to
C               Info (21) on a subsequent call, then a moderate
C               number of garbage collections (Info (15)) will
C               occur.
C
C           Info (22): memory needed in Index for the next call to
C               UMD2RF.
C
C           Info (23): memory needed in Value for the next call to
C               UMD2RF.
C
C           Info (24): zero.  Used by UMD2SO only.
C
C           Info (25 ... 40): reserved for future releases

C           ------------------------------------------------------------
C  Rinfo:   A double precision array of size 20.
C           Need not be set by caller on input.  Modified on output.
C           It contains information about the execution of UMD2FA.
C
C           Rinfo (1): total flop count in the BLAS
C
C           Rinfo (2): total assembly flop count
C
C           Rinfo (3): total flops during pivot search
C
C           Rinfo (4): Level-1 BLAS flops
C
C           Rinfo (5): Level-2 BLAS flops
C
C           Rinfo (6): Level-3 BLAS flops
C
C           Rinfo (7): zero.  Used by UMD2SO only.
C
C           Rinfo (8): zero.  Used by UMD2SO only.
C
C           Rinfo (9 ... 20): reserved for future releases

C=======================================================================
C  TO BE PRESERVED BETWEEN CALLS TO UMD2FA, UMD2RF, UMD2SO:
C=======================================================================
C
C  When calling UMD2SO to solve a linear system using the factors
C  computed by UMD2FA or UMD2RF, the following must be preserved:
C
C       n
C       Value (Keep (1)...Keep (2))
C       Index (Keep (3)...Keep (5))
C       Keep (1 ... 20)
C
C  When calling UMD2RF to factorize a subsequent matrix with a pattern
C  similar to that factorized by UMD2FA, the following must be
C  preserved:
C
C       n
C       Index (Keep (4)...Keep (5))
C       Keep (4 ... 20)
C
C  Note that the user may move the LU factors to a different position
C  in Value and/or Index, as long as Keep (1 ... 5) are modified
C  correspondingly.

C## End of user documentation for UMD2FA ###############################

C=======================================================================
C  CODING CONVENTIONS:
C=======================================================================
C
C  This package is written in ANSI Fortran 77.  To make the code more
C  understandable, the following coding conventions are followed for all
C  routines in this package:
C
C  1) Large code blocks are delimited with [...] comments.
C
C  2) GOTO usage:
C       a) Goto's used to return if an error condition is found are
C          written as "GO TO 9000" or "GO TO 9010".
C       b) Goto's used to exit loops prematurely are written as "GO TO",
C          and have a target label of 2000 or less.
C       c) Goto's used to jump to the next iteration of a do loop or
C          while loop (or to implement a while loop) are written as
C          "GOTO".
C       No other goto's are used in this package.
C
C  This package uses the following CRAY compiler directives to help
C  in the vectorization of loops.  Each of them operate on the
C  do-loop immediately following the directive.  Other compilers
C  normally treat these directives as ordinary comments.
C
C       CFPP$ NODEPCHK L        disables data dependency check, and
C                               asserts that no recursion exists.
C       CFPP$ NOLSTVAL L        disables the saving of last values of
C                               transformed scalars (indexes or promoted
C                               scalars, especially those in array
C                               subscripts).  Asserts that values do not
C                               need to be the same as in the scalar
C                               version (for later use of the scalars).
C       CDIR$ SHORTLOOP         asserts that the loop count is always
C                               64 or less.

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   user routine
C       subroutines called:     UMD2ER, UMD2P1, UMD2CO, UMD2F0
C       functions called:       MAX, MIN
        INTRINSIC MAX, MIN

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER I, NZ, LUX1, LUI1, IUSE, XUSE, LUIR1, NZOFF, NBLKS,
     $          MAXINT, NMAX
        LOGICAL PRESRV
        DOUBLE PRECISION
     $          IGNORE

C  Location of LU factors:
C  -----------------------
C  lux1:    real part of LU factors placed in Value (lux1 ... lvalue)
C  lui1:    integer part of LU factors placed in Index (lui1 ... lindex)
C  luir1:   Index (luir1 ... lindex) must be preserved for UMD2RF
C
C  Memory usage:
C  -------------
C  iuse:    current memory usage in Index
C  xuse:    current memory usage in Value
C
C  Matrix to factorize:
C  --------------------
C  nblks:   number of diagonal blocks (1 if BTF not used)
C  nzoff:   entries in off-diagonal part (0 if BTF not used)
C  nz:      entries in matrix after removing invalid/duplicate entries
C
C  Other:
C  ------
C  maxint:  largest representable positive integer
C  nmax:    largest permissible value of n
C  i:       general loop index
C  presrv:  true if original matrix to be preserved

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

C-----------------------------------------------------------------------
C  clear informational output, and Keep array (except Keep (6..8)):
C-----------------------------------------------------------------------

        DO 10 I = 1, 40 
           INFO (I) = 0
10      CONTINUE 
        DO 20 I = 1, 20 
           RINFO (I) = 0
20      CONTINUE 
        KEEP (1) = 0
        KEEP (2) = 0
        KEEP (3) = 0
        KEEP (4) = 0
        KEEP (5) = 0
        IGNORE = 0

C-----------------------------------------------------------------------
C  print input arguments if requested
C-----------------------------------------------------------------------

        CALL UMD2P1 (1, 1,
     $          N, NE, JOB, TRANSA, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO,
     $          IGNORE, IGNORE, 1, IGNORE, 1)

C-----------------------------------------------------------------------
C  initialize and check inputs
C-----------------------------------------------------------------------

        IUSE = 0
        XUSE = 0
        INFO (5) = NE
        INFO (6) = NE
        MAXINT = KEEP (6)
        NMAX = (MAXINT - 2) / 3
        IF (N .LT. 1) THEN 
C          n is too small
           CALL UMD2ER (1, ICNTL, INFO, -1, -1)
           GO TO 9000
        ENDIF 
        IF (NE .LT. 1) THEN 
C          ne is too small
           CALL UMD2ER (1, ICNTL, INFO, -2, -1)
           GO TO 9000
        ENDIF 

C-----------------------------------------------------------------------
C  get memory for conversion to column form
C-----------------------------------------------------------------------

        NZ = NE
        IUSE = 2*N+1 + MAX (2*NZ, N+1) + NZ
        XUSE = 2*NZ
        INFO (18) = IUSE
        INFO (20) = XUSE
        INFO (19) = IUSE
        INFO (21) = XUSE
        IF (LINDEX .LT. IUSE) THEN 
C          set error flag if out of integer memory:
           CALL UMD2ER (1, ICNTL, INFO, -3, IUSE)
        ENDIF 
        IF (LVALUE .LT. XUSE) THEN 
C          set error flag if out of real memory:
           CALL UMD2ER (1, ICNTL, INFO, -4, XUSE)
        ENDIF 
        IF (INFO (1) .LT. 0) THEN 
C          error return, if not enough integer and/or real memory:
           GO TO 9000
        ENDIF 

C-----------------------------------------------------------------------
C  convert to column-oriented form and remove duplicates
C-----------------------------------------------------------------------

        CALL UMD2CO (N, NZ, TRANSA, VALUE, LVALUE, INFO, ICNTL,
     $     INDEX, LINDEX-(2*N+1), INDEX(LINDEX-2*N), INDEX(LINDEX-N), 1)
        IF (INFO (1) .LT. 0) THEN 
C          error return, if all entries invalid (nz is now 0):
           GO TO 9000
        ENDIF 

C-----------------------------------------------------------------------
C  current memory usage:
C-----------------------------------------------------------------------

C       Index (1..n+1): column pointers.  input matrix is now in
C       Index (1..nz+n+1) and Value (1..nz)
C       col pattern: Index (n+1+ Index (col) ... n+1+ Index (col+1))
C       col values:  Value (     Index (col) ...      Index (col+1))
C       at this point, nz <= ne (nz = ne if there are no invalid or
C       duplicate entries; nz < ne otherwise).

        IUSE = NZ + (N+1)
        XUSE = NZ

C-----------------------------------------------------------------------
C  factorize
C-----------------------------------------------------------------------

        PRESRV = JOB .EQ. 1
        IF (PRESRV) THEN 

C          -------------------------------------------------------------
C          keep a copy of the original matrix in column-oriented form
C          -------------------------------------------------------------

C          copy column pointers (Cp (1..n+1) = Ap (1..n+1))
           IUSE = IUSE + (N+1)
CFPP$ NODEPCHK L
           DO 30 I = 1, N+1 
              INDEX (NZ+N+1+I) = INDEX (I)
30         CONTINUE 

           CALL UMD2F0 (N, NZ, INDEX (NZ+N+2),
     $          VALUE (NZ+1), LVALUE-NZ,
     $          INDEX (NZ+2*N+3), LINDEX-(NZ+2*N+2),
     $          LUX1, LUI1, IUSE, XUSE, NZOFF, NBLKS,
     $          ICNTL, CNTL, INFO, RINFO,
     $          PRESRV, INDEX, INDEX (N+2), VALUE, N, NZ, KEEP, NE)
           IF (INFO (1) .LT. 0) THEN 
C             error return, if UMD2F0 fails
              GO TO 9000
           ENDIF 
C          adjust pointers to reflect Index/Value, not II/XX:
           LUX1 = LUX1 + NZ
           LUI1 = LUI1 + (NZ+2*N+2)

C          move preserved copy of A to permanent place
           LUX1 = LUX1 - NZ
           LUI1 = LUI1 - (NZ+N+1)
           DO 40 I = NZ+N+1, 1, -1 
              INDEX (LUI1+I-1) = INDEX (I)
40         CONTINUE 
           DO 50 I = NZ, 1, -1 
              VALUE (LUX1+I-1) = VALUE (I)
50         CONTINUE 

        ELSE 

C          -------------------------------------------------------------
C          do not preserve the original matrix
C          -------------------------------------------------------------

           CALL UMD2F0 (N, NZ, INDEX,
     $          VALUE, LVALUE,
     $          INDEX (N+2), LINDEX-(N+1),
     $          LUX1, LUI1, IUSE, XUSE, NZOFF, NBLKS,
     $          ICNTL, CNTL, INFO, RINFO,
     $          PRESRV, 1, 1, IGNORE, 0, 1, KEEP, NE)
           IF (INFO (1) .LT. 0) THEN 
C             error return, if UMD2F0 fails
              GO TO 9000
           ENDIF 
C          adjust pointers to reflect Index/Value, not II/XX:
           LUI1 = LUI1 + (N+1)
        ENDIF 

C-----------------------------------------------------------------------
C  wrap-up
C-----------------------------------------------------------------------

        IF (TRANSA) THEN 
           INDEX (LINDEX-6) = 1
        ELSE 
           INDEX (LINDEX-6) = 0
        ENDIF 

        INDEX (LINDEX-5) = NZOFF
        INDEX (LINDEX-4) = NBLKS
        IF (PRESRV) THEN 
           INDEX (LINDEX-3) = 1
        ELSE 
           INDEX (LINDEX-3) = 0
        ENDIF 
        INDEX (LINDEX-2) = NZ
        INDEX (LINDEX-1) = N
        INDEX (LINDEX) = NE

C       do not need preserved matrix (n+1+nz), or off-diagonal entries
C       (nzoff) for UMD2RF:
        LUIR1 = LUI1
        IF (PRESRV) THEN 
C          do not need preserved matrix for UMD2RF
           LUIR1 = LUIR1 + N+1 + NZ
        ENDIF 
        IF (NBLKS .GT. 1) THEN 
C          do not need off-diagonal part for UMD2RF
           LUIR1 = LUIR1 + NZOFF
        ENDIF 

C       save location of LU factors
        KEEP (1) = LUX1
        KEEP (2) = LVALUE
        KEEP (3) = LUI1
        KEEP (4) = LUIR1
        KEEP (5) = LINDEX

C       update memory usage information
        IUSE = LINDEX - LUI1 + 1
        XUSE = LVALUE - LUX1 + 1
        INFO (22) = INFO (22) + (LINDEX - LUIR1 + 1)

C-----------------------------------------------------------------------
C  print the output arguments if requested, and return
C-----------------------------------------------------------------------

C       error return label:
9000    CONTINUE
        IF (INFO (1) .LT. 0) THEN 
           KEEP (1) = 0
           KEEP (2) = 0
           KEEP (3) = 0
           KEEP (4) = 0
           KEEP (5) = 0
        ENDIF 

        INFO (18) = MIN (LINDEX, MAX (INFO (18), IUSE))
        INFO (20) = MIN (LVALUE, MAX (INFO (20), XUSE))

        CALL UMD2P1 (1, 2,
     $          N, NE, JOB, TRANSA, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO,
     $          IGNORE, IGNORE, 1, IGNORE, 1)
        RETURN
        END 
