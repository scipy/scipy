        SUBROUTINE UMD2RF (N, NE, JOB, TRANSA, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO)
        INTEGER N, NE, JOB, LVALUE, LINDEX, INDEX (LINDEX), KEEP (20),
     $          ICNTL (20), INFO (40)
        DOUBLE PRECISION
     $          VALUE (LVALUE)
        DOUBLE PRECISION
     $          CNTL (10), RINFO (20)
        LOGICAL TRANSA
        
C=== UMD2RF ============================================================
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
C  HSL Compatibility:  this routine has the same arguments as MA38B/BD. 

C=======================================================================
C  USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Given a sparse matrix A, and a sparsity-preserving and numerically-
C  acceptable pivot order and symbolic factorization, compute the LU
C  factors, PAQ = LU.  Uses the sparsity pattern and permutations from
C  a prior factorization by UMD2FA or UMD2RF.  The matrix A should have
C  the same nonzero pattern as the matrix factorized by UMD2FA or
C  UMD2RF.  The matrix can have different numerical values.  No
C  variations are made in the pivot order computed by UMD2FA.  If a
C  zero pivot is encountered, an error flag is set and the
C  factorization terminates.
C
C  This routine can actually handle any matrix A such that (PAQ)_ij can
C  be nonzero only if (LU)_ij is be nonzero, where L and U are the LU
C  factors of the matrix factorized by UMD2FA.  If BTF (block triangular
C  form) is used, entries above the diagonal blocks of (PAQ)_ij can have
C  an arbitrary sparsity pattern.  Entries for which (LU)_ij is not
C  present, or those below the diagonal blocks are invalid and ignored
C  (a warning flag is set and the factorization proceeds without the
C  invalid entries).  A listing of the invalid entries can be printed.
C
C  This routine must be preceded by a call to UMD2FA or UMD2RF.
C  A call to UMD2RF can be followed by any number of calls to UMD2SO,
C  which solves a linear system using the LU factors computed by this
C  routine or by UMD2FA.  A call to UMD2RF can also be followed by any
C  number of calls to UMD2RF.

C=======================================================================
C  ARGUMENTS:
C=======================================================================

C           ------------------------------------------------------------
C  n:       An integer variable.
C           Must be set by caller on input (not modified).
C           Order of the matrix.  Must be identical to the value of n
C           in the last call to UMD2FA.

C           ------------------------------------------------------------
C  ne:      An integer variable.
C           Must be set by caller on input (not modified).
C           Number of entries in input matrix.  Normally not modified
C           since the last call to UMD2FA.
C           Restriction:  1 <= ne < (Keep (4)) / 2

C           ------------------------------------------------------------
C  job:     An integer variable.
C           Must be set by caller on input (not modified).
C           If job=1, then a column-oriented form of the input matrix
C           is preserved, otherwise, the input matrix is overwritten
C           with its LU factors.  If iterative refinement is to done
C           (Icntl (8) > 0), then job must be set to 1.  Can be
C           the same, or different, as the last call to UMD2FA.

C           ------------------------------------------------------------
C  transa:  A logical variable.
C           Must be set by caller on input (not modified).
C           If false then A is factorized: PAQ = LU.  Otherwise, A
C           transpose is factorized:  PA'Q = LU.  Normally the same as
C           the last call to UMD2FA.

C           ------------------------------------------------------------
C  lvalue:  An integer variable.
C           Must be set by caller on input (not modified).
C           Size of the Value array.  Restriction:  lvalue >= 2*ne,
C           although a larger will typically be required to complete
C           the factorization.  The exact value required is computed
C           by the last call to UMD2FA or UMD2RF (Info (23)).
C           This value assumes that the ne, job, and transa parameters
C           are the same as the last call.  Some garbage collection may
C           occur if lvalue is set to Info (23), but usually not
C           much.  We recommend lvalue => 1.2 * Info (23).  The
C           lvalue parameter is usually the same as in the last call to
C           UMD2FA, however.

C           ------------------------------------------------------------
C  lindex:  An integer variable.
C           Must be set by caller on input (not modified).
C           Size of the Index array.  Restriction:
C           lindex >= 3*ne+2*n+1 + (Keep (5) - Keep (4) + 1),
C           although a larger will typically be required to complete
C           the factorization.  The exact value required is computed
C           by the last call to UMD2FA or UMD2RF (Info (22)).
C           This value assumes that the ne, job, and transa parameters
C           are the same as the last call.  No garbage collection ever
C           occurs in the Index array, since UMD2RF does not create
C           external fragmentation in Index.  The lindex parameter is
C           usually the same as in the last call to UMD2FA, however.
C           Note that lindex >= Keep (5) is also required, since
C           the pattern of the prior LU factors reside in
C           Index (Keep (4) ... Keep (5)).

C           ------------------------------------------------------------
C  Value:   A double precision array of size lvalue.
C           Must be set by caller on input (normally from the last call
C           to UMD2FA or UMD2RF).  Modified on output.  On input,
C           Value (1..ne) holds the original matrix in triplet form.
C           On output, Value holds the LU factors, and (optionally) a
C           column-oriented form of the original matrix - otherwise
C           the input matrix is overwritten with the LU factors.

C           ------------------------------------------------------------
C  Index:   An integer array of size lindex.
C           Must be set by caller on input (normally from the last call
C           to UMD2FA or UMD2RF).  Modified on output.  On input,
C           Index (1..2*ne) holds the original matrix in triplet form,
C           and Index (Keep (4) ... Keep (5)) holds the pattern
C           of the prior LU factors.  On output, Index holds the LU
C           factors, and (optionally) a column-oriented form of the
C           original matrix - otherwise the input matrix is overwritten
C           with the LU factors.
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
C           On input, and the pattern of the prior LU factors is in
C               Index (Keep (4) ... Keep (5))
C
C           On output, the LU factors and the column-oriented form
C           of A (if preserved) are stored in:
C               Value (Keep (1)...Keep (2))
C               Index (Keep (3)...Keep (5))
C           where Keep (2) = lvalue, and Keep (5) = lindex.

C           ------------------------------------------------------------
C  Keep:    An integer array of size 20.
C
C           Keep (1 ... 3):  Need not be set by caller on input.
C               Modified on output.
C               Keep (1): new LU factors start here in Value
C               Keep (2) = lvalue: new LU factors end here in Value
C               Keep (3): new LU factors start here in Index
C
C           Keep (4 ... 5): Must be set by caller on input (normally
C               from the last call to UMD2FA or UMD2RF). Modified on
C               output.
C               Keep (4):  On input, the prior LU factors start here
C               in Index, not including the prior (optionally) preserved
C               input matrix, nor the off-diagonal pattern (if BTF was
C               used in the last call to UMD2FA).  On output, the new
C               LU factors needed for UMD2RF start here in Index.
C               Keep (5):  On input, the prior LU factors end here in
C               Index.  On output, Keep (5) is set to lindex, which
C               is where the new LU factors end in Index
C
C           Keep (6 ... 8):  Unused.  These are used by UMD2FA only.
C               Future releases may make use of them, however.
C
C           Keep (9 ... 20): Unused.  Reserved for future releases.

C           ------------------------------------------------------------
C  Cntl:    A double precision array of size 10.
C           Must be set by caller on input (not modified).
C           Control arguments, see UMD21I for a
C           description, which sets the default values.  The current
C           version of UMD2RF does not actually use Cntl.  It is
C           included to make the argument list of UMD2RF the same as
C           UMD2FA.  UMD2RF may use Cntl in future releases.

C           ------------------------------------------------------------
C  Icntl:   An integer array of size 20.
C           Must be set by caller on input (not modified).
C           Integer control arguments, see UMD21I for a description,
C           which sets the default values.  UMD2RF uses Icntl (1),
C           Icntl (2), Icntl (3), and Icntl (7).

C           ------------------------------------------------------------
C  Info:    An integer array of size 40.
C           Need not be set by caller on input.  Modified on output.
C           It contains information about the execution of UMD2RF.
C
C           Info (1): zero if no error occurred, negative if
C               an error occurred and the factorization was not
C               completed, positive if a warning occurred (the
C               factorization was completed).
C
C               These errors cause the factorization to terminate:
C
C               Error   Description
C               -1      n < 1 or n > maximum value
C               -2      ne < 1 or ne > maximum value
C               -3      lindex too small
C               -4      lvalue too small
C               -5      both lindex and lvalue are too small
C               -6      prior pivot ordering no longer acceptable
C               -7      LU factors are uncomputed, or are corrupted
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
C               These entries are ignored and a warning is set in
C               Info (1).
C
C           Info (4): invalid entries in A, not in prior LU
C               factors.  These entries are ignored and a warning is
C               set in Info (1).
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
C           Info (14): zero.  Used by UMD2FA only.
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
C               UMD2RF can be called if the value of Info (17) 
C               returned by UMD2FA was less than n, but the order
C               (n - Info (17)) active submatrix is still replaced
C               with the identity matrix.  Entries residing in this
C               submatrix are ignored, their number is included in
C               Info (4), and a warning is set in Info (1).
C
C           Info (18): memory used in Index.
C
C           Info (19): memory needed in Index (same as Info (18)).
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
C           It contains information about the execution of UMD2RF.
C
C           Rinfo (1): total flop count in the BLAS
C
C           Rinfo (2): total assembly flop count
C
C           Rinfo (3): zero.  Used by UMD2FA only.
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

C## End of user documentation for UMD2RF ###############################

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   user routine
C       subroutines called:     UMD2ER, UMD2P1, UMD2CO, UMD2R0
C       functions called:       MAX, MIN
        INTRINSIC MAX, MIN

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER I, NZ, LUX1, LUI1, IUSE, XUSE, N1, NZ1, NBLKS,
     $          LIND2, LUIR1, LUSIZ, LUI2, RPERMP, CPERMP,
     $          OFFPP, LUBLPP, BLKPP, ON, NZOFF, IP2, IO, PRL
        LOGICAL PRESRV, BADLU
        DOUBLE PRECISION
     $          IGNORE

C  Printing control:
C  -----------------
C  io:      I/O unit for diagnostic messages
C  prl:     printing level
C
C  Matrix to factorize:
C  --------------------
C  nz:      number of entries, after removing invalid/duplicate entries
C  presrv:  true if original matrix to be preserved
C
C  Memory usage:
C  -------------
C  iuse:    current memory usage in Index
C  xuse:    current memory usage in Value
C  lind2:   allocatable part of Index is (1..lind2)
C
C  Location and status of LU factors:
C  ----------------------------------
C  lui1:    integer part of LU factors start in Index (lui1...)
C  luir1:   Index (luir1 ... lui2) is needed for this call to UMD2RF
C  lusiz:   size of Index (luir1..lui2), needed from prior LU factors
C  lui2:    integer part of LU factors end in Index (..lui2)
C  lux1:    real part of LU factors in Value (lux1...lvalue)
C  ip2:     pointer into trailing part of LU factors in Index
C  badlu:   if true, then LU factors are corrupted or not computed
C
C  Arrays and scalars allocated in LU factors (in order):
C  ------------------------------------------------------
C  ...      LU factors of each diagonal block located here
C  lublpp:  LUblkp (1..nblks) array in Index (lublpp..lublpp+nblks-1)
C  blkpp:   Blkp (1..nblks+1) array loc. in Index (blkpp...blkpp+nblks)
C  offpp:   Offp (1..n+1) array located in Index (offpp...offpp+n)
C  on:      size of Offp array
C  cpermp:  Cperm (1..n) array located in Index (cpermp...cpermp+n-1)
C  rpermp:  Rperm (1..n) array located in Index (rpermp...rpermp+n-1)
C  nblks:   number of diagonal blocks
C  nz1:     number of entries when prior matrix factorize
C  n1:      N argument in UMD2FA or UMD2RF when prior matrix factorized 
C
C  Other:
C  ------
C  i:       loop index

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

        IGNORE = 0
        IO = ICNTL (2)
        PRL = ICNTL (3)

C-----------------------------------------------------------------------
C  clear informational output, and the unneeded part of the Keep array
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

C-----------------------------------------------------------------------
C  print input arguments if requested
C-----------------------------------------------------------------------

        CALL UMD2P1 (2, 1,
     $          N, NE, JOB, TRANSA, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO,
     $          IGNORE, IGNORE, 1, IGNORE, 1)

C-----------------------------------------------------------------------
C  check input arguments
C-----------------------------------------------------------------------

        IUSE = 0
        XUSE = 0
        INFO (5) = NE
        INFO (6) = NE
        IF (N .LT. 1) THEN 
C          n is too small
           CALL UMD2ER (2, ICNTL, INFO, -1, -1)
           GO TO 9000
        ENDIF 
        IF (NE .LT. 1) THEN 
C          ne is too small
           CALL UMD2ER (2, ICNTL, INFO, -2, -1)
           GO TO 9000
        ENDIF 

C-----------------------------------------------------------------------
C  get pointers to integer part of prior LU factors
C-----------------------------------------------------------------------

        LUIR1 = KEEP (4)
        LUI2 = KEEP (5)
        LUSIZ = LUI2 - LUIR1 + 1

        BADLU = LUIR1 .LE. 0 .OR. LUI2-6 .LT. LUIR1 .OR. LUI2.GT.LINDEX
        IF (BADLU) THEN 
           CALL UMD2ER (2, ICNTL, INFO, -7, 0)
C          error return, LU factors are corrupted:
           GO TO 9000
        ENDIF 
        IF (2*NE .GT. LUIR1) THEN 
           CALL UMD2ER (2, ICNTL, INFO, -2, LUIR1/2)
C          error return, ne is too large:
           GO TO 9000
        ENDIF 

C-----------------------------------------------------------------------
C  shift the prior LU factors down to the end of Index.  If Keep and
C  lindex are unmodified from the prior call to UMD2FA, then
C  Keep (5) = lindex, and this shift is not performed.
C-----------------------------------------------------------------------

        IF (LUI2 .LT. LINDEX) THEN 
           DO 30 I = LINDEX, LINDEX - LUSIZ + 1, -1 
              INDEX (I) = INDEX (I - LINDEX + LUI2)
30         CONTINUE 
           LUIR1 = LINDEX - LUSIZ + 1
           KEEP (5) = LINDEX
           KEEP (4) = LUIR1
        ENDIF 

C-----------------------------------------------------------------------
C  get seven scalars (transa, nzoff, nblks, presrv, nz, n, ne) from LU
C-----------------------------------------------------------------------

C       ne1 = Index (lindex), not required for UMD2RF
        N1 = INDEX (LINDEX-1)
        NZ1 = INDEX (LINDEX-2)
C       presr1 = Index (lindex-3) .ne. 0, not required for UMD2RF
        NBLKS = INDEX (LINDEX-4)
C       nzoff1 = Index (lindex-5), not required for UMD2RF
C       trans1 = Index (lindex-6) .ne. 0, not required for UMD2RF

C-----------------------------------------------------------------------
C  get pointers to permutation vectors
C-----------------------------------------------------------------------

        RPERMP = (LINDEX-6) - N
        CPERMP = RPERMP - N
        IP2 = CPERMP - 1

C-----------------------------------------------------------------------
C  get pointers to block-triangular information, if BTF was used
C-----------------------------------------------------------------------

        IF (NBLKS .GT. 1) THEN 

C          -------------------------------------------------------------
C          get pointers to BTF arrays
C          -------------------------------------------------------------

           OFFPP = CPERMP - (N+1)
           BLKPP = OFFPP - (NBLKS+1)
           LUBLPP = BLKPP - (NBLKS)
           IP2 = LUBLPP - 1
           ON = N

        ELSE 

C          -------------------------------------------------------------
C          matrix was factorized as a single block, pass dummy arg.
C          -------------------------------------------------------------

           OFFPP = 1
           BLKPP = 1
           LUBLPP = 1
           ON = 0

        ENDIF 

        BADLU = N .NE. N1 .OR. NZ1 .LE. 0 .OR. LUIR1 .GT. IP2 .OR.
     $          NBLKS .LE. 0 .OR. NBLKS .GT. N
        IF (BADLU) THEN 
           CALL UMD2ER (2, ICNTL, INFO, -7, 0)
C          error return, LU factors are corrupted:
           GO TO 9000
        ENDIF 

C-----------------------------------------------------------------------
C  get memory for conversion to column form
C-----------------------------------------------------------------------

        NZ = NE
        IUSE = 2*N+1 + MAX (2*NZ,N+1) + NZ + LUSIZ
        XUSE = 2*NZ
        INFO (18) = IUSE
        INFO (20) = XUSE
        INFO (19) = IUSE
        INFO (21) = XUSE
        INFO (23) = XUSE
        LIND2 = LUIR1 - 1
        IF (LINDEX .LT. IUSE) THEN 
C          set error flag if out of integer memory
           CALL UMD2ER (2, ICNTL, INFO, -3, IUSE)
        ENDIF 
        IF (LVALUE .LT. XUSE) THEN 
C          set error flag if out of real memory
           CALL UMD2ER (2, ICNTL, INFO, -4, XUSE)
        ENDIF 
        IF (INFO (1) .LT. 0) THEN 
C          error return, if not enough integer and/or real memory
           GO TO 9000
        ENDIF 

C-----------------------------------------------------------------------
C  convert to column-oriented form and remove duplicates
C-----------------------------------------------------------------------

        CALL UMD2CO (N, NZ, TRANSA, VALUE, LVALUE, INFO, ICNTL,
     $     INDEX, LIND2-(2*N+1), INDEX (LIND2-2*N), INDEX (LIND2-N), 2)
        IF (INFO (1) .LT. 0) THEN 
C          error return, if all entries are invalid (nz is now 0)
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
C       Pattern of prior LU factors and BTF arrays are in
C       Index (Keep (4) ... Keep (5))

        IUSE = NZ + (N+1) + LUSIZ
        XUSE = NZ

C-----------------------------------------------------------------------
C  refactorize
C-----------------------------------------------------------------------

        PRESRV = JOB .EQ. 1
        IF (PRESRV) THEN 

C          -------------------------------------------------------------
C          keep a copy of the original matrix in column-oriented form
C          -------------------------------------------------------------

C          copy column pointers (Cp (1..n+1) = Ap (1..n+1))
           IUSE = IUSE + (N+1)
CFPP$ NODEPCHK L
           DO 40 I = 1, N+1 
              INDEX (NZ+N+1+I) = INDEX (I)
40         CONTINUE 

           CALL UMD2R0 (N, NZ, INDEX (NZ+N+2),
     $          VALUE (NZ+1), LVALUE-NZ,
     $          INDEX (NZ+2*N+3), LIND2-(NZ+2*N+2),
     $          LUX1, LUI1, IUSE, XUSE, NZOFF, NBLKS,
     $          ICNTL, CNTL, INFO, RINFO,
     $          PRESRV, INDEX, INDEX (N+2), VALUE, N, NZ,
     $          INDEX (LUIR1), IP2 - LUIR1 + 1,
     $          INDEX (LUBLPP), INDEX (BLKPP), INDEX (OFFPP), ON,
     $          INDEX (CPERMP), INDEX (RPERMP), NE)
           IF (INFO (1) .LT. 0) THEN 
C             error return, if UMD2R0 fails
              GO TO 9000
           ENDIF 
C          adjust pointers to reflect Index/Value, not II/XX:
           LUX1 = LUX1 + NZ
           LUI1 = LUI1 + (NZ+2*N+2)

C          move preserved copy of A to permanent place
           LUX1 = LUX1 - (NZ)
           LUI1 = LUI1 - (NZ+N+1)
           DO 50 I = NZ+N+1, 1, -1 
              INDEX (LUI1+I-1) = INDEX (I)
50         CONTINUE 
           DO 60 I = NZ, 1, -1 
              VALUE (LUX1+I-1) = VALUE (I)
60         CONTINUE 

        ELSE 

C          -------------------------------------------------------------
C          do not preserve the original matrix
C          -------------------------------------------------------------

           CALL UMD2R0 (N, NZ, INDEX,
     $          VALUE, LVALUE,
     $          INDEX (N+2), LIND2-(N+1),
     $          LUX1, LUI1, IUSE, XUSE, NZOFF, NBLKS,
     $          ICNTL, CNTL, INFO, RINFO,
     $          PRESRV, 1, 1, IGNORE, 0, 1,
     $          INDEX (LUIR1), IP2 - LUIR1 + 1,
     $          INDEX (LUBLPP), INDEX (BLKPP), INDEX (OFFPP), ON,
     $          INDEX (CPERMP), INDEX (RPERMP), NE)
           IF (INFO (1) .LT. 0) THEN 
C             error return, if UMD2R0 fails
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

C       save location of LU factors
        KEEP (1) = LUX1
        KEEP (2) = LVALUE
        KEEP (3) = LUI1
        KEEP (4) = LUIR1
        KEEP (5) = LINDEX

C       update memory usage information
        IUSE = LINDEX - LUI1 + 1
        XUSE = LVALUE - LUX1 + 1

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
        INFO (19) = INFO (18)
        INFO (22) = INFO (19)
        INFO (20) = MIN (LVALUE, MAX (INFO (20), XUSE))

        CALL UMD2P1 (2, 2,
     $          N, NE, JOB, TRANSA, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO,
     $          IGNORE, IGNORE, 1, IGNORE, 1)
        RETURN
        END 
