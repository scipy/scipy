        SUBROUTINE UMC2F0 (N, NZ, CP, XX, XSIZE, II, ISIZE, XTAIL,
     $          ITAIL, IUSE, XUSE, NZOFF, NBLKS, ICNTL, CNTL, INFO,
     $          RINFO, PRESRV, AP, AI, AX, AN, ANZ, KEEP, NE)
        INTEGER N, NZ, ISIZE, II (ISIZE), ICNTL (20), INFO (40),
     $          CP (N+1), XSIZE, XTAIL, ITAIL, IUSE, XUSE, AN, ANZ,
     $          AP (AN+1), AI (ANZ), KEEP (20), NZOFF, NBLKS, NE
        LOGICAL PRESRV
        COMPLEX
     $          XX (XSIZE), AX (ANZ)
        REAL
     $          CNTL (10), RINFO (20)
        
C=== UMC2F0 ============================================================
C
C  Unsymmetric-pattern MultiFrontal Package (UMFPACK). Version 2.2c
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
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Factorize an unsymmetric sparse matrix in column-form, optionally
C  permuting the matrix to upper block triangular form and factorizing
C  each diagonal block.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       n:              order of matrix
C       nz:             entries in matrix, after removing duplicates
C                       and invalid entries.
C       ne:             number of triplets, unchanged from UMC2FA
C       Cp (1..n+1):    column pointers of input matrix
C       presrv:         if true then preserve original matrix
C       xsize:          size of XX
C       isize:          size of II
C       iuse:           memory usage in Index on input
C       xuse:           memory usage in Value on input
C       Icntl:          integer control parameters, see UMC21I
C       Cntl:           real control parameters, see UMC21I
C       Keep (6..8):    integer control parameters, see UMC21I
C
C       if presrv is true:
C           an:                 = n, order of preserved matrix
C           anz:                = anz, order of preserved matrix
C           Ap (1..an+1):       column pointers of preserved matrix
C           Ai (1..nz):         row indices of preserved matrix
C           Ax (1..nz):         values of preserved matrix
C           II:                 unused on input
C           XX:                 unused on input
C       else
C           an:                 1
C           anz:                1
C           Ap:                 unused
C           Ai:                 unused
C           Ax:                 unused
C           II (1..nz):         row indices of input matrix
C           XX (1..nz):         values of input matrix

C=======================================================================
C  OUTPUT: 
C=======================================================================
C
C       XX (xtail ... xsize), xtail:
C
C                       LU factors are located in XX (xtail ... xsize),
C                       including values in off-diagonal part if matrix
C                       was permuted to block triangular form.
C
C       II (itail ... isize), itail:
C
C                       LU factors are located in II (itail ... isize),
C                       including pattern, row and column permutations,
C                       block triangular information, etc.  See umf2fa
C                       for more information.
C
C       Info:           integer informational output, see UMC2FA
C       Rinfo:          real informational output, see UMC2FA
C
C       iuse:           memory usage in Index on output
C       xuse:           memory usage in Value on output
C
C       nzoff:          entries in off-diagonal part (0 if BTF not used)
C       nblks:          number of diagonal blocks (1 if BTF not used)

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   UMC2FA
C       subroutines called:     UMC2ER, UMC2FB, UMC2F1, UMC2OF
C       functions called:       MAX
        INTRINSIC MAX

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER KN, NZDIA, BLKPP, LUBLPP, P, OFFIP, XHEAD, ROW,
     $          OFFXP, OFFPP, IHEAD, K1, K2, BLK, PRP, P2, CPERMP,
     $          RPERMP, NSGLTN, NPIV, MNZ, NSYM, K, COL, RMAX, CMAX,
     $          TOTNLU, XRMAX, XRUSE
        LOGICAL TRYBTF, IOUT, XOUT
        COMPLEX
     $          A

C  Allocated array pointers:
C  -------------------------
C  blkpp:   Blkp (1..nblks+1) array located in II (blkpp..blkp+nblks)
C  lublpp:  LUblkp (1..nblks) array loc. in II (lublpp..lublpp+nblks-1)
C  offip:   Offi (1..nzoff) array located in II (offip..offip+nzoff-1)
C  offxp:   Offx (1..nzoff) array located in XX (offxp..offxp+nzoff-1)
C  offpp:   Offp (1..n+1) array located in II (offpp..offpp+n)
C  cpermp:  Cperm (1..n) array located in II (cpermp..cpermp+n-1)
C  rpermp:  Rperm (1..n) array located in II (rpermp..rpermp+n-1)
C  prp:     Pr (1..n) work array located in II (prp..prp+n-1)
C
C  BTF information:
C  ----------------
C  k1:      starting index of diagonal block being factorized
C  k2:      ending index of diagonal block being factorized
C  kn:      the order of the diagonal block being factorized
C  blk:     block number of diagonal block being factorized
C  trybtf:  true if BTF is to be attempted (= Icntl (4) .eq. 1)
C  nzdia:   number of entries in diagonal blocks (= nz if BTF not used)
C  nsgltn:  number of 1-by-1 diagonal blocks ("singletons")
C  npiv:    number of numerically valid singletons
C  a:       numerical value of a singleton
C  mnz:     nzoff
C
C  Memory usage:
C  -------------
C  xhead:   XX (1..xhead-1) is in use, XX (xhead..xtail-1) is free
C  ihead:   II (1..ihead-1) is in use, II (ihead..itail-1) is free
C  iout:    true if UMC2F1 ran out of integer memory, but did not
C           set error flag
C  xout:    true if UMC2F2 ran out of integer memory, but did not
C           set error flag
C
C  Estimated memory for UMC2RF:
C  ----------------------------
C  rmax:    largest contribution block is cmax-by-rmax
C  cmax:       "         "         "    "   "   "  "
C  totnlu:  total number of LU arrowheads in all diagonal blocks
C  xrmax:   estimated maximum real memory usage for UMC2RF
C  xruse:   estimated current real memory usage for UMC2RF
C
C  Other:
C  ------
C  k:       loop index (kth pivot)
C  row:     row index
C  col:     column index
C  p:       pointer
C  p2:      pointer
C  nsym:    number of symmetric pivots chosen

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

C-----------------------------------------------------------------------
C  get input parameters and initialize
C-----------------------------------------------------------------------

        NBLKS = 1
        NZOFF = 0
        NZDIA = NZ
        NSGLTN = 0
        NPIV = 0
        RMAX = 1
        CMAX = 1
        TOTNLU = 0
        IF (PRESRV) THEN 
C          original matrix is not in Cp/II/XX, but in Ap/Ai/Ax:
           IHEAD = 1
           XHEAD = 1
        ELSE 
           IHEAD = NZ + 1
           XHEAD = NZ + 1
        ENDIF 
        ITAIL = ISIZE + 1
        XTAIL = XSIZE + 1

C-----------------------------------------------------------------------
C  allocate permutation arrays: Cperm (1..n) and Rperm (1..n), and
C  seven scalars:  transa, nzoff, nblks, presrv, nz, n, ne
C  (in that order) at tail of II (in LU factors)
C-----------------------------------------------------------------------

        ITAIL = ITAIL - (2*N+7)
        IUSE = IUSE + (2*N+7)
        INFO (18) = MAX (INFO (18), IUSE)
        INFO (19) = INFO (18)
        CPERMP = ITAIL
        RPERMP = CPERMP + N
        IF (IHEAD .GT. ITAIL) THEN 
C          error return, if not enough integer memory:
           GO TO 9000
        ENDIF 

C-----------------------------------------------------------------------
C  Find permutations to block upper triangular form, if requested.
C-----------------------------------------------------------------------

        TRYBTF = ICNTL (4) .EQ. 1
        IF (TRYBTF) THEN 

C          -------------------------------------------------------------
C          get workspace at tail of II of size 6n+2
C          -------------------------------------------------------------

           ITAIL = ITAIL - (N+1)
           OFFPP = ITAIL
           ITAIL = ITAIL - (5*N+1)
           P = ITAIL
           IUSE = IUSE + (6*N+2)
           INFO (18) = MAX (INFO (18), IUSE)
           INFO (19) = INFO (18)

C          -------------------------------------------------------------
           IF (PRESRV) THEN 
C          find permutation, but do not convert to BTF form
C          -------------------------------------------------------------

              IF (IHEAD .GT. ITAIL) THEN 
C                error return, if not enough integer memory:
                 GO TO 9000
              ENDIF 
              CALL UMC2FB (AX, ANZ, AI, ANZ, N, NZ, NZDIA, NZOFF,
     $           NBLKS, CP, II (CPERMP), II (RPERMP), II(P), II(P+N),
     $           II (P+2*N), II (P+3*N), II (P+4*N), II (OFFPP),
     $           PRESRV, ICNTL)

C          -------------------------------------------------------------
           ELSE 
C          find permutation, convert to BTF form, and discard original
C          -------------------------------------------------------------

C             use additional size nz temporary workspace in II and XX
              IHEAD = IHEAD + NZ
              XHEAD = XHEAD + NZ
              IUSE = IUSE + NZ
              XUSE = XUSE + NZ
              INFO (18) = MAX (INFO (18), IUSE)
              INFO (20) = MAX (INFO (20), XUSE)
              INFO (19) = INFO (18)
              INFO (21) = INFO (20)
              IF (IHEAD .GT. ITAIL .OR. XHEAD .GT. XTAIL) THEN 
C                error return, if not enough integer and/or real memory:
                 GO TO 9000
              ENDIF 
              CALL UMC2FB (XX, 2*NZ, II, 2*NZ, N, NZ, NZDIA, NZOFF,
     $              NBLKS, CP, II (CPERMP), II (RPERMP), II(P), II(P+N),
     $              II (P+2*N), II (P+3*N), II (P+4*N), II (OFFPP),
     $              PRESRV, ICNTL)
C             deallocate extra workspace in II and XX
              IHEAD = IHEAD - NZ
              XHEAD = XHEAD - NZ
              IUSE = IUSE - NZ
              XUSE = XUSE - NZ
           ENDIF 

C          -------------------------------------------------------------
C          deallocate workspace, and allocate BTF arrays if required
C          -------------------------------------------------------------

           IF (NBLKS .GT. 1) THEN 
C             replace (6*n+2) workspace at tail of II with
C             Blkp (1..nblks+1) and LUblkp (1..nblks), Offp (1..n+1)
              BLKPP = OFFPP - (NBLKS+1)
              LUBLPP = BLKPP - (NBLKS)
              ITAIL = LUBLPP
              IUSE = IUSE - (6*N+2) + (2*NBLKS+N+2)
           ELSE 
C             The matrix is irreducible.  There is only one block. 
C             Remove everything at tail of II, except
C             for the 2*n permutation vectors and the 7 scalars.
C             (transa, nzoff, nblks, presrv, nz, n, ne).
              ITAIL = (ISIZE + 1) - (2*N+7)
              IUSE = IUSE - (6*N+2)
           ENDIF 

        ENDIF 

C-----------------------------------------------------------------------
C current memory usage:
C-----------------------------------------------------------------------

C       if .not. presrv then
C               input matrix is now in II (1..nz) and XX (1..nz)
C               off-diagonal part: in II/XX (1..nzoff)
C                       col pattern: II (Offp (col) ... Offp (col+1))
C                       col values:  XX (Offp (col) ... Offp (col+1))
C               diagonal blocks: in II/XX (nzoff+1..nz)
C                       col pattern: II (Cp (col) ... Cp (col+1))
C                       col values:  XX (Cp (col) ... Cp (col+1))
C               total: nz+n+1 integers, nz reals
C       else
C               input matrix is now in Ai (1..nz) and Ax (1..nz),
C               in original (non-BTF) order:
C                       col pattern: Ai (Ap (col) ... Ap (col+1))
C                       col values:  Ax (Ap (col) ... Ap (col+1))
C               Cp is a size n+1 integer workspace
C               total: nz+2*(n+1) integers, nz reals
C
C       if (nblks > 1) then
C               at tail of II (in order): 2*nblks+n+2
C                       LUblkp (1..nblks)
C                       Blkp (1..nblks+1)
C                       Offp (1..n+1)
C               total: (2*nblks+n+2) integers
C
C       remainder at tail of II:
C               Cperm (1..n)
C               Rperm (1..n)
C               seven scalars: transa, nzoff, nblks, presrv, nz, n, ne
C
C   Grand total current memory usage (including II,XX,Cp,Ai,Ap,Ax):
C
C       presrv  nblks>1         integers, iuse =
C       F       F               nz+  (n+1)+(2*n+7)
C       F       T               nz+  (n+1)+(2*n+7)+(2*nblks+n+2)
C       T       F               nz+2*(n+1)+(2*n+7)
C       T       T               nz+2*(n+1)+(2*n+7)+(2*nblks+n+2)
C
C   real usage is xuse = nz

C       ----------------------------------------------------------------
C       get memory usage for next call to UMC2RF
C       ----------------------------------------------------------------

        XRMAX = 2*NE
        XRUSE = NZ

C-----------------------------------------------------------------------
C factorization
C-----------------------------------------------------------------------

        IF (NBLKS .EQ. 1) THEN 

C          -------------------------------------------------------------
C          factorize the matrix as a single block
C          -------------------------------------------------------------

           CALL UMC2F1 (CP, N, II (CPERMP), II (RPERMP), NZOFF,
     $          ITAIL, XTAIL, XX, XSIZE, XUSE, II, ITAIL-1, IUSE,
     $          ICNTL, CNTL, INFO, RINFO, NBLKS,
     $          AP, AI, AX, PRESRV, 1, AN, ANZ, II, KEEP,
     $          RMAX, CMAX, TOTNLU, XRMAX, XRUSE, IOUT, XOUT)
           IF (IOUT .OR. XOUT) THEN 
C             error return, if not enough integer and/or real memory:
              GO TO 9000
           ENDIF 
           IF (INFO (1) .LT. 0) THEN 
C             error return, if error in UMC2F2:
              GO TO 9010
           ENDIF 
C          original matrix has been deallocated
           IHEAD = 1
           XHEAD = 1

C          -------------------------------------------------------------
C          make the index of the block relative to start of LU factors
C          -------------------------------------------------------------

           II (ITAIL) = 1

        ELSE 

C          -------------------------------------------------------------
C          factorize the block-upper-triangular form of the matrix
C          -------------------------------------------------------------

           PRP = OFFPP
           IF (PRESRV) THEN 
C             count the off-diagonal entries during factorization
              NZOFF = 0
C             compute temp inverse permutation in II (prp..prp+n-1)
CFPP$ NODEPCHK L
              DO 10 K = 1, N 
                 II (PRP + II (RPERMP+K-1) - 1) = K
10            CONTINUE 
           ENDIF 

           DO 30 BLK = NBLKS, 1, -1 

C             ----------------------------------------------------------
C             factorize the kn-by-kn block, A (k1..k2, k1..k2)
C             ----------------------------------------------------------

C             get k1 and k2, the start and end of this block
              K1 = II (BLKPP+BLK-1)
              K2 = II (BLKPP+BLK) - 1
              KN = K2-K1+1
              IF (.NOT. PRESRV) THEN 
                 P = CP (K1)
                 CP (K2+1) = IHEAD
              ENDIF 

              IF (KN .GT. 1) THEN 

C                -------------------------------------------------------
C                factor the block (the block is not a singleton)
C                -------------------------------------------------------

                 CALL UMC2F1 (CP (K1), KN,
     $              II (CPERMP+K1-1), II (RPERMP+K1-1), NZOFF,
     $              ITAIL, XTAIL, XX, XTAIL-1, XUSE, II, ITAIL-1,
     $              IUSE, ICNTL, CNTL, INFO, RINFO, NBLKS,
     $              AP, AI, AX, PRESRV, K1, AN, ANZ, II (PRP), KEEP,
     $              RMAX, CMAX, TOTNLU, XRMAX, XRUSE, IOUT, XOUT)
                 IF (IOUT .OR. XOUT) THEN 
C                   error return, if not enough int. and/or real memory:
                    GO TO 9000
                 ENDIF 
                 IF (INFO (1) .LT. 0) THEN 
C                   error return, if error in UMC2F2:
                    GO TO 9010
                 ENDIF 
                 IF (PRESRV) THEN 
                    IHEAD = 1
                    XHEAD = 1
                 ELSE 
                    IHEAD = P
                    XHEAD = P
                 ENDIF 

C                -------------------------------------------------------
C                save the location of the LU factors in LUbkp (blk)
C                -------------------------------------------------------

                 II (LUBLPP+BLK-1) = ITAIL

              ELSE 

C                -------------------------------------------------------
C                get the value of singleton at A (k1,k1), if it exists
C                -------------------------------------------------------

                 A = 0
                 IF (PRESRV) THEN 
C                   find the diagonal entry in the unpermuted matrix
                    COL = II (CPERMP + K1 - 1)
                    DO 20 P2 = AP (COL), AP (COL + 1) - 1 
                       ROW = II (PRP + AI (P2) - 1)
                       IF (ROW .LT. K1) THEN 
C                         entry in off-diagonal blocks
                          NZOFF = NZOFF + 1
                       ELSE 
                          A = AX (P2)
                       ENDIF 
20                  CONTINUE 
                    IHEAD = 1
                    XHEAD = 1
                 ELSE IF (P .NE. IHEAD) THEN 
                    A = XX (P)
                    IHEAD = P
                    XHEAD = P
                    IUSE = IUSE - 1
                    XUSE = XUSE - 1
                    XRUSE = XRUSE - 1
                 ENDIF 

C                -------------------------------------------------------
C                store the 1-by-1 LU factors of a singleton
C                -------------------------------------------------------

                 NSGLTN = NSGLTN + 1
                 IF (ABS (A) .EQ. 0) THEN 
C                   the diagonal entry is either not present, or present
C                   but numerically zero.  This is a singular matrix,
C                   replace with 1-by-1 identity matrix.
                    A = 1
                 ELSE 
C                   increment pivot count
                    NPIV = NPIV + 1
                 ENDIF 
                 XTAIL = XTAIL - 1
C                note: if the matrix is not preserved and nonsingular
C                then we will not run out of memory at this point.
                 XUSE = XUSE + 1
                 XRUSE = XRUSE + 1
                 XRMAX = MAX (XRMAX, XRUSE)
                 INFO (20) = MAX (INFO (20), XUSE)
                 INFO (21) = MAX (INFO (21), XUSE)
C                error return, if not enough real memory:
                 IF (XHEAD .GT. XTAIL) THEN 
                    GO TO 9000
                 ENDIF 
                 II (LUBLPP+BLK-1) = -XTAIL
                 XX (XTAIL) = A

              ENDIF 

30         CONTINUE 

C          -------------------------------------------------------------
C          make the index of each block relative to start of LU factors
C          -------------------------------------------------------------

CFPP$ NODEPCHK L
           DO 40 P = LUBLPP, LUBLPP + NBLKS - 1 
              IF (II (P) .GT. 0) THEN 
                 II (II (P)) = II (II (P)) - XTAIL + 1
                 II (P) = II (P) - ITAIL + 1
              ELSE 
C                this is a singleton
                 II (P) = (-II (P)) - XTAIL + 1
              ENDIF 
40         CONTINUE 

C          -------------------------------------------------------------
C          allocate temporary workspace for Pr (1..n) at head of II
C          -------------------------------------------------------------

           PRP = IHEAD
           IHEAD = IHEAD + N
           IUSE = IUSE + N

C          -------------------------------------------------------------
C          allocate a single entry in case the LU factors are empty
C          -------------------------------------------------------------

           IF (NBLKS .EQ. N) THEN 
C             otherwise, arrays in UMC2RF and UMC2SO would have
C             zero size, which can cause an address fault later on
              ITAIL = ITAIL - 1
              IUSE = IUSE + 1
              P2 = ITAIL
           ENDIF 

C          -------------------------------------------------------------
C          allocate permanent copy of off-diagonal blocks
C          -------------------------------------------------------------

           ITAIL = ITAIL - NZOFF
           OFFIP = ITAIL
           XTAIL = XTAIL - NZOFF
           OFFXP = XTAIL
           IUSE = IUSE + NZOFF
           XUSE = XUSE + NZOFF
           XRUSE = XRUSE + NZOFF
           XRMAX = MAX (XRMAX, XRUSE)
           INFO (18) = MAX (INFO (18), IUSE)
           INFO (19) = MAX (INFO (19), IUSE)
           INFO (20) = MAX (INFO (20), XUSE)
           INFO (21) = MAX (INFO (21), XUSE)
           IF (IHEAD .GT. ITAIL .OR. XHEAD .GT. XTAIL) THEN 
C             error return, if not enough integer and/or real memory:
              GO TO 9000
           ENDIF 

C          -------------------------------------------------------------
C          re-order the off-diagonal blocks according to pivot perm
C          -------------------------------------------------------------

C          use Cp as temporary work array:
           MNZ = NZOFF
           IF (PRESRV) THEN 
              CALL UMC2OF (CP, N, II (RPERMP), II (CPERMP), NZOFF,
     $          II (OFFPP), II (OFFIP), XX (OFFXP), II (PRP),
     $          ICNTL, AP, AI, AX, AN, ANZ, PRESRV, NBLKS, II (BLKPP),
     $          MNZ, 1, INFO, P)
           ELSE 
              CALL UMC2OF (CP, N, II (RPERMP), II (CPERMP), NZOFF,
     $          II (OFFPP), II (OFFIP), XX (OFFXP), II (PRP),
     $          ICNTL, 0, II, XX, 0, MNZ, PRESRV, 0, 0,
     $          MNZ, 1, INFO, P)
           ENDIF 
           IF (NBLKS .EQ. N) THEN 
C             zero the only entry in the integer part of the LU factors
              II (P2) = 0
           ENDIF 

C          -------------------------------------------------------------
C          deallocate Pr (1..n), and II/XX (1..nzoff) if present
C          -------------------------------------------------------------

           IHEAD = 1
           XHEAD = 1
           IUSE = IUSE - N
           IF (.NOT. PRESRV) THEN 
              IUSE = IUSE - NZOFF
              XUSE = XUSE - NZOFF
           ENDIF 

        ENDIF 

C-----------------------------------------------------------------------
C  normal and error return
C-----------------------------------------------------------------------

C       error return label:
9000    CONTINUE
        IF (IOUT .OR. IHEAD .GT. ITAIL) THEN 
C          set error flag if not enough integer memory
           CALL UMC2ER (1, ICNTL, INFO, -3, INFO (19))
        ENDIF 
        IF (XOUT .OR. XHEAD .GT. XTAIL) THEN 
C          set error flag if not enough real memory
           CALL UMC2ER (1, ICNTL, INFO, -4, INFO (21))
        ENDIF 

C       error return label, for error from UMC2F2:
9010    CONTINUE

        INFO (4) = 0
        NZDIA = NZ - NZOFF
        INFO (5) = NZ
        INFO (6) = NZDIA
        INFO (7) = NZOFF
        INFO (8) = NSGLTN
        INFO (9) = NBLKS
        INFO (12) = INFO (10) + INFO (11) + N + INFO (7)

C       Count the number of symmetric pivots chosen.  Note that some of
C       these may have been numerically unacceptable.
        NSYM = 0
        IF (INFO (1) .GE. 0) THEN 
           DO 50 K = 1, N 
              IF (II (CPERMP+K-1) .EQ. II (RPERMP+K-1)) THEN 
C                this kth pivot came from the diagonal of A
                 NSYM = NSYM + 1
              ENDIF 
50         CONTINUE 
        ENDIF 
        INFO (16) = NSYM

        INFO (17) = INFO (17) + NPIV
        RINFO (1) = RINFO (4) + RINFO (5) + RINFO (6)

        IF (INFO (1) .GE. 0 .AND. INFO (17) .LT. N) THEN 
C          set warning flag if matrix is singular
           CALL UMC2ER (1, ICNTL, INFO, 4, INFO (17))
        ENDIF 

C       ----------------------------------------------------------------
C       Determine an upper bound on the amount of integer memory needed
C       (LINDEX) for a subsequent call to UMC2RF.  If block-upper-
C       triangular-form is not in use (Info (9) is 1), then
C       this bound is exact.  If NE is higher in the call to UMC2RF
C       than in the call to UMC2FA, then add 3 integers for each
C       additional entry (including the 2 integers required for the
C       row and column indices of the additional triplet itself).
C       This estimate assumes that JOB and TRANSA are the same in
C       UMC2FA and UMC2RF.
C       ----------------------------------------------------------------

C       (Keep (5) - Keep (4) + 1), is added to Info (22)
C       in UMC2FA, to complete the computation of the estimate.

        IF (PRESRV) THEN 
           INFO (22) = MAX (3*NE+2*N+1, NE+3*N+2,
     $                           2*NZ+4*N+10+RMAX+3*CMAX+4*TOTNLU)
        ELSE 
           INFO (22) = MAX (3*NE+2*N+1, NE+3*N+2, 2*NZ+3*N+2,
     $                             NZ+3*N+ 9+RMAX+3*CMAX+4*TOTNLU)
        ENDIF 

C       ----------------------------------------------------------------
C       Approximate the amount of real memory needed (LVALUE) for a
C       subsequent call to UMC2RF.  The approximation is an upper bound
C       on the bare minimum amount needed.  Some garbage collection may
C       occur, but UMC2RF is guaranteed to finish if given an LVALUE of
C       size Info (23) and if the pattern is the same.  If NE is
C       higher in the call to UMC2RF than in the call to UMC2FA, then
C       add 2 reals for each additional entry (including the 1 real
C       required for the value of the additional triplet itself).
C       This estimate assumes that JOB and TRANSA are the same in
C       UMC2FA and UMC2RF.
C       ----------------------------------------------------------------

        INFO (23) = XRMAX
        RETURN
        END 
