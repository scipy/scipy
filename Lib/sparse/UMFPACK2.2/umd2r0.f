        SUBROUTINE UMD2R0 (N, NZ, CP, XX, XSIZE, II, ISIZE, XTAIL,
     $          ITAIL, IUSE, XUSE, NZOFF, NBLKS, ICNTL, CNTL, INFO,
     $          RINFO, PRESRV, AP, AI, AX, AN, ANZ, LUI, LUISIZ,
     $          LUBLKP, BLKP, OFFP, ON, CPERM, RPERM, NE)
        INTEGER N, NZ, ISIZE, II (ISIZE), ICNTL (20), INFO (40),
     $          CP (N+1), XSIZE, XTAIL, ITAIL, IUSE, XUSE, AN, ANZ,
     $          AP (AN+1), AI (ANZ), LUISIZ, LUI (LUISIZ), NBLKS,
     $          LUBLKP (NBLKS), BLKP (NBLKS+1), ON, OFFP (ON+1),
     $          CPERM (N), RPERM (N), NZOFF, NE
        LOGICAL PRESRV
        DOUBLE PRECISION
     $          XX (XSIZE), AX (ANZ)
        DOUBLE PRECISION
     $          CNTL (10), RINFO (20)
        
C=== UMD2R0 ============================================================
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
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Refactorize an unsymmetric sparse matrix in column-form, optionally
C  permuting the matrix to upper block triangular form and factorizing
C  each diagonal block.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       n:              order of matrix
C       nz:             entries in matrix
C       Cp (1..n+1):    column pointers of input matrix
C       presrv:         if true then preserve original matrix
C       xsize:          size of XX
C       isize:          size of II
C       iuse:           memory usage in Index on input
C       xuse:           memory usage in Value on input
C       Icntl:          integer control parameters, see UMD21I
C       Cntl:           real control parameters, see UMD21I
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
C           an:                 0
C           anz:                1
C           Ap:                 unused
C           Ai:                 unused
C           Ax:                 unused
C           II (1..nz):         row indices of input matrix
C           XX (1..nz):         values of input matrix
C
C       Information from prior LU factorization:
C
C       LUisiz:                 size of LUi
C       LUi (1..LUisiz):        patterns of LU factors, excluding
C                               prior preserved matrix (if it existed)
C                               and prior off-diagonal part (if it
C                               existed)
C       Cperm (1..n):           column permutations
C       Rperm (1..n):           row permutations
C       nblks:                  number of diagonal blocks for BTF
C       if nblks > 1:
C           LUblkp (1..nblks):  pointers to each diagonal LU factors
C           Blkp (1..nblks+1):  index range of diagonal blocks

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
C                       the off-diagonal nonzeros, if nblks > 1
C
C       Offp (1..n+1):  row pointers for off-diagonal part, if nblks > 1
C       Info:           integer informational output, see UMD2FA
C       Rinfo:          real informational output, see UMD2FA

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   UMD2RF
C       subroutines called:     UMD2ER, UMD2R2, UMD2RA, UMD2OF
C       functions called:       MAX
        INTRINSIC MAX

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER I, NZDIA, P, IHEAD, NSGLTN, NSYM, WP, ARIP, ARXP, NPIV,
     $          WRKSIZ, NLU, PRP, MC, MR, DUMMY1, DUMMY2, NZ2, K, BLK,
     $          K1, K2, KN, NZBLK, COL, ROW, PRL, IO, LUIP, MNZ, ARNZ,
     $          XHEAD, OFFIP, OFFXP, NOUTSD, NBELOW, NZORIG, XRMAX
        DOUBLE PRECISION
     $          A

C  Printing control:
C  -----------------
C  io:      I/O unit for diagnostic messages
C  prl:     printing level
C
C  Allocated array pointers:
C  -------------------------
C  wp:      W (1..n+1), or W (1..kn+1), work array located in II (wp...)
C  prp:     Pr (1..n) work array located in II (prp..prp+n-1)
C  arip:    Ari (1..nz) array located in II (arip..arip+nz-1)
C  arxp:    Arx (1..nz) array located in XX (arxp..arxp+nz-1)
C  offip:   Offi (1..nzoff) array located in II (offip..offip+nzoff-1)
C  offxp:   Offx (1..nzoff) array located in XX (offxp..offip+nzoff-1)
C
C  Arrowhead-form matrix:
C  ----------------------
C  nz2:     number of entries in arrowhead matrix
C  nzblk:   number of entries in arrowhead matrix of a single block
C  arnz:    arrowhead form of blocks located in II/XX (1..arnz)
C
C  BTF information:
C  ----------------
C  k1:      starting index of diagonal block being factorized
C  k2:      ending index of diagonal block being factorized
C  kn:      the order of the diagonal block being factorized
C  blk:     block number of diagonal block being factorized
C  nsgltn:  number of 1-by-1 diagonal blocks ("singletons")
C  a:       numerical value of a singleton
C  mnz:     nzoff
C  noutsd:  entries in diagonal blocks, but not in LU (invalid)
C  nbelow:  entries below diagonal blocks (invalid)
C  nzoff:   entries above diagonal blocks (valid)
C  nzdia:   entries in diagonal blocks (valid)
C  nzorig:  total number of original entries
C
C  Memory usage:
C  -------------
C  xhead:   XX (1..xhead-1) is in use, XX (xhead..xtail-1) is free
C  ihead:   II (1..ihead-1) is in use, II (ihead..itail-1) is free
C  wrksiz:  total size of work arrays need in II for call to UMD2R2
C  xrmax:   memory needed in Value for next call to UMD2RF
C
C  Symbolic information and pattern of prior LU factors:
C  -----------------------------------------------------
C  nlu:     number of elements in a diagonal block
C  luip:    integer part of LU factors located in LUi (luip...)
C  mr,mc:   largest frontal matrix for this diagonal block is mc-by-mr
C
C  Other:
C  ------
C  k:       loop index (kth pivot)
C  i:       loop index
C  row:     row index
C  col:     column index
C  p:       pointer
C  nsym:    number of symmetric pivots chosen
C  dummy1:  argument returned by UMD2RA, but not needed
C  dummy2:  argument returned by UMD2RA, but not needed

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

        IO = ICNTL (2)
        PRL = ICNTL (3)
        NZORIG = NZ

        IF (PRESRV) THEN 
C          original matrix is not in Cp/II/XX, but in Ap/Ai/Ax:
           IHEAD = 1
           XHEAD = 1
        ELSE 
           IHEAD = NZ + 1
           XHEAD = NZ + 1
        ENDIF 

        NZOFF = 0
        NZDIA = NZ
        NSGLTN = 0
        NPIV = 0
        NOUTSD = 0
        NBELOW = 0
        ITAIL = ISIZE + 1
        XTAIL = XSIZE + 1

C-----------------------------------------------------------------------
C current memory usage:
C-----------------------------------------------------------------------

C       if .not. presrv then
C               input matrix is now in II (1..nz) and XX (1..nz)
C                       col pattern: II (Cp (col) ... Cp (col+1))
C                       col values:  XX (Cp (col) ... Cp (col+1))
C               total: nz+n+1 integers, nz reals
C       else
C               input matrix is now in Ai (1..nz) and Ax (1..nz)
C                       col pattern: Ai (Ap (col) ... Ap (col+1))
C                       col values:  Ax (Ap (col) ... Ap (col+1))
C               Cp is a size n+1 integer workspace
C               total: nz+2*(n+1) integers, nz reals
C
C       if (nblks > 1) then
C                       LUblkp (1..nblks)
C                       Blkp (1..nblks+1)
C                       Offp (1..n+1)
C               total: (2*nblks+n+2) integers
C
C       Cperm (1..n) and Rperm (1..n)
C               total: 2*n integers
C
C   Grand total current memory usage (including II,XX,Cp,Ai,Ap,Ax
C       and LUi):
C
C       presrv  nblks>1 integers, iuse = 
C       F       F       LUisiz + nz+  (n+1)+(2*n+7)
C       F       T       LUisiz + nz+  (n+1)+(2*n+7)+(2*nblks+n+2)
C       T       F       LUisiz + nz+2*(n+1)+(2*n+7)
C       T       T       LUisiz + nz+2*(n+1)+(2*n+7)+(2*nblks+n+2)
C
C   real usage is xuse = nz

C-----------------------------------------------------------------------
C  get memory usage estimate for next call to UMD2RF
C-----------------------------------------------------------------------

        XRMAX = 2*NE

C-----------------------------------------------------------------------
C  convert matrix into arrowhead format (unless BTF and preserved)
C-----------------------------------------------------------------------

        IF (NBLKS .GT. 1 .AND. PRESRV) THEN 

C          -------------------------------------------------------------
C          BTF is to be used, and original matrix is to be preserved.
C          It is converted and factorized on a block-by-block basis,
C          using the inverse row permutation (computed and stored in
C          Offp (1..n)
C          -------------------------------------------------------------

           DO 10 K = 1, N 
              OFFP (RPERM (K)) = K
10         CONTINUE 

        ELSE 

C          -------------------------------------------------------------
C          convert the entire input matrix to arrowhead form
C          -------------------------------------------------------------

C          -------------------------------------------------------------
C          allocate workspace: W (n+1), Pr (n), Ari (nz), Arx (nz)
C          -------------------------------------------------------------

           ITAIL = ITAIL - (2*N+1)
           IUSE = IUSE + 2*N+1
           PRP = ITAIL
           WP = PRP + N
           IUSE = IUSE + NZ
           XUSE = XUSE + NZ
           ARXP = XHEAD
           ARIP = IHEAD
           IHEAD = IHEAD + NZ
           XHEAD = XHEAD + NZ
           INFO (18) = MAX (INFO (18), IUSE)
           INFO (20) = MAX (INFO (20), XUSE)
           INFO (21) = MAX (INFO (21), XUSE)
           IF (IHEAD .GT. ITAIL .OR. XHEAD .GT. XTAIL) THEN 
C             error return, if not enough integer and/or real memory:
              GO TO 9000
           ENDIF 

C          -------------------------------------------------------------
C          convert
C          -------------------------------------------------------------

           IF (NBLKS .EQ. 1) THEN 
              IF (PRESRV) THEN 
                 CALL UMD2RA (PRESRV, N, NZ, CPERM, RPERM, II (PRP),
     $              II (WP), NBLKS, XX (ARXP), II (ARIP), NZOFF, NZDIA,
     $              ICNTL, AP, BLKP, AI, AX, INFO, OFFP, ON, NZ,
     $              0, N, NZ2, I)
              ELSE 
                 CALL UMD2RA (PRESRV, N, NZ, CPERM, RPERM, II (PRP),
     $              II (WP), NBLKS, XX (ARXP), II (ARIP), NZOFF, NZDIA,
     $              ICNTL, CP, BLKP, II, XX, INFO, OFFP, ON, NZ,
     $              0, N, NZ2, I)
              ENDIF 
           ELSE 
C             note that presrv is false in this case
              CALL UMD2RA (PRESRV, N, NZ, CPERM, RPERM, II (PRP),
     $           II (WP), NBLKS, XX (ARXP), II (ARIP), NZOFF, NZDIA,
     $           ICNTL, CP, BLKP, II, XX, INFO, OFFP, ON, NZ,
     $           0, N, NZ2, NBELOW)
           ENDIF 

C          -------------------------------------------------------------
C          copy the arrowhead pointers from W (1..n+1) to Cp (1..n+1)
C          -------------------------------------------------------------

           DO 20 I = 1, N+1 
              CP (I) = II (WP+I-1)
20         CONTINUE 

C          -------------------------------------------------------------
C          deallocate W and Pr.  If not presrv deallocate Ari and Arx
C          -------------------------------------------------------------

           IUSE = IUSE - (2*N+1)
           IF (.NOT. PRESRV) THEN 
C             Ari and Arx have been deallocated.
              XUSE = XUSE - NZ
              IUSE = IUSE - NZ
           ENDIF 
           ITAIL = ISIZE + 1
           XTAIL = XSIZE + 1
           NZ = NZ2
           IHEAD = NZ + 1
           XHEAD = NZ + 1

        ENDIF 

        INFO (5) = NZ
        INFO (6) = NZDIA
        INFO (7) = NZOFF
        INFO (4) = NBELOW

C-----------------------------------------------------------------------
C  refactorization
C-----------------------------------------------------------------------

C       ----------------------------------------------------------------
C       if nblks=1
C          arrowhead form is now stored in II (1..nz) and XX (1..nz)
C          in reverse pivotal order (arrowhead n, n-1, ..., 2, 1).
C          The arrowhead form will be overwritten.
C       else if not presrv
C          off-diagonal part is in II (1..nzoff), XX (1..nzoff),
C          (with row pointers Offp (1..n+1)) followed by each diagonal
C          block (block 1, 2, ... nblks) in II/XX (nzoff+1...nz).
C          Each diagonal block is in arrowhead form, and in
C          reverse pivotal order (arrowhead k2, k2-1, ..., k1-1, k1).
C          The arrowhead form will be overwritten.
C       else (nblks > 1 and presrv)
C          II and XX are still empty.  Original matrix is in Ap, Ai,
C          and Ax.  Inverse row permutation (Pr) is in Offp (1..n).
C          The arrowhead form is not yet computed.
C       ----------------------------------------------------------------

        IF (NBLKS .EQ. 1) THEN 

C          -------------------------------------------------------------
C          refactorize the matrix as a single block
C          -------------------------------------------------------------

           NLU = LUI (2)
           MC = LUI (4)
           MR = LUI (5)
           WRKSIZ = 2*N + MR + 3*MC + 4*(NLU+2)
           ITAIL = ITAIL - WRKSIZ
           IUSE = IUSE + WRKSIZ
           P = ITAIL
           INFO (18) = MAX (INFO (18), IUSE)
           IF (IHEAD .GT. ITAIL) THEN 
C             error return, if not enough integer memory:
              GO TO 9000
           ENDIF 

           CALL UMD2R2 (CP, NZ, N, XTAIL,
     $          XX, XSIZE, XUSE, II, CPERM, RPERM,
     $          ICNTL, CNTL, INFO, RINFO, MC, MR,
     $          II (P), II (P+N), II (P+2*N), II (P+2*N+MR),
     $          II (P+2*N+MR+MC), II (P+2*N+MR+2*MC),
     $          II (P+2*N+MR+3*MC), II (P+2*N+MR+3*MC+(NLU+2)),
     $          II (P+2*N+MR+3*MC+2*(NLU+2)),
     $          II (P+2*N+MR+3*MC+3*(NLU+2)),
     $          NLU, LUI (6), LUI (NLU+6), NOUTSD,
     $          XRMAX)

           IF (INFO (1) .LT. 0) THEN 
C             error return, if not enough real memory or bad pivot found
              GO TO 9010
           ENDIF 

C          -------------------------------------------------------------
C          deallocate workspace and original matrix (reals already done)
C          -------------------------------------------------------------

           IUSE = IUSE - WRKSIZ - NZ
           ITAIL = ITAIL + WRKSIZ
           LUI (1) = 1
           IHEAD = 1
           XHEAD = 1

        ELSE 

C          -------------------------------------------------------------
C          refactorize the block-upper-triangular form of the matrix
C          -------------------------------------------------------------

           IF (PRESRV) THEN 
C             count the entries in off-diagonal part
              NZOFF = 0
           ENDIF 

           DO 70 BLK = NBLKS, 1, -1 

C             ----------------------------------------------------------
C             factorize the kn-by-kn block, A (k1..k2, k1..k2)
C             ----------------------------------------------------------

C             get k1 and k2, the start and end of this block
              K1 = BLKP (BLK)
              K2 = BLKP (BLK+1) - 1
              KN = K2-K1+1
              A = 0

C             ----------------------------------------------------------
C             get pointers to, or place the block in, arrowhead form
C             ----------------------------------------------------------

              IF (PRESRV) THEN 

                 IF (KN .GT. 1) THEN 

C                   ----------------------------------------------------
C                   convert a single block to arrowhead format, using
C                   the inverse row permutation stored in Offp
C                   ----------------------------------------------------

C                   ----------------------------------------------------
C                   compute nzblk, allocate II/XX (1..nzblk), W(1..kn+1)
C                   ----------------------------------------------------

                    NZBLK = 0
                    DO 40 K = K1, K2 
                       COL = CPERM (K)
CFPP$ NODEPCHK L
                       DO 30 P = AP (COL), AP (COL+1) - 1 
                          ROW = OFFP (AI (P))
                          IF (ROW .LT. K1) THEN 
C                            entry in off-diagonal part
                             NZOFF = NZOFF + 1
                          ELSE IF (ROW .LE. K2) THEN 
C                            entry in diagonal block
                             NZBLK = NZBLK + 1
                          ENDIF 
30                     CONTINUE 
40                  CONTINUE 

                    ITAIL = ITAIL - (KN+1)
                    WP = ITAIL
                    IHEAD = NZBLK + 1
                    XHEAD = NZBLK + 1
                    IUSE = IUSE + NZBLK + KN+1
                    XUSE = XUSE + NZBLK
                    XRMAX = MAX (XRMAX, XUSE)
                    INFO (18) = MAX (INFO (18), IUSE)
                    INFO (20) = MAX (INFO (20), XUSE)
                    INFO (21) = MAX (INFO (21), XUSE)
                    IF (IHEAD .GT. ITAIL .OR. XHEAD .GT. XTAIL) THEN 
C                      error return, if not enough integer
C                      and/or real memory:
                       GO TO 9000
                    ENDIF 

C                   ----------------------------------------------------
C                   convert blk from column-form in Ai/Ax to arrowhead
C                   form in II/XX (1..nzblk)
C                   ----------------------------------------------------

                    CALL UMD2RA (PRESRV, N, NZ, CPERM, RPERM, OFFP,
     $                 II (WP), NBLKS, XX, II, DUMMY1, DUMMY2,
     $                 ICNTL, AP, BLKP, AI, AX, INFO, 0, 0, NZBLK,
     $                 BLK, KN, NZ2, I)

C                   ----------------------------------------------------
C                   copy the arrowhead pointers from W (1..kn+1)
C                   to Cp (k1 ... k2+1)
C                   ----------------------------------------------------

                    DO 50 I = 0, KN 
                       CP (K1+I) = II (WP+I)
50                  CONTINUE 

C                   Cp (k1) is nzblk + 1 and Cp (k2+1) is 1

C                   ----------------------------------------------------
C                   deallocate W (1..kn+1)
C                   ----------------------------------------------------

                    IUSE = IUSE - (KN+1)
                    ITAIL = ITAIL + (KN+1)

                 ELSE 

C                   ----------------------------------------------------
C                   get the value of singleton at A (k1,k1) if it exists
C                   ----------------------------------------------------

C                   find the diagonal entry in the unpermuted matrix,
C                   and count the entries in the diagonal and
C                   off-diagonal blocks.
                    COL = CPERM (K1)
                    DO 60 P = AP (COL), AP (COL + 1) - 1 
C                      inverse row permutation is stored in Offp
                       ROW = OFFP (AI (P))
                       IF (ROW .LT. K1) THEN 
                          NZOFF = NZOFF + 1
                       ELSE IF (ROW .EQ. K1) THEN 
                          A = AX (P)
C                      else 
C                         this is an invalid entry, below the diagonal
C                         block.  It will be detected (and optionally
C                         printed) in the call to UMD2OF below.
                       ENDIF 
60                  CONTINUE 

                    IHEAD = 1
                    XHEAD = 1
                 ENDIF 

              ELSE 

C                -------------------------------------------------------
C                the block is located in II/XX (Cp (k2+1) ... Cp (k1)-1)
C                and has already been converted to arrowhead form
C                -------------------------------------------------------

                 IF (BLK .EQ. 1) THEN 
C                   this is the last block to factorize
                    CP (K2+1) = NZOFF + 1
                 ELSE 
                    CP (K2+1) = CP (BLKP (BLK-1))
                 ENDIF 

                 IHEAD = CP (K1)
                 XHEAD = IHEAD

                 IF (KN .EQ. 1) THEN 
C                   singleton block in II/XX (Cp (k1+1) ... Cp (k1)-1)
                    IF (CP (K1) .GT. CP (K1+1)) THEN 
                       A = XX (CP (K1) - 1)
                       IHEAD = IHEAD - 1
                       XHEAD = XHEAD - 1
                       IUSE = IUSE - 1
                       XUSE = XUSE - 1
                    ENDIF 
                 ENDIF 
                    
              ENDIF 

C             ----------------------------------------------------------
C             factor the block
C             ----------------------------------------------------------

              IF (KN .GT. 1) THEN 

C                -------------------------------------------------------
C                The A (k1..k2, k1..k2) block is not a singleton.
C                block is now in II/XX (Cp (k2+1) ... Cp (k1)-1), in
C                arrowhead form, and is to be overwritten with LU
C                -------------------------------------------------------

                 ARNZ = CP (K1) - 1

C                if (presrv) then 
C                   II/XX (1..arnz) holds just the current block, blk
C                else 
C                   II/XX (1..arnz) holds the off-diagonal part, and
C                   blocks 1..blk, in that order.
C                endif 

                 LUIP = LUBLKP (BLK)
C                luxp = LUi (luip), not needed for refactorization
                 NLU = LUI (LUIP+1)
C                npiv = LUi (luip+2), not needed for refactorization
                 MC = LUI (LUIP+3)
                 MR = LUI (LUIP+4)
                 WRKSIZ = 2*KN + MR + 3*MC + 4*(NLU+2)
                 ITAIL = ITAIL - WRKSIZ
                 IUSE = IUSE + WRKSIZ
                 P = ITAIL
                 INFO (18) = MAX (INFO (18), IUSE)
                 IF (IHEAD .GT. ITAIL) THEN 
C                   error return, if not enough integer memory:
                    GO TO 9000
                 ENDIF 

                 CALL UMD2R2 (CP (K1), ARNZ, KN, XTAIL,
     $                XX, XTAIL-1, XUSE, II, CPERM (K1), RPERM (K1),
     $                ICNTL, CNTL, INFO, RINFO, MC, MR,
     $                II (P), II (P+KN), II (P+2*KN), II (P+2*KN+MR),
     $                II (P+2*KN+MR+MC), II (P+2*KN+MR+2*MC),
     $                II (P+2*KN+MR+3*MC), II (P+2*KN+MR+3*MC+(NLU+2)),
     $                II (P+2*KN+MR+3*MC+2*(NLU+2)),
     $                II (P+2*KN+MR+3*MC+3*(NLU+2)),
     $                NLU, LUI (LUIP+5), LUI (LUIP+NLU+5), NOUTSD,
     $                XRMAX)

                 IF (INFO (1) .LT. 0) THEN 
C                   error return, if not enough real memory or bad pivot
                    GO TO 9010
                 ENDIF 

C                -------------------------------------------------------
C                deallocate workspace and original matrix (reals
C                already deallocated in UMD2R2)
C                -------------------------------------------------------

                 IUSE = IUSE - WRKSIZ
                 ITAIL = ITAIL + WRKSIZ
                 LUI (LUIP) = XTAIL
                 IUSE = IUSE - (IHEAD - CP (K2+1))
                 IHEAD = CP (K2+1)
                 XHEAD = IHEAD

              ELSE 

C                -------------------------------------------------------
C                factor the singleton A (k1,k1) block, in a
C                -------------------------------------------------------

                 NSGLTN = NSGLTN + 1
                 IF (ABS (A) .EQ. 0) THEN 
C                   this is a singular matrix, replace with 1-by-1
C                   identity matrix.
                    A = 1
                 ELSE 
C                   increment pivot count
                    NPIV = NPIV + 1
                 ENDIF 
                 XTAIL = XTAIL - 1
                 XUSE = XUSE + 1
                 XRMAX = MAX (XRMAX, XUSE)
                 INFO (20) = MAX (INFO (20), XUSE)
                 INFO (21) = MAX (INFO (21), XUSE)
C                note: if the matrix is not preserved and nonsingular
C                then we will not run out of memory
                 IF (XHEAD .GT. XTAIL) THEN 
C                   error return, if not enough real memory:
                    GO TO 9000
                 ENDIF 

C                -------------------------------------------------------
C                store the 1-by-1 LU factors
C                -------------------------------------------------------

                 XX (XTAIL) = A
                 LUBLKP (BLK) = -XTAIL

              ENDIF 
70         CONTINUE 

C          -------------------------------------------------------------
C          make the index of each block relative to start of LU factors
C          -------------------------------------------------------------
CFPP$ NODEPCHK L
           DO 80 BLK = 1, NBLKS 
              IF (LUBLKP (BLK) .GT. 0) THEN 
                 LUI (LUBLKP (BLK)) = LUI (LUBLKP (BLK)) - XTAIL + 1
              ELSE 
C                this is a singleton
                 LUBLKP (BLK) = (-LUBLKP (BLK)) - XTAIL + 1
              ENDIF 
80         CONTINUE 

C          -------------------------------------------------------------
C          store the off-diagonal blocks
C          -------------------------------------------------------------

           IF (PRESRV) THEN 

C             ----------------------------------------------------------
C             allocate temporary workspace for Pr (1..n) at head of II
C             ----------------------------------------------------------

              PRP = IHEAD
              IHEAD = IHEAD + N
              IUSE = IUSE + N

C             ----------------------------------------------------------
C             allocate permanent copy of off-diagonal blocks
C             ----------------------------------------------------------

              ITAIL = ITAIL - NZOFF
              OFFIP = ITAIL
              XTAIL = XTAIL - NZOFF
              OFFXP = XTAIL
              IUSE = IUSE + NZOFF
              XUSE = XUSE + NZOFF
              XRMAX = MAX (XRMAX, XUSE)
              INFO (18) = MAX (INFO (18), IUSE)
              INFO (20) = MAX (INFO (20), XUSE)
              INFO (21) = MAX (INFO (21), XUSE)
              IF (IHEAD .GT. ITAIL .OR. XHEAD .GT. XTAIL) THEN 
C                error return, if not enough integer and/or real memory:
                 GO TO 9000
              ENDIF 

C             ----------------------------------------------------------
C             re-order the off-diagonal blocks according to pivot perm
C             ----------------------------------------------------------

C             use Cp as temporary work array:
              MNZ = NZOFF
              IF (NZOFF .EQ. 0) THEN 
C                Offi and Offx are not accessed in UMD2OF.  Set offip
C                and offxp to 1 (since offip = itail = isize+1, which
C                can generate an address fault, otherwise).
                 OFFIP = 1
                 OFFXP = 1
              ENDIF 
              CALL UMD2OF (CP, N, RPERM, CPERM, NZOFF,
     $             OFFP, II (OFFIP), XX (OFFXP), II (PRP),
     $             ICNTL, AP, AI, AX, AN, ANZ, PRESRV, NBLKS, BLKP,
     $             MNZ, 2, INFO, NBELOW)

C             ----------------------------------------------------------
C             deallocate Pr (1..n)
C             ----------------------------------------------------------

              IHEAD = 1
              XHEAD = 1
              IUSE = IUSE - N

           ELSE 

C             off-diagonal entries are in II/XX (1..nzoff); shift down
C             to II/XX ( ... itail/xtail).  No extra memory needed.
              DO 90 I = NZOFF, 1, -1 
                 II (ITAIL+I-NZOFF-1) = II (I)
                 XX (XTAIL+I-NZOFF-1) = XX (I)
90            CONTINUE 
              IHEAD = 1
              XHEAD = 1
              ITAIL = ITAIL - NZOFF
              XTAIL = XTAIL - NZOFF
           ENDIF 

        ENDIF 

C       ----------------------------------------------------------------
C       clear the flags (negated row/col indices, and negated ludegr/c)
C       ----------------------------------------------------------------

        DO 100 I = 1, LUISIZ 
           LUI (I) = ABS (LUI (I))
100     CONTINUE 

C-----------------------------------------------------------------------
C  normal and error return
C-----------------------------------------------------------------------

C       error return label:
9000    CONTINUE
        IF (IHEAD .GT. ITAIL) THEN 
C          set error flag if not enough integer memory
           CALL UMD2ER (2, ICNTL, INFO, -3, INFO (18))
        ENDIF 
        IF (XHEAD .GT. XTAIL) THEN 
C          set error flag if not enough real memory
           CALL UMD2ER (2, ICNTL, INFO, -4, INFO (21))
        ENDIF 

C       error return label, for error return from UMD2R2:
9010    CONTINUE

C       Cp can now be deallocated in UMD2RF:
        IUSE = IUSE - (N+1)

        INFO (4) = NOUTSD + NBELOW
        NZDIA = NZORIG - NZOFF - NOUTSD - NBELOW
        INFO (5) = NZOFF + NZDIA
        INFO (6) = NZDIA
        INFO (7) = NZOFF
        INFO (8) = NSGLTN
        INFO (9) = NBLKS
        INFO (12) = INFO (10) + INFO (11) + N + INFO (7)

C       Count the number of symmetric pivots chosen.  Note that some
C       of these may have been numerically unacceptable.
        NSYM = 0
        DO 110 K = 1, N 
           IF (CPERM (K) .EQ. RPERM (K)) THEN 
C             this kth pivot came from the diagonal of A
              NSYM = NSYM + 1
           ENDIF 
110     CONTINUE 
        INFO (16) = NSYM

        INFO (17) = INFO (17) + NPIV
        RINFO (1) = RINFO (4) + RINFO (5) + RINFO (6)

C       set warning flag if entries outside prior pattern are present
        IF (INFO (4) .GT. 0) THEN 
           CALL UMD2ER (2, ICNTL, INFO, 1, -INFO (4))
        ENDIF 

C       set warning flag if matrix is singular
        IF (INFO (1) .GE. 0 .AND. INFO (17) .LT. N) THEN 
           CALL UMD2ER (2, ICNTL, INFO, 4, INFO (17))
        ENDIF 

C       ----------------------------------------------------------------
C       return memory usage estimate for next call to UMD2RF
C       ----------------------------------------------------------------

        INFO (23) = XRMAX

        RETURN
        END 
