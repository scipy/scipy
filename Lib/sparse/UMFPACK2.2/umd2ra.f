        SUBROUTINE UMD2RA (PRESRV, N, NZ, CPERM, RPERM, PR,
     $          W, NBLKS, ARX, ARI, NZOFF, NZDIA,
     $          ICNTL, MP, BLKP, MI, MX, INFO, OFFP, ON, NZBLK,
     $          CBLK, KN, NZ2, NBELOW)
        INTEGER N, NZ, CPERM (N), RPERM (N), PR (N), KN, W (KN+1),
     $          NBLKS, NZBLK, ARI (NZBLK), NZOFF, NZDIA, MP (N+1),
     $          MI (NZ), ON, ICNTL (20), BLKP (NBLKS+1), NZ2,
     $          INFO (40), OFFP (ON+1), CBLK, NBELOW
        LOGICAL PRESRV
        DOUBLE PRECISION
     $          ARX (NZBLK), MX (NZ)
        
C=== UMD2RA ============================================================
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
C  Convert a column-oriented matrix into an arrowhead format.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       n               size of entire matrix
C       Mi (1..nz):     row indices of column form of entire matrix
C       Mx (1..nz):     values of column form of entire matrix
C       Mp (1..n+1)     column pointers for entire matrix
C       Cperm (1..n):   column permutations
C       Rperm (1..n):   row permutations
C
C       if nblks > 1 and presrv
C           cblk:               the block to convert
C           kn:                 the size of the block to convert
C       else
C           cblk:               0
C           kn                  n, size of input matrix

C=======================================================================
C  OUTPUT: 
C=======================================================================
C
C       if nblks = 1 and not presrv
C
C           nzoff               0
C           nzdia               nz - (entries below in diagonal blocks)
C           nz2                 nzdia
C
C           Mi (1..nz2)         arrowheads for the diagonal block
C           Mx (1..nz2)
C           Ari, Arx            used as workspace
C           W (1..n+1)          pointer to each arrowhead in Mi/Mx
C
C           Offp                not accessed
C
C       if nblks = 1 and presrv
C
C           nzoff               0
C           nzdia               nz - (entries below in diagonal blocks)
C           nz2                 nzdia
C
C           Mi, Mx              not modified
C           Ari (1..nz2)        arrowheads for the diagonal block
C           Arx (1..nz2)
C           W (1..n+1)          pointer to each arrowhead in Ari/Arx
C
C           Offp                not accessed
C
C       else if nblks > 1 and not presrv
C
C           nzoff               number of entries in off-diagonal part
C           nzdia               number of entries in diagonal blocks
C                               (nz = nzoff + nzdia + entries below
C                               diagonal blocks)
C           nz2                 nzoff + nzdia
C
C           Mi (nzoff+1..nz2)   arrowheads for each diagonal block
C           Mx (nzoff+1..nz2)
C           Ari, Arx            used as workspace
C           W (1..n+1)          pointer to each arrowhead in Mi/Mx
C
C           Offp (1..n+1)       row pointers for off-diagonal part
C           Mi (1..nzoff)       col indices for off-diagonal part
C           Mx (1..nzoff)       values for off-diagonal part
C
C       else (nblks > 1 and presrv)
C
C           nzoff               0
C           nzdia               nonzeros in the diagonal block, cblk
C           nz2                 nzdia
C
C           Mi, Mx              not modified
C           Ari (1..nz2)        arrowheads for the diagonal block, cblk
C           Arx (1..nz2)
C           W (1..kn+1)         pointer to each arrowhead in Ari/Arx
C
C           Offp                not accessed

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   UMD2R0
C       subroutines called:     UMD2OF
C       functions called:       MIN
        INTRINSIC MIN

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER I, P, ROW, COL, BLK, BASE, K1, K2, K, B1, B2, K0

C  i:       loop index, arrowhead index
C  p:       pointer into column-form input matrix
C  row:     row index
C  col:     column index
C  blk:     current diagonal block
C  base:    where to start the construction of the arrowhead form
C  k1,k2:   current diagonal block is A (k1..k2, k1..k2)
C  k:       loop index, kth pivot
C  b1,b2:   convert blocks b1...b2 from column-form to arrowhead form
C  k0:      convert A (k0+1..., k0+1...) to arrowhead form

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

C-----------------------------------------------------------------------
C  if entire matrix is to be converted, then create the off-diagonal
C  part in row-oriented form in Ari (1..nzoff) and Arx (1..nzoff) and
C  compute inverse row permutation.  Otherwise, the inverse row
C  permutation has already been computed.
C-----------------------------------------------------------------------

        NZOFF = 0
        NBELOW = 0
        IF (NBLKS .EQ. 1) THEN 
           DO 10 K = 1, N 
              PR (RPERM (K)) = K
10         CONTINUE 
        ELSE IF (NBLKS .GT. 1 .AND. .NOT. PRESRV) THEN 
           CALL UMD2OF (W, N, RPERM, CPERM, NZOFF,
     $        OFFP, ARI, ARX, PR,
     $        ICNTL, MP, MI, MX, N, NZ, .TRUE., NBLKS, BLKP,
     $        NZ, 2, INFO, NBELOW)
        ENDIF 

C-----------------------------------------------------------------------
C  construct the arrowhead form for the diagonal block(s)
C-----------------------------------------------------------------------

        DO 20 I = 1, KN+1 
           W (I) = 0
20      CONTINUE 

        BASE = NZOFF + 1

        IF (CBLK .NE. 0) THEN 
C          convert just cblk
           K0 = BLKP (CBLK) - 1
           B1 = CBLK
           B2 = CBLK
        ELSE 
C          convert all the block(s)
           K0 = 0
           B1 = 1
           B2 = NBLKS
        ENDIF 

        DO 80 BLK = B1, B2 

C          -------------------------------------------------------------
C          get the starting and ending indices of this diagonal block
C          -------------------------------------------------------------

           IF (NBLKS .GT. 1) THEN 
              K1 = BLKP (BLK)
              K2 = BLKP (BLK+1) - 1
           ELSE 
              K1 = 1
              K2 = N
           ENDIF 

C          -------------------------------------------------------------
C          count the number of entries in each arrowhead
C          -------------------------------------------------------------

           DO 40 COL = K1, K2 
              DO 30 P = MP (CPERM (COL)), MP (CPERM (COL) + 1) - 1 
                 ROW = PR (MI (P))
                 IF (ROW .GE. K1 .AND. ROW .LE. K2) THEN 
C                   this is in a diagonal block, arrowhead i
                    I = MIN (ROW, COL) - K0
                    W (I) = W (I) + 1
                 ENDIF 
30            CONTINUE 
40         CONTINUE 

C          -------------------------------------------------------------
C          set pointers to point just past end of each arrowhead
C          -------------------------------------------------------------

           W (K2-K0+1) = W (K2-K0) + BASE
           DO 50 I = K2-K0, K1-K0+1, -1 
              W (I) = W (I+1) + W (I-1)
50         CONTINUE 
           W (K1-K0) = W (K1-K0+1)
C          W (i+1-k0) points just past end of arrowhead i in Ari/Arx

C          -------------------------------------------------------------
C          construct arrowhead form, leaving pointers in final state
C          -------------------------------------------------------------

           DO 70 COL = K1, K2 
              DO 60 P = MP (CPERM (COL)), MP (CPERM (COL) + 1) - 1 
                 ROW = PR (MI (P))
                 IF (ROW .GE. K1 .AND. ROW .LE. K2) THEN 
                    IF (ROW .GE. COL) THEN 
C                      diagonal, or lower triangular part
                       I = COL - K0 + 1
                       W (I) = W (I) - 1
                       ARI (W (I)) = ROW - K1 + 1
                       ARX (W (I)) = MX (P)
                    ELSE 
C                      upper triangular part, flag by negating col
                       I = ROW - K0 + 1
                       W (I) = W (I) - 1
                       ARI (W (I)) = -(COL - K1 + 1)
                       ARX (W (I)) = MX (P)
                    ENDIF 
                 ENDIF 
60            CONTINUE 
70         CONTINUE 

           BASE = W (K1-K0)
           W (K2-K0+1) = 0
80      CONTINUE 

        W (KN+1) = NZOFF + 1
        NZDIA = BASE - NZOFF - 1
        NZ2 = NZOFF + NZDIA

C       ----------------------------------------------------------------
C       if cblk = 0, the entire matrix has been converted:
C
C          W (i) now points just past end of arrowhead i in Ari/Arx
C          arrowhead i is located in Ari/Arx (W (i+1) ... W (i)-1),
C          except for the k2-th arrowhead in each block.  Those are
C          located in Ari/Arx (base ... W (k2) - 1), where base is
C          W (Blkp (blk-1)) if blk>1 or W (n+1) = nzoff + 1 otherwise.
C
C       otherwise, just one block has been converted:
C
C          W (i) now points just past end of arrowhead i in Ari/Arx,
C          where i = 1 is the first arrowhead of this block (not the
C          first arrowhead of the entire matrix).  Arrowhead i is
C          located in Ari/Arx (W (i+1) ... W (i)-1).
C          This option is used only if nblks>1 and presrv is true.
C       ----------------------------------------------------------------

C-----------------------------------------------------------------------
C  if not preserved, overwrite column-form with arrowhead form
C-----------------------------------------------------------------------

        IF (.NOT. PRESRV) THEN 
           DO 90 I = 1, NZ 
              MI (I) = ARI (I)
              MX (I) = ARX (I)
90         CONTINUE 
        ENDIF 

        RETURN
        END 
