        SUBROUTINE UMD2OF (W, N, RPERM, CPERM, NZOFF,
     $          OFFP, OFFI, OFFX, PR,
     $          ICNTL, MP, MI, MX, MN, MNZ, PRESRV, NBLKS, BLKP,
     $          ONZ, WHO, INFO, NBELOW)
        INTEGER N, NZOFF, W (N+1), RPERM (N), CPERM (N), ONZ,
     $          OFFP (N+1), OFFI (ONZ), PR (N), ICNTL (20), MN, MNZ,
     $          MP (MN+1), MI (MNZ), NBLKS, BLKP (NBLKS+1), WHO, NBELOW,
     $          INFO (40)
        LOGICAL PRESRV
        DOUBLE PRECISION
     $          OFFX (ONZ), MX (MNZ)
        
C=== UMD2OF ============================================================
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
C  Permute the off-diagonal blocks according to final pivot permutation.
C  This routine is called only if the block-triangular-form (BTF) is
C  used.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       n:              order of matrix
C       Rperm (1..n):   the final row permutations, including BTF
C                       If i is the k-th pivot row, then Rperm (k) = i
C       Cperm (1..n):   the final column permutations, including BTF
C                       If j is the k-th pivot col, then Cperm (k) = j
C       Icntl:          integer control parameters, see UMD21I
C       Info:           integer informational parameters
C       who:            who called (1: UMD2FA, 2: UMD2RF)
C
C       if presrv is true then
C           mn:                 order of preserved matrix
C           mnz:                number of entries in preserved matrix
C           Mp (1..mn+1):       column pointers of preserved matrix
C           Mi (1..mnz):        row indices of preserved matrix
C           Mx (1..mnz):        values of preserved matrix
C           Blkp (1..nblks+1):  the index range of the blocks
C           nblks:              the number of diagonal blocks
C       else
C           mn:                 0
C           mnz:                nzoff
C           Mp:                 unaccessed
C           Offp (1..n+1):      column pointers for off-diagonal entries
C                               in original order
C           Mi (1..mnz):        the row indices of off-diagonal entries,
C                               in original order
C           Mx (1..mnz):        the values of off-diagonal entries,
C                               in original order
C           nblks:              0
C           Blkp (1..nblks+1):  unaccessed
C           nzoff:              number of entries in off-diagonal blocks

C=======================================================================
C  WORKSPACE:
C=======================================================================
C
C       W (1..n)

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       Offp (1..n+1):          row pointers for off-diagonal part
C       Offi (1..nzoff):        column indices in off-diagonal part
C       Offx (1..nzoff):        values in off-diagonal part
C       nzoff:                  number of entries in off-diagonal blocks
C       Pr (1..n):              inverse row permutation
C       nbelow:                 entries that are below the diagonal
C                               blocks (can only occur if who = 2)

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   UMD2F0, UMD2RA, UMD2R0
C       subroutines called:     UMD2P2

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER ROW, COL, P, BLK, K, K1, K2, IO, PRL
        LOGICAL PR3

C  row:     row index
C  col:     column index
C  p:       pointer
C  blk:     current diagonal block
C  k:       kth pivot
C  k1,k2:   current diaogonal block is A (k1..k2, k1..k2)
C  io:      I/O unit for diagnostic messages
C  prl:     printing level
C  pr3:     true if printing entries below diagonal blocks (UMD2RF)

C=======================================================================
C  EXECUTABLE STATMENTS:
C=======================================================================

        IO = ICNTL (2)
        PRL = ICNTL (3)
        PR3 = PRL .GE. 3 .AND. IO .GE. 0

C-----------------------------------------------------------------------
C  compute inverse row permutation
C-----------------------------------------------------------------------

C       if original row i is the kth pivot row, then
C               Rperm (k) = i
C               Pr (i) = k
C       if original col j is the kth pivot col, then
C               Cperm (k) = j
CFPP$ NODEPCHK L
        DO 10 K = 1, N 
           PR (RPERM (K)) = K
10      CONTINUE 

C-----------------------------------------------------------------------
C  construct row-oriented pointers for permuted row-form
C-----------------------------------------------------------------------

        W (1) = 1
        DO 20 ROW = 2, N 
           W (ROW) = 0
20      CONTINUE 
        NBELOW = 0
        IF (PRESRV) THEN 
           DO 50 BLK = 1, NBLKS 
              K1 = BLKP (BLK)
              K2 = BLKP (BLK+1) - 1
              DO 40 COL = K1, K2 
CFPP$ NODEPCHK L
                 DO 30 P = MP (CPERM (COL)), MP (CPERM (COL)+1)-1 
                    ROW = PR (MI (P))
                    IF (ROW .LT. K1) THEN 
C                      offdiagonal entry
                       W (ROW) = W (ROW) + 1
                    ELSE IF (ROW .GT. K2 .AND. WHO .EQ. 2) THEN 
C                      This entry is below the diagonal block - invalid.
C                      This can only occur if who = 2 (UMD2RF).
                       IF (PR3) THEN 
C                         print the original row and column indices:
                          CALL UMD2P2 (2, 96, MI (P), COL, MX (P), IO)
                       ENDIF 
                       NBELOW = NBELOW + 1
                    ENDIF 
30               CONTINUE 
40            CONTINUE 
50         CONTINUE 
        ELSE 
           DO 70 COL = 1, N 
CFPP$ NODEPCHK L
              DO 60 P = OFFP (COL), OFFP (COL+1) - 1 
                 ROW = PR (MI (P))
                 W (ROW) = W (ROW) + 1
60            CONTINUE 
70         CONTINUE 
        ENDIF 
        DO 80 ROW = 2, N 
           W (ROW) = W (ROW) + W (ROW-1)
80      CONTINUE 
        W (N+1) = W (N)
C       W (row) now points just past end of row in Offi/x

C-----------------------------------------------------------------------
C  construct the row-oriented form of the off-diagonal values,
C  in the final pivot order.  The column indices in each row
C  are placed in ascending order (the access of Offi/Offx later on
C  does not require this, but it makes access more efficient).
C-----------------------------------------------------------------------

        IF (PRESRV) THEN 
           DO 110 BLK = NBLKS, 1, -1 
              K1 = BLKP (BLK)
              K2 = BLKP (BLK+1) - 1
              DO 100 COL = K2, K1, - 1 
CFPP$ NODEPCHK L
                 DO 90 P = MP (CPERM (COL)), MP (CPERM (COL)+1)-1 
                    ROW = PR (MI (P))
                    IF (ROW .LT. K1) THEN 
C                      offdiagonal entry
                       W (ROW) = W (ROW) - 1
                       OFFI (W (ROW)) = COL
                       OFFX (W (ROW)) = MX (P)
                    ENDIF 
90               CONTINUE 
100           CONTINUE 
110        CONTINUE 
        ELSE 
           DO 130 COL = N, 1, -1 
CFPP$ NODEPCHK L
              DO 120 P = OFFP (CPERM (COL)), OFFP (CPERM (COL) + 1) - 1
                 ROW = PR (MI (P))
                 W (ROW) = W (ROW) - 1
                 OFFI (W (ROW)) = COL
                 OFFX (W (ROW)) = MX (P)
120           CONTINUE 
130        CONTINUE 
        ENDIF 

C-----------------------------------------------------------------------
C  save the new row pointers
C-----------------------------------------------------------------------

        DO 140 ROW = 1, N+1 
           OFFP (ROW) = W (ROW)
140     CONTINUE 

        NZOFF = OFFP (N+1) - 1

        RETURN
        END 
