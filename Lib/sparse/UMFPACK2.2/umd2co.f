        SUBROUTINE UMD2CO (N, NZ, TRANSA, XX, XSIZE, INFO, ICNTL,
     $          II, ISIZE, W, WP, WHO)
        INTEGER ISIZE, II (ISIZE), N, NZ, W (N), WP (N+1), INFO (40),
     $          ICNTL (20), XSIZE, WHO
        DOUBLE PRECISION
     $          XX (XSIZE)
        LOGICAL TRANSA
        
C=== UMD2CO ============================================================
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
C  Convert input matrix (II,XX,n,nz) into from triplet form to column-
C  oriented form, optionally transposing the input matrix (single and
C  double precision versions only).  Remove invalid entries and
C  duplicate entries.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       n:              size of the matrix
C       nz:             number of nonzeros in the input matrix

C       transa:         if true then transpose input matrix

C       XX (1..nz):     values of triplet form
C       xsize:          size of XX, must be >= 2*nz
C       II (1..2*nz):   row and col indices of triplet form
C       isize:          size of II, must be >= max (2*nz,n+1) + nz
C       Icntl:          integer control parameters
C       who:            who called UMD2CO, 1: UMD2FA, 2: UMD2RF
C
C       II must be at least of size (nz + max (2*nz, n+1))
C       XX must be at least of size (nz + max (  nz, n+1))
C
C       input triplet matrix:

C          if (transa) is false:
C               II (p)          i, row index, for p = 1..nz
C               II (nz+p)       j, col index
C               XX (p)          a_ij
C          if (transa) is true:
C               II (p)          j, col index, for p = 1..nz
C               II (nz+p)       i, row index
C               XX (p)          a_ij

C=======================================================================
C  WORKSPACE:
C=======================================================================
C
C       W (1..n)

C=======================================================================
C  OUTPUT: 
C=======================================================================
C
C       nz:             number of nonzeros in the output matrix,
C                       after removing invalid entries, and summing up
C                       duplicate entries
C       II (n+2..nz+n+1): row indices in column-form
C       XX (1..nz):     values in column-form.
C       Info (1):       error flag
C       Info (3):       invalid entries
C       Info (2):       duplicate entries
C       Info (5):       remaining valid entries
C       Info (6):       remaining valid entries
C       Info (7):       0
C       Wp (1..n+1)     column pointers for column form
C       II (1..n+1)     column pointers for column form

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutines:  UMD2FA, UMD2RF
C       subroutines called:     UMD2ER, UMD2P2
C       functions called:       MAX
        INTRINSIC MAX

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER ROW, COL, PDEST, P, NZ1, PCOL, IP, XP, IO, PRL, NINVLD,
     $          NDUPL, I
        LOGICAL PR3

C  row:     row index
C  col:     column index
C  pdest:   location of an entry in the column-form, for dupl. removal
C  p:       pointer
C  nz1:     number of entries after removing invalid or duplic. entries
C  pcol:    column col starts here after duplicates removed
C  ip:      column-form copy of matrix placed in II (ip...ip+nz-1)
C  xp:      column-form copy of matrix placed in XX (xp...xp+nz-1)
C  ninvld:  number of invalid entries
C  ndupl:   number of duplicate entries
C  i:       a row index if transa is true, a column index otherwise
C  io:      I/O unit for warning messages (for invalid or dupl. entries)
C  prl:     printing level
C  pr3:     true if printing invalid and duplicate entries

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

C-----------------------------------------------------------------------
C  get arguments and check memory sizes
C-----------------------------------------------------------------------

        IO = ICNTL (2)
        PRL = ICNTL (3)
        PR3 = PRL .GE. 3 .AND. IO .GE. 0

C-----------------------------------------------------------------------
C  count nonzeros in columns and check for invalid entries
C-----------------------------------------------------------------------

        NINVLD = 0
        NDUPL = 0
        DO 10 COL = 1, N 
           W (COL) = 0
10      CONTINUE 
        NZ1 = NZ
        DO 20 P = NZ, 1, -1 
           ROW = II (P)
           COL = II (NZ+P)
           IF (ROW.LT.1.OR.ROW.GT.N.OR.COL.LT.1.OR.COL.GT.N) THEN 
C             this entry is invalid - delete it
              IF (PR3) THEN 
C                print the offending entry on the diagnostic I/O unit
                 CALL UMD2P2 (WHO, 99, ROW, COL, XX(P), IO)
              ENDIF 
              II (P)    = II (NZ1)
              II (NZ+P) = II (NZ+NZ1)
              XX (P)    = XX (NZ1)
              NZ1 = NZ1 - 1
           ELSE 

              IF (TRANSA) THEN 
C                factorizing A transpose
                 W (ROW) = W (ROW) + 1
              ELSE 
C                factorizing A
                 W (COL) = W (COL) + 1
              ENDIF 

           ENDIF 
20      CONTINUE 
        NINVLD = NZ - NZ1
        IF (NINVLD .NE. 0) THEN 
C          invalid entries found - set warning flag and continue
           CALL UMD2ER (WHO, ICNTL, INFO, 1, NINVLD)
        ENDIF 

C-----------------------------------------------------------------------
C  convert triplet form to column-form
C-----------------------------------------------------------------------

        WP (1) = 1
        DO 30 I = 1, N 
           WP (I+1) = WP (I) + W (I)
30      CONTINUE 
        DO 40 I = 1, N 
           W (I) = WP (I)
40      CONTINUE 

C       ----------------------------------------------------------------
C       construct column-form in II (2*nz+1..3*nz) and XX (nz+1..2*nz)
C       ----------------------------------------------------------------

        IP = MAX (2*NZ, N+1)
        XP = NZ

        IF (TRANSA) THEN 
           DO 50 P = 1, NZ1 
              ROW = II (P)
              COL = II (NZ+P)
              II (IP + W (ROW)) = COL
              XX (XP + W (ROW)) = XX (P)
              W (ROW) = W (ROW) + 1
50         CONTINUE 
        ELSE 

           DO 60 P = 1, NZ1 
              ROW = II (P)
              COL = II (NZ+P)
              II (IP + W (COL)) = ROW
              XX (XP + W (COL)) = XX (P)
              W (COL) = W (COL) + 1
60         CONTINUE 

        ENDIF 

C       ----------------------------------------------------------------
C       shift the matrix back to II (n+2..nz+n+1) and XX (n+2..nz+n+1)
C       ----------------------------------------------------------------

        NZ = NZ1
CFPP$ NODEPCHK L
        DO 70 P = 1, NZ 
           II (N+1+P) = II (IP+P)
           XX (P) = XX (XP+P)
70      CONTINUE 

C-----------------------------------------------------------------------
C  remove duplicate entries by adding them up
C-----------------------------------------------------------------------

        DO 80 ROW = 1, N 
           W (ROW) = 0
80      CONTINUE 
        PDEST = 1
        DO 100 COL = 1, N 
           PCOL = PDEST
           DO 90 P = WP (COL), WP (COL+1)-1 
              ROW = II (N+1+P)
              IF (W (ROW) .GE. PCOL) THEN 
C                this is a duplicate entry
                 XX (W (ROW)) = XX (W (ROW)) + XX (P)
                 IF (PR3) THEN 
C                   print the duplicate entry on the diagnostic I/O
C                   unit.  The row and column indices printed reflect
C                   the input matrix.

                    IF (TRANSA) THEN 
                       CALL UMD2P2 (WHO, 98, COL, ROW, XX (P), IO)
                    ELSE 
                       CALL UMD2P2 (WHO, 98, ROW, COL, XX (P), IO)
                    ENDIF 

                 ENDIF 
              ELSE 
C                this is a new entry, store and record where it is
                 W (ROW) = PDEST
                 IF (PDEST .NE. P) THEN 
                    II (N+1+PDEST) = ROW
                    XX (PDEST) = XX (P)
                 ENDIF 
                 PDEST = PDEST + 1
              ENDIF 
90         CONTINUE 
           WP (COL) = PCOL
100     CONTINUE 
        WP (N+1) = PDEST
        NZ1 = PDEST - 1
        NDUPL = NZ - NZ1
        IF (NDUPL .NE. 0) THEN 
C          duplicate entries found - set warning flag and continue
           CALL UMD2ER (WHO, ICNTL, INFO, 2, NDUPL)
        ENDIF 
        NZ = NZ1

C-----------------------------------------------------------------------
C  save column pointers in II (1..n+1)
C-----------------------------------------------------------------------

        DO 110 COL = 1, N+1 
           II (COL) = WP (COL)
110     CONTINUE 

        INFO (2) = NDUPL
        INFO (3) = NINVLD
        INFO (5) = NZ
        INFO (6) = NZ
        INFO (7) = 0
        IF (NZ .EQ. 0) THEN 
C          set error flag if all entries are invalid
           CALL UMD2ER (WHO, ICNTL, INFO, -2, -1)
        ENDIF 
        RETURN
        END 
