        SUBROUTINE UMD2FG (XX, XSIZE, XHEAD, XTAIL, XUSE,
     $          II, ISIZE, IHEAD, ITAIL, IUSE,
     $          CP, RP, DN, N, ICNTL, WIR, WIC, WR, WC,
     $          FFXP, FFSIZE, WXP, FFDIMC, DOSLOT,
     $          PFREE, XFREE, MHEAD, MTAIL, SLOTS)
        INTEGER N, DN, ISIZE, II (ISIZE), IHEAD, ITAIL, RP (N+DN),
     $          CP (N+1), ICNTL (20), WIR (N), WIC (N), XSIZE, XUSE,
     $          IUSE, XHEAD, XTAIL, FFXP, FFSIZE, WXP,
     $          FFDIMC, WR (N+DN), WC (N+DN), PFREE, XFREE, MHEAD,
     $          MTAIL, SLOTS
        LOGICAL DOSLOT
        DOUBLE PRECISION
     $          XX (XSIZE)
        
C=== UMD2FG ============================================================
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
C  Garbage collection for UMD2F2.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       II/XX:          integer/real workspace, containing matrix being
C                       factorized and partially-computed LU factors
C       isize:          size of II
C       xsize:          size of XX
C       xhead:          XX (1..xhead) is in use (matrix, frontal mtc's)
C       xtail:          XX (xtail..xsize) is in use (LU factors)
C       xuse:           memory usage in Value
C       ihead:          II (1..ihead) is in use (matrix, frontal mtc's)
C       itail:          II (itail..isize) is in use (LU factors)
C       iuse:           memory usage in Index
C       Cp (1..n+1):    pointers to columns
C       Rp (1..n+dn):   pointers to rows, frontal matrices, and LU
C                       arrowheads
C       dn:             number of dense columns
C       n:              order of matrix
C       Icntl:          integer control parameters, see UMD21I
C       Wr (1..n):      see UMD2F2
C       Wc (1..n):      see UMD2F2
C       ffxp:           pointer to current contribution block
C       ffsize:         size of current contribution block
C       mhead:          pointer to first block in memory list
C       mtail:          pointer to last block in memory list
C       doslot:         true if adding slots
C       if doslot:
C           WiR (1..n)  if WiR (row) >= 0 then add (or keep) an extra
C                       slot in the row's element list
C           WiC (1..n)  if WiR (col) >= 0 then add (or keep) an extra
C                       slot in the col's element list

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       II/XX:          external fragmentation is removed at head 
C       xhead:          XX (1..xhead) is in use, reduced in size
C       xuse:           memory usage in Value, reduced
C       ihead:          II (1..ihead) is in use, reduced in size
C       iuse:           memory usage in Index, reduced
C       pfree:          pointer to free block in memory list, set to 0
C       xfree:          size of free block in XX, set to -1
C       mhead:          pointer to first block in memory list
C       mtail:          pointer to last block in memory list
C       ffxp            current working array has been shifted
C       wxp             current work vector has been shifted

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   UMD2F2

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER WHAT, FSIZ, ROW, COL, P, IDP, XDP, I, E, EP, FDIMC,
     $          LUDEGR, LUDEGC, J, PC, CELN, CLEN, RELN, RLEN,
     $          CSIZ1, CSIZ2, RSIZ1, RSIZ2, FLUIP, CXP, FXP, RDEG,
     $          CDEG, CSCAL, RSCAL, FSCAL
        PARAMETER (CSCAL = 9, RSCAL = 2, FSCAL = 7)
        LOGICAL SLOT

C  Compression:
C  ------------
C  what:    what this block of memory is (a row, column, etc.)
C  idp:     int. destination pointer, current block moved to II (idp...)
C  xdp:     real destination pointer, current block moved to XX (xdp...)
C  slot:    true if adding, or keeping, a size-2 hole in an element list
C
C  Columns:
C  --------
C  cscal:   = 9, the number of scalars in column data structure
C  celn:    number of (e,f) tuples in element list of a column
C  clen:    number of unassembled original entries in a column
C  cdeg:    degree of a column (number of entries, including fill-in)
C  cxp:     a column is in XX (cxp...) prior to compression
C  pc:      column is in II (pc ...) prior to compression
C  csiz1:   size of a column in II, prior to compression
C  csiz2:   size of a column in II, after compression
C  col:     column index
C
C  Rows:
C  -----
C  rscal:   = 2, the number of scalars in row data structure
C  reln:    number of (e,f) tuples in element list of a row
C  rlen:    number of unassembled original entries in a row
C  rsiz1:   size of a row in II, prior to compression
C  rsiz2:   size of a row in II, after compression
C  rdeg:    degree of a row (number of entries, including fill-in)
C  row:     row index
C
C  Frontal matrices:
C  -----------------
C  fscal:   = 7, the number of scalars in element data structure
C  fluip:   element is in II (fluip...) prior to compression
C  fxp:     a frontal matrix is in XX (fxp...) prior to compression
C  e:       an element
C  fdimc:   column dimension (number of rows) of a frontal matrix
C  ludegr:  row degree (number of columns) of a contribution block
C  ludegc:  column degree (number of rows) of a contribution block
C  fsiz:    size of an artificial frontal matrix
C  ep:      an artificial frontal matrix is in II (ep ...) prior to comp
C
C  Other:
C  ------
C  p:       pointer
C  i:       general loop index
C  j:       general loop index

C=======================================================================
C  EXECUTABLE STATMENTS:
C=======================================================================

        SLOTS = 0

C-----------------------------------------------------------------------
C   prepare the non-pivotal rows/cols and unassembled elements
C-----------------------------------------------------------------------

C       place the size of each block of memory at the beginning,
C       and mark the 2nd entry in each block with what it is

C       ----------------------------------------------------------------
C       mark the columns
C       ----------------------------------------------------------------

CFPP$ NODEPCHK L
        DO 10 COL = 1, N 
           PC = CP (COL)
           IF (PC .NE. 0) THEN 
C             this is a non-pivotal, non-null column
              CDEG = II (PC+1)
              CP (COL) = CDEG
              II (PC+1) = COL+N
           ENDIF 
10      CONTINUE 

C       ----------------------------------------------------------------
C       mark the rows and frontal matrices
C       ----------------------------------------------------------------

CFPP$ NODEPCHK L
        DO 20 ROW = 1, N 
           P = RP (ROW)
           RLEN = WC (ROW)
           IF (P .EQ. 0) THEN 
C             a pivotal row
              CONTINUE
           ELSE IF (RLEN .GE. 0 .AND. RLEN .LE. N) THEN 
C             this is a non-pivotal, non-null row
              RDEG = II (P+1)
              RP (ROW) = RDEG
              II (P+1) = ROW+2*N
           ELSE IF (WR (ROW) .EQ. -(N+DN+2)) THEN 
C             a pivotal row, and an assembled element
              CONTINUE
           ELSE 
C             this is an unassembled frontal matrix
C             the size is implicitly fscal
              FDIMC = II (P+1)
              RP (ROW) = FDIMC
              II (P+1) = ROW
           ENDIF 
20      CONTINUE 

C       ----------------------------------------------------------------
C       mark the artificial frontal matrices
C       ----------------------------------------------------------------

CFPP$ NODEPCHK L
        DO 30 E = N+1, N+DN 
           EP = RP (E)
           IF (EP .NE. 0) THEN 
C             this is an unassembled artificial frontal matrix
C             the size is II (ep+1) + cscal
              FDIMC = II (EP+1)
              RP (E) = FDIMC
              II (EP+1) = E+2*N
           ENDIF 
30      CONTINUE 

C-----------------------------------------------------------------------
C  scan the link list and compress the reals
C-----------------------------------------------------------------------

        XDP = 1
        P = MHEAD
C       while (p .ne. 0) do
40      CONTINUE
        IF (P .NE. 0) THEN 

           WHAT = II (P+1)

C          -------------------------------------------------------------
           IF (WHAT .GT. 3*N) THEN 
C          -------------------------------------------------------------

C             this is an unassembled artificial frontal matrix
              E = WHAT - 2*N
              FXP = II (P+2)
              II (P+2) = XDP
CFPP$ NODEPCHK L
              DO 50 J = 0, RP (E) - 1 
                 XX (XDP+J) = XX (FXP+J)
50            CONTINUE 
              XDP = XDP + RP (E)

C          -------------------------------------------------------------
           ELSE IF (WHAT .EQ. -1 .OR. II (P+6) .EQ. 0) THEN 
C          -------------------------------------------------------------

C             this is a real hole - delete it from the link list
              IF (II (P+4) .NE. 0) THEN 
                 II (II (P+4)+3) = II (P+3)
              ELSE 
                 MHEAD = II (P+3)
              ENDIF 
              IF (II (P+3) .NE. 0) THEN 
                 II (II (P+3)+4) = II (P+4)
              ELSE 
                 MTAIL = II (P+4)
              ENDIF 

C          -------------------------------------------------------------
           ELSE IF (WHAT .LE. N) THEN 
C          -------------------------------------------------------------

C             this is an unassembled frontal matrix
              E = WHAT
              FXP = II (P+2)
              II (P+2) = XDP
              FLUIP = II (P)
              LUDEGR = II (FLUIP+2)
              LUDEGC = II (FLUIP+3)
              FDIMC = RP (E)
              IF (FDIMC .EQ. LUDEGC) THEN 
C                contribution block is already compressed
CFPP$ NODEPCHK L
                 DO 60 I = 0, (LUDEGR * LUDEGC) - 1 
                    XX (XDP+I) = XX (FXP+I)
60               CONTINUE 
              ELSE 
C                contribution block is not compressed
C                compress XX (fxp..) to XX (xdp..xdp+(ludegr*ludegc)-1)
                 DO 80 J = 0, LUDEGR - 1 
CFPP$ NODEPCHK L
                    DO 70 I = 0, LUDEGC - 1 
                       XX (XDP + J*LUDEGC + I) = XX (FXP + J*FDIMC + I)
70                  CONTINUE 
80               CONTINUE 
                 RP (E) = LUDEGC
              ENDIF 
              XDP = XDP + LUDEGR*LUDEGC

C          -------------------------------------------------------------
           ELSE IF (WHAT .LE. 2*N) THEN 
C          -------------------------------------------------------------

C             this is a column
              CXP = II (P+2)
              II (P+2) = XDP
              CLEN = II (P+6)
CFPP$ NODEPCHK L
              DO 90 J = 0, CLEN - 1 
                 XX (XDP+J) = XX (CXP+J)
90            CONTINUE 
              XDP = XDP + CLEN

C          -------------------------------------------------------------
           ENDIF 
C          -------------------------------------------------------------

C          -------------------------------------------------------------
C          get the next item in the link list
C          -------------------------------------------------------------

           P = II (P+3)

C       end while:
        GOTO 40
        ENDIF 

        PFREE = 0
        XFREE = -1

C       ----------------------------------------------------------------
C       shift the current working array (if it exists)
C       ----------------------------------------------------------------

        IF (FFXP .NE. 0) THEN 
CFPP$ NODEPCHK L
           DO 100 I = 0, FFSIZE - 1 
              XX (XDP+I) = XX (FFXP+I)
100        CONTINUE 
           FFXP = XDP
           XDP = XDP + FFSIZE
        ENDIF 

C       ----------------------------------------------------------------
C       shift the current work vector (if it exists)
C       ----------------------------------------------------------------

        IF (WXP .NE. 0) THEN 
           WXP = XDP
           XDP = XDP + FFDIMC
        ENDIF 

C-----------------------------------------------------------------------
C  scan from the top of integer memory (1) to bottom (ihead) and
C  compress the integers
C-----------------------------------------------------------------------

        P = 1
        IDP = P
C       while (p .lt. ihead) do:
110     CONTINUE
        IF (P .LT. IHEAD) THEN 

           WHAT = II (P+1)

C          -------------------------------------------------------------
           IF (WHAT .GT. 3*N) THEN 
C          -------------------------------------------------------------

C             this is an unassembled artificial frontal matrix
              E = WHAT - 2*N
              FSIZ = RP (E) + CSCAL
              II (P+1) = RP (E)
              RP (E) = IDP
CFPP$ NODEPCHK L
              DO 120 I = 0, FSIZ - 1 
                 II (IDP+I) = II (P+I)
120           CONTINUE 
C             shift pointers in the link list
              IF (II (IDP+4) .NE. 0) THEN 
                 II (II (IDP+4)+3) = IDP
              ELSE 
                 MHEAD = IDP
              ENDIF 
              IF (II (IDP+3) .NE. 0) THEN 
                 II (II (IDP+3)+4) = IDP
              ELSE 
                 MTAIL = IDP
              ENDIF 
              P = P + FSIZ
              IDP = IDP + FSIZ

C          -------------------------------------------------------------
           ELSE IF (WHAT .EQ. -1) THEN 
C          -------------------------------------------------------------

C             this is a integer hole
              P = P + II (P)

C          -------------------------------------------------------------
           ELSE IF (WHAT .GE. 1 .AND. WHAT .LE. N) THEN 
C          -------------------------------------------------------------

C             this is an unassembled frontal matrix (fscal integers)
              E = WHAT
              FDIMC = RP (E)
              II (P+1) = FDIMC
              RP (E) = IDP
CFPP$ NODEPCHK L
              DO 130 I = 0, FSCAL - 1 
                 II (IDP+I) = II (P+I)
130           CONTINUE 
C             shift pointers in the link list
              IF (II (IDP+4) .NE. 0) THEN 
                 II (II (IDP+4)+3) = IDP
              ELSE 
                 MHEAD = IDP
              ENDIF 
              IF (II (IDP+3) .NE. 0) THEN 
                 II (II (IDP+3)+4) = IDP
              ELSE 
                 MTAIL = IDP
              ENDIF 
              P = P + FSCAL
              IDP = IDP + FSCAL

C          -------------------------------------------------------------
           ELSE IF (WHAT .LE. 2*N) THEN 
C          -------------------------------------------------------------

C             this is a non-pivotal column
              CSIZ1 = II (P)
              COL = WHAT - N
              CELN = II (P+5)
              CLEN = II (P+6)
              CSIZ2 = 2*CELN + CLEN + CSCAL
              SLOT = DOSLOT .AND. WIC (COL) .GE. 0 .AND. P .GE. IDP+2
              IF (SLOT) THEN 
C                keep (or make) one extra slot for element list growth
                 CSIZ2 = CSIZ2 + 2
                 SLOTS = SLOTS + 2
              ENDIF 
              CDEG = CP (COL)
              II (P+1) = CDEG
              CP (COL) = IDP
              II (P) = CSIZ2
C             copy the cscal scalars and the celn (e,f) tuples
CFPP$ NODEPCHK L
              DO 140 I = 0, CSCAL + 2*CELN - 1 
                 II (IDP+I) = II (P+I)
140           CONTINUE 
              IF (CLEN .GT. 0) THEN 
C                shift pointers in the link list
                 IF (II (IDP+4) .NE. 0) THEN 
                    II (II (IDP+4)+3) = IDP
                 ELSE 
                    MHEAD = IDP
                 ENDIF 
                 IF (II (IDP+3) .NE. 0) THEN 
                    II (II (IDP+3)+4) = IDP
                 ELSE 
                    MTAIL = IDP
                 ENDIF 
              ENDIF 
              P = P + CSIZ1 - CLEN
              IDP = IDP + CSCAL + 2*CELN
              IF (SLOT) THEN 
C                skip past the slot
                 IDP = IDP + 2
              ENDIF 
C             copy the clen original row indices
CFPP$ NODEPCHK L
              DO 150 I = 0, CLEN - 1 
                 II (IDP+I) = II (P+I)
150           CONTINUE 
              P = P + CLEN
              IDP = IDP + CLEN

C          -------------------------------------------------------------
           ELSE 
C          -------------------------------------------------------------

C             this is a non-pivotal row
              RSIZ1 = II (P)
              ROW = WHAT - 2*N
              RELN = WR (ROW)
              RLEN = WC (ROW)
              RSIZ2 = 2*RELN + RLEN + RSCAL
              SLOT = DOSLOT .AND. WIR (ROW) .GE. 0 .AND. P .GE. IDP+2
              IF (SLOT) THEN 
C                keep (or make) one extra slot for element list growth
                 RSIZ2 = RSIZ2 + 2
                 SLOTS = SLOTS + 2
              ENDIF 
              RDEG = RP (ROW)
              II (P+1) = RDEG
              RP (ROW) = IDP
              II (P) = RSIZ2
C             copy the rscal scalars, and the reln (e,f) tuples
CFPP$ NODEPCHK L
              DO 160 I = 0, RSCAL + 2*RELN - 1 
                 II (IDP+I) = II (P+I)
160           CONTINUE 
              P = P + RSIZ1 - RLEN
              IDP = IDP + RSCAL + 2*RELN
              IF (SLOT) THEN 
C                skip past the slot
                 IDP = IDP + 2
              ENDIF 
C             copy the rlen original column indices
CFPP$ NODEPCHK L
              DO 170 I = 0, RLEN - 1 
                 II (IDP+I) = II (P+I)
170           CONTINUE 
              P = P + RLEN
              IDP = IDP + RLEN

C          -------------------------------------------------------------
           ENDIF 
C          -------------------------------------------------------------

C          -------------------------------------------------------------
C          move to the next block
C          -------------------------------------------------------------

C       end while:
        GOTO 110
        ENDIF 

C-----------------------------------------------------------------------
C  deallocate the unused space
C-----------------------------------------------------------------------

        IUSE = IUSE - (IHEAD - IDP)
        IHEAD = IDP
        XUSE = XUSE - (XHEAD - XDP)
        XHEAD = XDP
        RETURN
        END 
