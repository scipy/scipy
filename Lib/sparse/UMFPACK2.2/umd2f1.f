        SUBROUTINE UMD2F1 (CP, N, CPERM, RPERM, NZOFF,
     $          ITAIL, XTAIL, XX, XSIZE, XUSE, II, ISIZE, IUSE,
     $          ICNTL, CNTL, INFO, RINFO, NBLKS,
     $          AP, AI, AX, PRESRV, K1, AN, ANZ, PR, KEEP,
     $          RMAX, CMAX, TOTNLU, XRMAX, XRUSE, IOUT, XOUT)
        INTEGER XSIZE, ISIZE, N, ICNTL (20), INFO (40), XUSE, IUSE,
     $          ITAIL, XTAIL, II (ISIZE), CP (N+1), CPERM (N), NZOFF,
     $          AN, ANZ, RPERM (N), AI (ANZ), AP (AN+1), K1, PR (AN),
     $          NBLKS, KEEP (20), RMAX, CMAX, TOTNLU, XRMAX, XRUSE
        LOGICAL PRESRV, IOUT, XOUT
        DOUBLE PRECISION
     $          XX (XSIZE), AX (ANZ)
        DOUBLE PRECISION
     $          CNTL (10), RINFO (20)
        
C=== UMD2F1 ============================================================
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
C  UMD2F1 factorizes the n-by-n column-form matrix at the head of II/XX
C  or in Ap/Ai/Ax, and places its LU factors at the tail of II/XX.  The
C  input matrix overwritten if it is located in II/XX on input.  If
C  block-triangular-form (BTF) is in use, this routine factorizes a
C  single diagonal block.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       n:              order of matrix (or order of diagonal block
C                       if BTF is in use).
C       Cp (1..n+1):    column pointers for input matrix
C       nblks:          number of diagonal blocks in BTF form
C       isize:          size of II
C       xsize:          size of XX
C       k1:             first index of this matrix (1 if BTF not used)
C       Icntl:          integer control parameters, see UMD21I
C       Cntl:           real control parameters, see UMD21I
C       Keep (6..8):    integer control parameters, see UMD21I
C       iuse:           memory usage in Index
C       xuse:           memory usage in Value
C       rmax:           maximum ludegr seen so far (see UMD2F2 for info)
C       cmax:           maximum ludegc seen so far (see UMD2F2 for info)
C       totnlu:         total number of LU arrowheads constructed so far
C
C       if nblks>1 then:
C          Cperm (1..n):        col permutation to BTF
C          Rperm (1..n):        row permutation to BTF
C       else
C          Cperm (1..n):        undefined on input
C          Rperm (1..n):        undefined on input
C
C
C       presrv:         if true then input matrix is preserved
C
C       if presrv is true then:
C           an:                 order of preserved matrix (all blocks)
C           anz:                entries in preserved matrix
C           Ap (1..an+1):       column pointers for preserved matrix
C           Ai (1..anz):        row indices of preserved matrix
C           Ax (1..anz):        values of preserved matrix
C                               The preserved matrix is not in BTF form;
C                               it is in the orginal order.
C           if nblks > 1:
C               Pr (1..n):      inverse row permutations to BTF form
C               nzoff           entries in off-diagonal blocks
C                               seen so far
C           else
C               Pr (1..n):      undefined on input
C
C           II (1..isize):      undefined on input
C           XX (1..xsize):      undefined on input
C           Cp (1..n+1):        undefined on input
C
C       else, if presrv is false:
C           an:                         1
C           anz:                        1
C           II (1..Cp (1) - 1):         unused
C           II (Cp (1) ... Cp (n+1)-1): row indices of matrix to factor,
C                                       will be overwritten on output
C           II (Cp (n+1) ... isize):    unused on input
C
C           XX (1..Cp (1) - 1):         unused
C           XX (Cp (1) ... Cp (n+1)-1): values of matrix to factorize,
C                                       will be overwritten on output
C           XX (Cp (n+1) ... xsize):    unused on input
C                       If BTF is in use, then II and XX contain a
C                       single diagonal block.

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       XX (xtail ... xsize), xtail,  II (itail ... isize), itail:
C
C                       The LU factors of a single diagonal block.
C                       See UMD2F2 for a description.
C
C       II (cp1 ... itail-1):   undefined on output
C       XX (cp1 ... xtail-1):   undefined on output,
C                       where cp1 is equal to the value of Cp (1)
C                       if presrv is false, or cp1 = 1 if presrv is
C                       true.
C
C       Info:           integer informational output, see UMD2FA
C       Rinfo:          real informational output, see UMD2FA
C       Cperm (1..n):   the final col permutations, including BTF
C       Rperm (1..n):   the final row permutations, including BTF
C
C       iuse:           memory usage in Index
C       xuse:           memory usage in Value
C       rmax:           maximum ludegr seen so far (see UMD2F2 for info)
C       cmax:           maximum ludegc seen so far (see UMD2F2 for info)
C       totnlu:         total number of LU arrowheads constructed so far
C
C       if nblks>1 and presrv:
C           nzoff       entries in off-diagonal blocks seen so far
C
C       iout:           true if ran out of integer memory in UMD2F1
C       xout:           true if ran out of real memory in UMD2F1

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   UMD2F0
C       subroutines called:     UMD2F2
C       functions called:       MAX, SQRT
        INTRINSIC MAX, SQRT

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER CP1, PC, PEND, PCOL, CDEG, COL, CSIZ, NZ, XP, IP, IS, P,
     $          DN, DSIZ, WRKSIZ, I, CLEN, ROW, CSCAL
        PARAMETER (CSCAL = 9)

C  Original and expanded column-form:
C  ----------------------------------
C  cp1:     = value of Cp (1) on input
C  pc:      pointer to integer part of expanded column-form matrix
C  pend:    column col ends here in the input column-form matrix
C  pcol:    column col starts here in the input column-form matrix
C  cdeg:    degree (number of entries) in a column
C  clen:    number of original entries in a column (= degree, here,
C           but decreases in UMD2F2)
C  csiz:    size of the integer part of an expanded column (cdeg+cscal)
C  cscal:   = 9, the number of scalars in column data structure
C  nz:      number of entries in the diagonal block being factorized
C  xp:      pointer to real part of the expanded column
C  ip:      pointer to integer part of the expanded column
C
C  Memory usage:
C  -------------
C  wrksiz:  size of integer workspace needed by UMD2F2
C
C  "Dense" columns: (converted to a prior, or artificial, frontal mat.):
C  ----------------
C  dn:      number of "dense" columns
C  dsiz:    a column is "dense" if it has more than dsiz entries
C
C  Other:
C  ------
C  row:     row index
C  col:     a column index
C  i:       loop index
C  p:       pointer

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

        IOUT = .FALSE.
        XOUT = .FALSE.

C-----------------------------------------------------------------------
C  count "dense" columns (they are treated as a priori frontal matrices)
C-----------------------------------------------------------------------

C       a column is "dense" if it has more than dsiz entries
        DSIZ = MAX (0, KEEP (7),
     $                 KEEP (8) * INT (SQRT (REAL (N))))
        DN = 0
        IF (PRESRV) THEN 
           IF (NBLKS .EQ. 1) THEN 
              DO 10 COL = 1, N 
                 IF (AP (COL+1) - AP (COL) .GT. DSIZ) THEN 
C                   this is a "dense" column
                    DN = DN + 1
                 ENDIF 
10            CONTINUE 
           ELSE 
              DO 40 COL = 1, N 
C                if col might be dense, check more carefully:
                 CDEG = AP (CPERM (COL) + 1)- AP (CPERM (COL))
                 IF (CDEG .GT. DSIZ) THEN 
                    CDEG = 0
                    DO 20 P = AP (CPERM (COL)), AP (CPERM (COL) + 1) -1
                       ROW = PR (AI (P))
                       IF (ROW .GE. K1) THEN 
                          CDEG = CDEG + 1
                          IF (CDEG .GT. DSIZ) THEN 
C                            this is a "dense" column, exit out of loop
                             DN = DN + 1
                             GO TO 30
                          ENDIF 
                       ENDIF 
20                  CONTINUE 
C                   loop exit label:
30                  CONTINUE
                 ENDIF 
40            CONTINUE 
           ENDIF 
        ELSE 
           DO 50 COL = 1, N 
              IF (CP (COL+1) - CP (COL) .GT. DSIZ) THEN 
C                this is a "dense" column
                 DN = DN + 1
              ENDIF 
50         CONTINUE 
        ENDIF 

C-----------------------------------------------------------------------
C  get size of workspaces to allocate from II
C-----------------------------------------------------------------------

C       workspaces: WiR (n), WiC (n), WpR (n), WpC (n),
C       Wm (n), Head (n), Rp (n+dn), Wc (n+dn), Wr (n+dn), Wj (n)
        IF (NBLKS .EQ. 1) THEN 
C          Rperm (1..n) is used as WiR (1..n), and
C          Cperm (1..n) is used as WiC (1..n) in UMD2F2
           WRKSIZ = 8*N + 3*DN
        ELSE 
           WRKSIZ = 10*N + 3*DN
        ENDIF 

C-----------------------------------------------------------------------
C  construct the expanded column-form of the matrix or the diag. block
C-----------------------------------------------------------------------

        IF (PRESRV) THEN 

C          -------------------------------------------------------------
C          allocate space for wrksiz workspace and nz+cscal*n
C          integers and nz reals for the expanded column-form matrix.
C          -------------------------------------------------------------

           CP1 = 1
           XP = 1
           IP = 1 + WRKSIZ
           IF (NBLKS .EQ. 1) THEN 

C             ----------------------------------------------------------
C             construct copy of entire matrix
C             ----------------------------------------------------------

              NZ = ANZ
              IS = NZ + WRKSIZ + CSCAL*N
              IUSE = IUSE + IS
              XUSE = XUSE + NZ
              INFO (18) = MAX (INFO (18), IUSE)
              INFO (19) = MAX (INFO (19), IUSE)
              INFO (20) = MAX (INFO (20), XUSE)
              INFO (21) = MAX (INFO (21), XUSE)
              IOUT = IS .GT. ISIZE
              XOUT = NZ .GT. XSIZE
              IF (IOUT .OR. XOUT) THEN 
C                error return, if not enough integer and/or real memory:
                 GO TO 9000
              ENDIF 

              PC = IP
              DO 70 COL = 1, N 
                 CP (COL) = PC - WRKSIZ
                 CDEG = AP (COL+1) - AP (COL)
                 CLEN = CDEG
                 CSIZ = CDEG + CSCAL
                 II (PC) = CSIZ
                 II (PC+1) = CDEG
                 II (PC+5) = 0
                 II (PC+6) = CLEN
                 II (PC+7) = 0
                 II (PC+8) = 0
                 II (PC+2) = XP
                 XP = XP + CDEG
                 PC = PC + CSCAL
                 P = AP (COL)
                 DO 60 I = 0, CDEG - 1 
                    II (PC + I) = AI (P + I)
60               CONTINUE 
                 PC = PC + CDEG
70            CONTINUE 
              DO 80 P = 1, NZ 
                 XX (P) = AX (P)
80            CONTINUE 

           ELSE 

C             ----------------------------------------------------------
C             construct copy of a single block in BTF form
C             ----------------------------------------------------------

C             check for memory usage during construction of block
              DO 100 COL = 1, N 
                 PC = IP
                 CP (COL) = PC - WRKSIZ
                 IP = IP + CSCAL
                 IOUT = IP .GT. ISIZE
                 IF (IOUT) THEN 
C                   error return, if not enough integer memory:
                    GO TO 9000
                 ENDIF 
                 II (PC+2) = XP
                 CDEG = IP
                 DO 90 P = AP (CPERM (COL)), AP (CPERM (COL)+1)-1 
                    ROW = PR (AI (P))
                    IF (ROW .GE. K1) THEN 
                       IOUT = IP .GT. ISIZE
                       XOUT = XP .GT. XSIZE
                       IF (IOUT .OR. XOUT) THEN 
C                         error return, if not enough memory
                          GO TO 9000
                       ENDIF 
                       II (IP) = ROW - K1 + 1
                       XX (XP) = AX (P)
                       IP = IP + 1
                       XP = XP + 1
                    ELSE 
C                      entry in off-diagonal part
                       NZOFF = NZOFF + 1
                    ENDIF 
90               CONTINUE 
                 CDEG = IP - CDEG
                 CLEN = CDEG
                 CSIZ = CDEG + CSCAL
                 II (PC) = CSIZ
                 II (PC+1) = CDEG
                 II (PC+5) = 0
                 II (PC+6) = CLEN
                 II (PC+7) = 0
                 II (PC+8) = 0
100           CONTINUE 

              NZ = XP - 1
              IS = NZ + WRKSIZ + CSCAL*N
              IUSE = IUSE + IS
              XUSE = XUSE + NZ
              INFO (18) = MAX (INFO (18), IUSE)
              INFO (19) = MAX (INFO (19), IUSE)
              INFO (20) = MAX (INFO (20), XUSE)
              INFO (21) = MAX (INFO (21), XUSE)

           ENDIF 

C          -------------------------------------------------------------
C          get memory usage for next call to UMD2RF
C          -------------------------------------------------------------

           XRUSE = XRUSE + NZ
           XRMAX = MAX (XRMAX, XRUSE)

        ELSE 

C          -------------------------------------------------------------
C          allocate space for wrksiz workspace and additional cscal*n
C          space for the expanded column-form of the matrix.
C          -------------------------------------------------------------

           CP1 = CP (1)
           NZ = CP (N+1) - CP1
           PC = CP1 + WRKSIZ + (NZ+CSCAL*N)
           IUSE = IUSE + WRKSIZ + CSCAL*N
           INFO (18) = MAX (INFO (18), IUSE)
           INFO (19) = MAX (INFO (19), IUSE)
           IOUT = PC .GT. ISIZE+1
           IF (IOUT) THEN 
C             error return, if not enough integer memory:
              GO TO 9000
           ENDIF 

C          -------------------------------------------------------------
C          expand the column form in place and make space for workspace
C          -------------------------------------------------------------

           XP = NZ + 1
           IP = NZ + CSCAL*N + 1
           PEND = CP (N+1)
           DO 120 COL = N, 1, -1 
              PCOL = CP (COL)
              DO 110 P = PEND-1, PCOL, -1 
                 PC = PC - 1
                 II (PC) = II (P)
110           CONTINUE 
              PC = PC - CSCAL
              CDEG = PEND - PCOL
              CLEN = CDEG
              PEND = PCOL
              CSIZ = CDEG + CSCAL
              IP = IP - CSIZ
              CP (COL) = IP
              II (PC) = CSIZ
              II (PC+1) = CDEG
              II (PC+5) = 0
              II (PC+6) = CLEN
              II (PC+7) = 0
              II (PC+8) = 0
              XP = XP - CDEG
              II (PC+2) = XP
120        CONTINUE 
        ENDIF 

C-----------------------------------------------------------------------
C  factorize the expanded column-form, with allocated workspaces
C-----------------------------------------------------------------------

        XP = CP1
        IP = CP1 + WRKSIZ

        IF (NBLKS .EQ. 1) THEN 

C          pass Rperm and Cperm as the WiR and WiC arrays in UMD2F2:
           CALL UMD2F2 (CP, NZ, N, 1, 1, 1, ITAIL, XTAIL,
     $          XX (XP), XSIZE-XP+1, II (IP), ISIZE-IP+1, ICNTL, CNTL,
     $          INFO, RINFO, .FALSE., IUSE, XUSE,
     $          RPERM, CPERM, II (CP1), II (CP1+N),
     $          II (CP1+2*N), II (CP1+3*N), II (CP1+4*N), II (CP1+5*N),
     $          II (CP1+6*N+DN), II (CP1+7*N+2*DN),
     $          DN, DSIZ, KEEP, RMAX, CMAX, TOTNLU, XRMAX, XRUSE)

        ELSE 

C          pass Cperm, Rperm, WiC and WiR as separate arrays, and
C          change Cperm and Rperm from the BTF permutations to the
C          final permutations (including BTF and numerical pivoting).
           CALL UMD2F2 (CP, NZ, N, N, CPERM, RPERM, ITAIL, XTAIL,
     $          XX (XP), XSIZE-XP+1, II (IP), ISIZE-IP+1, ICNTL, CNTL,
     $          INFO, RINFO, .TRUE., IUSE, XUSE,
     $          II (CP1), II (CP1+N), II (CP1+2*N), II (CP1+3*N),
     $          II (CP1+4*N), II (CP1+5*N), II (CP1+6*N), II (CP1+7*N),
     $          II (CP1+8*N+DN), II (CP1+9*N+2*DN),
     $          DN, DSIZ, KEEP, RMAX, CMAX, TOTNLU, XRMAX, XRUSE)

        ENDIF 

        IF (INFO (1) .LT. 0) THEN 
C          error return, if error occured in UMD2F2:
           RETURN
        ENDIF 

C-----------------------------------------------------------------------
C  adjust tail pointers, and save pointer to numerical part of LU
C-----------------------------------------------------------------------

C       head = cp1
        IUSE = IUSE - WRKSIZ
        ITAIL = ITAIL + IP - 1
        XTAIL = XTAIL + XP - 1
        II (ITAIL) = XTAIL
        RETURN

C=======================================================================
C  Error return
C=======================================================================

C       error return label:
9000    CONTINUE
        RETURN
        END 
