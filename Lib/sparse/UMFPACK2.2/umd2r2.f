        SUBROUTINE UMD2R2 (CP, NZ, N, XTAIL, XX, XSIZE, XUSE, ARI,
     $          CPERM, RPERM, ICNTL, CNTL, INFO, RINFO, MC, MR,
     $          WIR, WIC, WPR, WPC, WM, WJ, FRDIMC, FRXP, FRNEXT,
     $          FRPREV, NLU, LUP, LUI, NOUTSD, XRMAX)
        INTEGER XSIZE, ICNTL (20), INFO (40), CPERM (N), RPERM (N),
     $          XTAIL, NZ, N, ARI (NZ), CP (N+1), MR, MC, NOUTSD,
     $          WIR (N), WIC (N), WPR (MR), XRMAX, WPC (MC), WM (MC),
     $          NLU, FRDIMC (NLU+2), FRXP (NLU+2), XUSE, WJ (MC),
     $          FRNEXT (NLU+2), FRPREV (NLU+2), LUP (NLU), LUI (*)
        DOUBLE PRECISION
     $          XX (XSIZE)
        DOUBLE PRECISION
     $          CNTL (10), RINFO (20)
        
C=== UMD2R2 ============================================================
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
C  UMD2R2 refactorizes the n-by-n input matrix at the head of XX
C  (in arrowhead form) and places its LU factors at the tail of
C  XX.  The input matrix is overwritten.   No BTF information is
C  used in this routine.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       Cp (1..n+1):    column pointers of arrowhead form
C       n:              order of input matrix
C       nz:             entries in input matrix
C       xsize:          size of XX
C       Icntl:          integer control parameters, see UMD21I
C       Cntl:           real control parameters, see UMD21I
C
C       Ari (1..nz):            arrowhead format of A
C       XX (1..nz):             arrowhead format of A, see below
C       XX (nz+1..xsize):       undefined on input, used as workspace
C
C       nlu:            number of LU arrowheads
C       LUp (1..nlu):   pointers to LU arrowheads in LUi
C       LUi (1.. ):     LU arrowheads
C
C       xuse:           memory usage in Value
C
C       noutsd:         entries not in prior LU pattern
C
C       Cperm (1..n):   column permutation
C       Rperm (1..n):   row permutation

C=======================================================================
C  WORKSPACE:
C=======================================================================
C
C       WiR (1..n)
C       WiC (1..n)
C
C       WpR (1.. max ludegr)
C       WpC (1.. max ludegc)
C       Wm  (1.. max ludegc)
C       Wj  (1.. max ludegc)
C
C       FRdimc (1..nlu+2)
C       FRxp   (1..nlu+2)
C       FRnext (1..nlu+2)
C       FRprev (1..nlu+2)

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       LUi (1..):              LU arrowheads, modified luxp pointers
C       XX (1..xtail-1):        undefined on output
C       XX (xtail..xsize):      LU factors of this matrix, see below
C
C       Info:           integer informational output, see UMD2FA
C       Rinfo:          real informational output, see UMD2FA
C
C       xuse:           memory usage in Value
C
C       noutsd:         entries not in prior LU pattern, incremented

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   UMD2R0
C       subroutines called:     UMD2ER, UMD2P2, UMD2RG, DGEMV,
C                               DGEMM, DTRSV, DTRSM
C       functions called:       ABS, MAX
        INTRINSIC ABS, MAX

C=======================================================================
C  DESCRIPTION OF DATA STRUCTURES:
C=======================================================================

C-----------------------------------------------------------------------
C  Matrix being factorized:
C-----------------------------------------------------------------------
C
C  The input matrix is held in an arrowhead format.  For the kth pivot,
C  the nonzeros in the pivot row (A (k, k...n)) and pivot column
C  (A (k...n, k)) are stored in the kth arrowhead.  The kth arrowhead
C  is located in:
C       Ari (Cp (k+1) ... Cp (k)-1):    pattern
C       XX  (Cp (k+1) ... Cp (k)-1):    values
C
C  Suppose p is in the range Cp (k+1) to Cp (k)-1.  If Ari (p) is
C  greater than zero, then the entry is in row Ari (p), column k,
C  with value XX (p).  If Ari (p) is less than zero, then the entry is
C  in row k, column -Ari (p), with value XX (p).  The arrowheads are
C  stored in reverse order (arrowhead n, n-1, ... 2, 1) in Ari and XX.
C  Note that Cp (n+1) = 1 unless BTF is in use and the original matrix
C  is not preserved.   In all cases, the real part of the arrowhead
C  format (XX (Cp (n+1) ... Cp (1)-1)) is overwritten with the LU
C  factors.  The integer part (Ari (Cp (n+1) ... Cp (1)-1)) is not
C  overwritten, since UMD2R2 does not require dynamic allocation of
C  integer memory.

C-----------------------------------------------------------------------
C  Frontal matrices
C-----------------------------------------------------------------------
C
C   Each unassembled frontal matrix (element) is stored as follows:
C       total size: fscal integers, (fdimr*fdimc) reals
C
C       if e is an unassembled element, and not the current frontal
C       matrix:
C
C       fluip = LUp (e) pointer to LU arrowhead in II
C       fdimc = FRdimc (e)      column dimension of contribution block
C       fxp   = FRxp (e)        pointer to contribution block in XX
C       next  = FRnext (e)      pointer to next block in XX
C       prev  = FRprev (e)      pointer to previous block in XX
C       fdegr = abs (LUi (fluip+2))
C       fdegc = abs (LUi (fluip+2))
C       XX (fxp ... )
C               a 2-dimensional array, C (1..fdimc, 1..fdimr), where
C               fdimr = fdegr if the contribution block is compressed,
C               or fdimr = LUi (fluip+5) if not.  Note, however, that
C               fdimr is not needed.  The contribution block is stored
C               in C (1..fdegc, 1..fdegr) in the C (1..fdimc,...) array.
C
C               If memory is limited, garbage collection will occur.
C               In this case, the C (1..fdimc, 1..fdimr) array is
C               compressed to be just large enough to hold the
C               unassembled contribution block,
C               C (1..fdegc, 1..fdegr).

C-----------------------------------------------------------------------
C  Current frontal matrix
C-----------------------------------------------------------------------
C
C  ffxp points to current frontal matrix (contribution block and LU
C  factors).  For example, if fflefc = 4, fflefr = 6, luk = 3,
C  ffdimc = 8, ffdimr = 12, then "x" is a term in the contribution
C  block, "l" in L1, "u" in U1, "L" in L2, "U" in U2, and "." is unused.
C  XX (fxp) is "X". The first 3 pivot values (diagonal entries in U1)
C  are labelled 1, 2, and 3.  The frontal matrix is ffdimc-by-ffdimr.
C
C                   |----------- col 1 of L1 and L2, etc.
C                   V
C       X x x x x x L L L . . .
C       x x x x x x L L L . . .
C       x x x x x x L L L . . .
C       x x x x x x L L L . . .
C       U U U U U U 3 l l . . .         <- row 3 of U1 and U2
C       U U U U U U u 2 l . . .         <- row 2 of U1 and U2
C       U U U U U U u u 1 . . .         <- row 1 of U1 and U2
C       . . . . . . . . . . . .

C-----------------------------------------------------------------------
C  LU factors
C-----------------------------------------------------------------------
C
C   The LU factors are placed at the tail of XX.  If this routine
C   is factorizing a single block, then this description is for the
C   factors of the single block:
C
C       LUi (1..):      integer info. for LU factors
C       XX (xtail..xsize):      real values in LU factors
C
C   Each LU arrowhead (or factorized element) is stored as follows:
C   ---------------------------------------------------------------
C
C       total size: (7 + ludegc + ludegr + lunson) integers,
C                   (luk**2 + ludegc*luk + luk*ludegc) reals
C
C       If e is an LU arrowhead, then luip = LUp (e).
C
C       luxp   = LUi (luip) pointer to numerical LU arrowhead
C       luk    = LUi (luip+1) number of pivots in LU arrowhead
C       ludegr = LUi (luip+2) degree of last row of U (excl. diag)
C       ludegc = LUi (luip+3) degree of last col of L (excl. diag)
C       lunson = LUi (luip+4) number of children in assembly DAG
C       ffdimr = LUi (luip+5)
C       ffdimc = LUi (luip+6)
C                       max front size for this LU arrowhead is
C                       ffdimr-by-ffdimc, or zero if this LU arrowhead
C                       factorized within the frontal matrix of a prior
C                       LU arrowhead.
C       lucp   = (luip + 7)
C                       pointer to pattern of column of L
C       lurp   = lucp + ludegc
C                       pointer to patter of row of U
C       lusonp = lurp + ludegr
C                       pointer to list of sons in the assembly DAG
C       LUi (lucp ... lucp + ludegc - 1)
C                       row indices of column of L
C       LUi (lurp ... lurp + ludegr - 1)
C                       column indices of row of U
C       LUi (lusonp ... lusonp + lunson - 1)
C                       list of sons
C       XX (luxp...luxp + luk**2 + ludegc*luk + luk*ludegr - 1)
C                       pivot block (luk-by-luk) and the L block
C                       (ludegc-by-luk) in a single (luk+ludegc)-by-luk
C                       array, followed by the U block in a
C                       luk-by-ludegr array.
C
C   Pivot column/row pattern (also columns/rows in contribution block):
C       If the column/row index is negated, the column/row has been
C       assembled out of the frontal matrix into a subsequent frontal
C       matrix.  After factorization, the negative flags are removed.
C
C   List of sons:
C       1 <= son <= n:           son an LUson
C       n+1 <= son <= 2n:        son-n is an Uson
C       2n+n <= son <= 3n:       son-2n is a Lson

C-----------------------------------------------------------------------
C  Workspaces:
C-----------------------------------------------------------------------
C
C  WpC (1..ludegr):     holds the pivot column pattern
C                       (excluding the pivot row indices)
C
C  WpR (1..ludegr):     holds the pivot row pattern
C                       (excluding the pivot column indices)
C
C  WiR (row) >= 0 for each row in pivot column pattern.
C               offset into pattern is given by:
C               WiR (row) == offset - 1
C               Otherwise, WiR (1..n) is < 0
C
C  WiC (col) >= 0 for each col in pivot row pattern.
C               WiC (col) == (offset - 1) * ffdimc
C               Otherwise, WiC (1..n) is < 0
C
C  Wm (1..degc) or Wm (1..fdegc):       a gathered copy of WiR
C  Wj (1..degc) or Wj (1..fdegc):       offset in pattern of a son 

C-----------------------------------------------------------------------
C  Memory allocation in XX:
C-----------------------------------------------------------------------
C
C   XX (1..xhead):      values of original entries in arrowheads of
C                       matrix, values of contribution blocks, followed
C                       by the current frontal matrix.
C
C   mtail = nlu+2
C   mhead = nlu+1:      FRnext (mhead) points to the first contribution
C                       block in the head of XX.  The FRnext and FRprev
C                       arrays form a doubly-linked list.  Traversing
C                       the list from mhead to mtail gives the
C                       contribution blocks in ascending ordering of
C                       address (FRxp).  A block is free if FRdimc <= 0.
C                       The largest known free block in XX is pfree,
C                       located in
C                       XX (FRxp (pfree) ... FRxp (pfree) + xfree -1),
C                       unless pfree = 0, in which case no largest free
C                       block is known.

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER SWPCOL, SWPROW, FDIMC, K0, COLPOS, ROWPOS, PIVOT, FFPP,
     $          P, I, J, LUDEGR, LUDEGC, KPOS, SP, FFRP, FFCP, TYPE,
     $          FXP, LURP, LUCP, NEXT, FFLEFR, PREV, XHEAD, FDEGR,
     $          FFLEFC, K, XCDP, XDP, XSP, S, FDEGC, FLURP, FLUCP,
     $          COL, E, ROW, MHEAD, MTAIL, UXP, LUK, IO, FLUIP, LUSONP,
     $          FFSIZE, FFXP, FFDIMR, FFDIMC, XRDP, NPIV, NB, LUNSON,
     $          XNEED, LDIMR, LDIMC, LXP, PRL, XP, LUIP, PFREE, XFREE,
     $          XS, LUXP, FSP, FLP, FDP, DEGC, NZU, NZL, XRUSE
        LOGICAL PR3, ALLCOL, ALLROW
        DOUBLE PRECISION
     $          ONE, PIV, TEMP
        DOUBLE PRECISION
     $          TMP

C  Printing control:
C  -----------------
C  prl:     invalid entries printed if prl >= 3
C  io:      I/O unit for warning messages (printing invalid entries)
C  pr3:     true if invalid entries are to be printed when found
C
C  Current working array:
C  ----------------------
C  ffxp:    current working array is in XX (ffxp ... ffxp+ffsize-1)
C  ffsize:  size of current working array in XX
C  ffdimr:  row degree (number of columns) of current working array
C  ffdimc:  column degree (number of rows) of current working array
C  fflefr:  row degree (number of columns) of current contribution block
C  fflefc:  column degree (number of rows) of current contribution block
C  ffrp:     U2 block is in XX (ffrp ...)
C  ffcp:     L2 block is in XX (ffcp ...)
C  ffpp:     location in XX of the current pivot value
C
C  Current element:
C  ----------------
C  s:       current element being factorized
C  luip:    current element is in LUi (luip ...)
C  luk:     number of pivots in current element
C  ludegc:  degree of pivot column (excluding pivots themselves)
C  ludegr:  degree of pivot row (excluding pivots themselves)
C  ldimr:   row degree (number of columns) of current element
C  ldimc:   column degree (number of row) of current element
C  lucp:    pattern of col(s) of current element in LUi (lucp...)
C  lurp:    pattern of row(s) of current element in LUi (lurp...)
C  lusonp:  list of sons of current element is in LUi (lusonp...)
C  lunson:  number of sons of current element
C  sp:      pointer into list of sons of current element
C  luxp:    numerical values of LU arrowhead stored in XX (luxp ...)
C  lxp:     L2 block is stored in XX (lxp ...) when computed
C  uxp:     U2 block is stored in XX (uxp ...) when computed
C  nzu:     nonzeros above diagonal in U in current LU arrowhead
C  nzl:     nonzeros below diagonal in L in current LU arrowhead
C  swpcol:  the non-pivotal column to be swapped with pivot column
C  swprow:  the non-pivotal row to be swapped with pivot row
C  colpos:  position in WpR of the pivot column
C  rowpos:  position in WpC of the pivot row
C  kpos:    position in C to place pivot row/column
C  k:       current pivot is kth pivot of current element, k=1..luk
C  k0:      contribution block, C, has been updated with pivots 1..k0
C  npiv:    number of pivots factorized so far, excl. current element
C  pivot:   current pivot entry is A (pivot, pivot)
C  xcdp:    current pivot column is in XX (xcdp ...)
C  xrdp:    current pivot row is in XX (xrdp ...)
C
C  Son, or element other than current element:
C  -------------------------------------------
C  e:       an element other than s (a son of s, for example)
C  fluip:   LU arrowhead of e is in LUi (fluip ...)
C  fxp:     contribution block of son is in XX (fxp ...)
C  fdimc:   leading dimension of contribution block of a son
C  fdegr:   row degree of contribution block of son (number of columns)
C  fdegc:   column degree of contribution block of son (number of rows)
C  allcol:  true if all columns are present in son
C  allrow:  true if all rows are present in son
C  flucp:   pattern of col(s) of son in LUi (flucp...)
C  flurp:   pattern of row(s) of son in LUi (flurp...)
C  type:    an LUson (type = 1), Uson (type = 2) or Lson (type = 3)
C  degc:    compressed column offset vector of son is in Wj/Wm (1..degc)
C
C  Memory allocation:
C  ------------------
C  mhead:   nlu+1, head pointer for contribution block link list
C  mtail:   nlu+2, tail pointer for contribution block link list
C  prev:    FRprev (e) of the element e
C  next:    FRnext (e) of the element e
C  pfree:   FRxp (pfree) is the largest known free block in XX
C  xfree:   size of largest known free block in XX
C  xneed:   bare minimum memory currently needed in XX
C  xhead:   XX (1..xhead-1) is in use, XX (xhead ..) is free
C  xruse:   estimated memory needed in XX for next call to UMD2RF,
C           assuming a modest number of garbage collections
C  xs:      size of a block of memory in XX
C
C  Other:
C  ------
C  xdp:     destination pointer, into XX
C  xsp:     source pointer, into XX
C  xp:      a pointer into XX
C  fsp:     source pointer, into XX
C  fsp:     destination pointer, into XX
C  flp:     last row/column in current contribution is in XX (flp...)
C  col:     a column index
C  row:     a row index
C  nb:      block size for tradeoff between Level-2 and Level-3 BLAS  
C  p, i, j, x:  various uses

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

C       ----------------------------------------------------------------
C       get control parameters and initialize various scalars
C       ----------------------------------------------------------------

        ONE = 1
        IO = ICNTL (2)
        PRL = ICNTL (3)
        NB = MAX (1, ICNTL (7))
        NPIV = 0
        XHEAD = CP (1)
        XTAIL = XSIZE + 1
        XNEED = XUSE
        XRUSE = XUSE
        XRMAX = MAX (XRMAX, XRUSE)
        MHEAD = NLU+1
        MTAIL = NLU+2
        XFREE = -1
        PFREE = 0
        PR3 = PRL .GE. 3 .AND. IO .GE. 0

C       ----------------------------------------------------------------
C       initialize workspaces
C       ----------------------------------------------------------------

        DO 10 I = 1, N 
           WIR (I) = -1
           WIC (I) = -1
10      CONTINUE 

        DO 20 E = 1, NLU+2 
           FRDIMC (E) = 0
           FRXP (E) = 0
           FRNEXT (E) = 0
           FRPREV (E) = 0
20      CONTINUE 
        FRNEXT (MHEAD) = MTAIL
        FRPREV (MTAIL) = MHEAD
        FRXP (MHEAD) = XHEAD
        FRXP (MTAIL) = XHEAD

C       count the numerical assembly of the original matrix
        RINFO (2) = RINFO (2) + (NZ)

C       current working array is empty:
        FFLEFR = 0
        FFLEFC = 0
        FFSIZE = 0
        FFXP = XHEAD

C=======================================================================
C  Factorization [
C=======================================================================

        DO 600 S = 1, NLU 

C=======================================================================
C  Get the next element to factorize
C=======================================================================

           LUIP = LUP (S)
           LUK = LUI (LUIP+1)
           LUDEGC = LUI (LUIP+3)
           LUDEGR = LUI (LUIP+2)
           LUNSON = LUI (LUIP+4)
           LUCP = (LUIP + 7)
           LURP = LUCP + LUDEGC
           LUSONP = LURP + LUDEGR
           LDIMC = LUK + LUDEGC
           LDIMR = LUK + LUDEGR

C=======================================================================
C  Start new frontal matrix or merge with prior contribution block [
C=======================================================================

C          =============================================================
           IF (LUI (LUIP+6) .NE. 0) THEN 
C          start new contribution block
C          =============================================================

C             ----------------------------------------------------------
C             clear the prior offsets
C             ----------------------------------------------------------

              DO 30 I = 1, FFLEFR 
                 WIC (WPR (I)) = -1
30            CONTINUE 
              DO 40 I = 1, FFLEFC 
                 WIR (WPC (I)) = -1
40            CONTINUE 

C             ----------------------------------------------------------
C             save prior contribution block (s-1), if it exists
C             ----------------------------------------------------------

              XS = FFLEFR * FFLEFC
              IF (FFSIZE .NE. 0) THEN 
C                one more frontal matrix is finished
                 XNEED = XNEED - (FFSIZE - XS)
                 XRUSE = XRUSE - (FFSIZE - XS)
                 INFO (13) = INFO (13) + 1
C             else 
C                prior contribution block does not exist
              ENDIF 

              IF (FFLEFR .LE. 0 .OR. FFLEFC .LE. 0) THEN 

C                -------------------------------------------------------
C                if prior contribution block nonexistent or empty
C                -------------------------------------------------------

                 XUSE = XUSE - (XHEAD - FRXP (MTAIL))
                 XHEAD = FRXP (MTAIL)

              ELSE 

C                -------------------------------------------------------
C                prepare the prior contribution block for later assembly
C                -------------------------------------------------------

                 E = S - 1

C                count the numerical assembly
                 RINFO (2) = RINFO (2) + (XS)

                 IF (XS .LE. XFREE) THEN 

C                   ----------------------------------------------------
C                   compress and store in a freed block
C                   ----------------------------------------------------

C                   place the new block in the list
                    XFREE = XFREE - XS
                    IF (PFREE .EQ. MTAIL) THEN 
C                      place the new block at start of tail block
                       PREV = FRPREV (MTAIL)
                       NEXT = MTAIL
                       XDP = FRXP (MTAIL)
                       FRXP (MTAIL) = XDP + XS
                    ELSE 
C                      place the new block at end of block
                       PREV = PFREE
                       NEXT = FRNEXT (PFREE)
                       XDP = FRXP (NEXT) - XS
                       IF (XFREE .EQ. 0 .AND. PFREE .NE. MHEAD) THEN 
C                         delete the free block if its size is zero
                          PREV = FRPREV (PREV)
                          PFREE = 0
                          XFREE = -1
                       ENDIF 
                    ENDIF 
                    DO 60 J = 0, FFLEFR - 1 
CFPP$ NODEPCHK L
                       DO 50 I = 0, FFLEFC - 1 
                          XX (XDP+J*FFLEFC+I) = XX (FFXP+J*FFDIMC+I)
50                     CONTINUE 
60                  CONTINUE 
                    XUSE = XUSE - (XHEAD - FRXP (MTAIL))
                    XHEAD = FRXP (MTAIL)
                    FRXP (E) = XDP
                    FRDIMC (E) = FFLEFC

                 ELSE 

C                   ----------------------------------------------------
C                   deallocate part of unused portion of frontal matrix
C                   ----------------------------------------------------

C                   leave the contribution block C (1:fflefc, 1:fflefr)
C                   at head of XX, with column dimension of ffdimc and
C                   space of size (fflefr-1)*ffdimc for the first
C                   fflefr columns, and fflefc for the last column.
                    XS = FFSIZE - (FFLEFC + (FFLEFR-1)*FFDIMC)
                    XHEAD = XHEAD - XS
                    XUSE = XUSE - XS
                    PREV = FRPREV (MTAIL)
                    NEXT = MTAIL
                    FRXP (MTAIL) = XHEAD
                    FRXP (E) = FFXP
                    FRDIMC (E) = FFDIMC
                 ENDIF 

                 FRNEXT (PREV) = E
                 FRPREV (NEXT) = E
                 FRNEXT (E) = NEXT
                 FRPREV (E) = PREV

              ENDIF 

              IF (PFREE .EQ. MTAIL) THEN 
                 PFREE = 0
                 XFREE = -1
              ENDIF 

C             ----------------------------------------------------------
C             allocate a new ffdimr-by-ffdimc frontal matrix
C             ----------------------------------------------------------

              FFDIMC = LUI (LUIP+6)
              FFDIMR = LUI (LUIP+5)
              FFSIZE = FFDIMR * FFDIMC
              FFXP = 0

C             ----------------------------------------------------------
C             allocate and zero the space, garbage collection if needed
C             ----------------------------------------------------------

              IF (FFSIZE .GT. XTAIL-XHEAD) THEN 
                 INFO (15) = INFO (15) + 1
                 CALL UMD2RG (XX, XSIZE, XHEAD, XTAIL, XUSE,
     $              LUI, FRDIMC, FRXP, FRNEXT, FRPREV, NLU, LUP,
     $              ICNTL, FFXP, FFSIZE, PFREE, XFREE)
              ENDIF 

              FFXP = XHEAD
              XHEAD = XHEAD + FFSIZE
              XUSE = XUSE + FFSIZE
              XNEED = XNEED + FFSIZE
              XRUSE = XRUSE + FFSIZE
              XRMAX = MAX (XRMAX, XRUSE)
              INFO (20) = MAX (INFO (20), XUSE)
              INFO (21) = MAX (INFO (21), XNEED)
              IF (XHEAD .GT. XTAIL) THEN 
C                error return, if not enough real memory:
                 GO TO 9000
              ENDIF 

C             ----------------------------------------------------------
C             zero the frontal matrix
C             ----------------------------------------------------------

              DO 70 P = FFXP, FFXP + FFSIZE - 1 
                 XX (P) = 0
70            CONTINUE 

C             ----------------------------------------------------------
C             place pivot rows and columns in correct position
C             ----------------------------------------------------------

              DO 80 K = 1, LUK 
                 WIC (NPIV + K) = (LDIMR - K) * FFDIMC
                 WIR (NPIV + K) =  LDIMC - K
80            CONTINUE 

C             ----------------------------------------------------------
C             get the pivot row pattern of the new LU arrowhead
C             ----------------------------------------------------------

              DO 90 I = 0, LUDEGR - 1 
                 COL = LUI (LURP+I)
                 WIC (COL) = I * FFDIMC
                 WPR (I+1) = COL
90            CONTINUE 

C             ----------------------------------------------------------
C             get the pivot column pattern of the new LU arrowhead
C             ----------------------------------------------------------

              DO 100 I = 0, LUDEGC - 1 
                 ROW = LUI (LUCP+I)
                 WIR (ROW) = I
                 WPC (I+1) = ROW
100           CONTINUE 

C          =============================================================
           ELSE 
C          merge with prior contribution block
C          =============================================================

C             ----------------------------------------------------------
C             prior block is located at XX (ffxp ... ffxp + ffsize - 1).
C             It holds a working array C (1..ffdimc, 1..ffdimr), with a
C             prior contribution block in C (1..fflefc, 1..fflefr).
C             The last pivot column pattern is WpC (1..fflefc), and
C             the last pivot row pattern is WpR (1..fflefr).  The
C             offsets WiR and WiC are:
C             WiR (WpC (i)) = i-1, for i = 1..fflefc, and -1 otherwise.
C             WiC (WpR (i)) = (i-1)*ffdimc, for i = 1..fflefr, else -1.
C             The prior LU arrowhead is an implicit LUson of the current
C             element (and is implicitly assembled into the same
C             frontal matrix).
C             ----------------------------------------------------------

C             ----------------------------------------------------------
C             zero the newly extended frontal matrix
C             ----------------------------------------------------------

C             zero the new columns in the contribution and LU blocks
C             C (1..ldimc, fflefr+1..ldimr) = 0
              DO 120 J = FFLEFR, LDIMR - 1 
                 DO 110 I = 0, LDIMC - 1 
                    XX (FFXP + J*FFDIMC + I) = 0
110              CONTINUE 
120           CONTINUE 

C             C (fflefc+1..ldimc, 1..fflefr) = 0
C             zero the new rows in the contribution and U blocks
              DO 140 I = FFLEFC, LDIMC - 1 
CFPP$ NODEPCHK L
                 DO 130 J = 0, FFLEFR - 1 
                    XX (FFXP + J*FFDIMC + I) = 0
130              CONTINUE 
140           CONTINUE 

C             ----------------------------------------------------------
C             move pivot rows and columns into correct position
C             ----------------------------------------------------------

              DO 220 K = 1, LUK 

C                -------------------------------------------------------
C                kth pivot of frontal matrix, (npiv+k)th pivot of LU
C                -------------------------------------------------------

                 PIVOT = NPIV + K

C                -------------------------------------------------------
C                move the kth pivot column into position
C                -------------------------------------------------------

                 XSP = WIC (PIVOT)
                 KPOS = LDIMR - K + 1
                 XDP = (KPOS - 1) * FFDIMC
                 WIC (PIVOT) = XDP

                 IF (XSP .GE. 0) THEN 
C                   pivot column is already in current frontal matrix,
C                   shift into proper position
                    COLPOS = (XSP / FFDIMC) + 1
                    FSP = FFXP + XSP
                    FDP = FFXP + XDP

                    IF (FFLEFR .LT. KPOS) THEN 

                       IF (FFLEFR .EQ. COLPOS) THEN 

C                         ----------------------------------------------
C                         move C(:,colpos) => C (:,kpos)
C                         C (:,colpos) = 0
C                         ----------------------------------------------
CFPP$ NODEPCHK L
                          DO 150 I = 0, LDIMC - 1 
                             XX (FDP+I) = XX (FSP+I)
                             XX (FSP+I) = 0
150                       CONTINUE 

                       ELSE 

C                         ----------------------------------------------
C                         move C(:,colpos) => C (:,kpos)
C                         move C(:,fflefr) => C (:,colpos)
C                         C (:,fflefr) = 0
C                         ----------------------------------------------

                          FLP = FFXP + (FFLEFR - 1) * FFDIMC
CFPP$ NODEPCHK L
                          DO 160 I = 0, LDIMC - 1 
                             XX (FDP+I) = XX (FSP+I)
                             XX (FSP+I) = XX (FLP+I)
                             XX (FLP+I) = 0
160                       CONTINUE 

                          SWPCOL = WPR (FFLEFR)
                          WPR (COLPOS) = SWPCOL
                          WIC (SWPCOL) = XSP
                       ENDIF 

                    ELSE IF (COLPOS .NE. KPOS) THEN 

C                      -------------------------------------------------
C                      swap C (:,colpos) <=> C (:,kpos)
C                      -------------------------------------------------
CFPP$ NODEPCHK L
                       DO 180 I = 0, LDIMC - 1 
                          TEMP = XX (FDP+I)
                          XX (FDP+I) = XX (FSP+I)
                          XX (FSP+I) = TEMP
180                    CONTINUE 

                       SWPCOL = WPR (KPOS)
                       WPR (COLPOS) = SWPCOL
                       WIC (SWPCOL) = XSP
                    ENDIF 

                    FFLEFR = FFLEFR - 1
                 ENDIF 

C                -------------------------------------------------------
C                move the kth pivot row into position
C                -------------------------------------------------------

                 XSP = WIR (PIVOT)
                 KPOS = LDIMC - K + 1
                 XDP = (KPOS - 1)
                 WIR (PIVOT) = XDP

                 IF (XSP .GE. 0) THEN 
C                   pivot row is already in current frontal matrix,
C                   shift into proper position
                    ROWPOS = XSP + 1
                    FSP = FFXP + XSP
                    FDP = FFXP + XDP

                    IF (FFLEFC .LT. KPOS) THEN 

                       IF (FFLEFC .EQ. ROWPOS) THEN 

C                         ----------------------------------------------
C                         move C(rowpos,:) => C (kpos,:)
C                         C (rowpos,:) = 0
C                         ----------------------------------------------
CFPP$ NODEPCHK L
                          DO 190 J = 0, (LDIMR - 1) * FFDIMC, FFDIMC 
                             XX (FDP+J) = XX (FSP+J)
                             XX (FSP+J) = 0
190                       CONTINUE 

                       ELSE 

C                         ----------------------------------------------
C                         move C(rowpos,:) => C (kpos,:)
C                         move C(fflefc,:) => C (rowpos,:)
C                         C (fflefc,:) = 0
C                         ----------------------------------------------

                          FLP = FFXP + (FFLEFC - 1)
CFPP$ NODEPCHK L
                          DO 200 J = 0, (LDIMR - 1) * FFDIMC, FFDIMC 
                             XX (FDP+J) = XX (FSP+J)
                             XX (FSP+J) = XX (FLP+J)
                             XX (FLP+J) = 0
200                       CONTINUE 

                          SWPROW = WPC (FFLEFC)
                          WPC (ROWPOS) = SWPROW
                          WIR (SWPROW) = XSP
                       ENDIF 

                    ELSE IF (ROWPOS .NE. KPOS) THEN 

C                      -------------------------------------------------
C                      swap C (rowpos,:) <=> C (kpos,:)
C                      -------------------------------------------------
CFPP$ NODEPCHK L
                       DO 210 J = 0, (LDIMR - 1) * FFDIMC, FFDIMC 
                          TEMP = XX (FDP+J)
                          XX (FDP+J) = XX (FSP+J)
                          XX (FSP+J) = TEMP
210                    CONTINUE 

                       SWPROW = WPC (KPOS)
                       WPC (ROWPOS) = SWPROW
                       WIR (SWPROW) = XSP
                    ENDIF 

                    FFLEFC = FFLEFC - 1
                 ENDIF 

220           CONTINUE 

C             ----------------------------------------------------------
C             merge with pivot row pattern of new LU arrowhead
C             ----------------------------------------------------------

              I = FFLEFR
              DO 230 P = LURP, LURP + LUDEGR - 1 
                 COL = LUI (P)
                 IF (WIC (COL) .LT. 0) THEN 
                    WIC (COL) = I * FFDIMC
                    I = I + 1
                    WPR (I) = COL
                 ENDIF 
230           CONTINUE 

C             ----------------------------------------------------------
C             merge with pivot column pattern of new LU arrowhead
C             ----------------------------------------------------------

              I = FFLEFC
              DO 240 P = LUCP, LUCP + LUDEGC - 1 
                 ROW = LUI (P)
                 IF (WIR (ROW) .LT. 0) THEN 
                    WIR (ROW) = I
                    I = I + 1
                    WPC (I) = ROW
                 ENDIF 
240           CONTINUE 

           ENDIF 

C=======================================================================
C  Done initializing frontal matrix ]
C=======================================================================

C=======================================================================
C  Assemble original arrowheads into the frontal matrix, and deallocate
C=======================================================================

C          -------------------------------------------------------------
C          current workspace usage:
C          -------------------------------------------------------------

C          WpC (1..ludegr):     holds the pivot column pattern
C                               (excluding the pivot row indices)
C
C          WpR (1..ludegr):     holds the pivot row pattern
C                               (excluding the pivot column indices)
C
C          C (1..ffdimr, 1..ffdimc):  space for the frontal matrix,
C               in XX (ffxp ... ffxp + ffsize - 1)
C
C          C (i,j) is located at XX (ffxp+((i)-1)+((j)-1)*ffdimc)
C
C          C (1..ludegc, 1..ludegr):            contribution block
C          C (ludegc+1..ludegc+luk, 1..ludegr):             U2 block
C          C (1..ludegc, ludegr+1..ludegr+luk):             L2 block
C          C (ludegc+1..ludegc+luk, ludegr+1..ludegr+luk):  L1\U1 block
C
C          WiR (row) >= 0 for each row in pivot column pattern.
C               offset into pattern is given by:
C               WiR (row) == offset - 1
C               Also, WiR (npiv+1 ... npiv+luk) is
C               ludegc+luk-1 ... ludegc, the offsets of the pivot rows.
C
C               Otherwise, WiR (1..n) is < 0
C
C          WiC (col) >= 0 for each col in pivot row pattern.
C               WiC (col) == (offset - 1) * ffdimc
C               Also, WiC (npiv+1 ... npiv+luk) is
C               ludegr+luk-1 ... ludegr, the offsets of the pivot rows.
C
C               Otherwise, WiC (1..n) is < 0

           DO 260 K = 1, LUK 
              I = NPIV + K
              XCDP = FFXP + WIC (I)
              XRDP = FFXP + WIR (I)
              DO 250 P = CP (I+1), CP (I) - 1 
                 J = ARI (P)
                 IF (J .GT. 0) THEN 
C                   a diagonal entry, or lower triangular entry
C                   row = j, col = i
                    XP = XCDP + WIR (J)
                    IF (XP .LT. XCDP) THEN 
C                      invalid entry - not in prior LU pattern
                       NOUTSD = NOUTSD + 1
                       IF (PR3) THEN 
C                         get original row and column index and print it
                          ROW = RPERM (J)
                          COL = CPERM (I)
                          CALL UMD2P2 (2, 97, ROW, COL, XX (P), IO)
                       ENDIF 
                    ELSE 
                       XX (XP) = XX (XP) + XX (P)
                    ENDIF 
                 ELSE 
C                   an upper triangular entry
C                   row = i, col = -j
                    XP = XRDP + WIC (-J)
                    IF (XP .LT. XRDP) THEN 
C                      invalid entry - not in prior LU pattern
                       NOUTSD = NOUTSD + 1
                       IF (PR3) THEN 
C                         get original row and column index and print it
                          ROW = RPERM (I)
                          COL = CPERM (-J)
                          CALL UMD2P2 (2, 97, ROW, COL, XX (P), IO)
                       ENDIF 
                    ELSE 
                       XX (XP) = XX (XP) + XX (P)
                    ENDIF 
                 ENDIF 
250           CONTINUE 
260        CONTINUE 

C          deallocate the original arrowheads
           P = CP (NPIV + LUK + 1)
           XS = CP (NPIV + 1) - P
           FRXP (MHEAD) = P
           XNEED = XNEED - XS
           IF (XS .GT. XFREE) THEN 
              XFREE = XS
              PFREE = MHEAD
           ENDIF 

C=======================================================================
C  Assemble LUsons, Usons, and Lsons into the frontal matrix [
C=======================================================================

           DO 480 SP = LUSONP, LUSONP + LUNSON - 1 

C             ----------------------------------------------------------
C             get the son and determine its type (LUson, Uson, or Lson)
C             ----------------------------------------------------------

              E = LUI (SP)
              IF (E .LE. N) THEN 
C                LUson
                 TYPE = 1
              ELSE IF (E .LE. 2*N) THEN 
C                Uson
                 E = E - N
                 TYPE = 2
              ELSE 
C                Lson
                 E = E - 2*N
                 TYPE = 3
              ENDIF 

C             ----------------------------------------------------------
C             if fdimc=0 this is the implicit LUson (already assembled)
C             ----------------------------------------------------------

              FDIMC = FRDIMC (E)
              IF (FDIMC .NE. 0) THEN 

C                -------------------------------------------------------
C                get scalar info of the son (it needs assembling)
C                -------------------------------------------------------

                 FXP = FRXP (E)
                 FLUIP = LUP (E)
                 FDEGR = LUI (FLUIP+2)
                 FDEGC = LUI (FLUIP+3)
                 ALLCOL = FDEGR .GT. 0
                 ALLROW = FDEGC .GT. 0
                 FDEGR = ABS (FDEGR)
                 FDEGC = ABS (FDEGC)
                 FLUCP = (FLUIP + 7)
                 FLURP = FLUCP + FDEGC

C                use Wm (1..fdegc) for offsets:

C                -------------------------------------------------------
                 IF (TYPE .EQ. 1) THEN 
C                this is an LUson - assemble an entire frontal matrix
C                -------------------------------------------------------

C                   ----------------------------------------------------
                    IF (ALLROW) THEN 
C                   no rows assembled out of this LUson yet
C                   ----------------------------------------------------

C                      compute the compressed column offset vector
                       DO 270 I = 0, FDEGC-1 
                          ROW = LUI (FLUCP+I)
                          WM (I+1) = WIR (ROW)
270                    CONTINUE 

C                      -------------------------------------------------
                       IF (ALLCOL) THEN 
C                      no rows or cols assembled out of LUson yet
C                      -------------------------------------------------

                          DO 290 J = 0, FDEGR-1 
                             COL = LUI (FLURP+J)
                             XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                             DO 280 I = 0, FDEGC-1 
                                XX (XDP + WM (I+1)) =
     $                          XX (XDP + WM (I+1)) +
     $                          XX (FXP + J*FDIMC + I)
280                          CONTINUE 
290                       CONTINUE 

C                      -------------------------------------------------
                       ELSE 
C                      some columns already assembled out of LUson
C                      -------------------------------------------------

                          DO 310 J = 0, FDEGR-1 
                             COL = LUI (FLURP+J)
                             IF (COL .GT. 0) THEN 
                                XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                                DO 300 I = 0, FDEGC-1 
                                   XX (XDP + WM (I+1)) =
     $                             XX (XDP + WM (I+1)) +
     $                             XX (FXP + J*FDIMC + I)
300                             CONTINUE 
                             ENDIF 
310                       CONTINUE 

                       ENDIF 

C                   ----------------------------------------------------
                    ELSE 
C                   some rows already assembled out of LUson
C                   ----------------------------------------------------

C                      compute the compressed column offset vector
                       DEGC = 0
                       DO 320 I = 0, FDEGC-1 
                          ROW = LUI (FLUCP+I)
                          IF (ROW .GT. 0) THEN 
                             DEGC = DEGC + 1
                             WJ (DEGC) = I
                             WM (DEGC) = WIR (ROW)
                          ENDIF 
320                    CONTINUE 

C                      -------------------------------------------------
                       IF (ALLCOL) THEN 
C                      some rows already assembled out of LUson
C                      -------------------------------------------------

                          DO 340 J = 0, FDEGR-1 
                             COL = LUI (FLURP+J)
                             XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                             DO 330 I = 1, DEGC 
                                XX (XDP + WM (I)) =
     $                          XX (XDP + WM (I)) +
     $                          XX (FXP + J*FDIMC + WJ (I))
330                          CONTINUE 
340                       CONTINUE 

C                      -------------------------------------------------
                       ELSE 
C                      rows and columns already assembled out of LUson
C                      -------------------------------------------------

                          DO 360 J = 0, FDEGR-1 
                             COL = LUI (FLURP+J)
                             IF (COL .GT. 0) THEN 
                                XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                                DO 350 I = 1, DEGC 
                                   XX (XDP + WM (I)) =
     $                             XX (XDP + WM (I)) +
     $                             XX (FXP + J*FDIMC + WJ (I))
350                             CONTINUE 
                             ENDIF 
360                       CONTINUE 

                       ENDIF 
                    ENDIF 

C                   ----------------------------------------------------
C                   deallocate the LUson frontal matrix
C                   ----------------------------------------------------

                    FRDIMC (E) = 0
                    PREV = FRPREV (E)
                    NEXT = FRNEXT (E)
                    XNEED = XNEED - FDEGR*FDEGC
                    XRUSE = XRUSE - FDEGR*FDEGC

                    IF (FRDIMC (PREV) .LE. 0) THEN 
C                      previous block is free - delete this block
                       FRNEXT (PREV) = NEXT
                       FRPREV (NEXT) = PREV
                       E = PREV
                       PREV = FRPREV (E)
                    ENDIF 

                    IF (FRDIMC (NEXT) .LE. 0) THEN 
C                      next block is free - delete this block
                       FRXP (NEXT) = FRXP (E)
                       IF (E .LE. NLU) THEN 
                          FRNEXT (PREV) = NEXT
                          FRPREV (NEXT) = PREV
                       ENDIF 
                       E = NEXT
                       NEXT = FRNEXT (E)
                       IF (FRNEXT (MHEAD) .EQ. MTAIL) THEN 
C                         no blocks left except mhead and mtail
                          FRXP (MTAIL) = FRXP (MHEAD)
                       ENDIF 
                    ENDIF 

C                   get the size of the freed block
                    IF (NEXT .EQ. 0) THEN 
C                      this is the mtail block
                       XS = FFXP - FRXP (E)
                    ELSE 
                       XS = FRXP (NEXT) - FRXP (E)
                    ENDIF 
                    IF (XS .GT. XFREE) THEN 
C                      keep track of the largest free block
                       XFREE = XS
                       PFREE = E
                    ENDIF 

C                -------------------------------------------------------
                 ELSE IF (TYPE .EQ. 2) THEN 
C                Uson:  assemble all possible columns
C                -------------------------------------------------------

C                   ----------------------------------------------------
                    IF (ALLROW) THEN 
C                   no rows assembled out of this Uson yet
C                   ----------------------------------------------------

C                      compute the compressed column offset vector
                       DO 370 I = 0, FDEGC-1 
                          ROW = LUI (FLUCP+I)
                          WM (I+1) = WIR (ROW)
370                    CONTINUE 

                       DO 390 J = 0, FDEGR-1 
                          COL = LUI (FLURP+J)
                          IF (COL .GT. 0) THEN 
                             IF (WIC (COL) .GE. 0) THEN 
                                XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                                DO 380 I = 0, FDEGC-1 
                                   XX (XDP + WM (I+1)) =
     $                             XX (XDP + WM (I+1)) +
     $                             XX (FXP + J*FDIMC + I)
380                             CONTINUE 
C                               flag this column as assembled
                                LUI (FLURP+J) = -COL
                             ENDIF 
                          ENDIF 
390                    CONTINUE 

C                   ----------------------------------------------------
                    ELSE 
C                   some rows already assembled out of this Uson
C                   ----------------------------------------------------

C                      compute the compressed column offset vector
                       DEGC = 0
                       DO 400 I = 0, FDEGC-1 
                          ROW = LUI (FLUCP+I)
                          IF (ROW .GT. 0) THEN 
                             DEGC = DEGC + 1
                             WJ (DEGC) = I
                             WM (DEGC) = WIR (ROW)
                          ENDIF 
400                    CONTINUE 

                       DO 420 J = 0, FDEGR-1 
                          COL = LUI (FLURP+J)
                          IF (COL .GT. 0) THEN 
                             IF (WIC (COL) .GE. 0) THEN 
                                XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                                DO 410 I = 1, DEGC 
                                   XX (XDP + WM (I)) =
     $                             XX (XDP + WM (I)) +
     $                             XX (FXP + J*FDIMC + WJ (I))
410                             CONTINUE 
C                               flag this column as assembled
                                LUI (FLURP+J) = -COL
                             ENDIF 
                          ENDIF 
420                    CONTINUE 

                    ENDIF 

C                   flag this element as missing some columns
                    LUI (FLUIP+2) = -FDEGR

C                -------------------------------------------------------
                 ELSE 
C                Lson:  assemble all possible rows
C                -------------------------------------------------------

C                   compute the compressed column offset vector
                    DEGC = 0
                    DO 430 I = 0, FDEGC-1 
                       ROW = LUI (FLUCP+I)
                       IF (ROW .GT. 0) THEN 
                          IF (WIR (ROW) .GE. 0) THEN 
C                            this row will be assembled in loop below
                             DEGC = DEGC + 1
                             WJ (DEGC) = I
                             WM (DEGC) = WIR (ROW)
C                            flag this row as assembled
                             LUI (FLUCP+I) = -ROW
                          ENDIF 
                       ENDIF 
430                 CONTINUE 

C                   ----------------------------------------------------
                    IF (ALLCOL) THEN 
C                   no columns assembled out of this Lson yet
C                   ----------------------------------------------------

                       DO 450 J = 0, FDEGR-1 
                          COL = LUI (FLURP+J)
                          XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                          DO 440 I = 1, DEGC 
                             XX (XDP + WM (I)) =
     $                       XX (XDP + WM (I)) +
     $                       XX (FXP + J*FDIMC + WJ (I))
440                       CONTINUE 
450                    CONTINUE 

C                   ----------------------------------------------------
                    ELSE 
C                   some columns already assembled out of this Lson
C                   ----------------------------------------------------

                       DO 470 J = 0, FDEGR-1 
                          COL = LUI (FLURP+J)
                          IF (COL .GT. 0) THEN 
                             XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                             DO 460 I = 1, DEGC 
                                XX (XDP + WM (I)) =
     $                          XX (XDP + WM (I)) +
     $                          XX (FXP + J*FDIMC + WJ (I))
460                          CONTINUE 
                          ENDIF 
470                    CONTINUE 

                    ENDIF 

C                   flag this element as missing some rows
                    LUI (FLUIP+3) = -FDEGC

                 ENDIF 

              ENDIF 

480        CONTINUE 

C=======================================================================
C  Done assemblying sons into the frontal matrix ]
C=======================================================================

C=======================================================================
C  Factorize the frontal matrix [
C=======================================================================

           K0 = 0
           FFLEFR = LDIMR
           FFLEFC = LDIMC
           FFCP = FFXP + FFLEFR * FFDIMC
           FFRP = FFXP + FFLEFC
           FFPP = FFXP + FFLEFC + FFLEFR * FFDIMC

           DO 500 K = 1, LUK 

C             ----------------------------------------------------------
C             compute kth column of U1, and update pivot column
C             ----------------------------------------------------------

              IF (K-K0-2 .GT. 0) THEN 
C                u1 = L1 \ u1.  Note that L1 transpose is stored, and
C                that u1 is stored with rows in reverse order.
                 CALL DTRSV ('U', 'N', 'U', K-K0-1,
     $                         XX (FFPP         ), FFDIMC,
     $                         XX (FFPP - FFDIMC), 1)
                 TMP = K-K0-2
                 TMP = TMP * (K-K0-1)
                 RINFO (5) = RINFO (5) + 2.0*(TMP) / 2.0
              ENDIF 
              IF (K-K0-1 .GT. 0) THEN 
C                l1 = l1 - L2*u1
                 CALL DGEMV ('N', FFLEFC, K-K0-1,
     $                   -ONE, XX (FFCP         ), FFDIMC,
     $                         XX (FFPP - FFDIMC), 1,
     $                    ONE, XX (FFCP - FFDIMC), 1)
                 TMP = FFLEFC
                 TMP = TMP * (K-K0-1)
                 RINFO (5) = RINFO (5) + 2.0*(TMP)
              ENDIF 

              FFCP = FFCP - FFDIMC
              FFRP = FFRP - 1
              FFPP = FFPP - FFDIMC - 1
              FFLEFR = FFLEFR - 1
              FFLEFC = FFLEFC - 1

C             ----------------------------------------------------------
C             divide pivot column by pivot
C             ----------------------------------------------------------

C             k-th pivot in frontal matrix located in XX (ffpp)
              PIV = XX (FFPP)
              IF (ABS (PIV) .EQ. 0) THEN 
C                error return, if pivot order from UMD2FA not acceptable
                 GO TO 9010
              ENDIF 
              PIV = 1 / PIV
              DO 490 P = FFCP, FFCP + FFLEFC - 1 
                 XX (P) = XX (P) * PIV
490           CONTINUE 
C             count this as a call to the Level-1 BLAS:
              RINFO (4) = RINFO (4) + (FFLEFC)
              INFO (17) = INFO (17) + 1

C             ----------------------------------------------------------
C             compute U1 (k0+1..k, k..ldimc) and
C             update contribution block: rank-nb, or if last pivot
C             ----------------------------------------------------------

              IF (K-K0 .GE. NB .OR. K .EQ. LUK) THEN 
                 CALL DTRSM ('L', 'U', 'N', 'U', K-K0, FFLEFR, ONE,
     $                      XX (FFPP), FFDIMC,
     $                      XX (FFRP), FFDIMC)
                 TMP = FFLEFR
                 TMP = TMP * (K-K0-1)
                 TMP = TMP * (K-K0)
                 RINFO (6) = RINFO (6) + 2.0*(TMP) / 2.0
                 CALL DGEMM ('N', 'N', FFLEFC, FFLEFR, K-K0,
     $                -ONE, XX (FFCP ), FFDIMC,
     $                      XX (FFRP ), FFDIMC,
     $                 ONE, XX (FFXP), FFDIMC)
                 TMP = FFLEFC
                 TMP = TMP * FFLEFR
                 TMP = TMP * (K-K0)
                 RINFO (6) = RINFO (6) + 2.0*(TMP)
                 K0 = K
              ENDIF 

500        CONTINUE 

C=======================================================================
C  Done factorizing the frontal matrix ]
C=======================================================================

C=======================================================================
C  Save the new LU arrowhead [
C=======================================================================

C          allocate permanent space for the LU arrowhead
           XS = LUK*LUDEGC + LUK*LUDEGR + LUK*LUK

           IF (XS .GT. XTAIL-XHEAD) THEN 
              INFO (15) = INFO (15) + 1
              CALL UMD2RG (XX, XSIZE, XHEAD, XTAIL, XUSE,
     $              LUI, FRDIMC, FRXP, FRNEXT, FRPREV, NLU, LUP,
     $              ICNTL, FFXP, FFSIZE, PFREE, XFREE)
           ENDIF 

           XTAIL = XTAIL - XS
           LUXP = XTAIL
           XUSE = XUSE + XS
           XNEED = XNEED + XS
           XRUSE = XRUSE + XS
           XRMAX = MAX (XRMAX, XRUSE)
           INFO (20) = MAX (INFO (20), XUSE)
           INFO (21) = MAX (INFO (21), XNEED)
           IF (XHEAD .GT. XTAIL) THEN 
C             error return, if not enough real memory:
              GO TO 9000
           ENDIF 

C          save the scalar data of the LU arrowhead
           LUI (LUIP) = LUXP

C          save column pattern (it may have been rearranged)
           DO 510 I = 0, LUDEGC-1 
              LUI (LUCP+I) = WPC (I+1)
510        CONTINUE 

C          save row pattern (it may have been rearranged)
           DO 520 I = 0, LUDEGR-1 
              LUI (LURP+I) = WPR (I+1)
520        CONTINUE 

C          move the L1,U1 matrix, compressing the dimension from
C          ffdimc to ldimc.  The LU arrowhead grows on top of stack.
           XP = FFXP + (LDIMR-1)*FFDIMC + LDIMC-1
           DO 540 J = 0, LUK-1 
CFPP$ NODEPCHK L
              DO 530 I = 0, LUK-1 
                 XX (LUXP + J*LDIMC + I) = XX (XP - J*FFDIMC - I)
530           CONTINUE 
540        CONTINUE 

C          move L2 matrix, compressing dimension from ffdimc to ldimc
           IF (LUDEGC .NE. 0) THEN 
              LXP = LUXP + LUK
              XP = FFXP + (LDIMR-1)*FFDIMC
              DO 560 J = 0, LUK-1 
CFPP$ NODEPCHK L
                 DO 550 I = 0, LUDEGC-1 
                    XX (LXP + J*LDIMC + I) = XX (XP - J*FFDIMC + I)
550              CONTINUE 
560           CONTINUE 
           ENDIF 

C          move the U2 block.
           IF (LUDEGR .NE. 0) THEN 
              UXP = LUXP + LUK * LDIMC
              XP = FFXP + LDIMC-1
              DO 580 J = 0, LUDEGR-1 
CFPP$ NODEPCHK L
                 DO 570 I = 0, LUK-1 
                    XX (UXP + J*LUK + I) = XX (XP + J*FFDIMC - I)
570              CONTINUE 
580           CONTINUE 
           ENDIF 

C          one more LU arrowhead has been refactorized
           NZU = (LUK*(LUK-1)/2) + LUK*LUDEGC
           NZL = (LUK*(LUK-1)/2) + LUK*LUDEGR
           INFO (10) = INFO (10) + NZL
           INFO (11) = INFO (11) + NZU

C          -------------------------------------------------------------
C          clear the pivot row and column offsets
C          -------------------------------------------------------------

           DO 590 PIVOT = NPIV + 1, NPIV + LUK 
              WIR (PIVOT) = -1
              WIC (PIVOT) = -1
590        CONTINUE 
           NPIV = NPIV + LUK

C=======================================================================
C  Done saving the new LU arrowhead ]
C=======================================================================

600     CONTINUE 

C=======================================================================
C  Factorization complete ]
C=======================================================================

C=======================================================================
C  Wrap-up:  store LU factors in their final form
C=======================================================================

C       ----------------------------------------------------------------
C       Flag remaining arrowheads as invalid entries, if prior matrix
C       was singular.  Print them if requested.
C       ----------------------------------------------------------------

        IF (NPIV .LT. N) THEN 
           IF (PR3) THEN 
              DO 620 I = NPIV+1, N 
                 DO 610 P = CP (I+1), CP (I) - 1 
                    J = ARI (P)
                    IF (J .GT. 0) THEN 
C                      a diagonal entry, or lower triangular entry
C                      get original row and column index
                       ROW = RPERM (J)
                       COL = CPERM (I)
                    ELSE 
C                      an upper triangular entry
C                      get original row and column index
                       ROW = RPERM (I)
                       COL = CPERM (-J)
                    ENDIF 
                    CALL UMD2P2 (2, 95, ROW, COL, XX(P), IO)
610              CONTINUE 
620           CONTINUE 
           ENDIF 
           NOUTSD = NOUTSD + (CP (NPIV+1) - CP (N+1))
        ENDIF 

C       ----------------------------------------------------------------
C       deallocate all remaining input arrowheads and frontal matrices
C       ----------------------------------------------------------------

        IF (FFSIZE .NE. 0) THEN 
           INFO (13) = INFO (13) + 1
        ENDIF 
        XUSE = XUSE - (XHEAD - CP (N+1))
        XNEED = XUSE
        XHEAD = CP (N+1)

        IF (NLU .EQ. 0) THEN 
C          LU factors are completely empty (A = 0).
C          Add one real, to simplify rest of code.
C          Otherwise, some arrays in UMD2RF or UMD2SO would have
C          zero size, which can cause an address fault.
           XTAIL = XSIZE
           XUSE = XUSE + 1
           XRUSE = XUSE
           XNEED = XUSE
           INFO (20) = MAX (INFO (20), XUSE)
           INFO (21) = MAX (INFO (21), XNEED)
        ENDIF 

        IF (XHEAD .LE. XTAIL) THEN 

C          -------------------------------------------------------------
C          sufficient memory to complete the factorization
C          -------------------------------------------------------------

           IF (NLU .EQ. 0) THEN 
C             zero the dummy entry, although it won't be accessed:
              XX (XTAIL) = 0
           ENDIF 

C          -------------------------------------------------------------
C          update pointers in LU factors
C          -------------------------------------------------------------

           DO 630 S = 1, NLU 
              LUIP = LUP (S)
              LUXP = LUI (LUIP)
              LUI (LUIP) = LUXP - XTAIL + 1
630        CONTINUE 

C          -------------------------------------------------------------
C          get memory usage estimate for next call to UMD2RF
C          -------------------------------------------------------------

           XRUSE = XUSE
           XRMAX = MAX (XRMAX, XRUSE)
           RETURN

        ENDIF 

C=======================================================================
C  Error conditions
C=======================================================================

C       error return label:
9000    CONTINUE
C       out of real memory
        CALL UMD2ER (2, ICNTL, INFO, -4, INFO (21))
        RETURN

C       error return label:
9010    CONTINUE
C       original pivot order computed by UMD2FA is no longer acceptable
        CALL UMD2ER (2, ICNTL, INFO, -6, 0)
        RETURN
        END 
