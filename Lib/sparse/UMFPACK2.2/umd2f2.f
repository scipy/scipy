        SUBROUTINE UMD2F2 (CP, NZ, N, PN, CPERM, RPERM, ITAIL, XTAIL,
     $          XX, XSIZE, II, ISIZE, ICNTL, CNTL, INFO, RINFO, PGIVEN,
     $          IUSE, XUSE, WIR, WIC, WPR, WPC, WM, HEAD,
     $          WJ, RP, WC, WR, DN, DSIZ, KEEP,
     $          RMAX, CMAX, TOTNLU, XRMAX, XRUSE)
        INTEGER XSIZE, ISIZE, ICNTL (20), INFO (40), PN,
     $          ITAIL, XTAIL, NZ, N, II (ISIZE), CP (N+1), DN, DSIZ,
     $          RPERM (PN), CPERM (PN), WIR (N), WIC (N), WPR (N),
     $          WPC (N), WM (N), HEAD (N), RP (N+DN), WC (N+DN),
     $          WR (N+DN), IUSE, XUSE, WJ (N), KEEP (20),
     $          RMAX, CMAX, TOTNLU, XRMAX, XRUSE
        LOGICAL PGIVEN
        DOUBLE PRECISION
     $          XX (XSIZE)
        DOUBLE PRECISION
     $          CNTL (10), RINFO (20)
        
C=== UMD2F2 ============================================================
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
C  UMD2F2 factorizes the n-by-n input matrix at the head of II/XX
C  (in expanded column-form) and places its LU factors at the tail of
C  II/XX.  The input matrix is overwritten.   No BTF information is
C  used in this routine, except that the BTF permutation arrays are
C  modified to include the final permutations.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       Cp (1..n+1):    column pointers of expanded column-form,
C                       undefined on output
C       n:              order of input matrix
C       nz:             entries in input matrix
C       isize:          size of II
C       xsize:          size of XX
C       iuse:           memory usage in Index
C       xuse:           memory usage in Value
C       Icntl:          integer control parameters, see UMD21I
C       Cntl:           real control parameters, see UMD21I
C       Keep (6)        integer control parameter, see UMD21I
C       dn:             number of dense columns
C       dsiz:           entries required for col to be treated as dense
C       rmax:           maximum ludegr seen so far (see below)
C       cmax:           maximum ludegc seen so far (see below)
C       totnlu:         total number of LU arrowheads constructed so far
C       xrmax:          maximum real memory usage for UMD2RF
C       xruse:          current real memory usage for UMD2RF
C
C       pgiven:         true if Cperm and Rperm are defined on input
C       if pgiven then:
C          Cperm (1..pn):       col permutation to BTF, n = pn
C          Rperm (1..pn):       row permutation to BTF
C       else
C          Cperm (1..pn):       unaccessed pn = 1
C          Rperm (1..pn):       unaccessed
C
C       II (1..nz+cscal*n):             expanded column-form, see below
C       II (nz+cscal*n+1..isize):       undefined on input
C       XX (1..nz):                     expanded column-form, see below
C       XX (nz+1..xsize):               undefined on input

C=======================================================================
C  WORKSPACE:
C=======================================================================
C
C       WiR (1..n)
C       WiC (1..n)
C       WpR (1..n)
C       WpC (1..n)
C       Wm (1..n)
C       Head (n)
C       Rp (1..n+dn)
C       Wr (1..n+dn)
C       Wc (1..n+dn)

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       II (1..itail-1):        undefined on output
C       II (itail..isize):      LU factors of this matrix, see below
C       XX (1..xtail-1):        undefined on output
C       XX (xtail..xsize):      LU factors of this matrix, see below
C
C       Info:           integer informational output, see UMD2FA
C       Rinfo:          real informational output, see UMD2FA
C       if pgiven:
C          Cperm (1..n): the final col permutations, including BTF
C          Rperm (1..n): the final row permutations, including BTF
C
C       WiC (1..n):     row permutations, not including BTF
C       WiR (1..n):     column permutations, not including BTF
C
C       iuse:           memory usage in Index
C       xuse:           memory usage in Value
C       rmax:           maximum ludegr seen so far (see below)
C       cmax:           maximum ludegc seen so far (see below)
C       totnlu:         total number of LU arrowheads constructed so far
C       xrmax:          maximum real memory usage for UMD2RF
C       xruse:          current real memory usage for UMD2RF

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   UMD2F1
C       subroutines called:     UMD2ER, UMD2FG, DGEMV, DGEMM 
C       functions called:       IDAMAX, ABS, MAX, MIN
        INTEGER IDAMAX
        INTRINSIC ABS, MAX, MIN

C=======================================================================
C  DESCRIPTION OF DATA STRUCTURES:
C=======================================================================

C-----------------------------------------------------------------------
C  Column/element/arrowhead pointers:
C-----------------------------------------------------------------------
C
C  The Cp (1..n) array contains information about non-pivotal columns
C
C       p = Cp (j)
C       if (p = 0) then j is a pivotal column
C       else i is a non-pivotal column
C
C  The Rp (1..n) array contains information about non-pivotal rows,
C  unassembled frontal matrices (elements), and the LU arrowheads
C
C       p = Rp (i)
C       if (i > n) then
C          i is an artificial frontal matrix (a dense column)
C          if (p = 0) then i is assembled, else unassembled
C       else if (p = 0) then i is pivotal but not element/arrowhead
C       else if (Wc (i) >= 0 and Wc (i) <= n) then
C          i is a non-pivotal row
C       else if (Wc (i) = -(n+dn+2)) then
C          i is a pivotal row, an assembled element, and an LU arrowhead
C       else i an unassembled element

C-----------------------------------------------------------------------
C  Matrix being factorized:
C-----------------------------------------------------------------------
C
C    Each column is stored in II and XX:
C    -----------------------------------
C
C       if j is a non-pivotal column, pc = Cp (j):
C
C       csiz = II (pc) size of the integer data structure for col j,
C                        including the cscal scalars
C       cdeg = II (pc+1) degree of column j
C       cxp  = II (pc+2) pointer into XX for numerical values
C       next = II (pc+3) pointer to next block of memory in XX
C       prev = II (pc+4) pointer to previous block of memory in XX
C       celn = II (pc+5) number of elements in column j element list
C       clen = II (pc+6) number of original entries in column j
C       cnxt = II (pc+7) next column with same degree as col j
C       cprv = II (pc+8) previous column with same degree as col j
C       cep = (pc+9) pointer to start of the element list
C       II (cep ... cep + 2*celn - 1)
C                       element list (e,f) for the column
C       II (cep + 2*celn ... pc + csiz - clen - 1)
C                       empty
C       II (pc + csiz - clen ... pc + csiz - 1)
C                       row indices of original nonzeros in the column
C       XX (xp ... xp + clen - 1)
C                       numerical values of original nonzeros in the col
C
C       if cdeg = II (pc+1) = -(n+2), then this is a singular column
C       if cdeg = -1, then this column is deallocated
C
C    Each row is stored in II only:
C    ------------------------------
C
C       if i is a non-pivotal row, pr = Rp (i)
C
C       rsiz = II (pr) size of the integer data structure for row i,
C                        including the rscal scalars
C       rdeg = II (pr+1) degree of row i
C       reln = Wr (i) number of elements in row i element list
C       rlen = Wc (i) number of original entries in row i
C       rep  = (pr+2) pointer to start of the element list
C       II (rep ... rep + 2*reln - 1)
C                       element list (e,f) for the row
C       II (rep + 2*reln ... pr + rsiz - rlen - 1)
C                       empty
C       II (pr + rsiz - rlen ... pr + rsiz - 1)
C                       column indices of original nonzeros in the row
C
C       if rdeg = -1, then this row is deallocated

C-----------------------------------------------------------------------
C  Frontal matrices
C-----------------------------------------------------------------------
C
C   Each unassembled frontal matrix (element) is stored as follows:
C       total size: fscal integers, (fdimr*fdimc) reals
C
C       if e is an unassembled element, ep = Rp (e), and e is also
C       the first pivot row in the frontal matrix.
C
C       fluip  = II (ep)        pointer to LU arrowhead in II
C       fdimc  = II (ep+1)      column dimension of contribution block
C       fxp    = II (ep+2)      pointer to contribution block in XX
C       next   = II (ep+3)      pointer to next block in XX
C       prev   = II (ep+4)      pointer to previous block in XX
C       fleftr = II (ep+5)      number of unassembled rows
C       fleftc = II (ep+6)      number of unassembled columns
C       fextr = Wr (e) - w0     external row degree of the frontal mtx
C       fextc = Wc (e) - w0     external col degree of the frontal mtx
C       XX (fxp ... )
C               a 2-dimensional array, C (1..fdimc, 1..fdimr).
C               note that fdimr is not kept (it is not needed,
C               except for the current frontal).  If this is not the
C               current frontal matrix, then luip points to the
C               corresponding LU arrowhead, and the contribution block
C               is stored in C (1..ludegc, 1..ludegr) in the
C               C (1..fdimc, ...) array.
C
C               If memory is limited, garbage collection will occur.
C               In this case, the C (1..fdimc, 1..fdimr) array is
C               compressed to be just large enough to hold the
C               unassembled contribution block,
C               C (1..ludegc, 1..ludegr).

C-----------------------------------------------------------------------
C  Artificial frontal matrices
C-----------------------------------------------------------------------
C
C   An artificial frontal matrix is an original column that is treated
C   as a c-by-1 frontal matrix, where c is the number of original
C   nonzeros in the column.  Dense columns (c > dsiz) are treated this
C   way.  An artificial frontal matrix is just the same as a frontal
C   matrix created by the elimination of one or more pivots, except
C   that there is no corresponding LU arrowhead.  The row and column
C   patterns are stored in:
C
C       ep = Rp (e), where e = n+1 .. n+dn, where there are dn
C                    artificial frontal matrices.
C               
C       lucp = (ep+9)   pointer to row pattern (just one column index)
C       lurp = (ep+8) pointer to column pattern (fdimc row indices)

C-----------------------------------------------------------------------
C  Current frontal matrix
C-----------------------------------------------------------------------
C
C  ffxp points to current frontal matrix (contribution block and LU
C  factors).  For example, if fflefc = 4, fflefr = 6, k = 3, and
C  gro = 2.0, then "x" is a term in the contribution block, "l" in L1,
C  "u" in U1, "L" in L2, "U" in U2, and "." is unused.  XX (fxp) is "X".
C  The first 3 pivot values (diagonal entries in U1) are 1,2, and 3.
C  For this frontal matrix, ffdimr = 12 (the number of columns), and
C  ffdimc = 8 (the number of rows).  The frontal matrix is
C  ffdimc-by-ffdimr
C
C                             |----------- col 1 of L1 and L2, etc.
C                             V
C       X x x x x x . . . L L L
C       x x x x x x . . . L L L
C       x x x x x x . . . L L L
C       x x x x x x . . . L L L
C       . . . . . . . . . . . .
C       U U U U U U . . . 3 l l         <- row 3 of U1 and U2
C       U U U U U U . . . u 2 l         <- row 2 of U1 and U2
C       U U U U U U . . . u u 1         <- row 1 of U1 and U2

C-----------------------------------------------------------------------
C  LU factors
C-----------------------------------------------------------------------
C
C   The LU factors are placed at the tail of II and XX.  If this routine
C   is factorizing a single block, then this decription is for the
C   factors of the single block:
C
C       II (itail):             xtail = start of LU factors in XX
C       II (itail+1):           nlu = number of LU arrowheads
C       II (itail+2):           npiv = number of pivots
C       II (itail+3):           maximum number of rows in any
C                               contribution block (max ludegc)
C       II (itail+4):           maximum number of columns in any
C                               contribution block (max ludegr)
C       II (itail+5..itail+nlu+4): LUp (1..nlu) array, pointers to each
C                               LU arrowhead, in order of their
C                               factorization
C       II (itail+nlu+5...isize):integer info. for LU factors
C       XX (xtail..xsize):      real values in LU factors
C
C   Each LU arrowhead is stored as follows:
C   ---------------------------------------
C
C       total size: (7 + ludegc + ludegr + nsons) integers,
C                   (luk**2 + ludegc*luk + luk*ludegc) reals
C
C       If e is an LU arrowhead, then luip = Rp (e), and luip >= itail.
C       When UMD2F2 returns, then luip is given by luip =
C       II (itail+s+1), where s = 1..nlu is the position of the LU
C       arrowhead in the LU factors (s=1,2,.. refers to the first,
C       second,.. LU arrowhead)
C
C       luxp   = II (luip) pointer to numerical LU arrowhead
C       luk    = II (luip+1) number of pivots in LU arrowhead
C       ludegr = II (luip+2) degree of last row of U (excl. diag)
C       ludegc = II (luip+3) degree of last col of L (excl. diag)
C       nsons  = II (luip+4) number of children in assembly DAG
C       ludimr = II (luip+5)
C       ludimc = II (luip+5) max front size is ludimr-by-ludimc,
C                       or zero if this LU arrowhead factorized within
C                       the frontal matrix of a prior LU arrowhead.
C       lucp   = (luip + 7)
C                       pointer to pattern of column of L
C       lurp   = lucp + ludegc
C                       pointer to patter of row of U
C       lusonp = lurp + ludegr
C                       pointer to list of sons in the assembly DAG
C       II (lucp ... lucp + ludegc - 1)
C                       row indices of column of L
C       II (lurp ... lurp + ludegr - 1)
C                       column indices of row of U
C       II (lusonp ... lusonp + nsons - 1)
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
C       matrix.  After factorization, the negative flags are removed,
C       and the row/col indices are replaced with their corresponding
C       index in the permuted LU factors.
C
C   List of sons:
C       1 <= son <= n:           son an LUson
C       n+1 <= son <= 2n:        son-n is an Uson
C       2n+n <= son <= 3n:       son-2n is a Lson
C       during factorzation, a son is referred to by its first
C       pivot column.  After factorization, they are numbered according
C       to their order in the LU factors.

C-----------------------------------------------------------------------
C  Workspaces:
C-----------------------------------------------------------------------
C
C   WiR (e):  link list of sons of the current element
C       WiR (e) = -1 means that e is not in the list.
C       WiR (e) = next+n+2 means that "next" is the element after e.
C       The end of the list is marked with WiR (e) = -(n+2).
C       sonlst points to the first element in the list, or 0 if
C       the sonlst is empty.
C
C   WiR (row), WiC (col):  used for pivot row/col offsets:
C
C       If WiR (row) >= 0 then the row is in the current
C       column pattern.  Similarly for WiC (col).
C
C       If WiR (row) is set to "empty" (<= -1), then
C       the row is not in the current pivot column pattern.
C
C       Similarly, if WiC (col) is set to -2, then the column is
C       not in the current pivot row pattern.
C
C       If WiC (col) = -1 then col is pivotal
C
C       After factorization, WiR/C holds the pivot permutations.
C
C   WpR/C (1..n):  the first part is used for the current frontal
C           matrix pattern.  During factorization, the last part holds
C           a stack of the row and column permutations (WpR/C (n-k+1)
C           is the k-th pivot row/column).
C
C   Head (1..n):        degree lists for columns.  Head (d) is the
C                       first column in list d with degree d.
C                       The cnxt and cprv pointers are stored in the
C                       column data structure itself.
C                       mindeg is the least non-empty list
C
C   Wm (1..n):          various uses
C   Wj (1..degc) or Wj (1..fdegc):      offset in pattern of a son 

C-----------------------------------------------------------------------
C  Memory allocation in II and XX:
C-----------------------------------------------------------------------
C
C   II (1..ihead):      rows and columns of active submatrix, and
C                       integer information for frontal matrices.
C   XX (1..xhead):      values of original entries in columns of
C                       matrix, values of contribution blocks, followed
C                       by the current frontal matrix.
C
C   mhead:              a pointer to the first block in the head of
C                       XX.  Each block (a column or frontal matrix)
C                       contains a next and prev pointer for this list.
C                       If the list is traversed starting at mhead,
C                       then the pointers to the reals (cxp or fxp)
C                       will appear in strictly increasing order.
C                       Note that the next, prev, and real pointers
C                       are in II.  next and prev point to the next
C                       and previous block in II, and the real pointer
C                       points to the real part in XX.
C
C   mtail:              the end of the memory list.

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER SWPCOL, SWPROW, FDIMC, K0, COLPOS, ROWPOS, ROW2, RDEG2,
     $          P, I, J, FFROW, PIVROW, PIVCOL, LUDEGR, LUDEGC, E1,
     $          FXP, LURP, LUCP, IP, NEXT, FFLEFR, PC, MNEXT, MPREV,
     $          FFLEFC, FEDEGR, FEDEGC, K, XUDP, XDP, XSP, XLP, S, COL2,
     $          BESTCO, COL, E, ROW, COST, SRCHED, PR, F1, RSCAN, REP,
     $          KLEFT1, FFSIZE, FFXP, W0, FFDIMR, FFDIMC, KLEFT, XLDP,
     $          EP, SCAN1, SCAN2, SCAN3, SCAN4, NZL, NZU, DEGC, CEP
        INTEGER MINDEG, NSRCH, NPIV, ESON, LUIP1, DNZ, IWORST, WXP,
     $          NB, LUPP, NLU, NSONS, INEED, XNEED, LDIMC, LXP, RLEN2,
     $          RSIZ, LSONS, SONLST, XHEAD, IHEAD, DELN, DLEN,
     $          SLIST, XP, LUIP, RDEG, CDEG1, PFREE, XFREE, CDEG2,
     $          F, CDEG, MTAIL, MHEAD, RSIZ2, CSIZ2, IP2, MAXDR, MAXDC,
     $          XS, IS, LUXP, FSP, FLP, FDP, JJ, USONS, NDN, P2,
     $          CSIZ, CELN, CLEN, RELN, RLEN, UXP, PC2, PR2
        INTEGER CNXT, CPRV, CXP, FLUIP, LUSONP, FLEFTR, FLEFTC, MAXINT,
     $          FMAXR, FMAXC, SLOTS, LIMIT, RSCAL, CSCAL, FSCAL, EXTRA,
     $          FMAX, W0BIG, MINMEM, DUMMY1, DUMMY2, DUMMY3, DUMMY4
        LOGICAL SYMSRC, PFOUND, MOVELU, OKCOL, OKROW, BETTER
        DOUBLE PRECISION
     $          TEMP, ONE
        DOUBLE PRECISION
     $          TOLER, MAXVAL, RELPT, GRO, APIV, TMP
        PARAMETER (RSCAL = 2, CSCAL = 9, FSCAL = 7,
     $          MINMEM = 24)

C  Current element and working array, C:
C  -------------------------------------
C  ffxp:    current working array is in XX (ffxp ... ffxp+ffsize-1)
C  ffsize:  size of current working array in XX
C  ffdimr:  row degree (number of columns) of current working array
C  ffdimc:  column degree (number of rows) of current working array
C  fflefr:  row degree (number of columns) of current contribution block
C  fflefc:  column degree (number of rows) of current contribution block
C  fmaxr:   max row degree (maximum front size is fmaxr-by-fmaxc)
C  fmaxc:   max col degree (maximum front size is fmaxr-by-fmaxc)
C  fedegr:  extended row degree
C  fedegc:  extended column degree
C  ffrow:   current element being factorized (a pivot row index)
C  pivrow:  current pivot row index
C  pivcol:  current pivot column index
C  e1:      first pivot row in the frontal matrix
C  gro:     frontal matrix amalgamation growth factor
C  usons:   pointer to a link list of Usons, in Wc, assembled this SCAN3
C  lsons:   pointer to a link list of Lsons, in Wr, assembled this SCAN4
C  sonlst:  pointer to a link list of sons, in WiR, of current element
C  swpcol:  the non-pivotal column to be swapped with pivot column
C  swprow:  the non-pivotal row to be swapped with pivot row
C  colpos:  position in WpR of the pivot column
C  rowpos:  position in WpC of the pivot row
C  k:       current pivot is kth pivot of current element
C  k0:      contribution block, C, has been updated with pivots 1..k0
C
C  LU arrowhead (a factorized element):
C  ------------------------------------
C  movelu:  true if a new LU arrowhead is to be created
C  luip:    current element is in II (luip ...)
C  luip1:   first element from current frontal matrix in II (luip1...) 
C  ludegc:  degree of pivot column (excluding pivots themselves)
C  ludegr:  degree of pivot row (excluding pivots themselves)
C  lucp:    pattern of col(s) of current element in II (lucp...)
C  lurp:    pattern of row(s) of current element in II (lurp...)
C  lusonp:  list of sons of current element is in II (lusonp...)
C  nsons:   number of sons of current element
C  ldimc:   column dimension (number of rows) of [L1\U1 L2] block
C  luxp:    numerical values of LU arrowhead stored in XX (luxp ...)
C  lxp:     L2 block is stored in XX (lxp ...) when computed
C  uxp:     U2 block is stored in XX (uxp ...) when computed
C  nzu:     nonzeros above diagonal in U in current LU arrowhead
C  nzl:     nonzeros below diagonal in L in current LU arrowhead
C
C  Son, or element other than current element:
C  -------------------------------------------
C  e:       an element
C  eson:    an element
C  s:       a renumbered element (1..nlu) for UMD2SO and UMD2RF
C  ep:      frontal matrix integer data struct. in II (ep...ep+fscal-1)
C  fscal:   = 7, size of frontal matrix data structure
C  fluip:   LU arrowhead of e is in II (fluip ...)
C  fxp:     contribution block of son is in XX (fxp ...)
C  fdimc:   leading dimension of contribution block of e
C  lucp:    pattern of col(s) of e in II (lucp...)
C  lurp:    pattern of row(s) of e in II (lurp...)
C  ludegr:  row degree of contribution block of e
C  ludegr:  column degree of contribution block of e
C  maxdr:   maximum ludegr for any LU arrowhead, for UMD2RF
C  maxdc:   maximum ludegc for any LU arrowhead, for UMD2RF
C  degc:    compressed column offset vector of son is in Wj/Wm (1..degc)
C  fleftr:  remaining row degree (number of columns) of a contrib. block
C  fleftc:  remaining column degree (number of rows) of a contrib. block
C  xudp:    pointer to a column of a prior contribution block
C  xldp:    pointer to a row of a prior contribution block
C
C  Memory allocation:
C  ------------------
C  mhead:   head pointer for link list of blocks in XX
C  mtail:   tail pointer for link list of blocks in XX
C  mprev:   previous block, II (p+4), of the block located at p
C  mnext:   next block, II (p+3), of the block located at p
C  pfree:   II (pfree+2) is the largest known free block in XX
C  xfree:   size of largest known free block in XX
C  xhead:   XX (1..xhead-1) is in use, XX (xhead ..xtail-1) is free
C  xtail:   XX (xtail..xsize) is in use, XX (xhead ..xtail-1) is free
C  xneed:   bare minimum memory currently needed in XX
C  ihead:   II (1..ihead-1) is in use, II (ihead ..itail-1) is free
C  itail:   II (itail..isize) is in use, II (ihead ..itail-1) is free
C  ineed:   bare minimum memory currently needed in II
C  iworst:  worst possible current integer memory required
C  xs:      size of a block of memory in XX
C  is:      size of a block of memory in II
C  wxp:     pointer to a temporary workspace in XX (wxp ... )
C  slots:   number of slots added to element lists during garbage coll.
C  minmem:  smallest isize allowed
C
C  Wr and Wc flag arrays:
C  ----------------------
C  w0:      marker value for Wr (1..n) and Wc (1..n) arrays 
C  w0big:   largest permissible value of w0 (w0+n must not overflow)
C  fmax:    largest row/col degree of an element seen so far
C
C  A column:
C  ---------
C  pc:      pointer to a column, in II (pc...)
C  pc2:     pointer to a column, in II (pc2...)
C  csiz:    size of integer data structure of a column
C  csiz2:   size of integer data structure of a column
C  cscal:   = 9, number of scalars in data structure of a column
C  cdeg:    degree of a column
C  cdeg1:   degree of a column
C  cdeg2:   degree of a column
C  celn:    number of elements in the element list of a column
C  clen:    number of original entries that remain in a column
C  cnxt:    next column with same degree as this column
C  cprv:    previous column with same degree as this column
C  cep:     pointer to the element list of a column
C  cxp:     pointer to the numerical values in a column
C  limit:   maximum size for row/col data structure (excl. scalars)
C
C  Dense columns:
C  --------------
C  dnz:     number of original entries that reside in "dense" columns
C  dn:      number of "dense" columns
C  ndn:     n + dn
C  extra:   number of extra slots to add to reconstructed "dense" cols
C
C  A row:
C  ------
C  pr:      pointer to a row, in II (pr...)
C  pr2:     pointer to a row, in II (pr2...)
C  rsiz:    size of integer data structure of a row
C  rsiz2:   size of integer data structure of a row
C  rscal:   = 2, number of scalars in data structure of a row
C  rdeg:    degree of a row
C  rdeg2:   degree of a row
C  reln:    number of elements in the element list of a row
C  rlen:    number of original entries that remain in a row
C  rlen2:   number of original entries that remain in a row
C  rep:     pointer to the element list of a row
C
C  Pivot search:
C  -------------
C  cost:    approximate Markowitz-cost of the current candidate pivot
C  bestco:  best approximate Markowitz-cost seen so far
C  srched:  number of non-singular candidates searched so far
C  mindeg:  minimum degree of columns in active submatrix
C  nsrch:   maximum number of columns to search
C  slist:   pointer to a link list of searched columns, in II
C  symsrc:  true if attempting to preserve symmetry
C  pfound:  true if pivot found during local search
C  okcol:   true if candidate pivot column is acceptable, so far
C  okrow:   true if candidate pivot row is acceptable, so far
C  toler:   pivot tolerance; abs(pivot) must be >= toler
C  maxval:  maximum absolute value in a candidate pivot column
C  relpt:   relative pivot tolerance (Cntl (1))
C  npiv:    number of pivots factorized so far, incl. current element
C  kleft:   number of rows/columns remaining in active submatrix
C  kleft1:  kleft - 1
C  better:  if true, then candidate is better than the prior candidate
C
C  Assembly:
C  ---------
C  f1:      degree prior to assembly next item
C  f:       offset into an element
C  rscan:   skip row assembly if more than rscan original entries
C  scan1:   start SCAN1 at WpC (scan1 ... fflefc) 
C  scan2:   start SCAN2 at WpR (scan2 ... fflefr) 
C  scan3:   start SCAN3 at WpR (scan3 ... fflefr) 
C  scan4:   start SCAN4 at WpC (scan4 ... fflefc) 
C  deln:    number of (e,f) tuples to delete from an element list
C  dlen:    number of original entries to delete from a row/col
C
C  Allocated arrays:
C  -----------------
C  lupp:    LUp (1..nlu) array located in II (lupp...lupp+nlu-1)
C  nlu:     number of LU arrowheads
C
C  Other:
C  ------
C  xdp:     destination pointer, into XX
C  xsp:     source pointer, into XX
C  xlp:     pointer into XX of location of last row/col in C 
C  xp:      pointer into XX
C  ip:      pointer into II
C  ip2:     pointer into II
C  p2:      pointer into II
C  fsp:     source pointer, into XX
C  fsp:     destination pointer, into XX
C  flp:     last row/column in current contribution is in XX (flp...)
C  col,col2: a column index
C  row,row2: a row index
C  nb:      block size for tradeoff between Level-2 and Level-3 BLAS  
C  p, i, j, k, x:  various uses
C  jj:      loop index
C  maxint:  largest representable positive integer
C  next:    next pointer, for a link list
C  dummy1:  dummy loop index for main factorization loop
C  dummy2:  dummy loop index for global pivot search loop
C  dummy3:  dummy loop index for outer frontal matrix factorization loop
C  dummy4:  dummy loop index for inner frontal matrix factorization loop

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

C       ----------------------------------------------------------------
C       get control parameters and initialize various scalars
C       ----------------------------------------------------------------

        ONE = 1
        NSRCH = MAX (1, ICNTL (5))
        SYMSRC = ICNTL (6) .NE. 0
        NB = MAX (1, ICNTL (7))
        RELPT = CNTL (1)
        IF (RELPT .LT. 0) RELPT = 0
        IF (RELPT .GT. 1) RELPT = 1
        GRO = CNTL (2)
        IF (GRO .LT. 1) GRO = 1
        MAXINT = KEEP (6)
        NDN = N + DN
        W0BIG = MAXINT - N
        W0 = NDN + 2
C       currently: w0 = n+dn+2 < 2n+2 < w0big = maxint - n
C       2n+2 < maxint - n must hold, so n < (maxint - 2) / 3 is the
C       largest that n can be.  This condition is checked in UMD2FA.
        KLEFT = N
        NPIV = 0
        NLU = 0
        MINDEG = 1
        FMAX = 1
        IHEAD = NZ + CSCAL*N + 1
        XHEAD = NZ + 1
        ITAIL = ISIZE + 1
        XTAIL = XSIZE + 1
C       Cp (1) must equal 1, the first block
        XFREE = -1
        PFREE = 0
C       make sure integer space is at least of size minmem (simplifies
C       link list management and memory management)
        INFO (19) = MAX (INFO (19), IUSE+MINMEM)
        IF (IHEAD.GT.ITAIL.OR.ISIZE.LT.MINMEM.OR.XHEAD.GT.XTAIL) THEN 
C          error return, if not enough integer and/or real memory:
           GO TO 9000
        ENDIF 
        BESTCO = 0
        LIMIT = N + 2*NDN
        LSONS = NDN + 1
        USONS = NDN + 1

C       ----------------------------------------------------------------
C       initialize workspaces
C       ----------------------------------------------------------------

        DO 10 I = 1, N 
           WIR (I) = -1
           WIC (I) = -2
           HEAD (I) = 0
           WC (I) = 0
           WR (I) = 0
10      CONTINUE 

C       ----------------------------------------------------------------
C       initialize the link list for keeping track of real memory usage
C       ----------------------------------------------------------------

        MHEAD = 0
        MTAIL = 0
        DO 20 COL = 1, N 
           PC = CP (COL)
           CLEN = II (PC+6)
           IF (CLEN .GT. 0) THEN 
C             place the column in the link list of blocks in XX
              IF (MHEAD .EQ. 0) THEN 
                 MHEAD = PC
              ENDIF 
              II (PC+4) = MTAIL
              II (PC+3) = 0
              IF (MTAIL .NE. 0) THEN 
                 II (MTAIL+3) = PC
              ENDIF 
              MTAIL = PC
           ELSE 
              II (PC+2) = 0
              II (PC+4) = 0
              II (PC+3) = 0
           ENDIF 
20      CONTINUE 

C       ----------------------------------------------------------------
C       convert dense columns to a-priori contribution blocks and
C       get the count of nonzeros in each row
C       ----------------------------------------------------------------

        E = N
        DNZ = 0
        DO 50 COL = 1, N 
           PC = CP (COL)
           CLEN = II (PC+6)
           CEP = (PC+9)
           IF (CLEN .GT. DSIZ) THEN 
C             this is a dense column - add to element list length
              DNZ = DNZ + CLEN
              DO 30 IP = CEP, CEP + CLEN - 1 
                 ROW = II (IP)
                 WR (ROW) = WR (ROW) + 1
30            CONTINUE 
C             convert dense column (in place) into a frontal matrix
              E = E + 1
              EP = PC
              RP (E) = EP
              FDIMC = CLEN
              FLEFTC = CLEN
              FLEFTR = 1
              II (EP+1) = FDIMC
              II (EP+5) = FLEFTR
              II (EP+6) = FLEFTC
              WR (E) = W0-1
              WC (E) = W0-1
              LURP = (EP+8)
              II (LURP) = COL
              FMAX = MAX (FMAX, FLEFTC)
           ELSE 
C             this is a sparse column - add to orig entry length
              DO 40 IP = CEP, CEP + CLEN - 1 
                 ROW = II (IP)
                 WC (ROW) = WC (ROW) + 1
40            CONTINUE 
           ENDIF 
50      CONTINUE 

C       ----------------------------------------------------------------
C       get memory for row-oriented form, and dense column element lists
C       ----------------------------------------------------------------

        PR = IHEAD
        CSIZ = CSCAL + 2
        IS = (NZ + RSCAL*N + DNZ) + (DN * CSIZ)
        IHEAD = IHEAD + IS
        IUSE = IUSE + IS
        INEED = IUSE
        XNEED = XUSE
        INFO (18) = MAX (INFO (18), IUSE)
        INFO (19) = MAX (INFO (19), INEED)
        IF (IHEAD .GT. ITAIL) THEN 
C          error return, if not enough integer memory:
           GO TO 9000
        ENDIF 

C       ----------------------------------------------------------------
C       if memory is available, add up to dsiz+6 extra slots in the
C       reconstructed dense columns to allow for element list growth
C       ----------------------------------------------------------------

        IF (DN .GT. 0) THEN 
           EXTRA = MIN ((ITAIL - IHEAD) / DN, DSIZ + 6)
           CSIZ = CSIZ + EXTRA
           IS = DN * EXTRA
           IHEAD = IHEAD + IS
           IUSE = IUSE + IS
           INFO (18) = MAX (INFO (18), IUSE)
        ENDIF 

C       ----------------------------------------------------------------
C       construct row pointers
C       ----------------------------------------------------------------

        DO 60 ROW = 1, N 
           RP (ROW) = PR
           REP  = (PR+2)
           RELN = WR (ROW)
           RLEN = WC (ROW)
           RSIZ = 2*RELN + RLEN + RSCAL
           II (PR) = RSIZ
           RDEG = RELN + RLEN
           II (PR+1) = RDEG
           WM (ROW) = REP
           PR = PR + RSIZ
60      CONTINUE 

C       ----------------------------------------------------------------
C       construct row element lists for dense columns
C       ----------------------------------------------------------------

        PC = PR
        DO 80 E = N+1, N+DN 
           EP = RP (E)
           LUCP = (EP+9)
           FDIMC = II (EP+1)
CFPP$ NODEPCHK L
           DO 70 F = 0, FDIMC - 1 
              ROW = II (LUCP+F)
              II (WM (ROW)    ) = E
              II (WM (ROW) + 1) = F
              WM (ROW) = WM (ROW) + 2
70         CONTINUE 
C          re-construct dense columns as just an element list,
C          containing a single element tuple (e,f), where f = 0
           LURP = (EP+8)
           COL = II (LURP)
           CP (COL) = PC
           II (PC) = CSIZ
           CDEG = FDIMC
           II (PC+1) = CDEG
           II (PC+2) = 0
           II (PC+4) = 0
           II (PC+3) = 0
           II (PC+5) = 1
           II (PC+6) = 0
           II (PC+7) = 0
           II (PC+8) = 0
C          store the (e,0) tuple:
           CEP = (PC+9)
           II (CEP  ) = E
           II (CEP+1) = 0
           PC = PC + CSIZ
80      CONTINUE 

C       ----------------------------------------------------------------
C       construct the nonzero pattern of the row-oriented form
C       ----------------------------------------------------------------

        DO 100 COL = 1, N 
           PC = CP (COL)
           CEP = (PC+9)
           CLEN = II (PC+6)
CFPP$ NODEPCHK L
           DO 90 P = CEP, CEP + CLEN - 1 
              ROW = II (P)
              II (WM (ROW)) = COL
              WM (ROW) = WM (ROW) + 1
90         CONTINUE 
100     CONTINUE 

C       count the numerical assembly of the original matrix
        RINFO (2) = RINFO (2) + (NZ)

C       ----------------------------------------------------------------
C       initialize the degree lists
C       ----------------------------------------------------------------

C       do so in reverse order to try to improve pivot tie-breaking
        DO 110 COL = N, 1, -1 
           PC = CP (COL)
           CDEG = II (PC+1)
           IF (CDEG .LE. 0) THEN 
C             empty column - remove from pivot search
              CDEG = -(N+2)
              II (PC+1) = CDEG
           ELSE 
              CNXT = HEAD (CDEG)
              II (PC+7) = CNXT
              II (PC+8) = 0
              IF (CNXT .NE. 0) THEN 
                 II (CP (CNXT)+8) = COL
              ENDIF 
              HEAD (CDEG) = COL
           ENDIF 
110     CONTINUE 

C=======================================================================
C=======================================================================
C  MAIN FACTORIZATION LOOP [
C=======================================================================
C=======================================================================

        DO 1540 DUMMY1 = 1, N 
C       (this loop is not indented due to its length)

C       ----------------------------------------------------------------
C       factorization is done if N pivots have been found
C       ----------------------------------------------------------------

        IF (NPIV .GE. N) THEN 
           GO TO 2000
        ENDIF 

C=======================================================================
C  Global pivot search, and initialization of a new frontal matrix [
C=======================================================================

        IF (MTAIL .NE. 0 .AND. II (MTAIL+6) .EQ. 0) THEN 
C          tail block is free, delete it
           XP = II (MTAIL+2)
           XUSE = XUSE - (XHEAD - XP)
           XHEAD = XP
           IF (MTAIL .EQ. PFREE) THEN 
              PFREE = 0
              XFREE = -1
           ENDIF 
           MTAIL = II (MTAIL+4)
           IF (MTAIL .NE. 0) THEN 
              II (MTAIL+3) = 0
           ELSE 
C             singular matrix.  No columns or contribution blocks left.
              MHEAD = 0
           ENDIF 
        ENDIF 

C=======================================================================
C  Global pivot search:  find pivot row and column
C=======================================================================

        NSONS = 0
        SONLST = 0
        SRCHED = 0
        PIVCOL = 0
        SLIST = 0

        DO 255 DUMMY2 = 1, N 

C          -------------------------------------------------------------
C          get col from column upper-bound degree list
C          -------------------------------------------------------------

           COL = 0
           DO 140 CDEG = MINDEG, N 
              COL = HEAD (CDEG)
              IF (COL .NE. 0) THEN 
C                exit out of loop if column found:
                 GO TO 150
              ENDIF 
140        CONTINUE 
           IF (COL .EQ. 0) THEN 
C             exit out of loop if column not found (singular matrix):
              GO TO 260
           ENDIF 
C          loop exit label:
150        CONTINUE
           PC = CP (COL)
           CNXT = II (PC+7)
           IF (CNXT .NE. 0) THEN 
              II (CP (CNXT)+8) = 0
           ENDIF 
           HEAD (CDEG) = CNXT
           MINDEG = CDEG

C          -------------------------------------------------------------
C          construct candidate column in Wm and XX (wxp..wxp+cdeg-1)
C          -------------------------------------------------------------

           XS = CDEG
C          use Wm (1..cdeg) for pattern [
C          use XX (wxp..wxp+xs-1) as workspace for values [

           IF (XS .GT. XTAIL-XHEAD) THEN 

              INFO (15) = INFO (15) + 1
              CALL UMD2FG (XX, XSIZE, XHEAD, XTAIL, XUSE,
     $                     II, ISIZE, IHEAD, ITAIL, IUSE,
     $                     CP, RP, DN, N, ICNTL, WIR, WIC, WR, WC,
     $                     0, 0, 0, 0, .FALSE.,
     $                     PFREE, XFREE, MHEAD, MTAIL, SLOTS)
C             at this point, iuse = ineed and xuse = xneed
              PC = CP (COL)
           ENDIF 

           WXP = XHEAD
           XHEAD = XHEAD + XS
           XUSE = XUSE + XS
           XNEED = XNEED + XS
           INFO (20) = MAX (INFO (20), XUSE)
           INFO (21) = MAX (INFO (21), XNEED)
           IF (XHEAD .GT. XTAIL) THEN 
C             error return, if not enough real memory:
              GO TO 9000
           ENDIF 

C          -------------------------------------------------------------
C          assemble the elements in the element list
C          -------------------------------------------------------------

           CDEG = 0
           CEP = (PC+9)
           CELN = II (PC+5)
           DO 190 IP = CEP, CEP + 2*CELN - 2, 2 
              E = II (IP)
              F = II (IP+1)
              EP = RP (E)
              FDIMC = II (EP+1)
              FXP = II (EP+2)
              IF (E .LE. N) THEN 
                 FLUIP = II (EP)
                 LUDEGC = II (FLUIP+3)
                 LUCP = (FLUIP + 7)
              ELSE 
                 LUDEGC = FDIMC
                 LUCP = (EP+9)
              ENDIF 
              XP = FXP + F * FDIMC
C             split into 3 loops so that they all vectorize on a CRAY
              CDEG1 = CDEG
              DO 160 P = LUCP, LUCP + LUDEGC - 1 
                 ROW = II (P)
                 IF (ROW .GT. 0) THEN 
                    IF (WIR (ROW) .LE. 0) THEN 
                       CDEG = CDEG + 1
                       WM (CDEG) = ROW
                    ENDIF 
                 ENDIF 
160           CONTINUE 
              DO 170 I = CDEG1+1, CDEG 
                 ROW = WM (I)
                 WIR (ROW) = I
                 XX (WXP+I-1) = 0
170           CONTINUE 
CFPP$ NODEPCHK L
              DO 180 J = 0, LUDEGC - 1 
                 ROW = II (LUCP+J)
                 IF (ROW .GT. 0) THEN 
                    XX (WXP + WIR (ROW) - 1) =
     $              XX (WXP + WIR (ROW) - 1) + XX (XP+J)
                 ENDIF 
180           CONTINUE 
190        CONTINUE 

C          -------------------------------------------------------------
C          assemble the original entries in the column
C          -------------------------------------------------------------

           CDEG1 = CDEG
           CLEN = II (PC+6)
           CSIZ = II (PC)
           IP = PC + CSIZ - CLEN
           CXP = II (PC+2)
CFPP$ NODEPCHK L
           DO 200 I = 0, CLEN - 1 
              ROW = II (IP+I)
              WM (CDEG+1+I) = ROW
              XX (WXP+CDEG+I) = XX (CXP+I)
200        CONTINUE 
           CDEG = CDEG + CLEN

C          -------------------------------------------------------------
C          update the degree of this column (exact, not upper bound)
C          -------------------------------------------------------------

           II (PC+1) = CDEG

C          Wm (1..cdeg) holds the pattern of col being searched.
C          XX (wxp..wxp+cdeg-1) holds the numerical values of col being
C          searched.  WiR (Wm (1..cdeg1)) is 1..cdeg1.

C          -------------------------------------------------------------
C          find the maximum absolute value in the column
C          -------------------------------------------------------------

           MAXVAL = ABS (XX (WXP - 1 + IDAMAX (CDEG, XX (WXP), 1)))
           RINFO (3) = RINFO (3) + (CDEG)
           TOLER = RELPT * MAXVAL
           RDEG = N+1

C          -------------------------------------------------------------
C          look for the best possible pivot row in this column
C          -------------------------------------------------------------

           IF (CDEG .NE. 0 .AND. MAXVAL .GT. 0) THEN 
              IF (SYMSRC) THEN 
C                prefer symmetric pivots, if numerically acceptable
                 ROW = COL
                 ROWPOS = WIR (ROW)
                 IF (ROWPOS .LE. 0) THEN 
C                   diagonal may be in original entries
                    DO 210 I = CDEG1 + 1, CDEG1 + CLEN 
                       IF (WM (I) .EQ. ROW) THEN 
                          ROWPOS = I
C                         exit out of loop if symmetric pivot found:
                          GO TO 220
                       ENDIF 
210                 CONTINUE 
C                   loop exit label:
220                 CONTINUE
                 ENDIF 
                 IF (ROWPOS .GT. 0) THEN 
C                   diagonal entry exists in the column pattern
                    APIV = ABS (XX (WXP-1+ROWPOS))
                    IF (APIV .GE. TOLER .AND. APIV .GT. 0) THEN 
C                      diagonal entry is numerically acceptable
                       PR = RP (ROW)
                       RDEG = II (PR+1)
                    ENDIF 
                 ENDIF 
              ENDIF 
              IF (RDEG .EQ. N+1) THEN 
C                Continue searching - no diagonal found or sought for.
C                Minimize row degree subject to abs(value) constraints.
                 ROW = N+1
                 DO 230 I = 1, CDEG 
                    ROW2 = WM (I)
                    PR = RP (ROW2)
                    RDEG2 = II (PR+1)
C                   among those numerically acceptable rows of least
C                   (upper bound) degree, select the row with the
C                   lowest row index
                    BETTER = RDEG2 .LT. RDEG .OR.
     $                      (RDEG2 .EQ. RDEG .AND. ROW2 .LT. ROW)
                    IF (BETTER) THEN 
                       APIV = ABS (XX (WXP-1+I))
                       IF (APIV .GE. TOLER .AND. APIV .GT. 0) THEN 
                          ROW = ROW2
                          RDEG = RDEG2
                          ROWPOS = I
                       ENDIF 
                    ENDIF 
230              CONTINUE 
              ENDIF 
           ENDIF 

C          -------------------------------------------------------------
C          deallocate workspace
C          -------------------------------------------------------------

           XHEAD = XHEAD - XS
           XUSE = XUSE - XS
           XNEED = XNEED - XS
C          done using XX (wxp...wxp+xs-1) ]

C          -------------------------------------------------------------
C          reset work vector
C          -------------------------------------------------------------

           DO 240 I = 1, CDEG1 
              WIR (WM (I)) = -1
240        CONTINUE 

C          -------------------------------------------------------------
C          check to see if a pivot column was found
C          -------------------------------------------------------------

           IF (RDEG .EQ. N+1) THEN 

C             ----------------------------------------------------------
C             no pivot found, column is zero
C             ----------------------------------------------------------

C             remove this singular column from any further pivot search
              CDEG = -(N+2)
              II (PC+1) = CDEG

           ELSE 

C             ----------------------------------------------------------
C             save a list of the columns searched (with nonzero degrees)
C             ----------------------------------------------------------

              SRCHED = SRCHED + 1
              II (PC+7) = SLIST
              SLIST = COL

C             ----------------------------------------------------------
C             check if this is the best pivot seen so far
C             ----------------------------------------------------------

C             compute the true Markowitz cost without scanning the row
C             Wm (1..cdeg) holds pivot column, including pivot row index
C             Wm (rowpos) contains the candidate pivot row index
              COST = (CDEG - 1) * (RDEG - 1)
              IF (PIVCOL .EQ. 0 .OR. COST .LT. BESTCO) THEN 
                 FFLEFC = CDEG
                 DO 250 I = 1, FFLEFC-1 
                    WPC (I) = WM (I)
250              CONTINUE 
C                remove the pivot row index from pivot column pattern
                 WPC (ROWPOS) = WM (FFLEFC)
                 PIVCOL = COL
                 PIVROW = ROW
                 BESTCO = COST
              ENDIF 
           ENDIF 

C          done using Wm (1..cdeg) for pattern ]
C          WpC (1..fflefc-1) holds pivot column (excl. pivot row index)

C          -------------------------------------------------------------
C          exit global pivot search if nsrch pivots have been searched
C          -------------------------------------------------------------

           IF (SRCHED .GE. NSRCH) THEN 
              GO TO 260
           ENDIF 

255     CONTINUE 
C       exit label for loop 255:
260     CONTINUE

C=======================================================================
C  Quit early if no pivot found (singular matrix detected)
C=======================================================================

        IF (PIVCOL .EQ. 0) THEN 
C          complete the column permutation vector in
C          WpC (n-npiv+1 ... n) in reverse order
           K = N - NPIV + 1
           DO 270 COL = 1, N 
              IF (CP (COL) .NE. 0) THEN 
C                this is a non-pivotal column
                 K = K - 1
                 WPC (K) = COL
                 CP (COL) = 0
              ENDIF 
270        CONTINUE 
C          complete the row permutation vector in
C          WpR (n-npiv+1 ... n) in reverse order
           K = N - NPIV + 1
           DO 280 ROW = 1, NDN 
              IF (ROW .GT. N) THEN 
C                this is an artificial frontal matrix
                 E = ROW
                 RP (E) = 0
              ELSE IF (RP (ROW) .NE. 0) THEN 
                 RLEN = WC (ROW)
                 IF (RLEN .GE. 0 .AND. RLEN .LE. N) THEN 
C                   this is a non-pivotal row
                    K = K - 1
                    WPR (K) = ROW
                    RP (ROW) = 0
                 ELSE IF (RLEN .NE. -(NDN+2)) THEN 
C                   this is an unassembled element: convert to LU
                    E = ROW
                    EP = RP (ROW)
                    WR (E) = -(NDN+2)
                    WC (E) = -(NDN+2)
                    FLUIP = II (EP)
                    RP (E) = FLUIP
                 ENDIF 
              ENDIF 
280        CONTINUE 
C          factorization is done, exit the main factorization loop:
           GO TO 2000
        ENDIF 

C=======================================================================
C  Place the non-pivotal columns searched back in degree lists
C=======================================================================

        DO 300 I = 1, SRCHED 
           COL = SLIST
           PC = CP (COL)
           SLIST = II (PC+7)
           IF (COL .NE. PIVCOL) THEN 
              CDEG = II (PC+1)
              CNXT = HEAD (CDEG)
              II (PC+7) = CNXT
              II (PC+8) = 0
              IF (CNXT .NE. 0) THEN 
                 II (CP (CNXT)+8) = COL
              ENDIF 
              HEAD (CDEG) = COL
              MINDEG = MIN (MINDEG, CDEG)
           ENDIF 
300     CONTINUE 

C=======================================================================
C  Construct pivot row pattern
C=======================================================================

C       At this point, WiR (1..n) = -1 and WiC (1..n) is -2 for
C       nonpivotal columns and -1 for pivotal columns.
C       WiC (WpR (1..fflefr+1)) is set to zero in the code below.  It
C       will be set to the proper offsets in do 775, once ffdimc is
C       known (offsets are dependent on ffdimc, which is dependent on
C       fflefr calculated below, and the memory allocation).

C       ----------------------------------------------------------------
C       assemble the elements in the element list
C       ----------------------------------------------------------------

        PR = RP (PIVROW)
        FFLEFR = 0
        REP = (PR+2)
        RELN = WR (PIVROW)
        DO 330 IP = REP, REP + 2*RELN - 2, 2 
           E = II (IP)
           EP = RP (E)
           IF (E .LE. N) THEN 
              FLUIP = II (EP)
              LUCP = (FLUIP + 7)
              LUDEGR = II (FLUIP+2)
              LUDEGC = II (FLUIP+3)
              LURP = LUCP + LUDEGC
C             split into two loops so that they both vectorize on a CRAY
              F1 = FFLEFR
              DO 310 P = LURP, LURP + LUDEGR - 1 
                 COL = II (P)
                 IF (COL .GT. 0) THEN 
                    IF (WIC (COL) .EQ. -2) THEN 
                       FFLEFR = FFLEFR + 1
                       WPR (FFLEFR) = COL
                    ENDIF 
                 ENDIF 
310           CONTINUE 
              DO 320 I = F1+1, FFLEFR 
                 WIC (WPR (I)) = 0
320           CONTINUE 
           ELSE 
C             this is an artifical element (a dense column)
              LURP = (EP+8)
              COL = II (LURP)
              IF (WIC (COL) .EQ. -2) THEN 
                 FFLEFR = FFLEFR + 1
                 WPR (FFLEFR) = COL
                 WIC (COL) = 0
              ENDIF 
           ENDIF 
330     CONTINUE 

C       ----------------------------------------------------------------
C       assemble the original entries in the pivot row
C       ----------------------------------------------------------------

        RSIZ = II (PR)
        RLEN = WC (PIVROW)
        DO 340 P = PR + RSIZ - RLEN, PR + RSIZ - 1 
           COL = II (P)
           IF (WIC (COL) .EQ. -2) THEN 
              FFLEFR = FFLEFR + 1
              WPR (FFLEFR) = COL
           ENDIF 
340     CONTINUE 
C       the exact degree of the pivot row is fflefr

C=======================================================================
C  Initialize the new frontal matrix
C=======================================================================

C       ffrow is the name of current frontal matrix
        FFROW = PIVROW
        E1 = PIVROW
        K = 1
        K0 = 0

        FFDIMR = MIN (KLEFT, INT (GRO * FFLEFR))
        FFDIMC = MIN (KLEFT, INT (GRO * FFLEFC))

        FMAXR = FFLEFR
        FMAXC = FFLEFC
        FFSIZE = FFDIMC * FFDIMR
        RSCAN = MAX (DSIZ, FFDIMR)

C       ----------------------------------------------------------------
C       compute the offsets for rows in the pivot column
C       and the offsets for columns in the pivot row
C       ----------------------------------------------------------------

        DO 350 I = 1, FFLEFC - 1 
           WIR (WPC (I)) = I - 1
350     CONTINUE 
        DO 360 I = 1, FFLEFR 
           WIC (WPR (I)) = (I - 1) * FFDIMC
360     CONTINUE 

C       ----------------------------------------------------------------
C       remove the pivot column index from the pivot row pattern
C       ----------------------------------------------------------------

        COL = WPR (FFLEFR)
        COLPOS = (WIC (PIVCOL)/FFDIMC)+1
        WPR (COLPOS) = COL
        WIC (COL) = WIC (PIVCOL)
        WIC (PIVCOL) = (FFDIMR - 1) * FFDIMC
        WIR (PIVROW) = FFDIMC - 1

C       ----------------------------------------------------------------
C       remove the pivot row/col from the nonzero count
C       ----------------------------------------------------------------

        FFLEFR = FFLEFR - 1
        FFLEFC = FFLEFC - 1

C       ----------------------------------------------------------------
C       allocate the working array, doing garbage collection if needed
C       also allocate space for a work vector of size ffdimc
C       ----------------------------------------------------------------

        IF (FFSIZE + FFDIMC .GT. XTAIL-XHEAD) THEN 
           INFO (15) = INFO (15) + 1
           CALL UMD2FG (XX, XSIZE, XHEAD, XTAIL, XUSE,
     $                  II, ISIZE, IHEAD, ITAIL, IUSE,
     $                  CP, RP, DN, N, ICNTL, WIR, WIC, WR, WC,
     $                  0, 0, 0, 0, .FALSE.,
     $                  PFREE, XFREE, MHEAD, MTAIL, SLOTS)
C          at this point, iuse = ineed and xuse = xneed
        ENDIF 

        FFXP = XHEAD
        XHEAD = XHEAD + FFSIZE
        WXP = XHEAD
        XHEAD = XHEAD + FFDIMC
        XUSE = XUSE + FFSIZE + FFDIMC
        XNEED = XNEED + FFSIZE + FFDIMC
        INFO (20) = MAX (INFO (20), XUSE)
        INFO (21) = MAX (INFO (21), XNEED)
        IF (XHEAD .GT. XTAIL) THEN 
C          error return, if not enough real memory:
           GO TO 9000
        ENDIF 

C       ----------------------------------------------------------------
C       get memory usage for next call to UMD2RF
C       ----------------------------------------------------------------

        XRUSE = XRUSE + FFSIZE
        XRMAX = MAX (XRMAX, XRUSE)

C       ----------------------------------------------------------------
C       zero the working array
C       ----------------------------------------------------------------

C       zero the contribution block:
        DO 380 J = 0, FFLEFR - 1 
           DO 370 I = 0, FFLEFC - 1 
              XX (FFXP + J*FFDIMC + I) = 0
370        CONTINUE 
380     CONTINUE 

C       zero the pivot row:
        DO 390 J = 0, FFLEFR - 1 
           XX (FFXP + J*FFDIMC + FFDIMC-1) = 0
390     CONTINUE 

C       zero the pivot column:
        DO 400 I = 0, FFLEFC - 1 
           XX (FFXP + (FFDIMR-1)*FFDIMC + I) = 0
400     CONTINUE 

C       zero the pivot entry itself:
        XX (FFXP + (FFDIMR-1)*FFDIMC + FFDIMC-1) = 0

C       ----------------------------------------------------------------
C       current workspace usage:
C       ----------------------------------------------------------------

C       WpC (1..fflefc):        holds the pivot column pattern
C                               (excluding the pivot row index)
C       WpC (fflefc+1 .. n-npiv):       unused
C       WpC (n-npiv+1 .. n):            pivot columns in reverse order
C
C       WpR (1..fflefr):        holds the pivot row pattern
C                               (excluding the pivot column index)
C       WpR (fflefr+1 .. n-npiv):       unused
C       WpR (n-npiv+1 .. n):            pivot rows in reverse order
C
C       C (1..ffdimr, 1..ffdimc):  space for the new frontal matrix.
C
C       C (i,j) is located at XX (ffxp+((i)-1)+((j)-1)*ffdimc)
C
C       WiR (row) >= 0 for each row in pivot column pattern.
C               offset into pattern is given by:
C               WiR (row) == offset - 1
C               Also, WiR (pivrow) is ffdimc-1, the offset in C of
C               the pivot row itself.
C               Otherwise, WiR (1..n) is -1
C
C       WiC (col) >= 0 for each col in pivot row pattern.
C               WiC (col) == (offset - 1) * ffdimc
C               Also, WiC (pivcol) is (ffdimr-1)*ffdimc,
C               the offset in C of the pivot column itself.
C               Otherwise, WiC (1..n) is -2 for nonpivotal columns,
C               and -1 for pivotal columns

C       ----------------------------------------------------------------
C       remove the columns affected by this element from degree lists
C       ----------------------------------------------------------------

        DO 410 J = 1, FFLEFR 
           PC = CP (WPR (J))
           CDEG = II (PC+1)
           IF (CDEG .GT. 0) THEN 
              CNXT = II (PC+7)
              CPRV = II (PC+8)
              IF (CNXT .NE. 0) THEN 
                 II (CP (CNXT)+8) = CPRV
              ENDIF 
              IF (CPRV .NE. 0) THEN 
                 II (CP (CPRV)+7) = CNXT
              ELSE 
                 HEAD (CDEG) = CNXT
              ENDIF 
           ENDIF 
410     CONTINUE 

C=======================================================================
C  Initialization of new frontal matrix is complete ]
C=======================================================================

C=======================================================================
C  Assemble and factorize the current frontal matrix [
C=======================================================================

C       for first pivot in frontal matrix, do all scans
        SCAN1 = 0
        SCAN2 = 0
        SCAN3 = 0
        SCAN4 = 0

        DO 1395 DUMMY3 = 1, N 
C       (this loop is not indented due to its length)

C=======================================================================
C  Degree update and numerical assembly [
C=======================================================================

        KLEFT1 = KLEFT - 1

C       ----------------------------------------------------------------
C       SCAN1:  scan the element lists of each row in the pivot col
C               and compute the external column degree for each frontal
C       ----------------------------------------------------------------

        ROW = PIVROW
        DO 440 J = SCAN1, FFLEFC 
           IF (J .NE. 0) THEN 
C             Get a row;  otherwise, scan the pivot row if j is zero.
              ROW = WPC (J)
           ENDIF 
           PR = RP (ROW)
           REP = (PR+2)
           RELN = WR (ROW)
CFPP$ NODEPCHK L
           DO 430 P = REP, REP + 2*RELN - 2, 2 
              E = II (P)
              IF (WC (E) .LT. W0) THEN 
C                this is the first time seen in either scan 1 or 2:
                 EP = RP (E)
                 FLEFTR = II (EP+5)
                 FLEFTC = II (EP+6)
                 WR (E) = FLEFTR + W0
                 WC (E) = FLEFTC + W0
              ENDIF 
              WC (E) = WC (E) - 1
430        CONTINUE 
440     CONTINUE 

C       ----------------------------------------------------------------
C       SCAN2:  scan the element lists of each col in the pivot row
C               and compute the external row degree for each frontal
C       ----------------------------------------------------------------

        COL = PIVCOL
        DO 460 J = SCAN2, FFLEFR 
           IF (J .NE. 0) THEN 
C             Get a col;  otherwise, scan the pivot col if j is zero.
              COL = WPR (J)
           ENDIF 
           PC = CP (COL)
           CELN = II (PC+5)
           CEP = (PC+9)
CFPP$ NODEPCHK L
           DO 450 P = CEP, CEP + 2*CELN - 2, 2 
              E = II (P)
              IF (WR (E) .LT. W0) THEN 
C                this is the first time seen in either scan 1 or 2:
                 EP = RP (E)
                 FLEFTR = II (EP+5)
                 FLEFTC = II (EP+6)
                 WR (E) = FLEFTR + W0
                 WC (E) = FLEFTC + W0
              ENDIF 
              WR (E) = WR (E) - 1
450        CONTINUE 
460     CONTINUE 

C       ----------------------------------------------------------------
C       SCAN3:  scan the element lists of each column in pivot row
C               do degree update for the columns
C               assemble effective Usons and LU-sons
C       ----------------------------------------------------------------

C       flag Usons in Wc (e) as scanned (all now unflagged) [
C       uses Wc (e) for the link list.  Wc (e) <= 0
C       means that e is in the list, the external column
C       degree is zero, and -(Wc (e)) is the next element in
C       the Uson list.

        COL = PIVCOL
        DO 700 JJ = SCAN3, FFLEFR 

C          -------------------------------------------------------------
C          assemble and update the degree of a column
C          -------------------------------------------------------------

           IF (JJ .NE. 0) THEN 
C             Get a col;  otherwise, scan the pivot col if jj is zero
              COL = WPR (JJ)
           ENDIF 

C          -------------------------------------------------------------
C          compute the degree, and partition the element list into
C          two parts.  The first part are not LUsons or Usons, and
C          are not assembled.  The second part is assembled.
C          -------------------------------------------------------------

           CDEG = 0
           DELN = 0
           PC = CP (COL)
           CEP = (PC+9)
           CELN = II (PC+5)
           IP2 = CEP + 2*CELN - 2
           XUDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
           DO 470 IP = CEP, IP2, 2 
              E = II (IP)
              IF (WC (E) .GT. W0) THEN 
C                this element cannot be assembled
                    CDEG = CDEG + (WC (E) - W0)
              ELSE 
C                delete this tuple from the element list
                 DELN = DELN + 1
                 WM (DELN) = IP
              ENDIF 
470       CONTINUE 

          IF (DELN .NE. 0) THEN 

C             ----------------------------------------------------------
C             move the deleted tuples to the end of the element list
C             ----------------------------------------------------------

              P2 = IP2
              DO 480 I = DELN, 1, -1 
                 E = II (WM (I)  )
                 F = II (WM (I)+1)
                 II (WM (I)  ) = II (P2  )
                 II (WM (I)+1) = II (P2+1)
                 II (P2  ) = E
                 II (P2+1) = F
                 P2 = P2 - 2
480           CONTINUE 

C             ----------------------------------------------------------
C             assemble from LUsons and Usons (the deleted tuples) 
C             ----------------------------------------------------------

              DO 670 IP = P2 + 2, IP2, 2 

C                -------------------------------------------------------
C                this is an LUson or Uson.  If fextc < 0 then this has
C                already been assembled.
C                -------------------------------------------------------

                 E = II (IP)
                 IF (WC (E) .LT. W0) THEN 
C                   go to next iteration if already assembled
                    GOTO 670
                 ENDIF 

C                -------------------------------------------------------
C                get scalar info, add son to list if not already there
C                -------------------------------------------------------

                 EP = RP (E)
                 FDIMC = II (EP+1)
                 FXP = II (EP+2)
                 FLEFTR = II (EP+5)
                 FLEFTC = II (EP+6)
                 IF (E .LE. N) THEN 
                    FLUIP = II (EP)
                    LUDEGR = II (FLUIP+2)
                    LUDEGC = II (FLUIP+3)
                    LUCP = (FLUIP + 7)
                    LURP = LUCP + LUDEGC
                    IF (WIR (E) .EQ. -1) THEN 
                       WIR (E) = SONLST - N - 2
                       SONLST = E
                       NSONS = NSONS + 1
                    ENDIF 
                 ELSE 
C                   an artificial frontal matrix
                    LUDEGR = 1
                    LUDEGC = FDIMC
                    LUCP = (EP+9)
                    LURP = (EP+8)
                 ENDIF 

C                -------------------------------------------------------
                 IF (WR (E) .EQ. W0) THEN 
C                this is an LUson - assemble an entire frontal matrix
C                -------------------------------------------------------

C                   ----------------------------------------------------
                    IF (LUDEGC .EQ. FLEFTC) THEN 
C                   no rows assembled out of this frontal yet
C                   ----------------------------------------------------

C                      compute the compressed column offset vector
C                      use Wm (1..ludegc for offsets) [
                       DO 490 I = 0, LUDEGC-1 
                          ROW2 = II (LUCP+I)
                          WM (I+1) = WIR (ROW2)
490                    CONTINUE 

C                      -------------------------------------------------
                       IF (LUDEGR .EQ. FLEFTR) THEN 
C                      no rows or cols assembled out of frontal yet
C                      -------------------------------------------------

                          DO 510 J = 0, LUDEGR-1 
                             COL2 = II (LURP+J)
                             XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                             DO 500 I = 0, LUDEGC-1 
                                XX (XDP + WM (I+1)) =
     $                          XX (XDP + WM (I+1)) +
     $                          XX (FXP + J*FDIMC + I)
500                          CONTINUE 
510                       CONTINUE 

C                      -------------------------------------------------
                       ELSE 
C                      only cols have been assembled out of frontal
C                      -------------------------------------------------

                          DO 530 J = 0, LUDEGR-1 
                             COL2 = II (LURP+J)
                             IF (COL2 .GT. 0) THEN 
                                XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                                DO 520 I = 0, LUDEGC-1 
                                   XX (XDP + WM (I+1)) =
     $                             XX (XDP + WM (I+1)) +
     $                             XX (FXP + J*FDIMC + I)
520                             CONTINUE 
                             ENDIF 
530                       CONTINUE 
                       ENDIF 
C                      done using Wm (1..ludegc for offsets) ]

C                   ----------------------------------------------------
                    ELSE 
C                   some rows have been assembled out of this frontal
C                   ----------------------------------------------------

C                      compute the compressed column offset vector
C                      use Wm (1..ludegc for offsets) [
                       DEGC = 0
                       DO 540 I = 0, LUDEGC-1 
                          ROW2 = II (LUCP+I)
                          IF (ROW2 .GT. 0) THEN 
                             DEGC = DEGC + 1
                             WJ (DEGC) = I
                             WM (DEGC) = WIR (ROW2)
                          ENDIF 
540                    CONTINUE 

C                      -------------------------------------------------
                       IF (LUDEGR .EQ. FLEFTR) THEN 
C                      only rows assembled out of this frontal
C                      -------------------------------------------------

                          DO 560 J = 0, LUDEGR-1 
                             COL2 = II (LURP+J)
                             XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                             DO 550 I = 1, DEGC 
                                XX (XDP + WM (I)) =
     $                          XX (XDP + WM (I)) +
     $                          XX (FXP + J*FDIMC + WJ (I))
550                          CONTINUE 
560                       CONTINUE 

C                      -------------------------------------------------
                       ELSE 
C                      both rows and columns assembled out of frontal
C                      -------------------------------------------------

                          DO 580 J = 0, LUDEGR-1 
                             COL2 = II (LURP+J)
                             IF (COL2 .GT. 0) THEN 
                                XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                                DO 570 I = 1, DEGC 
                                   XX (XDP + WM (I)) =
     $                             XX (XDP + WM (I)) +
     $                             XX (FXP + J*FDIMC + WJ (I))
570                             CONTINUE 
                             ENDIF 
580                       CONTINUE 
                       ENDIF 
C                      done using Wm (1..ludegc for offsets) ]
                    ENDIF 

C                   ----------------------------------------------------
C                   deallocate the LUson frontal matrix
C                   ----------------------------------------------------

                    WR (E) = -(NDN+2)
                    WC (E) = -(NDN+2)
                    IF (E .LE. N) THEN 
                       RP (E) = FLUIP
                       II (EP) = FSCAL
                       INEED = INEED - FSCAL
                    ELSE 
                       RP (E) = 0
                       II (EP) = FDIMC + CSCAL
                       INEED = INEED - (FDIMC + CSCAL)
                    ENDIF 
                    II (EP+1) = -1
                    II (EP+6) = 0
                    MPREV = II (EP+4)
                    MNEXT = II (EP+3)
                    XNEED = XNEED - LUDEGR*LUDEGC
                    IF (MNEXT .NE. 0 .AND. II (MNEXT+6) .EQ. 0) THEN 
C                      next block is free - delete it
                       MNEXT = II (MNEXT+3)
                       II (EP+3) = MNEXT
                       IF (MNEXT .NE. 0) THEN 
                          II (MNEXT+4) = EP
                       ELSE 
                          MTAIL = EP
                       ENDIF 
                    ENDIF 
                    IF (MPREV .NE. 0 .AND. II (MPREV+6) .EQ. 0) THEN 
C                      previous block is free - delete it
                       II (EP+2) = II (MPREV+2)
                       MPREV = II (MPREV+4)
                       II (EP+4) = MPREV
                       IF (MPREV .NE. 0) THEN 
                          II (MPREV+3) = EP
                       ELSE 
                          MHEAD = EP
                       ENDIF 
                    ENDIF 
C                   get the size of the freed block
                    IF (MNEXT .NE. 0) THEN 
                       XS = II (MNEXT+2) - II (EP+2)
                    ELSE 
                       XS = FFXP - II (EP+2)
                    ENDIF 
                    IF (XS .GT. XFREE) THEN 
C                      keep track of the largest free block
                       XFREE = XS
                       PFREE = EP
                    ENDIF 

C                   ----------------------------------------------------
C                   get memory usage for next call to UMD2RF
C                   ----------------------------------------------------

                    XRUSE = XRUSE - LUDEGR*LUDEGC

C                -------------------------------------------------------
                 ELSE IF (WR (E) - W0 .LE. FLEFTR/2) THEN 
C                this is a Uson - assemble all possible columns
C                -------------------------------------------------------

C                   ----------------------------------------------------
C                   add to Uson list - to be cleared just after scan 3
C                   ----------------------------------------------------

                    WC (E) = -USONS
                    USONS = E

C                   ----------------------------------------------------
                    IF (LUDEGC .EQ. FLEFTC) THEN 
C                   no rows assembled out of this Uson frontal yet
C                   ----------------------------------------------------

C                      compute the compressed column offset vector
C                      use Wm (1..ludegc for offsets)
                       DO 590 I = 0, LUDEGC-1 
                          ROW2 = II (LUCP+I)
                          WM (I+1) = WIR (ROW2)
590                    CONTINUE 

                       DO 610 J = 0, LUDEGR-1 
                          COL2 = II (LURP+J)
                          IF (COL2 .GT. 0) THEN 
                             IF (WIC (COL2) .GE. 0) THEN 
                                XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                                DO 600 I = 0, LUDEGC-1 
                                   XX (XDP + WM (I+1)) =
     $                             XX (XDP + WM (I+1)) +
     $                             XX (FXP + J*FDIMC + I)
600                             CONTINUE 
C                               flag this column as assembled from Uson
                                II (LURP+J) = -COL2
                             ENDIF 
                          ENDIF 
610                    CONTINUE 

C                   ----------------------------------------------------
                    ELSE 
C                   some rows already assembled out of this Uson frontal
C                   ----------------------------------------------------

C                      compute the compressed column offset vector
C                      use Wm (1..ludegc for offsets)
                       DEGC = 0
                       DO 620 I = 0, LUDEGC-1 
                          ROW2 = II (LUCP+I)
                          IF (ROW2 .GT. 0) THEN 
                             DEGC = DEGC + 1
                             WJ (DEGC) = I
                             WM (DEGC) = WIR (ROW2)
                          ENDIF 
620                    CONTINUE 

                       DO 640 J = 0, LUDEGR-1 
                          COL2 = II (LURP+J)
                          IF (COL2 .GT. 0) THEN 
                             IF (WIC (COL2) .GE. 0) THEN 
                                XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                                DO 630 I = 1, DEGC 
                                   XX (XDP + WM (I)) =
     $                             XX (XDP + WM (I)) +
     $                             XX (FXP + J*FDIMC + WJ (I))
630                             CONTINUE 
C                               flag this column as assembled from Uson
                                II (LURP+J) = -COL2
                             ENDIF 
                          ENDIF 
640                    CONTINUE 

                    ENDIF 

                    FLEFTR = WR (E) - W0
                    II (EP+5) = FLEFTR

C                -------------------------------------------------------
                 ELSE 
C                this is a Uson - assemble just one column
C                -------------------------------------------------------

C                   get the offset, f, from the (e,f) tuple
                    F = II (IP+1)

C                   ----------------------------------------------------
                    IF (LUDEGC .EQ. FLEFTC) THEN 
C                   no rows assembled out of this Uson yet
C                   ----------------------------------------------------

CFPP$ NODEPCHK L
                       DO 650 I = 0, LUDEGC-1 
                          ROW2 = II (LUCP+I)
                          XX (XUDP + WIR (ROW2)) =
     $                    XX (XUDP + WIR (ROW2)) + 
     $                    XX (FXP + F*FDIMC + I)
650                    CONTINUE 

C                   ----------------------------------------------------
                    ELSE 
C                   some rows already assembled out of this Uson
C                   ----------------------------------------------------

CFPP$ NODEPCHK L
                       DO 660 I = 0, LUDEGC-1 
                          ROW2 = II (LUCP+I)
                          IF (ROW2 .GT. 0) THEN 
                             XX (XUDP + WIR (ROW2)) =
     $                       XX (XUDP + WIR (ROW2)) +
     $                       XX (FXP + F*FDIMC + I)
                          ENDIF 
660                    CONTINUE 
                    ENDIF 

C                   ----------------------------------------------------
C                   decrement count of unassembled cols in frontal
C                   ----------------------------------------------------

                    II (EP+5) = FLEFTR - 1
C                   flag the column as assembled from the Uson
                    II (LURP+F) = -COL
                 ENDIF 

670           CONTINUE 

C             ----------------------------------------------------------
C             update the count of (e,f) tuples in the element list
C             ----------------------------------------------------------

              II (PC+5) = II (PC+5) - DELN
              INEED = INEED - 2*DELN
           ENDIF 

C          -------------------------------------------------------------
C          assemble the original column and update count of entries
C          -------------------------------------------------------------

           CLEN = II (PC+6)
           IF (CLEN .GT. 0) THEN 
              CSIZ = II (PC)
              IP = PC + CSIZ - CLEN
              DLEN = 0
CFPP$ NODEPCHK L
              DO 680 I = 0, CLEN - 1 
                 ROW = II (IP+I)
                 IF (WIR (ROW) .GE. 0) THEN 
C                   this entry can be assembled and deleted
                    DLEN = DLEN + 1
                    WM (DLEN) = I
                 ENDIF 
680           CONTINUE 
              IF (DLEN .NE. 0) THEN 
                 CXP = II (PC+2)
                 DO 690 J = 1, DLEN 
                    I = WM (J)
                    ROW = II (IP+I)
C                   assemble the entry
                    XX (XUDP + WIR (ROW)) =
     $              XX (XUDP + WIR (ROW)) + XX (CXP+I)
C                   and delete the entry
                    II (IP +I) = II (IP +J-1)
                    XX (CXP+I) = XX (CXP+J-1)
690              CONTINUE 
                 CLEN = CLEN - DLEN
                 CXP = CXP + DLEN
                 INEED = INEED - DLEN
                 XNEED = XNEED - DLEN
                 II (PC+6) = CLEN
                 IF (CLEN .NE. 0) THEN 
                    II (PC+2) = CXP
                 ELSE 
C                   deallocate the real portion of the column:
                    MPREV = II (PC+4)
                    MNEXT = II (PC+3)
                    IF (MNEXT .NE. 0 .AND. II (MNEXT+6) .EQ. 0) THEN
C                      next block is free - delete it
                       MNEXT = II (MNEXT+3)
                       II (PC+3) = MNEXT
                       IF (MNEXT .NE. 0) THEN 
                          II (MNEXT+4) = PC
                       ELSE 
                          MTAIL = PC
                       ENDIF 
                    ENDIF 
                    IF (MPREV .NE. 0 .AND. II (MPREV+6) .EQ. 0) THEN
C                      previous block is free - delete it
                       II (PC+2) = II (MPREV+2)
                       MPREV = II (MPREV+4)
                       II (PC+4) = MPREV
                       IF (MPREV .NE. 0) THEN 
                          II (MPREV+3) = PC
                       ELSE 
                          MHEAD = PC
                       ENDIF 
                    ENDIF 
                    IF (PC .EQ. MHEAD) THEN 
C                      adjust the start of the block if this is head
                       II (PC+2) = 1
                    ENDIF 
C                   get the size of the freed block
                    IF (MNEXT .NE. 0) THEN 
                       XS = II (MNEXT+2) - II (PC+2)
                    ELSE 
                       XS = FFXP - II (PC+2)
                    ENDIF 
                    IF (XS .GT. XFREE) THEN 
C                      keep track of the largest free block
                       XFREE = XS
                       PFREE = PC
                    ENDIF 
                 ENDIF 
              ENDIF 
              CDEG = CDEG + CLEN
           ENDIF 

C          -------------------------------------------------------------
C          compute the upper bound degree - excluding current front
C          -------------------------------------------------------------

           CDEG2 = II (PC+1)
           CDEG = MIN (KLEFT1 - FFLEFC, CDEG2, CDEG)
           II (PC+1) = CDEG

700     CONTINUE 

C       ----------------------------------------------------------------
C       SCAN-3 wrap-up:  remove flags from assembled Usons
C       ----------------------------------------------------------------

C       while (usons .ne. ndn+1) do
710     CONTINUE
        IF (USONS .NE. NDN+1) THEN 
           NEXT = -WC (USONS)
           WC (USONS) = W0
           USONS = NEXT
C       end while:
        GOTO 710
        ENDIF 
C       done un-flagging usons, all now unflagged in Wc (e) ]

C       ----------------------------------------------------------------
C       SCAN4:  scan element lists of each row in the pivot column
C               do degree update for the rows
C               assemble effective Lsons
C       ----------------------------------------------------------------

C       flag Lsons in Wr (e) (all are now unflagged) [
C       uses Wr (e) for the link list.  Wr (e) <= 0 means
C       that e is in the list, the external row degree is zero, and
C       -(Wr (e)) is the next element in the Lson list.

        ROW = PIVROW
        DO 840 JJ = SCAN4, FFLEFC 

C          -------------------------------------------------------------
C          assemble and update the degree of a row
C          -------------------------------------------------------------

           IF (JJ .NE. 0) THEN 
C             Get a row;  otherwise, scan the pivot row if jj is zero
              ROW = WPC (JJ)
           ENDIF 

C          -------------------------------------------------------------
C          compute the degree, and partition the element list into
C          two parts.  The first part are not LUsons or Lsons, and
C          are not assembled.  The second part is assembled.
C          -------------------------------------------------------------

           RDEG = 0
           DELN = 0
           PR = RP (ROW)
           REP = (PR+2)
           RELN = WR (ROW)
           IP2 = REP + 2*RELN - 2
CFPP$ NODEPCHK L
           DO 720 IP = REP, IP2, 2 
              E = II (IP)
              IF (WR (E) .GT. W0) THEN 
                 RDEG = RDEG + (WR (E) - W0)
              ELSE 
                 DELN = DELN + 1
                 WM (DELN) = IP
              ENDIF 
720        CONTINUE 
                    
           IF (DELN .NE. 0) THEN 

C             ----------------------------------------------------------
C             move the deleted tuples to the end of the element list
C             ----------------------------------------------------------

              P2 = IP2
              DO 730 I = DELN, 1, -1 
                 E = II (WM (I)  )
                 F = II (WM (I)+1)
                 II (WM (I)  ) = II (P2  )
                 II (WM (I)+1) = II (P2+1)
                 II (P2  ) = E
                 II (P2+1) = F
                 P2 = P2 - 2
730           CONTINUE 

C             ----------------------------------------------------------
C             assemble from Lsons (the deleted tuples) 
C             ----------------------------------------------------------

              DO 810 IP = P2 + 2, IP2, 2 

C                -------------------------------------------------------
C                this is an LUson or Lson.  If fextr < 0 then this has
C                already been assembled.  All LUsons have already been
C                assembled (in SCAN3, above).
C                -------------------------------------------------------

                 E = II (IP)
                 IF (WR (E) .LT. W0) THEN 
C                   go to next iteration if already assembled
                    GOTO 810
                 ENDIF 

C                -------------------------------------------------------
C                get scalar info, add to son list if not already there
C                -------------------------------------------------------

                 EP = RP (E)
                 FDIMC = II (EP+1)
                 FXP = II (EP+2)
                 FLEFTR = II (EP+5)
                 FLEFTC = II (EP+6)
                 IF (E .LE. N) THEN 
                    FLUIP = II (EP)
                    LUDEGR = II (FLUIP+2)
                    LUDEGC = II (FLUIP+3)
                    LUCP = (FLUIP + 7)
                    LURP = LUCP + LUDEGC
                    IF (WIR (E) .EQ. -1) THEN 
                       WIR (E) = SONLST - N - 2
                       SONLST = E
                       NSONS = NSONS + 1
                    ENDIF 
                 ELSE 
C                   an artificial frontal matrix
                    LUDEGR = 1
                    LUDEGC = FDIMC
                    LUCP = (EP+9)
                    LURP = (EP+8)
                 ENDIF 

C                -------------------------------------------------------
                 IF (WC (E) - W0 .LE. FLEFTC/2) THEN 
C                this is an Lson - assemble all possible rows
C                -------------------------------------------------------

C                   ----------------------------------------------------
C                   add to Lson list - to be cleared just after scan 4
C                   ----------------------------------------------------

                    WR (E) = -LSONS
                    LSONS = E

C                   compute the compressed column offset vector
C                   use Wm (1..ludegc for offsets) [
                    DEGC = 0
                    DO 740 I = 0, LUDEGC-1 
                       ROW2 = II (LUCP+I)
                       IF (ROW2 .GT. 0) THEN 
                          IF (WIR (ROW2) .GE. 0) THEN 
C                            this row will be assembled in loop below
                             DEGC = DEGC + 1
                             WJ (DEGC) = I
                             WM (DEGC) = WIR (ROW2)
C                            flag the row as assembled from the Lson
                             II (LUCP+I) = -ROW2
                          ENDIF 
                       ENDIF 
740                 CONTINUE 

C                   ----------------------------------------------------
                    IF (LUDEGR .EQ. FLEFTR) THEN 
C                   no columns assembled out this Lson yet
C                   ----------------------------------------------------

                       DO 760 J = 0, LUDEGR-1 
                          COL2 = II (LURP+J)
                          XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                          DO 750 I = 1, DEGC 
                             XX (XDP + WM (I)) =
     $                       XX (XDP + WM (I)) +
     $                       XX (FXP + J*FDIMC + WJ (I))
750                       CONTINUE 
760                    CONTINUE 

C                   ----------------------------------------------------
                    ELSE 
C                   some columns already assembled out of this Lson
C                   ----------------------------------------------------

                       DO 780 J = 0, LUDEGR-1 
                          COL2 = II (LURP+J)
                          IF (COL2 .GT. 0) THEN 
                             XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                             DO 770 I = 1, DEGC 
                                XX (XDP + WM (I)) =
     $                          XX (XDP + WM (I)) +
     $                          XX (FXP + J*FDIMC + WJ (I))
770                          CONTINUE 
                          ENDIF 
780                    CONTINUE 
                    ENDIF 

C                   done using Wm (1..ludegc for offsets) ]
                    FLEFTC = WC (E) - W0
                    II (EP+6) = FLEFTC

C                -------------------------------------------------------
                 ELSE 
C                this is an Lson - assemble just one row
C                -------------------------------------------------------

                    XLDP = FFXP + WIR (ROW)
C                   get the offset, f, from the (e,f) tuple
                    F = II (IP+1)

C                   ----------------------------------------------------
                    IF (LUDEGR .EQ. FLEFTR) THEN 
C                   no columns assembled out this Lson yet
C                   ----------------------------------------------------

CFPP$ NODEPCHK L
                       DO 790 J = 0, LUDEGR-1 
                          COL2 = II (LURP+J)
                          XX (XLDP + WIC (COL2)) =
     $                    XX (XLDP + WIC (COL2)) +
     $                    XX (FXP + J*FDIMC + F)
790                    CONTINUE 

C                   ----------------------------------------------------
                    ELSE 
C                   some columns already assembled out of this Lson
C                   ----------------------------------------------------

CFPP$ NODEPCHK L
                       DO 800 J = 0, LUDEGR-1 
                          COL2 = II (LURP+J)
                          IF (COL2 .GT. 0) THEN 
                             XX (XLDP + WIC (COL2)) =
     $                       XX (XLDP + WIC (COL2)) +
     $                       XX (FXP + J*FDIMC + F)
                          ENDIF 
800                    CONTINUE 
                    ENDIF 

                    II (EP+6) = FLEFTC - 1
C                   flag the row as assembled from the Lson
                    II (LUCP+F) = -ROW
                 ENDIF 

810           CONTINUE 

C             ----------------------------------------------------------
C             update the count of (e,f) tuples in the element list
C             ----------------------------------------------------------

              WR (ROW) = WR (ROW) - DELN
              INEED = INEED - 2*DELN
           ENDIF 

C          -------------------------------------------------------------
C          assemble the original row and update count of entries
C          -------------------------------------------------------------

           RLEN = WC (ROW)
           IF (RLEN .GT. 0) THEN 
C             do not scan a very long row:
              IF (RLEN .LE. RSCAN) THEN 
                 RSIZ = II (PR)
                 IP = PR + RSIZ - RLEN
                 DLEN = 0
CFPP$ NODEPCHK L
                 DO 820 P = IP, IP + RLEN - 1 
                    COL = II (P)
                    IF (WIC (COL) .NE. -2) THEN 
C                      this entry can be assembled and deleted
C                      if WiC (col) = -1, it is an older pivot col,
C                      otherwise (>=0) it is in the current element
                       DLEN = DLEN + 1
                       WM (DLEN) = P
                    ENDIF 
820              CONTINUE 
                 IF (DLEN .NE. 0) THEN 
                    DO 830 J = 1, DLEN 
C                      delete the entry
                       II (WM (J)) = II (IP+J-1)
830                 CONTINUE 
                    RLEN = RLEN - DLEN
                    INEED = INEED - DLEN
                    WC (ROW) = RLEN
                 ENDIF 
              ENDIF 
              RDEG = RDEG + RLEN
           ENDIF 

C          -------------------------------------------------------------
C          compute the upper bound degree - excluding current front
C          -------------------------------------------------------------

           RDEG2 = II (PR+1)
           RDEG = MIN (KLEFT1 - FFLEFR, RDEG2, RDEG)
           II (PR+1) = RDEG

840     CONTINUE 

C       ----------------------------------------------------------------
C       SCAN-4 wrap-up:  remove flags from assembled Lsons
C       ----------------------------------------------------------------

C       while (lsons .ne. ndn+1) do
850     CONTINUE
        IF (LSONS .NE. NDN+1) THEN 
           NEXT = -WR (LSONS)
           WR (LSONS) = W0
           LSONS = NEXT
C       end while:
        GOTO 850
        ENDIF 
C       done un-flagging Lsons, all now unflagged in Wr (e) ]

C=======================================================================
C  Degree update and numerical assemble is complete ]
C=======================================================================

C=======================================================================
C  Factorize frontal matrix until next pivot extends it [
C=======================================================================

        DO 1324 DUMMY4 = 1, N 
C       (this loop is not indented due to its length)

C       ----------------------------------------------------------------
C       Wc (e) = fextc+w0, where fextc is the external column
C               degree for each element (ep = Rp (e)) appearing in
C               the element lists for each row in the pivot column.
C               if Wc (e) < w0, then fextc is defined as II (ep+6)
C
C       Wr (e) = fextr+w0, where fextr is the external row
C               degree for each element (ep = Rp (e)) appearing in
C               the element lists for each column in the pivot row
C               if Wr (e) < w0, then fextr is defined as II (ep+5)
C
C       WiR (row) >= 0 for each row in pivot column pattern.
C               offset into pattern is given by:
C               WiR (row) == offset - 1
C               WiR (pivrow) is the offset of the latest pivot row
C
C       WiC (col) >= 0 for each col in pivot row pattern.
C               WiC (col) == (offset - 1) * ffdimc
C               WiC (pivcol) is the offset of the latest pivot column
C
C       WpR (1..fflefr) is the pivot row pattern (excl pivot cols)
C       WpC (1..fflefc) is the pivot col pattern (excl pivot rows)
C       ----------------------------------------------------------------

C=======================================================================
C  Divide pivot column by pivot
C=======================================================================

C       k-th pivot in frontal matrix located in C(ffdimc-k+1,ffdimr-k+1)
        XDP = FFXP + (FFDIMR - K) * FFDIMC
        TEMP = XX (XDP + FFDIMC - K)

C       divide C(1:fflefc,ffdimr-k+1) by pivot value
        TEMP = 1 / TEMP
        DO 870 P = XDP, XDP + FFLEFC-1 
           XX (P) = XX (P) * TEMP
870     CONTINUE 
C       count this as a call to the Level-1 BLAS:
        RINFO (4) = RINFO (4) + (FFLEFC)

C=======================================================================
C  A pivot step is complete
C=======================================================================

        KLEFT = KLEFT - 1
        NPIV = NPIV + 1
        INFO (17) = INFO (17) + 1

C       ----------------------------------------------------------------
C       the pivot column is fully assembled and scaled, and is now the
C       (npiv)-th column of L. The pivot row is the (npiv)-th row of U.
C       ----------------------------------------------------------------

        WPR (N-NPIV+1) = PIVROW
        WPC (N-NPIV+1) = PIVCOL
        WIR (PIVROW) = -1
        WIC (PIVCOL) = -1

C       ----------------------------------------------------------------
C       deallocate the pivot row and pivot column
C       ----------------------------------------------------------------

        RLEN = WC (PIVROW)
        INEED = INEED - CSCAL - RSCAL - RLEN
        PR = RP (PIVROW)
        PC = CP (PIVCOL)
        II (PR+1) = -1
        II (PC+1) = -1
        RP (PIVROW) = 0
        CP (PIVCOL) = 0

C=======================================================================
C  Local search for next pivot within current frontal matrix [
C=======================================================================

        FEDEGC = FFLEFC
        FEDEGR = FFLEFR
        PFOUND = .FALSE.
        OKCOL = FFLEFC .GT. 0
        OKROW = .FALSE.

C       ----------------------------------------------------------------
C       find column of minimum degree in current frontal row pattern
C       ----------------------------------------------------------------

C       among those columns of least (upper bound) degree, select the
C       column with the lowest column index
        IF (OKCOL) THEN 
           COLPOS = 0
           PIVCOL = N+1
           CDEG = N+1
C          can this be vectorized?  This is the most intensive
C          non-vector loop.
           DO 880 J = 1, FFLEFR 
              COL = WPR (J)
              PC = CP (COL)
              CDEG2 = II (PC+1)
              BETTER = CDEG2 .GE. 0 .AND.
     $                (CDEG2 .LT. CDEG .OR.
     $                (CDEG2 .EQ. CDEG .AND. COL .LT. PIVCOL))
              IF (BETTER) THEN 
                 CDEG = CDEG2
                 COLPOS = J
                 PIVCOL = COL
              ENDIF 
880        CONTINUE 
           OKCOL = COLPOS .NE. 0
        ENDIF 

C=======================================================================
C  Assemble candidate pivot column in temporary workspace
C=======================================================================

        IF (OKCOL) THEN 
           PC = CP (PIVCOL)
           CLEN = II (PC+6)
           OKCOL = FEDEGC + CLEN .LE. FFDIMC
        ENDIF 

        IF (OKCOL) THEN 

C          -------------------------------------------------------------
C          copy candidate column from current frontal matrix into
C          work vector XX (wxp ... wxp+ffdimc-1) [
C          -------------------------------------------------------------

           P = FFXP + (COLPOS - 1) * FFDIMC - 1
CFPP$ NODEPCHK L
           DO 890 I = 1, FFLEFC 
              XX (WXP-1+I) = XX (P+I)
890        CONTINUE 

C          -------------------------------------------------------------
C          update candidate column with previous pivots in this front
C          -------------------------------------------------------------

           IF (K-K0 .GT. 0 .AND. FFLEFC .NE. 0) THEN 
              CALL DGEMV ('N', FFLEFC, K-K0,
     $          -ONE, XX (FFXP + (FFDIMR - K) * FFDIMC)        ,FFDIMC,
     $                XX (FFXP + (COLPOS - 1) * FFDIMC + FFDIMC - K), 1,
     $           ONE, XX (WXP)                                      , 1)
              TMP = FFLEFC
              TMP = TMP * (K-K0)
              RINFO (3) = RINFO (3) + 2.0*(TMP)
           ENDIF 

C          -------------------------------------------------------------
C          Compute extended pivot column in XX (wxp..wxp-1+fedegc).
C          Pattern of pivot column is placed in WpC (1..fedegc)
C          -------------------------------------------------------------

C          assemble the elements in the element list
           CEP = (PC+9)
           CELN = II (PC+5)
           DO 930 IP = CEP, CEP + 2*CELN - 2, 2 
              E = II (IP)
              F = II (IP+1)
              EP = RP (E)
              FLEFTC = II (EP+6)
              FDIMC = II (EP+1)
              FXP = II (EP+2)
              IF (E .LE. N) THEN 
                 FLUIP = II (EP)
                 LUCP = (FLUIP + 7)
                 LUDEGC = II (FLUIP+3)
              ELSE 
                 LUCP = (EP+9)
                 LUDEGC = FDIMC
              ENDIF 
              XP = FXP + F * FDIMC
C             split into 3 loops so that they all vectorize on a CRAY
              F1 = FEDEGC
              DO 900 P = LUCP, LUCP + LUDEGC - 1 
                 ROW = II (P)
                 IF (ROW .GT. 0) THEN 
                    IF (WIR (ROW) .LT. 0) THEN 
                       F1 = F1 + 1
                       WPC (F1) = ROW
                    ENDIF 
                 ENDIF 
900           CONTINUE 
              OKCOL = F1 + CLEN .LE. FFDIMC
              IF (.NOT. OKCOL) THEN 
C                exit out of loop if column too long:
                 GO TO 940
              ENDIF 
              DO 910 I = FEDEGC+1, F1 
                 ROW = WPC (I)
                 WIR (ROW) = I - 1
                 XX (WXP-1+I) = 0
910           CONTINUE 
              FEDEGC = F1
CFPP$ NODEPCHK L
              DO 920 J = 0, LUDEGC - 1 
                 ROW = II (LUCP+J)
                 IF (ROW .GT. 0) THEN 
                    XX (WXP + WIR (ROW)) = 
     $              XX (WXP + WIR (ROW)) + XX (XP+J)
                 ENDIF 
920           CONTINUE 
930        CONTINUE 
C          loop exit label:
940        CONTINUE
        ENDIF 

C=======================================================================
C  Find candidate pivot row - unless candidate pivot column is too long
C=======================================================================

        IF (OKCOL) THEN 

C          -------------------------------------------------------------
C          assemble the original entries in the column
C          -------------------------------------------------------------

           CSIZ = II (PC)
           IP = PC + CSIZ - CLEN
           CXP = II (PC+2)
CFPP$ NODEPCHK L
           DO 950 I = 0, CLEN - 1 
              ROW = II (IP+I)
              WIR (ROW) = FEDEGC + I
              WPC (FEDEGC+1+I) = ROW
              XX  (WXP+FEDEGC+I) = XX (CXP+I)
950        CONTINUE 
           FEDEGC = FEDEGC + CLEN

C          -------------------------------------------------------------
C          update degree of candidate column - excluding current front
C          -------------------------------------------------------------

           CDEG = FEDEGC - FFLEFC
           II (PC+1) = CDEG

C          -------------------------------------------------------------
C          find the maximum absolute value in the column
C          -------------------------------------------------------------

           MAXVAL = ABS (XX (WXP-1 + IDAMAX (FEDEGC, XX (WXP), 1)))
           RINFO (3) = RINFO (3) + (FEDEGC)
           TOLER = RELPT * MAXVAL
           RDEG = N+1

C          -------------------------------------------------------------
C          look for the best possible pivot row in this column
C          -------------------------------------------------------------

           IF (MAXVAL .GT. 0) THEN 
              IF (SYMSRC) THEN 
C                prefer symmetric pivots, if numerically acceptable
                 PIVROW = PIVCOL
                 ROWPOS = WIR (PIVROW) + 1
                 IF (ROWPOS .GT. 0 .AND. ROWPOS .LE. FFLEFC) THEN 
C                   diagonal entry exists in the column pattern
C                   also within the current frontal matrix
                    APIV = ABS (XX (WXP-1+ROWPOS))
                    IF (APIV .GE. TOLER .AND. APIV .GT. 0) THEN 
C                      diagonal entry is numerically acceptable
                       PR = RP (PIVROW)
                       RDEG = II (PR+1)
                    ENDIF 
                 ENDIF 
              ENDIF 
              IF (RDEG .EQ. N+1) THEN 
C                Continue searching - no diagonal found or sought for.
C                Minimize row degree subject to abs(value) constraints.
                 PIVROW = N+1
                 DO 960 I = 1, FFLEFC 
                    ROW2 = WPC (I)
                    PR = RP (ROW2)
                    RDEG2 = II (PR+1)
C                   among those numerically acceptable rows of least
C                   (upper bound) degree, select the row with the
C                   lowest row index
                    BETTER = RDEG2 .LT. RDEG .OR.
     $                      (RDEG2 .EQ. RDEG .AND. ROW2 .LT. PIVROW)
                    IF (BETTER) THEN 
                       APIV = ABS (XX (WXP-1+I))
                       IF (APIV .GE. TOLER .AND. APIV .GT. 0) THEN 
                          PIVROW = ROW2
                          RDEG = RDEG2
                          ROWPOS = I
                       ENDIF 
                    ENDIF 
960              CONTINUE 
              ENDIF 
           ELSE 
C             remove this column from any further pivot search
              CDEG = -(N+2)
              II (PC+1) = CDEG
           ENDIF 
           OKROW = RDEG .NE. N+1
        ENDIF 

C       done using XX (wxp...wxp+ffdimc-1) ]

C=======================================================================
C  If found, construct candidate pivot row pattern
C=======================================================================

        IF (OKROW) THEN 

C          -------------------------------------------------------------
C          assemble the elements in the element list
C          -------------------------------------------------------------

           PR = RP (PIVROW)
           REP = (PR+2)
           RELN = WR (PIVROW)
           DO 990 IP = REP, REP + 2*RELN - 2, 2 
              E = II (IP)
              EP = RP (E)
              IF (E .LE. N) THEN 
                 FLUIP = II (EP)
                 LUCP = (FLUIP + 7)
                 LUDEGR = II (FLUIP+2)
                 LUDEGC = II (FLUIP+3)
                 LURP = LUCP + LUDEGC
                 FLEFTR = II (EP+5)
                 OKROW = FLEFTR .LE. FFDIMR
                 IF (.NOT. OKROW) THEN 
C                   exit out of loop if row too long:
                    GO TO 1000
                 ENDIF 
C                split into two loops so that both vectorize on a CRAY
                 F1 = FEDEGR
                 DO 970 P = LURP, LURP + LUDEGR - 1 
                    COL = II (P)
                    IF (COL .GT. 0) THEN 
                       IF (WIC (COL) .EQ. -2) THEN 
                          F1 = F1 + 1
                          WPR (F1) = COL
                       ENDIF 
                    ENDIF 
970              CONTINUE 
                 OKROW = F1 .LE. FFDIMR
                 IF (.NOT. OKROW) THEN 
C                   exit out of loop if row too long:
                    GO TO 1000
                 ENDIF 
                 DO 980 I = FEDEGR+1, F1 
                    WIC (WPR (I)) = (I - 1) * FFDIMC
980              CONTINUE 
                 FEDEGR = F1
              ELSE 
C                this is an artificial element (a dense column)
                 LURP = (EP+8)
                 COL = II (LURP)
                 IF (WIC (COL) .EQ. -2) THEN 
                    WIC (COL) = FEDEGR * FFDIMC
                    FEDEGR = FEDEGR + 1
                    WPR (FEDEGR) = COL
                    OKROW = FEDEGR .LE. FFDIMR
                    IF (.NOT. OKROW) THEN 
C                      exit out of loop if row too long:
                       GO TO 1000
                    ENDIF 
                 ENDIF 
              ENDIF 
990        CONTINUE 
C          loop exit label:
1000       CONTINUE
        ENDIF 

        IF (OKROW) THEN 

C          -------------------------------------------------------------
C          assemble the original entries in the row
C          -------------------------------------------------------------

           RLEN = WC (PIVROW)
           IF (RLEN .GT. 0) THEN 
              F1 = FEDEGR
              RSIZ = II (PR)
              P2 = PR + RSIZ
C             split into two loops so that they both vectorize on a CRAY
              DO 1010 P = P2 - RLEN, P2 - 1 
                 COL = II (P)
                 IF (WIC (COL) .EQ. -2) THEN 
C                   this entry cannot be assembled, do not delete
                    F1 = F1 + 1
                    WPR (F1) = COL
                 ENDIF 
1010          CONTINUE 
              RLEN2 = F1 - FEDEGR
              IF (RLEN2 .LT. RLEN) THEN 
C                delete one or more entries in the row
                 DO 1020 I = FEDEGR+1, F1 
                    II (P2 - F1 + I - 1) = WPR (I)
1020             CONTINUE 
                 INEED = INEED - (RLEN - RLEN2)
                 WC (PIVROW) = RLEN2
              ENDIF 

C             ----------------------------------------------------------
C             update the candidate row degree - excluding current front
C             ----------------------------------------------------------

              RDEG = F1 - FFLEFR
              II (PR+1) = RDEG

C             ----------------------------------------------------------
C             pivot is found if candidate pivot row is not too long
C             ----------------------------------------------------------

              OKROW = F1 .LE. FFDIMR
              IF (OKROW) THEN 
                 DO 1030 I = FEDEGR+1, F1 
                    WIC (WPR (I)) = (I - 1) * FFDIMC
1030             CONTINUE 
                 FEDEGR = F1
              ENDIF 

           ELSE 

C             ----------------------------------------------------------
C             update the candidate row degree - excluding current front
C             ----------------------------------------------------------

              RDEG = FEDEGR - FFLEFR
              II (PR+1) = RDEG
           ENDIF 
        ENDIF 

C       ----------------------------------------------------------------
C       if pivot not found: clear WiR and WiC
C       ----------------------------------------------------------------

        PFOUND = OKROW .AND. OKCOL
        IF (.NOT. PFOUND) THEN 
           MOVELU = K .GT. 0
           DO 1040 I = FFLEFR+1, FEDEGR 
              WIC (WPR (I)) = -2
1040       CONTINUE 
           FEDEGR = FFLEFR
           DO 1050 I = FFLEFC+1, FEDEGC 
              WIR (WPC (I)) = -1
1050       CONTINUE 
           FEDEGC = FFLEFC
        ELSE 
           MOVELU = FEDEGC .GT. FFDIMC - K .OR. FEDEGR .GT. FFDIMR - K
        ENDIF 

C       ----------------------------------------------------------------
C       WpR (1..fflefr)                 unextended pivot row pattern
C       WpR (fflefr+1 .. fedegr)        extended pattern, if pfound
C       WpR (fedegr+1 .. n-npiv)        empty space
C       WpR (n-npiv+1 .. n)             pivot row order
C
C       WpC (1..fflefc)                 unextended pivot column pattern
C       WpC (fflefc+1 .. fedegc)        extended pattern, if pfound
C       WpC (fedegc+1 .. n-npiv)        empty space
C       WpC (n-npiv+1 .. n)             pivot column order
C       ----------------------------------------------------------------

C=======================================================================
C  Local pivot search complete ]
C=======================================================================

C=======================================================================
C  Update contribution block: rank-nb, or if LU arrowhead to be moved
C=======================================================================

        IF (K-K0 .GE. NB .OR. MOVELU) THEN 
           CALL DGEMM ('N', 'N', FFLEFC, FFLEFR, K-K0,
     $          -ONE, XX (FFXP + (FFDIMR - K) * FFDIMC), FFDIMC,
     $                XX (FFXP +  FFDIMC - K)          , FFDIMC,
     $           ONE, XX (FFXP)                        , FFDIMC)
           TMP = FFLEFC
           TMP = TMP * FFLEFR
           TMP = TMP * (K-K0)
           RINFO (6) = RINFO (6) + 2.0*(TMP)
           K0 = K
        ENDIF 

C=======================================================================
C  Move the LU arrowhead if no pivot found, or pivot needs room
C=======================================================================

        IF (MOVELU) THEN 

C          allocate permanent space for the LU arrowhead
           LUDEGR = FFLEFR
           LUDEGC = FFLEFC
           XS = K*LUDEGC + K*LUDEGR + K*K
           IS = 7 + LUDEGC + LUDEGR + NSONS
           IF (IS .GT. ITAIL-IHEAD .OR. XS .GT. XTAIL-XHEAD) THEN 
              IF (IS .GT. ITAIL-IHEAD) THEN 
C                garbage collection because we ran out of integer mem
                 INFO (14) = INFO (14) + 1
              ENDIF 
              IF (XS .GT. XTAIL-XHEAD) THEN 
C                garbage collection because we ran out of real mem
                 INFO (15) = INFO (15) + 1
              ENDIF 
              CALL UMD2FG (XX, XSIZE, XHEAD, XTAIL, XUSE,
     $                     II, ISIZE, IHEAD, ITAIL, IUSE,
     $                     CP, RP, DN, N, ICNTL, WIR, WIC, WR, WC,
     $                     FFXP, FFSIZE, WXP, FFDIMC, .FALSE.,
     $                     PFREE, XFREE, MHEAD, MTAIL, SLOTS)
C             at this point, iuse = ineed and xuse = xneed
           ENDIF 

           ITAIL = ITAIL - IS
           LUIP = ITAIL
           IUSE = IUSE + IS
           INEED = INEED + IS
           XTAIL = XTAIL - XS
           LUXP = XTAIL
           XUSE = XUSE + XS
           XNEED = XNEED + XS
           INFO (18) = MAX (INFO (18), IUSE)
           INFO (19) = MAX (INFO (19), INEED)
           INFO (20) = MAX (INFO (20), XUSE)
           INFO (21) = MAX (INFO (21), XNEED)
           IF (IHEAD .GT. ITAIL .OR. XHEAD .GT. XTAIL) THEN 
C             error return, if not enough integer and/or real memory:
              GO TO 9000
           ENDIF 

C          -------------------------------------------------------------
C          get memory usage for next call to UMD2RF
C          -------------------------------------------------------------

           XRUSE = XRUSE + XS
           XRMAX = MAX (XRMAX, XRUSE)

C          -------------------------------------------------------------
C          save the new LU arrowhead
C          -------------------------------------------------------------

C          save the scalar data of the LU arrowhead
           II (LUIP) = LUXP
           II (LUIP+1) = K
           II (LUIP+2) = LUDEGR
           II (LUIP+3) = LUDEGC
           II (LUIP+4) = NSONS
           II (LUIP+5) = 0
           II (LUIP+6) = 0
           E = FFROW
           IF (E .EQ. E1) THEN 
C             this is the first LU arrowhead from this global pivot
              LUIP1 = LUIP
           ENDIF 
           WR (E) = -(NDN+2)
           WC (E) = -(NDN+2)

C          save column pattern
           LUCP = (LUIP + 7)
           DO 1060 I = 0, LUDEGC-1 
              II (LUCP+I) = WPC (I+1)
1060       CONTINUE 

C          save row pattern
           LURP = LUCP + LUDEGC
           DO 1070 I = 0, LUDEGR-1 
              II (LURP+I) = WPR (I+1)
1070       CONTINUE 

C          add list of sons after the end of the frontal matrix pattern
C          this list of sons is for the refactorization (UMD2RF) only.
           LUSONP = LURP + LUDEGR
           IP = LUSONP
           E = SONLST
C          while (e > 0) do
1080       CONTINUE
           IF (E .GT. 0) THEN 
              EP = RP (E)
              IF (WC (E) .EQ. -(NDN+2)) THEN 
C                LUson
                 II (IP) = E
              ELSE IF (WC (E) .EQ. W0) THEN 
C                Uson
                 II (IP) = E + N
              ELSE IF (WR (E) .EQ. W0) THEN 
C                Lson
                 II (IP) = E + 2*N
              ENDIF 
              NEXT = WIR (E) + N + 2
              WIR (E) = -1
              E = NEXT
              IP = IP + 1
C          end while:
           GOTO 1080
           ENDIF 
           NSONS = 0
           SONLST = 0

C          move the L1,U1 matrix, compressing the dimension from
C          ffdimc to ldimc.  The LU arrowhead grows on top of stack.
           LDIMC = K + LUDEGC
           XP = FFXP + (FFDIMR-1)*FFDIMC + FFDIMC-1
           DO 1100 J = 0, K-1 
CFPP$ NODEPCHK L
              DO 1090 I = 0, K-1 
                 XX (LUXP + J*LDIMC + I) = XX (XP - J*FFDIMC - I)
1090          CONTINUE 
1100       CONTINUE 

C          move L2 matrix, compressing dimension from ffdimc to ludegc+k
           IF (LUDEGC .NE. 0) THEN 
              LXP = LUXP + K
              XP = FFXP + (FFDIMR-1)*FFDIMC
              DO 1120 J = 0, K-1 
CFPP$ NODEPCHK L
                 DO 1110 I = 0, LUDEGC-1 
                    XX (LXP + J*LDIMC + I) = XX (XP - J*FFDIMC + I)
1110             CONTINUE 
1120          CONTINUE 
           ENDIF 

C          move the U2 block.
           IF (LUDEGR .NE. 0) THEN 
              UXP = LUXP + K * LDIMC
              XP = FFXP + FFDIMC-1
              DO 1140 J = 0, LUDEGR-1 
CFPP$ NODEPCHK L
                 DO 1130 I = 0, K-1 
                    XX (UXP + J*K + I) = XX (XP + J*FFDIMC - I)
1130             CONTINUE 
1140          CONTINUE 
           ENDIF 

C          one more LU arrowhead has been created
           NLU = NLU + 1
           NZU = (K*(K-1)/2) + K*LUDEGC
           NZL = (K*(K-1)/2) + K*LUDEGR
           INFO (10) = INFO (10) + NZL
           INFO (11) = INFO (11) + NZU

C          no more rows of U or columns of L in current frontal array
           K = 0
           K0 = 0

           IF (PFOUND) THEN 

C             ----------------------------------------------------------
C             Place the old frontal matrix as the only item in the son
C             list, since the next "implied" frontal matrix will have
C             this as its son.
C             ----------------------------------------------------------

              NSONS = 1
              E = FFROW
              WIR (E) = - N - 2
              SONLST = E

C             ----------------------------------------------------------
C             The contribution block of the old frontal matrix is still
C             stored in the current frontal matrix, and continues (in a
C             unifrontal sense) as a "new" frontal matrix (same array
C             but with a new name, and the LU arrowhead is removed and
C             placed in the LU factors).  Old name is "ffrow", new name
C             is "pivrow".
C             ----------------------------------------------------------

              RP (E) = LUIP
              FFROW = PIVROW
           ENDIF 
        ENDIF 

C=======================================================================
C  Stop the factorization of this frontal matrix if no pivot found
C=======================================================================

C       (this is the only way out of loop 1395)
        IF (.NOT. PFOUND) THEN 
C          exit out of loop 1395 if pivot not found:
           GO TO 1400
        ENDIF 

C=======================================================================
C  Update the pivot column, and move into position as (k+1)-st col of L
C=======================================================================

        XSP = (COLPOS - 1) * FFDIMC
        XDP = (FFDIMR - K - 1) * FFDIMC
        FSP = FFXP + XSP
        FDP = FFXP + XDP

        IF (K-K0 .GT. 0 .AND. FFLEFC .NE. 0) THEN 
           CALL DGEMV ('N', FFLEFC, K-K0,
     $          -ONE, XX (FDP + FFDIMC    ), FFDIMC,
     $                XX (FSP + FFDIMC - K), 1,
     $           ONE, XX (FSP             ), 1)
           TMP = FFLEFC
           TMP = TMP * (K-K0)
           RINFO (5) = RINFO (5) + 2.0*(TMP)
        ENDIF 

        IF (FFLEFR .LT. FFDIMR - K) THEN 

           XLP = (FFLEFR - 1) * FFDIMC
           IF (FFLEFR .EQ. COLPOS) THEN 

C             ----------------------------------------------------------
C             move C(:,colpos) => C(:,ffdimr-k)
C             ----------------------------------------------------------

C                column of the contribution block:
CFPP$ NODEPCHK L
                 DO 1160 I = 0, FFLEFC - 1 
                    XX (FDP+I) = XX (FSP+I)
1160             CONTINUE 
C                column of the U2 block
CFPP$ NODEPCHK L
                 DO 1170 I = FFDIMC - K, FFDIMC - 1 
                    XX (FDP+I) = XX (FSP+I)
1170             CONTINUE 

           ELSE 

C             ----------------------------------------------------------
C             move C(:,colpos) => C(:,ffdimr-k)
C             move C(:,fflefr) => C(:,colpos)
C             ----------------------------------------------------------

              FLP = FFXP + XLP

C                columns of the contribution block:
CFPP$ NODEPCHK L
                 DO 1190 I = 0, FFLEFC - 1 
                    XX (FDP+I) = XX (FSP+I)
                    XX (FSP+I) = XX (FLP+I)
1190             CONTINUE 
C                columns of the U2 block:
CFPP$ NODEPCHK L
                 DO 1200 I = FFDIMC - K, FFDIMC - 1 
                    XX (FDP+I) = XX (FSP+I)
                    XX (FSP+I) = XX (FLP+I)
1200             CONTINUE 

              SWPCOL = WPR (FFLEFR)
              WPR (COLPOS) = SWPCOL
              WIC (SWPCOL) = XSP
           ENDIF 

           IF (FEDEGR .NE. FFLEFR) THEN 
C             move column fedegr to column fflefr (pattern only)
              SWPCOL = WPR (FEDEGR)
              WPR (FFLEFR) = SWPCOL
              WIC (SWPCOL) = XLP
           ENDIF 

        ELSE IF (COLPOS .NE. FFDIMR - K) THEN 

C          -------------------------------------------------------------
C          swap C(:,colpos) <=> C (:,ffdimr-k)
C          -------------------------------------------------------------

C             swap only what needs to be swapped
C             columns of the contribution block:
CFPP$ NODEPCHK L
CFPP$ NOLSTVAL L
              DO 1220 I = 0, FFLEFC - 1 
                 TEMP = XX (FDP+I)
                 XX (FDP+I) = XX (FSP+I)
                 XX (FSP+I) = TEMP
1220          CONTINUE 
C             columns of the U2 block:
CFPP$ NODEPCHK L
CFPP$ NOLSTVAL L
              DO 1230 I = FFDIMC - K, FFDIMC - 1 
                 TEMP = XX (FDP+I)
                 XX (FDP+I) = XX (FSP+I)
                 XX (FSP+I) = TEMP
1230          CONTINUE 

           SWPCOL = WPR (FFDIMR - K)
           WPR (COLPOS) = SWPCOL
           WIC (SWPCOL) = XSP
        ENDIF 

        WIC (PIVCOL) = XDP
        FEDEGR = FEDEGR - 1
        SCAN2 = FFLEFR
        FFLEFR = FFLEFR - 1

C=======================================================================
C  Move pivot row into position as (k+1)-st row of U, and update
C=======================================================================

        XSP = ROWPOS - 1
        XDP = FFDIMC - K - 1
        FSP = FFXP + XSP
        FDP = FFXP + XDP

        IF (FFLEFC .LT. FFDIMC - K) THEN 

           XLP = FFLEFC - 1
           IF (FFLEFC .EQ. ROWPOS) THEN 

C             ----------------------------------------------------------
C             move C(rowpos,:) => C(ffdimc-k,:)
C             ----------------------------------------------------------

C                row of the contribution block:
CFPP$ NODEPCHK L
                 DO 1250 J = 0, (FFLEFR - 1) * FFDIMC, FFDIMC 
                    XX (FDP+J) = XX (FSP+J)
1250             CONTINUE 
C                row of the L2 block:
CFPP$ NODEPCHK L
                 DO 1260 J = (FFDIMR - K - 1) * FFDIMC,
     $                       (FFDIMR - 1) * FFDIMC, FFDIMC 
                    XX (FDP+J) = XX (FSP+J)
1260             CONTINUE 

           ELSE 

C             ----------------------------------------------------------
C             move C(rowpos,:) => C(ffdimc-k,:)
C             move C(fflefc,:) => C(rowpos,:)
C             ----------------------------------------------------------

              FLP = FFXP + XLP

C                rows of the contribution block:
CFPP$ NODEPCHK L
                 DO 1280 J = 0, (FFLEFR - 1) * FFDIMC, FFDIMC 
                    XX (FDP+J) = XX (FSP+J)
                    XX (FSP+J) = XX (FLP+J)
1280             CONTINUE 
C                rows of the L2 block:
CFPP$ NODEPCHK L
                 DO 1290 J = (FFDIMR - K - 1) * FFDIMC,
     $                       (FFDIMR - 1) * FFDIMC, FFDIMC 
                    XX (FDP+J) = XX (FSP+J)
                    XX (FSP+J) = XX (FLP+J)
1290             CONTINUE 

              SWPROW = WPC (FFLEFC)
              WPC (ROWPOS) = SWPROW
              WIR (SWPROW) = XSP
           ENDIF 

           IF (FEDEGC .NE. FFLEFC) THEN 
C             move row fedegc to row fflefc (pattern only)
              SWPROW = WPC (FEDEGC)
              WPC (FFLEFC) = SWPROW
              WIR (SWPROW) = XLP
           ENDIF 

        ELSE IF (ROWPOS .NE. FFDIMC - K) THEN 

C          -------------------------------------------------------------
C          swap C(rowpos,:) <=> C (ffdimc-k,:)
C          -------------------------------------------------------------

C             swap only what needs to be swapped
C             rows of the contribution block:
CFPP$ NODEPCHK L
CFPP$ NOLSTVAL L
              DO 1310 J = 0, (FFLEFR - 1) * FFDIMC, FFDIMC 
                 TEMP = XX (FDP+J)
                 XX (FDP+J) = XX (FSP+J)
                 XX (FSP+J) = TEMP
1310          CONTINUE 
C             rows of the L2 block:
CFPP$ NODEPCHK L
CFPP$ NOLSTVAL L
              DO 1320 J = (FFDIMR - K - 1) * FFDIMC,
     $                    (FFDIMR - 1) * FFDIMC, FFDIMC 
                 TEMP = XX (FDP+J)
                 XX (FDP+J) = XX (FSP+J)
                 XX (FSP+J) = TEMP
1320          CONTINUE 

           SWPROW = WPC (FFDIMC - K)
           WPC (ROWPOS) = SWPROW
           WIR (SWPROW) = XSP
        ENDIF 

        WIR (PIVROW) = XDP
        FEDEGC = FEDEGC - 1
        SCAN1 = FFLEFC
        FFLEFC = FFLEFC - 1

        IF (K-K0 .GT. 0 .AND. FFLEFR .GT. 0) THEN 
           CALL DGEMV ('T', K-K0, FFLEFR,
     $       -ONE, XX (FDP + 1)                    , FFDIMC,
     $             XX (FDP + (FFDIMR - K) * FFDIMC), FFDIMC,
     $        ONE, XX (FDP)                        , FFDIMC)
           TMP = K-K0
           TMP = TMP * FFLEFR
           RINFO (5) = RINFO (5) + 2.0*(TMP)
        ENDIF 

C=======================================================================
C  Prepare for degree update and next local pivot search
C=======================================================================

C       ----------------------------------------------------------------
C       if only column pattern has been extended:
C               scan1:  new rows only
C               scan2:  no columns scanned
C               scan3:  all columns
C               scan4:  new rows only
C
C       if only row pattern has been extended:
C               scan1:  no rows scanned
C               scan2:  new columns only
C               scan3:  new columns only
C               scan4:  all rows
C
C       if both row and column pattern have been extended:
C               scan1:  new rows only
C               scan2:  new columns only
C               scan3:  all columns
C               scan4:  all rows
C
C       if no patterns have been extended:
C               scan1-4: none
C       ----------------------------------------------------------------

        IF (FEDEGC .EQ. FFLEFC) THEN 
C          column pattern has not been extended
           SCAN3 = FFLEFR + 1
        ELSE 
C          column pattern has been extended.
           SCAN3 = 0
        ENDIF 

        IF (FEDEGR .EQ. FFLEFR) THEN 
C          row pattern has not been extended
           SCAN4 = FFLEFC + 1
        ELSE 
C          row pattern has been extended
           SCAN4 = 0
        ENDIF 

C=======================================================================
C  Finished with step k (except for assembly and scaling of pivot col)
C=======================================================================

        K = K + 1

C       ----------------------------------------------------------------
C       exit loop if frontal matrix has been extended
C       ----------------------------------------------------------------

        IF (FEDEGR .NE. FFLEFR .OR. FEDEGC .NE. FFLEFC) THEN 
           GO TO 1325
        ENDIF 

1324    CONTINUE 
C       exit label for loop 1324:
1325    CONTINUE

C=======================================================================
C  Finished factorizing while frontal matrix is not extended ]
C=======================================================================

C=======================================================================
C  Extend the frontal matrix [
C=======================================================================

C       ----------------------------------------------------------------
C       Zero the newly extended frontal matrix
C       ----------------------------------------------------------------

C       fill-in due to amalgamation caused by this step is
C       k*(fedegr-fflefr+fedegc-fflefc)

        DO 1350 J = FFLEFR, FEDEGR - 1 
C          zero the new columns in the contribution block:
           DO 1330 I = 0, FEDEGC - 1 
              XX (FFXP + J*FFDIMC + I) = 0
1330       CONTINUE 
C          zero the new columns in U block:
           DO 1340 I = FFDIMC - K, FFDIMC - 1 
              XX (FFXP + J*FFDIMC + I) = 0
1340       CONTINUE 
1350    CONTINUE 

CFPP$ NODEPCHK L
        DO 1380 I = FFLEFC, FEDEGC - 1 
C          zero the new rows in the contribution block:
CFPP$ NODEPCHK L
           DO 1360 J = 0, FFLEFR - 1 
              XX (FFXP + J*FFDIMC + I) = 0
1360       CONTINUE 
C          zero the new rows in L block:
CFPP$ NODEPCHK L
           DO 1370 J = FFDIMR - K, FFDIMR - 1 
              XX (FFXP + J*FFDIMC + I) = 0
1370       CONTINUE 
1380    CONTINUE 

C       ----------------------------------------------------------------
C       remove the new columns from the degree lists
C       ----------------------------------------------------------------

        DO 1390 J = FFLEFR+1, FEDEGR 
           PC = CP (WPR (J))
           CDEG = II (PC+1)
           IF (CDEG .GT. 0) THEN 
              CNXT = II (PC+7)
              CPRV = II (PC+8)
              IF (CNXT .NE. 0) THEN 
                 II (CP (CNXT)+8) = CPRV
              ENDIF 
              IF (CPRV .NE. 0) THEN 
                 II (CP (CPRV)+7) = CNXT
              ELSE 
                 HEAD (CDEG) = CNXT
              ENDIF 
           ENDIF 
1390    CONTINUE 

C       ----------------------------------------------------------------
C       finalize extended row and column pattern of the frontal matrix
C       ----------------------------------------------------------------

        FFLEFC = FEDEGC
        FFLEFR = FEDEGR
        FMAXR = MAX (FMAXR, FFLEFR + K)
        FMAXC = MAX (FMAXC, FFLEFC + K)

C=======================================================================
C  Done extending the current frontal matrix ]
C=======================================================================

1395    CONTINUE 
C       exit label for loop 1395:
1400    CONTINUE

C=======================================================================
C  Done assembling and factorizing the current frontal matrix ]
C=======================================================================

C=======================================================================
C  Wrap-up:  complete the current frontal matrix [
C=======================================================================

C       ----------------------------------------------------------------
C       store the maximum front size in the first LU arrowhead
C       ----------------------------------------------------------------

        II (LUIP1+5) = FMAXR
        II (LUIP1+6) = FMAXC

C       one more frontal matrix is finished
        INFO (13) = INFO (13) + 1

C       ----------------------------------------------------------------
C       add the current frontal matrix to the degrees of each column,
C       and place the modified columns back in the degree lists
C       ----------------------------------------------------------------

C       do so in reverse order to try to improve pivot tie-breaking
        DO 1410 J = FFLEFR, 1, -1 
           COL = WPR (J)
           PC = CP (COL)
C          add the current frontal matrix to the degree
           CDEG = II (PC+1)
           CDEG = MIN (KLEFT, CDEG + FFLEFC)
           IF (CDEG .GT. 0) THEN 
              II (PC+1) = CDEG
              CNXT = HEAD (CDEG)
              II (PC+7) = CNXT
              II (PC+8) = 0
              IF (CNXT .NE. 0) THEN 
                 II (CP (CNXT)+8) = COL
              ENDIF 
              HEAD (CDEG) = COL
              MINDEG = MIN (MINDEG, CDEG)
           ENDIF 
1410    CONTINUE 

C       ----------------------------------------------------------------
C       add the current frontal matrix to the degrees of each row
C       ----------------------------------------------------------------

CFPP$ NODEPCHK L
        DO 1420 I = 1, FFLEFC 
           ROW = WPC (I)
           PR = RP (ROW)
           RDEG = II (PR+1)
           RDEG = MIN (KLEFT, RDEG + FFLEFR)
           II (PR+1) = RDEG
1420    CONTINUE 

C       ----------------------------------------------------------------
C       Reset w0 so that Wr (1..n) < w0 and Wc (1..n) < w0.
C       Also ensure that w0 + n would not cause integer overflow
C       ----------------------------------------------------------------

        W0 = W0 + FMAX + 1
        IF (W0 .GE. W0BIG) THEN 
           W0 = NDN+2
           DO 1430 E = 1, N+DN 
              IF (WR (E) .GT. NDN) THEN 
C                this is a frontal matrix
                 WR (E) = W0-1
                 WC (E) = W0-1
              ENDIF 
1430       CONTINUE 
        ENDIF 

C       ----------------------------------------------------------------
C       deallocate work vector
C       ----------------------------------------------------------------

        XUSE = XUSE - FFDIMC
        XNEED = XNEED - FFDIMC
        XHEAD = XHEAD - FFDIMC

C       ----------------------------------------------------------------
C       get the name of this new frontal matrix, and size of
C       contribution block
C       ----------------------------------------------------------------

        E = FFROW
        XS = FFLEFR * FFLEFC
        FMAX = MAX (FMAX, FFLEFR, FFLEFC)

C       ----------------------------------------------------------------
C       get memory usage for next call to UMD2RF
C       ----------------------------------------------------------------

        XRUSE = XRUSE - FFSIZE + XS

C       ----------------------------------------------------------------
C       if contribution block empty, deallocate and continue next step
C       ----------------------------------------------------------------

        IF (FFLEFR .LE. 0 .OR. FFLEFC .LE. 0) THEN 
           RP (E) = LUIP
           XUSE = XUSE - FFSIZE
           XNEED = XNEED - FFSIZE
           XHEAD = FFXP
           DO 1440 I = 1, FFLEFR 
              WIC (WPR (I)) = -2
1440       CONTINUE 
           DO 1450 I = 1, FFLEFC 
              WIR (WPC (I)) = -1
1450       CONTINUE 
C          next iteration of main factorization loop 1540:
           GOTO 1540
        ENDIF 

C       ----------------------------------------------------------------
C       prepare the contribution block for later assembly
C       ----------------------------------------------------------------

        IF (FSCAL .GT. ITAIL-IHEAD) THEN 
           INFO (14) = INFO (14) + 1
           CALL UMD2FG (XX, XSIZE, XHEAD, XTAIL, XUSE,
     $                  II, ISIZE, IHEAD, ITAIL, IUSE,
     $                  CP, RP, DN, N, ICNTL, WIR, WIC, WR, WC,
     $                  FFXP, FFSIZE, 0, 0, .FALSE.,
     $                  PFREE, XFREE, MHEAD, MTAIL, SLOTS)
C          at this point, iuse = ineed and xuse = xneed
        ENDIF 

        EP = IHEAD
        IHEAD = IHEAD + FSCAL
        IUSE = IUSE + FSCAL
        INEED = INEED + FSCAL
        INFO (18) = MAX (INFO (18), IUSE)
        INFO (19) = MAX (INFO (19), INEED)
        IF (IHEAD .GT. ITAIL) THEN 
C          error return, if not enough integer memory:
C          (highly unlikely to run out of memory at this point)
           GO TO 9000
        ENDIF 

        RP (E) = EP
        II (EP) = LUIP
        II (EP+5) = FFLEFR
        II (EP+6) = FFLEFC
        WR (E) = W0-1
        WC (E) = W0-1

C       count the numerical assembly
        RINFO (2) = RINFO (2) + (XS)

        IF (XS .LE. XFREE) THEN 

C          -------------------------------------------------------------
C          compress and store the contribution block in a freed block
C          -------------------------------------------------------------

C          place the new block in the list in front of the free block
           XDP = II (PFREE+2)
           II (PFREE+2) = II (PFREE+2) + XS
           XFREE = XFREE - XS
           MPREV = II (PFREE+4)
           IF (XFREE .EQ. 0) THEN 
C             delete the free block if its size is zero
              MNEXT = II (PFREE+3)
              PFREE = 0
              XFREE = -1
           ELSE 
              MNEXT = PFREE
           ENDIF 
           IF (MNEXT .NE. 0) THEN 
              II (MNEXT+4) = EP
           ELSE 
              MTAIL = EP
           ENDIF 
           IF (MPREV .NE. 0) THEN 
              II (MPREV+3) = EP
           ELSE 
              MHEAD = EP
           ENDIF 
           DO 1470 J = 0, FFLEFR - 1 
CFPP$ NODEPCHK L
              DO 1460 I = 0, FFLEFC - 1 
                 XX (XDP + J*FFLEFC + I) = XX (FFXP + J*FFDIMC + I)
1460          CONTINUE 
1470       CONTINUE 
           XHEAD = FFXP
           XUSE = XUSE - FFSIZE
           XNEED = XNEED - FFSIZE + XS
           FFDIMC = FFLEFC
           II (EP+1) = FFDIMC
           II (EP+2) = XDP
           II (EP+3) = MNEXT
           II (EP+4) = MPREV

        ELSE 

C          -------------------------------------------------------------
C          deallocate part of the unused portion of the frontal matrix
C          -------------------------------------------------------------

C          leave the contribution block C (1..fflefc, 1..fflefr) at the
C          head of XX, with column dimension of ffdimc and in space
C          of size (fflefr-1)*ffdimc for the first fflefr columns, and
C          fflefc for the last column.
           XNEED = XNEED - FFSIZE + XS
           XS = FFSIZE - (FFLEFC + (FFLEFR-1)*FFDIMC)
           XHEAD = XHEAD - XS
           XUSE = XUSE - XS
           II (EP+1) = FFDIMC
           II (EP+2) = FFXP
           II (EP+3) = 0
           II (EP+4) = MTAIL
           IF (MTAIL .EQ. 0) THEN 
              MHEAD = EP
           ELSE 
              II (MTAIL+3) = EP
           ENDIF 
           MTAIL = EP
        ENDIF 

C       ----------------------------------------------------------------
C       add tuples to the amount of integer space needed - and add
C       limit+cscal to maximum need to account for worst-case possible
C       reallocation of rows/columns.  Required integer memory usage
C       is guaranteed not to exceed iworst during the placement of (e,f)
C       tuples in the two loops below.
C       ----------------------------------------------------------------

        INEED = INEED + 2*(FFLEFR+FFLEFC)
        IWORST = INEED + LIMIT + CSCAL
        INFO (19) = MAX (INFO (19), IWORST)
        INFO (18) = MAX (INFO (18), IWORST)

C       ----------------------------------------------------------------
C       place (e,f) in the element list of each column
C       ----------------------------------------------------------------

        DO 1500 I = 1, FFLEFR 
           COL = WPR (I)
           PC = CP (COL)
           CELN = II (PC+5)
           CSIZ = II (PC)
           CLEN = II (PC+6)
C          clear the column offset
           WIC (COL) = -2

C          -------------------------------------------------------------
C          make sure an empty slot exists - if not, create one
C          -------------------------------------------------------------

           IF (2*(CELN+1) + CLEN + CSCAL .GT. CSIZ) THEN 

C             ----------------------------------------------------------
C             no room exists - reallocate elsewhere
C             ----------------------------------------------------------

C             at least this much space is needed:
              IS = 2 * (CELN + 1) + CLEN
C             add some slots for growth: at least 8 tuples,
C             or double the size - whichever is larger (but with a total
C             size not larger than limit+cscal)
              IS = MIN (IS + MAX (16, IS), LIMIT)
              CSIZ2 = IS + CSCAL

C             ----------------------------------------------------------
C             make sure enough room exists: garbage collection if needed
C             ----------------------------------------------------------

              IF (CSIZ2 .GT. ITAIL-IHEAD) THEN 
C                garbage collection:
                 INFO (14) = INFO (14) + 1
                 CALL UMD2FG (XX, XSIZE, XHEAD, XTAIL, XUSE,
     $                        II, ISIZE, IHEAD, ITAIL, IUSE,
     $                        CP, RP, DN, N, ICNTL, WIR, WIC, WR, WC,
     $                        0, 0, 0, 0, .TRUE.,
     $                        PFREE, XFREE, MHEAD, MTAIL, SLOTS)
C                at this point, iuse+csiz2 <= iworst and xuse = xneed
                 PC = CP (COL)
                 CSIZ = II (PC)
              ENDIF 

C             ----------------------------------------------------------
C             get space for the new copy
C             ----------------------------------------------------------

              PC2 = IHEAD
              IHEAD = IHEAD + CSIZ2
              IUSE = IUSE + CSIZ2
              INFO (18) = MAX (INFO (18), IUSE)
              IF (IHEAD .GT. ITAIL) THEN 
C                error return, if not enough integer memory:
                 GO TO 9000
              ENDIF 

C             ----------------------------------------------------------
C             make the copy, leaving hole in middle for element list
C             ----------------------------------------------------------

C             copy the cscal scalars, and the element list
CFPP$ NODEPCHK L
              DO 1480 J = 0, CSCAL + 2*CELN - 1 
                 II (PC2+J) = II (PC+J)
1480          CONTINUE 

C             copy column indices of original entries (XX is unchanged)
CFPP$ NODEPCHK L
              DO 1490 J = 0, CLEN - 1 
                 II (PC2+CSIZ2-CLEN+J) = II (PC+CSIZ-CLEN+J)
1490          CONTINUE 

              IF (CLEN .GT. 0) THEN 
C                place the new block in the memory-list
                 MNEXT = II (PC2+3) 
                 MPREV = II (PC2+4)
                 IF (MNEXT .NE. 0) THEN 
                    II (MNEXT+4) = PC2
                 ELSE 
                    MTAIL = PC2
                 ENDIF 
                 IF (MPREV .NE. 0) THEN 
                    II (MPREV+3) = PC2
                 ELSE 
                    MHEAD = PC2
                 ENDIF 
              ENDIF 

              CP (COL) = PC2
              II (PC2) = CSIZ2

C             ----------------------------------------------------------
C             deallocate the old copy of the column in II (not in XX)
C             ----------------------------------------------------------

              II (PC+1) = -1
              II (PC+6) = 0
              PC = PC2
           ENDIF 

C          -------------------------------------------------------------
C          place the new (e,f) tuple in the element list of the column
C          -------------------------------------------------------------

           CEP = (PC+9)
           II (CEP + 2*CELN  ) = E
           II (CEP + 2*CELN+1) = I - 1
           II (PC+5) = CELN + 1
1500    CONTINUE 

C       ----------------------------------------------------------------
C       place (e,f) in the element list of each row
C       ----------------------------------------------------------------

        DO 1530 I = 1, FFLEFC 
           ROW = WPC (I)
           PR = RP (ROW)
           RSIZ = II (PR)
           RELN = WR (ROW)
           RLEN = WC (ROW)
C          clear the row offset
           WIR (ROW) = -1

C          -------------------------------------------------------------
C          make sure an empty slot exists - if not, create one
C          -------------------------------------------------------------

           IF (2*(RELN+1) + RLEN + RSCAL .GT. RSIZ) THEN 

C             ----------------------------------------------------------
C             no room exists - reallocate elsewhere
C             ----------------------------------------------------------

C             at least this much space is needed:
              IS = 2 * (RELN + 1) + RLEN
C             add some extra slots for growth - for at least 8
C             tuples, or double the size (but with a total size not
C             larger than limit+rscal)
              IS = MIN (IS + MAX (16, IS), LIMIT)
              RSIZ2 = IS + RSCAL

C             ----------------------------------------------------------
C             make sure enough room exists: garbage collection if needed
C             ----------------------------------------------------------

              IF (RSIZ2 .GT. ITAIL-IHEAD) THEN 
C                garbage collection:
                 INFO (14) = INFO (14) + 1
                 CALL UMD2FG (XX, XSIZE, XHEAD, XTAIL, XUSE,
     $                        II, ISIZE, IHEAD, ITAIL, IUSE,
     $                        CP, RP, DN, N, ICNTL, WIR, WIC, WR, WC,
     $                        0, 0, 0, 0, .TRUE.,
     $                        PFREE, XFREE, MHEAD, MTAIL, SLOTS)
C                at this point, iuse+rsiz2 <= iworst and xuse = xneed
                 PR = RP (ROW)
                 RSIZ = II (PR)
              ENDIF 

C             ----------------------------------------------------------
C             get space for the new copy
C             ----------------------------------------------------------

              PR2 = IHEAD
              IHEAD = IHEAD + RSIZ2
              IUSE = IUSE + RSIZ2
              INFO (18) = MAX (INFO (18), IUSE)
              IF (IHEAD .GT. ITAIL) THEN 
C                error return, if not enough integer memory:
                 GO TO 9000
              ENDIF 

C             ----------------------------------------------------------
C             make the copy, leaving hole in middle for element list
C             ----------------------------------------------------------

C             copy the rscal scalars, and the element list
CFPP$ NODEPCHK L
              DO 1510 J = 0, RSCAL + 2*RELN - 1 
                 II (PR2+J) = II (PR+J)
1510          CONTINUE 

C             copy the original entries
CFPP$ NODEPCHK L
              DO 1520 J = 0, RLEN - 1 
                 II (PR2+RSIZ2-RLEN+J) = II (PR+RSIZ-RLEN+J)
1520          CONTINUE 

              RP (ROW) = PR2
              II (PR2) = RSIZ2

C             ----------------------------------------------------------
C             deallocate the old copy of the row
C             ----------------------------------------------------------

              II (PR+1) = -1
              PR = PR2
           ENDIF 

C          -------------------------------------------------------------
C          place the new (e,f) tuple in the element list of the row
C          -------------------------------------------------------------

           REP = (PR+2)
           II (REP + 2*RELN  ) = E
           II (REP + 2*RELN+1) = I - 1
           WR (ROW) = RELN + 1
1530    CONTINUE 

C=======================================================================
C  Wrap-up of factorized frontal matrix is complete ]
C=======================================================================

1540    CONTINUE 
C       exit label for loop 1540:
2000    CONTINUE

C=======================================================================
C=======================================================================
C  END OF MAIN FACTORIZATION LOOP ]
C=======================================================================
C=======================================================================

C=======================================================================
C  Wrap-up:  store LU factors in their final form [
C=======================================================================

C       ----------------------------------------------------------------
C       deallocate all remaining columns, rows, and frontal matrices
C       ----------------------------------------------------------------

        IUSE = IUSE - (IHEAD - 1)
        XUSE = XUSE - (XHEAD - 1)
        INEED = IUSE
        XNEED = XUSE
        IHEAD = 1
        XHEAD = 1

        IF (NLU .EQ. 0) THEN 
C          LU factors are completely empty (A = 0).
C          Add one integer and one real, to simplify rest of code.
C          Otherwise, some arrays in UMD2RF or UMD2SO would have
C          zero size, which can cause an address fault.
           ITAIL = ISIZE
           XTAIL = XSIZE
           IUSE = IUSE + 1
           XUSE = XUSE + 1
           INEED = IUSE
           XNEED = XUSE
           IP = ITAIL
           XP = XTAIL
        ENDIF 

C       ----------------------------------------------------------------
C       compute permutation and inverse permutation vectors.
C       use WiR/C for the row/col permutation, and WpR/C for the
C       inverse row/col permutation.
C       ----------------------------------------------------------------

        DO 2010 K = 1, N 
C          the kth pivot row and column:
           ROW = WPR (N-K+1)
           COL = WPC (N-K+1)
           WIR (K) = ROW
           WIC (K) = COL
2010    CONTINUE 
C       replace WpR/C with the inversion permutations:
        DO 2020 K = 1, N 
           ROW = WIR (K)
           COL = WIC (K)
           WPR (ROW) = K
           WPC (COL) = K
2020    CONTINUE 

        IF (PGIVEN) THEN 
C          the input matrix had been permuted from the original ordering
C          according to Rperm and Cperm.  Combine the initial
C          permutations (now in Rperm and Cperm) and the pivoting
C          permutations, and place them back into Rperm and Cperm.
           DO 2030 ROW = 1, N 
              WM (WPR (ROW)) = RPERM (ROW)
2030       CONTINUE 
           DO 2040 ROW = 1, N 
              RPERM (ROW) = WM (ROW)
2040       CONTINUE 
           DO 2050 COL = 1, N 
              WM (WPC (COL)) = CPERM (COL)
2050       CONTINUE 
           DO 2060 COL = 1, N 
              CPERM (COL) = WM (COL)
2060       CONTINUE 
C       else 
C          the input matrix was not permuted on input.  Rperm and Cperm
C          in UMD2F1 have been passed to this routine as WiR and WiC,
C          which now contain the row and column permutations.  Rperm and
C          Cperm in this routine (UMD2F2) are not defined.
        ENDIF 

C       ----------------------------------------------------------------
C       allocate nlu+3 integers for xtail, nlu, npiv and LUp (1..nlu)
C       ----------------------------------------------------------------

        IS = NLU + 5
        LUIP1 = ITAIL
        ITAIL = ITAIL - IS
        IUSE = IUSE + IS
        INEED = IUSE
        INFO (18) = MAX (INFO (18), IUSE)
        INFO (19) = MAX (INFO (19), INEED)
        IF (IHEAD .LE. ITAIL) THEN 

C          -------------------------------------------------------------
C          sufficient memory exist to finish the factorization
C          -------------------------------------------------------------

           II (ITAIL+1) = NLU
           II (ITAIL+2) = NPIV
           LUPP = ITAIL+5
           IF (NLU .EQ. 0) THEN 
C             zero the dummy entries, if LU factors are empty
              II (IP) = 0
              XX (XP) = 0
           ENDIF 

C          -------------------------------------------------------------
C          convert the LU factors into the new pivot order
C          -------------------------------------------------------------

           S = 0
           MAXDR = 1
           MAXDC = 1
           DO 2100 K = 1, N 
              E = WIR (K)
              LUIP = RP (E)
              IF (LUIP .GT. 0) THEN 
C                this is an LU arrowhead - save a pointer in LUp:
                 S = S + 1
C                update pointers to LU arrowhead relative to start of LU
                 II (LUPP+S-1) = LUIP - LUIP1 + 1
                 LUXP = II (LUIP)
                 II (LUIP) = LUXP - XTAIL + 1
C                convert the row and column indices to their final order
C                pattern of a column of L:
                 P = (LUIP + 7)
                 LUDEGC = II (LUIP+3)
                 MAXDC = MAX (MAXDC, LUDEGC)
                 DO 2070 J = 1, LUDEGC 
                    II (P) = WPR (ABS (II (P)))
                    P = P + 1
2070             CONTINUE 
C                pattern of a row of U:
                 LUDEGR = II (LUIP+2)
                 MAXDR = MAX (MAXDR, LUDEGR)
                 DO 2080 J = 1, LUDEGR 
                    II (P) = WPC (ABS (II (P)))
                    P = P + 1
2080             CONTINUE 
C                convert the LUsons, Usons, and Lsons:
                 NSONS = II (LUIP+4)
                 DO 2090 J = 1, NSONS 
                    ESON = II (P)
                    IF (ESON .LE. N) THEN 
C                      an LUson
                       II (P) = WM (ESON)
                    ELSE IF (ESON .LE. 2*N) THEN 
C                      a Uson
                       II (P) = WM (ESON-N) + N
                    ELSE 
C                      an Lson
                       II (P) = WM (ESON-2*N) + 2*N
                    ENDIF 
                    P = P + 1
2090             CONTINUE 
C                renumber this LU arrowhead
                 WM (E) = S
              ENDIF 
2100       CONTINUE 

           CMAX = MAX (CMAX, MAXDC)
           RMAX = MAX (RMAX, MAXDR)
           TOTNLU = TOTNLU + NLU

           II (ITAIL+3) = MAXDC
           II (ITAIL+4) = MAXDR

C          -------------------------------------------------------------
C          get memory usage for next call to UMD2RF
C          -------------------------------------------------------------

           XRUSE = XRUSE - NZ
           RETURN
        ENDIF 

C=======================================================================
C  LU factors are now stored in their final form ]
C=======================================================================

C=======================================================================
C  Error conditions
C=======================================================================

C       error return label:
9000    CONTINUE
        IF (IHEAD .GT. ITAIL .OR. ISIZE .LT. MINMEM) THEN 
C          error return if out of integer memory
           CALL UMD2ER (1, ICNTL, INFO, -3, INFO (19))
        ENDIF 
        IF (XHEAD .GT. XTAIL) THEN 
C          error return if out of real memory
           CALL UMD2ER (1, ICNTL, INFO, -4, INFO (21))
        ENDIF 
        RETURN
        END 
