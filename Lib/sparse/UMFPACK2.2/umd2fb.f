        SUBROUTINE UMD2FB (XX, XSIZE, II, ISIZE, N, NZ, NZDIA, NZOFF,
     $          NBLKS, CP, CPERM, RPERM, PR, PC,
     $          W, ZPERM, BP, OFFP,
     $          PRESRV, ICNTL)
        INTEGER N, NZ, ISIZE, II (ISIZE), NZDIA, NZOFF, NBLKS, CP (N+1),
     $          CPERM (N), RPERM (N), PR (N), PC (N), W (N), ZPERM (N),
     $          BP (N+1), OFFP (N+1), ICNTL (20), XSIZE
        LOGICAL PRESRV
        DOUBLE PRECISION
     $          XX (XSIZE)
        
C=== UMD2FB ============================================================
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
C  Find permutations to block triangular form:
C       1) permute the matrix so that it has a zero-free diagonal.
C       2) find strongly-connected components of the corresponding
C          graph.  Each diagonal block corresponds to exactly one
C          strongly-connected component.
C       3) convert the matrix to block triangular form, unless it is
C          to be preserved in its original form.

C  Calls Harwell MA28 routines MC21B and MC13E, which can be obtained
C  separately from Netlib.  Send email to netlib@ornl.gov with the
C  message:
C       send mc13e.f mc21b.f from harwell

C=======================================================================
C  INSTALLATION NOTE:
C=======================================================================
C
C  If the MA28 Harwell Subroutine Library routines MC21B and MC13E
C  (which perform the permutation to block-triangular-form) are not
C  available, then you may comment out all executable code in this
C  routine, or place a "return" statement as the first executable
C  statement (see below).  If you do make this modification, please do
C  not delete any original code.  Add a comment and date to your
C  modifications.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       presrv:         true if original matrix is to be preserved
C       n:              order of matrix
C       nz:             entries in matrix
C       isize:          size of II
C       xsize:          size of XX
C       Cp (1..n+1):    column pointers
C       XX (1..nz):     values
C       II (1..nz):     row indices
C       Icntl:          integer control arguments
C
C          input matrix in column form is in:
C          XX (1..nz), II (1..nz), n, nz, Cp (1..n+1), where
C               II (Cp(col) ... Cp(col+1)-1): row indices
C               XX (Cp(col) ... Cp(col+1)-1): values
C          if presrv is false then xsize and isize must be >= 2*nz
C          otherwise, xsize and isize must be >= nz

C=======================================================================
C  WORKSPACE:
C=======================================================================
C
C       Pr (1..n), Pc (1..n), W (1..n), Zperm (1..n)

C======================================================================
C  OUTPUT: 
C=======================================================================
C
C       nblks: number of blocks
C       if (nblks > 1):
C
C           Cperm (1..n), Rperm (1..n): permutation to block form:
C               Rperm (newrow) = oldrow
C               Cperm (newcol) = oldcol
C       
C           Bp (n-nblks+1...n+1) holds the start/end of blocks 1..nblks
C
C           if (presrv is false) then
C
C              input matrix is converted to block-upper-tri. form,
C              using II/XX (nz+1..2*nz) as workspace.
C              nzdia: nonzeros in diagonal blocks
C              nzoff: nonzeros in off-diagonal blocks
C              (nz = nzdia + nzoff)
C
C              off-diagonal column-oriented form in XX/II (1..nzoff)
C              col is located in
C              XX/II (Offp (col) ... Offp (col+1)-1)
C
C              diagonal blocks now in XX/II (nzoff+1 .. nzoff+nzdia)
C              col is located in
C              XX/II (Cp (col) ... Cp (col+1)-1)
C
C       else, nblks=1: and no other output is generated.

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   UMD2F0
C       subroutines called:     MC21B, MC13E (in MA28 HSL package)

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER COL, NDIAG, I, PO, PB, BLK, P, ROW, K1, K

C  ndiag:   number of entries on the diagonal
C  po:      pointer into off-diagonal part
C  pb:      pointer into diagonal blocks
C  blk:     block number for current diagonal block
C  k1:      column col is in diagonal block A (k1.., k1...)
C  k:       kth row/col in BTF form is Rperm(k)/Cperm(k) in input matrix
C  p:       pointer
C  row:     row index
C  col:     column index
C  i:       general loop index

C=======================================================================
C  EXECUTABLE STATEMENTS:
C       if (MC21B and MC13E not available during installation) return
C=======================================================================

        NZDIA = NZ
        NZOFF = 0

C-----------------------------------------------------------------------
C compute the length of each column
C-----------------------------------------------------------------------

        DO 10 COL = 1, N 
           W (COL) = CP (COL+1) - CP (COL)
10      CONTINUE 

C-----------------------------------------------------------------------
C find a column permutation for a zero-free diagonal
C-----------------------------------------------------------------------

        CALL MC21B (N, II, NZ, CP, W, ZPERM, NDIAG, OFFP, CPERM, PR, PC)
C          MC21B calling interface:
C          input:       n, II (1..nz), nz, Cp (n), W (n):
C                       n-by-n matrix, col is of length W (col),
C                       and its pattern is located in
C                       II (Cp (col) ... Cp (col)+W(col)-1)
C          output:      Zperm (n), the permutation, such that
C                       colold = Zperm (col), and ndiag (number of 
C                       structural nonzeros on the diagonal.    
C                       matrix is structurally singular if ndiag < n
C          workspace:   Offp, Cperm, Pr, Pc

C-----------------------------------------------------------------------
C  permute the columns of the temporary matrix to get zero-free diagonal
C-----------------------------------------------------------------------

        DO 20 COL = 1, N 
           OFFP (COL) = CP (ZPERM (COL))
           W (COL) = CP (ZPERM (COL)+1) - CP (ZPERM (COL))
20      CONTINUE 

C-----------------------------------------------------------------------
C  find a symmetric permutation into upper block triangular form
C  (that is, find the strongly-connected components in the graph).
C-----------------------------------------------------------------------

        CALL MC13E (N, II, NZ, OFFP, W, RPERM, BP, NBLKS, CPERM, PR, PC)
C          MC13E calling interface:
C          input:       n, II (1..nz), nz, Offp (n), W (n)
C                       n-by-n matrix, col of length W(col),
C                       in II (Offp(col) ... Offp(col)+W(col)-1), where
C                       this permuted matrix has a zero-free diagonal
C                       (unless the matrix is structurally singular).
C          output:      Rperm (n), Bp (n+1), nblks
C                       old = Rperm (new) is the symmetric permutation,
C                       there are nblks diagonal blocks, Bp (i) is
C                       the position in new order of the ith block.
C          workspace:   Cperm, Pr, Pc

C-----------------------------------------------------------------------
C  if more than one block, get permutations and block pointers,
C  and convert to block-upper-triangular form (unless matrix preserved)
C-----------------------------------------------------------------------

        IF (NBLKS .NE. 1) THEN 

C          -------------------------------------------------------------
C          find the composite column permutation vector (Cperm):
C          -------------------------------------------------------------

           DO 30 COL = 1, N 
              CPERM (COL) = ZPERM (RPERM (COL))
30         CONTINUE 

C          -------------------------------------------------------------
C          convert to block-upper-triangular form, if not preserved
C          -------------------------------------------------------------

           IF (.NOT. PRESRV) THEN 

C             ----------------------------------------------------------
C             find the inverse permutation vectors, Pr and Pc
C             ----------------------------------------------------------

              DO 40 K = 1, N 
                 PC (CPERM (K)) = K
                 PR (RPERM (K)) = K
40            CONTINUE 

C             ----------------------------------------------------------
C             construct flag array to determine if entry in block or not
C             ----------------------------------------------------------

              BP (NBLKS+1) = N+1
              DO 60 BLK = 1, NBLKS 
                 DO 50 I = BP (BLK), BP (BLK+1)-1 
                    W (I) = BP (BLK)
50               CONTINUE 
60            CONTINUE 

C             ----------------------------------------------------------
C             construct block-diagonal form in XX/II (nz+1..nz+nzdia)
C             ----------------------------------------------------------

C             These blocks are in a permuted order (according to Rperm
C             and Cperm).  The row indices in each block range from 1
C             to the size of the block.

              PB = NZ + 1
              DO 80 COL = 1, N 
                 ZPERM (COL) = PB
                 K1 = W (COL)
CFPP$ NODEPCHK L
                 DO 70 P = CP (CPERM (COL)), CP (CPERM (COL)+1)-1 
                    ROW = PR (II (P)) 
                    IF (W (ROW) .EQ. K1) THEN 
C                      entry is in the diagonal block:
                       II (PB) = ROW - K1 + 1
                       XX (PB) = XX (P)
                       PB = PB + 1
                    ENDIF 
70               CONTINUE 
80            CONTINUE 
C             Zperm (n+1) == pb  ( but Zperm (n+1) does not exist )
              NZDIA = PB - (NZ + 1)
              NZOFF = NZ - NZDIA

C             diagonal blocks now in XX/II (nz+1..nz+nzdia)
C             col is located in XX/II (Zperm (col) ... Zperm (col+1)-1)

C             ----------------------------------------------------------
C             compress original matrix to off-diagonal part, in place
C             ----------------------------------------------------------

C             The rows/cols of off-diagonal form correspond to rows/cols
C             in the original, unpermuted matrix.  They are permuted to
C             the final pivot order and stored in a row-oriented form,
C             after the factorization is complete (by UMD2OF).

              PO = 1
              DO 100 COL = 1, N 
                 OFFP (COL) = PO
                 K1 = W (PC (COL))
CFPP$ NODEPCHK L
                 DO 90 P = CP (COL), CP (COL+1)-1 
                    ROW = PR (II (P))
                    IF (W (ROW) .NE. K1) THEN 
C                      offdiagonal entry
                       II (PO) = II (P)
                       XX (PO) = XX (P)
                       PO = PO + 1
                    ENDIF 
90               CONTINUE 
100           CONTINUE 
              OFFP (N+1) = PO

C             off-diagonal form now in XX/II (1..nzoff)
C             col is located in XX/II(Offp(col)..Offp(col+1)-1)

C             ----------------------------------------------------------
C             move block-diagonal part into place
C             ----------------------------------------------------------

              PB = NZ + 1
CFPP$ NODEPCHK L
              DO 110 I = 0, NZDIA - 1 
                 II (PO+I) = II (PB+I)
                 XX (PO+I) = XX (PB+I)
110           CONTINUE 
              DO 120 COL = 1, N 
                 CP (COL) = ZPERM (COL) - NZDIA
120           CONTINUE 
C             Cp (n+1) == nz+1  ( this is unchanged )

C             diagonal blocks now in XX/II (nzoff+1 .. nzoff+nzdia)
C             col is located in XX/II (Cp (col) ... Cp (col+1)-1)

           ENDIF 

C          -------------------------------------------------------------
C          shift Bp (1 .. nblks+1) down to Bp (1+n-nblks .. n+1), which
C          then becomes the Blkp (1 .. nblks+1) array.
C          -------------------------------------------------------------

           BP (NBLKS+1) = N+1
CFPP$ NODEPCHK L
           DO 130 BLK = NBLKS + 1, 1, -1 
              BP (BLK + (N-NBLKS)) = BP (BLK)
130        CONTINUE 
        ENDIF 

        RETURN
        END 
