        SUBROUTINE UMD2SU (NLU, NPIV, N, LUP, LUI, LUX, X, W)
        INTEGER NLU, NPIV, N, LUP (NLU), LUI (*)
        DOUBLE PRECISION
     $          LUX (*), X (N), W (N)
        
C=== UMD2SU ============================================================
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
C  solves Ux = b, where U is the upper triangular factor of a matrix
C  (if BTF not used) or a single diagonal block (if BTF is used).
C  B is overwritten with the solution X.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       nlu:            number of LU arrowheads in the LU factors
C       npiv:           number of pivots found (normally n)
C       n:              order of matrix
C       LUp (1..nlu):   pointer to LU arrowheads in LUi
C       LUi ( ... ):    integer values of LU arrowheads
C       LUx ( ... ):    real values of LU arroheads
C       X (1..n):       the right-hand-side

C=======================================================================
C  WORKSPACE:
C=======================================================================
C
C       W (1..n)

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       X (1..n):       the solution to Ux=b

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   UMD2S2
C       subroutines called:     DTRSV, DGEMV

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER J, K, S, LUIP, LUXP, LUK, LUDEGR, LUDEGC, LURP, UXP,
     $          LUCP, COL
        DOUBLE PRECISION
     $          ONE

C  s:       an element, or LU arrowhead
C  k:       kth pivot
C  j:       jth column in U2 array in element s
C  luip:    s is in LUi (luip...)
C  luxp:    real part of s is in LUx (luxp...)
C  luk:     number of pivots in s
C  ludegc:  column degree of non-pivotal part of s
C  ludegr:  row degree of non-pivotal part of s
C  lucp:    pattern of column of s in LUi (lucp...lucp+ludegc-1)
C  lurp:    pattern of row of s in LUi (lurp...lurp+ludegr-1)
C  uxp:     the luk-by-ludegr U2 block of s is in LUx (uxp...)
C  col:     column index

C=======================================================================
C  EXECUTABLE STATMENTS:
C=======================================================================

        ONE = 1
        K = NPIV
        DO 30 S = NLU, 1, -1 

C          -------------------------------------------------------------
C          get s-th LU arrowhead (s = nlu..1, in reverse pivotal order)
C          -------------------------------------------------------------

           LUIP   = LUP (S)
           LUXP   = LUI (LUIP)
           LUK    = LUI (LUIP+1)
           LUDEGR = LUI (LUIP+2)
           LUDEGC = LUI (LUIP+3)
           LUCP   = (LUIP + 7)
           LURP   = LUCP + LUDEGC
           UXP    = LUXP + LUK * (LUDEGC + LUK)

           IF (LUK .EQ. 1) THEN 

C             ----------------------------------------------------------
C             only one pivot, stride-1 sparse dot product
C             ----------------------------------------------------------

CFPP$ NODEPCHK L
              DO 10 J = 1, LUDEGR 
                 COL = LUI (LURP+J-1)
C                row: k, U (row,col): LUx (uxp+j-1)
                 X (K) = X (K) - LUX (UXP+J-1) * X (COL)
10            CONTINUE 
C             divide by pivot, U (k,k): LUx (luxp)
              X (K) = X (K) / LUX (LUXP)
              K = K - 1

           ELSE 

C             ----------------------------------------------------------
C             more than one pivot
C             ----------------------------------------------------------

              K = K - LUK
              DO 20 J = 1, LUDEGR 
                 COL = LUI (LURP+J-1)
                 W (J) = X (COL)
20            CONTINUE 
              CALL DGEMV ('N', LUK, LUDEGR, -ONE,
     $           LUX (UXP), LUK, W, 1, ONE, X (K+1), 1)
              CALL DTRSV ('U', 'N', 'N', LUK,
     $           LUX (LUXP), LUDEGC + LUK, X (K+1), 1)

           ENDIF 

30      CONTINUE 
        RETURN
        END 
