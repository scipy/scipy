        SUBROUTINE UMD2S2 (N, JOB, TRANSC, LUXSIZ, LUX,
     $          LUISIZ, LUI, B, X, R, Z, LY, Y, S, CNTL, ICNTL, INFO,
     $          RINFO, CPERM, RPERM, PRESRV, AN, ANZ, AP, AI, AX, ON,
     $          NZOFF, OFFP, OFFI, OFFX, NBLKS, LUBLKP, BLKP, IRSTEP)
        INTEGER N, JOB, LUXSIZ, LUISIZ, LUI (LUISIZ), LY, IRSTEP,
     $          ICNTL (20), INFO (40), CPERM (N), RPERM (N), AN,
     $          ANZ, AP (AN+1), AI (ANZ), ON, NZOFF, OFFP (ON+1),
     $          OFFI (NZOFF), NBLKS, LUBLKP (NBLKS), BLKP (NBLKS+1)
        LOGICAL TRANSC, PRESRV
        DOUBLE PRECISION
     $          LUX (LUXSIZ), B (N), X (N), R (N), Z (N), Y (LY),
     $          S (LY), AX (ANZ), OFFX (NZOFF)
        DOUBLE PRECISION
     $          CNTL (10), RINFO (20)
        
C=== UMD2S2 ============================================================
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
C  Solve a system, given LU factors, permutation arrays, original
C  matrix (if preserved), and off-diagonal blocks (if BTF was used).

C=======================================================================
C  INPUT:
C=======================================================================
C
C       n:              order of matrix
C       job:            0: solve Ax=b, 1: solve Lx=b, 2: solve Ux=b

C       transc:         if true, solve with transposed factors instead

C       LUxsiz:         size of LUx
C       LUx (1..LUxsiz) real values in LU factors for each block
C       LUisiz:         size of LUi
C       LUi (1..LUisiz) integers in LU factors for each block
C       B (1..n):       right-hand-side
C       ly:             size of Y and S, ly=n if Y and S are used
C       Cntl:           control parameters, see UMD21I
C       Icntl:          integer control parameters, see UMD21I
C       Cperm (1..n):   Q, column permutation array
C       Rperm (1..n):   P, row permutation array
C       presrv:         if true, then original matrix was preserved
C       nblks:          number of diagonoal blocks (1 if no BTF)
C       irstep:         maximum number of steps of iterative refinement
C
C       if presrv then
C           an:                 order of preserved matrix, n
C           anz:                number of entries in preserved matrix
C           Ap (1..an+1):       column pointers of preserved matrix
C           Ai (1..anz):        row indices of preserved matrix
C           Ax (1..anz):        values of preserved matrix
C           an, anz, Ap, Ai, Ax:        not accessed
C
C       if nblks > 1 then
C           on:                 n
C           nzoff:              number of off-diagonoal entries
C           Offp (1..n+1)       row pointers for off-diagonal part
C           Offi (1..nzoff):    column indices for off-diagonal part
C           Offx (1..nzoff):    values of off-diagonal part
C           LUBlkp (1..nblks):  pointers to LU factors of each block
C           Blkp (1..nblks+1):  index range of each block
C       else
C           on, nzoff, Offp, Offi, Offx, LUBlkp, Blkp:  not accessed

C=======================================================================
C  WORKSPACE:
C=======================================================================
C
C       R (1..n), Z (1..n)
C       Y (1..ly), S (1..ly):   unaccessed if no iterative refinement

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       X (1..n):       solution
C       Info:           integer informational output, see UMD21I
C       Rinfo:          real informational output, see UMD21I
C
C       if irsteps > 0 and presrv is true then
C           W (1..n):           residual
C           Rinfo (7):  sparse error estimate, omega1
C           Rinfo (8):  sparse error estimate, omega2

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   UMD2SO
C       subroutines called:     UMD2ER, UMD2SL, UMD2LT, UMD2SU, UMD2UT
C       functions called:       IDAMAX, ABS, MAX
        INTRINSIC ABS, MAX
        INTEGER IDAMAX

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER NLU, I, BLK, K1, K2, KN, P, STEP, NPIV, J
        DOUBLE PRECISION
     $          A, AXX, R2, X2, Y2, Z2
        DOUBLE PRECISION
     $          XNORM, TAU, NCTAU, OMEGA1, OMEGA2, D1,
     $          D2, OMEGA, OMLAST, OM1LST, OM2LST, EPS, MAXEPS
        PARAMETER (MAXEPS = 2.0 ** (-15))

C  LU factors:
C  -----------
C  blk:     current diagonal block
C  k1,k2:   current diagonal block is A (k1..k2, k1..k2)
C  kn:      size of diagonal block (= k2-k1+1)
C  nlu:     number of elements in the LU factors of a single diag block
C  npiv:    number of pivots in the LU factors of a single diag block
C
C  Iterative refinement and sparse backward error:
C  -----------------------------------------------
C  step:    number of steps of iterative refinement taken
C  xnorm:   ||x|| maxnorm of solution vector, x
C  tau:     threshold for selecting which estimate to use (1 or 2)
C  nctau:   1000*n*eps
C  eps:     largest positive value such that fl (1.0 + eps) = 1.0
C  omega1:  current sparse backward error estimate 1
C  omega2:  current sparse backward error estimate 2
C  d1:      divisor for omega1
C  d2:      divisor for omega2
C  omega:   omega1 + omega2
C  omlast:  value of omega from previous step
C  om1lst:  value of omega1 from previous step
C  om2lst:  value of omega2 from previous step
C  maxeps:  2**(-15), maximum value that eps is allowed to be
C  a:       value of an entry in A, A_ij
C  axx:     A_ij * x_j
C
C  Other:
C  ------
C  i,j:     loop indices
C  p:       pointer
C  r2:      R (i)
C  x2:      X (i)
C  y2:      Y (i)
C  z2:      Z (i)

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

C-----------------------------------------------------------------------
C  initializations for sparse backward error
C-----------------------------------------------------------------------

        OMEGA = 0
        OMEGA1 = 0
        OMEGA2 = 0
        EPS = CNTL (3)
        IF (EPS .LE. 0 .OR. EPS .GT. MAXEPS) THEN 
C          eps is too small or too big: set to a large default value
           EPS = MAXEPS
        ENDIF 
        NCTAU = 1000 * N * EPS

C-----------------------------------------------------------------------
C  get information on LU factorization if BTF was not used  
C-----------------------------------------------------------------------

        IF (NBLKS .EQ. 1) THEN 
C          p is 1, and LUi (p) is 1
           NLU = LUI (2)
           NPIV = LUI (3)
        ENDIF 

C-----------------------------------------------------------------------
        IF (JOB .EQ. 1) THEN 
C-----------------------------------------------------------------------

C          -------------------------------------------------------------
           IF (.NOT. TRANSC) THEN 
C          -------------------------------------------------------------

C             ----------------------------------------------------------
C             Solve P'Lx=b:  x = L \ Pb
C             ----------------------------------------------------------

              DO 10 I = 1, N 
                 X (I) = B (RPERM (I))
10            CONTINUE 
              IF (NBLKS .EQ. 1) THEN 
                 CALL UMD2SL (NLU, NPIV, N, LUI(6), LUI(6+NLU), LUX,X,Z)
              ELSE 
                 DO 20 BLK = 1, NBLKS 
                    K1 = BLKP (BLK)
                    K2 = BLKP (BLK+1) - 1
                    KN = K2-K1+1
                    IF (KN .GT. 1) THEN 
                       P = LUBLKP (BLK)
                       NLU = LUI (P+1)
                       NPIV = LUI (P+2)
                       CALL UMD2SL (NLU, NPIV, KN, LUI (P+5),
     $                    LUI (P+5+NLU), LUX (LUI (P)), X (K1), Z)
                    ENDIF 
20               CONTINUE 
              ENDIF 

C          -------------------------------------------------------------
           ELSE 
C          -------------------------------------------------------------

C             ----------------------------------------------------------
C             Solve L'Px=b:  x = P' (L' \ b)
C             ----------------------------------------------------------

              DO 30 I = 1, N 
                 R (I) = B (I)
30            CONTINUE 
              IF (NBLKS .EQ. 1) THEN 
                 CALL UMD2LT (NLU, NPIV, N, LUI(6), LUI(6+NLU), LUX,R,Z)
              ELSE 
                 DO 40 BLK = 1, NBLKS 
                    K1 = BLKP (BLK)
                    K2 = BLKP (BLK+1) - 1
                    KN = K2-K1+1
                    IF (KN .GT. 1) THEN 
                       P = LUBLKP (BLK)
                       NLU = LUI (P+1)
                       NPIV = LUI (P+2)
                       CALL UMD2LT (NLU, NPIV, KN, LUI (P+5),
     $                    LUI (P+5+NLU), LUX (LUI (P)), R (K1), Z)
                    ENDIF 
40               CONTINUE 
              ENDIF 
              DO 50 I = 1, N 
                 X (RPERM (I)) = R (I)
50            CONTINUE 

C          -------------------------------------------------------------
           ENDIF 
C          -------------------------------------------------------------

C-----------------------------------------------------------------------
        ELSE IF (JOB .EQ. 2) THEN 
C-----------------------------------------------------------------------

C          -------------------------------------------------------------
           IF (TRANSC) THEN 
C          -------------------------------------------------------------

C             ----------------------------------------------------------
C             Solve QU'x=b:  x = U' \ Q'b
C             ----------------------------------------------------------

              DO 60 I = 1, N 
                 X (I) = B (CPERM (I))
60            CONTINUE 
              IF (NBLKS .EQ. 1) THEN 
                 CALL UMD2UT (NLU, NPIV, N, LUI(6), LUI(6+NLU), LUX,X,Z)
              ELSE 
                 DO 100 BLK = 1, NBLKS 
                    K1 = BLKP (BLK)
                    K2 = BLKP (BLK+1) - 1
                    KN = K2-K1+1
                    IF (KN .EQ. 1) THEN 
                       X (K1) = X (K1) / LUX (LUBLKP (BLK))
                       R (K1) = X (K1)
                    ELSE 
                       P = LUBLKP (BLK)
                       NLU = LUI (P+1)
                       NPIV = LUI (P+2)
                       CALL UMD2UT (NLU, NPIV, KN, LUI (P+5),
     $                    LUI (P+5+NLU), LUX (LUI (P)), X (K1), Z)
                       DO 70 I = K1, K2 
                          R (I) = X (I)
70                     CONTINUE 
                       CALL UMD2LT (NLU, NPIV, KN, LUI (P+5),
     $                    LUI (P+5+NLU), LUX (LUI (P)), R (K1), Z)
                    ENDIF 
                    DO 90 I = K1, K2 
                       R2 = R (I)
                       DO 80 P = OFFP (I), OFFP (I+1)-1 
                          X (OFFI (P)) = X (OFFI (P)) - OFFX (P) * R2
80                     CONTINUE 
90                  CONTINUE 
100              CONTINUE 
              ENDIF 

C          -------------------------------------------------------------
           ELSE 
C          -------------------------------------------------------------

C             ----------------------------------------------------------
C             Solve UQ'x=b:  x = Q (U \ b)
C             ----------------------------------------------------------

              IF (NBLKS .EQ. 1) THEN 
                 DO 110 I = 1, N 
                    R (I) = B (I)
110              CONTINUE 
                 CALL UMD2SU (NLU, NPIV, N, LUI(6), LUI(6+NLU), LUX,R,Z)
              ELSE 
                 DO 150 BLK = NBLKS, 1, -1 
                    K1 = BLKP (BLK)
                    K2 = BLKP (BLK+1) - 1
                    KN = K2-K1+1
                    DO 130 I = K1, K2 
                       X2 = 0
                       DO 120 P = OFFP (I), OFFP (I+1)-1 
                          X2 = X2 + OFFX (P) * R (OFFI (P))
120                    CONTINUE 
                       X (I) = X2
130                 CONTINUE 
                    IF (KN .EQ. 1) THEN 
                       R (K1) = (B (K1) - X (K1)) / LUX (LUBLKP (BLK))
                    ELSE 
                       P = LUBLKP (BLK)
                       NLU = LUI (P+1)
                       NPIV = LUI (P+2)
                       CALL UMD2SL (NLU, NPIV, KN, LUI (P+5),
     $                    LUI (P+5+NLU), LUX (LUI (P)), X (K1), Z)
                       DO 140 I = K1, K2 
                          R (I) = B (I) - X (I)
140                    CONTINUE 
                       CALL UMD2SU (NLU, NPIV, KN, LUI (P+5),
     $                    LUI (P+5+NLU), LUX (LUI (P)), R (K1), Z)
                    ENDIF 
150              CONTINUE 
              ENDIF 
              DO 160 I = 1, N 
                 X (CPERM (I)) = R (I)
160           CONTINUE 

C          -------------------------------------------------------------
           ENDIF 
C          -------------------------------------------------------------

C-----------------------------------------------------------------------
        ELSE 
C-----------------------------------------------------------------------

           DO 450 STEP = 0, IRSTEP 

C             ----------------------------------------------------------
C             If transa was true in UMD2FA or UMD2RF, then C = A'.
C             Otherwise C = A.  In both cases, the factorization is
C             PCQ = LU, and C is stored in column-form in Ai,Ax,Ap if
C             it is preserved.
C             ----------------------------------------------------------

C             ----------------------------------------------------------
              IF (.NOT. TRANSC) THEN 
C             ----------------------------------------------------------

C                -------------------------------------------------------
C                Solve Cx=b (step 0):
C                   x = Q (U \ L \ Pb)
C                and then perform iterative refinement (step > 0):
C                   x = x + Q (U \ L \ P (b-Cx))
C                -------------------------------------------------------

                 IF (STEP .EQ. 0) THEN 
                    DO 170 I = 1, N 
                       R (I) = B (RPERM (I))
170                 CONTINUE 
                 ELSE 
                    DO 180 I = 1, N 
                       Z (I) = B (I)
180                 CONTINUE 
                    DO 200 I = 1, N 
                       X2 = X (I)
                       DO 190 P = AP (I), AP (I+1) - 1 
                          Z (AI (P)) = Z (AI (P)) - AX (P) * X2
190                    CONTINUE 
200                 CONTINUE 
                    DO 210 I = 1, N 
                       R (I) = Z (RPERM (I))
210                 CONTINUE 
                 ENDIF 
                 IF (NBLKS .EQ. 1) THEN 
                    CALL UMD2SL (NLU, NPIV, N,LUI(6),LUI(6+NLU),LUX,R,Z)
                    CALL UMD2SU (NLU, NPIV, N,LUI(6),LUI(6+NLU),LUX,R,Z)
                 ELSE 
                    DO 240 BLK = NBLKS, 1, -1 
                       K1 = BLKP (BLK)
                       K2 = BLKP (BLK+1) - 1
                       KN = K2-K1+1
                       DO 230 I = K1, K2 
                          R2 = R (I)
                          DO 220 P = OFFP (I), OFFP (I+1)-1 
                             R2 = R2 - OFFX (P) * R (OFFI (P))
220                       CONTINUE 
                          R (I) = R2
230                    CONTINUE 
                       IF (KN .EQ. 1) THEN 
                          R (K1) = R (K1) / LUX (LUBLKP (BLK))
                       ELSE 
                          P = LUBLKP (BLK)
                          NLU = LUI (P+1)
                          NPIV = LUI (P+2)
                          CALL UMD2SL (NLU, NPIV, KN, LUI (P+5),
     $                       LUI (P+5+NLU), LUX (LUI (P)), R (K1), Z)
                          CALL UMD2SU (NLU, NPIV, KN, LUI (P+5),
     $                       LUI (P+5+NLU), LUX (LUI (P)), R (K1), Z)
                       ENDIF 
240                 CONTINUE 
                 ENDIF 
                 IF (STEP .EQ. 0) THEN 
                    DO 250 I = 1, N 
                       X (CPERM (I)) = R (I)
250                 CONTINUE 
                 ELSE 
                    DO 260 I = 1, N 
                       X (CPERM (I)) = X (CPERM (I)) + R (I)
260                 CONTINUE 
                 ENDIF 

C             ----------------------------------------------------------
              ELSE 
C             ----------------------------------------------------------

C                -------------------------------------------------------
C                Solve C'x=b (step 0):
C                   x = P' (L' \ U' \ Q'b)
C                and then perform iterative refinement (step > 0):
C                   x = x + P' (L' \ U' \ Q' (b-C'x))
C                -------------------------------------------------------

                 IF (STEP .EQ. 0) THEN 
                    DO 270 I = 1, N 
                       R (I) = B (CPERM (I))
270                 CONTINUE 
                 ELSE 
                    DO 280 I = 1, N 
                       Z (I) = B (I)
280                 CONTINUE 
                    DO 300 I = 1, N 
                       Z2 = Z (I)
                       DO 290 P = AP (I), AP (I+1) - 1 
                          Z2 = Z2 - AX (P) * X (AI (P))
290                    CONTINUE 
                       Z (I) = Z2
300                 CONTINUE 
                    DO 310 I = 1, N 
                       R (I) = Z (CPERM (I))
310                 CONTINUE 
                 ENDIF 
                 IF (NBLKS .EQ. 1) THEN 
                    CALL UMD2UT (NLU, NPIV, N,LUI(6),LUI(6+NLU),LUX,R,Z)
                    CALL UMD2LT (NLU, NPIV, N,LUI(6),LUI(6+NLU),LUX,R,Z)
                 ELSE 
                    DO 340 BLK = 1, NBLKS 
                       K1 = BLKP (BLK)
                       K2 = BLKP (BLK+1) - 1
                       KN = K2-K1+1
                       IF (KN .EQ. 1) THEN 
                          R (K1) = R (K1) / LUX (LUBLKP (BLK))
                       ELSE 
                          P = LUBLKP (BLK)
                          NLU = LUI (P+1)
                          NPIV = LUI (P+2)
                          CALL UMD2UT (NLU, NPIV, KN, LUI (P+5),
     $                       LUI (P+5+NLU), LUX (LUI (P)), R (K1), Z)
                          CALL UMD2LT (NLU, NPIV, KN, LUI (P+5),
     $                       LUI (P+5+NLU), LUX (LUI (P)), R (K1), Z)
                       ENDIF 
                       DO 330 I = K1, K2 
                          R2 = R (I)
                          DO 320 P = OFFP (I), OFFP (I+1)-1 
                             R (OFFI (P)) = R (OFFI (P)) - OFFX (P) * R2
320                       CONTINUE 
330                    CONTINUE 
340                 CONTINUE 
                 ENDIF 
                 IF (STEP .EQ. 0) THEN 
                    DO 350 I = 1, N 
                       X (RPERM (I)) = R (I)
350                 CONTINUE 
                 ELSE 
                    DO 360 I = 1, N 
                       X (RPERM (I)) = X (RPERM (I)) + R (I)
360                 CONTINUE 
                 ENDIF 

C             ----------------------------------------------------------
              ENDIF 
C             ----------------------------------------------------------

C             ----------------------------------------------------------
C             sparse backward error estimate
C             ----------------------------------------------------------

              IF (IRSTEP .GT. 0) THEN 

C                xnorm = ||x|| maxnorm
                 XNORM = ABS (X (IDAMAX (N, X, 1)))

C                r (i) = (b-Ax)_i, residual (or A')
C                z (i) = (|A||x|)_i
C                y (i) = ||A_i||, maxnorm of row i of A (or A')
                 DO 370 I = 1, N 
                    R (I) = B (I)
                    Z (I) = 0
                    Y (I) = 0
370              CONTINUE 

                 IF (.NOT. TRANSC) THEN 

C                   ----------------------------------------------------
C                   sparse backward error for Cx=b, C stored by column
C                   ----------------------------------------------------

                    DO 390 J = 1, N 
                       X2 = X (J)
CFPP$ NODEPCHK L
                       DO 380 P = AP (J), AP (J+1) - 1 
                          I = AI (P)
                          A = AX (P)
                          AXX = A * X2
                          R (I) = R (I) -     (AXX)
                          Z (I) = Z (I) + ABS (AXX)
                          Y (I) = Y (I) + ABS (A)
380                    CONTINUE 
390                 CONTINUE 

                 ELSE 

C                   ----------------------------------------------------
C                   sparse backward error for C'x=b, C' stored by row
C                   ----------------------------------------------------

                    DO 410 I = 1, N 
                       R2 = R (I)
                       Z2 = Z (I)
                       Y2 = Y (I)
CFPP$ NODEPCHK L
                       DO 400 P = AP (I), AP (I+1) - 1 
                          J = AI (P)
                          A = AX (P)
                          AXX = A * X (J)
                          R2 = R2 -     (AXX)
                          Z2 = Z2 + ABS (AXX)
                          Y2 = Y2 + ABS (A)
400                    CONTINUE 
                       R (I) = R2
                       Z (I) = Z2
                       Y (I) = Y2
410                 CONTINUE 

                 ENDIF 

C                -------------------------------------------------------
C                save the last iteration in case we need to reinstate it
C                -------------------------------------------------------

                 OMLAST = OMEGA
                 OM1LST = OMEGA1
                 OM2LST = OMEGA2

C                -------------------------------------------------------
C                compute sparse backward errors: omega1 and omega2
C                -------------------------------------------------------

                 OMEGA1 = 0
                 OMEGA2 = 0
                 DO 420 I = 1, N 
                    TAU = (Y (I) * XNORM + ABS (B (I))) * NCTAU
                    D1 = Z (I) + ABS (B (I))
                    IF (D1 .GT. TAU) THEN 
                       OMEGA1 = MAX (OMEGA1, ABS (R (I)) / D1)
                    ELSE IF (TAU .GT. 0) THEN 
                       D2 = Z (I) + Y (I) * XNORM
                       OMEGA2 = MAX (OMEGA2, ABS (R (I)) / D2)
                    ENDIF 
420              CONTINUE 
                 OMEGA = OMEGA1 + OMEGA2
                 RINFO (7) = OMEGA1
                 RINFO (8) = OMEGA2

C                -------------------------------------------------------
C                stop the iterations if the backward error is small
C                -------------------------------------------------------

                 INFO (24) = STEP
                 IF (1 + OMEGA .LE. 1) THEN 
C                   further iterative refinement will no longer improve
C                   the solution
                    RETURN
                 ENDIF 

C                -------------------------------------------------------
C                stop if insufficient decrease in omega
C                -------------------------------------------------------

                 IF (STEP .GT. 0 .AND. OMEGA .GT. OMLAST / 2) THEN 
                    IF (OMEGA .GT. OMLAST) THEN 
C                      last iteration better than this one, reinstate it
                       DO 430 I = 1, N 
                          X (I) = S (I)
                          RINFO (7) = OM1LST
                          RINFO (8) = OM2LST
430                    CONTINUE 
                    ENDIF 
                    INFO (24) = STEP - 1
                    RETURN
                 ENDIF 

C                -------------------------------------------------------
C                save current solution in case we need to reinstate
C                -------------------------------------------------------

                 DO 440 I = 1, N 
                    S (I) = X (I)
440              CONTINUE 

              ENDIF 

450        CONTINUE 

C-----------------------------------------------------------------------
        ENDIF 
C-----------------------------------------------------------------------

        RETURN
        END 
