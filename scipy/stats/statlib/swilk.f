      SUBROUTINE SWILK (INIT, X, N, N1, N2, A, W, PW, IFAULT)
C
C        ALGORITHM AS R94 APPL. STATIST. (1995) VOL.44, NO.4
C
C        Calculates the Shapiro-Wilk W test and its significance level
C
C        IFAULT error code details from the R94 paper:
C        - 0 for no fault
C        - 1 if N1 < 3
C        - 2 if N > 5000 (a non-fatal error)
C        - 3 if N2 < N/2, so insufficient storage for A
C        - 4 if N1 > N or (N1 < N and N < 20)
C        - 5 if the proportion censored (N-N1)/N > 0.8
C        - 6 if the data have zero range (if sorted on input)
C
      INTEGER N, N1, N2, IFAULT
      REAL X(*), A(*), PW, W
      REAL C1(6), C2(6), C3(4), C4(4), C5(4), C6(3), C7(2)
      REAL C8(2), C9(2), G(2)
      REAL Z90, Z95, Z99, ZM, ZSS, BF1, XX90, XX95, ZERO, ONE, TWO
      REAL THREE, SQRTH, QTR, TH, SMALL, PI6, STQR
      REAL SUMM2, SSUMM2, FAC, RSN, AN, AN25, A1, A2, DELTA, RANGE
      REAL SA, SX, SSX, SSA, SAX, ASA, XSX, SSASSX, W1, Y, XX, XI
      REAL GAMMA, M, S, LD, BF, Z90F, Z95F, Z99F, ZFM, ZSD, ZBAR
C      
C        Auxiliary routines
C
      REAL PPND, POLY
      DOUBLE PRECISION ALNORM
C      
      INTEGER NCENS, NN2, I, I1, J
      LOGICAL INIT, UPPER
C
      DATA C1 /0.0E0, 0.221157E0, -0.147981E0, -0.207119E1,
     *     0.4434685E1, -0.2706056E1/
      DATA C2 /0.0E0, 0.42981E-1, -0.293762E0, -0.1752461E1,
     *     0.5682633E1, -0.3582633E1/
      DATA C3 /0.5440E0, -0.39978E0, 0.25054E-1, -0.6714E-3/
      DATA C4 /0.13822E1, -0.77857E0, 0.62767E-1, -0.20322E-2/
      DATA C5 /-0.15861E1, -0.31082E0, -0.83751E-1, 0.38915E-2/
      DATA C6 /-0.4803E0, -0.82676E-1, 0.30302E-2/
      DATA C7 /0.164E0, 0.533E0/
      DATA C8 /0.1736E0, 0.315E0/
      DATA C9 /0.256E0, -0.635E-2/
      DATA G  /-0.2273E1, 0.459E0/
      DATA Z90, Z95, Z99 /0.12816E1, 0.16449E1, 0.23263E1/
      DATA ZM, ZSS /0.17509E1, 0.56268E0/
      DATA BF1 /0.8378E0/, XX90, XX95 /0.556E0, 0.622E0/
      DATA ZERO /0.0E0/, ONE/1.0E0/, TWO/2.0E0/, THREE/3.0E0/
      DATA SQRTH /0.70711E0/, QTR/0.25E0/, TH/0.375E0/, SMALL/1E-19/
      DATA PI6 /0.1909859E1/, STQR/0.1047198E1/, UPPER/.TRUE./
C
      PW  =  ONE
      IF (W .GE. ZERO) W = ONE
      AN = N
      IFAULT = 3
      NN2 = N/2
      IF (N2 .LT. NN2) RETURN
      IFAULT = 1
      IF (N .LT. 3) RETURN
C
C        If INIT is false, calculates coefficients for the test
C
      IF (.NOT. INIT) THEN
        IF (N .EQ. 3) THEN
           A(1) = SQRTH
        ELSE
           AN25 = AN + QTR
           SUMM2 = ZERO
           DO 30 I = 1, N2
              A(I) = PPND((REAL(I) - TH)/AN25,IFAULT)
              SUMM2 = SUMM2 + A(I) ** 2
30          CONTINUE                
           SUMM2 = SUMM2 * TWO
           SSUMM2 = SQRT(SUMM2)
           RSN = ONE / SQRT(AN)
           A1 = POLY(C1, 6, RSN) - A(1) / SSUMM2
C
C        Normalize coefficients
C
           IF (N .GT. 5) THEN
              I1 = 3
              A2 = -A(2)/SSUMM2 + POLY(C2,6,RSN)
              FAC = SQRT((SUMM2 - TWO * A(1) ** 2 - TWO *
     *               A(2) ** 2)/(ONE - TWO * A1 ** 2 - TWO * A2 ** 2))
              A(1) = A1
              A(2) = A2
           ELSE
              I1 = 2
              FAC = SQRT((SUMM2 - TWO * A(1) ** 2)/
     *                   (ONE - TWO * A1 ** 2))
              A(1) = A1
           END IF
           DO 40 I = I1, NN2
              A(I) = -A(I)/FAC
   40       CONTINUE
        END IF
        INIT = .TRUE.
      END IF
      IF (N1 .LT. 3) RETURN
      NCENS = N - N1
      IFAULT = 4
      IF (NCENS .LT. 0 .OR. (NCENS .GT. 0 .AND. N .LT. 20)) RETURN
      IFAULT = 5
      DELTA = FLOAT(NCENS)/AN
      IF (DELTA .GT. 0.8) RETURN
C
C        If W input as negative, calculate significance level of -W
C
      IF (W .LT. ZERO) THEN
        W1 = ONE + W
        IFAULT = 0
        GOTO 70
      END IF
C
C        Check for zero range
C
      IFAULT = 6
      RANGE = X(N1) - X(1)
      IF (RANGE .LT. SMALL) RETURN
C
C        Check for correct sort order on range - scaled X
C
      IFAULT = 7
      XX = X(1)/RANGE
      SX = XX
      SA = -A(1)
      J = N - 1
      DO 50 I = 2, N1
        XI = X(I)/RANGE
CCCCC   IF (XX-XI .GT. SMALL) PRINT *,' ANYTHING'
        SX = SX + XI
        IF (I .NE. J) SA = SA + SIGN(1, I - J) * A(MIN(I, J))
        XX = XI
        J = J - 1
50    CONTINUE
      IFAULT = 0
      IF (N .GT. 5000) IFAULT = 2
C
C        Calculate W statistic as squared correlation
C        between data and coefficients
C
      SA = SA/N1
      SX = SX/N1
      SSA = ZERO
      SSX = ZERO
      SAX = ZERO
      J = N
      DO 60 I = 1, N1
        IF (I .NE. J) THEN
           ASA = SIGN(1, I - J) * A(MIN(I, J)) - SA
        ELSE
           ASA = -SA
        END IF
        XSX = X(I)/RANGE - SX
        SSA = SSA + ASA * ASA
        SSX = SSX + XSX * XSX
        SAX = SAX + ASA * XSX
        J = J - 1
   60 CONTINUE
C
C        W1 equals (1-W) calculated to avoid excessive rounding error
C        for W very near 1 (a potential problem in very large samples)
C
      SSASSX = SQRT(SSA * SSX)
      W1 = (SSASSX - SAX) * (SSASSX + SAX)/(SSA * SSX)
   70 W = ONE - W1
C
C        Calculate significance level for W (exact for N=3)
C
      IF (N .EQ. 3) THEN
         PW = PI6 * (ASIN(SQRT(W)) - STQR)
         RETURN
      END IF
      Y = LOG(W1)
      XX = LOG(AN)
      M = ZERO
      S = ONE
      IF (N .LE. 11) THEN
        GAMMA = POLY(G, 2, AN)
        IF (Y .GE. GAMMA) THEN
           PW = SMALL
           RETURN
        END IF
        Y = -LOG(GAMMA - Y)
        M = POLY(C3, 4, AN)
        S = EXP(POLY(C4, 4, AN))
      ELSE
        M = POLY(C5, 4, XX)
        S = EXP(POLY(C6, 3, XX))
      END IF
      IF (NCENS .GT. 0) THEN
C
C        Censoring by proportion NCENS/N.  Calculate mean and sd
C        of normal equivalent deviate of W.
C
        LD = -LOG(DELTA)
        BF = ONE + XX * BF1
        Z90F = Z90 + BF * POLY(C7, 2, XX90 ** XX) ** LD
        Z95F = Z95 + BF * POLY(C8, 2, XX95 ** XX) ** LD
        Z99F = Z99 + BF * POLY(C9, 2, XX) ** LD
C
C        Regress Z90F,...,Z99F on normal deviates Z90,...,Z99 to get
C        pseudo-mean and pseudo-sd of z as the slope and intercept
C
        ZFM = (Z90F + Z95F + Z99F)/THREE
        ZSD = (Z90*(Z90F-ZFM)+Z95*(Z95F-ZFM)+Z99*(Z99F-ZFM))/ZSS
        ZBAR = ZFM - ZSD * ZM
        M = M + ZBAR * S
        S = S * ZSD
      END IF
      PW = REAL(ALNORM(DBLE((Y - M)/S), UPPER))
C
      RETURN
      END

      DOUBLE PRECISION FUNCTION ALNORM(X, UPPER)
C
C       EVALUATES THE TAIL AREA OF THE STANDARDIZED NORMAL CURVE FROM
C       X TO INFINITY IF UPPER IS .TRUE. OR FROM MINUS INFINITY TO X
C       IF UPPER IS .FALSE.
C
C  NOTE NOVEMBER 2001: MODIFY UTZERO.  ALTHOUGH NOT NECESSARY
C  WHEN USING ALNORM FOR SIMPLY COMPUTING PERCENT POINTS,
C  EXTENDING RANGE IS HELPFUL FOR USE WITH FUNCTIONS THAT
C  USE ALNORM IN INTERMEDIATE COMPUTATIONS.
C
      DOUBLE PRECISION LTONE,UTZERO,ZERO,HALF,ONE,CON,
     $ A1,A2,A3,A4,A5,A6,A7,B1,B2,
     $ B3,B4,B5,B6,B7,B8,B9,B10,B11,B12,X,Y,Z,ZEXP
      LOGICAL UPPER,UP
C
C       LTONE AND UTZERO MUST BE SET TO SUIT THE PARTICULAR COMPUTER
C
CCCCC DATA LTONE, UTZERO /7.0D0, 18.66D0/
      DATA LTONE, UTZERO /7.0D0, 38.00D0/
      DATA ZERO,HALF,ONE,CON /0.0D0,0.5D0,1.0D0,1.28D0/
      DATA          A1,             A2,            A3,
     $              A4,             A5,            A6,
     $              A7
     $ /0.398942280444D0, 0.399903438504D0, 5.75885480458D0,
     $   29.8213557808D0,  2.62433121679D0, 48.6959930692D0,
     $   5.92885724438D0/
      DATA          B1,             B2,             B3,
     $              B4,             B5,             B6,
     $              B7,             B8,             B9,
     $             B10,            B11,            B12
     $ /0.398942280385D0,      3.8052D-8,    1.00000615302D0,
     $   3.98064794D-4,     1.98615381364D0, 0.151679116635D0,
     $   5.29330324926D0,   4.8385912808D0,  15.1508972451D0,
     $  0.742380924027D0,   30.789933034D0,  3.99019417011D0/
C
      ZEXP(Z) = DEXP(Z)
C
      UP = UPPER
      Z = X
      IF (Z .GE. ZERO) GOTO 10
      UP = .NOT. UP
      Z = -Z
  10  IF (Z .LE. LTONE .OR. UP .AND. Z .LE. UTZERO) GOTO 20
      ALNORM = ZERO
      GOTO 40
  20  Y = HALF * Z * Z
      IF (Z .GT. CON) GOTO 30
C
      ALNORM = HALF - Z * (A1- A2 * Y / (Y + A3- A4 / (Y + A5 + A6 /
     $ (Y + A7))))
      GOTO 40
C
  30  ALNORM = B1* ZEXP(-Y)/(Z - B2 + B3/ (Z +B4 +B5/(Z -B6 +B7/
     $ (Z +B8 -B9/ (Z +B10 +B11/ (Z + B12))))))
C
  40  IF (.NOT. UP) ALNORM = ONE - ALNORM
      RETURN
      END

      REAL FUNCTION PPND(P, IFAULT)
C
C  ALGORITHM AS 111  APPL. STATIST. (1977), VOL.26, NO.1
C
C  PRODUCES NORMAL DEVIATE CORRESPONDING TO LOWER TAIL AREA OF P
C  REAL VERSION FOR EPS = 2 **(-31)
C  THE HASH SUMS ARE THE SUMS OF THE MODULI OF THE COEFFICIENTS
C  THEY HAVE NO INHERENT MEANINGS BUT ARE INCLUDED FOR USE IN
C  CHECKING TRANSCRIPTIONS
C  STANDARD FUNCTIONS ABS, ALOG AND SQRT ARE USED
C
C  NOTE: WE COULD USE DATAPLOT NORPPF, BUT VARIOUS APPLIED
C        STATISTICS ALGORITHMS USE THIS.  SO WE PROVIDE IT TO
C        MAKE USE OF APPLIED STATISTICS ALGORITHMS EASIER.
C
      REAL ZERO, SPLIT, HALF, ONE
      REAL A0, A1, A2, A3, B1, B2, B3, B4, C0, C1, C2, C3, D1, D2
      REAL P, Q, R
      INTEGER IFAULT
      DATA ZERO /0.0E0/, HALF/0.5E0/, ONE/1.0E0/
      DATA SPLIT /0.42E0/
      DATA A0 / 2.50662823884E0/
      DATA A1 / -18.61500062529E0/
      DATA A2 / 41.39119773534E0/
      DATA A3 / -25.44106049637E0/
      DATA B1 / -8.47351093090E0/
      DATA B2 / 23.08336743743E0/
      DATA B3 / -21.06224101826E0/
      DATA B4 / 3.13082909833E0/
      DATA C0 / -2.78718931138E0/
      DATA C1 / -2.29796479134E0/
      DATA C2 / 4.85014127135E0/
      DATA C3 / 2.32121276858E0/
      DATA D1 / 3.54388924762E0/
      DATA D2 / 1.63706781897E0/
C
      IFAULT = 0
      Q = P - HALF
      IF (ABS(Q) .GT. SPLIT) GOTO 1
      R = Q*Q
      PPND = Q * (((A3*R + A2)*R + A1) * R + A0) /
     *  ((((B4*R + B3)*R + B2) * R + B1) * R + ONE)
      RETURN
1     R = P
      IF (Q .GT. ZERO)R = ONE - P
      IF (R .LE. ZERO) GOTO 2
      R = SQRT(-ALOG(R))
      PPND = (((C3 * R + C2) * R + C1) * R + C0)/
     *  ((D2*R + D1) * R + ONE)
      IF (Q .LT. ZERO) PPND = -PPND
      RETURN
2     IFAULT = 1
      PPND = ZERO
      RETURN
      END

      REAL FUNCTION POLY(C, NORD, X)
C
C
C        ALGORITHM AS 181.2   APPL. STATIST.  (1982) VOL. 31, NO. 2
C
C        CALCULATES THE ALGEBRAIC POLYNOMIAL OF ORDER NORED-1 WITH
C        ARRAY OF COEFFICIENTS C.  ZERO ORDER COEFFICIENT IS C(1)
C
      REAL C(NORD)
      POLY = C(1)
      IF(NORD.EQ.1) RETURN
      P = X*C(NORD)
      IF(NORD.EQ.2) GOTO 20
      N2 = NORD-2
      J = N2+1
      DO 10 I = 1,N2
      P = (P+C(J))*X
      J = J-1
   10 CONTINUE
   20 POLY = POLY+P
      RETURN
      END
