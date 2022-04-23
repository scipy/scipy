      SUBROUTINE cdfchi(which,p,q,x,df,status,bound)
C**********************************************************************
C
C      SUBROUTINE CDFCHI( WHICH, P, Q, X, DF, STATUS, BOUND )
C               Cumulative Distribution Function
C               CHI-Square distribution
C
C
C                              Function
C
C
C     Calculates any one parameter of the chi-square
C     distribution given values for the others.
C
C
C                              Arguments
C
C
C     WHICH --> Integer indicating which of the next three argument
C               values is to be calculated from the others.
C               Legal range: 1..3
C               iwhich = 1 : Calculate P and Q from X and DF
C               iwhich = 2 : Calculate X from P,Q and DF
C               iwhich = 3 : Calculate DF from P,Q and X
C                    INTEGER WHICH
C
C     P <--> The integral from 0 to X of the chi-square
C            distribution.
C            Input range: [0, 1].
C                    DOUBLE PRECISION P
C
C     Q <--> 1-P.
C            Input range: (0, 1].
C            P + Q = 1.0.
C                    DOUBLE PRECISION Q
C
C     X <--> Upper limit of integration of the non-central
C            chi-square distribution.
C            Input range: [0, +infinity).
C            Search range: [0,1E100]
C                    DOUBLE PRECISION X
C
C     DF <--> Degrees of freedom of the
C             chi-square distribution.
C             Input range: (0, +infinity).
C             Search range: [ 1E-100, 1E100]
C                    DOUBLE PRECISION DF
C
C     STATUS <-- 0 if calculation completed correctly
C               -I if input parameter number I is out of range
C                1 if answer appears to be lower than lowest
C                  search bound
C                2 if answer appears to be higher than greatest
C                  search bound
C                3 if P + Q .ne. 1
C               10 indicates error returned from cumgam.  See
C                  references in cdfgam
C                    INTEGER STATUS
C
C     BOUND <-- Undefined if STATUS is 0
C
C               Bound exceeded by parameter number I if STATUS
C               is negative.
C
C               Lower search bound if STATUS is 1.
C
C               Upper search bound if STATUS is 2.
C
C
C                              Method
C
C
C     Formula    26.4.19   of Abramowitz  and     Stegun, Handbook  of
C     Mathematical Functions   (1966) is used   to reduce the chisqure
C     distribution to the incomplete distribution.
C
C     Computation of other parameters involve a search for a value that
C     produces  the desired  value  of P.   The search relies  on  the
C     monotinicity of P with the other parameter.
C
C**********************************************************************
C     .. Parameters ..
      DOUBLE PRECISION tol
      PARAMETER (tol=1.0D-8)
      DOUBLE PRECISION atol
      PARAMETER (atol=1.0D-50)
      DOUBLE PRECISION zero,inf
      PARAMETER (zero=1.0D-100,inf=1.0D100)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION bound,df,p,q,x
      INTEGER status,which
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ccum,cum,fx,porq,pq
      LOGICAL qhi,qleft,qporq
C     ..
C     .. External Functions ..
      DOUBLE PRECISION spmpar
      EXTERNAL spmpar
C     ..
C     .. External Subroutines ..
      EXTERNAL cumchi,dinvr,dstinv
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..
      IF (.NOT. ((which.LT.1).OR. (which.GT.3))) GO TO 30
      IF (.NOT. (which.LT.1)) GO TO 10
      bound = 1.0D0
      GO TO 20

   10 bound = 3.0D0
   20 status = -1
      RETURN

   30 IF (which.EQ.1) GO TO 70
      IF (.NOT. ((p.LT.0.0D0).OR. (p.GT.1.0D0))) GO TO 60
      IF (.NOT. (p.LT.0.0D0)) GO TO 40
      bound = 0.0D0
      GO TO 50

   40 bound = 1.0D0
   50 status = -2
      RETURN

   60 CONTINUE
   70 IF (which.EQ.1) GO TO 110
      IF (.NOT. ((q.LE.0.0D0).OR. (q.GT.1.0D0))) GO TO 100
      IF (.NOT. (q.LE.0.0D0)) GO TO 80
      bound = 0.0D0
      GO TO 90

   80 bound = 1.0D0
   90 status = -3
      RETURN

  100 CONTINUE
  110 IF (which.EQ.2) GO TO 130
      IF (.NOT. (x.LT.0.0D0)) GO TO 120
      bound = 0.0D0
      status = -4
      RETURN

  120 CONTINUE
  130 IF (which.EQ.3) GO TO 150
      IF (.NOT. (df.LE.0.0D0)) GO TO 140
      bound = 0.0D0
      status = -5
      RETURN

  140 CONTINUE
  150 IF (which.EQ.1) GO TO 190
      pq = p + q
      IF (.NOT. (abs(((pq)-0.5D0)-0.5D0).GT.
     +    (3.0D0*spmpar(1)))) GO TO 180
      IF (.NOT. (pq.LT.0.0D0)) GO TO 160
      bound = 0.0D0
      GO TO 170

  160 bound = 1.0D0
  170 status = 3
      RETURN

  180 CONTINUE
  190 IF (which.EQ.1) GO TO 220
      qporq = p .LE. q
      IF (.NOT. (qporq)) GO TO 200
      porq = p
      GO TO 210

  200 porq = q
  210 CONTINUE
  220 IF ((1).EQ. (which)) THEN
          status = 0
          CALL cumchi(x,df,p,q)
          IF (porq.GT.1.5D0) THEN
              status = 10
              RETURN

          END IF

      ELSE IF ((2).EQ. (which)) THEN
          x = 5.0D0
          CALL dstinv(0.0D0,inf,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,x,fx,qleft,qhi)
  230     IF (.NOT. (status.EQ.1)) GO TO 270
          CALL cumchi(x,df,cum,ccum)
          IF (.NOT. (qporq)) GO TO 240
          fx = cum - p
          GO TO 250

  240     fx = ccum - q
  250     IF (.NOT. ((fx+porq).GT.1.5D0)) GO TO 260
          status = 10
          RETURN

  260     CALL dinvr(status,x,fx,qleft,qhi)
          GO TO 230

  270     IF (.NOT. (status.EQ.-1)) GO TO 300
          IF (.NOT. (qleft)) GO TO 280
          status = 1
          bound = 0.0D0
          GO TO 290

  280     status = 2
          bound = inf
  290     CONTINUE
  300     CONTINUE

      ELSE IF ((3).EQ. (which)) THEN
          df = 5.0D0
          CALL dstinv(zero,inf,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,df,fx,qleft,qhi)
  310     IF (.NOT. (status.EQ.1)) GO TO 350
          CALL cumchi(x,df,cum,ccum)
          IF (.NOT. (qporq)) GO TO 320
          fx = cum - p
          GO TO 330

  320     fx = ccum - q
  330     IF (.NOT. ((fx+porq).GT.1.5D0)) GO TO 340
          status = 10
          RETURN

  340     CALL dinvr(status,df,fx,qleft,qhi)
          GO TO 310

  350     IF (.NOT. (status.EQ.-1)) GO TO 380
          IF (.NOT. (qleft)) GO TO 360
          status = 1
          bound = zero
          GO TO 370

  360     status = 2
          bound = inf
  370     CONTINUE
  380 END IF

      RETURN

      END
