      SUBROUTINE cdft(which,p,q,t,df,status,bound)
C**********************************************************************
C
C      SUBROUTINE CDFT( WHICH, P, Q, T, DF, STATUS, BOUND )
C               Cumulative Distribution Function
C                         T distribution
C
C
C                              Function
C
C
C     Calculates any one parameter of the t distribution given
C     values for the others.
C
C
C                              Arguments
C
C
C     WHICH --> Integer indicating which  argument
C               values is to be calculated from the others.
C               Legal range: 1..3
C               iwhich = 1 : Calculate P and Q from T and DF
C               iwhich = 2 : Calculate T from P,Q and DF
C               iwhich = 3 : Calculate DF from P,Q and T
C                    INTEGER WHICH
C
C        P <--> The integral from -infinity to t of the t-density.
C              Input range: (0,1].
C                    DOUBLE PRECISION P
C
C     Q <--> 1-P.
C            Input range: (0, 1].
C            P + Q = 1.0.
C                    DOUBLE PRECISION Q
C
C        T <--> Upper limit of integration of the t-density.
C               Input range: ( -infinity, +infinity).
C               Search range: [ -1E100, 1E100 ]
C                    DOUBLE PRECISION T
C
C        DF <--> Degrees of freedom of the t-distribution.
C                Input range: (0 , +infinity).
C                Search range: [1e-100, 1E10]
C                    DOUBLE PRECISION DF
C
C     STATUS <-- 0 if calculation completed correctly
C               -I if input parameter number I is out of range
C                1 if answer appears to be lower than lowest
C                  search bound
C                2 if answer appears to be higher than greatest
C                  search bound
C                3 if P + Q .ne. 1
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
C     Formula  26.5.27  of   Abramowitz   and  Stegun,   Handbook   of
C     Mathematical Functions  (1966) is used to reduce the computation
C     of the cumulative distribution function to that of an incomplete
C     beta.
C
C     Computation of other parameters involve a seach for a value that
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
      DOUBLE PRECISION rtinf
      PARAMETER (rtinf=1.0D100)
      DOUBLE PRECISION maxdf
      PARAMETER (maxdf=1.0d10)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION bound,df,p,q,t
      INTEGER status,which
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ccum,cum,fx,pq
      LOGICAL qhi,qleft,qporq
C     ..
C     .. External Functions ..
      DOUBLE PRECISION dt1,spmpar
      EXTERNAL dt1,spmpar
C     ..
C     .. External Subroutines ..
      EXTERNAL cumt,dinvr,dstinv
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
      IF (.NOT. ((p.LE.0.0D0).OR. (p.GT.1.0D0))) GO TO 60
      IF (.NOT. (p.LE.0.0D0)) GO TO 40
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
  110 IF (which.EQ.3) GO TO 130
      IF (.NOT. (df.LE.0.0D0)) GO TO 120
      bound = 0.0D0
      status = -5
      RETURN

  120 CONTINUE
  130 IF (which.EQ.1) GO TO 170
      pq = p + q
      IF (.NOT. (abs(((pq)-0.5D0)-0.5D0).GT.
     +    (3.0D0*spmpar(1)))) GO TO 160
      IF (.NOT. (pq.LT.0.0D0)) GO TO 140
      bound = 0.0D0
      GO TO 150

  140 bound = 1.0D0
  150 status = 3
      RETURN

  160 CONTINUE
  170 IF (.NOT. (which.EQ.1)) qporq = p .LE. q
      IF ((1).EQ. (which)) THEN
          CALL cumt(t,df,p,q)
          status = 0

      ELSE IF ((2).EQ. (which)) THEN
          t = dt1(p,q,df)
          CALL dstinv(-rtinf,rtinf,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,t,fx,qleft,qhi)
  180     IF (.NOT. (status.EQ.1)) GO TO 210
          CALL cumt(t,df,cum,ccum)
          IF (.NOT. (qporq)) GO TO 190
          fx = cum - p
          GO TO 200

  190     fx = ccum - q
  200     CALL dinvr(status,t,fx,qleft,qhi)
          GO TO 180

  210     IF (.NOT. (status.EQ.-1)) GO TO 240
          IF (.NOT. (qleft)) GO TO 220
          status = 1
          bound = -rtinf
          GO TO 230

  220     status = 2
          bound = rtinf
  230     CONTINUE
  240     CONTINUE

      ELSE IF ((3).EQ. (which)) THEN
          df = 5.0D0
          CALL dstinv(zero,maxdf,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,df,fx,qleft,qhi)
  250     IF (.NOT. (status.EQ.1)) GO TO 280
          CALL cumt(t,df,cum,ccum)
          IF (.NOT. (qporq)) GO TO 260
          fx = cum - p
          GO TO 270

  260     fx = ccum - q
  270     CALL dinvr(status,df,fx,qleft,qhi)
          GO TO 250

  280     IF (.NOT. (status.EQ.-1)) GO TO 310
          IF (.NOT. (qleft)) GO TO 290
          status = 1
          bound = zero
          GO TO 300

  290     status = 2
          bound = maxdf
  300     CONTINUE
  310 END IF

      RETURN

      END
