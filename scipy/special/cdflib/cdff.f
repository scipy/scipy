      SUBROUTINE cdff(which,p,q,f,dfn,dfd,status,bound)
C**********************************************************************
C
C      SUBROUTINE CDFF( WHICH, P, Q, F, DFN, DFD, STATUS, BOUND )
C               Cumulative Distribution Function
C               F distribution
C
C
C                              Function
C
C
C     Calculates any one parameter of the F distribution
C     given values for the others.
C
C
C                              Arguments
C
C
C     WHICH --> Integer indicating which of the next four argument
C               values is to be calculated from the others.
C               Legal range: 1..4
C               iwhich = 1 : Calculate P and Q from F,DFN and DFD
C               iwhich = 2 : Calculate F from P,Q,DFN and DFD
C               iwhich = 3 : Calculate DFN from P,Q,F and DFD
C               iwhich = 4 : Calculate DFD from P,Q,F and DFN
C                    INTEGER WHICH
C
C       P <--> The integral from 0 to F of the f-density.
C              Input range: [0,1].
C                    DOUBLE PRECISION P
C
C       Q <--> 1-P.
C              Input range: (0, 1].
C              P + Q = 1.0.
C                    DOUBLE PRECISION Q
C
C       F <--> Upper limit of integration of the f-density.
C              Input range: [0, +infinity).
C              Search range: [0,1E100]
C                    DOUBLE PRECISION F
C
C     DFN < --> Degrees of freedom of the numerator sum of squares.
C               Input range: (0, +infinity).
C               Search range: [ 1E-100, 1E100]
C                    DOUBLE PRECISION DFN
C
C     DFD < --> Degrees of freedom of the denominator sum of squares.
C               Input range: (0, +infinity).
C               Search range: [ 1E-100, 1E100]
C                    DOUBLE PRECISION DFD
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
C     Formula   26.6.2   of   Abramowitz   and   Stegun,  Handbook  of
C     Mathematical  Functions (1966) is used to reduce the computation
C     of the  cumulative  distribution function for the  F  variate to
C     that of an incomplete beta.
C
C     Computation of other parameters involve a search for a value that
C     produces  the desired  value  of P.   The search relies  on  the
C     monotinicity of P with the other parameter.
C
C                              WARNING
C
C     The value of the  cumulative  F distribution is  not necessarily
C     monotone in  either degrees of freedom.  There  thus may  be two
C     values  that  provide a given CDF  value.   This routine assumes
C     monotonicity and will find an arbitrary one of the two values.
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
      DOUBLE PRECISION bound,dfd,dfn,f,p,q
      INTEGER status,which
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ccum,cum,fx,pq
      LOGICAL qhi,qleft,qporq
C     ..
C     .. External Functions ..
      DOUBLE PRECISION spmpar
      EXTERNAL spmpar
C     ..
C     .. External Subroutines ..
      EXTERNAL cumf,dinvr,dstinv
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..
      IF (.NOT. ((which.LT.1).OR. (which.GT.4))) GO TO 30
      IF (.NOT. (which.LT.1)) GO TO 10
      bound = 1.0D0
      GO TO 20

   10 bound = 4.0D0
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
      IF (.NOT. (f.LT.0.0D0)) GO TO 120
      bound = 0.0D0
      status = -4
      RETURN

  120 CONTINUE
  130 IF (which.EQ.3) GO TO 150
      IF (.NOT. (dfn.LE.0.0D0)) GO TO 140
      bound = 0.0D0
      status = -5
      RETURN

  140 CONTINUE
  150 IF (which.EQ.4) GO TO 170
      IF (.NOT. (dfd.LE.0.0D0)) GO TO 160
      bound = 0.0D0
      status = -6
      RETURN

  160 CONTINUE
  170 IF (which.EQ.1) GO TO 210
      pq = p + q
      IF (.NOT. (abs(((pq)-0.5D0)-0.5D0).GT.
     +    (3.0D0*spmpar(1)))) GO TO 200
      IF (.NOT. (pq.LT.0.0D0)) GO TO 180
      bound = 0.0D0
      GO TO 190

  180 bound = 1.0D0
  190 status = 3
      RETURN

  200 CONTINUE
  210 IF (.NOT. (which.EQ.1)) qporq = p .LE. q
      IF ((1).EQ. (which)) THEN
          CALL cumf(f,dfn,dfd,p,q)
          status = 0

      ELSE IF ((2).EQ. (which)) THEN
          f = 5.0D0
          CALL dstinv(0.0D0,inf,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,f,fx,qleft,qhi)
  220     IF (.NOT. (status.EQ.1)) GO TO 250
          CALL cumf(f,dfn,dfd,cum,ccum)
          IF (.NOT. (qporq)) GO TO 230
          fx = cum - p
          GO TO 240

  230     fx = ccum - q
  240     CALL dinvr(status,f,fx,qleft,qhi)
          GO TO 220

  250     IF (.NOT. (status.EQ.-1)) GO TO 280
          IF (.NOT. (qleft)) GO TO 260
          status = 1
          bound = 0.0D0
          GO TO 270

  260     status = 2
          bound = inf
  270     CONTINUE
  280     CONTINUE

      ELSE IF ((3).EQ. (which)) THEN
          dfn = 5.0D0
          CALL dstinv(zero,inf,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,dfn,fx,qleft,qhi)
  290     IF (.NOT. (status.EQ.1)) GO TO 320
          CALL cumf(f,dfn,dfd,cum,ccum)
          IF (.NOT. (qporq)) GO TO 300
          fx = cum - p
          GO TO 310

  300     fx = ccum - q
  310     CALL dinvr(status,dfn,fx,qleft,qhi)
          GO TO 290

  320     IF (.NOT. (status.EQ.-1)) GO TO 350
          IF (.NOT. (qleft)) GO TO 330
          status = 1
          bound = zero
          GO TO 340

  330     status = 2
          bound = inf
  340     CONTINUE
  350     CONTINUE

      ELSE IF ((4).EQ. (which)) THEN
          dfd = 5.0D0
          CALL dstinv(zero,inf,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,dfd,fx,qleft,qhi)
  360     IF (.NOT. (status.EQ.1)) GO TO 390
          CALL cumf(f,dfn,dfd,cum,ccum)
          IF (.NOT. (qporq)) GO TO 370
          fx = cum - p
          GO TO 380

  370     fx = ccum - q
  380     CALL dinvr(status,dfd,fx,qleft,qhi)
          GO TO 360

  390     IF (.NOT. (status.EQ.-1)) GO TO 420
          IF (.NOT. (qleft)) GO TO 400
          status = 1
          bound = zero
          GO TO 410

  400     status = 2
          bound = inf
  410     CONTINUE
  420 END IF

      RETURN

      END
