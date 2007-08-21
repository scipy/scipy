      SUBROUTINE cdfbet(which,p,q,x,y,a,b,status,bound)
C**********************************************************************
C
C      SUBROUTINE CDFBET( WHICH, P, Q, X, Y, A, B, STATUS, BOUND )
C               Cumulative Distribution Function
C                         BETa Distribution
C
C
C                              Function
C
C
C     Calculates any one parameter of the beta distribution given
C     values for the others.
C
C
C                              Arguments
C
C
C     WHICH --> Integer indicating which of the next four argument
C               values is to be calculated from the others.
C               Legal range: 1..4
C               iwhich = 1 : Calculate P and Q from X,Y,A and B
C               iwhich = 2 : Calculate X and Y from P,Q,A and B
C               iwhich = 3 : Calculate A from P,Q,X,Y and B
C               iwhich = 4 : Calculate B from P,Q,X,Y and A
C
C                    INTEGER WHICH
C
C     P <--> The integral from 0 to X of the chi-square
C            distribution.
C            Input range: [0, 1].
C                    DOUBLE PRECISION P
C
C     Q <--> 1-P.
C            Input range: [0, 1].
C            P + Q = 1.0.
C                    DOUBLE PRECISION Q
C
C     X <--> Upper limit of integration of beta density.
C            Input range: [0,1].
C            Search range: [0,1]
C                    DOUBLE PRECISION X
C
C     Y <--> 1-X.
C            Input range: [0,1].
C            Search range: [0,1]
C            X + Y = 1.0.
C                    DOUBLE PRECISION Y
C
C     A <--> The first parameter of the beta density.
C            Input range: (0, +infinity).
C            Search range: [1D-100,1D100]
C                    DOUBLE PRECISION A
C
C     B <--> The second parameter of the beta density.
C            Input range: (0, +infinity).
C            Search range: [1D-100,1D100]
C                    DOUBLE PRECISION B
C
C     STATUS <-- 0 if calculation completed correctly
C               -I if input parameter number I is out of range
C                1 if answer appears to be lower than lowest
C                  search bound
C                2 if answer appears to be higher than greatest
C                  search bound
C                3 if P + Q .ne. 1
C                4 if X + Y .ne. 1
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
C     Cumulative distribution function  (P)  is calculated directly by
C     code associated with the following reference.
C
C     DiDinato, A. R. and Morris,  A.   H.  Algorithm 708: Significant
C     Digit Computation of the Incomplete  Beta  Function Ratios.  ACM
C     Trans. Math.  Softw. 18 (1993), 360-373.
C
C     Computation of other parameters involve a seach for a value that
C     produces  the desired  value  of P.   The search relies  on  the
C     monotinicity of P with the other parameter.
C
C
C                              Note
C
C
C     The beta density is proportional to
C               t^(A-1) * (1-t)^(B-1)
C
C**********************************************************************
C     .. Parameters ..
      DOUBLE PRECISION tol
      PARAMETER (tol=1.0D-8)
      DOUBLE PRECISION atol
      PARAMETER (atol=1.0D-50)
      DOUBLE PRECISION zero,inf
      PARAMETER (zero=1.0D-100,inf=1.0D100)
      DOUBLE PRECISION one
      PARAMETER (one=1.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,bound,p,q,x,y
      INTEGER status,which
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ccum,cum,fx,pq,xhi,xlo,xy
      LOGICAL qhi,qleft,qporq
C     ..
C     .. External Functions ..
      DOUBLE PRECISION spmpar
      EXTERNAL spmpar
C     ..
C     .. External Subroutines ..
      EXTERNAL cumbet,dinvr,dstinv,dstzr,dzror
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
      IF (.NOT. ((q.LT.0.0D0).OR. (q.GT.1.0D0))) GO TO 100
      IF (.NOT. (q.LT.0.0D0)) GO TO 80
      bound = 0.0D0
      GO TO 90

   80 bound = 1.0D0
   90 status = -3
      RETURN

  100 CONTINUE
  110 IF (which.EQ.2) GO TO 150
      IF (.NOT. ((x.LT.0.0D0).OR. (x.GT.1.0D0))) GO TO 140
      IF (.NOT. (x.LT.0.0D0)) GO TO 120
      bound = 0.0D0
      GO TO 130

  120 bound = 1.0D0
  130 status = -4
      RETURN

  140 CONTINUE
  150 IF (which.EQ.2) GO TO 190
      IF (.NOT. ((y.LT.0.0D0).OR. (y.GT.1.0D0))) GO TO 180
      IF (.NOT. (y.LT.0.0D0)) GO TO 160
      bound = 0.0D0
      GO TO 170

  160 bound = 1.0D0
  170 status = -5
      RETURN

  180 CONTINUE
  190 IF (which.EQ.3) GO TO 210
      IF (.NOT. (a.LE.0.0D0)) GO TO 200
      bound = 0.0D0
      status = -6
      RETURN

  200 CONTINUE
  210 IF (which.EQ.4) GO TO 230
      IF (.NOT. (b.LE.0.0D0)) GO TO 220
      bound = 0.0D0
      status = -7
      RETURN

  220 CONTINUE
  230 IF (which.EQ.1) GO TO 270
      pq = p + q
      IF (.NOT. (abs(((pq)-0.5D0)-0.5D0).GT.
     +    (3.0D0*spmpar(1)))) GO TO 260
      IF (.NOT. (pq.LT.0.0D0)) GO TO 240
      bound = 0.0D0
      GO TO 250

  240 bound = 1.0D0
  250 status = 3
      RETURN

  260 CONTINUE
  270 IF (which.EQ.2) GO TO 310
      xy = x + y
      IF (.NOT. (abs(((xy)-0.5D0)-0.5D0).GT.
     +    (3.0D0*spmpar(1)))) GO TO 300
      IF (.NOT. (xy.LT.0.0D0)) GO TO 280
      bound = 0.0D0
      GO TO 290

  280 bound = 1.0D0
  290 status = 4
      RETURN

  300 CONTINUE
  310 IF (.NOT. (which.EQ.1)) qporq = p .LE. q
      IF ((1).EQ. (which)) THEN
          CALL cumbet(x,y,a,b,p,q)
          status = 0

      ELSE IF ((2).EQ. (which)) THEN
          CALL dstzr(0.0D0,1.0D0,atol,tol)
          IF (.NOT. (qporq)) GO TO 340
          status = 0
          CALL dzror(status,x,fx,xlo,xhi,qleft,qhi)
          y = one - x
  320     IF (.NOT. (status.EQ.1)) GO TO 330
          CALL cumbet(x,y,a,b,cum,ccum)
          fx = cum - p
          CALL dzror(status,x,fx,xlo,xhi,qleft,qhi)
          y = one - x
          GO TO 320

  330     GO TO 370

  340     status = 0
          CALL dzror(status,y,fx,xlo,xhi,qleft,qhi)
          x = one - y
  350     IF (.NOT. (status.EQ.1)) GO TO 360
          CALL cumbet(x,y,a,b,cum,ccum)
          fx = ccum - q
          CALL dzror(status,y,fx,xlo,xhi,qleft,qhi)
          x = one - y
          GO TO 350

  360     CONTINUE
  370     IF (.NOT. (status.EQ.-1)) GO TO 400
          IF (.NOT. (qleft)) GO TO 380
          status = 1
          bound = 0.0D0
          GO TO 390

  380     status = 2
          bound = 1.0D0
  390     CONTINUE
  400     CONTINUE

      ELSE IF ((3).EQ. (which)) THEN
          a = 5.0D0
          CALL dstinv(zero,inf,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,a,fx,qleft,qhi)
  410     IF (.NOT. (status.EQ.1)) GO TO 440
          CALL cumbet(x,y,a,b,cum,ccum)
          IF (.NOT. (qporq)) GO TO 420
          fx = cum - p
          GO TO 430

  420     fx = ccum - q
  430     CALL dinvr(status,a,fx,qleft,qhi)
          GO TO 410

  440     IF (.NOT. (status.EQ.-1)) GO TO 470
          IF (.NOT. (qleft)) GO TO 450
          status = 1
          bound = zero
          GO TO 460

  450     status = 2
          bound = inf
  460     CONTINUE
  470     CONTINUE

      ELSE IF ((4).EQ. (which)) THEN
          b = 5.0D0
          CALL dstinv(zero,inf,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,b,fx,qleft,qhi)
  480     IF (.NOT. (status.EQ.1)) GO TO 510
          CALL cumbet(x,y,a,b,cum,ccum)
          IF (.NOT. (qporq)) GO TO 490
          fx = cum - p
          GO TO 500

  490     fx = ccum - q
  500     CALL dinvr(status,b,fx,qleft,qhi)
          GO TO 480

  510     IF (.NOT. (status.EQ.-1)) GO TO 540
          IF (.NOT. (qleft)) GO TO 520
          status = 1
          bound = zero
          GO TO 530

  520     status = 2
          bound = inf
  530     CONTINUE
  540 END IF

      RETURN

      END
