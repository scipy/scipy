      SUBROUTINE cdfbin(which,p,q,s,xn,pr,ompr,status,bound)
C**********************************************************************
C
C      SUBROUTINE CDFBIN ( WHICH, P, Q, S, XN, PR, OMPR, STATUS, BOUND )
C               Cumulative Distribution Function
C                         BINomial distribution
C
C
C                              Function
C
C
C     Calculates any one parameter of the binomial
C     distribution given values for the others.
C
C
C                              Arguments
C
C
C     WHICH --> Integer indicating which of the next four argument
C               values is to be calculated from the others.
C               Legal range: 1..4
C               iwhich = 1 : Calculate P and Q from S,XN,PR and OMPR
C               iwhich = 2 : Calculate S from P,Q,XN,PR and OMPR
C               iwhich = 3 : Calculate XN from P,Q,S,PR and OMPR
C               iwhich = 4 : Calculate PR and OMPR from P,Q,S and XN
C                    INTEGER WHICH
C
C     P <--> The cumulation from 0 to S of the binomial distribution.
C            (Probablility of S or fewer successes in XN trials each
C            with probability of success PR.)
C            Input range: [0,1].
C                    DOUBLE PRECISION P
C
C     Q <--> 1-P.
C            Input range: [0, 1].
C            P + Q = 1.0.
C                    DOUBLE PRECISION Q
C
C     S <--> The number of successes observed.
C            Input range: [0, XN]
C            Search range: [0, XN]
C                    DOUBLE PRECISION S
C
C     XN  <--> The number of binomial trials.
C              Input range: (0, +infinity).
C              Search range: [1E-100, 1E100]
C                    DOUBLE PRECISION XN
C
C     PR  <--> The probability of success in each binomial trial.
C              Input range: [0,1].
C              Search range: [0,1]
C                    DOUBLE PRECISION PR
C
C     OMPR  <--> 1-PR
C              Input range: [0,1].
C              Search range: [0,1]
C              PR + OMPR = 1.0
C                    DOUBLE PRECISION OMPR
C
C     STATUS <-- 0 if calculation completed correctly
C               -I if input parameter number I is out of range
C                1 if answer appears to be lower than lowest
C                  search bound
C                2 if answer appears to be higher than greatest
C                  search bound
C                3 if P + Q .ne. 1
C                4 if PR + OMPR .ne. 1
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
C     Formula  26.5.24    of   Abramowitz  and    Stegun,  Handbook   of
C     Mathematical   Functions (1966) is   used  to reduce the  binomial
C     distribution  to  the  cumulative incomplete    beta distribution.
C
C     Computation of other parameters involve a search for a value that
C     produces  the desired  value  of P.   The search relies  on  the
C     monotinicity of P with the other parameter.
C
C
C**********************************************************************
C     .. Parameters ..
      DOUBLE PRECISION atol
      PARAMETER (atol=1.0D-50)
      DOUBLE PRECISION tol
      PARAMETER (tol=1.0D-8)
      DOUBLE PRECISION zero,inf
      PARAMETER (zero=1.0D-100,inf=1.0D100)
      DOUBLE PRECISION one
      PARAMETER (one=1.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION bound,ompr,p,pr,q,s,xn
      INTEGER status,which
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ccum,cum,fx,pq,prompr,xhi,xlo
      LOGICAL qhi,qleft,qporq
C     ..
C     .. External Functions ..
      DOUBLE PRECISION spmpar
      EXTERNAL spmpar
C     ..
C     .. External Subroutines ..
      EXTERNAL cumbin,dinvr,dstinv,dstzr,dzror
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..
      IF (.NOT. ((which.LT.1).AND. (which.GT.4))) GO TO 30
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
  110 IF (which.EQ.3) GO TO 130
      IF (.NOT. (xn.LE.0.0D0)) GO TO 120
      bound = 0.0D0
      status = -5
      RETURN

  120 CONTINUE
  130 IF (which.EQ.2) GO TO 170
      IF (.NOT. ((s.LT.0.0D0).OR. ((which.NE.3).AND.
     +    (s.GT.xn)))) GO TO 160
      IF (.NOT. (s.LT.0.0D0)) GO TO 140
      bound = 0.0D0
      GO TO 150

  140 bound = xn
  150 status = -4
      RETURN

  160 CONTINUE
  170 IF (which.EQ.4) GO TO 210
      IF (.NOT. ((pr.LT.0.0D0).OR. (pr.GT.1.0D0))) GO TO 200
      IF (.NOT. (pr.LT.0.0D0)) GO TO 180
      bound = 0.0D0
      GO TO 190

  180 bound = 1.0D0
  190 status = -6
      RETURN

  200 CONTINUE
  210 IF (which.EQ.4) GO TO 250
      IF (.NOT. ((ompr.LT.0.0D0).OR. (ompr.GT.1.0D0))) GO TO 240
      IF (.NOT. (ompr.LT.0.0D0)) GO TO 220
      bound = 0.0D0
      GO TO 230

  220 bound = 1.0D0
  230 status = -7
      RETURN

  240 CONTINUE
  250 IF (which.EQ.1) GO TO 290
      pq = p + q
      IF (.NOT. (abs(((pq)-0.5D0)-0.5D0).GT.
     +    (3.0D0*spmpar(1)))) GO TO 280
      IF (.NOT. (pq.LT.0.0D0)) GO TO 260
      bound = 0.0D0
      GO TO 270

  260 bound = 1.0D0
  270 status = 3
      RETURN

  280 CONTINUE
  290 IF (which.EQ.4) GO TO 330
      prompr = pr + ompr
      IF (.NOT. (abs(((prompr)-0.5D0)-0.5D0).GT.
     +    (3.0D0*spmpar(1)))) GO TO 320
      IF (.NOT. (prompr.LT.0.0D0)) GO TO 300
      bound = 0.0D0
      GO TO 310

  300 bound = 1.0D0
  310 status = 4
      RETURN

  320 CONTINUE
  330 IF (.NOT. (which.EQ.1)) qporq = p .LE. q
      IF ((1).EQ. (which)) THEN
          CALL cumbin(s,xn,pr,ompr,p,q)
          status = 0

      ELSE IF ((2).EQ. (which)) THEN
          s = xn/2.0D0
          CALL dstinv(0.0D0,xn,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,s,fx,qleft,qhi)
  340     IF (.NOT. (status.EQ.1)) GO TO 370
          CALL cumbin(s,xn,pr,ompr,cum,ccum)
          IF (.NOT. (qporq)) GO TO 350
          fx = cum - p
          GO TO 360

  350     fx = ccum - q
  360     CALL dinvr(status,s,fx,qleft,qhi)
          GO TO 340

  370     IF (.NOT. (status.EQ.-1)) GO TO 400
          IF (.NOT. (qleft)) GO TO 380
          status = 1
          bound = 0.0D0
          GO TO 390

  380     status = 2
          bound = xn
  390     CONTINUE
  400     CONTINUE

      ELSE IF ((3).EQ. (which)) THEN
          xn = 5.0D0
          CALL dstinv(zero,inf,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,xn,fx,qleft,qhi)
  410     IF (.NOT. (status.EQ.1)) GO TO 440
          CALL cumbin(s,xn,pr,ompr,cum,ccum)
          IF (.NOT. (qporq)) GO TO 420
          fx = cum - p
          GO TO 430

  420     fx = ccum - q
  430     CALL dinvr(status,xn,fx,qleft,qhi)
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
          CALL dstzr(0.0D0,1.0D0,atol,tol)
          IF (.NOT. (qporq)) GO TO 500
          status = 0
          CALL dzror(status,pr,fx,xlo,xhi,qleft,qhi)
          ompr = one - pr
  480     IF (.NOT. (status.EQ.1)) GO TO 490
          CALL cumbin(s,xn,pr,ompr,cum,ccum)
          fx = cum - p
          CALL dzror(status,pr,fx,xlo,xhi,qleft,qhi)
          ompr = one - pr
          GO TO 480

  490     GO TO 530

  500     status = 0
          CALL dzror(status,ompr,fx,xlo,xhi,qleft,qhi)
          pr = one - ompr
  510     IF (.NOT. (status.EQ.1)) GO TO 520
          CALL cumbin(s,xn,pr,ompr,cum,ccum)
          fx = ccum - q
          CALL dzror(status,ompr,fx,xlo,xhi,qleft,qhi)
          pr = one - ompr
          GO TO 510

  520     CONTINUE
  530     IF (.NOT. (status.EQ.-1)) GO TO 560
          IF (.NOT. (qleft)) GO TO 540
          status = 1
          bound = 0.0D0
          GO TO 550

  540     status = 2
          bound = 1.0D0
  550     CONTINUE
  560 END IF

      RETURN

      END
