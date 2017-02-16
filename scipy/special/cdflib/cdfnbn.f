      SUBROUTINE cdfnbn(which,p,q,s,xn,pr,ompr,status,bound)
C**********************************************************************
C
C      SUBROUTINE CDFNBN ( WHICH, P, S, XN, PR, STATUS, BOUND )
C               Cumulative Distribution Function
C               Negative BiNomial distribution
C
C
C                              Function
C
C
C     Calculates any one parameter of the negative binomial
C     distribution given values for the others.
C
C     The  cumulative  negative   binomial  distribution  returns  the
C     probability that there  will be  F or fewer failures before  the
C     XNth success in binomial trials each of which has probability of
C     success PR.
C
C     The individual term of the negative binomial is the probability of
C     S failures before XN successes and is
C          Choose( S, XN+S-1 ) * PR^(XN) * (1-PR)^S
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
C     P <--> The cumulation from 0 to S of the  negative
C            binomial distribution.
C            Input range: [0,1].
C                    DOUBLE PRECISION P
C
C     Q <--> 1-P.
C            Input range: (0, 1].
C            P + Q = 1.0.
C                    DOUBLE PRECISION Q
C
C     S <--> The upper limit of cumulation of the binomial distribution.
C            There are F or fewer failures before the XNth success.
C            Input range: [0, +infinity).
C            Search range: [0, 1E100]
C                    DOUBLE PRECISION S
C
C     XN  <--> The number of successes.
C              Input range: [0, +infinity).
C              Search range: [0, 1E100]
C                    DOUBLE PRECISION XN
C
C     PR  <--> The probability of success in each binomial trial.
C              Input range: [0,1].
C              Search range: [0,1].
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
C     Formula   26.5.26   of   Abramowitz  and  Stegun,  Handbook   of
C     Mathematical Functions (1966) is used  to  reduce calculation of
C     the cumulative distribution  function to that of  an  incomplete
C     beta.
C
C     Computation of other parameters involve a seach for a value that
C     produces  the desired  value  of P.   The search relies  on  the
C     monotinicity of P with the other parameter.
C
C
C**********************************************************************
C     .. Parameters ..
      DOUBLE PRECISION tol
      PARAMETER (tol=1.0D-8)
      DOUBLE PRECISION atol
      PARAMETER (atol=1.0D-50)
      DOUBLE PRECISION inf
      PARAMETER (inf=1.0D100)
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
      EXTERNAL cumnbn,dinvr,dstinv,dstzr,dzror
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
      IF (.NOT. (s.LT.0.0D0)) GO TO 120
      bound = 0.0D0
      status = -4
      RETURN

  120 CONTINUE
  130 IF (which.EQ.3) GO TO 150
      IF (.NOT. (xn.LT.0.0D0)) GO TO 140
      bound = 0.0D0
      status = -5
      RETURN

  140 CONTINUE
  150 IF (which.EQ.4) GO TO 190
      IF (.NOT. ((pr.LT.0.0D0).OR. (pr.GT.1.0D0))) GO TO 180
      IF (.NOT. (pr.LT.0.0D0)) GO TO 160
      bound = 0.0D0
      GO TO 170

  160 bound = 1.0D0
  170 status = -6
      RETURN

  180 CONTINUE
  190 IF (which.EQ.4) GO TO 230
      IF (.NOT. ((ompr.LT.0.0D0).OR. (ompr.GT.1.0D0))) GO TO 220
      IF (.NOT. (ompr.LT.0.0D0)) GO TO 200
      bound = 0.0D0
      GO TO 210

  200 bound = 1.0D0
  210 status = -7
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
  270 IF (which.EQ.4) GO TO 310
      prompr = pr + ompr
      IF (.NOT. (abs(((prompr)-0.5D0)-0.5D0).GT.
     +    (3.0D0*spmpar(1)))) GO TO 300
      IF (.NOT. (prompr.LT.0.0D0)) GO TO 280
      bound = 0.0D0
      GO TO 290

  280 bound = 1.0D0
  290 status = 4
      RETURN

  300 CONTINUE
  310 IF (.NOT. (which.EQ.1)) qporq = p .LE. q
      IF ((1).EQ. (which)) THEN
          CALL cumnbn(s,xn,pr,ompr,p,q)
          status = 0

      ELSE IF ((2).EQ. (which)) THEN
          s = 5.0D0
          CALL dstinv(0.0D0,inf,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,s,fx,qleft,qhi)
  320     IF (.NOT. (status.EQ.1)) GO TO 350
          CALL cumnbn(s,xn,pr,ompr,cum,ccum)
          IF (.NOT. (qporq)) GO TO 330
          fx = cum - p
          GO TO 340

  330     fx = ccum - q
  340     CALL dinvr(status,s,fx,qleft,qhi)
          GO TO 320

  350     IF (.NOT. (status.EQ.-1)) GO TO 380
          IF (.NOT. (qleft)) GO TO 360
          status = 1
          bound = 0.0D0
          GO TO 370

  360     status = 2
          bound = inf
  370     CONTINUE
  380     CONTINUE

      ELSE IF ((3).EQ. (which)) THEN
          xn = 5.0D0
          CALL dstinv(0.0D0,inf,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,xn,fx,qleft,qhi)
  390     IF (.NOT. (status.EQ.1)) GO TO 420
          CALL cumnbn(s,xn,pr,ompr,cum,ccum)
          IF (.NOT. (qporq)) GO TO 400
          fx = cum - p
          GO TO 410

  400     fx = ccum - q
  410     CALL dinvr(status,xn,fx,qleft,qhi)
          GO TO 390

  420     IF (.NOT. (status.EQ.-1)) GO TO 450
          IF (.NOT. (qleft)) GO TO 430
          status = 1
          bound = 0.0D0
          GO TO 440

  430     status = 2
          bound = inf
  440     CONTINUE
  450     CONTINUE

      ELSE IF ((4).EQ. (which)) THEN
          CALL dstzr(0.0D0,1.0D0,atol,tol)
          IF (.NOT. (qporq)) GO TO 480
          status = 0
          CALL dzror(status,pr,fx,xlo,xhi,qleft,qhi)
          ompr = one - pr
  460     IF (.NOT. (status.EQ.1)) GO TO 470
          CALL cumnbn(s,xn,pr,ompr,cum,ccum)
          fx = cum - p
          CALL dzror(status,pr,fx,xlo,xhi,qleft,qhi)
          ompr = one - pr
          GO TO 460

  470     GO TO 510

  480     status = 0
          CALL dzror(status,ompr,fx,xlo,xhi,qleft,qhi)
          pr = one - ompr
  490     IF (.NOT. (status.EQ.1)) GO TO 500
          CALL cumnbn(s,xn,pr,ompr,cum,ccum)
          fx = ccum - q
          CALL dzror(status,ompr,fx,xlo,xhi,qleft,qhi)
          pr = one - ompr
          GO TO 490

  500     CONTINUE
  510     IF (.NOT. (status.EQ.-1)) GO TO 540
          IF (.NOT. (qleft)) GO TO 520
          status = 1
          bound = 0.0D0
          GO TO 530

  520     status = 2
          bound = 1.0D0
  530     CONTINUE
  540 END IF

      RETURN

      END
