      SUBROUTINE cdfgam(which,p,q,x,shape,scale,status,bound)
C**********************************************************************
C
C      SUBROUTINE CDFGAM( WHICH, P, Q, X, SHAPE, SCALE, STATUS, BOUND )
C               Cumulative Distribution Function
C                         GAMma Distribution
C
C
C                              Function
C
C
C     Calculates any one parameter of the gamma
C     distribution given values for the others.
C
C
C                              Arguments
C
C
C     WHICH --> Integer indicating which of the next four argument
C               values is to be calculated from the others.
C               Legal range: 1..4
C               iwhich = 1 : Calculate P and Q from X,SHAPE and SCALE
C               iwhich = 2 : Calculate X from P,Q,SHAPE and SCALE
C               iwhich = 3 : Calculate SHAPE from P,Q,X and SCALE
C               iwhich = 4 : Calculate SCALE from P,Q,X and SHAPE
C                    INTEGER WHICH
C
C     P <--> The integral from 0 to X of the gamma density.
C            Input range: [0,1].
C                    DOUBLE PRECISION P
C
C     Q <--> 1-P.
C            Input range: (0, 1].
C            P + Q = 1.0.
C                    DOUBLE PRECISION Q
C
C
C     X <--> The upper limit of integration of the gamma density.
C            Input range: [0, +infinity).
C            Search range: [0,1E100]
C                    DOUBLE PRECISION X
C
C     SHAPE <--> The shape parameter of the gamma density.
C                Input range: (0, +infinity).
C                Search range: [1E-100,1E100]
C                  DOUBLE PRECISION SHAPE
C
C
C     SCALE <--> The scale parameter of the gamma density.
C                Input range: (0, +infinity).
C                Search range: (1E-100,1E100]
C                   DOUBLE PRECISION SCALE
C
C     STATUS <-- 0 if calculation completed correctly
C               -I if input parameter number I is out of range
C                1 if answer appears to be lower than lowest
C                  search bound
C                2 if answer appears to be higher than greatest
C                  search bound
C                3 if P + Q .ne. 1
C                10 if the gamma or inverse gamma routine cannot
C                   compute the answer.  Usually happens only for
C                   X and SHAPE very large (gt 1E10 or more)
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
C     Cumulative distribution function (P) is calculated directly by
C     the code associated with:
C
C     DiDinato, A. R. and Morris, A. H. Computation of the  incomplete
C     gamma function  ratios  and their  inverse.   ACM  Trans.  Math.
C     Softw. 12 (1986), 377-393.
C
C     Computation of other parameters involve a seach for a value that
C     produces  the desired  value  of P.   The search relies  on  the
C     monotinicity of P with the other parameter.
C
C
C                              Note
C
C
C
C     The gamma density is proportional to
C       T**(SHAPE - 1) * EXP(- SCALE * T)
C
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
      DOUBLE PRECISION bound,p,q,scale,shape,x
      INTEGER status,which
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ccum,cum,fx,porq,pq,xscale,xx
      INTEGER ierr
      LOGICAL qhi,qleft,qporq
C     ..
C     .. External Functions ..
      DOUBLE PRECISION spmpar
      EXTERNAL spmpar
C     ..
C     .. External Subroutines ..
      EXTERNAL cumgam,dinvr,dstinv,gaminv
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

   40 bound = 1.0d0
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
      IF (.NOT. (shape.LE.0.0D0)) GO TO 140
      bound = 0.0D0
      status = -5
      RETURN

  140 CONTINUE
  150 IF (which.EQ.4) GO TO 170
      IF (.NOT. (scale.LE.0.0D0)) GO TO 160
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
  210 IF (which.EQ.1) GO TO 240
      qporq = p .LE. q
      IF (.NOT. (qporq)) GO TO 220
      porq = p
      GO TO 230

  220 porq = q
  230 CONTINUE
  240 IF ((1).EQ. (which)) THEN
          status = 0
          xscale = x*scale
          CALL cumgam(xscale,shape,p,q)
          IF (p.GT.1.5D0) status = 10

      ELSE IF ((2).EQ. (which)) THEN
          CALL gaminv(shape,xx,-1.0D0,p,q,ierr)
          IF (ierr.LT.0.0D0) THEN
              status = 10
              RETURN

          ELSE
              x = xx/scale
              status = 0
          END IF

      ELSE IF ((3).EQ. (which)) THEN
          shape = 5.0D0
          xscale = x*scale
          CALL dstinv(zero,inf,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,shape,fx,qleft,qhi)
  250     IF (.NOT. (status.EQ.1)) GO TO 290
          CALL cumgam(xscale,shape,cum,ccum)
          IF (.NOT. (qporq)) GO TO 260
          fx = cum - p
          GO TO 270

  260     fx = ccum - q
  270     IF (.NOT. ((qporq.AND. (cum.GT.1.5D0)).OR.
     +        ((.NOT.qporq).AND. (ccum.GT.1.5D0)))) GO TO 280
          status = 10
          RETURN

  280     CALL dinvr(status,shape,fx,qleft,qhi)
          GO TO 250

  290     IF (.NOT. (status.EQ.-1)) GO TO 320
          IF (.NOT. (qleft)) GO TO 300
          status = 1
          bound = zero
          GO TO 310

  300     status = 2
          bound = inf
  310     CONTINUE
  320     CONTINUE

      ELSE IF ((4).EQ. (which)) THEN
          CALL gaminv(shape,xx,-1.0D0,p,q,ierr)
          IF (ierr.LT.0.0D0) THEN
              status = 10
              RETURN

          ELSE
              scale = xx/x
              status = 0
          END IF

      END IF

      RETURN

      END
