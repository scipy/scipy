      SUBROUTINE cdfchn(which,p,q,x,df,pnonc,status,bound)
C**********************************************************************
C
C      SUBROUTINE CDFCHN( WHICH, P, Q, X, DF, PNONC, STATUS, BOUND )
C               Cumulative Distribution Function
C               Non-central Chi-Square
C
C
C                              Function
C
C
C     Calculates any one parameter of the non-central chi-square
C     distribution given values for the others.
C
C
C                              Arguments
C
C
C     WHICH --> Integer indicating which of the next three argument
C               values is to be calculated from the others.
C               Input range: 1..4
C               iwhich = 1 : Calculate P and Q from X and DF
C               iwhich = 2 : Calculate X from P,DF and PNONC
C               iwhich = 3 : Calculate DF from P,X and PNONC
C               iwhich = 3 : Calculate PNONC from P,X and DF
C                    INTEGER WHICH
C
C     P <--> The integral from 0 to X of the non-central chi-square
C            distribution.
C            Input range: [0, 1-1E-16).
C                    DOUBLE PRECISION P
C
C     Q <--> 1-P.
C            Q is not used by this subroutine and is only included
C            for similarity with other cdf* routines.
C                    DOUBLE PRECISION Q
C
C     X <--> Upper limit of integration of the non-central
C            chi-square distribution.
C            Input range: [0, +infinity).
C            Search range: [0,1E300]
C                    DOUBLE PRECISION X
C
C     DF <--> Degrees of freedom of the non-central
C             chi-square distribution.
C             Input range: (0, +infinity).
C             Search range: [ 1E-300, 1E300]
C                    DOUBLE PRECISION DF
C
C     PNONC <--> Non-centrality parameter of the non-central
C                chi-square distribution.
C                Input range: [0, +infinity).
C                Search range: [0,1E4]
C                    DOUBLE PRECISION PNONC
C
C     STATUS <-- 0 if calculation completed correctly
C               -I if input parameter number I is out of range
C                1 if answer appears to be lower than lowest
C                  search bound
C                2 if answer appears to be higher than greatest
C                  search bound
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
C     Formula  26.4.25   of   Abramowitz   and   Stegun,  Handbook  of
C     Mathematical  Functions (1966) is used to compute the cumulative
C     distribution function.
C
C     Computation of other parameters involve a search for a value that
C     produces  the desired  value  of P.   The search relies  on  the
C     monotinicity of P with the other parameter.
C
C
C                            WARNING
C
C     The computation time  required for this  routine is proportional
C     to the noncentrality  parameter  (PNONC).  Very large  values of
C     this parameter can consume immense  computer resources.  This is
C     why the search range is bounded by 1e9.
C
C**********************************************************************
C     .. Parameters ..
      DOUBLE PRECISION tent9
      PARAMETER (tent9=1.0D9)
      DOUBLE PRECISION tol
      PARAMETER (tol=1.0D-8)
      DOUBLE PRECISION atol
      PARAMETER (atol=1.0D-50)
      DOUBLE PRECISION zero,one,inf
      PARAMETER (zero=1.0D-300,one=1.0D0-1.0D-16,inf=1.0D300)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION bound,df,p,pnonc,q,x
      INTEGER status,which
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ccum,cum,fx
      LOGICAL qhi,qleft
C     ..
C     .. External Subroutines ..
      EXTERNAL cumchn,dinvr,dstinv
C     ..
      IF (x.GT.inf) THEN
          x = inf
      END IF
      IF (df.GT.inf) THEN
          df = inf
      END IF
      IF (pnonc.GT.tent9) THEN
          pnonc = tent9
      END IF

      IF (.NOT. ((which.LT.1).OR. (which.GT.4))) GO TO 30
      IF (.NOT. (which.LT.1)) GO TO 10
      bound = 1.0D0
      GO TO 20

   10 bound = 4.0D0
   20 status = -1
      RETURN

   30 IF (which.EQ.1) GO TO 70
      IF (.NOT. ((p.LT.0.0D0).OR. (p.GT.one))) GO TO 60
      IF (.NOT. (p.LT.0.0D0)) GO TO 40
      bound = 0.0D0
      GO TO 50

   40 bound = one
   50 status = -2
      RETURN

   60 CONTINUE
   70 IF (which.EQ.2) GO TO 90
      IF (x.GE.0.0D0) GO TO 80
      bound = 0.0D0
      status = -4
      RETURN

   80 CONTINUE
   90 IF (which.EQ.3) GO TO 110
      IF (df.GT.0.0D0) GO TO 100
      bound = 0.0D0
      status = -5
      RETURN

  100 CONTINUE
  110 IF (which.EQ.4) GO TO 130
      IF (pnonc.GE.0.0D0) GO TO 120
      bound = 0.0D0
      status = -6
      RETURN

  120 CONTINUE
  130 IF ((1).EQ. (which)) THEN
          CALL cumchn(x,df,pnonc,p,q)
          status = 0

      ELSE IF ((2).EQ. (which)) THEN
          x = 5.0D0
          CALL dstinv(0.0D0,inf,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,x,fx,qleft,qhi)
  140     IF (.NOT. (status.EQ.1)) GO TO 150
          CALL cumchn(x,df,pnonc,cum,ccum)
          fx = cum - p
          CALL dinvr(status,x,fx,qleft,qhi)
          GO TO 140

  150     IF (.NOT. (status.EQ.-1)) GO TO 180
          IF (.NOT. (qleft)) GO TO 160
          status = 1
          bound = 0.0D0
          GO TO 170

  160     status = 2
          bound = inf
  170     CONTINUE
  180     CONTINUE

      ELSE IF ((3).EQ. (which)) THEN
          df = 5.0D0
          CALL dstinv(zero,inf,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,df,fx,qleft,qhi)
  190     IF (.NOT. (status.EQ.1)) GO TO 200
          CALL cumchn(x,df,pnonc,cum,ccum)
          fx = cum - p
          CALL dinvr(status,df,fx,qleft,qhi)
          GO TO 190

  200     IF (.NOT. (status.EQ.-1)) GO TO 230
          IF (.NOT. (qleft)) GO TO 210
          status = 1
          bound = zero
          GO TO 220

  210     status = 2
          bound = inf
  220     CONTINUE
  230     CONTINUE

      ELSE IF ((4).EQ. (which)) THEN
          pnonc = 5.0D0
          CALL dstinv(0.0D0,tent9,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,pnonc,fx,qleft,qhi)
  240     IF (.NOT. (status.EQ.1)) GO TO 250
          CALL cumchn(x,df,pnonc,cum,ccum)
          fx = cum - p
          CALL dinvr(status,pnonc,fx,qleft,qhi)
          GO TO 240

  250     IF (.NOT. (status.EQ.-1)) GO TO 280
          IF (.NOT. (qleft)) GO TO 260
          status = 1
          bound = zero
          GO TO 270

  260     status = 2
          bound = tent9
  270     CONTINUE
  280 END IF

      RETURN

      END
