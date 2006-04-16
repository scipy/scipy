      SUBROUTINE cdftnc(which,p,q,t,df,pnonc,status,bound)
C***********************************************************************
C
C      SUBROUTINE CDFTNC( WHICH, P, Q, T, DF, PNONC, STATUS, BOUND )
C               Cumulative Distribution Function
C                  Non-Central T distribution
C
C                               Function
C
C     Calculates any one parameter of the noncentral t distribution give
C     values for the others.
C
C                               Arguments
C
C     WHICH --> Integer indicating which  argument
C               values is to be calculated from the others.
C               Legal range: 1..3
C               iwhich = 1 : Calculate P and Q from T,DF,PNONC
C               iwhich = 2 : Calculate T from P,Q,DF,PNONC
C               iwhich = 3 : Calculate DF from P,Q,T
C               iwhich = 4 : Calculate PNONC from P,Q,DF,T
C                    INTEGER WHICH
C
C        P <--> The integral from -infinity to t of the noncentral t-den
C              Input range: (0,1].
C                    DOUBLE PRECISION P
C
C     Q <--> 1-P.
C            Input range: (0, 1].
C            P + Q = 1.0.
C                    DOUBLE PRECISION Q
C
C        T <--> Upper limit of integration of the noncentral t-density.
C               Input range: ( -infinity, +infinity).
C               Search range: [ -1E100, 1E100 ]
C                    DOUBLE PRECISION T
C
C        DF <--> Degrees of freedom of the noncentral t-distribution.
C                Input range: (0 , +infinity).
C                Search range: [1e-100, 1E10]
C                    DOUBLE PRECISION DF
C
C     PNONC <--> Noncentrality parameter of the noncentral t-distributio
C                Input range: [-infinity , +infinity).
C                Search range: [-1e4, 1E4]
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
C                                Method
C
C     Upper tail    of  the  cumulative  noncentral t is calculated usin
C     formulae  from page 532  of Johnson, Kotz,  Balakrishnan, Coninuou
C     Univariate Distributions, Vol 2, 2nd Edition.  Wiley (1995)
C
C     Computation of other parameters involve a seach for a value that
C     produces  the desired  value  of P.   The search relies  on  the
C     monotinicity of P with the other parameter.
C
C***********************************************************************
C     .. Parameters ..
      DOUBLE PRECISION tent4
      PARAMETER (tent4=1.0D4)
      DOUBLE PRECISION tol
      PARAMETER (tol=1.0D-8)
      DOUBLE PRECISION atol
      PARAMETER (atol=1.0D-50)
      DOUBLE PRECISION zero,one,inf
      PARAMETER (zero=1.0D-100,one=1.0D0-1.0D-16,inf=1.0D100)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION bound,df,p,pnonc,q,t
      INTEGER status,which
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ccum,cum,fx
      LOGICAL qhi,qleft
C     ..
C     .. External Subroutines ..
      EXTERNAL cumtnc,dinvr,dstinv
C     ..
      IF (.NOT. ((which.LT.1).OR. (which.GT.4))) GO TO 30
      IF (.NOT. (which.LT.1)) GO TO 10
      bound = 1.0D0
      GO TO 20

   10 bound = 5.0D0
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
   70 IF (which.EQ.3) GO TO 90
      IF (.NOT. (df.LE.0.0D0)) GO TO 80
      bound = 0.0D0
      status = -5
      RETURN

   80 CONTINUE
   90 IF (which.EQ.4) GO TO 100
  100 IF ((1).EQ. (which)) THEN
          CALL cumtnc(t,df,pnonc,p,q)
          status = 0

      ELSE IF ((2).EQ. (which)) THEN
          t = 5.0D0
          CALL dstinv(-inf,inf,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,t,fx,qleft,qhi)
  110     IF (.NOT. (status.EQ.1)) GO TO 120
          CALL cumtnc(t,df,pnonc,cum,ccum)
          fx = cum - p
          CALL dinvr(status,t,fx,qleft,qhi)
          GO TO 110

  120     IF (.NOT. (status.EQ.-1)) GO TO 150
          IF (.NOT. (qleft)) GO TO 130
          status = 1
          bound = -inf
          GO TO 140

  130     status = 2
          bound = inf
  140     CONTINUE
  150     CONTINUE

      ELSE IF ((3).EQ. (which)) THEN
          df = 5.0D0
          CALL dstinv(zero,tent4,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,df,fx,qleft,qhi)
  160     IF (.NOT. (status.EQ.1)) GO TO 170
          CALL cumtnc(t,df,pnonc,cum,ccum)
          fx = cum - p
          CALL dinvr(status,df,fx,qleft,qhi)
          GO TO 160

  170     IF (.NOT. (status.EQ.-1)) GO TO 200
          IF (.NOT. (qleft)) GO TO 180
          status = 1
          bound = zero
          GO TO 190

  180     status = 2
          bound = inf
  190     CONTINUE
  200     CONTINUE

      ELSE IF ((4).EQ. (which)) THEN
          pnonc = 5.0D0
          CALL dstinv(-tent4,tent4,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,pnonc,fx,qleft,qhi)
  210     IF (.NOT. (status.EQ.1)) GO TO 220
          CALL cumtnc(t,df,pnonc,cum,ccum)
          fx = cum - p
          CALL dinvr(status,pnonc,fx,qleft,qhi)
          GO TO 210

  220     IF (.NOT. (status.EQ.-1)) GO TO 250
          IF (.NOT. (qleft)) GO TO 230
          status = 1
          bound = 0.0D0
          GO TO 240

  230     status = 2
          bound = tent4
  240     CONTINUE
  250 END IF

      RETURN

      END
