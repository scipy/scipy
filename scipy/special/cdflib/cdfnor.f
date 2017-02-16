      SUBROUTINE cdfnor(which,p,q,x,mean,sd,status,bound)
C**********************************************************************
C
C      SUBROUTINE CDFNOR( WHICH, P, Q, X, MEAN, SD, STATUS, BOUND )
C               Cumulative Distribution Function
C               NORmal distribution
C
C
C                              Function
C
C
C     Calculates any one parameter of the normal
C     distribution given values for the others.
C
C
C                              Arguments
C
C
C     WHICH  --> Integer indicating  which of the  next  parameter
C     values is to be calculated using values  of the others.
C     Legal range: 1..4
C               iwhich = 1 : Calculate P and Q from X,MEAN and SD
C               iwhich = 2 : Calculate X from P,Q,MEAN and SD
C               iwhich = 3 : Calculate MEAN from P,Q,X and SD
C               iwhich = 4 : Calculate SD from P,Q,X and MEAN
C                    INTEGER WHICH
C
C     P <--> The integral from -infinity to X of the normal density.
C            Input range: (0,1].
C                    DOUBLE PRECISION P
C
C     Q <--> 1-P.
C            Input range: (0, 1].
C            P + Q = 1.0.
C                    DOUBLE PRECISION Q
C
C     X < --> Upper limit of integration of the normal-density.
C             Input range: ( -infinity, +infinity)
C                    DOUBLE PRECISION X
C
C     MEAN <--> The mean of the normal density.
C               Input range: (-infinity, +infinity)
C                    DOUBLE PRECISION MEAN
C
C     SD <--> Standard Deviation of the normal density.
C             Input range: (0, +infinity).
C                    DOUBLE PRECISION SD
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
C
C
C     A slightly modified version of ANORM from
C
C     Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
C     Package of Special Function Routines and Test Drivers"
C     acm Transactions on Mathematical Software. 19, 22-32.
C
C     is used to calulate the  cumulative standard normal distribution.
C
C     The rational functions from pages  90-95  of Kennedy and Gentle,
C     Statistical  Computing,  Marcel  Dekker, NY,  1980 are  used  as
C     starting values to Newton's Iterations which compute the inverse
C     standard normal.  Therefore no  searches  are necessary for  any
C     parameter.
C
C     For X < -15, the asymptotic expansion for the normal is used  as
C     the starting value in finding the inverse standard normal.
C     This is formula 26.2.12 of Abramowitz and Stegun.
C
C
C                              Note
C
C
C      The normal density is proportional to
C      exp( - 0.5 * (( X - MEAN)/SD)**2)
C
C
C**********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION bound,mean,p,q,sd,x
      INTEGER status,which
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION pq,z
C     ..
C     .. External Functions ..
      DOUBLE PRECISION dinvnr,spmpar
      EXTERNAL dinvnr,spmpar
C     ..
C     .. External Subroutines ..
      EXTERNAL cumnor
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..
      status = 0
      IF (.NOT. ((which.LT.1).OR. (which.GT.4))) GO TO 30
      IF (.NOT. (which.LT.1)) GO TO 10
      bound = 1.0D0
      GO TO 20

   10 bound = 4.0D0
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
  110 IF (which.EQ.1) GO TO 150
      pq = p + q
      IF (.NOT. (abs(((pq)-0.5D0)-0.5D0).GT.
     +    (3.0D0*spmpar(1)))) GO TO 140
      IF (.NOT. (pq.LT.0.0D0)) GO TO 120
      bound = 0.0D0
      GO TO 130

  120 bound = 1.0D0
  130 status = 3
      RETURN

  140 CONTINUE
  150 IF (which.EQ.4) GO TO 170
      IF (.NOT. (sd.LE.0.0D0)) GO TO 160
      bound = 0.0D0
      status = -6
      RETURN

  160 CONTINUE
  170 IF ((1).EQ. (which)) THEN
          z = (x-mean)/sd
          CALL cumnor(z,p,q)

      ELSE IF ((2).EQ. (which)) THEN
          z = dinvnr(p,q)
          x = sd*z + mean

      ELSE IF ((3).EQ. (which)) THEN
          z = dinvnr(p,q)
          mean = x - sd*z

      ELSE IF ((4).EQ. (which)) THEN
          z = dinvnr(p,q)
          sd = (x-mean)/z
      END IF

      RETURN

      END
