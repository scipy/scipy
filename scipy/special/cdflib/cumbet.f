      SUBROUTINE cumbet(x,y,a,b,cum,ccum)
C**********************************************************************
C
C     SUBROUTINE CUMBET(X,Y,A,B,CUM,CCUM)
C          Double precision cUMulative incomplete BETa distribution
C
C
C                              Function
C
C
C     Calculates the cdf to X of the incomplete beta distribution
C     with parameters a and b.  This is the integral from 0 to x
C     of (1/B(a,b))*f(t)) where f(t) = t**(a-1) * (1-t)**(b-1)
C
C
C                              Arguments
C
C
C     X --> Upper limit of integration.
C                                        X is DOUBLE PRECISION
C
C     Y --> 1 - X.
C                                        Y is DOUBLE PRECISION
C
C     A --> First parameter of the beta distribution.
C                                        A is DOUBLE PRECISION
C
C     B --> Second parameter of the beta distribution.
C                                        B is DOUBLE PRECISION
C
C     CUM <-- Cumulative incomplete beta distribution.
C                                        CUM is DOUBLE PRECISION
C
C     CCUM <-- Compliment of Cumulative incomplete beta distribution.
C                                        CCUM is DOUBLE PRECISION
C
C
C                              Method
C
C
C     Calls the routine BRATIO.
C
C                                   References
C
C     Didonato, Armido R. and Morris, Alfred H. Jr. (1992) Algorithim
C     708 Significant Digit Computation of the Incomplete Beta Function
C     Ratios. ACM ToMS, Vol.18, No. 3, Sept. 1992, 360-373.
C
C**********************************************************************

C     .. Scalar Arguments ..
      DOUBLE PRECISION x,y,a,b,cum,ccum
C     ..
C     .. Local Scalars ..
      INTEGER ierr
C     ..
C     .. External Routines ..
      EXTERNAL bratio
C     ..
C     .. Executable Statements ..
      IF (.NOT. (x.LE.0.0D0)) GO TO 10
      cum = 0.0D0
      ccum = 1.0D0
      RETURN

   10 IF (.NOT. (y.LE.0.0D0)) GO TO 20
      cum = 1.0D0
      ccum = 0.0D0
      RETURN

   20 CALL bratio(a,b,x,y,cum,ccum,ierr)

C     Call bratio routine


      RETURN

      END
