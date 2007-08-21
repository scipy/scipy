      DOUBLE PRECISION FUNCTION stvaln(p)
C
C**********************************************************************
C
C     DOUBLE PRECISION FUNCTION STVALN(P)
C                    STarting VALue for Neton-Raphon
C                calculation of Normal distribution Inverse
C
C
C                              Function
C
C
C     Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from -
C     infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P
C
C
C                              Arguments
C
C
C     P --> The probability whose normal deviate is sought.
C                    P is DOUBLE PRECISION
C
C
C                              Method
C
C
C     The  rational   function   on  page 95    of Kennedy  and  Gentle,
C     Statistical Computing, Marcel Dekker, NY , 1980.
C
C**********************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION p
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION sign,y,z
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION xden(5),xnum(5)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION devlpl
      EXTERNAL devlpl
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dble,log,sqrt
C     ..
C     .. Data statements ..
      DATA xnum/-0.322232431088D0,-1.000000000000D0,-0.342242088547D0,
     +     -0.204231210245D-1,-0.453642210148D-4/
      DATA xden/0.993484626060D-1,0.588581570495D0,0.531103462366D0,
     +     0.103537752850D0,0.38560700634D-2/
C     ..
C     .. Executable Statements ..
      IF (.NOT. (p.LE.0.5D0)) GO TO 10
      sign = -1.0D0
      z = p
      GO TO 20

   10 sign = 1.0D0
      z = 1.0D0 - p
   20 y = sqrt(-2.0D0*log(z))
      stvaln = y + devlpl(xnum,5,y)/devlpl(xden,5,y)
      stvaln = sign*stvaln
      RETURN

      END
