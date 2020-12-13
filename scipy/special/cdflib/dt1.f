      DOUBLE PRECISION FUNCTION dt1(p,q,df)
C**********************************************************************
C
C     DOUBLE PRECISION FUNCTION DT1(P,Q,DF)
C     Double precision Initialize Approximation to
C           INVerse of the cumulative T distribution
C
C
C                              Function
C
C
C     Returns  the  inverse   of  the T   distribution   function, i.e.,
C     the integral from 0 to INVT of the T density is P. This is an
C     initial approximation
C
C
C                              Arguments
C
C
C     P --> The p-value whose inverse from the T distribution is
C          desired.
C                    P is DOUBLE PRECISION
C
C     Q --> 1-P.
C                    Q is DOUBLE PRECISION
C
C     DF --> Degrees of freedom of the T distribution.
C                    DF is DOUBLE PRECISION
C
C**********************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION df,p,q
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION denpow,sum,term,x,xp,xx
      INTEGER i
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION coef(5,4),denom(4)
      INTEGER ideg(4)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION dinvnr,devlpl
      EXTERNAL dinvnr,devlpl
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..
C     .. Data statements ..
      DATA (coef(i,1),i=1,5)/1.0D0,1.0D0,3*0.0D0/
      DATA (coef(i,2),i=1,5)/3.0D0,16.0D0,5.0D0,2*0.0D0/
      DATA (coef(i,3),i=1,5)/-15.0D0,17.0D0,19.0D0,3.0D0,0.0D0/
      DATA (coef(i,4),i=1,5)/-945.0D0,-1920.0D0,1482.0D0,776.0D0,79.0D0/
      DATA ideg/2,3,4,5/
      DATA denom/4.0D0,96.0D0,384.0D0,92160.0D0/
C     ..
C     .. Executable Statements ..
      x = abs(dinvnr(p,q))
      xx = x*x
      sum = x
      denpow = 1.0D0
      DO 10,i = 1,4
          term = devlpl(coef(1,i),ideg(i),xx)*x
          denpow = denpow*df
          sum = sum + term/ (denpow*denom(i))
   10 CONTINUE
      IF (.NOT. (p.GE.0.5D0)) GO TO 20
      xp = sum
      GO TO 30

   20 xp = -sum
   30 dt1 = xp
      RETURN

      END
