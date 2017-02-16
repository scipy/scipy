      DOUBLE PRECISION FUNCTION erfc1(ind,x)
C-----------------------------------------------------------------------
C         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION
C
C          ERFC1(IND,X) = ERFC(X)            IF IND = 0
C          ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
      INTEGER ind
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ax,bot,c,e,t,top,w
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION a(5),b(3),p(8),q(8),r(5),s(4)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION exparg
      EXTERNAL exparg
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dble,exp
C     ..
C     .. Data statements ..
C-------------------------
C-------------------------
C-------------------------
C-------------------------
      DATA c/.564189583547756D0/
      DATA a(1)/.771058495001320D-04/,a(2)/-.133733772997339D-02/,
     +     a(3)/.323076579225834D-01/,a(4)/.479137145607681D-01/,
     +     a(5)/.128379167095513D+00/
      DATA b(1)/.301048631703895D-02/,b(2)/.538971687740286D-01/,
     +     b(3)/.375795757275549D+00/
      DATA p(1)/-1.36864857382717D-07/,p(2)/5.64195517478974D-01/,
     +     p(3)/7.21175825088309D+00/,p(4)/4.31622272220567D+01/,
     +     p(5)/1.52989285046940D+02/,p(6)/3.39320816734344D+02/,
     +     p(7)/4.51918953711873D+02/,p(8)/3.00459261020162D+02/
      DATA q(1)/1.00000000000000D+00/,q(2)/1.27827273196294D+01/,
     +     q(3)/7.70001529352295D+01/,q(4)/2.77585444743988D+02/,
     +     q(5)/6.38980264465631D+02/,q(6)/9.31354094850610D+02/,
     +     q(7)/7.90950925327898D+02/,q(8)/3.00459260956983D+02/
      DATA r(1)/2.10144126479064D+00/,r(2)/2.62370141675169D+01/,
     +     r(3)/2.13688200555087D+01/,r(4)/4.65807828718470D+00/,
     +     r(5)/2.82094791773523D-01/
      DATA s(1)/9.41537750555460D+01/,s(2)/1.87114811799590D+02/,
     +     s(3)/9.90191814623914D+01/,s(4)/1.80124575948747D+01/
C     ..
C     .. Executable Statements ..
C-------------------------
C
C                     ABS(X) .LE. 0.5
C
      ax = abs(x)
      IF (ax.GT.0.5D0) GO TO 10
      t = x*x
      top = ((((a(1)*t+a(2))*t+a(3))*t+a(4))*t+a(5)) + 1.0D0
      bot = ((b(1)*t+b(2))*t+b(3))*t + 1.0D0
      erfc1 = 0.5D0 + (0.5D0-x* (top/bot))
      IF (ind.NE.0) erfc1 = exp(t)*erfc1
      RETURN
C
C                  0.5 .LT. ABS(X) .LE. 4
C
   10 IF (ax.GT.4.0D0) GO TO 20
      top = ((((((p(1)*ax+p(2))*ax+p(3))*ax+p(4))*ax+p(5))*ax+p(6))*ax+
     +      p(7))*ax + p(8)
      bot = ((((((q(1)*ax+q(2))*ax+q(3))*ax+q(4))*ax+q(5))*ax+q(6))*ax+
     +      q(7))*ax + q(8)
      erfc1 = top/bot
      GO TO 40
C
C                      ABS(X) .GT. 4
C
   20 IF (x.LE.-5.6D0) GO TO 60
      IF (ind.NE.0) GO TO 30
      IF (x.GT.100.0D0) GO TO 70
      IF (x*x.GT.-exparg(1)) GO TO 70
C
   30 t = (1.0D0/x)**2
      top = (((r(1)*t+r(2))*t+r(3))*t+r(4))*t + r(5)
      bot = (((s(1)*t+s(2))*t+s(3))*t+s(4))*t + 1.0D0
      erfc1 = (c-t*top/bot)/ax
C
C                      FINAL ASSEMBLY
C
   40 IF (ind.EQ.0) GO TO 50
      IF (x.LT.0.0D0) erfc1 = 2.0D0*exp(x*x) - erfc1
      RETURN

   50 w = dble(x)*dble(x)
      t = w
      e = w - dble(t)
      erfc1 = ((0.5D0+ (0.5D0-e))*exp(-t))*erfc1
      IF (x.LT.0.0D0) erfc1 = 2.0D0 - erfc1
      RETURN
C
C             LIMIT VALUE FOR LARGE NEGATIVE X
C
   60 erfc1 = 2.0D0
      IF (ind.NE.0) erfc1 = 2.0D0*exp(x*x)
      RETURN
C
C             LIMIT VALUE FOR LARGE POSITIVE X
C                       WHEN IND = 0
C
   70 erfc1 = 0.0D0
      RETURN

      END
