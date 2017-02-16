      DOUBLE PRECISION FUNCTION rcomp(a,x)
C     -------------------
C     EVALUATION OF EXP(-X)*X**A/GAMMA(A)
C     -------------------
C     RT2PIN = 1/SQRT(2*PI)
C     -------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION rt2pin,t,t1,u
C     ..
C     .. External Functions ..
      DOUBLE PRECISION gam1,gamma,rlog
      EXTERNAL gam1,gamma,rlog
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dlog,exp,sqrt
C     ..
C     .. Data statements ..
      DATA rt2pin/.398942280401433D0/
C     ..
C     .. Executable Statements ..
C     -------------------
      rcomp = 0.0D0
      IF (a.GE.20.0D0) GO TO 20
      t = a*dlog(x) - x
      IF (a.GE.1.0D0) GO TO 10
      rcomp = (a*exp(t))* (1.0D0+gam1(a))
      RETURN

   10 rcomp = exp(t)/gamma(a)
      RETURN
C
   20 u = x/a
      IF (u.EQ.0.0D0) RETURN
      t = (1.0D0/a)**2
      t1 = (((0.75D0*t-1.0D0)*t+3.5D0)*t-105.0D0)/ (a*1260.0D0)
      t1 = t1 - a*rlog(u)
      rcomp = rt2pin*sqrt(a)*exp(t1)
      RETURN

      END
