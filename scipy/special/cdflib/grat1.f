      SUBROUTINE grat1(a,x,r,p,q,eps)
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,eps,p,q,r,x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a2n,a2nm1,am0,an,an0,b2n,b2nm1,c,cma,g,h,j,l,sum,
     +                 t,tol,w,z
C     ..
C     .. External Functions ..
      DOUBLE PRECISION erf,erfc1,gam1,rexp
      EXTERNAL erf,erfc1,gam1,rexp
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dlog,exp,sqrt
C     ..
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C        EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS
C                      P(A,X) AND Q(A,X)
C
C     IT IS ASSUMED THAT A .LE. 1.  EPS IS THE TOLERANCE TO BE USED.
C     THE INPUT ARGUMENT R HAS THE VALUE E**(-X)*X**A/GAMMA(A).
C-----------------------------------------------------------------------
      IF (a*x.EQ.0.0D0) GO TO 120
      IF (a.EQ.0.5D0) GO TO 100
      IF (x.LT.1.1D0) GO TO 10
      GO TO 60
C
C             TAYLOR SERIES FOR P(A,X)/X**A
C
   10 an = 3.0D0
      c = x
      sum = x/ (a+3.0D0)
      tol = 0.1D0*eps/ (a+1.0D0)
   20 an = an + 1.0D0
      c = -c* (x/an)
      t = c/ (a+an)
      sum = sum + t
      IF (abs(t).GT.tol) GO TO 20
      j = a*x* ((sum/6.0D0-0.5D0/ (a+2.0D0))*x+1.0D0/ (a+1.0D0))
C
      z = a*dlog(x)
      h = gam1(a)
      g = 1.0D0 + h
      IF (x.LT.0.25D0) GO TO 30
      IF (a.LT.x/2.59D0) GO TO 50
      GO TO 40

   30 IF (z.GT.-.13394D0) GO TO 50
C
   40 w = exp(z)
      p = w*g* (0.5D0+ (0.5D0-j))
      q = 0.5D0 + (0.5D0-p)
      RETURN
C
   50 l = rexp(z)
      w = 0.5D0 + (0.5D0+l)
      q = (w*j-l)*g - h
      IF (q.LT.0.0D0) GO TO 90
      p = 0.5D0 + (0.5D0-q)
      RETURN
C
C              CONTINUED FRACTION EXPANSION
C
   60 a2nm1 = 1.0D0
      a2n = 1.0D0
      b2nm1 = x
      b2n = x + (1.0D0-a)
      c = 1.0D0
   70 a2nm1 = x*a2n + c*a2nm1
      b2nm1 = x*b2n + c*b2nm1
      am0 = a2nm1/b2nm1
      c = c + 1.0D0
      cma = c - a
      a2n = a2nm1 + cma*a2n
      b2n = b2nm1 + cma*b2n
      an0 = a2n/b2n
      IF (abs(an0-am0).GE.eps*an0) GO TO 70
      q = r*an0
      p = 0.5D0 + (0.5D0-q)
      RETURN
C
C                SPECIAL CASES
C
   80 p = 0.0D0
      q = 1.0D0
      RETURN
C
   90 p = 1.0D0
      q = 0.0D0
      RETURN
C
  100 IF (x.GE.0.25D0) GO TO 110
      p = erf(sqrt(x))
      q = 0.5D0 + (0.5D0-p)
      RETURN

  110 q = erfc1(0,sqrt(x))
      p = 0.5D0 + (0.5D0-q)
      RETURN
C
  120 IF (x.LE.a) GO TO 80
      GO TO 90

      END
