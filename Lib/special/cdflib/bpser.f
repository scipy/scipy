      DOUBLE PRECISION FUNCTION bpser(a,b,x,eps)
C-----------------------------------------------------------------------
C     POWER SERIES EXPANSION FOR EVALUATING IX(A,B) WHEN B .LE. 1
C     OR B*X .LE. 0.7.  EPS IS THE TOLERANCE USED.
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,eps,x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a0,apb,b0,c,n,sum,t,tol,u,w,z
      INTEGER i,m
C     ..
C     .. External Functions ..
      DOUBLE PRECISION algdiv,betaln,gam1,gamln1
      EXTERNAL algdiv,betaln,gam1,gamln1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dble,dlog,dmax1,dmin1,exp
C     ..
C     .. Executable Statements ..
C
      bpser = 0.0D0
      IF (x.EQ.0.0D0) RETURN
C-----------------------------------------------------------------------
C            COMPUTE THE FACTOR X**A/(A*BETA(A,B))
C-----------------------------------------------------------------------
      a0 = dmin1(a,b)
      IF (a0.LT.1.0D0) GO TO 10
      z = a*dlog(x) - betaln(a,b)
      bpser = exp(z)/a
      GO TO 100

   10 b0 = dmax1(a,b)
      IF (b0.GE.8.0D0) GO TO 90
      IF (b0.GT.1.0D0) GO TO 40
C
C            PROCEDURE FOR A0 .LT. 1 AND B0 .LE. 1
C
      bpser = x**a
      IF (bpser.EQ.0.0D0) RETURN
C
      apb = a + b
      IF (apb.GT.1.0D0) GO TO 20
      z = 1.0D0 + gam1(apb)
      GO TO 30

   20 u = dble(a) + dble(b) - 1.D0
      z = (1.0D0+gam1(u))/apb
C
   30 c = (1.0D0+gam1(a))* (1.0D0+gam1(b))/z
      bpser = bpser*c* (b/apb)
      GO TO 100
C
C         PROCEDURE FOR A0 .LT. 1 AND 1 .LT. B0 .LT. 8
C
   40 u = gamln1(a0)
      m = b0 - 1.0D0
      IF (m.LT.1) GO TO 60
      c = 1.0D0
      DO 50 i = 1,m
          b0 = b0 - 1.0D0
          c = c* (b0/ (a0+b0))
   50 CONTINUE
      u = dlog(c) + u
C
   60 z = a*dlog(x) - u
      b0 = b0 - 1.0D0
      apb = a0 + b0
      IF (apb.GT.1.0D0) GO TO 70
      t = 1.0D0 + gam1(apb)
      GO TO 80

   70 u = dble(a0) + dble(b0) - 1.D0
      t = (1.0D0+gam1(u))/apb
   80 bpser = exp(z)* (a0/a)* (1.0D0+gam1(b0))/t
      GO TO 100
C
C            PROCEDURE FOR A0 .LT. 1 AND B0 .GE. 8
C
   90 u = gamln1(a0) + algdiv(a0,b0)
      z = a*dlog(x) - u
      bpser = (a0/a)*exp(z)
  100 IF (bpser.EQ.0.0D0 .OR. a.LE.0.1D0*eps) RETURN
C-----------------------------------------------------------------------
C                     COMPUTE THE SERIES
C-----------------------------------------------------------------------
      sum = 0.0D0
      n = 0.0D0
      c = 1.0D0
      tol = eps/a
  110 n = n + 1.0D0
      c = c* (0.5D0+ (0.5D0-b/n))*x
      w = c/ (a+n)
      sum = sum + w
      IF (abs(w).GT.tol) GO TO 110
      bpser = bpser* (1.0D0+a*sum)
      RETURN

      END
