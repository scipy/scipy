      DOUBLE PRECISION FUNCTION brcmp1(mu,a,b,x,y)
C-----------------------------------------------------------------------
C          EVALUATION OF  EXP(MU) * (X**A*Y**B/BETA(A,B))
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,x,y
      INTEGER mu
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a0,apb,b0,c,const,e,h,lambda,lnx,lny,t,u,v,x0,y0,
     +                 z
      INTEGER i,n
C     ..
C     .. External Functions ..
      DOUBLE PRECISION algdiv,alnrel,bcorr,betaln,esum,gam1,gamln1,rlog1
      EXTERNAL algdiv,alnrel,bcorr,betaln,esum,gam1,gamln1,rlog1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dble,dlog,dmax1,dmin1,exp,sqrt
C     ..
C     .. Data statements ..
C-----------------
C     CONST = 1/SQRT(2*PI)
C-----------------
      DATA const/.398942280401433D0/
C     ..
C     .. Executable Statements ..
C
      a0 = dmin1(a,b)
      IF (a0.GE.8.0D0) GO TO 130
C
      IF (x.GT.0.375D0) GO TO 10
      lnx = dlog(x)
      lny = alnrel(-x)
      GO TO 30

   10 IF (y.GT.0.375D0) GO TO 20
      lnx = alnrel(-y)
      lny = dlog(y)
      GO TO 30

   20 lnx = dlog(x)
      lny = dlog(y)
C
   30 z = a*lnx + b*lny
      IF (a0.LT.1.0D0) GO TO 40
      z = z - betaln(a,b)
      brcmp1 = esum(mu,z)
      RETURN
C-----------------------------------------------------------------------
C              PROCEDURE FOR A .LT. 1 OR B .LT. 1
C-----------------------------------------------------------------------
   40 b0 = dmax1(a,b)
      IF (b0.GE.8.0D0) GO TO 120
      IF (b0.GT.1.0D0) GO TO 70
C
C                   ALGORITHM FOR B0 .LE. 1
C
      brcmp1 = esum(mu,z)
      IF (brcmp1.EQ.0.0D0) RETURN
C
      apb = a + b
      IF (apb.GT.1.0D0) GO TO 50
      z = 1.0D0 + gam1(apb)
      GO TO 60

   50 u = dble(a) + dble(b) - 1.D0
      z = (1.0D0+gam1(u))/apb
C
   60 c = (1.0D0+gam1(a))* (1.0D0+gam1(b))/z
      brcmp1 = brcmp1* (a0*c)/ (1.0D0+a0/b0)
      RETURN
C
C                ALGORITHM FOR 1 .LT. B0 .LT. 8
C
   70 u = gamln1(a0)
      n = b0 - 1.0D0
      IF (n.LT.1) GO TO 90
      c = 1.0D0
      DO 80 i = 1,n
          b0 = b0 - 1.0D0
          c = c* (b0/ (a0+b0))
   80 CONTINUE
      u = dlog(c) + u
C
   90 z = z - u
      b0 = b0 - 1.0D0
      apb = a0 + b0
      IF (apb.GT.1.0D0) GO TO 100
      t = 1.0D0 + gam1(apb)
      GO TO 110

  100 u = dble(a0) + dble(b0) - 1.D0
      t = (1.0D0+gam1(u))/apb
  110 brcmp1 = a0*esum(mu,z)* (1.0D0+gam1(b0))/t
      RETURN
C
C                   ALGORITHM FOR B0 .GE. 8
C
  120 u = gamln1(a0) + algdiv(a0,b0)
      brcmp1 = a0*esum(mu,z-u)
      RETURN
C-----------------------------------------------------------------------
C              PROCEDURE FOR A .GE. 8 AND B .GE. 8
C-----------------------------------------------------------------------
  130 IF (a.GT.b) GO TO 140
      h = a/b
      x0 = h/ (1.0D0+h)
      y0 = 1.0D0/ (1.0D0+h)
      lambda = a - (a+b)*x
      GO TO 150

  140 h = b/a
      x0 = 1.0D0/ (1.0D0+h)
      y0 = h/ (1.0D0+h)
      lambda = (a+b)*y - b
C
  150 e = -lambda/a
      IF (abs(e).GT.0.6D0) GO TO 160
      u = rlog1(e)
      GO TO 170

  160 u = e - dlog(x/x0)
C
  170 e = lambda/b
      IF (abs(e).GT.0.6D0) GO TO 180
      v = rlog1(e)
      GO TO 190

  180 v = e - dlog(y/y0)
C
  190 z = esum(mu,- (a*u+b*v))
      brcmp1 = const*sqrt(b*x0)*z*exp(-bcorr(a,b))
      RETURN

      END
