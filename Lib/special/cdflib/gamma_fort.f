      DOUBLE PRECISION FUNCTION gamma(a)
C-----------------------------------------------------------------------
C
C         EVALUATION OF THE GAMMA FUNCTION FOR REAL ARGUMENTS
C
C                           -----------
C
C     GAMMA(A) IS ASSIGNED THE VALUE 0 WHEN THE GAMMA FUNCTION CANNOT
C     BE COMPUTED.
C
C-----------------------------------------------------------------------
C     WRITTEN BY ALFRED H. MORRIS, JR.
C          NAVAL SURFACE WEAPONS CENTER
C          DAHLGREN, VIRGINIA
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION bot,d,g,lnx,pi,r1,r2,r3,r4,r5,s,t,top,w,x,z
      INTEGER i,j,m,n
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION p(7),q(7)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION exparg,spmpar
      EXTERNAL exparg,spmpar
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dble,dlog,exp,int,mod,sin
C     ..
C     .. Data statements ..
C--------------------------
C     D = 0.5*(LN(2*PI) - 1)
C--------------------------
C--------------------------
C--------------------------
      DATA pi/3.1415926535898D0/
      DATA d/.41893853320467274178D0/
      DATA p(1)/.539637273585445D-03/,p(2)/.261939260042690D-02/,
     +     p(3)/.204493667594920D-01/,p(4)/.730981088720487D-01/,
     +     p(5)/.279648642639792D+00/,p(6)/.553413866010467D+00/,
     +     p(7)/1.0D0/
      DATA q(1)/-.832979206704073D-03/,q(2)/.470059485860584D-02/,
     +     q(3)/.225211131035340D-01/,q(4)/-.170458969313360D+00/,
     +     q(5)/-.567902761974940D-01/,q(6)/.113062953091122D+01/,
     +     q(7)/1.0D0/
      DATA r1/.820756370353826D-03/,r2/-.595156336428591D-03/,
     +     r3/.793650663183693D-03/,r4/-.277777777770481D-02/,
     +     r5/.833333333333333D-01/
C     ..
C     .. Executable Statements ..
C--------------------------
      gamma = 0.0D0
      x = a
      IF (abs(a).GE.15.0D0) GO TO 110
C-----------------------------------------------------------------------
C            EVALUATION OF GAMMA(A) FOR ABS(A) .LT. 15
C-----------------------------------------------------------------------
      t = 1.0D0
      m = int(a) - 1
C
C     LET T BE THE PRODUCT OF A-J WHEN A .GE. 2
C
      IF (m) 40,30,10
   10 DO 20 j = 1,m
          x = x - 1.0D0
          t = x*t
   20 CONTINUE
   30 x = x - 1.0D0
      GO TO 80
C
C     LET T BE THE PRODUCT OF A+J WHEN A .LT. 1
C
   40 t = a
      IF (a.GT.0.0D0) GO TO 70
      m = -m - 1
      IF (m.EQ.0) GO TO 60
      DO 50 j = 1,m
          x = x + 1.0D0
          t = x*t
   50 CONTINUE
   60 x = (x+0.5D0) + 0.5D0
      t = x*t
      IF (t.EQ.0.0D0) RETURN
C
   70 CONTINUE
C
C     THE FOLLOWING CODE CHECKS IF 1/T CAN OVERFLOW. THIS
C     CODE MAY BE OMITTED IF DESIRED.
C
      IF (abs(t).GE.1.D-30) GO TO 80
      IF (abs(t)*spmpar(3).LE.1.0001D0) RETURN
      gamma = 1.0D0/t
      RETURN
C
C     COMPUTE GAMMA(1 + X) FOR  0 .LE. X .LT. 1
C
   80 top = p(1)
      bot = q(1)
      DO 90 i = 2,7
          top = p(i) + x*top
          bot = q(i) + x*bot
   90 CONTINUE
      gamma = top/bot
C
C     TERMINATION
C
      IF (a.LT.1.0D0) GO TO 100
      gamma = gamma*t
      RETURN

  100 gamma = gamma/t
      RETURN
C-----------------------------------------------------------------------
C            EVALUATION OF GAMMA(A) FOR ABS(A) .GE. 15
C-----------------------------------------------------------------------
  110 IF (abs(a).GE.1.D3) RETURN
      IF (a.GT.0.0D0) GO TO 120
      x = -a
      n = x
      t = x - n
      IF (t.GT.0.9D0) t = 1.0D0 - t
      s = sin(pi*t)/pi
      IF (mod(n,2).EQ.0) s = -s
      IF (s.EQ.0.0D0) RETURN
C
C     COMPUTE THE MODIFIED ASYMPTOTIC SUM
C
  120 t = 1.0D0/ (x*x)
      g = ((((r1*t+r2)*t+r3)*t+r4)*t+r5)/x
C
C     ONE MAY REPLACE THE NEXT STATEMENT WITH  LNX = ALOG(X)
C     BUT LESS ACCURACY WILL NORMALLY BE OBTAINED.
C
      lnx = dlog(x)
C
C     FINAL ASSEMBLY
C
      z = x
      g = (d+g) + (z-0.5D0)* (lnx-1.D0)
      w = g
      t = g - dble(w)
      IF (w.GT.0.99999D0*exparg(0)) RETURN
      gamma = exp(w)* (1.0D0+t)
      IF (a.LT.0.0D0) gamma = (1.0D0/ (gamma*s))/x
      RETURN

      END
