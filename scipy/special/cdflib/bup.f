      DOUBLE PRECISION FUNCTION bup(a,b,x,y,n,eps)
C-----------------------------------------------------------------------
C     EVALUATION OF IX(A,B) - IX(A+N,B) WHERE N IS A POSITIVE INTEGER.
C     EPS IS THE TOLERANCE USED.
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,eps,x,y
      INTEGER n
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ap1,apb,d,l,r,t,w
      INTEGER i,k,kp1,mu,nm1
C     ..
C     .. External Functions ..
      DOUBLE PRECISION brcmp1,exparg
      EXTERNAL brcmp1,exparg
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp
C     ..
C     .. Executable Statements ..
C
C          OBTAIN THE SCALING FACTOR EXP(-MU) AND
C             EXP(MU)*(X**A*Y**B/BETA(A,B))/A
C
      apb = a + b
      ap1 = a + 1.0D0
      mu = 0
      d = 1.0D0
      IF (n.EQ.1 .OR. a.LT.1.0D0) GO TO 10
      IF (apb.LT.1.1D0*ap1) GO TO 10
      mu = abs(exparg(1))
      k = exparg(0)
      IF (k.LT.mu) mu = k
      t = mu
      d = exp(-t)
C
   10 bup = brcmp1(mu,a,b,x,y)/a
      IF (n.EQ.1 .OR. bup.EQ.0.0D0) RETURN
      nm1 = n - 1
      w = d
C
C          LET K BE THE INDEX OF THE MAXIMUM TERM
C
      k = 0
      IF (b.LE.1.0D0) GO TO 50
      IF (y.GT.1.D-4) GO TO 20
      k = nm1
      GO TO 30

   20 r = (b-1.0D0)*x/y - a
      IF (r.LT.1.0D0) GO TO 50
      k = nm1
      t = nm1
      IF (r.LT.t) k = r
C
C          ADD THE INCREASING TERMS OF THE SERIES
C
   30 DO 40 i = 1,k
          l = i - 1
          d = ((apb+l)/ (ap1+l))*x*d
          w = w + d
   40 CONTINUE
      IF (k.EQ.nm1) GO TO 70
C
C          ADD THE REMAINING TERMS OF THE SERIES
C
   50 kp1 = k + 1
      DO 60 i = kp1,nm1
          l = i - 1
          d = ((apb+l)/ (ap1+l))*x*d
          w = w + d
          IF (d.LE.eps*w) GO TO 70
   60 CONTINUE
C
C               TERMINATE THE PROCEDURE
C
   70 bup = bup*w
      RETURN

      END
