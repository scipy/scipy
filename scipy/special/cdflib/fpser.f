      DOUBLE PRECISION FUNCTION fpser(a,b,x,eps)
C-----------------------------------------------------------------------
C
C                 EVALUATION OF I (A,B)
C                                X
C
C          FOR B .LT. MIN(EPS,EPS*A) AND X .LE. 0.5.
C
C-----------------------------------------------------------------------
C
C                  SET  FPSER = X**A
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,eps,x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION an,c,s,t,tol
C     ..
C     .. External Functions ..
      DOUBLE PRECISION exparg
      EXTERNAL exparg
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dlog,exp
C     ..
C     .. Executable Statements ..

      fpser = 1.0D0
      IF (a.LE.1.D-3*eps) GO TO 10
      fpser = 0.0D0
      t = a*dlog(x)
      IF (t.LT.exparg(1)) RETURN
      fpser = exp(t)
C
C                NOTE THAT 1/B(A,B) = B
C
   10 fpser = (b/a)*fpser
      tol = eps/a
      an = a + 1.0D0
      t = x
      s = t/an
   20 an = an + 1.0D0
      t = x*t
      c = t/an
      s = s + c
      IF (abs(c).GT.tol) GO TO 20
C
      fpser = fpser* (1.0D0+a*s)
      RETURN

      END
