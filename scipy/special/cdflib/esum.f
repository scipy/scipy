      DOUBLE PRECISION FUNCTION esum(mu,x)
C-----------------------------------------------------------------------
C                    EVALUATION OF EXP(MU + X)
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
      INTEGER mu
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION w
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC exp
C     ..
C     .. Executable Statements ..

      IF (x.GT.0.0D0) GO TO 10
C
      IF (mu.LT.0) GO TO 20
      w = mu + x
      IF (w.GT.0.0D0) GO TO 20
      esum = exp(w)
      RETURN
C
   10 IF (mu.GT.0) GO TO 20
      w = mu + x
      IF (w.LT.0.0D0) GO TO 20
      esum = exp(w)
      RETURN
C
   20 w = mu
      esum = exp(w)*exp(x)
      RETURN

      END
