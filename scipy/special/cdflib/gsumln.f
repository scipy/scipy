      DOUBLE PRECISION FUNCTION gsumln(a,b)
C-----------------------------------------------------------------------
C          EVALUATION OF THE FUNCTION LN(GAMMA(A + B))
C          FOR 1 .LE. A .LE. 2  AND  1 .LE. B .LE. 2
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION x
C     ..
C     .. External Functions ..
      DOUBLE PRECISION alnrel,gamln1
      EXTERNAL alnrel,gamln1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dble,dlog
C     ..
C     .. Executable Statements ..
      x = dble(a) + dble(b) - 2.D0
      IF (x.GT.0.25D0) GO TO 10
      gsumln = gamln1(1.0D0+x)
      RETURN

   10 IF (x.GT.1.25D0) GO TO 20
      gsumln = gamln1(x) + alnrel(x)
      RETURN

   20 gsumln = gamln1(x-1.0D0) + dlog(x* (1.0D0+x))
      RETURN

      END
