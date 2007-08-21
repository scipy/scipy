      DOUBLE PRECISION FUNCTION apser(a,b,x,eps)
C-----------------------------------------------------------------------
C     APSER YIELDS THE INCOMPLETE BETA RATIO I(SUB(1-X))(B,A) FOR
C     A .LE. MIN(EPS,EPS*B), B*X .LE. 1, AND X .LE. 0.5. USED WHEN
C     A IS VERY SMALL. USE ONLY IF ABOVE INEQUALITIES ARE SATISFIED.
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,eps,x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION aj,bx,c,g,j,s,t,tol
C     ..
C     .. External Functions ..
      DOUBLE PRECISION psi
      EXTERNAL psi
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dlog
C     ..
C     .. Data statements ..
C--------------------
      DATA g/.577215664901533D0/
C     ..
C     .. Executable Statements ..
C--------------------
      bx = b*x
      t = x - bx
      IF (b*eps.GT.2.D-2) GO TO 10
      c = dlog(x) + psi(b) + g + t
      GO TO 20

   10 c = dlog(bx) + g + t
C
   20 tol = 5.0D0*eps*abs(c)
      j = 1.0D0
      s = 0.0D0
   30 j = j + 1.0D0
      t = t* (x-bx/j)
      aj = t/j
      s = s + aj
      IF (abs(aj).GT.tol) GO TO 30
C
      apser = -a* (c+s)
      RETURN

      END
