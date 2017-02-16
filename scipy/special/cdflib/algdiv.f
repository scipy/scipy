      DOUBLE PRECISION FUNCTION algdiv(a,b)
C-----------------------------------------------------------------------
C
C     COMPUTATION OF LN(GAMMA(B)/GAMMA(A+B)) WHEN B .GE. 8
C
C                         --------
C
C     IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY
C     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X).
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION c,c0,c1,c2,c3,c4,c5,d,h,s11,s3,s5,s7,s9,t,u,v,w,
     +                 x,x2
C     ..
C     .. External Functions ..
      DOUBLE PRECISION alnrel
      EXTERNAL alnrel
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dlog
C     ..
C     .. Data statements ..
      DATA c0/.833333333333333D-01/,c1/-.277777777760991D-02/,
     +     c2/.793650666825390D-03/,c3/-.595202931351870D-03/,
     +     c4/.837308034031215D-03/,c5/-.165322962780713D-02/
C     ..
C     .. Executable Statements ..
C------------------------
      IF (a.LE.b) GO TO 10
      h = b/a
      c = 1.0D0/ (1.0D0+h)
      x = h/ (1.0D0+h)
      d = a + (b-0.5D0)
      GO TO 20

   10 h = a/b
      c = h/ (1.0D0+h)
      x = 1.0D0/ (1.0D0+h)
      d = b + (a-0.5D0)
C
C                SET SN = (1 - X**N)/(1 - X)
C
   20 x2 = x*x
      s3 = 1.0D0 + (x+x2)
      s5 = 1.0D0 + (x+x2*s3)
      s7 = 1.0D0 + (x+x2*s5)
      s9 = 1.0D0 + (x+x2*s7)
      s11 = 1.0D0 + (x+x2*s9)
C
C                SET W = DEL(B) - DEL(A + B)
C
      t = (1.0D0/b)**2
      w = ((((c5*s11*t+c4*s9)*t+c3*s7)*t+c2*s5)*t+c1*s3)*t + c0
      w = w* (c/b)
C
C                    COMBINE THE RESULTS
C
      u = d*alnrel(a/b)
      v = a* (dlog(b)-1.0D0)
      IF (u.LE.v) GO TO 30
      algdiv = (w-v) - u
      RETURN

   30 algdiv = (w-u) - v
      RETURN

      END
