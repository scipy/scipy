      DOUBLE PRECISION FUNCTION bcorr(a0,b0)
C-----------------------------------------------------------------------
C
C     EVALUATION OF  DEL(A0) + DEL(B0) - DEL(A0 + B0)  WHERE
C     LN(GAMMA(A)) = (A - 0.5)*LN(A) - A + 0.5*LN(2*PI) + DEL(A).
C     IT IS ASSUMED THAT A0 .GE. 8 AND B0 .GE. 8.
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a0,b0
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a,b,c,c0,c1,c2,c3,c4,c5,h,s11,s3,s5,s7,s9,t,w,x,
     +                 x2
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dmax1,dmin1
C     ..
C     .. Data statements ..
      DATA c0/.833333333333333D-01/,c1/-.277777777760991D-02/,
     +     c2/.793650666825390D-03/,c3/-.595202931351870D-03/,
     +     c4/.837308034031215D-03/,c5/-.165322962780713D-02/
C     ..
C     .. Executable Statements ..
C------------------------
      a = dmin1(a0,b0)
      b = dmax1(a0,b0)
C
      h = a/b
      c = h/ (1.0D0+h)
      x = 1.0D0/ (1.0D0+h)
      x2 = x*x
C
C                SET SN = (1 - X**N)/(1 - X)
C
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
C                   COMPUTE  DEL(A) + W
C
      t = (1.0D0/a)**2
      bcorr = (((((c5*t+c4)*t+c3)*t+c2)*t+c1)*t+c0)/a + w
      RETURN

      END
