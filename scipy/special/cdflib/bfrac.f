      DOUBLE PRECISION FUNCTION bfrac(a,b,x,y,lambda,eps)
C-----------------------------------------------------------------------
C     CONTINUED FRACTION EXPANSION FOR IX(A,B) WHEN A,B .GT. 1.
C     IT IS ASSUMED THAT  LAMBDA = (A + B)*Y - B.
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,eps,lambda,x,y
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION alpha,an,anp1,beta,bn,bnp1,c,c0,c1,e,n,p,r,r0,s,
     +                 t,w,yp1
C     ..
C     .. External Functions ..
      DOUBLE PRECISION brcomp
      EXTERNAL brcomp
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..
C     .. Executable Statements ..
C--------------------
      bfrac = brcomp(a,b,x,y)
      IF (bfrac.EQ.0.0D0) RETURN
C
      c = 1.0D0 + lambda
      c0 = b/a
      c1 = 1.0D0 + 1.0D0/a
      yp1 = y + 1.0D0
C
      n = 0.0D0
      p = 1.0D0
      s = a + 1.0D0
      an = 0.0D0
      bn = 1.0D0
      anp1 = 1.0D0
      bnp1 = c/c1
      r = c1/c
C
C        CONTINUED FRACTION CALCULATION
C
   10 n = n + 1.0D0
      t = n/a
      w = n* (b-n)*x
      e = a/s
      alpha = (p* (p+c0)*e*e)* (w*x)
      e = (1.0D0+t)/ (c1+t+t)
      beta = n + w/s + e* (c+n*yp1)
      p = 1.0D0 + t
      s = s + 2.0D0
C
C        UPDATE AN, BN, ANP1, AND BNP1
C
      t = alpha*an + beta*anp1
      an = anp1
      anp1 = t
      t = alpha*bn + beta*bnp1
      bn = bnp1
      bnp1 = t
C
      r0 = r
      r = anp1/bnp1
      IF (.NOT.(abs(r-r0).GT.eps*r)) GO TO 20
C
C        RESCALE AN, BN, ANP1, AND BNP1
C
      an = an/bnp1
      bn = bn/bnp1
      anp1 = r
      bnp1 = 1.0D0
      GO TO 10
C
C                 TERMINATION
C
   20 bfrac = bfrac*r
      RETURN

      END
