      DOUBLE PRECISION FUNCTION basym(a,b,lambda,eps)
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR IX(A,B) FOR LARGE A AND B.
C     LAMBDA = (A + B)*Y - B  AND EPS IS THE TOLERANCE USED.
C     IT IS ASSUMED THAT LAMBDA IS NONNEGATIVE AND THAT
C     A AND B ARE GREATER THAN OR EQUAL TO 15.
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,eps,lambda
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION bsum,dsum,e0,e1,f,h,h2,hn,j0,j1,r,r0,r1,s,sum,t,
     +                 t0,t1,u,w,w0,z,z0,z2,zn,znm1
      INTEGER i,im1,imj,j,m,mm1,mmj,n,np1,num
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION a0(21),b0(21),c(21),d(21)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION bcorr,erfc1,rlog1
      EXTERNAL bcorr,erfc1,rlog1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp,sqrt
C     ..
C     .. Data statements ..
C------------------------
C     ****** NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP
C            ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN.
C            THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1.
C
C------------------------
C     E0 = 2/SQRT(PI)
C     E1 = 2**(-3/2)
C------------------------
      DATA num/20/
      DATA e0/1.12837916709551D0/,e1/.353553390593274D0/
C     ..
C     .. Executable Statements ..
C------------------------
      basym = 0.0D0
      IF (a.GE.b) GO TO 10
      h = a/b
      r0 = 1.0D0/ (1.0D0+h)
      r1 = (b-a)/b
      w0 = 1.0D0/sqrt(a* (1.0D0+h))
      GO TO 20

   10 h = b/a
      r0 = 1.0D0/ (1.0D0+h)
      r1 = (b-a)/a
      w0 = 1.0D0/sqrt(b* (1.0D0+h))
C
   20 f = a*rlog1(-lambda/a) + b*rlog1(lambda/b)
      t = exp(-f)
      IF (t.EQ.0.0D0) RETURN
      z0 = sqrt(f)
      z = 0.5D0* (z0/e1)
      z2 = f + f
C
      a0(1) = (2.0D0/3.0D0)*r1
      c(1) = -0.5D0*a0(1)
      d(1) = -c(1)
      j0 = (0.5D0/e0)*erfc1(1,z0)
      j1 = e1
      sum = j0 + d(1)*w0*j1
C
      s = 1.0D0
      h2 = h*h
      hn = 1.0D0
      w = w0
      znm1 = z
      zn = z2
      DO 70 n = 2,num,2
          hn = h2*hn
          a0(n) = 2.0D0*r0* (1.0D0+h*hn)/ (n+2.0D0)
          np1 = n + 1
          s = s + hn
          a0(np1) = 2.0D0*r1*s/ (n+3.0D0)
C
          DO 60 i = n,np1
              r = -0.5D0* (i+1.0D0)
              b0(1) = r*a0(1)
              DO 40 m = 2,i
                  bsum = 0.0D0
                  mm1 = m - 1
                  DO 30 j = 1,mm1
                      mmj = m - j
                      bsum = bsum + (j*r-mmj)*a0(j)*b0(mmj)
   30             CONTINUE
                  b0(m) = r*a0(m) + bsum/m
   40         CONTINUE
              c(i) = b0(i)/ (i+1.0D0)
C
              dsum = 0.0D0
              im1 = i - 1
              DO 50 j = 1,im1
                  imj = i - j
                  dsum = dsum + d(imj)*c(j)
   50         CONTINUE
              d(i) = - (dsum+c(i))
   60     CONTINUE
C
          j0 = e1*znm1 + (n-1.0D0)*j0
          j1 = e1*zn + n*j1
          znm1 = z2*znm1
          zn = z2*zn
          w = w0*w
          t0 = d(n)*w*j0
          w = w0*w
          t1 = d(np1)*w*j1
          sum = sum + (t0+t1)
          IF ((abs(t0)+abs(t1)).LE.eps*sum) GO TO 80
   70 CONTINUE
C
   80 u = exp(-bcorr(a,b))
      basym = e0*t*u*sum
      RETURN

      END
