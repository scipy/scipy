      SUBROUTINE bgrat(a,b,x,y,w,eps,ierr)
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR IX(A,B) WHEN A IS LARGER THAN B.
C     THE RESULT OF THE EXPANSION IS ADDED TO W. IT IS ASSUMED
C     THAT A .GE. 15 AND B .LE. 1.  EPS IS THE TOLERANCE USED.
C     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,eps,w,x,y
      INTEGER ierr
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION bm1,bp2n,cn,coef,dj,j,l,lnx,n2,nu,p,q,r,s,sum,t,
     +                 t2,u,v,z
      INTEGER i,n,nm1
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION c(30),d(30)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION algdiv,alnrel,gam1
      EXTERNAL algdiv,alnrel,gam1
C     ..
C     .. External Subroutines ..
      EXTERNAL grat1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dlog,exp
C     ..
C     .. Executable Statements ..
C
      bm1 = (b-0.5D0) - 0.5D0
      nu = a + 0.5D0*bm1
      IF (y.GT.0.375D0) GO TO 10
      lnx = alnrel(-y)
      GO TO 20

   10 lnx = dlog(x)
   20 z = -nu*lnx
      IF (b*z.EQ.0.0D0) GO TO 70
C
C                 COMPUTATION OF THE EXPANSION
C                 SET R = EXP(-Z)*Z**B/GAMMA(B)
C
      r = b* (1.0D0+gam1(b))*exp(b*dlog(z))
      r = r*exp(a*lnx)*exp(0.5D0*bm1*lnx)
      u = algdiv(b,a) + b*dlog(nu)
      u = r*exp(-u)
      IF (u.EQ.0.0D0) GO TO 70
      CALL grat1(b,z,r,p,q,eps)
C
      v = 0.25D0* (1.0D0/nu)**2
      t2 = 0.25D0*lnx*lnx
      l = w/u
      j = q/r
      sum = j
      t = 1.0D0
      cn = 1.0D0
      n2 = 0.0D0
      DO 50 n = 1,30
          bp2n = b + n2
          j = (bp2n* (bp2n+1.0D0)*j+ (z+bp2n+1.0D0)*t)*v
          n2 = n2 + 2.0D0
          t = t*t2
          cn = cn/ (n2* (n2+1.0D0))
          c(n) = cn
          s = 0.0D0
          IF (n.EQ.1) GO TO 40
          nm1 = n - 1
          coef = b - n
          DO 30 i = 1,nm1
              s = s + coef*c(i)*d(n-i)
              coef = coef + b
   30     CONTINUE
   40     d(n) = bm1*cn + s/n
          dj = d(n)*j
          sum = sum + dj
          IF (sum.LE.0.0D0) GO TO 70
          IF (abs(dj).LE.eps* (sum+l)) GO TO 60
   50 CONTINUE
C
C                    ADD THE RESULTS TO W
C
   60 ierr = 0
      w = w + u*sum
      RETURN
C
C               THE EXPANSION CANNOT BE COMPUTED
C
   70 ierr = 1
      RETURN

      END
