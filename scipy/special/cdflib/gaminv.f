      SUBROUTINE gaminv(a,x,x0,p,q,ierr)
C ----------------------------------------------------------------------
C            INVERSE INCOMPLETE GAMMA RATIO FUNCTION
C
C     GIVEN POSITIVE A, AND NONEGATIVE P AND Q WHERE P + Q = 1.
C     THEN X IS COMPUTED WHERE P(A,X) = P AND Q(A,X) = Q. SCHRODER
C     ITERATION IS EMPLOYED. THE ROUTINE ATTEMPTS TO COMPUTE X
C     TO 10 SIGNIFICANT DIGITS IF THIS IS POSSIBLE FOR THE
C     PARTICULAR COMPUTER ARITHMETIC BEING USED.
C
C                      ------------
C
C     X IS A VARIABLE. IF P = 0 THEN X IS ASSIGNED THE VALUE 0,
C     AND IF Q = 0 THEN X IS SET TO THE LARGEST FLOATING POINT
C     NUMBER AVAILABLE. OTHERWISE, GAMINV ATTEMPTS TO OBTAIN
C     A SOLUTION FOR P(A,X) = P AND Q(A,X) = Q. IF THE ROUTINE
C     IS SUCCESSFUL THEN THE SOLUTION IS STORED IN X.
C
C     X0 IS AN OPTIONAL INITIAL APPROXIMATION FOR X. IF THE USER
C     DOES NOT WISH TO SUPPLY AN INITIAL APPROXIMATION, THEN SET
C     X0 .LE. 0.
C
C     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
C     WHEN THE ROUTINE TERMINATES, IERR HAS ONE OF THE FOLLOWING
C     VALUES ...
C
C       IERR =  0    THE SOLUTION WAS OBTAINED. ITERATION WAS
C                    NOT USED.
C       IERR.GT.0    THE SOLUTION WAS OBTAINED. IERR ITERATIONS
C                    WERE PERFORMED.
C       IERR = -2    (INPUT ERROR) A .LE. 0
C       IERR = -3    NO SOLUTION WAS OBTAINED. THE RATIO Q/A
C                    IS TOO LARGE.
C       IERR = -4    (INPUT ERROR) P + Q .NE. 1
C       IERR = -6    20 ITERATIONS WERE PERFORMED. THE MOST
C                    RECENT VALUE OBTAINED FOR X IS GIVEN.
C                    THIS CANNOT OCCUR IF X0 .LE. 0.
C       IERR = -7    ITERATION FAILED. NO VALUE IS GIVEN FOR X.
C                    THIS MAY OCCUR WHEN X IS APPROXIMATELY 0.
C       IERR = -8    A VALUE FOR X HAS BEEN OBTAINED, BUT THE
C                    ROUTINE IS NOT CERTAIN OF ITS ACCURACY.
C                    ITERATION CANNOT BE PERFORMED IN THIS
C                    CASE. IF X0 .LE. 0, THIS CAN OCCUR ONLY
C                    WHEN P OR Q IS APPROXIMATELY 0. IF X0 IS
C                    POSITIVE THEN THIS CAN OCCUR WHEN A IS
C                    EXCEEDINGLY CLOSE TO X AND A IS EXTREMELY
C                    LARGE (SAY A .GE. 1.E20).
C ----------------------------------------------------------------------
C     WRITTEN BY ALFRED H. MORRIS, JR.
C        NAVAL SURFACE WEAPONS CENTER
C        DAHLGREN, VIRGINIA
C     -------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,p,q,x,x0
      INTEGER ierr
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a0,a1,a2,a3,am1,amax,ap1,ap2,ap3,apn,b,b1,b2,b3,
     +                 b4,c,c1,c2,c3,c4,c5,d,e,e2,eps,g,h,ln10,pn,qg,qn,
     +                 r,rta,s,s2,sum,t,tol,u,w,xmax,xmin,xn,y,z
      INTEGER iop
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION amin(2),bmin(2),dmin(2),emin(2),eps0(2)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION alnrel,gamln,gamln1,gamma,rcomp,spmpar
      EXTERNAL alnrel,gamln,gamln1,gamma,rcomp,spmpar
C     ..
C     .. External Subroutines ..
      EXTERNAL gratio
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dble,dlog,dmax1,exp,sqrt
C     ..
C     .. Data statements ..
C     -------------------
C     LN10 = LN(10)
C     C = EULER CONSTANT
C     -------------------
C     -------------------
C     -------------------
C     -------------------
      DATA ln10/2.302585D0/
      DATA c/.577215664901533D0/
      DATA a0/3.31125922108741D0/,a1/11.6616720288968D0/,
     +     a2/4.28342155967104D0/,a3/.213623493715853D0/
      DATA b1/6.61053765625462D0/,b2/6.40691597760039D0/,
     +     b3/1.27364489782223D0/,b4/.036117081018842D0/
      DATA eps0(1)/1.D-10/,eps0(2)/1.D-08/
      DATA amin(1)/500.0D0/,amin(2)/100.0D0/
      DATA bmin(1)/1.D-28/,bmin(2)/1.D-13/
      DATA dmin(1)/1.D-06/,dmin(2)/1.D-04/
      DATA emin(1)/2.D-03/,emin(2)/6.D-03/
      DATA tol/1.D-5/
C     ..
C     .. Executable Statements ..
C     -------------------
C     ****** E, XMIN, AND XMAX ARE MACHINE DEPENDENT CONSTANTS.
C            E IS THE SMALLEST NUMBER FOR WHICH 1.0 + E .GT. 1.0.
C            XMIN IS THE SMALLEST POSITIVE NUMBER AND XMAX IS THE
C            LARGEST POSITIVE NUMBER.
C
      e = spmpar(1)
      xmin = spmpar(2)
      xmax = spmpar(3)
C     -------------------
      x = 0.0D0
      IF (a.LE.0.0D0) GO TO 300
      t = dble(p) + dble(q) - 1.D0
      IF (abs(t).GT.e) GO TO 320
C
      ierr = 0
      IF (p.EQ.0.0D0) RETURN
      IF (q.EQ.0.0D0) GO TO 270
      IF (a.EQ.1.0D0) GO TO 280
C
      e2 = 2.0D0*e
      amax = 0.4D-10/ (e*e)
      iop = 1
      IF (e.GT.1.D-10) iop = 2
      eps = eps0(iop)
      xn = x0
      IF (x0.GT.0.0D0) GO TO 160
C
C        SELECTION OF THE INITIAL APPROXIMATION XN OF X
C                       WHEN A .LT. 1
C
      IF (a.GT.1.0D0) GO TO 80
      g = gamma(a+1.0D0)
      qg = q*g
      IF (qg.EQ.0.0D0) GO TO 360
      b = qg/a
      IF (qg.GT.0.6D0*a) GO TO 40
      IF (a.GE.0.30D0 .OR. b.LT.0.35D0) GO TO 10
      t = exp(- (b+c))
      u = t*exp(t)
      xn = t*exp(u)
      GO TO 160
C
   10 IF (b.GE.0.45D0) GO TO 40
      IF (b.EQ.0.0D0) GO TO 360
      y = -dlog(b)
      s = 0.5D0 + (0.5D0-a)
      z = dlog(y)
      t = y - s*z
      IF (b.LT.0.15D0) GO TO 20
      xn = y - s*dlog(t) - dlog(1.0D0+s/ (t+1.0D0))
      GO TO 220

   20 IF (b.LE.0.01D0) GO TO 30
      u = ((t+2.0D0* (3.0D0-a))*t+ (2.0D0-a)* (3.0D0-a))/
     +    ((t+ (5.0D0-a))*t+2.0D0)
      xn = y - s*dlog(t) - dlog(u)
      GO TO 220

   30 c1 = -s*z
      c2 = -s* (1.0D0+c1)
      c3 = s* ((0.5D0*c1+ (2.0D0-a))*c1+ (2.5D0-1.5D0*a))
      c4 = -s* (((c1/3.0D0+ (2.5D0-1.5D0*a))*c1+ ((a-6.0D0)*a+7.0D0))*
     +     c1+ ((11.0D0*a-46)*a+47.0D0)/6.0D0)
      c5 = -s* ((((-c1/4.0D0+ (11.0D0*a-17.0D0)/6.0D0)*c1+ ((-3.0D0*a+
     +     13.0D0)*a-13.0D0))*c1+0.5D0* (((2.0D0*a-25.0D0)*a+72.0D0)*a-
     +     61.0D0))*c1+ (((25.0D0*a-195.0D0)*a+477.0D0)*a-379.0D0)/
     +     12.0D0)
      xn = ((((c5/y+c4)/y+c3)/y+c2)/y+c1) + y
      IF (a.GT.1.0D0) GO TO 220
      IF (b.GT.bmin(iop)) GO TO 220
      x = xn
      RETURN
C
   40 IF (b*q.GT.1.D-8) GO TO 50
      xn = exp(- (q/a+c))
      GO TO 70

   50 IF (p.LE.0.9D0) GO TO 60
      xn = exp((alnrel(-q)+gamln1(a))/a)
      GO TO 70

   60 xn = exp(dlog(p*g)/a)
   70 IF (xn.EQ.0.0D0) GO TO 310
      t = 0.5D0 + (0.5D0-xn/ (a+1.0D0))
      xn = xn/t
      GO TO 160
C
C        SELECTION OF THE INITIAL APPROXIMATION XN OF X
C                       WHEN A .GT. 1
C
   80 IF (q.LE.0.5D0) GO TO 90
      w = dlog(p)
      GO TO 100

   90 w = dlog(q)
  100 t = sqrt(-2.0D0*w)
      s = t - (((a3*t+a2)*t+a1)*t+a0)/ ((((b4*t+b3)*t+b2)*t+b1)*t+1.0D0)
      IF (q.GT.0.5D0) s = -s
C
      rta = sqrt(a)
      s2 = s*s
      xn = a + s*rta + (s2-1.0D0)/3.0D0 + s* (s2-7.0D0)/ (36.0D0*rta) -
     +     ((3.0D0*s2+7.0D0)*s2-16.0D0)/ (810.0D0*a) +
     +     s* ((9.0D0*s2+256.0D0)*s2-433.0D0)/ (38880.0D0*a*rta)
      xn = dmax1(xn,0.0D0)
      IF (a.LT.amin(iop)) GO TO 110
      x = xn
      d = 0.5D0 + (0.5D0-x/a)
      IF (abs(d).LE.dmin(iop)) RETURN
C
  110 IF (p.LE.0.5D0) GO TO 130
      IF (xn.LT.3.0D0*a) GO TO 220
      y = - (w+gamln(a))
      d = dmax1(2.0D0,a* (a-1.0D0))
      IF (y.LT.ln10*d) GO TO 120
      s = 1.0D0 - a
      z = dlog(y)
      GO TO 30

  120 t = a - 1.0D0
      xn = y + t*dlog(xn) - alnrel(-t/ (xn+1.0D0))
      xn = y + t*dlog(xn) - alnrel(-t/ (xn+1.0D0))
      GO TO 220
C
  130 ap1 = a + 1.0D0
      IF (xn.GT.0.70D0*ap1) GO TO 170
      w = w + gamln(ap1)
      IF (xn.GT.0.15D0*ap1) GO TO 140
      ap2 = a + 2.0D0
      ap3 = a + 3.0D0
      x = exp((w+x)/a)
      x = exp((w+x-dlog(1.0D0+ (x/ap1)* (1.0D0+x/ap2)))/a)
      x = exp((w+x-dlog(1.0D0+ (x/ap1)* (1.0D0+x/ap2)))/a)
      x = exp((w+x-dlog(1.0D0+ (x/ap1)* (1.0D0+ (x/ap2)* (1.0D0+
     +    x/ap3))))/a)
      xn = x
      IF (xn.GT.1.D-2*ap1) GO TO 140
      IF (xn.LE.emin(iop)*ap1) RETURN
      GO TO 170
C
  140 apn = ap1
      t = xn/apn
      sum = 1.0D0 + t
  150 apn = apn + 1.0D0
      t = t* (xn/apn)
      sum = sum + t
      IF (t.GT.1.D-4) GO TO 150
      t = w - dlog(sum)
      xn = exp((xn+t)/a)
      xn = xn* (1.0D0- (a*dlog(xn)-xn-t)/ (a-xn))
      GO TO 170
C
C                 SCHRODER ITERATION USING P
C
  160 IF (p.GT.0.5D0) GO TO 220
  170 IF (p.LE.1.D10*xmin) GO TO 350
      am1 = (a-0.5D0) - 0.5D0
  180 IF (a.LE.amax) GO TO 190
      d = 0.5D0 + (0.5D0-xn/a)
      IF (abs(d).LE.e2) GO TO 350
C
  190 IF (ierr.GE.20) GO TO 330
      ierr = ierr + 1
      CALL gratio(a,xn,pn,qn,0)
      IF (pn.EQ.0.0D0 .OR. qn.EQ.0.0D0) GO TO 350
      r = rcomp(a,xn)
      IF (r.EQ.0.0D0) GO TO 350
      t = (pn-p)/r
      w = 0.5D0* (am1-xn)
      IF (abs(t).LE.0.1D0 .AND. abs(w*t).LE.0.1D0) GO TO 200
      x = xn* (1.0D0-t)
      IF (x.LE.0.0D0) GO TO 340
      d = abs(t)
      GO TO 210
C
  200 h = t* (1.0D0+w*t)
      x = xn* (1.0D0-h)
      IF (x.LE.0.0D0) GO TO 340
      IF (abs(w).GE.1.0D0 .AND. abs(w)*t*t.LE.eps) RETURN
      d = abs(h)
  210 xn = x
      IF (d.GT.tol) GO TO 180
      IF (d.LE.eps) RETURN
      IF (abs(p-pn).LE.tol*p) RETURN
      GO TO 180
C
C                 SCHRODER ITERATION USING Q
C
  220 IF (q.LE.1.D10*xmin) GO TO 350
      am1 = (a-0.5D0) - 0.5D0
  230 IF (a.LE.amax) GO TO 240
      d = 0.5D0 + (0.5D0-xn/a)
      IF (abs(d).LE.e2) GO TO 350
C
  240 IF (ierr.GE.20) GO TO 330
      ierr = ierr + 1
      CALL gratio(a,xn,pn,qn,0)
      IF (pn.EQ.0.0D0 .OR. qn.EQ.0.0D0) GO TO 350
      r = rcomp(a,xn)
      IF (r.EQ.0.0D0) GO TO 350
      t = (q-qn)/r
      w = 0.5D0* (am1-xn)
      IF (abs(t).LE.0.1D0 .AND. abs(w*t).LE.0.1D0) GO TO 250
      x = xn* (1.0D0-t)
      IF (x.LE.0.0D0) GO TO 340
      d = abs(t)
      GO TO 260
C
  250 h = t* (1.0D0+w*t)
      x = xn* (1.0D0-h)
      IF (x.LE.0.0D0) GO TO 340
      IF (abs(w).GE.1.0D0 .AND. abs(w)*t*t.LE.eps) RETURN
      d = abs(h)
  260 xn = x
      IF (d.GT.tol) GO TO 230
      IF (d.LE.eps) RETURN
      IF (abs(q-qn).LE.tol*q) RETURN
      GO TO 230
C
C                       SPECIAL CASES
C
  270 x = xmax
      RETURN
C
  280 IF (q.LT.0.9D0) GO TO 290
      x = -alnrel(-p)
      RETURN

  290 x = -dlog(q)
      RETURN
C
C                       ERROR RETURN
C
  300 ierr = -2
      RETURN
C
  310 ierr = -3
      RETURN
C
  320 ierr = -4
      RETURN
C
  330 ierr = -6
      RETURN
C
  340 ierr = -7
      RETURN
C
  350 x = xn
      ierr = -8
      RETURN
C
  360 x = xmax
      ierr = -8
      RETURN

      END
