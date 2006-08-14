      SUBROUTINE gratio(a,x,ans,qans,ind)
C ----------------------------------------------------------------------
C        EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS
C                      P(A,X) AND Q(A,X)
C
C                        ----------
C
C     IT IS ASSUMED THAT A AND X ARE NONNEGATIVE, WHERE A AND X
C     ARE NOT BOTH 0.
C
C     ANS AND QANS ARE VARIABLES. GRATIO ASSIGNS ANS THE VALUE
C     P(A,X) AND QANS THE VALUE Q(A,X). IND MAY BE ANY INTEGER.
C     IF IND = 0 THEN THE USER IS REQUESTING AS MUCH ACCURACY AS
C     POSSIBLE (UP TO 14 SIGNIFICANT DIGITS). OTHERWISE, IF
C     IND = 1 THEN ACCURACY IS REQUESTED TO WITHIN 1 UNIT OF THE
C     6-TH SIGNIFICANT DIGIT, AND IF IND .NE. 0,1 THEN ACCURACY
C     IS REQUESTED TO WITHIN 1 UNIT OF THE 3RD SIGNIFICANT DIGIT.
C
C     ERROR RETURN ...
C        ANS IS ASSIGNED THE VALUE 2 WHEN A OR X IS NEGATIVE,
C     WHEN A*X = 0, OR WHEN P(A,X) AND Q(A,X) ARE INDETERMINANT.
C     P(A,X) AND Q(A,X) ARE COMPUTATIONALLY INDETERMINANT WHEN
C     X IS EXCEEDINGLY CLOSE TO A AND A IS EXTREMELY LARGE.
C ----------------------------------------------------------------------
C     WRITTEN BY ALFRED H. MORRIS, JR.
C        NAVAL SURFACE WEAPONS CENTER
C        DAHLGREN, VIRGINIA
C     --------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,ans,qans,x
      INTEGER ind
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a2n,a2nm1,acc,alog10,am0,amn,an,an0,apn,b2n,
     +                 b2nm1,c,c0,c1,c2,c3,c4,c5,c6,cma,d10,d20,d30,d40,
     +                 d50,d60,d70,e,e0,g,h,j,l,r,rt2pin,rta,rtpi,rtx,s,
     +                 sum,t,t1,third,tol,twoa,u,w,x0,y,z
      INTEGER i,iop,m,max,n
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION acc0(3),big(3),d0(13),d1(12),d2(10),d3(8),d4(6),
     +                 d5(4),d6(2),e00(3),wk(20),x00(3)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION erf,erfc1,gam1,gamma,rexp,rlog,spmpar
      EXTERNAL erf,erfc1,gam1,gamma,rexp,rlog,spmpar
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dble,dlog,dmax1,exp,int,sqrt
C     ..
C     .. Data statements ..
C     --------------------
C     --------------------
C     ALOG10 = LN(10)
C     RT2PIN = 1/SQRT(2*PI)
C     RTPI   = SQRT(PI)
C     --------------------
C     --------------------
C     --------------------
C     --------------------
C     --------------------
C     --------------------
C     --------------------
C     --------------------
C     --------------------
      DATA acc0(1)/5.D-15/,acc0(2)/5.D-7/,acc0(3)/5.D-4/
      DATA big(1)/20.0D0/,big(2)/14.0D0/,big(3)/10.0D0/
      DATA e00(1)/.25D-3/,e00(2)/.25D-1/,e00(3)/.14D0/
      DATA x00(1)/31.0D0/,x00(2)/17.0D0/,x00(3)/9.7D0/
      DATA alog10/2.30258509299405D0/
      DATA rt2pin/.398942280401433D0/
      DATA rtpi/1.77245385090552D0/
      DATA third/.333333333333333D0/
      DATA d0(1)/.833333333333333D-01/,d0(2)/-.148148148148148D-01/,
     +     d0(3)/.115740740740741D-02/,d0(4)/.352733686067019D-03/,
     +     d0(5)/-.178755144032922D-03/,d0(6)/.391926317852244D-04/,
     +     d0(7)/-.218544851067999D-05/,d0(8)/-.185406221071516D-05/,
     +     d0(9)/.829671134095309D-06/,d0(10)/-.176659527368261D-06/,
     +     d0(11)/.670785354340150D-08/,d0(12)/.102618097842403D-07/,
     +     d0(13)/-.438203601845335D-08/
      DATA d10/-.185185185185185D-02/,d1(1)/-.347222222222222D-02/,
     +     d1(2)/.264550264550265D-02/,d1(3)/-.990226337448560D-03/,
     +     d1(4)/.205761316872428D-03/,d1(5)/-.401877572016461D-06/,
     +     d1(6)/-.180985503344900D-04/,d1(7)/.764916091608111D-05/,
     +     d1(8)/-.161209008945634D-05/,d1(9)/.464712780280743D-08/,
     +     d1(10)/.137863344691572D-06/,d1(11)/-.575254560351770D-07/,
     +     d1(12)/.119516285997781D-07/
      DATA d20/.413359788359788D-02/,d2(1)/-.268132716049383D-02/,
     +     d2(2)/.771604938271605D-03/,d2(3)/.200938786008230D-05/,
     +     d2(4)/-.107366532263652D-03/,d2(5)/.529234488291201D-04/,
     +     d2(6)/-.127606351886187D-04/,d2(7)/.342357873409614D-07/,
     +     d2(8)/.137219573090629D-05/,d2(9)/-.629899213838006D-06/,
     +     d2(10)/.142806142060642D-06/
      DATA d30/.649434156378601D-03/,d3(1)/.229472093621399D-03/,
     +     d3(2)/-.469189494395256D-03/,d3(3)/.267720632062839D-03/,
     +     d3(4)/-.756180167188398D-04/,d3(5)/-.239650511386730D-06/,
     +     d3(6)/.110826541153473D-04/,d3(7)/-.567495282699160D-05/,
     +     d3(8)/.142309007324359D-05/
      DATA d40/-.861888290916712D-03/,d4(1)/.784039221720067D-03/,
     +     d4(2)/-.299072480303190D-03/,d4(3)/-.146384525788434D-05/,
     +     d4(4)/.664149821546512D-04/,d4(5)/-.396836504717943D-04/,
     +     d4(6)/.113757269706784D-04/
      DATA d50/-.336798553366358D-03/,d5(1)/-.697281375836586D-04/,
     +     d5(2)/.277275324495939D-03/,d5(3)/-.199325705161888D-03/,
     +     d5(4)/.679778047793721D-04/
      DATA d60/.531307936463992D-03/,d6(1)/-.592166437353694D-03/,
     +     d6(2)/.270878209671804D-03/
      DATA d70/.344367606892378D-03/
C     ..
C     .. Executable Statements ..
C     --------------------
C     ****** E IS A MACHINE DEPENDENT CONSTANT. E IS THE SMALLEST
C            FLOATING POINT NUMBER FOR WHICH 1.0 + E .GT. 1.0 .
C
      e = spmpar(1)
C
C     --------------------
      IF (a.LT.0.0D0 .OR. x.LT.0.0D0) GO TO 430
      IF (a.EQ.0.0D0 .AND. x.EQ.0.0D0) GO TO 430
      IF (a*x.EQ.0.0D0) GO TO 420
C
      iop = ind + 1
      IF (iop.NE.1 .AND. iop.NE.2) iop = 3
      acc = dmax1(acc0(iop),e)
      e0 = e00(iop)
      x0 = x00(iop)
C
C            SELECT THE APPROPRIATE ALGORITHM
C
      IF (a.GE.1.0D0) GO TO 10
      IF (a.EQ.0.5D0) GO TO 390
      IF (x.LT.1.1D0) GO TO 160
      t1 = a*dlog(x) - x
      u = a*exp(t1)
      IF (u.EQ.0.0D0) GO TO 380
      r = u* (1.0D0+gam1(a))
      GO TO 250
C
   10 IF (a.GE.big(iop)) GO TO 30
      IF (a.GT.x .OR. x.GE.x0) GO TO 20
      twoa = a + a
      m = int(twoa)
      IF (twoa.NE.dble(m)) GO TO 20
      i = m/2
      IF (a.EQ.dble(i)) GO TO 210
      GO TO 220

   20 t1 = a*dlog(x) - x
      r = exp(t1)/gamma(a)
      GO TO 40
C
   30 l = x/a
      IF (l.EQ.0.0D0) GO TO 370
      s = 0.5D0 + (0.5D0-l)
      z = rlog(l)
      IF (z.GE.700.0D0/a) GO TO 410
      y = a*z
      rta = sqrt(a)
      IF (abs(s).LE.e0/rta) GO TO 330
      IF (abs(s).LE.0.4D0) GO TO 270
C
      t = (1.0D0/a)**2
      t1 = (((0.75D0*t-1.0D0)*t+3.5D0)*t-105.0D0)/ (a*1260.0D0)
      t1 = t1 - y
      r = rt2pin*rta*exp(t1)
C
   40 IF (r.EQ.0.0D0) GO TO 420
      IF (x.LE.dmax1(a,alog10)) GO TO 50
      IF (x.LT.x0) GO TO 250
      GO TO 100
C
C                 TAYLOR SERIES FOR P/R
C
   50 apn = a + 1.0D0
      t = x/apn
      wk(1) = t
      DO 60 n = 2,20
          apn = apn + 1.0D0
          t = t* (x/apn)
          IF (t.LE.1.D-3) GO TO 70
          wk(n) = t
   60 CONTINUE
      n = 20
C
   70 sum = t
      tol = 0.5D0*acc
   80 apn = apn + 1.0D0
      t = t* (x/apn)
      sum = sum + t
      IF (t.GT.tol) GO TO 80
C
      max = n - 1
      DO 90 m = 1,max
          n = n - 1
          sum = sum + wk(n)
   90 CONTINUE
      ans = (r/a)* (1.0D0+sum)
      qans = 0.5D0 + (0.5D0-ans)
      RETURN
C
C                 ASYMPTOTIC EXPANSION
C
  100 amn = a - 1.0D0
      t = amn/x
      wk(1) = t
      DO 110 n = 2,20
          amn = amn - 1.0D0
          t = t* (amn/x)
          IF (abs(t).LE.1.D-3) GO TO 120
          wk(n) = t
  110 CONTINUE
      n = 20
C
  120 sum = t
  130 IF (abs(t).LE.acc) GO TO 140
      amn = amn - 1.0D0
      t = t* (amn/x)
      sum = sum + t
      GO TO 130
C
  140 max = n - 1
      DO 150 m = 1,max
          n = n - 1
          sum = sum + wk(n)
  150 CONTINUE
      qans = (r/x)* (1.0D0+sum)
      ans = 0.5D0 + (0.5D0-qans)
      RETURN
C
C             TAYLOR SERIES FOR P(A,X)/X**A
C
  160 an = 3.0D0
      c = x
      sum = x/ (a+3.0D0)
      tol = 3.0D0*acc/ (a+1.0D0)
  170 an = an + 1.0D0
      c = -c* (x/an)
      t = c/ (a+an)
      sum = sum + t
      IF (abs(t).GT.tol) GO TO 170
      j = a*x* ((sum/6.0D0-0.5D0/ (a+2.0D0))*x+1.0D0/ (a+1.0D0))
C
      z = a*dlog(x)
      h = gam1(a)
      g = 1.0D0 + h
      IF (x.LT.0.25D0) GO TO 180
      IF (a.LT.x/2.59D0) GO TO 200
      GO TO 190

  180 IF (z.GT.-.13394D0) GO TO 200
C
  190 w = exp(z)
      ans = w*g* (0.5D0+ (0.5D0-j))
      qans = 0.5D0 + (0.5D0-ans)
      RETURN
C
  200 l = rexp(z)
      w = 0.5D0 + (0.5D0+l)
      qans = (w*j-l)*g - h
      IF (qans.LT.0.0D0) GO TO 380
      ans = 0.5D0 + (0.5D0-qans)
      RETURN
C
C             FINITE SUMS FOR Q WHEN A .GE. 1
C                 AND 2*A IS AN INTEGER
C
  210 sum = exp(-x)
      t = sum
      n = 1
      c = 0.0D0
      GO TO 230
C
  220 rtx = sqrt(x)
      sum = erfc1(0,rtx)
      t = exp(-x)/ (rtpi*rtx)
      n = 0
      c = -0.5D0
C
  230 IF (n.EQ.i) GO TO 240
      n = n + 1
      c = c + 1.0D0
      t = (x*t)/c
      sum = sum + t
      GO TO 230

  240 qans = sum
      ans = 0.5D0 + (0.5D0-qans)
      RETURN
C
C              CONTINUED FRACTION EXPANSION
C
  250 tol = dmax1(5.0D0*e,acc)
      a2nm1 = 1.0D0
      a2n = 1.0D0
      b2nm1 = x
      b2n = x + (1.0D0-a)
      c = 1.0D0
  260 a2nm1 = x*a2n + c*a2nm1
      b2nm1 = x*b2n + c*b2nm1
      am0 = a2nm1/b2nm1
      c = c + 1.0D0
      cma = c - a
      a2n = a2nm1 + cma*a2n
      b2n = b2nm1 + cma*b2n
      an0 = a2n/b2n
      IF (abs(an0-am0).GE.tol*an0) GO TO 260
C
      qans = r*an0
      ans = 0.5D0 + (0.5D0-qans)
      RETURN
C
C                GENERAL TEMME EXPANSION
C
  270 IF (abs(s).LE.2.0D0*e .AND. a*e*e.GT.3.28D-3) GO TO 430
      c = exp(-y)
      w = 0.5D0*erfc1(1,sqrt(y))
      u = 1.0D0/a
      z = sqrt(z+z)
      IF (l.LT.1.0D0) z = -z
      IF (iop.lt.2) GO TO 280
      IF (iop.eq.2) GO TO 290
      GO TO 300
C
  280 IF (abs(s).LE.1.D-3) GO TO 340
      c0 = ((((((((((((d0(13)*z+d0(12))*z+d0(11))*z+d0(10))*z+d0(9))*z+
     +     d0(8))*z+d0(7))*z+d0(6))*z+d0(5))*z+d0(4))*z+d0(3))*z+d0(2))*
     +     z+d0(1))*z - third
      c1 = (((((((((((d1(12)*z+d1(11))*z+d1(10))*z+d1(9))*z+d1(8))*z+
     +     d1(7))*z+d1(6))*z+d1(5))*z+d1(4))*z+d1(3))*z+d1(2))*z+d1(1))*
     +     z + d10
      c2 = (((((((((d2(10)*z+d2(9))*z+d2(8))*z+d2(7))*z+d2(6))*z+
     +     d2(5))*z+d2(4))*z+d2(3))*z+d2(2))*z+d2(1))*z + d20
      c3 = (((((((d3(8)*z+d3(7))*z+d3(6))*z+d3(5))*z+d3(4))*z+d3(3))*z+
     +     d3(2))*z+d3(1))*z + d30
      c4 = (((((d4(6)*z+d4(5))*z+d4(4))*z+d4(3))*z+d4(2))*z+d4(1))*z +
     +     d40
      c5 = (((d5(4)*z+d5(3))*z+d5(2))*z+d5(1))*z + d50
      c6 = (d6(2)*z+d6(1))*z + d60
      t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u + c0
      GO TO 310
C
  290 c0 = (((((d0(6)*z+d0(5))*z+d0(4))*z+d0(3))*z+d0(2))*z+d0(1))*z -
     +     third
      c1 = (((d1(4)*z+d1(3))*z+d1(2))*z+d1(1))*z + d10
      c2 = d2(1)*z + d20
      t = (c2*u+c1)*u + c0
      GO TO 310
C
  300 t = ((d0(3)*z+d0(2))*z+d0(1))*z - third
C
  310 IF (l.LT.1.0D0) GO TO 320
      qans = c* (w+rt2pin*t/rta)
      ans = 0.5D0 + (0.5D0-qans)
      RETURN

  320 ans = c* (w-rt2pin*t/rta)
      qans = 0.5D0 + (0.5D0-ans)
      RETURN
C
C               TEMME EXPANSION FOR L = 1
C
  330 IF (a*e*e.GT.3.28D-3) GO TO 430
      c = 0.5D0 + (0.5D0-y)
      w = (0.5D0-sqrt(y)* (0.5D0+ (0.5D0-y/3.0D0))/rtpi)/c
      u = 1.0D0/a
      z = sqrt(z+z)
      IF (l.LT.1.0D0) z = -z
      IF (iop.lt.2) GO TO 340
      IF (iop.eq.2) GO TO 350
      GO TO 360
C
  340 c0 = ((((((d0(7)*z+d0(6))*z+d0(5))*z+d0(4))*z+d0(3))*z+d0(2))*z+
     +     d0(1))*z - third
      c1 = (((((d1(6)*z+d1(5))*z+d1(4))*z+d1(3))*z+d1(2))*z+d1(1))*z +
     +     d10
      c2 = ((((d2(5)*z+d2(4))*z+d2(3))*z+d2(2))*z+d2(1))*z + d20
      c3 = (((d3(4)*z+d3(3))*z+d3(2))*z+d3(1))*z + d30
      c4 = (d4(2)*z+d4(1))*z + d40
      c5 = (d5(2)*z+d5(1))*z + d50
      c6 = d6(1)*z + d60
      t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u + c0
      GO TO 310
C
  350 c0 = (d0(2)*z+d0(1))*z - third
      c1 = d1(1)*z + d10
      t = (d20*u+c1)*u + c0
      GO TO 310
C
  360 t = d0(1)*z - third
      GO TO 310
C
C                     SPECIAL CASES
C
  370 ans = 0.0D0
      qans = 1.0D0
      RETURN
C
  380 ans = 1.0D0
      qans = 0.0D0
      RETURN
C
  390 IF (x.GE.0.25D0) GO TO 400
      ans = erf(sqrt(x))
      qans = 0.5D0 + (0.5D0-ans)
      RETURN

  400 qans = erfc1(0,sqrt(x))
      ans = 0.5D0 + (0.5D0-qans)
      RETURN
C
  410 IF (abs(s).LE.2.0D0*e) GO TO 430
  420 IF (x.LE.a) GO TO 370
      GO TO 380
C
C                     ERROR RETURN
C
  430 ans = 2.0D0
      RETURN

      END
