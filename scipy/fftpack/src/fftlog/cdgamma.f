c-----------------------------------------------------------------------
c This code was modified from a subroutine from the gamerf package at
c http://momonga.t.u-tokyo.ac.jp/~ooura/gamerf.html .
c The original copyright statement states:
c   Copyright(C) 1996 Takuya OOURA (email: ooura@mmm.t.u-tokyo.ac.jp).
c   You may use, copy, modify this code for any purpose and 
c   without fee. You may distribute this ORIGINAL package.
c
c Permission to distribute this modified gamma function code
c with the FFTLog package has been granted
c (email from Takuya Ooura to Andrew Hamilton dated 16 March 1999).
c
c Original gamerf2a.doc documentation states:
c
c   Gamma(z)=sqrt(2*pi)*(z+r)^(z-1/2)*exp(-z-r)
c            *(a_0
c            +a_1*(z-1)/z
c            +a_2*(z-1)*(z-2)/z/(z+1)
c            +a_3*(z-1)*(z-2)*(z-3)/z/(z+1)/(z+2)
c            +...)
c        a_n= f_n*(2*n)*(2*n-1)*(2*n-2)*...*(n+1)/1/2/3/.../n
c            -a_0*(2*n)*(2*n-1)*(2*n-2)*...*(n+1)/1/2/3/.../n
c            -a_1*(2*n)*(2*n-1)*(2*n-2)*...*(n+2)/1/2/3/.../(n-1)
c            -...
c            -a_(n-1)*(2*n)/1
c        f_n=1/sqrt(2*pi)*(1*2*3*...*n)*(n+1+r)^(-n-1/2)*exp(n+1+r)

c   C.Lanczos,A Precision Approximation of the Gamma Function,
c   J.SIAM Numer.Anal.Ser.B,Vol.1,1964
c
c Modified 28 Oct 98 by Andrew J S Hamilton
c http://casa.colorado.edu/~ajsh/
c (1) to return ln[Gamma(x)] with the correct phase,
c     as well as Gamma(x), and
c (2) to remain accurate for large absolute values of input x,
c     and for x near 0 or negative integers.
c-----------------------------------------------------------------------
      complex*16 function cdgamma(x,l)
      integer l
      complex*16 x
c *
c * Complex Gamma function in double precision.
c *
c * l = 0: Gamma(x)
c *     1: ln[Gamma(x)]
c *
      real*8 pi,pv,pu,pr,p1,p2,p3,p4,p5,p6,q1,q2,q3,q4,q5,q6
      parameter (
     &    pi = 3.14159265358979324d+00, 
     &    pv = 7.31790632447016203d+00, 
     &    pu = 3.48064577727581257d+00, 
     &    pr = 3.27673720261526849d-02, 
     &    p1 = 1.05400280458730808d+01, 
     &    p2 = 4.73821439163096063d+01, 
     &    p3 = 9.11395751189899762d+01, 
     &    p4 = 6.62756400966213521d+01, 
     &    p5 = 1.32280130755055088d+01, 
     &    p6 = 2.93729529320536228d-01)
      parameter (
     &    q1 = 9.99999999999975753d-01, 
     &    q2 = 2.00000000000603851d+00, 
     &    q3 = 2.99999999944915534d+00, 
     &    q4 = 4.00000003016801681d+00, 
     &    q5 = 4.99999857982434025d+00, 
     &    q6 = 6.00009857740312429d+00)
      real*8 big
      parameter (big=1.d20)
      real*8 t,ui,ur,vi,vr,wr,wi,xi,xr,yi,yr,zero
      data zero /0.d0/
c
      xr = dble(x)
      xi = dimag(x)
c---x = 0, -1, -2
      if (xr .eq. aint(xr) .and. xr .le. 0.d0 .and. xi .eq. 0.d0) then
c...Gamma
        if (l .eq. 0) then
          wr = xr / 2.d0
c   +Infinity at even negative integers
          if (wr .eq. aint(wr)) then
            yr = 1.d0 / zero
c   -Infinity at odd negative integers
          else
            yr = -1.d0 / zero
          endif
          yi = 0.d0
c...lnGamma
        elseif (l .eq. 1) then
c   real part is +Infinity
          yr = 1.d0 / zero
          yi = pi * aint(xr)
        endif
        goto 200
      endif
c---Re(x) < 1/2 : use reflection formula
      if (xr .lt. .5d0) then
          wr = 1.d0 - xr
          wi = -xi
      else
          wr = xr
          wi = xi
      end if
c---large |x|: keep only leading term of rational function
      t = wr * wr + wi * wi
      if (t .gt. big * big) then
c   Rational function v
        vr = wr / t + pr
        vi = - wi / t
c   ln(overall factor)
c   u = ln(x + pv) - 1
        yr = wr + pv
        ur = yr
        if (ur .lt. 0.d0) ur = - ur
        ui = wi
        if (ui .lt. 0.d0) ui = - ui
        if (ur.ge.ui) then
          t = wi / yr
          ur = log(ur) + log(1.d0 + t * t) / 2.d0 - 1.d0
        else
          t = yr / wi
          ur = log(ui) + log(1.d0 + t * t) / 2.d0 - 1.d0
        endif
        ui = atan2(wi, yr)
c---not large |x|
      else
c   u = u(x) = x + q6 = O(x)
        ur = wr + q6
c   v = v(x,u) = (x + q5) u = (x + q5)(x + q6) = O(x^2)
        vr = ur * (wr + q5) - wi * wi
        vi = wi * (wr + q5) + ur * wi
c   y = y(x,u,v) = p6 + p5 u + p4 v = p4 x^2 + ... = O(x^2)
        yr = p6 + (p5 * ur + p4 * vr)
        yi = p5 * wi + p4 * vi
c   u = u(x,v) = (x + q4) v = (x + q4)(x + q5)(x + q6) = O(x^3)
        ur = vr * (wr + q4) - vi * wi
        ui = vi * (wr + q4) + vr * wi
c   v = v(x,u) = (x + q3) u = (x + q3)(x + q4)...(x + q6) = O(x^4)
        vr = ur * (wr + q3) - ui * wi
        vi = ui * (wr + q3) + ur * wi
c   y = y(y,u,v) = y + p3 u + p2 v = p2 x^4 + ... = O(x^4)
        yr = yr + (p3 * ur + p2 * vr)
        yi = yi + (p3 * ui + p2 * vi)
c   u = u(x,v) = (x + q2) v = (x + q2)(x + q3)...(x + q6) = O(x^5)
        ur = vr * (wr + q2) - vi * wi
        ui = vi * (wr + q2) + vr * wi
c   v = v(x,u) = (x + q1) u = (x + q1)(x + q2)...(x + q6) = O(x^6)
        vr = ur * (wr + q1) - ui * wi
        vi = ui * (wr + q1) + ur * wi
c   Numerator
c   y = y(y,u,v) = y + p1 u + v = x^6 + ... = O(x^6)
c     = (x+q1)...(x+q6) + p1 (x+q2)...(x+q6) + ... + p5 (x+q6) + p6
        yr = yr + (p1 * ur + vr)
        yi = yi + (p1 * ui + vi)
c   Denominator
c   u = x v = x(x + q1)(x + q2)...(x + q6) = O(x^7)
        ur = vr * wr - vi * wi
        ui = vi * wr + vr * wi
c   t = |u|^2
        t = ur * ur + ui * ui
c   Rational function v = y u*/|u|^2 + pr = y/u + pr = pr + 1/x + ...
c     = pr + 1/x ( 1 + 1/(x+q1) ( p1 + 1/(x+q2) ( p2 + ...
        vr = (yr * ur + yi * ui) / t + pr
        vi = (yi * ur - yr * ui) / t
c   Overall factor
c   u = ln(x + pv) - 1
        yr = wr + pv
        ur = log(yr * yr + wi * wi) / 2.d0 - 1.d0
        ui = atan2(wi, yr)
      endif
c---lnGamma
c   y = u(x - .5) - pu
c     = (x - .5) [ln(x + pv) - 1] - pu
      yr = ur * (wr - 0.5d0) - ui * wi - pu
      yi = ui * (wr - 0.5d0) + ur * wi
c   y = y + ln(v)
c     = (x - .5) [ln(x + pv) - 1] - pu + ln(Rational)
c     = lnGamma(x)
      yr = yr + log(vr * vr + vi * vi) / 2.d0
      yi = yi + atan2(vi, vr)
c---Reflection formula Gamma(x) Gamma(1-x) = pi/sin(pi x)
c   sign of Gamma
      t = 1.d0
      if (xr .lt. .5d0) then
        wi = anint(xr)
        wr = xr - wi
        if (wi .gt. xr) wi = wi - 1.d0
c   case of real x
        if (xi .eq. 0.d0) then
c   w = ln[sin(pi x)]
          wr = log(sin(pi * abs(wr)))
          if (l .eq. 0) then
            if (wi .ne. 2.d0 * aint(wi / 2)) t = -1.d0
            wi = 0.d0
          elseif (l .eq. 1) then
            wi = - pi * wi
          endif
c   case where imaginary part of x is < 1 in absolute value
        elseif (abs(xi) .lt. 1.d0) then
          if (l .eq. 0) then
            if (wi .ne. 2.d0 * aint(wi / 2.d0)) t = -1.d0
            ui = 0.d0
          elseif (l .eq. 1) then
            ui = -pi * wi
            if (xi .lt. 0.d0) ui = -ui
          endif
          wr = pi * wr
          wi = pi * xi
          vr = sin(wr) * cosh(wi)
          vi = cos(wr) * sinh(wi)
          if (wr .lt. 0.d0) then
            vr = -vr
            vi = -vi
          endif
c   w = ln[sin(pi x)]
          wr = log(vr * vr + vi * vi) / 2.d0
          wi = ui + atan2(vi, vr)
c   case where imaginary part of x is >= 1 in absolute value
        else
          if (l .eq. 0) then
            if (wi .ne. 2.d0 * aint(wi / 2)) t = -1.d0
            if (wr .ge. 0.d0) then
              ui = pi * (.5d0 - wr)
            else
              ui = pi * (- .5d0 - wr)
            endif
          elseif (l .eq. 1) then
            ui = pi * (.5d0 - xr)
          endif
          wi = exp(- 2.d0 * pi * abs(xi))
          wr = 2.d0 * pi * wr
          vr = (1.d0 - cos(wr) * wi) / 2.d0
          vi = - sin(wr) * wi / 2.d0
          ur = pi * xi
c   w = ln[sin(pi x)]
          if (xi .gt. 0.d0) then
            wr = ur + log(vr * vr + vi * vi) / 2.d0
            wi = ui + atan2(vi, vr)
          elseif (xi .lt. 0.d0) then
            wr = - ur + log(vr * vr + vi * vi) / 2.d0
            wi = - ui - atan2(vi, vr)
          endif
        endif
c   y = ln[Gamma(x)]
        yr = log(pi) - yr - wr
        yi = - yi - wi
      endif
c---Gamma
      if (l .eq. 0) then
        ur = exp(yr)
        if (xi .eq. 0.d0) then
          yr = t * ur
          yi = 0.d0
        else
          yr = t * ur * cos(yi)
          yi = ur * sin(yi)
        endif
      endif
c---finish
  200 cdgamma = dcmplx(yr, yi)
      end
c
