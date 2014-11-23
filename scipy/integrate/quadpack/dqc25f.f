      subroutine dqc25f(f,a,b,omega,integr,nrmom,maxp1,ksave,result,
     *   abserr,neval,resabs,resasc,momcom,chebmo)
c***begin prologue  dqc25f
c***date written   810101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a2
c***keywords  integration rules for functions with cos or sin
c             factor, clenshaw-curtis, gauss-kronrod
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute the integral i=integral of f(x) over (a,b)
c            where w(x) = cos(omega*x) or w(x)=sin(omega*x) and to
c            compute j = integral of abs(f) over (a,b). for small value
c            of omega or small intervals (a,b) the 15-point gauss-kronro
c            rule is used. otherwise a generalized clenshaw-curtis
c            method is used.
c***description
c
c        integration rules for functions with cos or sin factor
c        standard fortran subroutine
c        double precision version
c
c        parameters
c         on entry
c           f      - double precision
c                    function subprogram defining the integrand
c                    function f(x). the actual name for f needs to
c                    be declared e x t e r n a l in the calling program.
c
c           a      - double precision
c                    lower limit of integration
c
c           b      - double precision
c                    upper limit of integration
c
c           omega  - double precision
c                    parameter in the weight function
c
c           integr - integer
c                    indicates which weight function is to be used
c                       integr = 1   w(x) = cos(omega*x)
c                       integr = 2   w(x) = sin(omega*x)
c
c           nrmom  - integer
c                    the length of interval (a,b) is equal to the length
c                    of the original integration interval divided by
c                    2**nrmom (we suppose that the routine is used in an
c                    adaptive integration process, otherwise set
c                    nrmom = 0). nrmom must be zero at the first call.
c
c           maxp1  - integer
c                    gives an upper bound on the number of chebyshev
c                    moments which can be stored, i.e. for the
c                    intervals of lengths abs(bb-aa)*2**(-l),
c                    l = 0,1,2, ..., maxp1-2.
c
c           ksave  - integer
c                    key which is one when the moments for the
c                    current interval have been computed
c
c         on return
c           result - double precision
c                    approximation to the integral i
c
c           abserr - double precision
c                    estimate of the modulus of the absolute
c                    error, which should equal or exceed abs(i-result)
c
c           neval  - integer
c                    number of integrand evaluations
c
c           resabs - double precision
c                    approximation to the integral j
c
c           resasc - double precision
c                    approximation to the integral of abs(f-i/(b-a))
c
c         on entry and return
c           momcom - integer
c                    for each interval length we need to compute the
c                    chebyshev moments. momcom counts the number of
c                    intervals for which these moments have already been
c                    computed. if nrmom.lt.momcom or ksave = 1, the
c                    chebyshev moments for the interval (a,b) have
c                    already been computed and stored, otherwise we
c                    compute them and we increase momcom.
c
c           chebmo - double precision
c                    array of dimension at least (maxp1,25) containing
c                    the modified chebyshev moments for the first momcom
c                    momcom interval lengths
c
c ......................................................................
c***references  (none)
c***routines called  d1mach,dgtsl,dqcheb,dqk15w,dqwgtf
c***end prologue  dqc25f
c
      double precision a,abserr,ac,an,an2,as,asap,ass,b,centr,chebmo,
     *  cheb12,cheb24,conc,cons,cospar,d,dabs,dcos,dsin,dqwgtf,d1,
     *  d1mach,d2,estc,ests,f,fval,hlgth,oflow,omega,parint,par2,par22,
     *  p2,p3,p4,resabs,resasc,resc12,resc24,ress12,ress24,result,
     *  sinpar,v,x
      integer i,iers,integr,isym,j,k,ksave,m,momcom,neval,maxp1,
     *  noequ,noeq1,nrmom
c
      dimension chebmo(maxp1,25),cheb12(13),cheb24(25),d(25),d1(25),
     *  d2(25),fval(25),v(28),x(11)
c
      external f,dqwgtf
c
c           the vector x contains the values cos(k*pi/24)
c           k = 1, ...,11, to be used for the chebyshev expansion of f
c
      data x(1) / 0.9914448613 7381041114 4557526928 563d0 /
      data x(2) / 0.9659258262 8906828674 9743199728 897d0 /
      data x(3) / 0.9238795325 1128675612 8183189396 788d0 /
      data x(4) / 0.8660254037 8443864676 3723170752 936d0 /
      data x(5) / 0.7933533402 9123516457 9776961501 299d0 /
      data x(6) / 0.7071067811 8654752440 0844362104 849d0 /
      data x(7) / 0.6087614290 0872063941 6097542898 164d0 /
      data x(8) / 0.5000000000 0000000000 0000000000 000d0 /
      data x(9) / 0.3826834323 6508977172 8459984030 399d0 /
      data x(10) / 0.2588190451 0252076234 8898837624 048d0 /
      data x(11) / 0.1305261922 2005159154 8406227895 489d0 /
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the integration interval
c           hlgth  - half-length of the integration interval
c           fval   - value of the function f at the points
c                    (b-a)*0.5*cos(k*pi/12) + (b+a)*0.5, k = 0, ..., 24
c           cheb12 - coefficients of the chebyshev series expansion
c                    of degree 12, for the function f, in the
c                    interval (a,b)
c           cheb24 - coefficients of the chebyshev series expansion
c                    of degree 24, for the function f, in the
c                    interval (a,b)
c           resc12 - approximation to the integral of
c                    cos(0.5*(b-a)*omega*x)*f(0.5*(b-a)*x+0.5*(b+a))
c                    over (-1,+1), using the chebyshev series
c                    expansion of degree 12
c           resc24 - approximation to the same integral, using the
c                    chebyshev series expansion of degree 24
c           ress12 - the analogue of resc12 for the sine
c           ress24 - the analogue of resc24 for the sine
c
c
c           machine dependent constant
c           --------------------------
c
c           oflow is the largest positive magnitude.
c
c***first executable statement  dqc25f
      oflow = d1mach(2)
c
      centr = 0.5d+00*(b+a)
      hlgth = 0.5d+00*(b-a)
      parint = omega*hlgth
c
c           compute the integral using the 15-point gauss-kronrod
c           formula if the value of the parameter in the integrand
c           is small.
c
      if(dabs(parint).gt.0.2d+01) go to 10
      call dqk15w(f,dqwgtf,omega,p2,p3,p4,integr,a,b,result,
     *  abserr,resabs,resasc)
      neval = 15
      go to 170
c
c           compute the integral using the generalized clenshaw-
c           curtis method.
c
   10 conc = hlgth*dcos(centr*omega)
      cons = hlgth*dsin(centr*omega)
      resasc = oflow
      neval = 25
c
c           check whether the chebyshev moments for this interval
c           have already been computed.
c
      if(nrmom.lt.momcom.or.ksave.eq.1) go to 120
c
c           compute a new set of chebyshev moments.
c
      m = momcom+1
      par2 = parint*parint
      par22 = par2+0.2d+01
      sinpar = dsin(parint)
      cospar = dcos(parint)
c
c           compute the chebyshev moments with respect to cosine.
c
      v(1) = 0.2d+01*sinpar/parint
      v(2) = (0.8d+01*cospar+(par2+par2-0.8d+01)*sinpar/parint)/par2
      v(3) = (0.32d+02*(par2-0.12d+02)*cospar+(0.2d+01*
     *  ((par2-0.80d+02)*par2+0.192d+03)*sinpar)/parint)/(par2*par2)
      ac = 0.8d+01*cospar
      as = 0.24d+02*parint*sinpar
      if(dabs(parint).gt.0.24d+02) go to 30
c
c           compute the chebyshev moments as the solutions of a
c           boundary value problem with 1 initial value (v(3)) and 1
c           end value (computed using an asymptotic formula).
c
      noequ = 25
      noeq1 = noequ-1
      an = 0.6d+01
      do 20 k = 1,noeq1
        an2 = an*an
        d(k) = -0.2d+01*(an2-0.4d+01)*(par22-an2-an2)
        d2(k) = (an-0.1d+01)*(an-0.2d+01)*par2
        d1(k+1) = (an+0.3d+01)*(an+0.4d+01)*par2
        v(k+3) = as-(an2-0.4d+01)*ac
        an = an+0.2d+01
   20 continue
      an2 = an*an
      d(noequ) = -0.2d+01*(an2-0.4d+01)*(par22-an2-an2)
      v(noequ+3) = as-(an2-0.4d+01)*ac
      v(4) = v(4)-0.56d+02*par2*v(3)
      ass = parint*sinpar
      asap = (((((0.210d+03*par2-0.1d+01)*cospar-(0.105d+03*par2
     *  -0.63d+02)*ass)/an2-(0.1d+01-0.15d+02*par2)*cospar
     *  +0.15d+02*ass)/an2-cospar+0.3d+01*ass)/an2-cospar)/an2
      v(noequ+3) = v(noequ+3)-0.2d+01*asap*par2*(an-0.1d+01)*
     *   (an-0.2d+01)
c
c           solve the tridiagonal system by means of gaussian
c           elimination with partial pivoting.
c
c***        call to dgtsl has been replaced by call to
c***        lapack routine dgtsv
c
c      call dgtsl(noequ,d1,d,d2,v(4),iers)
      call dgtsv(noequ,1,d1(2),d,d2,v(4),noequ,iers)
      go to 50
c
c           compute the chebyshev moments by means of forward
c           recursion.
c
   30 an = 0.4d+01
      do 40 i = 4,13
        an2 = an*an
        v(i) = ((an2-0.4d+01)*(0.2d+01*(par22-an2-an2)*v(i-1)-ac)
     *  +as-par2*(an+0.1d+01)*(an+0.2d+01)*v(i-2))/
     *  (par2*(an-0.1d+01)*(an-0.2d+01))
        an = an+0.2d+01
   40 continue
   50 do 60 j = 1,13
        chebmo(m,2*j-1) = v(j)
   60 continue
c
c           compute the chebyshev moments with respect to sine.
c
      v(1) = 0.2d+01*(sinpar-parint*cospar)/par2
      v(2) = (0.18d+02-0.48d+02/par2)*sinpar/par2
     *  +(-0.2d+01+0.48d+02/par2)*cospar/parint
      ac = -0.24d+02*parint*cospar
      as = -0.8d+01*sinpar
      if(dabs(parint).gt.0.24d+02) go to 80
c
c           compute the chebyshev moments as the solutions of a boundary
c           value problem with 1 initial value (v(2)) and 1 end value
c           (computed using an asymptotic formula).
c
      an = 0.5d+01
      do 70 k = 1,noeq1
        an2 = an*an
        d(k) = -0.2d+01*(an2-0.4d+01)*(par22-an2-an2)
        d2(k) = (an-0.1d+01)*(an-0.2d+01)*par2
        d1(k+1) = (an+0.3d+01)*(an+0.4d+01)*par2
        v(k+2) = ac+(an2-0.4d+01)*as
        an = an+0.2d+01
   70 continue
      an2 = an*an
      d(noequ) = -0.2d+01*(an2-0.4d+01)*(par22-an2-an2)
      v(noequ+2) = ac+(an2-0.4d+01)*as
      v(3) = v(3)-0.42d+02*par2*v(2)
      ass = parint*cospar
      asap = (((((0.105d+03*par2-0.63d+02)*ass+(0.210d+03*par2
     *  -0.1d+01)*sinpar)/an2+(0.15d+02*par2-0.1d+01)*sinpar-
     *  0.15d+02*ass)/an2-0.3d+01*ass-sinpar)/an2-sinpar)/an2
      v(noequ+2) = v(noequ+2)-0.2d+01*asap*par2*(an-0.1d+01)
     *  *(an-0.2d+01)
c
c           solve the tridiagonal system by means of gaussian
c           elimination with partial pivoting.
c
c***        call to dgtsl has been replaced by call to
c***        lapack routine dgtsv
c
c      call dgtsl(noequ,d1,d,d2,v(3),iers)
      call dgtsv(noequ,1,d1(2),d,d2,v(3),noequ,iers)
      go to 100
c
c           compute the chebyshev moments by means of forward recursion.
c
   80 an = 0.3d+01
      do 90 i = 3,12
        an2 = an*an
        v(i) = ((an2-0.4d+01)*(0.2d+01*(par22-an2-an2)*v(i-1)+as)
     *  +ac-par2*(an+0.1d+01)*(an+0.2d+01)*v(i-2))
     *  /(par2*(an-0.1d+01)*(an-0.2d+01))
        an = an+0.2d+01
   90 continue
  100 do 110 j = 1,12
        chebmo(m,2*j) = v(j)
  110 continue
  120 if (nrmom.lt.momcom) m = nrmom+1
       if (momcom.lt.(maxp1-1).and.nrmom.ge.momcom) momcom = momcom+1
c
c           compute the coefficients of the chebyshev expansions
c           of degrees 12 and 24 of the function f.
c
      fval(1) = 0.5d+00*f(centr+hlgth)
      fval(13) = f(centr)
      fval(25) = 0.5d+00*f(centr-hlgth)
      do 130 i = 2,12
        isym = 26-i
        fval(i) = f(hlgth*x(i-1)+centr)
        fval(isym) = f(centr-hlgth*x(i-1))
  130 continue
      call dqcheb(x,fval,cheb12,cheb24)
c
c           compute the integral and error estimates.
c
      resc12 = cheb12(13)*chebmo(m,13)
      ress12 = 0.0d+00
      k = 11
      do 140 j = 1,6
        resc12 = resc12+cheb12(k)*chebmo(m,k)
        ress12 = ress12+cheb12(k+1)*chebmo(m,k+1)
        k = k-2
  140 continue
      resc24 = cheb24(25)*chebmo(m,25)
      ress24 = 0.0d+00
      resabs = dabs(cheb24(25))
      k = 23
      do 150 j = 1,12
        resc24 = resc24+cheb24(k)*chebmo(m,k)
        ress24 = ress24+cheb24(k+1)*chebmo(m,k+1)
        resabs = dabs(cheb24(k))+dabs(cheb24(k+1))
        k = k-2
  150 continue
      estc = dabs(resc24-resc12)
      ests = dabs(ress24-ress12)
      resabs = resabs*dabs(hlgth)
      if(integr.eq.2) go to 160
      result = conc*resc24-cons*ress24
      abserr = dabs(conc*estc)+dabs(cons*ests)
      go to 170
  160 result = conc*ress24+cons*resc24
      abserr = dabs(conc*ests)+dabs(cons*estc)
  170 return
      end
