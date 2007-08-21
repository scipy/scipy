      subroutine dqc25c(f,a,b,c,result,abserr,krul,neval)
c***begin prologue  dqc25c
c***date written   810101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a2,j4
c***keywords  25-point clenshaw-curtis integration
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f*w over (a,b) with
c            error estimate, where w(x) = 1/(x-c)
c***description
c
c        integration rules for the computation of cauchy
c        principal value integrals
c        standard fortran subroutine
c        double precision version
c
c        parameters
c           f      - double precision
c                    function subprogram defining the integrand function
c                    f(x). the actual name for f needs to be declared
c                    e x t e r n a l  in the driver program.
c
c           a      - double precision
c                    left end point of the integration interval
c
c           b      - double precision
c                    right end point of the integration interval, b.gt.a
c
c           c      - double precision
c                    parameter in the weight function
c
c           result - double precision
c                    approximation to the integral
c                    result is computed by using a generalized
c                    clenshaw-curtis method if c lies within ten percent
c                    of the integration interval. in the other case the
c                    15-point kronrod rule obtained by optimal addition
c                    of abscissae to the 7-point gauss rule, is applied.
c
c           abserr - double precision
c                    estimate of the modulus of the absolute error,
c                    which should equal or exceed abs(i-result)
c
c           krul   - integer
c                    key which is decreased by 1 if the 15-point
c                    gauss-kronrod scheme has been used
c
c           neval  - integer
c                    number of integrand evaluations
c
c.......................................................................
c***references  (none)
c***routines called  dqcheb,dqk15w,dqwgtc
c***end prologue  dqc25c
c
      double precision a,abserr,ak22,amom0,amom1,amom2,b,c,cc,centr,
     *  cheb12,cheb24,dabs,dlog,dqwgtc,f,fval,hlgth,p2,p3,p4,resabs,
     *  resasc,result,res12,res24,u,x
      integer i,isym,k,kp,krul,neval
c
      dimension x(11),fval(25),cheb12(13),cheb24(25)
c
      external f,dqwgtc
c
c           the vector x contains the values cos(k*pi/24),
c           k = 1, ..., 11, to be used for the chebyshev series
c           expansion of f
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
c           ----------------------
c           fval   - value of the function f at the points
c                    cos(k*pi/24),  k = 0, ..., 24
c           cheb12 - chebyshev series expansion coefficients,
c                    for the function f, of degree 12
c           cheb24 - chebyshev series expansion coefficients,
c                    for the function f, of degree 24
c           res12  - approximation to the integral corresponding
c                    to the use of cheb12
c           res24  - approximation to the integral corresponding
c                    to the use of cheb24
c           dqwgtc - external function subprogram defining
c                    the weight function
c           hlgth  - half-length of the interval
c           centr  - mid point of the interval
c
c
c           check the position of c.
c
c***first executable statement  dqc25c
      cc = (0.2d+01*c-b-a)/(b-a)
      if(dabs(cc).lt.0.11d+01) go to 10
c
c           apply the 15-point gauss-kronrod scheme.
c
      krul = krul-1
      call dqk15w(f,dqwgtc,c,p2,p3,p4,kp,a,b,result,abserr,
     *  resabs,resasc)
      neval = 15
      if (resasc.eq.abserr) krul = krul+1
      go to 50
c
c           use the generalized clenshaw-curtis method.
c
   10 hlgth = 0.5d+00*(b-a)
      centr = 0.5d+00*(b+a)
      neval = 25
      fval(1) = 0.5d+00*f(hlgth+centr)
      fval(13) = f(centr)
      fval(25) = 0.5d+00*f(centr-hlgth)
      do 20 i=2,12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = f(u+centr)
        fval(isym) = f(centr-u)
   20 continue
c
c           compute the chebyshev series expansion.
c
      call dqcheb(x,fval,cheb12,cheb24)
c
c           the modified chebyshev moments are computed by forward
c           recursion, using amom0 and amom1 as starting values.
c
      amom0 = dlog(dabs((0.1d+01-cc)/(0.1d+01+cc)))
      amom1 = 0.2d+01+cc*amom0
      res12 = cheb12(1)*amom0+cheb12(2)*amom1
      res24 = cheb24(1)*amom0+cheb24(2)*amom1
      do 30 k=3,13
        amom2 = 0.2d+01*cc*amom1-amom0
        ak22 = (k-2)*(k-2)
        if((k/2)*2.eq.k) amom2 = amom2-0.4d+01/(ak22-0.1d+01)
        res12 = res12+cheb12(k)*amom2
        res24 = res24+cheb24(k)*amom2
        amom0 = amom1
        amom1 = amom2
   30 continue
      do 40 k=14,25
        amom2 = 0.2d+01*cc*amom1-amom0
        ak22 = (k-2)*(k-2)
        if((k/2)*2.eq.k) amom2 = amom2-0.4d+01/(ak22-0.1d+01)
        res24 = res24+cheb24(k)*amom2
        amom0 = amom1
        amom1 = amom2
   40 continue
      result = res24
      abserr = dabs(res24-res12)
   50 return
      end
