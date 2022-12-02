      recursive subroutine dqk31(f,a,b,result,abserr,resabs,resasc)
c***begin prologue  dqk31
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  31-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f over (a,b) with error
c                           estimate
c                       j = integral of abs(f) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           double precision version
c
c           parameters
c            on entry
c              f      - double precision
c                       function subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the calling program.
c
c              a      - double precision
c                       lower limit of integration
c
c              b      - double precision
c                       upper limit of integration
c
c            on return
c              result - double precision
c                       approximation to the integral i
c                       result is computed by applying the 31-point
c                       gauss-kronrod rule (resk), obtained by optimal
c                       addition of abscissae to the 15-point gauss
c                       rule (resg).
c
c              abserr - double precision
c                       estimate of the modulus of the modulus,
c                       which should not exceed abs(i-result)
c
c              resabs - double precision
c                       approximation to the integral j
c
c              resasc - double precision
c                       approximation to the integral of abs(f-i/(b-a))
c                       over (a,b)
c
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk31
      double precision a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,
     *  d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc,
     *  resg,resk,reskh,result,uflow,wg,wgk,xgk
      integer j,jtw,jtwm1
      external f
c
      dimension fv1(15),fv2(15),xgk(16),wgk(16),wg(8)
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 31-point kronrod rule
c                    xgk(2), xgk(4), ...  abscissae of the 15-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 15-point gauss rule
c
c           wgk    - weights of the 31-point kronrod rule
c
c           wg     - weights of the 15-point gauss rule
c
c
c gauss quadrature weights and kronron quadrature abscissae and weights
c as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
c bell labs, nov. 1981.
c
      data wg  (  1) / 0.0307532419 9611726835 4628393577 204 d0 /
      data wg  (  2) / 0.0703660474 8810812470 9267416450 667 d0 /
      data wg  (  3) / 0.1071592204 6717193501 1869546685 869 d0 /
      data wg  (  4) / 0.1395706779 2615431444 7804794511 028 d0 /
      data wg  (  5) / 0.1662692058 1699393355 3200860481 209 d0 /
      data wg  (  6) / 0.1861610000 1556221102 6800561866 423 d0 /
      data wg  (  7) / 0.1984314853 2711157645 6118326443 839 d0 /
      data wg  (  8) / 0.2025782419 2556127288 0620199967 519 d0 /
c
      data xgk (  1) / 0.9980022986 9339706028 5172840152 271 d0 /
      data xgk (  2) / 0.9879925180 2048542848 9565718586 613 d0 /
      data xgk (  3) / 0.9677390756 7913913425 7347978784 337 d0 /
      data xgk (  4) / 0.9372733924 0070590430 7758947710 209 d0 /
      data xgk (  5) / 0.8972645323 4408190088 2509656454 496 d0 /
      data xgk (  6) / 0.8482065834 1042721620 0648320774 217 d0 /
      data xgk (  7) / 0.7904185014 4246593296 7649294817 947 d0 /
      data xgk (  8) / 0.7244177313 6017004741 6186054613 938 d0 /
      data xgk (  9) / 0.6509967412 9741697053 3735895313 275 d0 /
      data xgk ( 10) / 0.5709721726 0853884753 7226737253 911 d0 /
      data xgk ( 11) / 0.4850818636 4023968069 3655740232 351 d0 /
      data xgk ( 12) / 0.3941513470 7756336989 7207370981 045 d0 /
      data xgk ( 13) / 0.2991800071 5316881216 6780024266 389 d0 /
      data xgk ( 14) / 0.2011940939 9743452230 0628303394 596 d0 /
      data xgk ( 15) / 0.1011420669 1871749902 7074231447 392 d0 /
      data xgk ( 16) / 0.0000000000 0000000000 0000000000 000 d0 /
c
      data wgk (  1) / 0.0053774798 7292334898 7792051430 128 d0 /
      data wgk (  2) / 0.0150079473 2931612253 8374763075 807 d0 /
      data wgk (  3) / 0.0254608473 2671532018 6874001019 653 d0 /
      data wgk (  4) / 0.0353463607 9137584622 2037948478 360 d0 /
      data wgk (  5) / 0.0445897513 2476487660 8227299373 280 d0 /
      data wgk (  6) / 0.0534815246 9092808726 5343147239 430 d0 /
      data wgk (  7) / 0.0620095678 0067064028 5139230960 803 d0 /
      data wgk (  8) / 0.0698541213 1872825870 9520077099 147 d0 /
      data wgk (  9) / 0.0768496807 5772037889 4432777482 659 d0 /
      data wgk ( 10) / 0.0830805028 2313302103 8289247286 104 d0 /
      data wgk ( 11) / 0.0885644430 5621177064 7275443693 774 d0 /
      data wgk ( 12) / 0.0931265981 7082532122 5486872747 346 d0 /
      data wgk ( 13) / 0.0966427269 8362367850 5179907627 589 d0 /
      data wgk ( 14) / 0.0991735987 2179195933 2393173484 603 d0 /
      data wgk ( 15) / 0.1007698455 2387559504 4946662617 570 d0 /
      data wgk ( 16) / 0.1013300070 1479154901 7374792767 493 d0 /
c
c
c           list of major variables
c           -----------------------
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc   - abscissa
c           fval*  - function value
c           resg   - result of the 15-point gauss formula
c           resk   - result of the 31-point kronrod formula
c           reskh  - approximation to the mean value of f over (a,b),
c                    i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c***first executable statement  dqk31
      epmach = d1mach(4)
      uflow = d1mach(1)
c
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
c
c           compute the 31-point kronrod approximation to
c           the integral, and estimate the absolute error.
c
      fc = f(centr)
      resg = wg(8)*fc
      resk = wgk(16)*fc
      resabs = dabs(resk)
      do 10 j=1,7
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j = 1,8
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(16)*dabs(fc-reskh)
      do 20 j=1,15
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)
     *  abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1
     *  ((epmach*0.5d+02)*resabs,abserr)
      return
      end
