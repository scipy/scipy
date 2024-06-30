      recursive subroutine dqk61(f,a,b,result,abserr,resabs,resasc)
c***begin prologue  dqk61
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  61-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f over (a,b) with error
c                           estimate
c                       j = integral of dabs(f) over (a,b)
c***description
c
c        integration rule
c        standard fortran subroutine
c        double precision version
c
c
c        parameters
c         on entry
c           f      - double precision
c                    function subprogram defining the integrand
c                    function f(x). the actual name for f needs to be
c                    declared e x t e r n a l in the calling program.
c
c           a      - double precision
c                    lower limit of integration
c
c           b      - double precision
c                    upper limit of integration
c
c         on return
c           result - double precision
c                    approximation to the integral i
c                    result is computed by applying the 61-point
c                    kronrod rule (resk) obtained by optimal addition of
c                    abscissae to the 30-point gauss rule (resg).
c
c           abserr - double precision
c                    estimate of the modulus of the absolute error,
c                    which should equal or exceed dabs(i-result)
c
c           resabs - double precision
c                    approximation to the integral j
c
c           resasc - double precision
c                    approximation to the integral of dabs(f-i/(b-a))
c
c
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk61
c
      double precision a,dabsc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,
     *  d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc,
     *  resg,resk,reskh,result,uflow,wg,wgk,xgk
      integer j,jtw,jtwm1
      external f
c
      dimension fv1(30),fv2(30),xgk(31),wgk(31),wg(15)
c
c           the abscissae and weights are given for the
c           interval (-1,1). because of symmetry only the positive
c           abscissae and their corresponding weights are given.
c
c           xgk   - abscissae of the 61-point kronrod rule
c                   xgk(2), xgk(4)  ... abscissae of the 30-point
c                   gauss rule
c                   xgk(1), xgk(3)  ... optimally added abscissae
c                   to the 30-point gauss rule
c
c           wgk   - weights of the 61-point kronrod rule
c
c           wg    - weigths of the 30-point gauss rule
c
c
c gauss quadrature weights and kronron quadrature abscissae and weights
c as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
c bell labs, nov. 1981.
c
      data wg  (  1) / 0.0079681924 9616660561 5465883474 674 d0 /
      data wg  (  2) / 0.0184664683 1109095914 2302131912 047 d0 /
      data wg  (  3) / 0.0287847078 8332336934 9719179611 292 d0 /
      data wg  (  4) / 0.0387991925 6962704959 6801936446 348 d0 /
      data wg  (  5) / 0.0484026728 3059405290 2938140422 808 d0 /
      data wg  (  6) / 0.0574931562 1761906648 1721689402 056 d0 /
      data wg  (  7) / 0.0659742298 8218049512 8128515115 962 d0 /
      data wg  (  8) / 0.0737559747 3770520626 8243850022 191 d0 /
      data wg  (  9) / 0.0807558952 2942021535 4694938460 530 d0 /
      data wg  ( 10) / 0.0868997872 0108297980 2387530715 126 d0 /
      data wg  ( 11) / 0.0921225222 3778612871 7632707087 619 d0 /
      data wg  ( 12) / 0.0963687371 7464425963 9468626351 810 d0 /
      data wg  ( 13) / 0.0995934205 8679526706 2780282103 569 d0 /
      data wg  ( 14) / 0.1017623897 4840550459 6428952168 554 d0 /
      data wg  ( 15) / 0.1028526528 9355884034 1285636705 415 d0 /
c
      data xgk (  1) / 0.9994844100 5049063757 1325895705 811 d0 /
      data xgk (  2) / 0.9968934840 7464954027 1630050918 695 d0 /
      data xgk (  3) / 0.9916309968 7040459485 8628366109 486 d0 /
      data xgk (  4) / 0.9836681232 7974720997 0032581605 663 d0 /
      data xgk (  5) / 0.9731163225 0112626837 4693868423 707 d0 /
      data xgk (  6) / 0.9600218649 6830751221 6871025581 798 d0 /
      data xgk (  7) / 0.9443744447 4855997941 5831324037 439 d0 /
      data xgk (  8) / 0.9262000474 2927432587 9324277080 474 d0 /
      data xgk (  9) / 0.9055733076 9990779854 6522558925 958 d0 /
      data xgk ( 10) / 0.8825605357 9205268154 3116462530 226 d0 /
      data xgk ( 11) / 0.8572052335 4606109895 8658510658 944 d0 /
      data xgk ( 12) / 0.8295657623 8276839744 2898119732 502 d0 /
      data xgk ( 13) / 0.7997278358 2183908301 3668942322 683 d0 /
      data xgk ( 14) / 0.7677774321 0482619491 7977340974 503 d0 /
      data xgk ( 15) / 0.7337900624 5322680472 6171131369 528 d0 /
      data xgk ( 16) / 0.6978504947 9331579693 2292388026 640 d0 /
      data xgk ( 17) / 0.6600610641 2662696137 0053668149 271 d0 /
      data xgk ( 18) / 0.6205261829 8924286114 0477556431 189 d0 /
      data xgk ( 19) / 0.5793452358 2636169175 6024932172 540 d0 /
      data xgk ( 20) / 0.5366241481 4201989926 4169793311 073 d0 /
      data xgk ( 21) / 0.4924804678 6177857499 3693061207 709 d0 /
      data xgk ( 22) / 0.4470337695 3808917678 0609900322 854 d0 /
      data xgk ( 23) / 0.4004012548 3039439253 5476211542 661 d0 /
      data xgk ( 24) / 0.3527047255 3087811347 1037207089 374 d0 /
      data xgk ( 25) / 0.3040732022 7362507737 2677107199 257 d0 /
      data xgk ( 26) / 0.2546369261 6788984643 9805129817 805 d0 /
      data xgk ( 27) / 0.2045251166 8230989143 8957671002 025 d0 /
      data xgk ( 28) / 0.1538699136 0858354696 3794672743 256 d0 /
      data xgk ( 29) / 0.1028069379 6673703014 7096751318 001 d0 /
      data xgk ( 30) / 0.0514718425 5531769583 3025213166 723 d0 /
      data xgk ( 31) / 0.0000000000 0000000000 0000000000 000 d0 /
c
      data wgk (  1) / 0.0013890136 9867700762 4551591226 760 d0 /
      data wgk (  2) / 0.0038904611 2709988405 1267201844 516 d0 /
      data wgk (  3) / 0.0066307039 1593129217 3319826369 750 d0 /
      data wgk (  4) / 0.0092732796 5951776342 8441146892 024 d0 /
      data wgk (  5) / 0.0118230152 5349634174 2232898853 251 d0 /
      data wgk (  6) / 0.0143697295 0704580481 2451432443 580 d0 /
      data wgk (  7) / 0.0169208891 8905327262 7572289420 322 d0 /
      data wgk (  8) / 0.0194141411 9394238117 3408951050 128 d0 /
      data wgk (  9) / 0.0218280358 2160919229 7167485738 339 d0 /
      data wgk ( 10) / 0.0241911620 7808060136 5686370725 232 d0 /
      data wgk ( 11) / 0.0265099548 8233310161 0601709335 075 d0 /
      data wgk ( 12) / 0.0287540487 6504129284 3978785354 334 d0 /
      data wgk ( 13) / 0.0309072575 6238776247 2884252943 092 d0 /
      data wgk ( 14) / 0.0329814470 5748372603 1814191016 854 d0 /
      data wgk ( 15) / 0.0349793380 2806002413 7499670731 468 d0 /
      data wgk ( 16) / 0.0368823646 5182122922 3911065617 136 d0 /
      data wgk ( 17) / 0.0386789456 2472759295 0348651532 281 d0 /
      data wgk ( 18) / 0.0403745389 5153595911 1995279752 468 d0 /
      data wgk ( 19) / 0.0419698102 1516424614 7147541285 970 d0 /
      data wgk ( 20) / 0.0434525397 0135606931 6831728117 073 d0 /
      data wgk ( 21) / 0.0448148001 3316266319 2355551616 723 d0 /
      data wgk ( 22) / 0.0460592382 7100698811 6271735559 374 d0 /
      data wgk ( 23) / 0.0471855465 6929915394 5261478181 099 d0 /
      data wgk ( 24) / 0.0481858617 5708712914 0779492298 305 d0 /
      data wgk ( 25) / 0.0490554345 5502977888 7528165367 238 d0 /
      data wgk ( 26) / 0.0497956834 2707420635 7811569379 942 d0 /
      data wgk ( 27) / 0.0504059214 0278234684 0893085653 585 d0 /
      data wgk ( 28) / 0.0508817958 9874960649 2297473049 805 d0 /
      data wgk ( 29) / 0.0512215478 4925877217 0656282604 944 d0 /
      data wgk ( 30) / 0.0514261285 3745902593 3862879215 781 d0 /
      data wgk ( 31) / 0.0514947294 2945156755 8340433647 099 d0 /
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           dabsc  - abscissa
c           fval*  - function value
c           resg   - result of the 30-point gauss rule
c           resk   - result of the 61-point kronrod rule
c           reskh  - approximation to the mean value of f
c                    over (a,b), i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
      epmach = d1mach(4)
      uflow = d1mach(1)
c
      centr = 0.5d+00*(b+a)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
c
c           compute the 61-point kronrod approximation to the
c           integral, and estimate the absolute error.
c
c***first executable statement  dqk61
      resg = 0.0d+00
      fc = f(centr)
      resk = wgk(31)*fc
      resabs = dabs(resk)
      do 10 j=1,15
        jtw = j*2
        dabsc = hlgth*xgk(jtw)
        fval1 = f(centr-dabsc)
        fval2 = f(centr+dabsc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j=1,15
        jtwm1 = j*2-1
        dabsc = hlgth*xgk(jtwm1)
        fval1 = f(centr-dabsc)
        fval2 = f(centr+dabsc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
  15    continue
      reskh = resk*0.5d+00
      resasc = wgk(31)*dabs(fc-reskh)
      do 20 j=1,30
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
