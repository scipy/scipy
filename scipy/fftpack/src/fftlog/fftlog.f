c=======================================================================
c FFTLog
c    computes the discrete Fast Fourier Transform
c    or Fast Hankel Transform (of arbitrary real index)
c    of a periodic logarithmic sequence.
c
c-----------------------------------------------------------------------
c
c Permission to distribute this modified fftlog.f under the BSD-3-Clause
c license has been granted (email from Andrew Hamilton to Dieter
c Werthmüller dated 07 October 2016).
c
c Modifications (all marked with comments: %DW):
c   - drfft* -> dfft*
c   - kr: remove kropt from fhti (see getkr)
c         => kr has to be defined exactly before calling fhti
c         => remove parameters : kropt, ok, lnblnk, stblnk, go, temp
c         => remove function   : stblnk(string)
c   - krgood: change from function krgood to subroutine getkr
c   - split wsave into wsave and xsave
c         => wsave is as in regular fft
c         => xsave is the fftlog-addition to wsave
c         => moved dffti from fhti to fhtq
c
c---Start of revision history-------------------------------------------
c
c Original version 27 March 1999.
c
c 5 July 1999:
c Bug fix in singular transform.
c Store amp, arg of last element of wsave separately, for even n,
c enlarging wsave by 1, for even n and q != 0.

c 17 July 1999:
c Reorganized wsave more logically.
c drffti now gets first, not last, 2*n+15 elements of wsave.
c fhti now treats last element of wsave with others, not by itself;
c wsave is enlarged from original, for even n,
c by 1 for q = 0, and by 2 for q != 0.
c
c 18 July 1999:
c A backward transform (dir = -1) is the same as
c a forward transform with q -> -q [and rk -> 1/rk]
c for any kr if n is odd,
c but only for low-ringing kr if n is even.
c Original falsely stated that transforms were same irrespective of kr.
c
c 19 Dec 1999:
c Slightly spiffed up interactive option to change kr to low-ringing kr.
c
c 13 Mar 2000:
c Made g77-safe:
c 1. Checked for any automatic variables that needed explicitly
c    `save'ing - actually there weren't any.
c    g77 is unusual in that it does NOT save variables not in
c    common or data statements.
c 2. All double precision constants explicitly suffixed d0.
c    g77 is unusual in that it does NOT automatically promote
c    real constants such as 1.1 to double precision 1.1d0.
c 
c---End of revision history---------------------------------------------
c
c For more information about FFTLog, see
c
c http://casa.colorado.edu/~ajsh/FFTLog/
c
c Andrew J S Hamilton March 1999
c email: Andrew.Hamilton@Colorado.EDU
c
c Refs:	Talman J. D., 1978, J. Comp. Phys., 29, 35
c	Hamilton A. J. S., 2000, MNRAS, 312, 257
c	( http://xxx.lanl.gov/abs/astro-ph/9905191 )
c
c FFTLog uses the NCAR suite of FFT routines,
c and a modified version of the complex Gamma function
c from the gamerf package at
c http://momonga.t.u-tokyo.ac.jp/~ooura/gamerf.html .
c The original gamerf copyright statement states:
c   Copyright(C) 1996 Takuya OOURA (email: ooura@mmm.t.u-tokyo.ac.jp).
c   You may use, copy, modify this code for any purpose and
c   without fee. You may distribute this ORIGINAL package.
c
c Permission to distribute the modified gamma function code
c with the FFTLog package has been granted
c (email from Takuya Ooura to Andrew Hamilton dated 16 March 1999).
c
c-----------------------------------------------------------------------
c Note: email programs may zap the 8-bit characters in the program
c       comments by converting them into 7-bit characters.
c       If this happens, get an unzapped copy from
c       http://casa.colorado.edu/~ajsh/FFTLog/ .
c       The fortran code itself contains only 7-bit characters;
c       the 8-bit characters occur only in comments.
c
c       If you are on a UNIX machine and 8-bit characters are not
c       showing up properly (they may appear as \nnn where n are
c       integers), then try putting
c          setenv LC_CTYPE iso_8859_1
c       in your .login file.
c
c-----------------------------------------------------------------------
c FFTLog computes a discrete version of
c the Hankel Transform (= Fourier-Bessel Transform)
c with a power law bias (k r)^q
c
c          infinity
c           /           q
c    ã(k) = | a(r) (k r)  J  (k r) k dr
c           /              mu
c          0
c
c          infinity
c           /           -q
c    a(r) = | ã(k) (k r)   J  (k r) r dk
c           /               mu
c          0
c
c where J_mu is the Bessel function of order mu.
c The index mu may be any real number, positive or negative.
c
c The input array a_j is a periodic sequence of length n,
c uniformly logarithmically spaced with spacing dlnr
c    a_j = a(r_j)   at   r_j = r_c exp[(j-j_c) dlnr]
c centred about the point r_c.  The central index j_c = (n+1)/2
c is 1/2 integral if n is even.
c Similarly, the output array ã_j is a periodic sequence of length n,
c also uniformly logarithmically spaced with spacing dlnr
c    ã_j = ã(k_j)   at   k_j = k_c exp[(j-j_c) dlnr]
c centred about the point k_c.
c
c The centre points r_c and k_c of the periodic intervals may be
c chosen arbitrarily; but it would be normal to choose the product
c kr = k_c r_c = k_j r_(n+1-j) = k_(n+1-j) r_j
c to be about 1 (or 2, or pi, to taste).
c
c The FFTLog algorithm is (see Hamilton 2000):
c 1. FFT the input array a_j to obtain the Fourier coefficients c_m ;
c 2. Multiply c_m by 
c       u_m = (kr)^[- i 2 m pi/(n dlnr)] U_mu[q + i 2 m pi/(n dlnr)]
c    where
c       U_mu(x) = 2^x Gamma[(mu+1+x)/2] / Gamma[(mu+1-x)/2]
c    to obtain c_m u_m ;
c 3. FFT c_m u_m back to obtain the discrete Hankel transform ã_j .
c
c-----------------------------------------------------------------------
c The Fourier sine and cosine transforms
c
c                     infinity
c                      /
c    Ã(k) = sqrt(2/pi) | A(r) sin(k r) dr
c                      /
c                     0
c
c                     infinity
c                      /
c    Ã(k) = sqrt(2/pi) | A(r) cos(k r) dr
c                      /
c                     0
c
c may be regarded as special cases of the Hankel transform
c with mu = 1/2 and -1/2 since
c
c    sqrt(2/pi) sin(x) = sqrt(x) J   (x)
c                                 1/2
c
c    sqrt(2/pi) cos(x) = sqrt(x) J    (x)
c                                 -1/2
c
c The Fourier transforms may be done by making the substitutions
c                 q-(1/2)                      -q-(1/2)
c    A(r) = a(r) r          and   Ã(k) = ã(k) k    
c
c and Hankel transforming a(r) with a power law bias (k r)^q
c
c          infinity
c           /           q
c    ã(k) = | a(r) (k r)  J    (k r) k dr
c           /              ±1/2
c          0
c
c Different choices of power law bias q lead to different discrete
c Fourier transforms of A(r), because the assumption of periodicity of
c a(r) = A(r) r^[-q+(1/2)] is different for different q.
c
c If A(r) is a power law, A(r) proportional to r^[q-(1/2)],
c then applying a bias q yields a discrete Fourier transform Ã(k)
c that is exactly equal to the continuous Fourier transform,
c because then a(r) is a constant, which is a periodic function.
c
c-----------------------------------------------------------------------
c The Hankel transform
c
c          infinity
c           /
c    Ã(k) = | A(r) J  (k r) k dr
c           /       mu
c          0
c
c may be done by making the substitutions
c                 q                      -q
c    A(r) = a(r) r    and   Ã(k) = ã(k) k    
c
c and Hankel transforming a(r) with a power law bias (k r)^q
c
c          infinity
c           /           q
c    ã(k) = | a(r) (k r)  J  (k r) k dr
c           /              mu
c          0
c
c Different choices of power law bias q lead to different discrete
c Hankel transforms of A(r), because the assumption of periodicity of
c a(r) = A(r) r^-q is different for different q.
c
c If A(r) is a power law, A(r) proportional to r^q,
c then applying a bias q yields a discrete Hankel transform Ã(k)
c that is exactly equal to the continuous Hankel transform,
c because then a(r) is a constant, which is a periodic function.
c
c-----------------------------------------------------------------------
c ------------------------
c There are five routines:
c ------------------------
c Comments in the subroutines contain further details.
c
c (1) subroutine fhti(n,mu,q,dlnr,kr,kropt,wsave,ok)
c     is an initialization routine.
c
c (2) subroutine fftl(n,a,rk,dir,wsave)
c     computes the discrete Fourier sine or cosine transform
c     of a logarithmically spaced periodic sequence.
c     This is a driver routine that calls fhtq.
c
c (3) subroutine fht(n,a,dir,wsave)
c     computes the discrete Hankel transform
c     of a logarithmically spaced periodic sequence.
c     This is a driver routine that calls fhtq.
c
c (4) subroutine fhtq(n,a,dir,wsave)
c     computes the biased discrete Hankel transform
c     of a logarithmically spaced periodic sequence.
c     This is the basic FFTLog routine.
c
c (5) real*8 function krgood(mu,q,dlnr,kr)
c     takes an input kr and returns the nearest low-ringing kr.
c     This is an optional routine called by fhti.
c
c=======================================================================
c THE FFTLog CODE
c=======================================================================
!      subroutine fhti(n,mu,q,dlnr,kr,kropt,wsave,ok)
      subroutine fhti(n,mu,q,dlnr,kr,xsave) !%DW remove kropt,ok
      integer n  !%DW remove kropt          !    wsave-> xsave
!      logical ok !%DW remove ok
      real*8 mu,q,dlnr,kr,xsave(*) !%DW wsave->xsave
c
c fhti initializes the working array wsave
c used by fftl, fht, and fhtq, and by NCAR routines drfftf and drfftb.
c fhti need be called once, whereafter fftl, fht, or fhtq may be
c called many times, as long as n, mu, q, dlnr, and kr remain unchanged.
c fhti should be called each time n, mu, q, dlnr, or kr is changed.
c The work array wsave should not be changed between calls to fftl,
c fht, or fhtq.
c
c If you are using the g77 fortran compiler, which by default does not
c save locally declared variables (unlike most other fortran compilers),
c then you should explicitly
c       save wsave
c in the calling program.
c
c  Input: n = number of points in the array to be transformed;
c             n may be any positive integer, but the NCAR FFT routines
c             run fastest if n is a product of small primes 2, 3, 5.
c	  mu = index of J_mu in Hankel transform;
c	       mu may be any real number, positive or negative.
c	  q = exponent of power law bias;
c	      q may be any real number, positive or negative.
c             If in doubt, use q = 0, for which case the Hankel
c             transform is orthogonal, i.e. self-inverse,
c             provided also that, for n even, kr is low-ringing.
c             Non-zero q may yield better approximations to the
c             continuous Hankel transform for some functions.
c         dlnr = separation between natural log of points;
c		 dlnr may be positive or negative.
c	  kr = k_c r_c where c is central point of array
c              = k_j r_(n+1-j) = k_(n+1-j) r_j .
c	       Normally one would choose kr to be about 1
c	       (or 2, or pi, to taste).
c         kropt = 0 to use input kr as is;
c                 1 to change kr to nearest low-ringing kr, quietly;
c                 2 to change kr to nearest low-ringing kr, verbosely;
c                 3 for option to change kr interactively.
c Output: wsave = working array used by fftl, fht, and fhtq,
c                 and by NCAR FFT routines drfftf and drfftb.
c                 wsave should be dimensioned at least:
c                 for q = 0 (unbiased transform):
c                    2*n+2*(n/2)+18
c                       [= 3*n+18 if n even, 3*n+17 if n odd]
c                 for q != 0 (biased transform):
c                    2*n+3*(n/2)+19
c                       [= (7*n)/2+19 if n even, (7*n)/2+18 if n odd]
c                 The first 2*n+15 elements of wsave are used by
c                 the NCAR FFT routines.
c         ok = .true. if all went ok;
c              .false. if not; currently only occurs if interactive
c                      response causes exit.
c
c        parameters
      real*8 PI
      parameter (PI=3.141592653589793238462643383279502884197d0)
      real*8 ZERO,ONE,TWO
      parameter (ZERO=0.d0,ONE=1.d0,TWO=2.d0)
      real*8 ROUND
      parameter (ROUND=1.d-15)
c        externals
!      integer lnblnk,stblnk  %DW removed, unused
      real*8 krgood
      complex*16 cdgamma
c        local (automatic) variables
!      character*1 go  %DW removed, unused
!      character*64 temp  %DW removed, unused
      integer l,m
      real*8 amp,arg,d,ln2,ln2kr,xm,xp,y
      complex*16 zm,zp
c
!-----%DW removed kropt, kr has to be provided exatly => subroutine getkr
!c--------adjust kr
!c        keep kr as is
!      if (kropt.eq.0) then
!        continue
!c        change kr to low-ringing kr quietly
!      elseif (kropt.eq.1) then
!        kr=krgood(mu,q,dlnr,kr)
!c        change kr to low-ringing kr verbosely
!      elseif (kropt.eq.2) then
!        d=krgood(mu,q,dlnr,kr)
!        if (abs(kr/d-ONE).gt.ROUND) then
!          kr=d
!          write (*,'(" kr changed to",g24.16)') kr
!        endif
!c        option to change kr to low-ringing kr interactively
!      else
!        d=krgood(mu,q,dlnr,kr)
!        if (abs(kr/d-ONE).gt.ROUND) then
!c        fortran demonstrates its inferiority to C
!          write (*,'(" change kr = ",$)')
!          write (temp,*) kr
!          write (*,'(a,$)') temp(stblnk(temp):lnblnk(temp))
!          write (*,'(" to low-ringing kr = ",$)')
!          write (temp,*) d
!          write (*,'(a,$)') temp(stblnk(temp):lnblnk(temp))
!          write (*,'("? [CR,y=yes, n=no, x=exit]: ",$)')
!          read (*,'(a1)') go
!          if (go.eq.' '.or.go.eq.'y'.or.go.eq.'Y') then
!            kr=d
!            write (*,'(" kr changed to",g24.16)') kr
!          elseif (go.eq.'n'.or.go.eq.'N') then
!            print *,'kr left unchanged at',kr
!          else
!            print *,'exit'
!            goto 300
!          endif
!        endif
!      endif
c--------return if n is <= 0
      if (n.le.0) goto 200
c--------initialize normal FFT  !%DW move to fhtq
!     call dffti(n,wsave)       !%DW Change drffti->dffti
c        drfft uses first 2*n+15 elements of wsave
      l=0 !%DW adjust for xsave 2*n+15->0
c--------put q, dlnr, kr in next 3 elements of wsave
      xsave(l+1)=q    !%DW wsave->xsave
      xsave(l+2)=dlnr !    "
      xsave(l+3)=kr   !    "
c        so far used 2*n+18 elements of wsave
      l=l+3
c--------rest of wsave used by fhtq - unbiased case (q = 0)
      if (q.eq.ZERO) then
        ln2kr=log(TWO/kr)
        xp=(mu+ONE)/TWO
        d=PI/(n*dlnr)
        do m=1,n/2
c        y = m pi/(n dlnr)
          y=m*d
          zp=dcmplx(xp,y)
          zp=cdgamma(zp,1)
c        Argument of kr^(-2 i y) U_mu(2 i y)
          arg=TWO*(ln2kr*y+dimag(zp))
          xsave(l+2*m-1)=cos(arg) !%DW wsave->xsave
          xsave(l+2*m)=sin(arg)   !    "
        enddo
c Altogether 2*n+18 + 2*(n/2) = 2*((3*n)/2)+18 elements used for q = 0,
c which is 3*n+18 for even n, 3*n+17 for odd n.
c--------rest of wsave used by fhtq - biased case (q != 0)
      else
        ln2=log(TWO)
        ln2kr=log(TWO/kr)
        xp=(mu+ONE+q)/TWO
        xm=(mu+ONE-q)/TWO
c........first element of rest of wsave
        y=ZERO
c        case where xp or xm is a negative integer
        if ((anint(xp).eq.xp.and.anint(xp).le.ZERO)
     *    .or.(anint(xm).eq.xm.and.anint(xm).le.ZERO)) then
c        case where xp and xm are both negative integers
c U_mu(q) = 2^q Gamma[xp]/Gamma[xm] is finite in this case
          if ((anint(xp).eq.xp.and.anint(xp).le.ZERO)
     *      .and.(anint(xm).eq.xm.and.anint(xm).le.ZERO)) then
c        Amplitude and Argument of U_mu(q)
            amp=exp(ln2*q)
            if (xp.gt.xm) then
              do m=1,nint(xp-xm)
                amp=amp*(xm+m-ONE)
              enddo
            elseif (xp.lt.xm) then
              do m=1,nint(xm-xp)
                amp=amp/(xp+m-ONE)
              enddo
            endif
            arg=anint(xp+xm)*PI
c        one of xp or xm is a negative integer
          else
c Transformation is singular if xp is -ve integer,
c and inverse transformation is singular if xm is -ve integer,
c but transformation may be well-defined if sum_j a_j = 0,
c as may well occur in physical cases.
c Policy is to drop the potentially infinite constant in the transform.
            if (anint(xp).eq.xp.and.anint(xp).le.ZERO) then
              print *,'fhti: (mu+1+q)/2 =',nint(xp),
     &          ' is -ve integer, yields singular transform:'
              print *,'   transform will omit additive constant that is 
     &generically infinite,'
              print *,'   but that may be finite or zero if'
              print *,'   the sum of the elements of the input array a_j
     & is zero.'
            else
              print *,'fhti: (mu+1-q)/2 =',nint(xm),
     &          ' is -ve integer, yields singular inverse transform:'
              print *,'   inverse transform will omit additive constant 
     &that is generically'
              print *,'   infinite, but that may be finite or zero if'
              print *,'   the sum of the elements of the input array a_j
     & is zero.'
            endif
            amp=ZERO
            arg=ZERO
          endif
c        neither xp nor xm is a negative integer
        else
          zp=dcmplx(xp,y)
          zm=dcmplx(xm,y)
          zp=cdgamma(zp,1)
          zm=cdgamma(zm,1)
c        Amplitude and Argument of U_mu(q)
          amp=exp(ln2*q+dble(zp)-dble(zm))
c        note +Im(zm) to get conjugate value below real axis
          arg=dimag(zp)+dimag(zm)
        endif
c        cos(arg) = ±1, sin(arg) = 0
        xsave(l+1)=amp*cos(arg) !%DW wsave->xsave
c........remaining elements of wsave
        d=PI/(n*dlnr)
        do m=1,n/2
c        y = m pi/(n dlnr)
          y=m*d
          zp=dcmplx(xp,y)
          zm=dcmplx(xm,y)
          zp=cdgamma(zp,1)
          zm=cdgamma(zm,1)
c        Amplitude and Argument of kr^(-2 i y) U_mu(q + 2 i y)
          amp=exp(ln2*q+dble(zp)-dble(zm))
          arg=2*ln2kr*y+dimag(zp)+dimag(zm)
          xsave(l+3*m-1)=amp      !%DW wsave->xsave
          xsave(l+3*m)=cos(arg)   !    "
          xsave(l+3*m+1)=sin(arg) !    "
        enddo
c Altogether 2*n+18 + 3*(n/2)+1 elements used for q != 0,
c which is (7*n)/2+19 for even n, (7*n)/2+18 for odd n.
c For even n, the very last element of wsave
c [i.e. wsave(l+3*m+1)=sin(arg) for m=n/2] is not used within FFTLog;
c if a low-ringing kr is used, this element should be zero.
c The last element is computed in case somebody wants it.
      endif
  200 continue  !%DW ok=.true. -> continue
      return
c
c--------error returns
!  300 ok=.false.  %DW remove ok
!      return      %DW "
      end
c
c=======================================================================
      subroutine fftl(n,a,rk,dir,wsave,xsave) !%DW add xsave
      integer n,dir
      real*8 a(n),rk,wsave(*),xsave(*) !%DW add xsave
c
c This is a driver routine that calls fhtq.
c
c fftl computes a discrete version of the Fourier
c sine (if mu = 1/2) or cosine (if mu = -1/2) transform
c
c                     infinity
c                      /
c    Ã(k) = sqrt(2/pi) | A(r) sin(k r) dr
c                      /
c                     0
c
c                     infinity
c                      /
c    Ã(k) = sqrt(2/pi) | A(r) cos(k r) dr
c                      /
c                     0
c
c by making the substitutions
c                 q-(1/2)                      -q-(1/2)
c    A(r) = a(r) r          and   Ã(k) = ã(k) k    
c
c and applying a biased Hankel transform to a(r).
c
c The steps are:
c 1. a(r) = A(r) r^[-dir*(q-.5)]
c 2. call fhtq to transform a(r) -> ã(k)
c 3. Ã(k) = ã(k) k^[-dir*(q+.5)]
c
c fhti must be called before the first call to fftl,
c with mu = 1/2 for a sine transform,
c or mu = -1/2 for a cosine transform.
c
c A call to fftl with dir=1 followed by
c a call to fftl with dir=-1 (and rk unchanged), or vice versa,
c leaves the array a unchanged.
c
c  Input: n = length of a array.
c         rk = r_c/k_c
c              = r_j/k_j (a constant, the same constant for any j);
c              rk is not (necessarily) the same quantity as kr.
c              rk is used only to multiply the output array by
c              sqrt(rk)^dir, so if you want to do the normalization
c              later, or you don't care about the normalization,
c              you can set rk = 1.
c	  dir = 1 for forward transform,
c		-1 for backward transform.
c               A backward transform (dir = -1) is the same as
c               a forward transform with q -> -q and rk -> 1/rk,
c               for any kr if n is odd,
c               for low-ringing kr if n is even.
c	  wsave = working array set up by fhti.
c Input/Output:
c         a on  input is the array A(r) to transform:
c             a(j) is A(r_j) at r_j = r_c exp[(j-jc) dlnr]
c                where jc = (n+1)/2 = central index of array.
c         a on output is the transformed array Ã(k):
c             a(j) is Ã(k_j) at k_j = k_c exp[(j-jc) dlnr].
c
c        parameters
      real*8 ONE,TWO,HALF
      parameter (ONE=1.d0,TWO=2.d0,HALF=ONE/TWO)
c        local (automatic) variables
      integer j,l
      real*8 dlnr,jc,kr,lnkr,lnrk,q
c
      l=0 !%DW adjust for xsave 2*n+15->0
      q=xsave(l+1)    !%DW wsave->xsave
      dlnr=xsave(l+2) !    "
      kr=xsave(l+3)   !    "
c........a(r) = A(r) (r/rc)^[-dir*(q-.5)]
c        centre point of array
      jc=dble(n+1)/TWO
      do j=1,n
        a(j)=a(j)*exp(-dir*(q-HALF)*(j-jc)*dlnr)
      enddo
c........transform a(r) -> ã(k)
      call fhtq(n,a,dir,wsave,xsave) !%DW add xsave
c........Ã(k) = ã(k) k^[-dir*(q+.5)] rc^[-dir*(q-.5)]
c        = ã(k) (k/kc)^[-dir*(q+.5)] (kc rc)^(-dir*q) (rc/kc)^(dir*.5)
      lnkr=log(kr)
      lnrk=log(rk)
      do j=1,n
        a(j)=a(j)*exp(-dir*((q+HALF)*(j-jc)*dlnr+q*lnkr-HALF*lnrk))
      enddo
      return
      end
c
c=======================================================================
      subroutine fht(n,a,dir,wsave,xsave) !%DW add xsave
      integer n,dir
      real*8 a(n),wsave(*),xsave(*) !%DW add xsave
c
c This is a driver routine that calls fhtq.
c
c fht computes a discrete version of the Hankel transform
c
c          infinity
c           /
c    Ã(k) = | A(r) J  (k r) k dr
c           /       mu
c          0
c
c by making the substitutions
c                 q                      -q
c    A(r) = a(r) r    and   Ã(k) = ã(k) k    
c
c and applying a biased Hankel transform to a(r).
c
c The steps are:
c 1. a(r) = A(r) r^(-dir*q)
c 2. call fhtq to transform a(r) -> ã(k)
c 3. Ã(k) = ã(k) k^(-dir*q)
c
c fhti must be called before the first call to fht.
c
c A call to fht with dir=1 followed by
c a call to fht with dir=-1, or vice versa,
c leaves the array a unchanged.
c
c  Input: n = length of a array.
c	  dir = 1 for forward transform,
c		-1 for backward transform.
c               A backward transform (dir = -1) is the same as
c               a forward transform with q -> -q,
c               for any kr if n is odd,
c               for low-ringing kr if n is even.
c	  wsave = working array set up by fhti.
c Input/Output:
c         a on  input is the array A(r) to transform:
c             a(j) is A(r_j) at r_j = r_c exp[(j-jc) dlnr]
c                where jc = (n+1)/2 = central index of array.
c         a on output is the transformed array Ã(k):
c             a(j) is Ã(k_j) at k_j = k_c exp[(j-jc) dlnr].
c
c        parameters
      real*8 ZERO,TWO
      parameter (ZERO=0.d0,TWO=2.d0)
c        local (automatic) variables
      integer j,l
      real*8 dlnr,jc,kr,lnkr,q
c
      l=0 !%DW adjust for xsave 2*n+15->0
      q=xsave(l+1)    !%DW wsave -> xsave
      dlnr=xsave(l+2) !    "
      kr=xsave(l+3)   !    "
c........a(r) = A(r) (r/rc)^(-dir*q)
      if (q.ne.ZERO) then
c        centre point of array
        jc=dble(n+1)/TWO
        do j=1,n
          a(j)=a(j)*exp(-dir*q*(j-jc)*dlnr)
        enddo
      endif
c........transform a(r) -> ã(k)
      call fhtq(n,a,dir,wsave,xsave) !%DW add xsave
c........Ã(k) = ã(k) (k rc)^(-dir*q)
c             = ã(k) (k/kc)^(-dir*q) (kc rc)^(-dir*q)
      if (q.ne.ZERO) then
        lnkr=log(kr)
        do j=1,n
          a(j)=a(j)*exp(-dir*q*((j-jc)*dlnr+lnkr))
        enddo
      endif
      return
      end
c
c=======================================================================
      subroutine fhtq(n,a,dir,wsave,xsave) !%DW add xsave
      integer n,dir
      real*8 a(n),wsave(*),xsave(*) !%DW add xsave
c
c This is the basic FFTLog routine.
c
c fhtq computes a discrete version of the biased Hankel transform
c
c          infinity
c           /           q
c    ã(k) = | a(r) (k r)  J  (k r) k dr
c           /              mu
c          0
c
c fhti must be called before the first call to fhtq.
c
c A call to fhtq with dir=1 followed by
c a call to fhtq with dir=-1, or vice versa,
c leaves the array a unchanged.
c
c  Input: n = length of a array.
c	  dir = 1 for forward transform,
c		-1 for backward transform.
c               A backward transform (dir = -1) is the same as
c               a forward transform with q -> -q,
c               for any kr if n is odd,
c               for low-ringing kr if n is even.
c Input/Output:
c         a on  input is the periodic array a(r) to transform:
c             a(j) is a(r_j) at r_j = r_c exp[(j-jc) dlnr]
c                 where jc = (n+1)/2 = central index of array.
c         a on output is the transformed periodic array ã(k):
c             a(j) is ã(k_j) at k_j = k_c exp[(j-jc) dlnr].
c
c        parameters
      real*8 ZERO
      parameter (ZERO=0.d0)
c        local (automatic) variables
      integer l,m
      real*8 ai,ar,q
c
c--------initialize normal FFT !%DW moved from fhti
      call dffti(n,wsave)      !%DW change drffti->dffti

      l=0 !%DW adjust for xsave 2*n+15->0
      q=xsave(l+1) !%DW wsave -> xsave
      l=l+3
c--------normal FFT
      call dfftf(n,a,wsave) !%DW change drfftf->dfftf
c--------unbiased (q = 0) transform
      if (q.eq.ZERO) then
c........multiply by
c        (kr)^[- i 2 m pi/(n dlnr)] U_mu[i 2 m pi/(n dlnr)]
        do m=1,(n-1)/2
          ar=a(2*m)
          ai=a(2*m+1)
          a(2*m)=ar*xsave(l+2*m-1)-ai*xsave(l+2*m)   !%DW wsave -> xsave
          a(2*m+1)=ar*xsave(l+2*m)+ai*xsave(l+2*m-1) !    "
        enddo
c        problematical last element, for even n
        if (mod(n,2).eq.0) then
          ar=xsave(l+n-1) !%DW wsave -> xsave
c        forward transform: multiply by real part
c Why?  See http://casa.colorado.edu/~ajsh/FFTLog/index.html#ure
          if (dir.eq.1) then
            a(n)=a(n)*ar
c        backward transform: divide by real part
          elseif (dir.eq.-1) then
c Real part ar can be zero for maximally bad choice of kr.
c This is unlikely to happen by chance, but if it does,
c policy is to let it happen.
c For low-ringing kr, imaginary part ai is zero by construction,
c and real part ar is guaranteed nonzero.
            a(n)=a(n)/ar
          endif
        endif
c--------biased (q != 0) transform
      else
c........multiply by
c        (kr)^[- i 2 m pi/(n dlnr)] U_mu[q + i 2 m pi/(n dlnr)]
c        phase
        do m=1,(n-1)/2
          ar=a(2*m)
          ai=a(2*m+1)
          a(2*m)=ar*xsave(l+3*m)-ai*xsave(l+3*m+1)   !%DW wsave -> xsave
          a(2*m+1)=ar*xsave(l+3*m+1)+ai*xsave(l+3*m) !    "
        enddo
c        forward transform: multiply by amplitude
        if (dir.eq.1) then
          a(1)=a(1)*xsave(l+1) !%DW wsave -> xsave
          do m=1,(n-1)/2
            a(2*m)=a(2*m)*xsave(l+3*m-1)     !%DW wsave -> xsave
            a(2*m+1)=a(2*m+1)*xsave(l+3*m-1) !    "
          enddo
c        backward transform: divide by amplitude
        elseif (dir.eq.-1) then
c        amplitude of m=0 element
          ar=xsave(l+1) !%DW wsave -> xsave
          if (ar.eq.ZERO) then
c Amplitude of m=0 element can be zero for some mu, q combinations
c (singular inverse); policy is to drop potentially infinite constant.
            a(1)=ZERO
          else
            a(1)=a(1)/ar
          endif
c        remaining amplitudes should never be zero
          do m=1,(n-1)/2
            a(2*m)=a(2*m)/xsave(l+3*m-1)     !%DW wsave -> xsave
            a(2*m+1)=a(2*m+1)/xsave(l+3*m-1) !    "
          enddo
        endif
c        problematical last element, for even n
        if (mod(n,2).eq.0) then
          m=n/2
          ar=xsave(l+3*m)*xsave(l+3*m-1) !%DW wsave -> xsave
c        forward transform: multiply by real part
          if (dir.eq.1) then
            a(n)=a(n)*ar
c        backward transform: divide by real part
          elseif (dir.eq.-1) then
c Real part ar can be zero for maximally bad choice of kr.
c This is unlikely to happen by chance, but if it does,
c policy is to let it happen.
c For low-ringing kr, imaginary part ai is zero by construction,
c and real part ar is guaranteed nonzero.
            a(n)=a(n)/ar
          endif
        endif
      endif
c--------normal FFT back
      call dfftb(n,a,wsave) !%DW change drfftb->dfftb
c--------reverse the array
c        and at the same time undo the FFTs' multiplication by n 
      do m=1,n/2
        ar=a(m)
        a(m)=a(n+1-m)/n
        a(n+1-m)=ar/n
      enddo
      if (mod(n,2).eq.1) then
        m=(n+1)/2
        a(m)=a(m)/n
      endif
      return
      end
c
c=======================================================================
!      real*8 function krgood(mu,q,dlnr,kr) !%DW replace function with
      subroutine getkr(mu,q,dlnr,kr)        !    subroutine
      real*8 mu,q,dlnr,kr
c
c Use of this routine is optional.
c
c Choosing kr so that
c     (kr)^(- i pi/dlnr) U_mu(q + i pi/dlnr)
c is real may reduce ringing of the discrete Hankel transform,
c because it makes the transition of this function across the period
c boundary smoother.
c
c  Input: mu = index of J_mu in Hankel transform.
c	  q = exponent of power law bias.
c         dlnr = separation between natural log of points.
c	  kr = suggested value of kr.
c Output: krgood = low-ringing value of kr nearest to input kr.
c                  ln(krgood) is always within dlnr/2 of ln(kr).
c
c        parameters
      real*8 PI
      parameter (PI=3.141592653589793238462643383279502884197d0)
      real*8 ZERO,ONE,TWO
      parameter (ZERO=0.d0,ONE=1.d0,TWO=2.d0)
c        externals
      complex*16 cdgamma
c        local (automatic) variables
      real*8 arg,iarg,xm,xp,y
      complex*16 zm,zp
c
!      krgood=kr !%DW no krgood, adjust kr in place
      if (dlnr.eq.ZERO) return
      xp=(mu+ONE+q)/TWO
      xm=(mu+ONE-q)/TWO
      y=PI/(TWO*dlnr)
      zp=dcmplx(xp,y)
      zm=dcmplx(xm,y)
      zp=cdgamma(zp,1)
      zm=cdgamma(zm,1)
c        low-ringing condition is that following should be integral
      arg=log(TWO/kr)/dlnr+(dimag(zp)+dimag(zm))/PI
      iarg=anint(arg)
c        if should ensure arg = +-Infinity dealt with correctly
      if (arg.ne.iarg) then
c        low-ringing kr
        kr=kr*exp((arg-iarg)*dlnr) !%DW no krgood, adjust kr in place
      endif
      return
      end
c
c-----------------------------------------------------------------------
!      integer function stblnk(string)  %DW comment stblnk, not used
!      character*(*) string
!c
!c        externals
!      integer lnblnk
!c        local (automatic) variables
!      integer i
!c *
!c * Return index of first non-blank character in string.
!c * If string is all blank, returned index is 1+lnblnk(string).
!c *
!      do i=1,lnblnk(string)
!        if (string(i:i).ne.' ') goto 120
!      enddo
!  120 stblnk=i
!      return
!      end
!c
