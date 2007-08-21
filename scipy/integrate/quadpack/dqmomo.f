      subroutine dqmomo(alfa,beta,ri,rj,rg,rh,integr)
c***begin prologue  dqmomo
c***date written   820101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1,c3a2
c***keywords  modified chebyshev moments
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  this routine computes modified chebsyshev moments. the k-th
c            modified chebyshev moment is defined as the integral over
c            (-1,1) of w(x)*t(k,x), where t(k,x) is the chebyshev
c            polynomial of degree k.
c***description
c
c        modified chebyshev moments
c        standard fortran subroutine
c        double precision version
c
c        parameters
c           alfa   - double precision
c                    parameter in the weight function w(x), alfa.gt.(-1)
c
c           beta   - double precision
c                    parameter in the weight function w(x), beta.gt.(-1)
c
c           ri     - double precision
c                    vector of dimension 25
c                    ri(k) is the integral over (-1,1) of
c                    (1+x)**alfa*t(k-1,x), k = 1, ..., 25.
c
c           rj     - double precision
c                    vector of dimension 25
c                    rj(k) is the integral over (-1,1) of
c                    (1-x)**beta*t(k-1,x), k = 1, ..., 25.
c
c           rg     - double precision
c                    vector of dimension 25
c                    rg(k) is the integral over (-1,1) of
c                    (1+x)**alfa*log((1+x)/2)*t(k-1,x), k = 1, ..., 25.
c
c           rh     - double precision
c                    vector of dimension 25
c                    rh(k) is the integral over (-1,1) of
c                    (1-x)**beta*log((1-x)/2)*t(k-1,x), k = 1, ..., 25.
c
c           integr - integer
c                    input parameter indicating the modified
c                    moments to be computed
c                    integr = 1 compute ri, rj
c                           = 2 compute ri, rj, rg
c                           = 3 compute ri, rj, rh
c                           = 4 compute ri, rj, rg, rh
c
c***references  (none)
c***routines called  (none)
c***end prologue  dqmomo
c
      double precision alfa,alfp1,alfp2,an,anm1,beta,betp1,betp2,ralf,
     *  rbet,rg,rh,ri,rj
      integer i,im1,integr
c
      dimension rg(25),rh(25),ri(25),rj(25)
c
c
c***first executable statement  dqmomo
      alfp1 = alfa+0.1d+01
      betp1 = beta+0.1d+01
      alfp2 = alfa+0.2d+01
      betp2 = beta+0.2d+01
      ralf = 0.2d+01**alfp1
      rbet = 0.2d+01**betp1
c
c           compute ri, rj using a forward recurrence relation.
c
      ri(1) = ralf/alfp1
      rj(1) = rbet/betp1
      ri(2) = ri(1)*alfa/alfp2
      rj(2) = rj(1)*beta/betp2
      an = 0.2d+01
      anm1 = 0.1d+01
      do 20 i=3,25
        ri(i) = -(ralf+an*(an-alfp2)*ri(i-1))/(anm1*(an+alfp1))
        rj(i) = -(rbet+an*(an-betp2)*rj(i-1))/(anm1*(an+betp1))
        anm1 = an
        an = an+0.1d+01
   20 continue
      if(integr.eq.1) go to 70
      if(integr.eq.3) go to 40
c
c           compute rg using a forward recurrence relation.
c
      rg(1) = -ri(1)/alfp1
      rg(2) = -(ralf+ralf)/(alfp2*alfp2)-rg(1)
      an = 0.2d+01
      anm1 = 0.1d+01
      im1 = 2
      do 30 i=3,25
        rg(i) = -(an*(an-alfp2)*rg(im1)-an*ri(im1)+anm1*ri(i))/
     *  (anm1*(an+alfp1))
        anm1 = an
        an = an+0.1d+01
        im1 = i
   30 continue
      if(integr.eq.2) go to 70
c
c           compute rh using a forward recurrence relation.
c
   40 rh(1) = -rj(1)/betp1
      rh(2) = -(rbet+rbet)/(betp2*betp2)-rh(1)
      an = 0.2d+01
      anm1 = 0.1d+01
      im1 = 2
      do 50 i=3,25
        rh(i) = -(an*(an-betp2)*rh(im1)-an*rj(im1)+
     *  anm1*rj(i))/(anm1*(an+betp1))
        anm1 = an
        an = an+0.1d+01
        im1 = i
   50 continue
      do 60 i=2,25,2
        rh(i) = -rh(i)
   60 continue
   70 do 80 i=2,25,2
        rj(i) = -rj(i)
   80 continue
   90 return
      end
