      subroutine surev(idim,tu,nu,tv,nv,c,u,mu,v,mv,f,mf,wrk,lwrk,
     * iwrk,kwrk,ier)
c  subroutine surev evaluates on a grid (u(i),v(j)),i=1,...,mu; j=1,...
c  ,mv a bicubic spline surface of dimension idim, given in the
c  b-spline representation.
c
c  calling sequence:
c     call surev(idim,tu,nu,tv,nv,c,u,mu,v,mv,f,mf,wrk,lwrk,
c    * iwrk,kwrk,ier)
c
c  input parameters:
c   idim  : integer, specifying the dimension of the spline surface.
c   tu    : real array, length nu, which contains the position of the
c           knots in the u-direction.
c   nu    : integer, giving the total number of knots in the u-direction
c   tv    : real array, length nv, which contains the position of the
c           knots in the v-direction.
c   nv    : integer, giving the total number of knots in the v-direction
c   c     : real array, length (nu-4)*(nv-4)*idim, which contains the
c           b-spline coefficients.
c   u     : real array of dimension (mu).
c           before entry u(i) must be set to the u co-ordinate of the
c           i-th grid point along the u-axis.
c           tu(4)<=u(i-1)<=u(i)<=tu(nu-3), i=2,...,mu.
c   mu    : on entry mu must specify the number of grid points along
c           the u-axis. mu >=1.
c   v     : real array of dimension (mv).
c           before entry v(j) must be set to the v co-ordinate of the
c           j-th grid point along the v-axis.
c           tv(4)<=v(j-1)<=v(j)<=tv(nv-3), j=2,...,mv.
c   mv    : on entry mv must specify the number of grid points along
c           the v-axis. mv >=1.
c   mf    : on entry, mf must specify the dimension of the array f.
c           mf >= mu*mv*idim
c   wrk   : real array of dimension lwrk. used as workspace.
c   lwrk  : integer, specifying the dimension of wrk.
c           lwrk >= 4*(mu+mv)
c   iwrk  : integer array of dimension kwrk. used as workspace.
c   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mu+mv.
c
c  output parameters:
c   f     : real array of dimension (mf).
c           on succesful exit f(mu*mv*(l-1)+mv*(i-1)+j) contains the
c           l-th co-ordinate of the bicubic spline surface at the
c           point (u(i),v(j)),l=1,...,idim,i=1,...,mu;j=1,...,mv.
c   ier   : integer error flag
c    ier=0 : normal return
c    ier=10: invalid input data (see restrictions)
c
c  restrictions:
c   mu >=1, mv >=1, lwrk>=4*(mu+mv), kwrk>=mu+mv , mf>=mu*mv*idim
c   tu(4) <= u(i-1) <= u(i) <= tu(nu-3), i=2,...,mu
c   tv(4) <= v(j-1) <= v(j) <= tv(nv-3), j=2,...,mv
c
c  other subroutines required:
c    fpsuev,fpbspl
c
c  references :
c    de boor c : on calculating with b-splines, j. approximation theory
c                6 (1972) 50-62.
c    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
c                applics 10 (1972) 134-149.
c    dierckx p. : curve and surface fitting with splines, monographs on
c                 numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      integer idim,nu,nv,mu,mv,mf,lwrk,kwrk,ier
c  ..array arguments..
      integer iwrk(kwrk)
      real*8 tu(nu),tv(nv),c((nu-4)*(nv-4)*idim),u(mu),v(mv),f(mf),
     * wrk(lwrk)
c  ..local scalars..
      integer i,muv
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(mf.lt.mu*mv*idim) go to 100
      muv = mu+mv
      if(lwrk.lt.4*muv) go to 100
      if(kwrk.lt.muv) go to 100
      if (mu.lt.1) go to 100
      if (mu.eq.1) go to 30
      go to 10
  10  do 20 i=2,mu
        if(u(i).lt.u(i-1)) go to 100
  20  continue
  30  if (mv.lt.1) go to 100
      if (mv.eq.1) go to 60
      go to 40
  40  do 50 i=2,mv
        if(v(i).lt.v(i-1)) go to 100
  50  continue
  60  ier = 0
      call fpsuev(idim,tu,nu,tv,nv,c,u,mu,v,mv,f,wrk(1),wrk(4*mu+1),
     * iwrk(1),iwrk(mu+1))
 100  return
      end
