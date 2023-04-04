      recursive subroutine fpsuev(idim,tu,nu,tv,nv,c,u,mu,v,mv,f,
     *   wu,wv,lu,lv)
      implicit none
c  ..scalar arguments..
      integer idim,nu,nv,mu,mv
c  ..array arguments..
      integer lu(mu),lv(mv)
      real*8 tu(nu),tv(nv),c((nu-4)*(nv-4)*idim),u(mu),v(mv),
     * f(mu*mv*idim),wu(mu,4),wv(mv,4)
c  ..local scalars..
      integer i,i1,j,j1,k,l,l1,l2,l3,m,nuv,nu4,nv4
      real*8 arg,sp,tb,te
c  ..local arrays..
      real*8 h(4)
c  ..subroutine references..
c    fpbspl
c  ..
      nu4 = nu-4
      tb = tu(4)
      te = tu(nu4+1)
      l = 4
      l1 = l+1
      do 40 i=1,mu
        arg = u(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
  10    if(arg.lt.tu(l1) .or. l.eq.nu4) go to 20
        l = l1
        l1 = l+1
        go to 10
  20    call fpbspl(tu,nu,3,arg,l,h)
        lu(i) = l-4
        do 30 j=1,4
          wu(i,j) = h(j)
  30    continue
  40  continue
      nv4 = nv-4
      tb = tv(4)
      te = tv(nv4+1)
      l = 4
      l1 = l+1
      do 80 i=1,mv
        arg = v(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
  50    if(arg.lt.tv(l1) .or. l.eq.nv4) go to 60
        l = l1
        l1 = l+1
        go to 50
  60    call fpbspl(tv,nv,3,arg,l,h)
        lv(i) = l-4
        do 70 j=1,4
          wv(i,j) = h(j)
  70    continue
  80  continue
      m = 0
      nuv = nu4*nv4
      do 140 k=1,idim
        l3 = (k-1)*nuv
        do 130 i=1,mu
          l = lu(i)*nv4+l3
          do 90 i1=1,4
            h(i1) = wu(i,i1)
  90      continue
          do 120 j=1,mv
            l1 = l+lv(j)
            sp = 0.
            do 110 i1=1,4
              l2 = l1
              do 100 j1=1,4
                l2 = l2+1
                sp = sp+c(l2)*h(i1)*wv(j,j1)
 100          continue
              l1 = l1+nv4
 110        continue
            m = m+1
            f(m) = sp
 120      continue
 130    continue
 140  continue
      return
      end
