      subroutine fpgrre(ifsx,ifsy,ifbx,ifby,x,mx,y,my,z,mz,kx,ky,tx,nx,
     * ty,ny,p,c,nc,fp,fpx,fpy,mm,mynx,kx1,kx2,ky1,ky2,spx,spy,right,q,
     * ax,ay,bx,by,nrx,nry)
c  ..
c  ..scalar arguments..
      real*8 p,fp
      integer ifsx,ifsy,ifbx,ifby,mx,my,mz,kx,ky,nx,ny,nc,mm,mynx,
     * kx1,kx2,ky1,ky2
c  ..array arguments..
      real*8 x(mx),y(my),z(mz),tx(nx),ty(ny),c(nc),spx(mx,kx1),spy(my,ky
     *1)
     * ,right(mm),q(mynx),ax(nx,kx2),bx(nx,kx2),ay(ny,ky2),by(ny,ky2),
     * fpx(nx),fpy(ny)
      integer nrx(mx),nry(my)
c  ..local scalars..
      real*8 arg,cos,fac,pinv,piv,sin,term,one,half
      integer i,ibandx,ibandy,ic,iq,irot,it,iz,i1,i2,i3,j,k,k1,k2,l,
     * l1,l2,ncof,nk1x,nk1y,nrold,nroldx,nroldy,number,numx,numx1,
     * numy,numy1,n1
c  ..local arrays..
      real*8 h(7)
c  ..subroutine references..
c    fpback,fpbspl,fpgivs,fpdisc,fprota
c  ..
c  the b-spline coefficients of the smoothing spline are calculated as
c  the least-squares solution of the over-determined linear system of
c  equations  (ay) c (ax)' = q       where
c
c               |   (spx)    |            |   (spy)    |
c        (ax) = | ---------- |     (ay) = | ---------- |
c               | (1/p) (bx) |            | (1/p) (by) |
c
c                                | z  ' 0 |
c                            q = | ------ |
c                                | 0  ' 0 |
c
c  with c      : the (ny-ky-1) x (nx-kx-1) matrix which contains the
c                b-spline coefficients.
c       z      : the my x mx matrix which contains the function values.
c       spx,spy: the mx x (nx-kx-1) and  my x (ny-ky-1) observation
c                matrices according to the least-squares problems in
c                the x- and y-direction.
c       bx,by  : the (nx-2*kx-1) x (nx-kx-1) and (ny-2*ky-1) x (ny-ky-1)
c                matrices which contain the discontinuity jumps of the
c                derivatives of the b-splines in the x- and y-direction.
      one = 1
      half = 0.5
      nk1x = nx-kx1
      nk1y = ny-ky1
      if(p.gt.0.) pinv = one/p
c  it depends on the value of the flags ifsx,ifsy,ifbx and ifby and on
c  the value of p whether the matrices (spx),(spy),(bx) and (by) still
c  must be determined.
      if(ifsx.ne.0) go to 50
c  calculate the non-zero elements of the matrix (spx) which is the
c  observation matrix according to the least-squares spline approximat-
c  ion problem in the x-direction.
      l = kx1
      l1 = kx2
      number = 0
      do 40 it=1,mx
        arg = x(it)
  10    if(arg.lt.tx(l1) .or. l.eq.nk1x) go to 20
        l = l1
        l1 = l+1
        number = number+1
        go to 10
  20    call fpbspl(tx,nx,kx,arg,l,h)
        do 30 i=1,kx1
          spx(it,i) = h(i)
  30    continue
        nrx(it) = number
  40  continue
      ifsx = 1
  50  if(ifsy.ne.0) go to 100
c  calculate the non-zero elements of the matrix (spy) which is the
c  observation matrix according to the least-squares spline approximat-
c  ion problem in the y-direction.
      l = ky1
      l1 = ky2
      number = 0
      do 90 it=1,my
        arg = y(it)
  60    if(arg.lt.ty(l1) .or. l.eq.nk1y) go to 70
        l = l1
        l1 = l+1
        number = number+1
        go to 60
  70    call fpbspl(ty,ny,ky,arg,l,h)
        do 80 i=1,ky1
          spy(it,i) = h(i)
  80    continue
        nry(it) = number
  90  continue
      ifsy = 1
 100  if(p.le.0.) go to 120
c  calculate the non-zero elements of the matrix (bx).
      if(ifbx.ne.0 .or. nx.eq.2*kx1) go to 110
      call fpdisc(tx,nx,kx2,bx,nx)
      ifbx = 1
c  calculate the non-zero elements of the matrix (by).
 110  if(ifby.ne.0 .or. ny.eq.2*ky1) go to 120
      call fpdisc(ty,ny,ky2,by,ny)
      ifby = 1
c  reduce the matrix (ax) to upper triangular form (rx) using givens
c  rotations. apply the same transformations to the rows of matrix q
c  to obtain the my x (nx-kx-1) matrix g.
c  store matrix (rx) into (ax) and g into q.
 120  l = my*nk1x
c  initialization.
      do 130 i=1,l
        q(i) = 0.
 130  continue
      do 140 i=1,nk1x
        do 140 j=1,kx2
          ax(i,j) = 0.
 140  continue
      l = 0
      nrold = 0
c  ibandx denotes the bandwidth of the matrices (ax) and (rx).
      ibandx = kx1
      do 270 it=1,mx
        number = nrx(it)
 150    if(nrold.eq.number) go to 180
        if(p.le.0.) go to 260
        ibandx = kx2
c  fetch a new row of matrix (bx).
        n1 = nrold+1
        do 160 j=1,kx2
          h(j) = bx(n1,j)*pinv
 160    continue
c  find the appropriate column of q.
        do 170 j=1,my
          right(j) = 0.
 170    continue
        irot = nrold
        go to 210
c  fetch a new row of matrix (spx).
 180    h(ibandx) = 0.
        do 190 j=1,kx1
          h(j) = spx(it,j)
 190    continue
c  find the appropriate column of q.
        do 200 j=1,my
          l = l+1
          right(j) = z(l)
 200    continue
        irot = number
c  rotate the new row of matrix (ax) into triangle.
 210    do 240 i=1,ibandx
          irot = irot+1
          piv = h(i)
          if(piv.eq.0.) go to 240
c  calculate the parameters of the givens transformation.
          call fpgivs(piv,ax(irot,1),cos,sin)
c  apply that transformation to the rows of matrix q.
          iq = (irot-1)*my
          do 220 j=1,my
            iq = iq+1
            call fprota(cos,sin,right(j),q(iq))
 220      continue
c  apply that transformation to the columns of (ax).
          if(i.eq.ibandx) go to 250
          i2 = 1
          i3 = i+1
          do 230 j=i3,ibandx
            i2 = i2+1
            call fprota(cos,sin,h(j),ax(irot,i2))
 230      continue
 240    continue
 250    if(nrold.eq.number) go to 270
 260    nrold = nrold+1
        go to 150
 270  continue
c  reduce the matrix (ay) to upper triangular form (ry) using givens
c  rotations. apply the same transformations to the columns of matrix g
c  to obtain the (ny-ky-1) x (nx-kx-1) matrix h.
c  store matrix (ry) into (ay) and h into c.
      ncof = nk1x*nk1y
c  initialization.
      do 280 i=1,ncof
        c(i) = 0.
 280  continue
      do 290 i=1,nk1y
        do 290 j=1,ky2
          ay(i,j) = 0.
 290  continue
      nrold = 0
c  ibandy denotes the bandwidth of the matrices (ay) and (ry).
      ibandy = ky1
      do 420 it=1,my
        number = nry(it)
 300    if(nrold.eq.number) go to 330
        if(p.le.0.) go to 410
        ibandy = ky2
c  fetch a new row of matrix (by).
        n1 = nrold+1
        do 310 j=1,ky2
          h(j) = by(n1,j)*pinv
 310    continue
c  find the appropriate row of g.
        do 320 j=1,nk1x
          right(j) = 0.
 320    continue
        irot = nrold
        go to 360
c  fetch a new row of matrix (spy)
 330    h(ibandy) = 0.
        do 340 j=1,ky1
          h(j) = spy(it,j)
 340    continue
c  find the appropriate row of g.
        l = it
        do 350 j=1,nk1x
          right(j) = q(l)
          l = l+my
 350    continue
        irot = number
c  rotate the new row of matrix (ay) into triangle.
 360    do 390 i=1,ibandy
          irot = irot+1
          piv = h(i)
          if(piv.eq.0.) go to 390
c  calculate the parameters of the givens transformation.
          call fpgivs(piv,ay(irot,1),cos,sin)
c  apply that transformation to the columns of matrix g.
          ic = irot
          do 370 j=1,nk1x
            call fprota(cos,sin,right(j),c(ic))
            ic = ic+nk1y
 370      continue
c  apply that transformation to the columns of matrix (ay).
          if(i.eq.ibandy) go to 400
          i2 = 1
          i3 = i+1
          do 380 j=i3,ibandy
            i2 = i2+1
            call fprota(cos,sin,h(j),ay(irot,i2))
 380      continue
 390    continue
 400    if(nrold.eq.number) go to 420
 410    nrold = nrold+1
        go to 300
 420  continue
c  backward substitution to obtain the b-spline coefficients as the
c  solution of the linear system    (ry) c (rx)' = h.
c  first step: solve the system  (ry) (c1) = h.
      k = 1
      do 450 i=1,nk1x
        call fpback(ay,c(k),nk1y,ibandy,c(k),ny)
        k = k+nk1y
 450  continue
c  second step: solve the system  c (rx)' = (c1).
      k = 0
      do 480 j=1,nk1y
        k = k+1
        l = k
        do 460 i=1,nk1x
          right(i) = c(l)
          l = l+nk1y
 460    continue
        call fpback(ax,right,nk1x,ibandx,right,nx)
        l = k
        do 470 i=1,nk1x
          c(l) = right(i)
          l = l+nk1y
 470    continue
 480  continue
c  calculate the quantities
c    res(i,j) = (z(i,j) - s(x(i),y(j)))**2 , i=1,2,..,mx;j=1,2,..,my
c    fp = sumi=1,mx(sumj=1,my(res(i,j)))
c    fpx(r) = sum''i(sumj=1,my(res(i,j))) , r=1,2,...,nx-2*kx-1
c                  tx(r+kx) <= x(i) <= tx(r+kx+1)
c    fpy(r) = sumi=1,mx(sum''j(res(i,j))) , r=1,2,...,ny-2*ky-1
c                  ty(r+ky) <= y(j) <= ty(r+ky+1)
      fp = 0.
      do 490 i=1,nx
        fpx(i) = 0.
 490  continue
      do 500 i=1,ny
        fpy(i) = 0.
 500  continue
      nk1y = ny-ky1
      iz = 0
      nroldx = 0
c  main loop for the different grid points.
      do 550 i1=1,mx
        numx = nrx(i1)
        numx1 = numx+1
        nroldy = 0
        do 540 i2=1,my
          numy = nry(i2)
          numy1 = numy+1
          iz = iz+1
c  evaluate s(x,y) at the current grid point by making the sum of the
c  cross products of the non-zero b-splines at (x,y), multiplied with
c  the appropriate b-spline coefficients.
          term = 0.
          k1 = numx*nk1y+numy
          do 520 l1=1,kx1
            k2 = k1
            fac = spx(i1,l1)
            do 510 l2=1,ky1
              k2 = k2+1
              term = term+fac*spy(i2,l2)*c(k2)
 510        continue
            k1 = k1+nk1y
 520      continue
c  calculate the squared residual at the current grid point.
          term = (z(iz)-term)**2
c  adjust the different parameters.
          fp = fp+term
          fpx(numx1) = fpx(numx1)+term
          fpy(numy1) = fpy(numy1)+term
          fac = term*half
          if(numy.eq.nroldy) go to 530
          fpy(numy1) = fpy(numy1)-fac
          fpy(numy) = fpy(numy)+fac
 530      nroldy = numy
          if(numx.eq.nroldx) go to 540
          fpx(numx1) = fpx(numx1)-fac
          fpx(numx) = fpx(numx)+fac
 540    continue
        nroldx = numx
 550  continue
      return
      end

