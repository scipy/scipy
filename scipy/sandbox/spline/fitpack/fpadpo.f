      subroutine fpadpo(idim,t,n,c,nc,k,cp,np,cc,t1,t2)
c  given a idim-dimensional spline curve of degree k, in its b-spline
c  representation ( knots t(j),j=1,...,n , b-spline coefficients c(j),
c  j=1,...,nc) and given also a polynomial curve in its b-spline
c  representation ( coefficients cp(j), j=1,...,np), subroutine fpadpo
c  calculates the b-spline representation (coefficients c(j),j=1,...,nc)
c  of the sum of the two curves.
c
c  other subroutine required : fpinst
c
c  ..
c  ..scalar arguments..
      integer idim,k,n,nc,np
c  ..array arguments..
      real*8 t(n),c(nc),cp(np),cc(nc),t1(n),t2(n)
c  ..local scalars..
      integer i,ii,j,jj,k1,l,l1,n1,n2,nk1,nk2
c  ..
      k1 = k+1
      nk1 = n-k1
c  initialization
      j = 1
      l = 1
      do 20 jj=1,idim
        l1 = j
        do 10 ii=1,k1
          cc(l1) = cp(l)
          l1 = l1+1
          l = l+1
  10    continue
        j = j+n
        l = l+k1
  20  continue
      if(nk1.eq.k1) go to 70
      n1 = k1*2
      j = n
      l = n1
      do 30 i=1,k1
        t1(i) = t(i)
        t1(l) = t(j)
        l = l-1
        j = j-1
  30  continue
c  find the b-spline representation of the given polynomial curve
c  according to the given set of knots.
      nk2 = nk1-1
      do 60 l=k1,nk2
        l1 = l+1
        j = 1
        do 40 i=1,idim
          call fpinst(0,t1,n1,cc(j),k,t(l1),l,t2,n2,cc(j),n)
          j = j+n
  40    continue
        do 50 i=1,n2
          t1(i) = t2(i)
  50    continue
        n1 = n2
  60  continue
c  find the b-spline representation of the resulting curve.
  70  j = 1
      do 90 jj=1,idim
        l = j
        do 80 i=1,nk1
          c(l) = cc(l)+c(l)
          l = l+1
  80    continue
        j = j+n
  90  continue
      return
      end
