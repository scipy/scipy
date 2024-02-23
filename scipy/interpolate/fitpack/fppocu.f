      recursive subroutine fppocu(idim,k,a,b,ib,db,nb,ie,de,ne,cp,np)
      implicit none
c  subroutine fppocu finds a idim-dimensional polynomial curve p(u) =
c  (p1(u),p2(u),...,pidim(u)) of degree k, satisfying certain derivative
c  constraints at the end points a and b, i.e.
c                  (l)
c    if ib > 0 : pj   (a) = db(idim*l+j), l=0,1,...,ib-1
c                  (l)
c    if ie > 0 : pj   (b) = de(idim*l+j), l=0,1,...,ie-1
c
c  the polynomial curve is returned in its b-spline representation
c  ( coefficients cp(j), j=1,2,...,np )
c  ..
c  ..scalar arguments..
      integer idim,k,ib,nb,ie,ne,np
      real*8 a,b
c  ..array arguments..
      real*8 db(nb),de(ne),cp(np)
c  ..local scalars..
      real*8 ab,aki
      integer i,id,j,jj,l,ll,k1,k2
c  ..local array..
      real*8 work(6,6)
c  ..
      k1 = k+1
      k2 = 2*k1
      ab = b-a
      do 110 id=1,idim
        do 10 j=1,k1
          work(j,1) = 0.
  10    continue
        if(ib.eq.0) go to 50
        l = id
        do 20 i=1,ib
          work(1,i) = db(l)
          l = l+idim
  20    continue
        if(ib.eq.1) go to 50
        ll = ib
        do 40 j=2,ib
          ll =  ll-1
          do 30 i=1,ll
            aki = k1-i
            work(j,i) = ab*work(j-1,i+1)/aki + work(j-1,i)
  30      continue
  40    continue
  50    if(ie.eq.0) go to 90
        l = id
        j = k1
        do 60 i=1,ie
          work(j,i) = de(l)
          l = l+idim
          j = j-1
  60    continue
        if(ie.eq.1) go to 90
        ll = ie
        do 80 jj=2,ie
          ll =  ll-1
          j = k1+1-jj
          do 70 i=1,ll
            aki = k1-i
            work(j,i) = work(j+1,i) - ab*work(j,i+1)/aki
            j = j-1
  70      continue
  80    continue
  90    l = (id-1)*k2
        do 100 j=1,k1
          l = l+1
          cp(l) = work(j,1)
 100    continue
 110  continue
      return
      end
