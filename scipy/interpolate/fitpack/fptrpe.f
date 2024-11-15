      recursive subroutine fptrpe(m,mm,idim,n,nr,sp,p,b,z,a,aa,q,right)
      implicit none
c  subroutine fptrpe reduces the (m+n-7) x (n-7) cyclic bandmatrix a
c  to upper triangular form and applies the same givens transformations
c  to the (m) x (mm) x (idim) matrix z to obtain the (n-7) x (mm) x
c  (idim) matrix q.
c  ..
c  ..scalar arguments..
      real*8 p
      integer m,mm,idim,n
c  ..array arguments..
      real*8 sp(m,4),b(n,5),z(m*mm*idim),a(n,5),aa(n,4),q((n-7)*mm*idim)
     *,
     * right(mm*idim)
      integer nr(m)
c  ..local scalars..
      real*8 co,pinv,piv,si,one
      integer i,irot,it,ii,i2,i3,j,jj,l,mid,nmd,m2,m3,
     * nrold,n4,number,n1,n7,n11,m1
      integer i1, ij,j1,jk,jper,l0,l1, ik
c  ..local arrays..
      real*8 h(5),h1(5),h2(4)
c  ..subroutine references..
c    fpgivs,fprota
c  ..
      one = 1
      if(p.gt.0.) pinv = one/p
      n4 = n-4
      n7 = n-7
      n11 = n-11
      mid = mm*idim
      m2 = m*mm
      m3 = n7*mm
      m1 = m-1
c  we determine the matrix (a) and then we reduce her to
c  upper triangular form (r) using givens rotations.
c  we apply the same transformations to the rows of matrix
c  z to obtain the (mm) x (n-7) matrix g.
c  we store matrix (r) into a and aa, g into q.
c  the n7 x n7 upper triangular matrix (r) has the form
c             | a1 '     |
c       (r) = |    ' a2  |
c             |  0 '     |
c  with (a2) a n7 x 4 matrix and (a1) a n11 x n11 upper
c  triangular matrix of bandwidth 5.
c  initialization.
      nmd = n7*mid
      do 50 i=1,nmd
        q(i) = 0.
  50  continue
      do 100 i=1,n4
        a(i,5) = 0.
        do 100 j=1,4
          a(i,j) = 0.
          aa(i,j) = 0.
 100  continue
      jper = 0
      nrold = 0
      do 760 it=1,m1
        number = nr(it)
 120    if(nrold.eq.number) go to 180
        if(p.le.0.) go to 740
c  fetch a new row of matrix (b).
        n1 = nrold+1
        do 140 j=1,5
          h(j) = b(n1,j)*pinv
 140    continue
c  find the appropriate row of q.
        do 160 j=1,mid
          right(j) = 0.
 160    continue
        go to 240
c  fetch a new row of matrix (sp)
 180    h(5) = 0.
        do 200 j=1,4
          h(j) = sp(it,j)
 200    continue
c  find the appropriate row of q.
        j = 0
        do 220 ii=1,idim
          l = (ii-1)*m2+(it-1)*mm
          do 220 jj=1,mm
            j = j+1
            l = l+1
            right(j) = z(l)
 220    continue
c  test whether there are non-zero values in the new row of (a)
c  corresponding to the b-splines n(j,*),j=n7+1,...,n4.
 240     if(nrold.lt.n11) go to 640
         if(jper.ne.0) go to 320
c  initialize the matrix (aa).
         jk = n11+1
         do 300 i=1,4
            ik = jk
            do 260 j=1,5
               if(ik.le.0) go to 280
               aa(ik,i) = a(ik,j)
               ik = ik-1
 260        continue
 280        jk = jk+1
 300     continue
         jper = 1
c  if one of the non-zero elements of the new row corresponds to one of
c  the b-splines n(j;*),j=n7+1,...,n4,we take account of the periodicity
c  conditions for setting up this row of (a).
 320     do 340 i=1,4
            h1(i) = 0.
            h2(i) = 0.
 340     continue
         h1(5) = 0.
         j = nrold-n11
         do 420 i=1,5
            j = j+1
            l0 = j
 360        l1 = l0-4
            if(l1.le.0) go to 400
            if(l1.le.n11) go to 380
            l0 = l1-n11
            go to 360
 380        h1(l1) = h(i)
            go to 420
 400        h2(l0) = h2(l0) + h(i)
 420     continue
c  rotate the new row of (a) into triangle.
         if(n11.le.0) go to 560
c  rotations with the rows 1,2,...,n11 of (a).
         do 540 irot=1,n11
            piv = h1(1)
            i2 = min0(n11-irot,4)
            if(piv.eq.0.) go to 500
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(irot,1),co,si)
c  apply that transformation to the columns of matrix q.
            j = 0
            do 440 ii=1,idim
               l = (ii-1)*m3+irot
               do 440 jj=1,mm
                 j = j+1
                 call fprota(co,si,right(j),q(l))
                 l = l+n7
 440        continue
c  apply that transformation to the rows of (a) with respect to aa.
            do 460 i=1,4
               call fprota(co,si,h2(i),aa(irot,i))
 460        continue
c  apply that transformation to the rows of (a) with respect to a.
            if(i2.eq.0) go to 560
            do 480 i=1,i2
               i1 = i+1
               call fprota(co,si,h1(i1),a(irot,i1))
 480        continue
 500        do 520 i=1,i2
               h1(i) = h1(i+1)
 520        continue
            h1(i2+1) = 0.
 540     continue
c  rotations with the rows n11+1,...,n7 of a.
 560     do 620 irot=1,4
            ij = n11+irot
            if(ij.le.0) go to 620
            piv = h2(irot)
            if(piv.eq.0.) go to 620
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,aa(ij,irot),co,si)
c  apply that transformation to the columns of matrix q.
            j = 0
            do 580 ii=1,idim
               l = (ii-1)*m3+ij
               do 580 jj=1,mm
                 j = j+1
                 call fprota(co,si,right(j),q(l))
                 l = l+n7
 580        continue
            if(irot.eq.4) go to 620
c  apply that transformation to the rows of (a) with respect to aa.
            j1 = irot+1
            do 600 i=j1,4
               call fprota(co,si,h2(i),aa(ij,i))
 600        continue
 620     continue
         go to 720
c  rotation into triangle of the new row of (a), in case the elements
c  corresponding to the b-splines n(j;*),j=n7+1,...,n4 are all zero.
 640     irot =nrold
         do 700 i=1,5
            irot = irot+1
            piv = h(i)
            if(piv.eq.0.) go to 700
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(irot,1),co,si)
c  apply that transformation to the columns of matrix g.
            j = 0
            do 660 ii=1,idim
               l = (ii-1)*m3+irot
               do 660 jj=1,mm
                 j = j+1
                 call fprota(co,si,right(j),q(l))
                 l = l+n7
 660        continue
c  apply that transformation to the rows of (a).
            if(i.eq.5) go to 700
            i2 = 1
            i3 = i+1
            do 680 j=i3,5
               i2 = i2+1
               call fprota(co,si,h(j),a(irot,i2))
 680        continue
 700     continue
 720     if(nrold.eq.number) go to 760
 740     nrold = nrold+1
         go to 120
 760  continue
      return
      end
