      subroutine fprank(a,f,n,m,na,tol,c,sq,rank,aa,ff,h)
c  subroutine fprank finds the minimum norm solution of a least-
c  squares problem in case of rank deficiency.
c
c  input parameters:
c    a : array, which contains the non-zero elements of the observation
c        matrix after triangularization by givens transformations.
c    f : array, which contains the transformed right hand side.
c    n : integer,which contains the dimension of a.
c    m : integer, which denotes the bandwidth of a.
c  tol : real value, giving a threshold to determine the rank of a.
c
c  output parameters:
c    c : array, which contains the minimum norm solution.
c   sq : real value, giving the contribution of reducing the rank
c        to the sum of squared residuals.
c rank : integer, which contains the rank of matrix a.
c
c  ..scalar arguments..
      integer n,m,na,rank
      real*8 tol,sq
c  ..array arguments..
      real*8 a(na,m),f(n),c(n),aa(n,m),ff(n),h(m)
c  ..local scalars..
      integer i,ii,ij,i1,i2,j,jj,j1,j2,j3,k,kk,m1,nl
      real*8 cos,fac,piv,sin,yi
      double precision store,stor1,stor2,stor3
c  ..function references..
      integer min0
c  ..subroutine references..
c    fpgivs,fprota
c  ..
      m1 = m-1
c  the rank deficiency nl is considered to be the number of sufficient
c  small diagonal elements of a.
      nl = 0
      sq = 0.
      do 90 i=1,n
        if(a(i,1).gt.tol) go to 90
c  if a sufficient small diagonal element is found, we put it to
c  zero. the remainder of the row corresponding to that zero diagonal
c  element is then rotated into triangle by givens rotations .
c  the rank deficiency is increased by one.
        nl = nl+1
        if(i.eq.n) go to 90
        yi = f(i)
        do 10 j=1,m1
          h(j) = a(i,j+1)
  10    continue
        h(m) = 0.
        i1 = i+1
        do 60 ii=i1,n
          i2 = min0(n-ii,m1)
          piv = h(1)
          if(piv.eq.0.) go to 30
          call fpgivs(piv,a(ii,1),cos,sin)
          call fprota(cos,sin,yi,f(ii))
          if(i2.eq.0) go to 70
          do 20 j=1,i2
            j1 = j+1
            call fprota(cos,sin,h(j1),a(ii,j1))
            h(j) = h(j1)
  20      continue
          go to 50
  30      if(i2.eq.0) go to 70
          do 40 j=1,i2
            h(j) = h(j+1)
  40      continue
  50      h(i2+1) = 0.
  60    continue
c  add to the sum of squared residuals the contribution of deleting
c  the row with small diagonal element.
  70    sq = sq+yi**2
  90  continue
c  rank denotes the rank of a.
      rank = n-nl
c  let b denote the (rank*n) upper trapezoidal matrix which can be
c  obtained from the (n*n) upper triangular matrix a by deleting
c  the rows and interchanging the columns corresponding to a zero
c  diagonal element. if this matrix is factorized using givens
c  transformations as  b = (r) (u)  where
c    r is a (rank*rank) upper triangular matrix,
c    u is a (rank*n) orthonormal matrix
c  then the minimal least-squares solution c is given by c = b' v,
c  where v is the solution of the system  (r) (r)' v = g  and
c  g denotes the vector obtained from the old right hand side f, by
c  removing the elements corresponding to a zero diagonal element of a.
c  initialization.
      do 100 i=1,rank
        do 100 j=1,m
          aa(i,j) = 0.
 100  continue
c  form in aa the upper triangular matrix obtained from a by
c  removing rows and columns with zero diagonal elements. form in ff
c  the new right hand side by removing the elements of the old right
c  hand side corresponding to a deleted row.
      ii = 0
      do 120 i=1,n
        if(a(i,1).le.tol) go to 120
        ii = ii+1
        ff(ii) = f(i)
        aa(ii,1) = a(i,1)
        jj = ii
        kk = 1
        j = i
        j1 = min0(j-1,m1)
        if(j1.eq.0) go to 120
        do 110 k=1,j1
          j = j-1
          if(a(j,1).le.tol) go to 110
          kk = kk+1
          jj = jj-1
          aa(jj,kk) = a(j,k+1)
 110    continue
 120  continue
c  form successively in h the columns of a with a zero diagonal element.
      ii = 0
      do 200 i=1,n
        ii = ii+1
        if(a(i,1).gt.tol) go to 200
        ii = ii-1
        if(ii.eq.0) go to 200
        jj = 1
        j = i
        j1 = min0(j-1,m1)
        do 130 k=1,j1
          j = j-1
          if(a(j,1).le.tol) go to 130
          h(jj) = a(j,k+1)
          jj = jj+1
 130    continue
        do 140 kk=jj,m
          h(kk) = 0.
 140    continue
c  rotate this column into aa by givens transformations.
        jj = ii
        do 190 i1=1,ii
          j1 = min0(jj-1,m1)
          piv = h(1)
          if(piv.ne.0.) go to 160
          if(j1.eq.0) go to 200
          do 150 j2=1,j1
            j3 = j2+1
            h(j2) = h(j3)
 150      continue
          go to 180
 160      call fpgivs(piv,aa(jj,1),cos,sin)
          if(j1.eq.0) go to 200
          kk = jj
          do 170 j2=1,j1
            j3 = j2+1
            kk = kk-1
            call fprota(cos,sin,h(j3),aa(kk,j3))
            h(j2) = h(j3)
 170      continue
 180      jj = jj-1
          h(j3) = 0.
 190    continue
 200  continue
c  solve the system (aa) (f1) = ff
      ff(rank) = ff(rank)/aa(rank,1)
      i = rank-1
      if(i.eq.0) go to 230
      do 220 j=2,rank
        store = ff(i)
        i1 = min0(j-1,m1)
        k = i
        do 210 ii=1,i1
          k = k+1
          stor1 = ff(k)
          stor2 = aa(i,ii+1)
          store = store-stor1*stor2
 210    continue
        stor1 = aa(i,1)
        ff(i) = store/stor1
        i = i-1
 220  continue
c  solve the system  (aa)' (f2) = f1
 230  ff(1) = ff(1)/aa(1,1)
      if(rank.eq.1) go to 260
      do 250 j=2,rank
        store = ff(j)
        i1 = min0(j-1,m1)
        k = j
        do 240 ii=1,i1
          k = k-1
          stor1 = ff(k)
          stor2 = aa(k,ii+1)
          store = store-stor1*stor2
 240    continue
        stor1 = aa(j,1)
        ff(j) = store/stor1
 250  continue
c  premultiply f2 by the transpoze of a.
 260  k = 0
      do 280 i=1,n
        store = 0.
        if(a(i,1).gt.tol) k = k+1
        j1 = min0(i,m)
        kk = k
        ij = i+1
        do 270 j=1,j1
          ij = ij-1
          if(a(ij,1).le.tol) go to 270
          stor1 = a(ij,j)
          stor2 = ff(kk)
          store = store+stor1*stor2
          kk = kk-1
 270    continue
        c(i) = store
 280  continue
c  add to the sum of squared residuals the contribution of putting
c  to zero the small diagonal elements of matrix (a).
      stor3 = 0.
      do 310 i=1,n
        if(a(i,1).gt.tol) go to 310
        store = f(i)
        i1 = min0(n-i,m1)
        if(i1.eq.0) go to 300
        do 290 j=1,i1
          ij = i+j
          stor1 = c(ij)
          stor2 = a(i,j+1)
          store = store-stor1*stor2
 290    continue
 300    fac = a(i,1)*c(i)
        stor1 = a(i,1)
        stor2 = c(i)
        stor1 = stor1*stor2
        stor3 = stor3+stor1*(stor1-store-store)
 310  continue
      fac = stor3
      sq = sq+fac
      return
      end

