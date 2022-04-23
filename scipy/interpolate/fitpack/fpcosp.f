      subroutine fpcosp(m,x,y,w,n,t,e,maxtr,maxbin,c,sq,sx,bind,nm,mb,a,
     *
     * b,const,z,zz,u,q,info,up,left,right,jbind,ibind,ier)
c  ..
c  ..scalar arguments..
      real*8 sq
      integer m,n,maxtr,maxbin,nm,mb,ier
c  ..array arguments..
      real*8 x(m),y(m),w(m),t(n),e(n),c(n),sx(m),a(n,4),b(nm,maxbin),
     * const(n),z(n),zz(n),u(maxbin),q(m,4)
      integer info(maxtr),up(maxtr),left(maxtr),right(maxtr),jbind(mb),
     * ibind(mb)
      logical bind(n)
c  ..local scalars..
      integer count,i,i1,j,j1,j2,j3,k,kdim,k1,k2,k3,k4,k5,k6,
     * l,lp1,l1,l2,l3,merk,nbind,number,n1,n4,n6
      real*8 f,wi,xi
c  ..local array..
      real*8 h(4)
c  ..subroutine references..
c    fpbspl,fpadno,fpdeno,fpfrno,fpseno
c  ..
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  if we use the b-spline representation of s(x) our approximation     c
c  problem results in a quadratic programming problem:                 c
c    find the b-spline coefficients c(j),j=1,2,...n-4 such that        c
c        (1) sumi((wi*(yi-sumj(cj*nj(xi))))**2),i=1,2,...m is minimal  c
c        (2) sumj(cj*n''j(t(l+3)))*e(l) <= 0, l=1,2,...n-6.            c
c  to solve this problem we use the theil-van de panne procedure.      c
c  if the inequality constraints (2) are numbered from 1 to n-6,       c
c  this algorithm finds a subset of constraints ibind(1)..ibind(nbind) c
c  such that the solution of the minimization problem (1) with these   c
c  constraints in equality form, satisfies all constraints. such a     c
c  feasible solution is optimal if the lagrange parameters associated  c
c  with that problem with equality constraints, are all positive.      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  determine n6, the number of inequality constraints.
      n6 = n-6
c  fix the parameters which determine these constraints.
      do 10 i=1,n6
        const(i) = e(i)*(t(i+4)-t(i+1))/(t(i+5)-t(i+2))
  10  continue
c  initialize the triply linked tree which is used to find the subset
c  of constraints ibind(1),...ibind(nbind).
      count = 1
      info(1) = 0
      left(1) = 0
      right(1) = 0
      up(1) = 1
      merk = 1
c  set up the normal equations  n'nc=n'y  where n denotes the m x (n-4)
c  observation matrix with elements ni,j = wi*nj(xi)  and y is the
c  column vector with elements yi*wi.
c  from the properties of the b-splines nj(x),j=1,2,...n-4, it follows
c  that  n'n  is a (n-4) x (n-4)  positive definite bandmatrix of
c  bandwidth 7. the matrices n'n and n'y are built up in a and z.
      n4 = n-4
c  initialization
      do 20 i=1,n4
        z(i) = 0.
        do 20 j=1,4
          a(i,j) = 0.
  20  continue
      l = 4
      lp1 = l+1
      do 70 i=1,m
c  fetch the current row of the observation matrix.
        xi = x(i)
        wi = w(i)**2
c  search for knot interval  t(l) <= xi < t(l+1)
  30    if(xi.lt.t(lp1) .or. l.eq.n4) go to 40
        l = lp1
        lp1 = l+1
        go to 30
c  evaluate the four non-zero cubic b-splines nj(xi),j=l-3,...l.
  40    call fpbspl(t,n,3,xi,l,h)
c  store in q these values h(1),h(2),...h(4).
        do 50 j=1,4
          q(i,j) = h(j)
  50    continue
c  add the contribution of the current row of the observation matrix
c  n to the normal equations.
        l3 = l-3
        k1 = 0
        do 60 j1 = l3,l
          k1 = k1+1
          f = h(k1)
          z(j1) = z(j1)+f*wi*y(i)
          k2 = k1
          j2 = 4
          do 60 j3 = j1,l
            a(j3,j2) = a(j3,j2)+f*wi*h(k2)
            k2 = k2+1
            j2 = j2-1
  60    continue
  70  continue
c  since n'n is a symmetric matrix it can be factorized as
c        (3)  n'n = (r1)'(d1)(r1)
c  with d1 a diagonal matrix and r1 an (n-4) x (n-4)  unit upper
c  triangular matrix of bandwidth 4. the matrices r1 and d1 are built
c  up in a. at the same time we solve the systems of equations
c        (4)  (r1)'(z2) = n'y
c        (5)  (d1) (z1) = (z2)
c  the vectors z2 and z1 are kept in zz and z.
      do 140 i=1,n4
        k1 = 1
        if(i.lt.4) k1 = 5-i
        k2 = i-4+k1
        k3 = k2
        do 100 j=k1,4
          k4 = j-1
          k5 = 4-j+k1
          f = a(i,j)
          if(k1.gt.k4) go to 90
          k6 = k2
          do 80 k=k1,k4
            f = f-a(i,k)*a(k3,k5)*a(k6,4)
            k5 = k5+1
            k6 = k6+1
  80      continue
  90      if(j.eq.4) go to 110
          a(i,j) = f/a(k3,4)
          k3 = k3+1
 100    continue
 110    a(i,4) = f
        f = z(i)
        if(i.eq.1) go to 130
        k4 = i
        do 120 j=k1,3
          k = k1+3-j
          k4 = k4-1
          f = f-a(i,k)*z(k4)*a(k4,4)
 120    continue
 130    z(i) = f/a(i,4)
        zz(i) = f
 140  continue
c  start computing the least-squares cubic spline without taking account
c  of any constraint.
      nbind = 0
      n1 = 1
      ibind(1) = 0
c  main loop for the least-squares problems with different subsets of
c  the constraints (2) in equality form. the resulting b-spline coeff.
c  c and lagrange parameters u are the solution of the system
c            ! n'n  b' ! ! c !   ! n'y !
c        (6) !         ! !   ! = !     !
c            !  b   0  ! ! u !   !  0  !
c  z1 is stored into array c.
 150  do 160 i=1,n4
        c(i) = z(i)
 160  continue
c  if there are no equality constraints, compute the coeff. c directly.
      if(nbind.eq.0) go to 370
c  initialization
      kdim = n4+nbind
      do 170 i=1,nbind
        do 170 j=1,kdim
          b(j,i) = 0.
 170  continue
c  matrix b is built up,expressing that the constraints nrs ibind(1),...
c  ibind(nbind) must be satisfied in equality form.
      do 180 i=1,nbind
        l = ibind(i)
        b(l,i) = e(l)
        b(l+1,i) = -(e(l)+const(l))
        b(l+2,i) = const(l)
 180  continue
c  find the matrix (b1) as the solution of the system of equations
c        (7)  (r1)'(d1)(b1) = b'
c  (b1) is built up in the upper part of the array b(rows 1,...n-4).
      do 220 k1=1,nbind
        l = ibind(k1)
        do 210 i=l,n4
          f = b(i,k1)
          if(i.eq.1) go to 200
          k2 = 3
          if(i.lt.4) k2 = i-1
          do 190 k3=1,k2
            l1 = i-k3
            l2 = 4-k3
            f = f-b(l1,k1)*a(i,l2)*a(l1,4)
 190      continue
 200      b(i,k1) = f/a(i,4)
 210    continue
 220  continue
c  factorization of the symmetric matrix  -(b1)'(d1)(b1)
c        (8)  -(b1)'(d1)(b1) = (r2)'(d2)(r2)
c  with (d2) a diagonal matrix and (r2) an nbind x nbind unit upper
c  triangular matrix. the matrices r2 and d2 are built up in the lower
c  part of the array b (rows n-3,n-2,...n-4+nbind).
      do 270 i=1,nbind
        i1 = i-1
        do 260 j=i,nbind
          f = 0.
          do 230 k=1,n4
            f = f+b(k,i)*b(k,j)*a(k,4)
 230      continue
          k1 = n4+1
          if(i1.eq.0) go to 250
          do 240 k=1,i1
            f = f+b(k1,i)*b(k1,j)*b(k1,k)
            k1 = k1+1
 240      continue
 250      b(k1,j) = -f
          if(j.eq.i) go to 260
          b(k1,j) = b(k1,j)/b(k1,i)
 260    continue
 270  continue
c  according to (3),(7) and (8) the system of equations (6) becomes
c         ! (r1)'    0  ! ! (d1)    0  ! ! (r1)  (b1) ! ! c !   ! n'y !
c    (9)  !             ! !            ! !            ! !   ! = !     !
c         ! (b1)'  (r2)'! !   0   (d2) ! !   0   (r2) ! ! u !   !  0  !
c  backward substitution to obtain the b-spline coefficients c(j),j=1,..
c  n-4 and the lagrange parameters u(j),j=1,2,...nbind.
c  first step of the backward substitution: solve the system
c             ! (r1)'(d1)      0     ! ! (c1) !   ! n'y !
c        (10) !                      ! !      ! = !     !
c             ! (b1)'(d1)  (r2)'(d2) ! ! (u1) !   !  0  !
c  from (4) and (5) we know that this is equivalent to
c        (11)  (c1) = (z1)
c        (12)  (r2)'(d2)(u1) = -(b1)'(z2)
      do 310 i=1,nbind
        f = 0.
        do 280 j=1,n4
          f = f+b(j,i)*zz(j)
 280    continue
        i1 = i-1
        k1 = n4+1
        if(i1.eq.0) go to 300
        do 290 j=1,i1
          f = f+u(j)*b(k1,i)*b(k1,j)
          k1 = k1+1
 290    continue
 300    u(i) = -f/b(k1,i)
 310  continue
c  second step of the backward substitution: solve the system
c             ! (r1)  (b1) ! ! c !   ! c1 !
c        (13) !            ! !   ! = !    !
c             !   0   (r2) ! ! u !   ! u1 !
      k1 = nbind
      k2 = kdim
c  find the lagrange parameters u.
      do 340 i=1,nbind
        f = u(k1)
        if(i.eq.1) go to 330
        k3 = k1+1
        do 320 j=k3,nbind
          f = f-u(j)*b(k2,j)
 320    continue
 330    u(k1) = f
        k1 = k1-1
        k2 = k2-1
 340  continue
c  find the b-spline coefficients c.
      do 360 i=1,n4
        f = c(i)
        do 350 j=1,nbind
          f = f-u(j)*b(i,j)
 350    continue
        c(i) = f
 360  continue
 370  k1 = n4
      do 390 i=2,n4
        k1 = k1-1
        f = c(k1)
        k2 = 1
        if(i.lt.5) k2 = 5-i
        k3 = k1
        l = 3
        do 380 j=k2,3
          k3 = k3+1
          f = f-a(k3,l)*c(k3)
          l = l-1
 380    continue
        c(k1) = f
 390  continue
c  test whether the solution of the least-squares problem with the
c  constraints ibind(1),...ibind(nbind) in equality form, satisfies
c  all of the constraints (2).
      k = 1
c  number counts the number of violated inequality constraints.
      number = 0
      do 440 j=1,n6
        l = ibind(k)
        k = k+1
        if(j.eq.l) go to 440
        k = k-1
c  test whether constraint j is satisfied
        f = e(j)*(c(j)-c(j+1))+const(j)*(c(j+2)-c(j+1))
        if(f.le.0.) go to 440
c  if constraint j is not satisfied, add a branch of length nbind+1
c  to the tree. the nodes of this branch contain in their information
c  field the number of the constraints ibind(1),...ibind(nbind) and j,
c  arranged in increasing order.
        number = number+1
        k1 = k-1
        if(k1.eq.0) go to 410
        do 400 i=1,k1
          jbind(i) = ibind(i)
 400    continue
 410    jbind(k) = j
        if(l.eq.0) go to 430
        do 420 i=k,nbind
          jbind(i+1) = ibind(i)
 420    continue
 430    call fpadno(maxtr,up,left,right,info,count,merk,jbind,n1,ier)
c  test whether the storage space which is required for the tree,exceeds
c  the available storage space.
        if(ier.ne.0) go to 560
 440  continue
c  test whether the solution of the least-squares problem with equality
c  constraints is a feasible solution.
      if(number.eq.0) go to 470
c  test whether there are still cases with nbind constraints in
c  equality form to be considered.
 450  if(merk.gt.1) go to 460
      nbind = n1
c  test whether the number of knots where s''(x)=0 exceeds maxbin.
      if(nbind.gt.maxbin) go to 550
      n1 = n1+1
      ibind(n1) = 0
c  search which cases with nbind constraints in equality form
c  are going to be considered.
      call fpdeno(maxtr,up,left,right,nbind,merk)
c  test whether the quadratic programming problem has a solution.
      if(merk.eq.1) go to 570
c  find a new case with nbind constraints in equality form.
 460  call fpseno(maxtr,up,left,right,info,merk,ibind,nbind)
      go to 150
c  test whether the feasible solution is optimal.
 470  ier = 0
      do 480 i=1,n6
        bind(i) = .false.
 480  continue
      if(nbind.eq.0) go to 500
      do 490 i=1,nbind
        if(u(i).le.0.) go to 450
        j = ibind(i)
        bind(j) = .true.
 490  continue
c  evaluate s(x) at the data points x(i) and calculate the weighted
c  sum of squared residual right hand sides sq.
 500  sq = 0.
      l = 4
      lp1 = 5
      do 530 i=1,m
 510    if(x(i).lt.t(lp1) .or. l.eq.n4) go to 520
        l = lp1
        lp1 = l+1
        go to 510
 520    sx(i) = c(l-3)*q(i,1)+c(l-2)*q(i,2)+c(l-1)*q(i,3)+c(l)*q(i,4)
        sq = sq+(w(i)*(y(i)-sx(i)))**2
 530  continue
      go to 600
c  error codes and messages.
 550  ier = 1
      go to 600
 560  ier = 2
      go to 600
 570  ier = 3
 600  return
      end
