      recursive subroutine sproot(t,n,c,zero,mest,m,ier)
      implicit none
c  subroutine sproot finds the zeros of a cubic spline s(x),which is
c  given in its normalized b-spline representation.
c
c  calling sequence:
c     call sproot(t,n,c,zero,mest,m,ier)
c
c  input parameters:
c    t    : real array,length n, containing the knots of s(x).
c    n    : integer, containing the number of knots.  n>=8
c    c    : real array,length n, containing the b-spline coefficients.
c    mest : integer, specifying the dimension of array zero.
c
c  output parameters:
c    zero : real array,length mest, containing the zeros of s(x).
c    m    : integer,giving the number of zeros.
c    ier  : error flag:
c      ier = 0: normal return.
c      ier = 1: the number of zeros exceeds mest.
c      ier =10: invalid input data (see restrictions).
c
c  other subroutines required: fpcuro
c
c  restrictions:
c    1) n>= 8.
c    2) t(4) < t(5) < ... < t(n-4) < t(n-3).
c       t(1) <= t(2) <= t(3) <= t(4)
c       t(n-3) <= t(n-2) <= t(n-1) <= t(n)
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c ..
c ..scalar arguments..
      integer n,mest,m,ier
c  ..array arguments..
      real*8 t(n),c(n),zero(mest)
c  ..local scalars..
      integer i,j,j1,l,n4
      real*8 ah,a0,a1,a2,a3,bh,b0,b1,c1,c2,c3,c4,c5,d4,d5,h1,h2,
     * three,two,t1,t2,t3,t4,t5,zz
      logical z0,z1,z2,z3,z4,nz0,nz1,nz2,nz3,nz4
c  ..local array..
      real*8 y(3)
c  ..
c  set some constants
      two = 0.2d+01
      three = 0.3d+01
c  before starting computations a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      n4 = n-4
      ier = 10
      if(n.lt.8) go to 800
      j = n
      do 10 i=1,3
        if(t(i).gt.t(i+1)) go to 800
        if(t(j).lt.t(j-1)) go to 800
        j = j-1
  10  continue
      do 20 i=4,n4
        if(t(i).ge.t(i+1)) go to 800
  20  continue
c  the problem considered reduces to finding the zeros of the cubic
c  polynomials pl(x) which define the cubic spline in each knot
c  interval t(l)<=x<=t(l+1). a zero of pl(x) is also a zero of s(x) on
c  the condition that it belongs to the knot interval.
c  the cubic polynomial pl(x) is determined by computing s(t(l)),
c  s'(t(l)),s(t(l+1)) and s'(t(l+1)). in fact we only have to compute
c  s(t(l+1)) and s'(t(l+1)); because of the continuity conditions of
c  splines and their derivatives, the value of s(t(l)) and s'(t(l))
c  is already known from the foregoing knot interval.
      ier = 0
c  evaluate some constants for the first knot interval
      h1 = t(4)-t(3)
      h2 = t(5)-t(4)
      t1 = t(4)-t(2)
      t2 = t(5)-t(3)
      t3 = t(6)-t(4)
      t4 = t(5)-t(2)
      t5 = t(6)-t(3)
c  calculate a0 = s(t(4)) and ah = s'(t(4)).
      c1 = c(1)
      c2 = c(2)
      c3 = c(3)
      c4 = (c2-c1)/t4
      c5 = (c3-c2)/t5
      d4 = (h2*c1+t1*c2)/t4
      d5 = (t3*c2+h1*c3)/t5
      a0 = (h2*d4+h1*d5)/t2
      ah = three*(h2*c4+h1*c5)/t2
      z1 = .true.
      if(ah.lt.0.0d0) z1 = .false.
      nz1 = .not.z1
      m = 0
c  main loop for the different knot intervals.
      do 300 l=4,n4
c  evaluate some constants for the knot interval t(l) <= x <= t(l+1).
        h1 = h2
        h2 = t(l+2)-t(l+1)
        t1 = t2
        t2 = t3
        t3 = t(l+3)-t(l+1)
        t4 = t5
        t5 = t(l+3)-t(l)
c  find a0 = s(t(l)), ah = s'(t(l)), b0 = s(t(l+1)) and bh = s'(t(l+1)).
        c1 = c2
        c2 = c3
        c3 = c(l)
        c4 = c5
        c5 = (c3-c2)/t5
        d4 = (h2*c1+t1*c2)/t4
        d5 = (h1*c3+t3*c2)/t5
        b0 = (h2*d4+h1*d5)/t2
        bh = three*(h2*c4+h1*c5)/t2
c  calculate the coefficients a0,a1,a2 and a3 of the cubic polynomial
c  pl(x) = ql(y) = a0+a1*y+a2*y**2+a3*y**3 ; y = (x-t(l))/(t(l+1)-t(l)).
        a1 = ah*h1
        b1 = bh*h1
        a2 = three*(b0-a0)-b1-two*a1
        a3 = two*(a0-b0)+b1+a1
c  test whether or not pl(x) could have a zero in the range
c  t(l) <= x <= t(l+1).
        z3 = .true.
        if(b1.lt.0.0d0) z3 = .false.
        nz3 = .not.z3
        if(a0*b0.le.0.0d0) go to 100
        z0 = .true.
        if(a0.lt.0.0d0) z0 = .false.
        nz0 = .not.z0
        z2 = .true.
        if(a2.lt.0.) z2 = .false.
        nz2 = .not.z2
        z4 = .true.
        if(3.0d0*a3+a2.lt.0.0d0) z4 = .false.
        nz4 = .not.z4
        if(.not.((z0.and.(nz1.and.(z3.or.z2.and.nz4).or.nz2.and.
     * z3.and.z4).or.nz0.and.(z1.and.(nz3.or.nz2.and.z4).or.z2.and.
     * nz3.and.nz4))))go to 200
c  find the zeros of ql(y).
 100    call fpcuro(a3,a2,a1,a0,y,j)
        if(j.eq.0) go to 200
c  find which zeros of pl(x) are zeros of s(x).
        do 150 i=1,j
          if(y(i).lt.0.0d0 .or. y(i).gt.1.0d0) go to 150
c  test whether the number of zeros of s(x) exceeds mest.
          if(m.ge.mest) go to 700
          m = m+1
          zero(m) = t(l)+h1*y(i)
 150    continue
 200    a0 = b0
        ah = bh
        z1 = z3
        nz1 = nz3
 300  continue
c  the zeros of s(x) are arranged in increasing order.
      if(m.lt.2) go to 800
      do 400 i=2,m
        j = i
 350    j1 = j-1
        if(j1.eq.0) go to 400
        if(zero(j).ge.zero(j1)) go to 400
        zz = zero(j)
        zero(j) = zero(j1)
        zero(j1) = zz
        j = j1
        go to 350
 400  continue
      j = m
      m = 1
      do 500 i=2,j
        if(zero(i).eq.zero(m)) go to 500
        m = m+1
        zero(m) = zero(i)
 500  continue
      go to 800
 700  ier = 1
 800  return
      end
