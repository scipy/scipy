      subroutine nnfc
     *     (n, r,c,ic, ia,ja,a, z, b,
     *      lmax,il,jl,ijl,l, d, umax,iu,ju,iju,u,
     *      row, tmp, irl,jrl, flag)
clll. optimize
c*** subroutine nnfc
c*** numerical ldu-factorization of sparse nonsymmetric matrix and
c      solution of system of linear equations (compressed pointer
c      storage)
c
c
c       input variables..  n, r, c, ic, ia, ja, a, b,
c                          il, jl, ijl, lmax, iu, ju, iju, umax
c       output variables.. z, l, d, u, flag
c
c       parameters used internally..
c nia   - irl,  - vectors used to find the rows of  l.  at the kth step
c nia   - jrl       of the factorization,  jrl(k)  points to the head
c       -           of a linked list in  jrl  of column indices j
c       -           such j .lt. k and  l(k,j)  is nonzero.  zero
c       -           indicates the end of the list.  irl(j)  (j.lt.k)
c       -           points to the smallest i such that i .ge. k and
c       -           l(i,j)  is nonzero.
c       -           size of each = n.
c fia   - row   - holds intermediate values in calculation of  u and l.
c       -           size = n.
c fia   - tmp   - holds new right-hand side  b*  for solution of the
c       -           equation ux = b*.
c       -           size = n.
c
c  internal variables..
c    jmin, jmax - indices of the first and last positions in a row to
c      be examined.
c    sum - used in calculating  tmp.
c
      integer rk,umax
      integer  r(1), c(1), ic(1), ia(1), ja(1), il(1), jl(1), ijl(1)
      integer  iu(1), ju(1), iju(1), irl(1), jrl(1), flag
      double precision  a(1), l(1), d(1), u(1), z(1), b(1), row(1)
      double precision  tmp(1), lki, sum, dk
c
c  ******  initialize pointers and test storage  ***********************
      if(il(n+1)-1 .gt. lmax) go to 104
      if(iu(n+1)-1 .gt. umax) go to 107
      do 1 k=1,n
        irl(k) = il(k)
        jrl(k) = 0
   1    continue
c
c  ******  for each row  ***********************************************
      do 19 k=1,n
c  ******  reverse jrl and zero row where kth row of l will fill in  ***
        row(k) = 0
        i1 = 0
        if (jrl(k) .eq. 0) go to 3
        i = jrl(k)
   2    i2 = jrl(i)
        jrl(i) = i1
        i1 = i
        row(i) = 0
        i = i2
        if (i .ne. 0) go to 2
c  ******  set row to zero where u will fill in  ***********************
   3    jmin = iju(k)
        jmax = jmin + iu(k+1) - iu(k) - 1
        if (jmin .gt. jmax) go to 5
        do 4 j=jmin,jmax
   4      row(ju(j)) = 0
c  ******  place kth row of a in row  **********************************
   5    rk = r(k)
        jmin = ia(rk)
        jmax = ia(rk+1) - 1
        do 6 j=jmin,jmax
          row(ic(ja(j))) = a(j)
   6      continue
c  ******  initialize sum, and link through jrl  ***********************
        sum = b(rk)
        i = i1
        if (i .eq. 0) go to 10
c  ******  assign the kth row of l and adjust row, sum  ****************
   7      lki = -row(i)
c  ******  if l is not required, then comment out the following line  **
          l(irl(i)) = -lki
          sum = sum + lki * tmp(i)
          jmin = iu(i)
          jmax = iu(i+1) - 1
          if (jmin .gt. jmax) go to 9
          mu = iju(i) - jmin
          do 8 j=jmin,jmax
   8        row(ju(mu+j)) = row(ju(mu+j)) + lki * u(j)
   9      i = jrl(i)
          if (i .ne. 0) go to 7
c
c  ******  assign kth row of u and diagonal d, set tmp(k)  *************
  10    if (row(k) .eq. 0.0d0) go to 108
        dk = 1.0d0 / row(k)
        d(k) = dk
        tmp(k) = sum * dk
        if (k .eq. n) go to 19
        jmin = iu(k)
        jmax = iu(k+1) - 1
        if (jmin .gt. jmax)  go to 12
        mu = iju(k) - jmin
        do 11 j=jmin,jmax
  11      u(j) = row(ju(mu+j)) * dk
  12    continue
c
c  ******  update irl and jrl, keeping jrl in decreasing order  ********
        i = i1
        if (i .eq. 0) go to 18
  14    irl(i) = irl(i) + 1
        i1 = jrl(i)
        if (irl(i) .ge. il(i+1)) go to 17
        ijlb = irl(i) - il(i) + ijl(i)
        j = jl(ijlb)
  15    if (i .gt. jrl(j)) go to 16
          j = jrl(j)
          go to 15
  16    jrl(i) = jrl(j)
        jrl(j) = i
  17    i = i1
        if (i .ne. 0) go to 14
  18    if (irl(k) .ge. il(k+1)) go to 19
        j = jl(ijl(k))
        jrl(k) = jrl(j)
        jrl(j) = k
  19    continue
c
c  ******  solve  ux = tmp  by back substitution  **********************
      k = n
      do 22 i=1,n
        sum =  tmp(k)
        jmin = iu(k)
        jmax = iu(k+1) - 1
        if (jmin .gt. jmax)  go to 21
        mu = iju(k) - jmin
        do 20 j=jmin,jmax
  20      sum = sum - u(j) * tmp(ju(mu+j))
  21    tmp(k) =  sum
        z(c(k)) =  sum
  22    k = k-1
      flag = 0
      return
c
c ** error.. insufficient storage for l
 104  flag = 4*n + 1
      return
c ** error.. insufficient storage for u
 107  flag = 7*n + 1
      return
c ** error.. zero pivot
 108  flag = 8*n + k
      return
      end
