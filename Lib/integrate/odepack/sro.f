      subroutine sro
     *     (n, ip, ia,ja,a, q, r, dflag)
clll. optimize
c***********************************************************************
c  sro -- symmetric reordering of sparse symmetric matrix
c***********************************************************************
c
c  description
c
c    the nonzero entries of the matrix m are assumed to be stored
c    symmetrically in (ia,ja,a) format (i.e., not both m(i,j) and m(j,i)
c    are stored if i ne j).
c
c    sro does not rearrange the order of the rows, but does move
c    nonzeroes from one row to another to ensure that if m(i,j) will be
c    in the upper triangle of m with respect to the new ordering, then
c    m(i,j) is stored in row i (and thus m(j,i) is not stored),  whereas
c    if m(i,j) will be in the strict lower triangle of m, then m(j,i) is
c    stored in row j (and thus m(i,j) is not stored).
c
c
c  additional parameters
c
c    q     - integer one-dimensional work array.  dimension = n
c
c    r     - integer one-dimensional work array.  dimension = number of
c            nonzero entries in the upper triangle of m
c
c    dflag - logical variable.  if dflag = .true., then store nonzero
c            diagonal elements at the beginning of the row
c
c-----------------------------------------------------------------------
c
      integer  ip(1),  ia(1), ja(1),  q(1), r(1)
      double precision  a(1),  ak
      logical  dflag
c
c
c--phase 1 -- find row in which to store each nonzero
c----initialize count of nonzeroes to be stored in each row
      do 1 i=1,n
  1     q(i) = 0
c
c----for each nonzero element a(j)
      do 3 i=1,n
        jmin = ia(i)
        jmax = ia(i+1) - 1
        if (jmin.gt.jmax)  go to 3
        do 2 j=jmin,jmax
c
c--------find row (=r(j)) and column (=ja(j)) in which to store a(j) ...
          k = ja(j)
          if (ip(k).lt.ip(i))  ja(j) = i
          if (ip(k).ge.ip(i))  k = i
          r(j) = k
c
c--------... and increment count of nonzeroes (=q(r(j)) in that row
  2       q(k) = q(k) + 1
  3     continue
c
c
c--phase 2 -- find new ia and permutation to apply to (ja,a)
c----determine pointers to delimit rows in permuted (ja,a)
      do 4 i=1,n
        ia(i+1) = ia(i) + q(i)
  4     q(i) = ia(i+1)
c
c----determine where each (ja(j),a(j)) is stored in permuted (ja,a)
c----for each nonzero element (in reverse order)
      ilast = 0
      jmin = ia(1)
      jmax = ia(n+1) - 1
      j = jmax
      do 6 jdummy=jmin,jmax
        i = r(j)
        if (.not.dflag .or. ja(j).ne.i .or. i.eq.ilast)  go to 5
c
c------if dflag, then put diagonal nonzero at beginning of row
          r(j) = ia(i)
          ilast = i
          go to 6
c
c------put (off-diagonal) nonzero in last unused location in row
  5       q(i) = q(i) - 1
          r(j) = q(i)
c
  6     j = j-1
c
c
c--phase 3 -- permute (ja,a) to upper triangular form (wrt new ordering)
      do 8 j=jmin,jmax
  7     if (r(j).eq.j)  go to 8
          k = r(j)
          r(j) = r(k)
          r(k) = k
          jak = ja(k)
          ja(k) = ja(j)
          ja(j) = jak
          ak = a(k)
          a(k) = a(j)
          a(j) = ak
          go to 7
  8     continue
c
      return
      end
