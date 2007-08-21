      subroutine odrv
     *     (n, ia,ja,a, p,ip, nsp,isp, path, flag)
clll. optimize
c                                                                 5/2/83
c***********************************************************************
c  odrv -- driver for sparse matrix reordering routines
c***********************************************************************
c
c  description
c
c    odrv finds a minimum degree ordering of the rows and columns
c    of a matrix m stored in (ia,ja,a) format (see below).  for the
c    reordered matrix, the work and storage required to perform
c    gaussian elimination is (usually) significantly less.
c
c    note.. odrv and its subordinate routines have been modified to
c    compute orderings for general matrices, not necessarily having any
c    symmetry.  the miminum degree ordering is computed for the
c    structure of the symmetric matrix  m + m-transpose.
c    modifications to the original odrv module have been made in
c    the coding in subroutine mdi, and in the initial comments in
c    subroutines odrv and md.
c
c    if only the nonzero entries in the upper triangle of m are being
c    stored, then odrv symmetrically reorders (ia,ja,a), (optionally)
c    with the diagonal entries placed first in each row.  this is to
c    ensure that if m(i,j) will be in the upper triangle of m with
c    respect to the new ordering, then m(i,j) is stored in row i (and
c    thus m(j,i) is not stored),  whereas if m(i,j) will be in the
c    strict lower triangle of m, then m(j,i) is stored in row j (and
c    thus m(i,j) is not stored).
c
c
c  storage of sparse matrices
c
c    the nonzero entries of the matrix m are stored row-by-row in the
c    array a.  to identify the individual nonzero entries in each row,
c    we need to know in which column each entry lies.  these column
c    indices are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
c    ja(k) = j.  to identify the individual rows, we need to know where
c    each row starts.  these row pointers are stored in the array ia.
c    i.e., if m(i,j) is the first nonzero entry (stored) in the i-th row
c    and  a(k) = m(i,j),  then  ia(i) = k.  moreover, ia(n+1) points to
c    the first location following the last element in the last row.
c    thus, the number of entries in the i-th row is  ia(i+1) - ia(i),
c    the nonzero entries in the i-th row are stored consecutively in
c
c            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
c
c    and the corresponding column indices are stored consecutively in
c
c            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
c
c    since the coefficient matrix is symmetric, only the nonzero entries
c    in the upper triangle need be stored.  for example, the matrix
c
c             ( 1  0  2  3  0 )
c             ( 0  4  0  0  0 )
c         m = ( 2  0  5  6  0 )
c             ( 3  0  6  7  8 )
c             ( 0  0  0  8  9 )
c
c    could be stored as
c
c            - 1  2  3  4  5  6  7  8  9 10 11 12 13
c         ---+--------------------------------------
c         ia - 1  4  5  8 12 14
c         ja - 1  3  4  2  1  3  4  1  3  4  5  4  5
c          a - 1  2  3  4  2  5  6  3  6  7  8  8  9
c
c    or (symmetrically) as
c
c            - 1  2  3  4  5  6  7  8  9
c         ---+--------------------------
c         ia - 1  4  5  7  9 10
c         ja - 1  3  4  2  3  4  4  5  5
c          a - 1  2  3  4  5  6  7  8  9          .
c
c
c  parameters
c
c    n    - order of the matrix
c
c    ia   - integer one-dimensional array containing pointers to delimit
c           rows in ja and a.  dimension = n+1
c
c    ja   - integer one-dimensional array containing the column indices
c           corresponding to the elements of a.  dimension = number of
c           nonzero entries in (the upper triangle of) m
c
c    a    - real one-dimensional array containing the nonzero entries in
c           (the upper triangle of) m, stored by rows.  dimension =
c           number of nonzero entries in (the upper triangle of) m
c
c    p    - integer one-dimensional array used to return the permutation
c           of the rows and columns of m corresponding to the minimum
c           degree ordering.  dimension = n
c
c    ip   - integer one-dimensional array used to return the inverse of
c           the permutation returned in p.  dimension = n
c
c    nsp  - declared dimension of the one-dimensional array isp.  nsp
c           must be at least  3n+4k,  where k is the number of nonzeroes
c           in the strict upper triangle of m
c
c    isp  - integer one-dimensional array used for working storage.
c           dimension = nsp
c
c    path - integer path specification.  values and their meanings are -
c             1  find minimum degree ordering only
c             2  find minimum degree ordering and reorder symmetrically
c                  stored matrix (used when only the nonzero entries in
c                  the upper triangle of m are being stored)
c             3  reorder symmetrically stored matrix as specified by
c                  input permutation (used when an ordering has already
c                  been determined and only the nonzero entries in the
c                  upper triangle of m are being stored)
c             4  same as 2 but put diagonal entries at start of each row
c             5  same as 3 but put diagonal entries at start of each row
c
c    flag - integer error flag.  values and their meanings are -
c               0    no errors detected
c              9n+k  insufficient storage in md
c             10n+1  insufficient storage in odrv
c             11n+1  illegal path specification
c
c
c  conversion from real to double precision
c
c    change the real declarations in odrv and sro to double precision
c    declarations.
c
c-----------------------------------------------------------------------
c
      integer  ia(1), ja(1),  p(1), ip(1),  isp(1),  path,  flag,
     *   v, l, head,  tmp, q
      double precision  a(1)
      logical  dflag
c
c----initialize error flag and validate path specification
      flag = 0
      if (path.lt.1 .or. 5.lt.path)  go to 111
c
c----allocate storage and find minimum degree ordering
      if ((path-1) * (path-2) * (path-4) .ne. 0)  go to 1
        max = (nsp-n)/2
        v    = 1
        l    = v     +  max
        head = l     +  max
        next = head  +  n
        if (max.lt.n)  go to 110
c
        call  md
     *     (n, ia,ja, max,isp(v),isp(l), isp(head),p,ip, isp(v), flag)
        if (flag.ne.0)  go to 100
c
c----allocate storage and symmetrically reorder matrix
   1  if ((path-2) * (path-3) * (path-4) * (path-5) .ne. 0)  go to 2
        tmp = (nsp+1) -      n
        q   = tmp     - (ia(n+1)-1)
        if (q.lt.1)  go to 110
c
        dflag = path.eq.4 .or. path.eq.5
        call sro
     *     (n,  ip,  ia, ja, a,  isp(tmp),  isp(q),  dflag)
c
   2  return
c
c ** error -- error detected in md
 100  return
c ** error -- insufficient storage
 110  flag = 10*n + 1
      return
c ** error -- illegal path specified
 111  flag = 11*n + 1
      return
      end
