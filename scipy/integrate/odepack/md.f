      subroutine md
     *     (n, ia,ja, max, v,l, head,last,next, mark, flag)
clll. optimize
c***********************************************************************
c  md -- minimum degree algorithm (based on element model)
c***********************************************************************
c
c  description
c
c    md finds a minimum degree ordering of the rows and columns of a
c    general sparse matrix m stored in (ia,ja,a) format.
c    when the structure of m is nonsymmetric, the ordering is that
c    obtained for the symmetric matrix  m + m-transpose.
c
c
c  additional parameters
c
c    max  - declared dimension of the one-dimensional arrays v and l.
c           max must be at least  n+2k,  where k is the number of
c           nonzeroes in the strict upper triangle of m + m-transpose
c
c    v    - integer one-dimensional work array.  dimension = max
c
c    l    - integer one-dimensional work array.  dimension = max
c
c    head - integer one-dimensional work array.  dimension = n
c
c    last - integer one-dimensional array used to return the permutation
c           of the rows and columns of m corresponding to the minimum
c           degree ordering.  dimension = n
c
c    next - integer one-dimensional array used to return the inverse of
c           the permutation returned in last.  dimension = n
c
c    mark - integer one-dimensional work array (may be the same as v).
c           dimension = n
c
c    flag - integer error flag.  values and their meanings are -
c             0     no errors detected
c             9n+k  insufficient storage in md
c
c
c  definitions of internal parameters
c
c    ---------+---------------------------------------------------------
c    v(s)     - value field of list entry
c    ---------+---------------------------------------------------------
c    l(s)     - link field of list entry  (0 =) end of list)
c    ---------+---------------------------------------------------------
c    l(vi)    - pointer to element list of uneliminated vertex vi
c    ---------+---------------------------------------------------------
c    l(ej)    - pointer to boundary list of active element ej
c    ---------+---------------------------------------------------------
c    head(d)  - vj =) vj head of d-list d
c             -  0 =) no vertex in d-list d
c
c
c             -                  vi uneliminated vertex
c             -          vi in ek           -       vi not in ek
c    ---------+-----------------------------+---------------------------
c    next(vi) - undefined but nonnegative   - vj =) vj next in d-list
c             -                             -  0 =) vi tail of d-list
c    ---------+-----------------------------+---------------------------
c    last(vi) - (not set until mdp)         - -d =) vi head of d-list d
c             --vk =) compute degree        - vj =) vj last in d-list
c             - ej =) vi prototype of ej    -  0 =) vi not in any d-list
c             -  0 =) do not compute degree -
c    ---------+-----------------------------+---------------------------
c    mark(vi) - mark(vk)                    - nonneg. tag .lt. mark(vk)
c
c
c             -                   vi eliminated vertex
c             -      ei active element      -           otherwise
c    ---------+-----------------------------+---------------------------
c    next(vi) - -j =) vi was j-th vertex    - -j =) vi was j-th vertex
c             -       to be eliminated      -       to be eliminated
c    ---------+-----------------------------+---------------------------
c    last(vi) -  m =) size of ei = m        - undefined
c    ---------+-----------------------------+---------------------------
c    mark(vi) - -m =) overlap count of ei   - undefined
c             -       with ek = m           -
c             - otherwise nonnegative tag   -
c             -       .lt. mark(vk)         -
c
c-----------------------------------------------------------------------
c
      integer  ia(1), ja(1),  v(1), l(1),  head(1), last(1), next(1),
     *   mark(1),  flag,  tag, dmin, vk,ek, tail
      equivalence  (vk,ek)
c
c----initialization
      tag = 0
      call  mdi
     *   (n, ia,ja, max,v,l, head,last,next, mark,tag, flag)
      if (flag.ne.0)  return
c
      k = 0
      dmin = 1
c
c----while  k .lt. n  do
   1  if (k.ge.n)  go to 4
c
c------search for vertex of minimum degree
   2    if (head(dmin).gt.0)  go to 3
          dmin = dmin + 1
          go to 2
c
c------remove vertex vk of minimum degree from degree list
   3    vk = head(dmin)
        head(dmin) = next(vk)
        if (head(dmin).gt.0)  last(head(dmin)) = -dmin
c
c------number vertex vk, adjust tag, and tag vk
        k = k+1
        next(vk) = -k
        last(ek) = dmin - 1
        tag = tag + last(ek)
        mark(vk) = tag
c
c------form element ek from uneliminated neighbors of vk
        call  mdm
     *     (vk,tail, v,l, last,next, mark)
c
c------purge inactive elements and do mass elimination
        call  mdp
     *     (k,ek,tail, v,l, head,last,next, mark)
c
c------update degrees of uneliminated vertices in ek
        call  mdu
     *     (ek,dmin, v,l, head,last,next, mark)
c
        go to 1
c
c----generate inverse permutation from permutation
   4  do 5 k=1,n
        next(k) = -next(k)
   5    last(next(k)) = k
c
      return
      end
