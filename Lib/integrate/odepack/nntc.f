      subroutine nntc
     *     (n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, z, b, tmp)
clll. optimize
c*** subroutine nntc
c*** numeric solution of the transpose of a sparse nonsymmetric system
c      of linear equations given lu-factorization (compressed pointer
c      storage)
c
c
c       input variables..  n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, b
c       output variables.. z
c
c       parameters used internally..
c fia   - tmp   - temporary vector which gets result of solving ut y = b
c       -           size = n.
c
c  internal variables..
c    jmin, jmax - indices of the first and last positions in a row of
c      u or l  to be used.
c
      integer r(1), c(1), il(1), jl(1), ijl(1), iu(1), ju(1), iju(1)
      double precision l(1), d(1), u(1), b(1), z(1), tmp(1), tmpk,sum
c
c  ******  set tmp to reordered b  *************************************
      do 1 k=1,n
   1    tmp(k) = b(c(k))
c  ******  solve  ut y = b  by forward substitution  *******************
      do 3 k=1,n
        jmin = iu(k)
        jmax = iu(k+1) - 1
        tmpk = -tmp(k)
        if (jmin .gt. jmax) go to 3
        mu = iju(k) - jmin
        do 2 j=jmin,jmax
   2      tmp(ju(mu+j)) = tmp(ju(mu+j)) + tmpk * u(j)
   3    continue
c  ******  solve  lt x = y  by back substitution  **********************
      k = n
      do 6 i=1,n
        sum = -tmp(k)
        jmin = il(k)
        jmax = il(k+1) - 1
        if (jmin .gt. jmax) go to 5
        ml = ijl(k) - jmin
        do 4 j=jmin,jmax
   4      sum = sum + l(j) * tmp(jl(ml+j))
   5    tmp(k) = -sum * d(k)
        z(r(k)) = tmp(k)
        k = k - 1
   6    continue
      return
      end
