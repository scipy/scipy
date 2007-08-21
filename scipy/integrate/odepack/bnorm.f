      double precision function bnorm (n, a, nra, ml, mu, w)
clll. optimize
c-----------------------------------------------------------------------
c this function computes the norm of a banded n by n matrix,
c stored in the array a, that is consistent with the weighted max-norm
c on vectors, with weights stored in the array w.
c ml and mu are the lower and upper half-bandwidths of the matrix.
c nra is the first dimension of the a array, nra .ge. ml+mu+1.
c in terms of the matrix elements a(i,j), the norm is given by..
c   bnorm = max(i=1,...,n) ( w(i) * sum(j=1,...,n) abs(a(i,j))/w(j) )
c-----------------------------------------------------------------------
      integer n, nra, ml, mu
      integer i, i1, jlo, jhi, j
      double precision a, w
      double precision an, sum
      dimension a(nra,n), w(n)
      an = 0.0d0
      do 20 i = 1,n
        sum = 0.0d0
        i1 = i + mu + 1
        jlo = max0(i-ml,1)
        jhi = min0(i+mu,n)
        do 10 j = jlo,jhi
 10       sum = sum + dabs(a(i1-j,j))/w(j)
        an = dmax1(an,sum*w(i))
 20     continue
      bnorm = an
      return
c----------------------- end of function bnorm -------------------------
      end
