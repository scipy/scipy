      double precision function vmnorm (n, v, w)
clll. optimize
c-----------------------------------------------------------------------
c this function routine computes the weighted max-norm
c of the vector of length n contained in the array v, with weights
c contained in the array w of length n..
c   vmnorm = max(i=1,...,n) abs(v(i))*w(i)
c-----------------------------------------------------------------------
      integer n,   i
      double precision v, w,   vm
      dimension v(n), w(n)
      vm = 0.0d0
      do 10 i = 1,n
 10     vm = dmax1(vm,dabs(v(i))*w(i))
      vmnorm = vm
      return
c----------------------- end of function vmnorm ------------------------
      end
