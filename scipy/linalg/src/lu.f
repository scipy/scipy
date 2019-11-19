
c     Calculate LU decomposition of a matrix
c     Author: Pearu Peterson, March 2002
c
c     prefixes: d,z,s,c   (double,complex double,float,complex float)
c     suffixes: _c,_r     (column major order,row major order)

      subroutine dlu_c(p,l,u,a,m,n,k,piv,info,permute_l,m1)
      integer m,n,piv(k),i,j,k,permute_l,m1,info
      double precision l(m,k),u(k,n),a(m,n)
      double precision p(m1,m1)

cf2py intent(in,copy) :: a
cf2py intent(out) :: info
cf2py integer intent(hide,cache),depend(k),dimension(k) :: piv
cf2py integer intent(hide),depend(a) :: m = shape(a,0)
cf2py integer intent(hide),depend(a) :: n = shape(a,1)
cf2py integer intent(hide),depend(m,n) :: k = (m<n?m:n)
cf2py intent(out) :: p,l,u
cf2py integer optional,intent(in):: permute_l = 0
cf2py integer intent(hide),depend(permute_l,m) :: m1 = (permute_l?1:m)
cf2py depend(m1) :: p

cf2py callprotoargument double*,double*,double*,double*,int*,int*,int*,int*,int*,int*,int*

      external dgetrf,dlaswp
      call dgetrf(m,n,a,m,piv,info)
      if (info.lt.0) then
         return
      endif
      do 20 j=1,k
         l(j,j) = 1d0
         do 10, i=j+1,m
            l(i,j) = a(i,j)
 10      continue
 20   continue
      do 40 j=1,k
         do 30, i=1,j
            u(i,j) = a(i,j)
 30      continue
 40   continue
      do 60 j=k+1,n
         do 50, i=1,k
            u(i,j) = a(i,j)
 50      continue
 60   continue
      if (permute_l.ne.0) then
         call dlaswp(k,l,m,1,k,piv,-1)
      else
         do 65 i=1,m
            p(i,i)=1d0
 65       continue
         call dlaswp(m,p,m,1,k,piv,-1)
      endif
      end

      subroutine zlu_c(p,l,u,a,m,n,k,piv,info,permute_l,m1)
      integer m,n,piv(k),i,j,k,permute_l,m1,info
      complex*16 l(m,k),u(k,n),a(m,n)
      double precision p(m1,m1)

cf2py intent(in,copy) :: a
cf2py intent(out) :: info
cf2py integer intent(hide,cache),depend(k),dimension(k) :: piv
cf2py integer intent(hide),depend(a) :: m = shape(a,0)
cf2py integer intent(hide),depend(a) :: n = shape(a,1)
cf2py integer intent(hide),depend(m,n) :: k = (m<n?m:n)
cf2py intent(out) :: p,l,u
cf2py integer optional,intent(in):: permute_l = 0
cf2py integer intent(hide),depend(permute_l,m) :: m1 = (permute_l?1:m)
cf2py depend(m1) :: p

cf2py callprotoargument double*,complex_double*,complex_double*,complex_double*,int*,int*,int*,int*,int*,int*,int*

      external zgetrf,zlaswp,dlaswp
      call zgetrf(m,n,a,m,piv,info)
      if (info.lt.0) then
         return
      endif
      do 20 j=1,k
         l(j,j) = 1d0
         do 10, i=j+1,m
            l(i,j) = a(i,j)
 10      continue
 20   continue
      do 40 j=1,k
         do 30, i=1,j
            u(i,j) = a(i,j)
 30      continue
 40   continue
      do 60 j=k+1,n
         do 50, i=1,k
            u(i,j) = a(i,j)
 50      continue
 60   continue
      if (permute_l.ne.0) then
         call zlaswp(k,l,m,1,k,piv,-1)
      else
         do 65 i=1,m
            p(i,i)=1d0
 65       continue
         call dlaswp(m,p,m,1,k,piv,-1)
      endif
      end

      subroutine slu_c(p,l,u,a,m,n,k,piv,info,permute_l,m1)
      integer m,n,piv(k),i,j,k,permute_l,m1,info
      real l(m,k),u(k,n),a(m,n)
      real p(m1,m1)

cf2py intent(in,copy) :: a
cf2py intent(out) :: info
cf2py integer intent(hide,cache),depend(k),dimension(k) :: piv
cf2py integer intent(hide),depend(a) :: m = shape(a,0)
cf2py integer intent(hide),depend(a) :: n = shape(a,1)
cf2py integer intent(hide),depend(m,n) :: k = (m<n?m:n)
cf2py intent(out) :: p,l,u
cf2py integer optional,intent(in):: permute_l = 0
cf2py integer intent(hide),depend(permute_l,m) :: m1 = (permute_l?1:m)
cf2py depend(m1) :: p

cf2py callprotoargument float*,float*,float*,float*,int*,int*,int*,int*,int*,int*,int*

      external sgetrf,slaswp
      call sgetrf(m,n,a,m,piv,info)
      if (info.lt.0) then
         return
      endif
      do 20 j=1,k
         l(j,j) = 1d0
         do 10, i=j+1,m
            l(i,j) = a(i,j)
 10      continue
 20   continue
      do 40 j=1,k
         do 30, i=1,j
            u(i,j) = a(i,j)
 30      continue
 40   continue
      do 60 j=k+1,n
         do 50, i=1,k
            u(i,j) = a(i,j)
 50      continue
 60   continue
      if (permute_l.ne.0) then
         call slaswp(k,l,m,1,k,piv,-1)
      else
         do 65 i=1,m
            p(i,i)=1e0
 65       continue
         call slaswp(m,p,m,1,k,piv,-1)
      endif
      end

      subroutine clu_c(p,l,u,a,m,n,k,piv,info,permute_l,m1)
      integer m,n,piv(k),i,j,k,permute_l,m1,info
      complex l(m,k),u(k,n),a(m,n)
      real p(m1,m1)

cf2py intent(in,copy) :: a
cf2py intent(out) :: info
cf2py integer intent(hide,cache),depend(k),dimension(k) :: piv
cf2py integer intent(hide),depend(a) :: m = shape(a,0)
cf2py integer intent(hide),depend(a) :: n = shape(a,1)
cf2py integer intent(hide),depend(m,n) :: k = (m<n?m:n)
cf2py intent(out) :: p,l,u
cf2py integer optional,intent(in):: permute_l = 0
cf2py integer intent(hide),depend(permute_l,m) :: m1 = (permute_l?1:m)
cf2py depend(m1) :: p

cf2py callprotoargument float*,complex_float*,complex_float*,complex_float*,int*,int*,int*,int*,int*,int*,int*

      external cgetrf,claswp,slaswp
      call cgetrf(m,n,a,m,piv,info)
      if (info.lt.0) then
         return
      endif
      do 20 j=1,k
         l(j,j) = 1d0
         do 10, i=j+1,m
            l(i,j) = a(i,j)
 10      continue
 20   continue
      do 40 j=1,k
         do 30, i=1,j
            u(i,j) = a(i,j)
 30      continue
 40   continue
      do 60 j=k+1,n
         do 50, i=1,k
            u(i,j) = a(i,j)
 50      continue
 60   continue
      if (permute_l.ne.0) then
         call claswp(k,l,m,1,k,piv,-1)
      else
         do 65 i=1,m
            p(i,i)=1e0
 65       continue
         call slaswp(m,p,m,1,k,piv,-1)
      endif
      end
