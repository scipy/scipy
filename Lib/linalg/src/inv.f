c
c     Calculate inverse of square matrix
c     Author: Pearu Peterson, March 2002
c
c     prefixes: d,z,s,c   (double,complex double,float,complex float)
c     suffixes: _c,_r     (column major order,row major order)

      subroutine dinv_c(a,n,piv,work,lwork,info)
      integer n,piv(n),lwork

      double precision a(n,n),work(lwork)
cf2py callprotoargument double*,int*,int*,double*,int*,int*

cf2py intent(in,copy,out,out=inv_a) :: a
cf2py intent(out) :: info
cf2py integer intent(hide,cache),depend(n),dimension(n) :: piv
cf2py integer intent(hide),depend(a) :: n = shape(a,0)
cf2py check(shape(a,0)==shape(a,1)) :: a
cf2py intent(hide,cache) :: work
cf2py depend(lwork) :: work
cf2py integer intent(hide),depend(n) :: lwork = 30*n
      external dgetrf,dgetri
      call dgetrf(n,n,a,n,piv,info)
      if (info.ne.0) then
         return
      endif
      call dgetri(n,a,n,piv,work,lwork,info)
      end

      subroutine dinv_r(a,n,piv,info)
      integer n,piv(n)
      double precision a(n,n)
cf2py callprotoargument double*,int*,int*,int*

cf2py intent(c,in,copy,out,out=inv_a) :: a
cf2py intent(out) :: info
cf2py integer intent(hide,cache),depend(n),dimension(n) :: piv
cf2py integer intent(hide),depend(a) :: n = shape(a,0)
cf2py check(shape(a,0)==shape(a,1)) :: a
      external flinalg_dinv_r
      call flinalg_dinv_r(a,n,piv,info)
      end


