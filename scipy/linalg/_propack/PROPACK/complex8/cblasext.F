c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c

#ifdef _OPENMP
      real function pscnrm2(n, x, incx)
      implicit none
      integer n, incx
      complex x(*)
      
      integer i
      real sum

      if ((n.gt.0).and.(incx.ne.0)) then    
         sum = 0.0
         if (incx.eq.1) then
c$OMP PARALLEL DO  reduction(+:sum) schedule(static)
            do i=1,n
               sum = sum + conjg(x(i))*x(i)
            enddo
         else
c$OMP PARALLEL DO firstprivate(incx) reduction(+:sum) schedule(static)
            do i=1,n
               sum = sum + conjg(x(1+(i-1)*incx))*x(1+(i-1)*incx)
            enddo
         endif
         pscnrm2 = sqrt(sum)
      else
         pscnrm2 = 0.0
      endif   
      return
      end
c
c****************************************************************************
c 
      subroutine pcscal(n, alpha, x , incx)
      implicit none
      integer n, incx
      complex alpha,x(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha) schedule(static)
            do i=1,n
               x(i) = alpha*x(i)
            enddo
         else
c$OMP PARALLEL DO firstprivate(alpha, incx) schedule(static)
            do i=1,n
               x(1+(i-1)*incx) = alpha*x(1+(i-1)*incx)
            enddo
         endif
      endif
      return
      end


      subroutine pcsscal(n, alpha, x , incx)
      implicit none
      integer n, incx
      real alpha
      complex x(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha) schedule(static)
            do i=1,n
               x(i) = alpha*x(i)
            enddo
         else
c$OMP PARALLEL DO firstprivate(alpha, incx) schedule(static)
            do i=1,n
               x(1+(i-1)*incx) = alpha*x(1+(i-1)*incx)
            enddo
         endif
      endif
      return
      end
      
c
c****************************************************************************
c 
      subroutine pccopy(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex x(*),y(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0).and.(incy.ne.0)) then         
         if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO  schedule(static)
            do i=1,n
               y(i) = x(i)
            enddo
         else
c$OMP PARALLEL DO firstprivate(incx, incy) schedule(static)
            do i=1,n
               y(1+(i-1)*incy) = x(1+(i-1)*incx)
            enddo
         endif
      endif
      return
      end

c
c****************************************************************************
c 
      subroutine pcaxpy(n, alpha, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex alpha,x(*),y(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0).and.(incy.ne.0)) then         
         if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha)  schedule(static)
           do i=1,n
               y(i) = alpha*x(i) + y(i)
            enddo
         else
c$OMP PARALLEL DO firstprivate(alpha,incx,incy) schedule(static)
            do i=1,n
               y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx) + 
     c              y(1+(i-1)*incy)
            enddo
         endif
      endif
      return
      end


      subroutine pcsaxpy(n, alpha, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      real alpha
      complex x(*),y(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0).and.(incy.ne.0)) then         
         if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha)  schedule(static)
           do i=1,n
               y(i) = alpha*x(i) + y(i)
            enddo
         else
c$OMP PARALLEL DO firstprivate(alpha,incx,incy) schedule(static)
            do i=1,n
               y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx) + 
     c              y(1+(i-1)*incy)
            enddo
         endif
      endif
      return
      end

      
c
c****************************************************************************
c 
      complex function pcdotc(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex x(*),y(*)
      
      integer i
      complex sum

      if ((n.gt.0).and.(incx.ne.0).and.(incy.ne.0)) then    
         if (incx.eq.1 .and. incy.eq.1) then
            sum = (0.0,0.0)
c$OMP PARALLEL DO reduction(+:sum) schedule(static)
            do i=1,n
               sum = sum + conjg(x(i)) * y(i)
            enddo
         else
            sum = (0.0,0.0)
c$OMP PARALLEL DO firstprivate(incx, incy) reduction(+:sum) 
c$OMP& schedule(static) 
            do i=1,n
               sum = sum + conjg(x(1+(i-1)*incx)) * y(1+(i-1)*incy)
            enddo
         endif
         pcdotc = sum
      else
         pcdotc = (0.0,0.0)
      endif   
      return
      end

      complex function pcdotu(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex x(*),y(*)
      
      integer i
      complex sum

      if ((n.gt.0).and.(incx.ne.0).and.(incy.ne.0)) then    
         if (incx.eq.1 .and. incy.eq.1) then
            sum = (0.0,0.0)
c$OMP PARALLEL DO reduction(+:sum) schedule(static)
            do i=1,n
               sum = sum + x(i) * y(i)
            enddo
         else
            sum = (0.0,0.0)
c$OMP PARALLEL DO firstprivate(incx, incy) reduction(+:sum) 
c$OMP& schedule(static) 
            do i=1,n
               sum = sum + x(1+(i-1)*incx) * y(1+(i-1)*incy)
            enddo
         endif
         pcdotu = sum
      else
         pcdotu = (0.0,0.0)
      endif   
      return
      end

c
c****************************************************************************
c 

#else

      real function pscnrm2(n, x, incx)
      implicit none
      integer n, incx
      complex x(*)
      real scnrm2
      external scnrm2

      pscnrm2 = scnrm2(n, x, incx)
      end
c
c****************************************************************************
c 


      subroutine pcscal(n, alpha, x , incx)
      implicit none
      integer n, incx
      complex alpha,x(*)

      call cscal(n, alpha, x , incx)
      end

c
c****************************************************************************
c 

      subroutine pcsscal(n, alpha, x , incx)
      implicit none
      integer n, incx
      real alpha
      complex x(*)

      call csscal(n, alpha, x , incx)
      end

c
c****************************************************************************
c 

      subroutine pccopy(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex x(*),y(*)

      call ccopy(n, x , incx, y, incy)
      end
c
c****************************************************************************
c 

      subroutine pcaxpy(n, alpha, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex alpha,x(*),y(*)
     
      call caxpy(n, alpha, x , incx, y, incy)
      end

      subroutine pcsaxpy(n, alpha, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      real alpha
      complex x(*),y(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0).and.(incy.ne.0)) then         
         if (incx.eq.1 .and. incy.eq.1) then
           do i=1,n
               y(i) = alpha*x(i) + y(i)
            enddo
         else
            do i=1,n
               y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx) + 
     c              y(1+(i-1)*incy)
            enddo
         endif
      endif
      return
      end

c
c****************************************************************************
c 
      
      complex function pcdotc(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex x(*),y(*)
#ifndef __APPLE__
      complex cdotc
      external cdotc
      pcdotc = cdotc(n, x , incx, y, incy)
#else
      integer i
      pcdotc = (0.0,0.0)
      if (incx.eq.1 .and. incy.eq.1) then
         do i=1,n
            pcdotc = pcdotc + conjg(x(i))*y(i)
         enddo
      else
         do i=1,n
            pcdotc = pcdotc + conjg(x(1+(i-1)*incx))*y(1+(i-1)*incy)
         enddo
      endif
#endif
      end

      complex function pcdotu(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex x(*),y(*)
#ifndef __APPLE__
      complex cdotu
      external cdotu
      pcdotu = cdotu(n, x , incx, y, incy)
#else
      integer i
      pcdotu = (0.0,0.0)
      if (incx.eq.1 .and. incy.eq.1) then
         do i=1,n
            pcdotu = pcdotu + x(i)*y(i)
         enddo
      else
         do i=1,n
            pcdotu = pcdotu + x(1+(i-1)*incx)*y(1+(i-1)*incy)
         enddo
      endif
#endif
      end
#endif

c
c****************************************************************************
c 
         
      subroutine pczero(n, x , incx)
      implicit none
      integer n, incx
      complex x(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
c$OMP PARALLEL DO 
            do i=1,n
               x(i) = (0.0,0.0)
            enddo
         else
c$OMP PARALLEL DO firstprivate(incx) schedule(static)
            do i=1,n
               x(1+(i-1)*incx) = (0.0,0.0)
            enddo
         endif
      endif
      return
      end

c
c****************************************************************************
c 

      subroutine pcset(n, alpha, x , incx)
      implicit none
      integer n, incx
      complex alpha,x(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha) schedule(static)
            do i=1,n
               x(i) = alpha
            enddo
         else
c$OMP PARALLEL DO firstprivate(alpha, incx) schedule(static)
            do i=1,n
               x(1+(i-1)*incx) = alpha
            enddo
         endif
      endif
      return
      end

c
c****************************************************************************
c 

  
      subroutine pcsaxpby(n,alpha,x,incx,beta,y,incy)
c
c     Y = alpha*X + beta*Y
c     

      implicit none
      real one,zero
      parameter(one = 1.0,zero = 0.0)
      integer n,incx,incy,i
      real alpha,beta
      complex x(n),y(n)

      if (n.le.0 .or. incy.eq.0 .or. incx.eq.0) return
      if (alpha.eq.zero .and. beta.eq.zero) then
         if (incy.eq.1) then
c$OMP PARALLEL DO  schedule(static)
            do i=1,n
               y(i) = cmplx(zero,zero)
            enddo
         else
c$OMP PARALLEL DO firstprivate(incy) schedule(static)
            do i=1,n
               y(1+(i-1)*incy) = cmplx(zero,zero)
            enddo
         endif
         
      else if (alpha.eq.zero .and. beta.ne.zero) then
         
         call pcsscal(n,beta,y,incy)

      else if (alpha.ne.zero .and. beta.eq.zero) then

         if (alpha.eq.one) then
            call pccopy(n,x,incx,y,incy)
         else
            if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha) schedule(static)
               do i=1,n
                  y(i) = alpha*x(i)
               enddo
            else
c$OMP PARALLEL DO firstprivate(incx, incy, alpha) 
c$OMP& schedule(static) 
               do i=1,n
                  y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx)
               enddo
            endif
         endif

      else

         if (beta.eq.one) then
c DAXPY
            call pcsaxpy(n,alpha,x,incx,y,incy)
         else
            if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha,beta) 
c$OMP& schedule(static) 
               do i=1,n
                  y(i) = alpha*x(i) + beta*y(i)
               enddo
            else
c$OMP PARALLEL DO firstprivate(alpha,beta,incx,incy)  
c$OMP& schedule(static) 
               do i=1,n
                  y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx) + 
     c                 beta*y(1+(i-1)*incy)
               enddo
            endif
         endif
      endif
      return
      end


      subroutine pcaxpby(n,alpha,x,incx,beta,y,incy)
c
c     Y = alpha*X + beta*Y
c     

      implicit none
      real one,zero
      parameter(one = 1.0,zero = 0.0)
      integer n,incx,incy,i
      complex alpha,beta,x(n),y(n)

      if (n.le.0 .or. incy.eq.0 .or. incx.eq.0) return
      if (alpha.eq.cmplx(zero,zero) .and. 
     c     beta.eq.cmplx(zero,zero)) then
         if (incy.eq.1) then
c$OMP PARALLEL DO  schedule(static)
            do i=1,n
               y(i) = cmplx(zero,zero)
            enddo
         else
c$OMP PARALLEL DO firstprivate(incy) schedule(static)
            do i=1,n
               y(1+(i-1)*incy) = cmplx(zero,zero)
            enddo
         endif
         
      else if (alpha.eq.cmplx(zero,zero) .and. 
     c        beta.ne.cmplx(zero,zero)) then
         
         call pcscal(n,beta,y,incy)

      else if (alpha.ne.cmplx(zero,zero) .and. 
     c        beta.eq.cmplx(zero,zero)) then

         if (alpha.eq.one) then
            call pccopy(n,x,incx,y,incy)
         else
            if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha) schedule(static)
               do i=1,n
                  y(i) = alpha*x(i)
               enddo
            else
c$OMP PARALLEL DO firstprivate(incx, incy, alpha) 
c$OMP& schedule(static) 
               do i=1,n
                  y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx)
               enddo
            endif
         endif

      else

         if (beta.eq.one) then
c DAXPY
            call pcaxpy(n,alpha,x,incx,y,incy)
         else
            if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha,beta) 
c$OMP& schedule(static) 
               do i=1,n
                  y(i) = alpha*x(i) + beta*y(i)
               enddo
            else
c$OMP PARALLEL DO firstprivate(alpha,beta,incx,incy)  
c$OMP& schedule(static) 
               do i=1,n
                  y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx) + 
     c                 beta*y(1+(i-1)*incy)
               enddo
            endif
         endif
      endif
      return
      end


c
c****************************************************************************
c 

      subroutine pcaxty(n,alpha,x,incx,y,incy)
c
c     Y = alpha*X*Y
c     

      implicit none
      real one,zero
      parameter(one = 1.0,zero = 0.0)
      integer n,incx,incy,i
      complex alpha,x(n),y(n)

      if (n.le.0 .or. incy.eq.0 .or. incx.eq.0) return
      if (alpha.eq.cmplx(zero,zero)) then
         if (incy.eq.1) then
c$OMP PARALLEL DO schedule(static)
           do i=1,n
               y(i) = cmplx(zero,zero)
            enddo
         else
c$OMP PARALLEL DO firstprivate(incy) schedule(static)
           do i=1,n
               y(1+(i-1)*incy) = cmplx(zero,zero)
            enddo
         endif
         
      else if (alpha.ne.cmplx(zero,zero)) then

         if (alpha.eq.one) then
            if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO  schedule(static)
               do i=1,n
                  y(i) = x(i)*y(i)
               enddo
            else
c$OMP PARALLEL DO shared(y,x) schedule(static)
               do i=1,n
                  y(1+(i-1)*incy) = x(1+(i-1)*incx)*y(1+(i-1)*incy)
               enddo
            endif

         else
            if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha) schedule(static)
               do i=1,n
                  y(i) = alpha*x(i)*y(i)
               enddo
            else
c$OMP PARALLEL DO firstprivate(alpha,incx,incy) 
c$OMP& schedule(static) 
               do i=1,n
                  y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx)*
     c                 y(1+(i-1)*incy)
               enddo
            endif
         endif
      endif
      return
      end


      subroutine pcsaxty(n,alpha,x,incx,y,incy)
c
c     Y = alpha*X*Y
c     

      implicit none
      real one,zero
      parameter(one = 1.0,zero = 0.0)
      integer n,incx,incy,i
      real alpha,x(n),y(n)

      if (n.le.0 .or. incy.eq.0 .or. incx.eq.0) return
      if (alpha.eq.zero) then
         if (incy.eq.1) then
c$OMP PARALLEL DO schedule(static)
           do i=1,n
               y(i) = cmplx(zero,zero)
            enddo
         else
c$OMP PARALLEL DO firstprivate(incy) schedule(static)
           do i=1,n
               y(1+(i-1)*incy) = cmplx(zero,zero)
            enddo
         endif
         
      else if (alpha.ne.zero) then

         if (alpha.eq.one) then
            if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO  schedule(static)
               do i=1,n
                  y(i) = x(i)*y(i)
               enddo
            else
c$OMP PARALLEL DO firstprivate(incx,incy) schedule(static)
               do i=1,n
                  y(1+(i-1)*incy) = x(1+(i-1)*incx)*y(1+(i-1)*incy)
               enddo
            endif

         else
            if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha) schedule(static)
               do i=1,n
                  y(i) = alpha*x(i)*y(i)
               enddo
            else
c$OMP PARALLEL DO firstprivate(alpha,incx,incy) 
c$OMP& schedule(static) 
               do i=1,n
                  y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx)*
     c                 y(1+(i-1)*incy)
               enddo
            endif
         endif
      endif
      return
      end


c
c****************************************************************************
c 
                 
      subroutine czero(n, x , incx)
      implicit none
      integer n, incx
      complex x(*)
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
            do i=1,n
               x(i) = cmplx(0.0,0.0)
            enddo
         else
            do i=1,n
               x(1+(i-1)*incx) = cmplx(0.0,0.0)
            enddo
         endif
      endif
      return
      end

#ifdef BLAH
      complex function cdotc(n, cx, incx, cy, incy)
      complex cx(*), cy(*), c
      integer n, incx, incy

      call cblas_cdotc_sub(%val(n), cx, %val(incx), cy, %val(incy), c)

      cdotc = c
      return
      end

      complex function cdotu(n, cx, incx, cy, incy)
      complex cx(*), cy(*), c
      integer n, incx, incy

      call cblas_cdotu_sub(%val(n), cx, %val(incx), cy, %val(incy), c)

      cdotu = c
      return
      end
#endif
