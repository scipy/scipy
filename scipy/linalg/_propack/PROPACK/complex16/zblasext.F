c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c

#ifdef _OPENMP
      double precision function pdznrm2(n, x, incx)
      implicit none
      integer n, incx
      complex*16 x(*)
      
      integer i
      double precision sum

      if ((n.gt.0).and.(incx.ne.0)) then    
         sum = 0d0
         if (incx.eq.1) then
c$OMP PARALLEL DO  reduction(+:sum) schedule(static)
            do i=1,n
               sum = sum + dconjg(x(i))*x(i)
            enddo
         else
c$OMP PARALLEL DO firstprivate(incx) reduction(+:sum) schedule(static)
            do i=1,n
               sum = sum + dconjg(x(1+(i-1)*incx))*x(1+(i-1)*incx)
            enddo
         endif
         pdznrm2 = sqrt(sum)
      else
         pdznrm2 = 0d0
      endif   
      return
      end
c
c****************************************************************************
c 
      subroutine pzscal(n, alpha, x , incx)
      implicit none
      integer n, incx
      complex*16 alpha,x(*)
      
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


      subroutine pzdscal(n, alpha, x , incx)
      implicit none
      integer n, incx
      double precision alpha
      complex*16 x(*)
      
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
      subroutine pzcopy(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex*16 x(*),y(*)
      
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
      subroutine pzaxpy(n, alpha, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex*16 alpha,x(*),y(*)
      
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


      subroutine pzdaxpy(n, alpha, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      double precision alpha
      complex*16 x(*),y(*)
      
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
      complex*16 function pzdotc(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex*16 x(*),y(*)
      
      integer i
      complex*16 sum

      if ((n.gt.0).and.(incx.ne.0).and.(incy.ne.0)) then    
         if (incx.eq.1 .and. incy.eq.1) then
            sum = (0d0,0d0)
c$OMP PARALLEL DO reduction(+:sum) schedule(static)
            do i=1,n
               sum = sum + dconjg(x(i)) * y(i)
            enddo
         else
            sum = (0d0,0d0)
c$OMP PARALLEL DO firstprivate(incx, incy) reduction(+:sum) 
c$OMP& schedule(static) 
            do i=1,n
               sum = sum + dconjg(x(1+(i-1)*incx)) * y(1+(i-1)*incy)
            enddo
         endif
         pzdotc = sum
      else
         pzdotc = (0d0,0d0)
      endif   
      return
      end

      complex*16 function pzdotu(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex*16 x(*),y(*)
      
      integer i
      complex*16 sum

      if ((n.gt.0).and.(incx.ne.0).and.(incy.ne.0)) then    
         if (incx.eq.1 .and. incy.eq.1) then
            sum = (0d0,0d0)
c$OMP PARALLEL DO reduction(+:sum) schedule(static)
            do i=1,n
               sum = sum + x(i) * y(i)
            enddo
         else
            sum = (0d0,0d0)
c$OMP PARALLEL DO firstprivate(incx, incy) reduction(+:sum) 
c$OMP& schedule(static) 
            do i=1,n
               sum = sum + x(1+(i-1)*incx) * y(1+(i-1)*incy)
            enddo
         endif
         pzdotu = sum
      else
         pzdotu = (0d0,0d0)
      endif   
      return
      end

c
c****************************************************************************
c 

#else

      double precision function pdznrm2(n, x, incx)
      implicit none
      integer n, incx
      complex*16 x(*)
      double precision dznrm2
      external dznrm2

      pdznrm2 = dznrm2(n, x, incx)
      end
c
c****************************************************************************
c 


      subroutine pzscal(n, alpha, x , incx)
      implicit none
      integer n, incx
      complex*16 alpha,x(*)

      call zscal(n, alpha, x , incx)
      end

c
c****************************************************************************
c 

      subroutine pzdscal(n, alpha, x , incx)
      implicit none
      integer n, incx
      double precision alpha
      complex*16 x(*)

      call zdscal(n, alpha, x , incx)
      end

c
c****************************************************************************
c 

      subroutine pzcopy(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex*16 x(*),y(*)

      call zcopy(n, x , incx, y, incy)
      end
c
c****************************************************************************
c 

      subroutine pzaxpy(n, alpha, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex*16 alpha,x(*),y(*)
     
      call zaxpy(n, alpha, x , incx, y, incy)
      end

      subroutine pzdaxpy(n, alpha, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      double precision alpha
      complex*16 x(*),y(*)
      
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
      
      complex*16 function pzdotc(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex*16 x(*),y(*)
#ifndef __APPLE__
      complex*16 zdotc
      external zdotc

      pzdotc = zdotc(n, x , incx, y, incy)
#else
      integer i
      pzdotc = (0.0,0.0)
      if (incx.eq.1 .and. incy.eq.1) then
         do i=1,n
            pzdotc = pzdotc + conjg(x(i))*y(i)
         enddo
      else
         do i=1,n
            pzdotc = pzdotc + conjg(x(1+(i-1)*incx))*y(1+(i-1)*incy)
         enddo
      endif
#endif      
      end

      complex*16 function pzdotu(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex*16 x(*),y(*)
#ifndef __APPLE__
      complex*16 zdotu
      external zdotu
      
      pzdotu = zdotu(n, x , incx, y, incy)
#else
      integer i
      pzdotu = (0.0,0.0)
      if (incx.eq.1 .and. incy.eq.1) then
         do i=1,n
            pzdotu = pzdotu + x(i)*y(i)
         enddo
      else
         do i=1,n
            pzdotu = pzdotu + x(1+(i-1)*incx)*y(1+(i-1)*incy)
         enddo
      endif
#endif
      end
#endif

c
c****************************************************************************
c 
         
      subroutine pzzero(n, x , incx)
      implicit none
      integer n, incx
      complex*16 x(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
c$OMP PARALLEL DO 
            do i=1,n
               x(i) = (0d0,0d0)
            enddo
         else
c$OMP PARALLEL DO firstprivate(incx) schedule(static)
            do i=1,n
               x(1+(i-1)*incx) = (0d0,0d0)
            enddo
         endif
      endif
      return
      end

c
c****************************************************************************
c 

      subroutine pzset(n, alpha, x , incx)
      implicit none
      integer n, incx
      complex*16 alpha,x(*)
      
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

  
      subroutine pzdaxpby(n,alpha,x,incx,beta,y,incy)
c
c     Y = alpha*X + beta*Y
c     

      implicit none
      double precision one,zero
      parameter(one = 1d0,zero = 0d0)
      integer n,incx,incy,i
      double precision alpha,beta
      complex*16 x(n),y(n)

      if (n.le.0 .or. incy.eq.0 .or. incx.eq.0) return
      if (alpha.eq.zero .and. beta.eq.zero) then
         if (incy.eq.1) then
c$OMP PARALLEL DO  schedule(static)
            do i=1,n
               y(i) = dcmplx(zero,zero)
            enddo
         else
c$OMP PARALLEL DO firstprivate(incy) schedule(static)
            do i=1,n
               y(1+(i-1)*incy) = dcmplx(zero,zero)
            enddo
         endif
         
      else if (alpha.eq.zero .and. beta.ne.zero) then
         
         call pzdscal(n,beta,y,incy)

      else if (alpha.ne.zero .and. beta.eq.zero) then

         if (alpha.eq.one) then
            call pzcopy(n,x,incx,y,incy)
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
            call pzdaxpy(n,alpha,x,incx,y,incy)
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


      subroutine pzaxpby(n,alpha,x,incx,beta,y,incy)
c
c     Y = alpha*X + beta*Y
c     

      implicit none
      double precision one,zero
      parameter(one = 1d0,zero = 0d0)
      integer n,incx,incy,i
      complex*16 alpha,beta,x(n),y(n)

      if (n.le.0 .or. incy.eq.0 .or. incx.eq.0) return
      if (alpha.eq.dcmplx(zero,zero) .and. 
     c     beta.eq.dcmplx(zero,zero)) then
         if (incy.eq.1) then
c$OMP PARALLEL DO  schedule(static)
            do i=1,n
               y(i) = dcmplx(zero,zero)
            enddo
         else
c$OMP PARALLEL DO firstprivate(incy) schedule(static)
            do i=1,n
               y(1+(i-1)*incy) = dcmplx(zero,zero)
            enddo
         endif
         
      else if (alpha.eq.dcmplx(zero,zero) .and. 
     c        beta.ne.dcmplx(zero,zero)) then
         
         call pzscal(n,beta,y,incy)

      else if (alpha.ne.dcmplx(zero,zero) .and. 
     c        beta.eq.dcmplx(zero,zero)) then

         if (alpha.eq.one) then
            call pzcopy(n,x,incx,y,incy)
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
            call pzaxpy(n,alpha,x,incx,y,incy)
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

      subroutine pzaxty(n,alpha,x,incx,y,incy)
c
c     Y = alpha*X*Y
c     

      implicit none
      double precision one,zero
      parameter(one = 1d0,zero = 0d0)
      integer n,incx,incy,i
      complex*16 alpha,x(n),y(n)

      if (n.le.0 .or. incy.eq.0 .or. incx.eq.0) return
      if (alpha.eq.dcmplx(zero,zero)) then
         if (incy.eq.1) then
c$OMP PARALLEL DO schedule(static)
           do i=1,n
               y(i) = dcmplx(zero,zero)
            enddo
         else
c$OMP PARALLEL DO firstprivate(incy) schedule(static)
           do i=1,n
               y(1+(i-1)*incy) = dcmplx(zero,zero)
            enddo
         endif
         
      else if (alpha.ne.dcmplx(zero,zero)) then

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


      subroutine pzdaxty(n,alpha,x,incx,y,incy)
c
c     Y = alpha*X*Y
c     

      implicit none
      double precision one,zero
      parameter(one = 1d0,zero = 0d0)
      integer n,incx,incy,i
      double precision alpha,x(n),y(n)

      if (n.le.0 .or. incy.eq.0 .or. incx.eq.0) return
      if (alpha.eq.zero) then
         if (incy.eq.1) then
c$OMP PARALLEL DO schedule(static)
           do i=1,n
               y(i) = dcmplx(zero,zero)
            enddo
         else
c$OMP PARALLEL DO firstprivate(incy) schedule(static)
           do i=1,n
               y(1+(i-1)*incy) = dcmplx(zero,zero)
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
                 
      subroutine zzero(n, x , incx)
      implicit none
      integer n, incx
      complex*16 x(*)
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
            do i=1,n
               x(i) = dcmplx(0d0,0d0)
            enddo
         else
            do i=1,n
               x(1+(i-1)*incx) = dcmplx(0d0,0d0)
            enddo
         endif
      endif
      return
      end

