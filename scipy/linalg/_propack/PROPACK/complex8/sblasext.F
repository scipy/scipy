c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c

#ifdef _OPENMP
      real function psnrm2(n, x, incx)
      implicit none
      integer n, incx
      real x(*)
      
      integer i
      real sum

      if ((n.gt.0).and.(incx.ne.0)) then    
         sum = 0.0
         if (incx.eq.1) then
c$OMP PARALLEL DO  reduction(+:sum) schedule(static)
            do i=1,n
               sum = sum + x(i)**2
            enddo
         else
c$OMP PARALLEL DO firstprivate(incx) reduction(+:sum) schedule(static)
            do i=1,n
               sum = sum + x(1+(i-1)*incx)**2
            enddo
         endif
         psnrm2 = sqrt(sum)
      else
         psnrm2 = 0.0
      endif   
      return
      end
c
c****************************************************************************
c 
      

      subroutine psscal(n, alpha, x , incx)
      implicit none
      integer n, incx
      real alpha,x(*)
      
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
         

      subroutine pscopy(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      real x(*),y(*)
      
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
      subroutine psaxpy(n, alpha, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      real alpha,x(*),y(*)
      
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
     
         
      real function psdot(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      real x(*),y(*)
      
      integer i
      real sum

      if ((n.gt.0).and.(incx.ne.0).and.(incy.ne.0)) then    
         if (incx.eq.1 .and. incy.eq.1) then
            sum = 0.0
c$OMP PARALLEL DO reduction(+:sum) schedule(static)
            do i=1,n
               sum = sum + x(i) * y(i)
            enddo
         else
            sum = 0.0
c$OMP PARALLEL DO firstprivate(incx, incy) reduction(+:sum) 
c$OMP& schedule(static) 
            do i=1,n
               sum = sum + x(1+(i-1)*incx) * y(1+(i-1)*incy)
            enddo
         endif
         psdot = sum
      else
         psdot = 0.0
      endif   
      return
      end


#else

      real function psnrm2(n, x, incx)
      implicit none
      integer n, incx
      real x(*), snrm2
      external snrm2

      psnrm2 = snrm2(n, x, incx)
      end

      subroutine psscal(n, alpha, x , incx)
      implicit none
      integer n, incx
      real alpha,x(*)

      call sscal(n, alpha, x , incx)
      end

      subroutine pscopy(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      real x(*),y(*)

      call scopy(n, x , incx, y, incy)
      end

      subroutine psaxpy(n, alpha, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      real alpha,x(*),y(*)
     
      call saxpy(n, alpha, x , incx, y, incy)
      end

      
      real function psdot(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      real x(*),y(*),sdot
      external sdot
      
      psdot = sdot(n, x , incx, y, incy)
      end
#endif

c
c****************************************************************************
c 
         
      subroutine pszero(n, x , incx)
      implicit none
      integer n, incx
      real x(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
c$OMP PARALLEL DO 
            do i=1,n
               x(i) = 0.0
            enddo
         else
c$OMP PARALLEL DO firstprivate(incx) schedule(static)
            do i=1,n
               x(1+(i-1)*incx) = 0.0
            enddo
         endif
      endif
      return
      end


      subroutine pizero(n, x , incx)
      implicit none
      integer n, incx
      integer x(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
c$OMP PARALLEL DO  schedule(static)
            do i=1,n
               x(i) = 0
            enddo
         else
c$OMP PARALLEL DO firstprivate(incx) schedule(static)
            do i=1,n
               x(1+(i-1)*incx) = 0
            enddo
         endif
      endif
      return
      end


      subroutine psset(n, alpha, x , incx)
      implicit none
      integer n, incx
      real alpha,x(*)
      
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


  
      subroutine psaxpby(n,alpha,x,incx,beta,y,incy)
c
c     Y = alpha*X + beta*Y
c     

      implicit none
      real one,zero
      parameter(one = 1.0,zero = 0.0)
      integer n,incx,incy,i
      real alpha,beta,x(n),y(n)

      if (n.le.0 .or. incy.eq.0 .or. incx.eq.0) return
      if (alpha.eq.zero .and. beta.eq.zero) then
         if (incy.eq.1) then
c$OMP PARALLEL DO  schedule(static)
            do i=1,n
               y(i) = zero
            enddo
         else
c$OMP PARALLEL DO firstprivate(incy) schedule(static)
            do i=1,n
               y(1+(i-1)*incy) = zero
            enddo
         endif
         
      else if (alpha.eq.zero .and. beta.ne.zero) then
         
         call psscal(n,beta,y,incy)

      else if (alpha.ne.zero .and. beta.eq.zero) then

         if (alpha.eq.one) then
            call pscopy(n,x,incx,y,incy)
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
            call psaxpy(n,alpha,x,incx,y,incy)
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



      subroutine psaxty(n,alpha,x,incx,y,incy)
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
               y(i) = zero
            enddo
         else
c$OMP PARALLEL DO firstprivate(incy) schedule(static)
           do i=1,n
               y(1+(i-1)*incy) = zero
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




      subroutine szero(n, x , incx)
      implicit none
      integer n, incx
      real x(*),zero
      parameter (zero = 0.0)          
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
            do i=1,n
               x(i) = zero
            enddo
         else
            do i=1,n
               x(1+(i-1)*incx) = zero
            enddo
         endif
      endif
      return
      end

                 
      subroutine dzero(n, x , incx)
      implicit none
      integer n, incx
      double precision x(*),zero
      parameter (zero = 0.0)          
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
            do i=1,n
               x(i) = zero
            enddo
         else
            do i=1,n
               x(1+(i-1)*incx) = zero
            enddo
         endif
      endif
      return
      end


      subroutine izero(n, x , incx)
      implicit none
      integer n, incx
      integer x(*),zero
      parameter (zero = 0)          
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
            do i=1,n
               x(i) = zero
            enddo
         else
            do i=1,n
               x(1+(i-1)*incx) = zero
            enddo
         endif
      endif
      return
      end

         
