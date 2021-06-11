c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c

#ifdef _OPENMP
      double precision function pdnrm2(n, x, incx)
      implicit none
      integer n, incx
      double precision x(*)
      
      integer i
      double precision sum

      if ((n.gt.0).and.(incx.ne.0)) then    
         sum = 0d0
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
         pdnrm2 = sqrt(sum)
      else
         pdnrm2 = 0d0
      endif   
      return
      end
c
c****************************************************************************
c 
      

      subroutine pdscal(n, alpha, x , incx)
      implicit none
      integer n, incx
      double precision alpha,x(*)
      
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
         

      subroutine pdcopy(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      double precision x(*),y(*)
      
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
      subroutine pdaxpy(n, alpha, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      double precision alpha,x(*),y(*)
      
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
     
         
      double precision function pddot(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      double precision x(*),y(*)
      
      integer i
      double precision sum

      if ((n.gt.0).and.(incx.ne.0).and.(incy.ne.0)) then    
         if (incx.eq.1 .and. incy.eq.1) then
            sum = 0d0
c$OMP PARALLEL DO reduction(+:sum) schedule(static)
            do i=1,n
               sum = sum + x(i) * y(i)
            enddo
         else
            sum = 0d0
c$OMP PARALLEL DO firstprivate(incx, incy) reduction(+:sum) 
c$OMP& schedule(static) 
            do i=1,n
               sum = sum + x(1+(i-1)*incx) * y(1+(i-1)*incy)
            enddo
         endif
         pddot = sum
      else
         pddot = 0d0
      endif   
      return
      end


#else

      double precision function pdnrm2(n, x, incx)
      implicit none
      integer n, incx
      double precision x(*), dnrm2
      external dnrm2

      pdnrm2 = dnrm2(n, x, incx)
      end

      subroutine pdscal(n, alpha, x , incx)
      implicit none
      integer n, incx
      double precision alpha,x(*)

      call dscal(n, alpha, x , incx)
      end

      subroutine pdcopy(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      double precision x(*),y(*)

      call dcopy(n, x , incx, y, incy)
      end

      subroutine pdaxpy(n, alpha, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      double precision alpha,x(*),y(*)
     
      call daxpy(n, alpha, x , incx, y, incy)
      end

      
      double precision function pddot(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      double precision x(*),y(*),ddot
      external ddot
      
      pddot = ddot(n, x , incx, y, incy)
      end
#endif

c
c****************************************************************************
c 
         
      subroutine pdzero(n, x , incx)
      implicit none
      integer n, incx
      double precision x(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
c$OMP PARALLEL DO 
            do i=1,n
               x(i) = 0d0
            enddo
         else
c$OMP PARALLEL DO firstprivate(incx) schedule(static)
            do i=1,n
               x(1+(i-1)*incx) = 0d0
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


      subroutine pdset(n, alpha, x , incx)
      implicit none
      integer n, incx
      double precision alpha,x(*)
      
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


  
      subroutine pdaxpby(n,alpha,x,incx,beta,y,incy)
c
c     Y = alpha*X + beta*Y
c     

      implicit none
      double precision one,zero
      parameter(one = 1d0,zero = 0d0)
      integer n,incx,incy,i
      double precision alpha,beta,x(n),y(n)

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
         
         call pdscal(n,beta,y,incy)

      else if (alpha.ne.zero .and. beta.eq.zero) then

         if (alpha.eq.one) then
            call pdcopy(n,x,incx,y,incy)
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
            call pdaxpy(n,alpha,x,incx,y,incy)
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



      subroutine pdaxty(n,alpha,x,incx,y,incy)
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

         
