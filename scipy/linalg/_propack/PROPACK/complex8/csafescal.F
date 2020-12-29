c
c     (C) Rasmus Munk Larsen, Stanford University, 2000,2004
c
      subroutine csafescal(n,alpha,x)
c
c     Scale the vector x by 1/alpha avoiding unnecessary under- and overflow.
c

c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      integer n
      real alpha
      complex x(*)

c     %------------%
c     | Parameters |
c     %------------%
      real one, zero
      parameter(one = 1.0, zero = 0.0)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,info
      real sfmin

c     %----------------------%
c     | External Subroutines |
c     %----------------------%
      external pcsscal,clascl

c     %--------------------%
c     | External Functions |
c     %--------------------%
      real slamch
      external slamch

c     %-----------------%
c     | Data statements |
c     %-----------------%
      save
      data sfmin /-1d0/
      
      if (sfmin.eq.-1d0) then         
         sfmin = slamch('s')
      endif

      if (abs(alpha).ge.sfmin) then
         call pcsscal(n,one/alpha, x, 1)
      else
         call clascl('General',i,i,alpha,one,n,1,x,n,info)
      endif
      
      end
