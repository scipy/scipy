c
c     (C) Rasmus Munk Larsen, Stanford University, 2000,2004
c
      subroutine dsafescal(n,alpha,x)
c
c     Scale the vector x by 1/alpha avoiding unnecessary under- and overflow.
c

c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      integer n
      double precision alpha, x(*)

c     %------------%
c     | Parameters |
c     %------------%
      double precision one, zero
      parameter(one = 1.0, zero = 0.0)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,info
      double precision sfmin

c     %----------------------%
c     | External Subroutines |
c     %----------------------%
      external dscal,dlascl

c     %--------------------%
c     | External Functions |
c     %--------------------%
      double precision dlamch
      external dlamch

c     %-----------------%
c     | Data statements |
c     %-----------------%
      save
      data sfmin /-1d0/
      
      if (sfmin.eq.-1d0) then         
         sfmin = dlamch('s')
      endif

      if (abs(alpha).ge.sfmin) then
         call pdscal(n,one/alpha, x, 1)
      else
         call dlascl('General',i,i,alpha,one,n,1,x,n,info)
      endif
      
      end
