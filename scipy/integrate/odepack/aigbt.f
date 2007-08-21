      subroutine aigbt (res, adda, neq, t, y, ydot,
     1                   mb, nb, pw, ipvt, ier )
clll. optimize
      external res, adda
      integer neq, mb, nb, ipvt, ier
      integer i, lenpw, lblox, lpb, lpc
      double precision t, y, ydot, pw
      dimension y(1), ydot(1), pw(1), ipvt(1), neq(1)
c-----------------------------------------------------------------------
c this subroutine computes the initial value
c of the vector ydot satisfying
c     a * ydot = g(t,y)
c when a is nonsingular.  it is called by lsoibt for
c initialization only, when istate = 0 .
c aigbt returns an error flag ier..
c   ier  =  0  means aigbt was successful.
c   ier .ge. 2 means res returned an error flag ires = ier.
c   ier .lt. 0 means the a matrix was found to have a singular
c              diagonal block (hence ydot could not be solved for).
c-----------------------------------------------------------------------
      lblox = mb*mb*nb
      lpb = 1 + lblox
      lpc = lpb + lblox
      lenpw = 3*lblox
      do 10 i = 1,lenpw
 10     pw(i) = 0.0d0
      ier = 1
      call res (neq, t, y, pw, ydot, ier)
      if (ier .gt. 1) return
      call adda (neq, t, y, mb, nb, pw(1), pw(lpb), pw(lpc) )
      call decbt (mb, nb, pw, pw(lpb), pw(lpc), ipvt, ier)
      if (ier .eq. 0) go to 20
      ier = -ier
      return
 20   call solbt (mb, nb, pw, pw(lpb), pw(lpc), ydot, ipvt)
      return
c-------------------- end of subroutine aigbt --------------------------
      end
