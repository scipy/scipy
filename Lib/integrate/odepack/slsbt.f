      subroutine slsbt (wm, iwm, x, tem)
clll. optimize
      integer iwm
      integer lblox, lpb, lpc, mb, nb
      double precision wm, x, tem
      dimension wm(*), iwm(*), x(1), tem(1)
c-----------------------------------------------------------------------
c this routine acts as an interface between the core integrator
c routine and the solbt routine for the solution of the linear system
c arising from chord iteration.
c communication with slsbt uses the following variables..
c wm    = real work space containing the lu decomposition,
c         starting at wm(3).
c iwm   = integer work space containing pivot information, starting at
c         iwm(21).  iwm also contains block structure parameters
c         mb = iwm(1) and nb = iwm(2).
c x     = the right-hand side vector on input, and the solution vector
c         on output, of length n.
c tem   = vector of work space of length n, not used in this version.
c-----------------------------------------------------------------------
      mb = iwm(1)
      nb = iwm(2)
      lblox = mb*mb*nb
      lpb = 3 + lblox
      lpc = lpb + lblox
      call solbt (mb, nb, wm(3), wm(lpb), wm(lpc), x, iwm(21))
      return
c----------------------- end of subroutine slsbt -----------------------
      end
