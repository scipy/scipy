      block data
c-----------------------------------------------------------------------
c this data subprogram loads variables into the internal common
c blocks used by the odepack solvers.  the variables are
c defined as follows..
c   illin  = counter for the number of consecutive times the package
c            was called with illegal input.  the run is stopped when
c            illin reaches 5.
c   ntrep  = counter for the number of consecutive times the package
c            was called with istate = 1 and tout = t.  the run is
c            stopped when ntrep reaches 5.
c   mesflg = flag to control printing of error messages.  1 means print,
c            0 means no printing.
c   lunit  = default value of logical unit number for printing of error
c            messages.
c-----------------------------------------------------------------------
      integer illin, iduma, ntrep, idumb, iowns, icomm, mesflg, lunit
      double precision rowns, rcomm
      common /ls0001/ rowns(209), rcomm(9),
     1   illin, iduma(10), ntrep, idumb(2), iowns(6), icomm(19)
      common /eh0001/ mesflg, lunit
      data illin/0/, ntrep/0/
      data mesflg/1/, lunit/6/
c
c----------------------- end of block data -----------------------------
      end
