      subroutine srcom (rsav, isav, job)
c-----------------------------------------------------------------------
c this routine saves or restores (depending on job) the contents of
c the common blocks ls0001 and eh0001, which are used internally
c by one or more odepack solvers.
c
c rsav = real array of length 218 or more.
c isav = integer array of length 41 or more.
c job  = flag indicating to save or restore the common blocks..
c        job  = 1 if common is to be saved (written to rsav/isav)
c        job  = 2 if common is to be restored (read from rsav/isav)
c        a call with job = 2 presumes a prior call with job = 1.
c-----------------------------------------------------------------------
      integer isav, job
      integer ieh, ils
      integer i, lenils, lenrls
      double precision rsav,   rls
      dimension rsav(1), isav(1)
      common /ls0001/ rls(218), ils(39)
      common /eh0001/ ieh(2)
      data lenrls/218/, lenils/39/
c
      if (job .eq. 2) go to 100
c
      do 10 i = 1,lenrls
 10     rsav(i) = rls(i)
      do 20 i = 1,lenils
 20     isav(i) = ils(i)
      isav(lenils+1) = ieh(1)
      isav(lenils+2) = ieh(2)
      return
c
 100  continue
      do 110 i = 1,lenrls
 110     rls(i) = rsav(i)
      do 120 i = 1,lenils
 120     ils(i) = isav(i)
      ieh(1) = isav(lenils+1)
      ieh(2) = isav(lenils+2)
      return
c----------------------- end of subroutine srcom -----------------------
      end
