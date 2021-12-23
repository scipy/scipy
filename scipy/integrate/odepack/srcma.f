      subroutine srcma (rsav, isav, job)
c-----------------------------------------------------------------------
c this routine saves or restores (depending on job) the contents of
c the common blocks ls0001, lsa001, and eh0001, which are used
c internally by one or more odepack solvers.
c
c rsav = real array of length 240 or more.
c isav = integer array of length 50 or more.
c job  = flag indicating to save or restore the common blocks..
c        job  = 1 if common is to be saved (written to rsav/isav)
c        job  = 2 if common is to be restored (read from rsav/isav)
c        a call with job = 2 presumes a prior call with job = 1.
c-----------------------------------------------------------------------
      integer isav(*), job
      integer ieh, ils, ilsa
      integer i, lenrls, lenils, lenrla, lenila
      double precision rsav(*)
      double precision rls, rlsa
      common /ls0001/ rls(218), ils(39)
      common /lsa001/ rlsa(22), ilsa(9)
      common /eh0001/ ieh(2)
      data lenrls/218/, lenils/39/, lenrla/22/, lenila/9/
c
      if (job .eq. 2) go to 100
      do 10 i = 1,lenrls
 10     rsav(i) = rls(i)
      do 15 i = 1,lenrla
 15     rsav(lenrls+i) = rlsa(i)
c
      do 20 i = 1,lenils
 20     isav(i) = ils(i)
      do 25 i = 1,lenila
 25     isav(lenils+i) = ilsa(i)
c
      isav(lenils+lenila+1) = ieh(1)
      isav(lenils+lenila+2) = ieh(2)
      return
c
 100  continue
      do 110 i = 1,lenrls
 110     rls(i) = rsav(i)
      do 115 i = 1,lenrla
 115     rlsa(i) = rsav(lenrls+i)
c
      do 120 i = 1,lenils
 120     ils(i) = isav(i)
      do 125 i = 1,lenila
 125     ilsa(i) = isav(lenils+i)
c
      ieh(1) = isav(lenils+lenila+1)
      ieh(2) = isav(lenils+lenila+2)
      return
c----------------------- end of subroutine srcma -----------------------
      end
