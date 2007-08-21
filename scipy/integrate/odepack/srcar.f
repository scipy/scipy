      subroutine srcar (rsav, isav, job)
c-----------------------------------------------------------------------
c this routine saves or restores (depending on job) the contents of
c the common blocks ls0001, lsa001, lar001, and eh0001, which are used
c internally by one or more odepack solvers.
c
c rsav = real array of length 245 or more.
c isav = integer array of length 59 or more.
c job  = flag indicating to save or restore the common blocks..
c        job  = 1 if common is to be saved (written to rsav/isav)
c        job  = 2 if common is to be restored (read from rsav/isav)
c        a call with job = 2 presumes a prior call with job = 1.
c-----------------------------------------------------------------------
      integer isav, job
      integer ieh, ils, ilsa, ilsr
      integer i, ioff, lenrls, lenils, lenrla, lenila, lenrlr, lenilr
      double precision rsav
      double precision rls, rlsa, rlsr
      dimension rsav(1), isav(1)
      common /ls0001/ rls(218), ils(39)
      common /lsa001/ rlsa(22), ilsa(9)
      common /lsr001/ rlsr(5), ilsr(9)
      common /eh0001/ ieh(2)
      data lenrls/218/, lenils/39/, lenrla/22/, lenila/9/
      data lenrlr/5/, lenilr/9/
c
      if (job .eq. 2) go to 100
      do 10 i = 1,lenrls
 10     rsav(i) = rls(i)
      do 15 i = 1,lenrla
 15     rsav(lenrls+i) = rlsa(i)
      ioff = lenrls + lenrla
      do 20 i = 1,lenrlr
 20     rsav(ioff+i) = rlsr(i)
c
      do 30 i = 1,lenils
 30     isav(i) = ils(i)
      do 35 i = 1,lenila
 35     isav(lenils+i) = ilsa(i)
      ioff = lenils + lenila
      do 40 i = 1,lenilr
 40     isav(ioff+i) = ilsr(i)
c
      ioff = ioff + lenilr
      isav(ioff+1) = ieh(1)
      isav(ioff+2) = ieh(2)
      return
c
 100  continue
      do 110 i = 1,lenrls
 110     rls(i) = rsav(i)
      do 115 i = 1,lenrla
 115     rlsa(i) = rsav(lenrls+i)
      ioff = lenrls + lenrla
      do 120 i = 1,lenrlr
 120     rlsr(i) = rsav(ioff+i)
c
      do 130 i = 1,lenils
 130     ils(i) = isav(i)
      do 135 i = 1,lenila
 135     ilsa(i) = isav(lenils+i)
      ioff = lenils + lenila
      do 140 i = 1,lenilr
 140     ilsr(i) = isav(ioff+i)
c
      ioff = ioff + lenilr
      ieh(1) = isav(ioff+1)
      ieh(2) = isav(ioff+2)
      return
c----------------------- end of subroutine srcar -----------------------
      end
