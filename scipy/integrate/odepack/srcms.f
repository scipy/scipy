      subroutine srcms (rsav, isav, job)
c-----------------------------------------------------------------------
c this routine saves or restores (depending on job) the contents of
c the common blocks ls0001, lss001, and eh0001, which are used
c internally by one or more odepack solvers.
c
c rsav = real array of length 224 or more.
c isav = integer array of length 75 or more.
c job  = flag indicating to save or restore the common blocks..
c        job  = 1 if common is to be saved (written to rsav/isav)
c        job  = 2 if common is to be restored (read from rsav/isav)
c        a call with job = 2 presumes a prior call with job = 1.
c-----------------------------------------------------------------------
      integer isav, job
      integer ieh, ils, ilss
      integer i, lenil, leniss, lenrl, lenrss
      double precision rsav,   rls, rlss
      dimension rsav(1), isav(1)
      common /ls0001/ rls(218), ils(39)
      common /lss001/ rlss(6), ilss(34)
      common /eh0001/ ieh(2)
      data lenrl/218/, lenil/39/, lenrss/6/, leniss/34/
c
      if (job .eq. 2) go to 100
      do 10 i = 1,lenrl
 10     rsav(i) = rls(i)
      do 15 i = 1,lenrss
 15     rsav(lenrl+i) = rlss(i)
c
      do 20 i = 1,lenil
 20     isav(i) = ils(i)
      do 25 i = 1,leniss
 25     isav(lenil+i) = ilss(i)
c
      isav(lenil+leniss+1) = ieh(1)
      isav(lenil+leniss+2) = ieh(2)
      return
c
 100  continue
      do 110 i = 1,lenrl
 110     rls(i) = rsav(i)
      do 115 i = 1,lenrss
 115     rlss(i) = rsav(lenrl+i)
c
      do 120 i = 1,lenil
 120     ils(i) = isav(i)
      do 125 i = 1,leniss
 125     ilss(i) = isav(lenil+i)
c
      ieh(1) = isav(lenil+leniss+1)
      ieh(2) = isav(lenil+leniss+2)
      return
c----------------------- end of subroutine srcms -----------------------
      end
