c
c\SCCS Information: @(#) 
c FILE: stats.F   SID: 2.1   DATE OF SID: 4/19/96   RELEASE: 2
c     %---------------------------------------------%
c     | Initialize statistic and timing information |
c     | for symmetric Arnoldi code.                 |
c     %---------------------------------------------%
 
      subroutine dstats

c     %--------------------------------%
c     | See stat.doc for documentation |
c     %--------------------------------%
      include   'stat.h'
 
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%

      nopx   = 0
      nbx    = 0
      nrorth = 0
      nitref = 0
      nrstrt = 0
 
      tsaupd = 0.0D+0
      tsaup2 = 0.0D+0
      tsaitr = 0.0D+0
      tseigt = 0.0D+0
      tsgets = 0.0D+0
      tsapps = 0.0D+0
      tsconv = 0.0D+0
      titref = 0.0D+0
      tgetv0 = 0.0D+0
      trvec  = 0.0D+0
 
c     %----------------------------------------------------%
c     | User time including reverse communication overhead |
c     %----------------------------------------------------%
      tmvopx = 0.0D+0
      tmvbx  = 0.0D+0
 
      return
c
c     End of dstats
c
      end
