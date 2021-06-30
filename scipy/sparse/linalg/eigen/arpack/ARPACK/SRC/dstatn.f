c
c     %---------------------------------------------%
c     | Initialize statistic and timing information |
c     | for nonsymmetric Arnoldi code.              |
c     %---------------------------------------------%
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\SCCS Information: @(#)
c FILE: statn.F   SID: 2.4   DATE OF SID: 4/20/96   RELEASE: 2
c
      subroutine dstatn
c
c     %--------------------------------%
c     | See stat.doc for documentation |
c     %--------------------------------%
c
      include   'stat.h'
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      nopx   = 0
      nbx    = 0
      nrorth = 0
      nitref = 0
      nrstrt = 0
c
      tnaupd = 0.0D+0
      tnaup2 = 0.0D+0
      tnaitr = 0.0D+0
      tneigh = 0.0D+0
      tngets = 0.0D+0
      tnapps = 0.0D+0
      tnconv = 0.0D+0
      titref = 0.0D+0
      tgetv0 = 0.0D+0
      trvec  = 0.0D+0
c
c     %----------------------------------------------------%
c     | User time including reverse communication overhead |
c     %----------------------------------------------------%
c
      tmvopx = 0.0D+0
      tmvbx  = 0.0D+0
c
      return
c
c
c     %---------------%
c     | End of dstatn |
c     %---------------%
c
      end
