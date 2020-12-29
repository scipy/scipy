c
c     (C) Rasmus Munk Larsen, Stanford University, 2000
c
      subroutine clearstat
      implicit none
      include 'stat.h'
      nopx = 0
      nreorth = 0
      ndot = 0
      nitref = 0
      nbsvd = 0
      nrestart = 0
      tmvopx = 0
      tgetu0 = 0
      tupdmu = 0
      tupdnu = 0
      tintv = 0
      tlanbpro = 0
      treorth = 0      
      treorthu = 0
      treorthv = 0
      telru = 0
      telrv = 0
      tbsvd = 0
      tnorm2 = 0
      tdot = 0
      tlansvd = 0
      nlandim = 0
      nsing = 0
      tritzvec = 0
      trestart = 0
      end

      
      subroutine printstat
      implicit none
      include 'stat.h'

      print *,'+---------------------------------------------------',
     c     '--------+'
      print *,'Dimension of Lanczos basis                  = ',nlandim
      print *,'Number of singular values requested         = ',nsing
      print *,'Number of restarts                          = ',nrestart
      print *,'Number of matrix-vector multiplications     = ',nopx
      print *,'Number of reorthogonalizations              = ',nreorth
      print *,'Number of inner products in reorth.         = ',ndot
c      print *,'Number of iterative refinement steps        = ',nitref
      print *,'Number of bidiagonal SVDs calculated        = ',nbsvd
      print *

      print *
      print *,'  Time spent doing matrix-vector multiply   = ',tmvopx
      print *,'  Time spent generating starting vectors    = ',tgetu0
      print *,'    Time spent reorthogonalizing U_{j+1}    = ',treorthu
      print *,'    Time spent reorthogonalizing V_{j}      = ',treorthv
      print *,'  Time spent reorthogonalizing              = ',treorth
      print *,'Total Time spent in LANBPRO                 = ',tlanbpro

c      print *
c      print *,'Time spent updating mu-recurrence           = ',tupdmu
c      print *,'Time spent updating nu-recurrence           = ',tupdnu
c      print *,'Time spent on local reorth. on U_{j+1}      = ',telru
c      print *,'Time spent on local reorth. on V_{j+1}      = ',telrv
c      print *,'Time spent in PDNORM2                       = ',tnorm2
c      print *,'Time spent in PDDOT                         = ',tdot
      print *
      print *,'  Time spent in LANBPRO                     = ',tlanbpro
      print *,'  Time spent computing bidiagonal SVDs      = ',tbsvd
      print *,'  Time spent doing implicit restarts        = ',trestart
      print *,'  Time spent computing Ritz vectors         = ',tritzvec
      print *
      print *,'Total Time spent in LANSVD                  = ',tlansvd
      print *,'+----------------------------------------------------',
     c     '-------+'
      end

      
