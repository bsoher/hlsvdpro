      subroutine clearstat
      implicit none
      include 'stat.h'
      nopx = 0
      nreorth = 0
      ndot = 0
      nitref = 0
      nrstrt = 0
      nbsvd = 0
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
      end

      
      subroutine printstat
      implicit none
      include 'stat.h'

      print *,'+-----------------------------------------------------+'
      print *,'Dimension of Lanczos basis                  = ',nlandim
      print *,'Number of matrix-vector multiplications     = ',nopx
      print *,'Number of reorthogonalizations              = ',nreorth
      print *,'Number of inner products in reorth.         = ',ndot
      print *,'Number of iterative refinement steps        = ',nitref
      print *,'Number of restarts                          = ',nrstrt
      print *,'Number of bidiagonal SVDs calculated        = ',nbsvd
      print *
      print *,'Time spent doing matrix-vector multiply     = ',tmvopx
      print *,'Time spent generating (re)starting vectors  = ',tgetu0
      print *,'Time spent updating mu-recurrence           = ',tupdmu
      print *,'Time spent updating nu-recurrence           = ',tupdnu
      print *,'Time spent in the body of lanbpro           = ',tlanbpro
      print *,'Time spent reorthogonalizing                = ',treorth
      print *,'Time spent reorthogonalizing U_{j+1}        = ',treorthu
      print *,'Time spent reorthogonalizing V_{j}          = ',treorthv
      print *,'Time spent on ext. local reorth. on U_{j+1} = ',telru
      print *,'Time spent on ext. local reorth. on V_{j+1} = ',telrv
      print *,'Time spent computing bidiagonal SVDs        = ',tbsvd
      print *,'Time spent in PDNORM2                       = ',tnorm2
      print *,'Time spent in PDDOT                         = ',tdot
      print *
      print *,'Total time in lansvd                        = ',tlansvd
      print *,'+-----------------------------------------------------+'
      end

      
