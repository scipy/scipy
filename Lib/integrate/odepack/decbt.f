      subroutine decbt (m, n, a, b, c, ip, ier)
      integer m, n, ip(m,n), ier
      double precision a(m,m,n), b(m,m,n), c(m,m,n)
c-----------------------------------------------------------------------
c the following line is for optimized compilation on llnl compilers.
clll. optimize
c-----------------------------------------------------------------------
c block-tridiagonal matrix decomposition routine.
c written by a. c. hindmarsh.
c latest revision.. november 10, 1983 (ach)
c reference.. ucid-30150
c             solution of block-tridiagonal systems of linear
c             algebraic equations
c             a.c. hindmarsh
c             february 1977
c the input matrix contains three blocks of elements in each block-row,
c including blocks in the (1,3) and (n,n-2) block positions.
c decbt uses block gauss elimination and subroutines dgefa and dgesl
c for solution of blocks.  partial pivoting is done within
c block-rows only.
c
c note.. this version uses linpack routines dgefa/dgesl instead of
c of dec/sol for solution of blocks, and it uses the bla routine ddot
c for dot product calculations.
c
c input..
c     m = order of each block.
c     n = number of blocks in each direction of the matrix.
c         n must be 4 or more.  the complete matrix has order m*n.
c     a = m by m by n array containing diagonal blocks.
c         a(i,j,k) contains the (i,j) element of the k-th block.
c     b = m by m by n array containing the super-diagonal blocks
c         (in b(*,*,k) for k = 1,...,n-1) and the block in the (n,n-2)
c         block position (in b(*,*,n)).
c     c = m by m by n array containing the subdiagonal blocks
c         (in c(*,*,k) for k = 2,3,...,n) and the block in the
c         (1,3) block position (in c(*,*,1)).
c    ip = integer array of length m*n for working storage.
c output..
c a,b,c = m by m by n arrays containing the block lu decomposition
c         of the input matrix.
c    ip = m by n array of pivot information.  ip(*,k) contains
c         information for the k-th digonal block.
c   ier = 0  if no trouble occurred, or
c       = -1 if the input value of m or n was illegal, or
c       = k  if a singular matrix was found in the k-th diagonal block.
c use solbt to solve the associated linear system.
c
c external routines required.. dgefa and dgesl (from linpack) and
c ddot (from the blas, or basic linear algebra package).
c-----------------------------------------------------------------------
      integer nm1, nm2, km1, i, j, k
      double precision dp, ddot
      if (m .lt. 1 .or. n .lt. 4) go to 210
      nm1 = n - 1
      nm2 = n - 2
c process the first block-row. -----------------------------------------
      call dgefa (a, m, m, ip, ier)
      k = 1
      if (ier .ne. 0) go to 200
      do 10 j = 1,m
        call dgesl (a, m, m, ip, b(1,j,1), 0)
        call dgesl (a, m, m, ip, c(1,j,1), 0)
 10     continue
c adjust b(*,*,2). -----------------------------------------------------
      do 40 j = 1,m
        do 30 i = 1,m
          dp = ddot (m, c(i,1,2), m, c(1,j,1), 1)
          b(i,j,2) = b(i,j,2) - dp
 30       continue
 40     continue
c main loop.  process block-rows 2 to n-1. -----------------------------
      do 100 k = 2,nm1
        km1 = k - 1
        do 70 j = 1,m
          do 60 i = 1,m
            dp = ddot (m, c(i,1,k), m, b(1,j,km1), 1)
            a(i,j,k) = a(i,j,k) - dp
 60         continue
 70       continue
        call dgefa (a(1,1,k), m, m, ip(1,k), ier)
        if (ier .ne. 0) go to 200
        do 80 j = 1,m
 80       call dgesl (a(1,1,k), m, m, ip(1,k), b(1,j,k), 0)
 100    continue
c process last block-row and return. -----------------------------------
      do 130 j = 1,m
        do 120 i = 1,m
          dp = ddot (m, b(i,1,n), m, b(1,j,nm2), 1)
          c(i,j,n) = c(i,j,n) - dp
 120      continue
 130    continue
      do 160 j = 1,m
        do 150 i = 1,m
          dp = ddot (m, c(i,1,n), m, b(1,j,nm1), 1)
          a(i,j,n) = a(i,j,n) - dp
 150      continue
 160    continue
      call dgefa (a(1,1,n), m, m, ip(1,n), ier)
      k = n
      if (ier .ne. 0) go to 200
      return
c error returns. -------------------------------------------------------
 200  ier = k
      return
 210  ier = -1
      return
c-----------------------  end of subroutine decbt  ---------------------
      end
