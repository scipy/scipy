c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c          BASIC MATRIX-VECTOR OPERATIONS - MATVEC MODULE              c
c         Matrix-vector Mulitiplications and Triang. Solves            c
c----------------------------------------------------------------------c
c contents: (as of Nov 18, 1991)                                       c
c----------                                                            c
c 1) Matrix-vector products:                                           c
c---------------------------                                           c
c amux  : A times a vector. Compressed Sparse Row (CSR) format.        c
c amuxms: A times a vector. Modified Compress Sparse Row format.       c
c atmux : Transp(A) times a vector. CSR format.                        c
c atmuxr: Transp(A) times a vector. CSR format. A rectangular.         c
c amuxe : A times a vector. Ellpack/Itpack (ELL) format.               c
c amuxd : A times a vector. Diagonal (DIA) format.                     c
c amuxj : A times a vector. Jagged Diagonal (JAD) format.              c
c vbrmv : Sparse matrix-full vector product, in VBR format             c
c                                                                      c
c 2) Triangular system solutions:                                      c
c-------------------------------                                       c
c lsol  : Unit Lower Triang. solve. Compressed Sparse Row (CSR) format.c
c ldsol : Lower Triang. solve.  Modified Sparse Row (MSR) format.      c
c lsolc : Unit Lower Triang. solve. Comp. Sparse Column (CSC) format.  c
c ldsolc: Lower Triang. solve. Modified Sparse Column (MSC) format.    c
c ldsoll: Lower Triang. solve with level scheduling. MSR format.       c
c usol  : Unit Upper Triang. solve. Compressed Sparse Row (CSR) format.c
c udsol : Upper Triang. solve.  Modified Sparse Row (MSR) format.      c
c usolc : Unit Upper Triang. solve. Comp. Sparse Column (CSC) format.  c
c udsolc: Upper Triang. solve.  Modified Sparse Column (MSC) format.   c
c----------------------------------------------------------------------c
c 1)     M A T R I X    B Y    V E C T O R     P R O D U C T S         c
c----------------------------------------------------------------------c
      subroutine zamux (n, x, y, a,ja,ia) 
      complex*16  x(*), y(*), a(*) 
      integer n, ja(*), ia(*)
c-----------------------------------------------------------------------
c         A times a vector
c----------------------------------------------------------------------- 
c multiplies a matrix by a vector using the dot product form
c Matrix A is stored in compressed sparse row storage.
c
c on entry:
c----------
c n     = row dimension of A
c x     = real array of length equal to the column dimension of
c         the A matrix.
c a, ja,
c    ia = input matrix in compressed sparse row format.
c
c on return:
c-----------
c y     = real array of length n, containing the product y=Ax
c
c-----------------------------------------------------------------------
c local variables
c
      complex*16 t
      integer i, k
c-----------------------------------------------------------------------
      do 100 i = 1,n
c
c     compute the inner product of row i with vector x
c 
         t = 0.0d0
         do 99 k=ia(i), ia(i+1)-1 
            t = t + a(k)*x(ja(k))
 99      continue
c
c     store result in y(i) 
c
         y(i) = t
 100  continue
c
      return
c---------end-of-amux---------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine zamuxms (n, x, y, a,ja)
      complex*16  x(*), y(*), a(*)
      integer n, ja(*)
c-----------------------------------------------------------------------
c         A times a vector in MSR format
c-----------------------------------------------------------------------
c multiplies a matrix by a vector using the dot product form
c Matrix A is stored in Modified Sparse Row storage.
c
c on entry:
c----------
c n     = row dimension of A
c x     = real array of length equal to the column dimension of
c         the A matrix.
c a, ja,= input matrix in modified compressed sparse row format.
c
c on return:
c-----------
c y     = real array of length n, containing the product y=Ax
c
c-----------------------------------------------------------------------
c local variables
c
      integer i, k
c-----------------------------------------------------------------------
        do 10 i=1, n
        y(i) = a(i)*x(i)
 10     continue
      do 100 i = 1,n
c
c     compute the inner product of row i with vector x
c
         do 99 k=ja(i), ja(i+1)-1
            y(i) = y(i) + a(k) *x(ja(k))
 99      continue
 100  continue
c
      return
c---------end-of-amuxm--------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine zatmux (n, x, y, a, ja, ia)
      complex*16 x(*), y(*), a(*) 
      integer n, ia(*), ja(*)
c-----------------------------------------------------------------------
c         transp( A ) times a vector
c----------------------------------------------------------------------- 
c multiplies the transpose of a matrix by a vector when the original
c matrix is stored in compressed sparse row storage. Can also be
c viewed as the product of a matrix by a vector when the original
c matrix is stored in the compressed sparse column format.
c-----------------------------------------------------------------------
c
c on entry:
c----------
c n     = row dimension of A
c x     = real array of length equal to the column dimension of
c         the A matrix.
c a, ja,
c    ia = input matrix in compressed sparse row format.
c
c on return:
c-----------
c y     = real array of length n, containing the product y=transp(A)*x
c
c-----------------------------------------------------------------------
c     local variables 
c
      integer i, k 
c-----------------------------------------------------------------------
c
c     zero out output vector
c 
      do 1 i=1,n
         y(i) = 0.0
 1    continue
c
c loop over the rows
c
      do 100 i = 1,n
         do 99 k=ia(i), ia(i+1)-1 
            y(ja(k)) = y(ja(k)) + x(i)*a(k)
 99      continue
 100  continue
c
      return
c-------------end-of-atmux---------------------------------------------- 
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine zatmuxr (m, n, x, y, a, ja, ia)
      complex*16 x(*), y(*), a(*) 
      integer m, n, ia(*), ja(*)
c-----------------------------------------------------------------------
c         transp( A ) times a vector, A can be rectangular
c----------------------------------------------------------------------- 
c See also atmux.  The essential difference is how the solution vector
c is initially zeroed.  If using this to multiply rectangular CSC 
c matrices by a vector, m number of rows, n is number of columns.
c-----------------------------------------------------------------------
c
c on entry:
c----------
c m     = column dimension of A
c n     = row dimension of A
c x     = real array of length equal to the column dimension of
c         the A matrix.
c a, ja,
c    ia = input matrix in compressed sparse row format.
c
c on return:
c-----------
c y     = real array of length n, containing the product y=transp(A)*x
c
c-----------------------------------------------------------------------
c     local variables 
c
      integer i, k 
c-----------------------------------------------------------------------
c
c     zero out output vector
c 
      do 1 i=1,m
         y(i) = 0.0
 1    continue
c
c loop over the rows
c
      do 100 i = 1,n
         do 99 k=ia(i), ia(i+1)-1 
            y(ja(k)) = y(ja(k)) + x(i)*a(k)
 99      continue
 100  continue
c
      return
c-------------end-of-atmuxr--------------------------------------------- 
      end
c-----------------------------------------------------------------------
      subroutine zactmux (n, x, y, a, ja, ia)
      complex*16 x(*), y(*), a(*) 
      integer n, ia(*), ja(*)
c-----------------------------------------------------------------------
c         conj(transp( A )) times a vector
c----------------------------------------------------------------------- 
c multiplies the transpose of a matrix by a vector when the original
c matrix is stored in compressed sparse row storage. Can also be
c viewed as the product of a matrix by a vector when the original
c matrix is stored in the compressed sparse column format.
c-----------------------------------------------------------------------
c
c on entry:
c----------
c n     = row dimension of A
c x     = real array of length equal to the column dimension of
c         the A matrix.
c a, ja,
c    ia = input matrix in compressed sparse row format.
c
c on return:
c-----------
c y     = real array of length n, containing the product y=transp(A)*x
c
c-----------------------------------------------------------------------
c     local variables 
c
      integer i, k 
c-----------------------------------------------------------------------
c
c     zero out output vector
c 
      do 1 i=1,n
         y(i) = 0.0
 1    continue
c
c loop over the rows
c
      do 100 i = 1,n
         do 99 k=ia(i), ia(i+1)-1 
            y(ja(k)) = y(ja(k)) + x(i)*conjg(a(k))
 99      continue
 100  continue
c
      return
c-------------end-of-actmux---------------------------------------------- 
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine zactmuxr (m, n, x, y, a, ja, ia)
      complex*16 x(*), y(*), a(*) 
      integer m, n, ia(*), ja(*)
c-----------------------------------------------------------------------
c         conj(transp( A )) times a vector, A can be rectangular
c----------------------------------------------------------------------- 
c See also atmux.  The essential difference is how the solution vector
c is initially zeroed.  If using this to multiply rectangular CSC 
c matrices by a vector, m number of rows, n is number of columns.
c-----------------------------------------------------------------------
c
c on entry:
c----------
c m     = column dimension of A
c n     = row dimension of A
c x     = real array of length equal to the column dimension of
c         the A matrix.
c a, ja,
c    ia = input matrix in compressed sparse row format.
c
c on return:
c-----------
c y     = real array of length n, containing the product y=transp(A)*x
c
c-----------------------------------------------------------------------
c     local variables 
c
      integer i, k 
c-----------------------------------------------------------------------
c
c     zero out output vector
c 
      do 1 i=1,m
         y(i) = 0.0
 1    continue
c
c loop over the rows
c
      do 100 i = 1,n
         do 99 k=ia(i), ia(i+1)-1 
            y(ja(k)) = y(ja(k)) + x(i)*conjg(a(k))
 99      continue
 100  continue
c
      return
c-------------end-of-actmuxr--------------------------------------------- 
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine zamuxe (n,x,y,na,ncol,a,ja) 
      complex*16 x(n), y(n), a(na,*)  
      integer  n, na, ncol, ja(na,*)
c-----------------------------------------------------------------------
c        A times a vector in Ellpack Itpack format (ELL)               
c----------------------------------------------------------------------- 
c multiplies a matrix by a vector when the original matrix is stored 
c in the ellpack-itpack sparse format.
c-----------------------------------------------------------------------
c
c on entry:
c----------
c n     = row dimension of A
c x     = real array of length equal to the column dimension of
c         the A matrix.
c na    = integer. The first dimension of arrays a and ja
c         as declared by the calling program.
c ncol  = integer. The number of active columns in array a.
c         (i.e., the number of generalized diagonals in matrix.)
c a, ja = the real and integer arrays of the itpack format
c         (a(i,k),k=1,ncol contains the elements of row i in matrix
c          ja(i,k),k=1,ncol contains their column numbers) 
c
c on return:
c-----------
c y     = real array of length n, containing the product y=y=A*x
c
c-----------------------------------------------------------------------
c local variables
c
      integer i, j 
c-----------------------------------------------------------------------
      do 1 i=1, n
         y(i) = 0.0 
 1    continue
      do 10 j=1,ncol
         do 25 i = 1,n
            y(i) = y(i)+a(i,j)*x(ja(i,j))
 25      continue
 10   continue
c
      return
c--------end-of-amuxe--------------------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine zamuxd (n,x,y,diag,ndiag,idiag,ioff) 
      integer n, ndiag, idiag, ioff(idiag) 
      complex*16 x(n), y(n), diag(ndiag,idiag)
c-----------------------------------------------------------------------
c        A times a vector in Diagonal storage format (DIA) 
c----------------------------------------------------------------------- 
c multiplies a matrix by a vector when the original matrix is stored 
c in the diagonal storage format.
c-----------------------------------------------------------------------
c
c on entry:
c----------
c n     = row dimension of A
c x     = real array of length equal to the column dimension of
c         the A matrix.
c ndiag  = integer. The first dimension of array adiag as declared in
c         the calling program.
c idiag  = integer. The number of diagonals in the matrix.
c diag   = real array containing the diagonals stored of A.
c idiag  = number of diagonals in matrix.
c diag   = real array of size (ndiag x idiag) containing the diagonals
c          
c ioff   = integer array of length idiag, containing the offsets of the
c   	   diagonals of the matrix:
c          diag(i,k) contains the element a(i,i+ioff(k)) of the matrix.
c
c on return:
c-----------
c y     = real array of length n, containing the product y=A*x
c
c-----------------------------------------------------------------------
c local variables 
c
      integer j, k, io, i1, i2 
c-----------------------------------------------------------------------
      do 1 j=1, n
         y(j) = 0.0d0
 1    continue	
      do 10 j=1, idiag
         io = ioff(j)
         i1 = max0(1,1-io)
         i2 = min0(n,n-io)
         do 9 k=i1, i2	
            y(k) = y(k)+diag(k,j)*x(k+io)
 9       continue
 10   continue
c 
      return
c----------end-of-amuxd-------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine zamuxj (n, x, y, jdiag, a, ja, ia)
      integer n, jdiag, ja(*), ia(*)
      complex*16 x(n), y(n), a(*)  
c-----------------------------------------------------------------------
c        A times a vector in Jagged-Diagonal storage format (JAD) 
c----------------------------------------------------------------------- 
c multiplies a matrix by a vector when the original matrix is stored 
c in the jagged diagonal storage format.
c-----------------------------------------------------------------------
c
c on entry:
c----------
c n      = row dimension of A
c x      = real array of length equal to the column dimension of
c         the A matrix.
c jdiag  = integer. The number of jadded-diagonals in the data-structure.
c a      = real array containing the jadded diagonals of A stored
c          in succession (in decreasing lengths) 
c j      = integer array containing the colum indices of the 
c          corresponding elements in a.
c ia     = integer array containing the lengths of the  jagged diagonals
c
c on return:
c-----------
c y      = real array of length n, containing the product y=A*x
c
c Note:
c------- 
c Permutation related to the JAD format is not performed.
c this can be done by:
c     call permvec (n,y,y,iperm) 
c after the call to amuxj, where iperm is the permutation produced
c by csrjad.
c-----------------------------------------------------------------------
c local variables 
c
      integer i, ii, k1, len, j 
c-----------------------------------------------------------------------
      do 1 i=1, n
         y(i) = 0.0d0
 1    continue
      do 70 ii=1, jdiag
         k1 = ia(ii)-1
         len = ia(ii+1)-k1-1
         do 60 j=1,len
            y(j)= y(j)+a(k1+j)*x(ja(k1+j)) 
 60      continue
 70   continue
c
      return
c----------end-of-amuxj------------------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine zvbrmv(nr, nc, ia, ja, ka, a, kvstr, kvstc, x, b)
c-----------------------------------------------------------------------
      integer nr, nc, ia(nr+1), ja(*), ka(*), kvstr(nr+1), kvstc(*)
      complex*16  a(*), x(*), b(*)
c-----------------------------------------------------------------------
c     Sparse matrix-full vector product, in VBR format.
c-----------------------------------------------------------------------
c     On entry:
c--------------
c     nr, nc  = number of block rows and columns in matrix A
c     ia,ja,ka,a,kvstr,kvstc = matrix A in variable block row format
c     x       = multiplier vector in full format
c
c     On return:
c---------------
c     b = product of matrix A times vector x in full format
c
c     Algorithm:
c---------------
c     Perform multiplication by traversing a in order.
c
c-----------------------------------------------------------------------
c-----local variables
      integer n, i, j, ii, jj, k, istart, istop
      complex*16  xjj
c---------------------------------
      n = kvstc(nc+1)-1
      do i = 1, n
         b(i) = 0.d0
      enddo
c---------------------------------
      k = 1
      do i = 1, nr
         istart = kvstr(i)
         istop  = kvstr(i+1)-1
         do j = ia(i), ia(i+1)-1
            do jj = kvstc(ja(j)), kvstc(ja(j)+1)-1
               xjj = x(jj)
               do ii = istart, istop
                  b(ii) = b(ii) + xjj*a(k)
                  k = k + 1
               enddo
            enddo
         enddo
      enddo
c---------------------------------
      return
      end
c-----------------------------------------------------------------------
c----------------------end-of-vbrmv-------------------------------------
c-----------------------------------------------------------------------
c----------------------------------------------------------------------c
c 2)     T R I A N G U L A R    S Y S T E M    S O L U T I O N S       c
c----------------------------------------------------------------------c
      subroutine zlsol (n,x,y,al,jal,ial)
      integer n, jal(*),ial(n+1) 
      complex*16  x(n), y(n), al(*) 
c-----------------------------------------------------------------------
c   solves    L x = y ; L = lower unit triang. /  CSR format
c----------------------------------------------------------------------- 
c solves a unit lower triangular system by standard (sequential )
c forward elimination - matrix stored in CSR format. 
c-----------------------------------------------------------------------
c
c On entry:
c---------- 
c n      = integer. dimension of problem.
c y      = real array containg the right side.
c
c al,
c jal,
c ial,    = Lower triangular matrix stored in compressed sparse row
c          format. 
c
c On return:
c----------- 
c	x  = The solution of  L x  = y.
c--------------------------------------------------------------------
c local variables 
c
      integer k, j 
      complex*16  t
c-----------------------------------------------------------------------
      x(1) = y(1) 
      do 150 k = 2, n
         t = y(k) 
         do 100 j = ial(k), ial(k+1)-1
            t = t-al(j)*x(jal(j))
 100     continue
         x(k) = t 
 150  continue
c
      return
c----------end-of-lsol-------------------------------------------------- 
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine zldsol (n,x,y,al,jal) 
      integer n, jal(*) 
      complex*16 x(n), y(n), al(*) 
c----------------------------------------------------------------------- 
c     Solves L x = y    L = triangular. MSR format 
c-----------------------------------------------------------------------
c solves a (non-unit) lower triangular system by standard (sequential) 
c forward elimination - matrix stored in MSR format 
c with diagonal elements already inverted (otherwise do inversion,
c al(1:n) = 1.0/al(1:n),  before calling ldsol).
c-----------------------------------------------------------------------
c
c On entry:
c---------- 
c n      = integer. dimension of problem.
c y      = real array containg the right hand side.
c
c al,
c jal,   = Lower triangular matrix stored in Modified Sparse Row 
c          format. 
c
c On return:
c----------- 
c	x = The solution of  L x = y .
c--------------------------------------------------------------------
c local variables 
c
      integer k, j 
      complex*16 t 
c-----------------------------------------------------------------------
      x(1) = y(1)*al(1) 
      do 150 k = 2, n
         t = y(k) 
         do 100 j = jal(k), jal(k+1)-1
            t = t - al(j)*x(jal(j))
 100     continue
         x(k) = al(k)*t 
 150  continue
      return
c----------end-of-ldsol-------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine zlsolc (n,x,y,al,jal,ial)
      integer n, jal(*),ial(*) 
      complex*16  x(n), y(n), al(*) 
c-----------------------------------------------------------------------
c       SOLVES     L x = y ;    where L = unit lower trang. CSC format
c-----------------------------------------------------------------------
c solves a unit lower triangular system by standard (sequential )
c forward elimination - matrix stored in CSC format. 
c-----------------------------------------------------------------------
c
c On entry:
c---------- 
c n      = integer. dimension of problem.
c y      = complex*16 array containg the right side.
c
c al,
c jal,
c ial,    = Lower triangular matrix stored in compressed sparse column 
c          format. 
c
c On return:
c----------- 
c	x  = The solution of  L x  = y.
c-----------------------------------------------------------------------
c local variables 
c
      integer k, j
      complex*16 t
c-----------------------------------------------------------------------
      do 140 k=1,n
         x(k) = y(k) 
 140  continue
      do 150 k = 1, n-1
         t = x(k) 
         do 100 j = ial(k), ial(k+1)-1
            x(jal(j)) = x(jal(j)) - t*al(j) 
 100     continue
 150  continue
c
      return
c----------end-of-lsolc------------------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine zldsolc (n,x,y,al,jal) 
      integer n, jal(*)
      complex*16 x(n), y(n), al(*)
c-----------------------------------------------------------------------
c    Solves     L x = y ;    L = nonunit Low. Triang. MSC format 
c----------------------------------------------------------------------- 
c solves a (non-unit) lower triangular system by standard (sequential) 
c forward elimination - matrix stored in Modified Sparse Column format 
c with diagonal elements already inverted (otherwise do inversion,
c al(1:n) = 1.0/al(1:n),  before calling ldsol).
c-----------------------------------------------------------------------
c
c On entry:
c---------- 
c n      = integer. dimension of problem.
c y      = real array containg the right hand side.
c
c al,
c jal,
c ial,    = Lower triangular matrix stored in Modified Sparse Column
c           format.
c
c On return:
c----------- 
c	x = The solution of  L x = y .
c--------------------------------------------------------------------
c local variables
c
      integer k, j
      complex*16 t 
c-----------------------------------------------------------------------
      do 140 k=1,n
         x(k) = y(k) 
 140  continue
      do 150 k = 1, n 
         x(k) = x(k)*al(k) 
         t = x(k) 
         do 100 j = jal(k), jal(k+1)-1
            x(jal(j)) = x(jal(j)) - t*al(j) 
 100     continue
 150  continue
c
      return
c----------end-of-lsolc------------------------------------------------ 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------  
      subroutine zldsoll (n,x,y,al,jal,nlev,lev,ilev) 
      integer n, nlev, jal(*), ilev(nlev+1), lev(n)
      complex*16 x(n), y(n), al(*)
c-----------------------------------------------------------------------
c    Solves L x = y    L = triangular. Uses LEVEL SCHEDULING/MSR format 
c-----------------------------------------------------------------------
c
c On entry:
c---------- 
c n      = integer. dimension of problem.
c y      = real array containg the right hand side.
c
c al,
c jal,   = Lower triangular matrix stored in Modified Sparse Row 
c          format. 
c nlev   = number of levels in matrix
c lev    = integer array of length n, containing the permutation
c          that defines the levels in the level scheduling ordering.
c ilev   = pointer to beginning of levels in lev.
c          the numbers lev(i) to lev(i+1)-1 contain the row numbers
c          that belong to level number i, in the level shcheduling
c          ordering.
c
c On return:
c----------- 
c	x = The solution of  L x = y .
c--------------------------------------------------------------------
      integer ii, jrow, i 
      complex*16 t 
c     
c     outer loop goes through the levels. (SEQUENTIAL loop)
c     
      do 150 ii=1, nlev
c     
c     next loop executes within the same level. PARALLEL loop
c     
         do 100 i=ilev(ii), ilev(ii+1)-1 
            jrow = lev(i)
c
c compute inner product of row jrow with x
c 
            t = y(jrow) 
            do 130 k=jal(jrow), jal(jrow+1)-1 
               t = t - al(k)*x(jal(k))
 130        continue
            x(jrow) = t*al(jrow) 
 100     continue
 150  continue
      return
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine zusol (n,x,y,au,jau,iau)
      integer n, jau(*),iau(n+1) 
      complex*16  x(n), y(n), au(*) 
c----------------------------------------------------------------------- 
c             Solves   U x = y    U = unit upper triangular. 
c-----------------------------------------------------------------------
c solves a unit upper triangular system by standard (sequential )
c backward elimination - matrix stored in CSR format. 
c-----------------------------------------------------------------------
c
c On entry:
c---------- 
c n      = integer. dimension of problem.
c y      = real array containg the right side.
c
c au,
c jau,
c iau,    = Lower triangular matrix stored in compressed sparse row
c          format. 
c
c On return:
c----------- 
c	x = The solution of  U x = y . 
c-------------------------------------------------------------------- 
c local variables 
c
      integer k, j 
      complex*16  t
c-----------------------------------------------------------------------
      x(n) = y(n) 
      do 150 k = n-1,1,-1 
         t = y(k) 
         do 100 j = iau(k), iau(k+1)-1
            t = t - au(j)*x(jau(j))
 100     continue
         x(k) = t 
 150  continue
c
      return
c----------end-of-usol-------------------------------------------------- 
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine zudsol (n,x,y,au,jau) 
      integer n, jau(*) 
      complex*16  x(n), y(n),au(*) 
c----------------------------------------------------------------------- 
c             Solves   U x = y  ;   U = upper triangular in MSR format
c-----------------------------------------------------------------------
c solves a non-unit upper triangular matrix by standard (sequential )
c backward elimination - matrix stored in MSR format. 
c with diagonal elements already inverted (otherwise do inversion,
c au(1:n) = 1.0/au(1:n),  before calling).
c-----------------------------------------------------------------------
c
c On entry:
c---------- 
c n      = integer. dimension of problem.
c y      = real array containg the right side.
c
c au,
c jau,    = Lower triangular matrix stored in modified sparse row
c          format. 
c
c On return:
c----------- 
c	x = The solution of  U x = y .
c--------------------------------------------------------------------
c local variables 
c
      integer k, j
      complex*16 t
c-----------------------------------------------------------------------
      x(n) = y(n)*au(n)
      do 150 k = n-1,1,-1
         t = y(k) 
         do 100 j = jau(k), jau(k+1)-1
            t = t - au(j)*x(jau(j))
 100     continue
         x(k) = au(k)*t 
 150  continue
c
      return
c----------end-of-udsol-------------------------------------------------
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine zusolc (n,x,y,au,jau,iau)
      complex*16  x(*), y(*), au(*) 
      integer n, jau(*),iau(*)
c-----------------------------------------------------------------------
c       SOUVES     U x = y ;    where U = unit upper trang. CSC format
c-----------------------------------------------------------------------
c solves a unit upper triangular system by standard (sequential )
c forward elimination - matrix stored in CSC format. 
c-----------------------------------------------------------------------
c
c On entry:
c---------- 
c n      = integer. dimension of problem.
c y      = complex*16 array containg the right side.
c
c au,
c jau,
c iau,    = Uower triangular matrix stored in compressed sparse column 
c          format. 
c
c On return:
c----------- 
c	x  = The solution of  U x  = y.
c-----------------------------------------------------------------------
c local variables 
c     
      integer k, j
      complex*16 t
c-----------------------------------------------------------------------
      do 140 k=1,n
         x(k) = y(k) 
 140  continue
      do 150 k = n,1,-1
         t = x(k) 
         do 100 j = iau(k), iau(k+1)-1
            x(jau(j)) = x(jau(j)) - t*au(j) 
 100     continue
 150  continue
c
      return
c----------end-of-usolc------------------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine zudsolc (n,x,y,au,jau)   
      integer n, jau(*) 
      complex*16 x(n), y(n), au(*)  
c-----------------------------------------------------------------------
c    Solves     U x = y ;    U = nonunit Up. Triang. MSC format 
c----------------------------------------------------------------------- 
c solves a (non-unit) upper triangular system by standard (sequential) 
c forward elimination - matrix stored in Modified Sparse Column format 
c with diagonal elements already inverted (otherwise do inversion,
c auuuul(1:n) = 1.0/au(1:n),  before calling ldsol).
c-----------------------------------------------------------------------
c
c On entry:
c---------- 
c n      = integer. dimension of problem.
c y      = complex*16 array containg the right hand side.
c
c au,
c jau,   = Upper triangular matrix stored in Modified Sparse Column
c          format.
c
c On return:
c----------- 
c	x = The solution of  U x = y .
c--------------------------------------------------------------------
c local variables 
c 
      integer k, j
      complex*16 t
c----------------------------------------------------------------------- 
      do 140 k=1,n
         x(k) = y(k) 
 140  continue
      do 150 k = n,1,-1
         x(k) = x(k)*au(k) 
         t = x(k) 
         do 100 j = jau(k), jau(k+1)-1
            x(jau(j)) = x(jau(j)) - t*au(j) 
 100     continue
 150  continue
c
      return
c----------end-of-udsolc------------------------------------------------ 
c-----------------------------------------------------------------------
      end
