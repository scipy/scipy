c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c                     UNARY SUBROUTINES MODULE                         c
c----------------------------------------------------------------------c
c contents:                                                            c
c----------                                                            c
c submat : extracts a submatrix from a sparse matrix.                  c
c filter : filters elements from a matrix according to their magnitude.c
c filterm: same as above, but for the MSR format                       c
c csort  : sorts the elements in increasing order of columns           c
c clncsr : clean up the CSR format matrix, remove duplicate entry, etc c
c transp : in-place transposition routine (see also csrcsc in formats) c
c copmat : copy of a matrix into another matrix (both stored csr)      c
c msrcop : copies a matrix in MSR format into a matrix in MSR format   c
c getelm : returns a(i,j) for any (i,j) from a CSR-stored matrix.      c
c getdia : extracts a specified diagonal from a matrix.                c
c getl   : extracts lower triangular part                              c
c getu   : extracts upper triangular part                              c
c levels : gets the level scheduling structure for lower triangular    c
c          matrices.                                                   c
c amask  : extracts     C = A mask M                                   c
c rperm  : permutes the rows of a matrix (B = P A)                     c
c cperm  : permutes the columns of a matrix (B = A Q)                  c
c dperm  : permutes both the rows and columns of a matrix (B = P A Q ) c
c dperm1 : general extractiob routine (extracts arbitrary rows)        c
c dperm2 : general submatrix permutation/extraction routine            c
c dmperm : symmetric permutation of row and column (B=PAP') in MSR fmt c
c dvperm : permutes a real vector (in-place)                           c
c ivperm : permutes an integer vector (in-place)                       c
c retmx  : returns the max absolute value in each row of the matrix    c
c diapos : returns the positions of the diagonal elements in A.        c
c extbdg : extracts the main diagonal blocks of a matrix.              c
c getbwd : returns the bandwidth information on a matrix.              c
c blkfnd : finds the block-size of a matrix.                           c
c blkchk : checks whether a given integer is the block size of A.      c
c infdia : obtains information on the diagonals of A.                  c
c amubdg : gets number of nonzeros in each row of A*B (as well as NNZ) c 
c aplbdg : gets number of nonzeros in each row of A+B (as well as NNZ) c
c rnrms  : computes the norms of the rows of A                         c
c cnrms  : computes the norms of the columns of A                      c
c roscal : scales the rows of a matrix by their norms.                 c
c coscal : scales the columns of a matrix by their norms.              c
c addblk : Adds a matrix B into a block of A.                          c
c get1up : Collects the first elements of each row of the upper        c
c          triangular portion of the matrix.                           c
c xtrows : extracts given rows from a matrix in CSR format.            c
c csrkvstr:  Finds block row partitioning of matrix in CSR format      c
c csrkvstc:  Finds block column partitioning of matrix in CSR format   c
c kvstmerge: Merges block partitionings, for conformal row/col pattern c
c----------------------------------------------------------------------c
      subroutine ssubmat (n,job,i1,i2,j1,j2,a,ja,ia,nr,nc,ao,jao,iao)
      integer n,job,i1,i2,j1,j2,nr,nc,ia(*),ja(*),jao(*),iao(*)
      real*4 a(*),ao(*) 
c-----------------------------------------------------------------------
c extracts the submatrix A(i1:i2,j1:j2) and puts the result in 
c matrix ao,iao,jao
c---- In place: ao,jao,iao may be the same as a,ja,ia.
c-------------- 
c on input
c---------
c n	= row dimension of the matrix 
c i1,i2 = two integers with i2 .ge. i1 indicating the range of rows to be
c          extracted. 
c j1,j2 = two integers with j2 .ge. j1 indicating the range of columns 
c         to be extracted.
c         * There is no checking whether the input values for i1, i2, j1,
c           j2 are between 1 and n. 
c a,
c ja,
c ia    = matrix in compressed sparse row format. 
c
c job	= job indicator: if job .ne. 1 then the real values in a are NOT
c         extracted, only the column indices (i.e. data structure) are.
c         otherwise values as well as column indices are extracted...
c         
c on output
c-------------- 
c nr	= number of rows of submatrix 
c nc	= number of columns of submatrix 
c	  * if either of nr or nc is nonpositive the code will quit.
c
c ao,
c jao,iao = extracted matrix in general sparse format with jao containing
c	the column indices,and iao being the pointer to the beginning 
c	of the row,in arrays a,ja.
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c
      nr = i2-i1+1
      nc = j2-j1+1
c     
      if ( nr .le. 0 .or. nc .le. 0) return
c     
      klen = 0
c     
c     simple procedure. proceeds row-wise...
c     
      do 100 i = 1,nr
         ii = i1+i-1
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         iao(i) = klen+1
c-----------------------------------------------------------------------
         do 60 k=k1,k2
            j = ja(k)
            if (j .ge. j1 .and. j .le. j2) then
               klen = klen+1
               if (job .eq. 1) ao(klen) = a(k)
               jao(klen) = j - j1+1
            endif
 60      continue
 100  continue
      iao(nr+1) = klen+1
      return
c------------end-of submat---------------------------------------------- 
c----------------------------------------------------------------------- 
      end
c-----------------------------------------------------------------------
      subroutine sfilter (n,job,drptol,a,ja,ia,b,jb,ib,len,ierr)
      real*4 a(*),b(*)
      real*4 drptol
      integer ja(*),jb(*),ia(*),ib(*),n,job,len,ierr
c-----------------------------------------------------------------------
c     This module removes any elements whose absolute value
c     is small from an input matrix A and puts the resulting
c     matrix in B.  The input parameter job selects a definition
c     of small.
c-----------------------------------------------------------------------
c on entry:
c---------
c  n	 = integer. row dimension of matrix
c  job   = integer. used to determine strategy chosen by caller to
c         drop elements from matrix A. 
c          job = 1  
c              Elements whose absolute value is less than the
c              drop tolerance are removed.
c          job = 2
c              Elements whose absolute value is less than the 
c              product of the drop tolerance and the Euclidean
c              norm of the row are removed. 
c          job = 3
c              Elements whose absolute value is less that the
c              product of the drop tolerance and the largest
c              element in the row are removed.
c 
c drptol = real. drop tolerance used for dropping strategy.
c a	
c ja
c ia     = input matrix in compressed sparse format
c len	 = integer. the amount of space available in arrays b and jb.
c
c on return:
c---------- 
c b	
c jb
c ib    = resulting matrix in compressed sparse format. 
c 
c ierr	= integer. containing error message.
c         ierr .eq. 0 indicates normal return
c         ierr .gt. 0 indicates that there is'nt enough
c         space is a and ja to store the resulting matrix.
c         ierr then contains the row number where filter stopped.
c note:
c------ This module is in place. (b,jb,ib can ne the same as 
c       a, ja, ia in which case the result will be overwritten).
c----------------------------------------------------------------------c
c           contributed by David Day,  Sep 19, 1989.                   c
c----------------------------------------------------------------------c
c local variables
      real*4 norm,loctol
      integer index,row,k,k1,k2 
c
      index = 1
      do 10 row= 1,n
         k1 = ia(row)
         k2 = ia(row+1) - 1
         ib(row) = index
	 goto (100,200,300) job
 100     norm = 1.0d0
         goto 400
 200     norm = 0.0d0
         do 22 k = k1,k2
            norm = norm + a(k) * a(k)
 22      continue
         norm = sqrt(norm)
         goto 400
 300     norm = 0.0d0
         do 23 k = k1,k2
            if( abs(a(k))  .gt. norm) then
               norm = abs(a(k))
            endif
 23      continue
 400     loctol = drptol * norm
	 do 30 k = k1,k2
	    if( abs(a(k)) .gt. loctol)then 
               if (index .gt. len) then
               ierr = row 
               return
            endif
            b(index) =  a(k)
            jb(index) = ja(k)
            index = index + 1
         endif
 30   continue
 10   continue
      ib(n+1) = index
      return
c--------------------end-of-filter -------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine sfilterm (n,job,drop,a,ja,b,jb,len,ierr)
      real*4 a(*),b(*)
      real*4 drop
      integer ja(*),jb(*),n,job,len,ierr
c-----------------------------------------------------------------------
c     This subroutine removes any elements whose absolute value
c     is small from an input matrix A. Same as filter but
c     uses the MSR format.
c-----------------------------------------------------------------------
c on entry:
c---------
c  n	 = integer. row dimension of matrix
c  job   = integer. used to determine strategy chosen by caller to
c         drop elements from matrix A. 
c          job = 1  
c              Elements whose absolute value is less than the
c              drop tolerance are removed.
c          job = 2
c              Elements whose absolute value is less than the 
c              product of the drop tolerance and the Euclidean
c              norm of the row are removed. 
c          job = 3
c              Elements whose absolute value is less that the
c              product of the drop tolerance and the largest
c              element in the row are removed.
c 
c drop = real. drop tolerance used for dropping strategy.
c a	
c ja     = input matrix in Modifief Sparse Row format
c len	 = integer. the amount of space in arrays b and jb.
c
c on return:
c---------- 
c
c b, jb = resulting matrix in Modifief Sparse Row format
c 
c ierr	= integer. containing error message.
c         ierr .eq. 0 indicates normal return
c         ierr .gt. 0 indicates that there is'nt enough
c         space is a and ja to store the resulting matrix.
c         ierr then contains the row number where filter stopped.
c note:
c------ This module is in place. (b,jb can ne the same as 
c       a, ja in which case the result will be overwritten).
c----------------------------------------------------------------------c
c           contributed by David Day,  Sep 19, 1989.                   c
c----------------------------------------------------------------------c
c local variables
c
      real*4 norm,loctol
      integer index,row,k,k1,k2 
c
      index = n+2
      do 10 row= 1,n
         k1 = ja(row)
         k2 = ja(row+1) - 1
         jb(row) = index
	 goto (100,200,300) job
 100     norm = 1.0d0
         goto 400
 200     norm = a(row)**2 
         do 22 k = k1,k2
            norm = norm + a(k) * a(k)
 22      continue
         norm = sqrt(norm)
         goto 400
 300     norm = abs(a(row)) 
         do 23 k = k1,k2
            norm = max(abs(a(k)),norm) 
 23      continue
 400     loctol = drop * norm
	 do 30 k = k1,k2
	    if( abs(a(k)) .gt. loctol)then 
               if (index .gt. len) then
                  ierr = row 
                  return
               endif
               b(index) =  a(k)
               jb(index) = ja(k)
               index = index + 1
            endif
 30      continue
 10   continue
      jb(n+1) = index
      return
c--------------------end-of-filterm-------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine scsort (n,a,ja,ia,iwork,values) 
      logical values
      integer n, ja(*), ia(n+1), iwork(*) 
      real*4 a(*) 
c-----------------------------------------------------------------------
c This routine sorts the elements of  a matrix (stored in Compressed
c Sparse Row Format) in increasing order of their column indices within 
c each row. It uses a form of bucket sort with a cost of O(nnz) where
c nnz = number of nonzero elements. 
c requires an integer work array of length 2*nnz.  
c-----------------------------------------------------------------------
c on entry:
c--------- 
c n     = the row dimension of the matrix
c a     = the matrix A in compressed sparse row format.
c ja    = the array of column indices of the elements in array a.
c ia    = the array of pointers to the rows. 
c iwork = integer work array of length max ( n+1, 2*nnz ) 
c         where nnz = (ia(n+1)-ia(1))  ) .
c values= logical indicating whether or not the real values a(*) must 
c         also be permuted. if (.not. values) then the array a is not
c         touched by csort and can be a dummy array. 
c 
c on return:
c----------
c the matrix stored in the structure a, ja, ia is permuted in such a
c way that the column indices are in increasing order within each row.
c iwork(1:nnz) contains the permutation used  to rearrange the elements.
c----------------------------------------------------------------------- 
c Y. Saad - Feb. 1, 1991.
c-----------------------------------------------------------------------
c local variables
      integer i, k, j, ifirst, nnz, next  
c
c count the number of elements in each column
c
      do 1 i=1,n+1
         iwork(i) = 0
 1    continue
      do 3 i=1, n
         do 2 k=ia(i), ia(i+1)-1 
            j = ja(k)+1
            iwork(j) = iwork(j)+1
 2       continue 
 3    continue
c
c compute pointers from lengths. 
c
      iwork(1) = 1
      do 4 i=1,n
         iwork(i+1) = iwork(i) + iwork(i+1)
 4    continue
c 
c get the positions of the nonzero elements in order of columns.
c
      ifirst = ia(1) 
      nnz = ia(n+1)-ifirst
      do 5 i=1,n
         do 51 k=ia(i),ia(i+1)-1 
            j = ja(k) 
            next = iwork(j) 
            iwork(nnz+next) = k
            iwork(j) = next+1
 51      continue
 5    continue
c
c convert to coordinate format
c 
      do 6 i=1, n
         do 61 k=ia(i), ia(i+1)-1 
            iwork(k) = i
 61      continue
 6    continue
c
c loop to find permutation: for each element find the correct 
c position in (sorted) arrays a, ja. Record this in iwork. 
c 
      do 7 k=1, nnz
         ko = iwork(nnz+k) 
         irow = iwork(ko)
         next = ia(irow)
c
c the current element should go in next position in row. iwork
c records this position. 
c 
         iwork(ko) = next
         ia(irow)  = next+1
 7       continue
c
c perform an in-place permutation of the  arrays.
c 
         call ivperm (nnz, ja(ifirst), iwork) 
         if (values) call sdvperm (nnz, a(ifirst), iwork) 
c
c reshift the pointers of the original matrix back.
c 
      do 8 i=n,1,-1
         ia(i+1) = ia(i)
 8    continue
      ia(1) = ifirst 
c
      return 
c---------------end-of-csort-------------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine sclncsr (job,value2,nrow,a,ja,ia,indu,iwk)
c     .. Scalar Arguments ..
      integer job, nrow, value2
c     ..
c     .. Array Arguments ..
      integer ia(nrow+1),indu(nrow),iwk(nrow+1),ja(*)
      real*4  a(*)
c     ..
c
c     This routine performs two tasks to clean up a CSR matrix
c     -- remove duplicate/zero entries,
c     -- perform a partial ordering, new order lower triangular part,
c        main diagonal, upper triangular part.
c
c     On entry:
c
c     job   = options
c         0 -- nothing is done
c         1 -- eliminate duplicate entries, zero entries.
c         2 -- eliminate duplicate entries and perform partial ordering.
c         3 -- eliminate duplicate entries, sort the entries in the
c              increasing order of clumn indices.
c
c     value2  -- 0 the matrix is pattern only (a is not touched)
c                1 matrix has values too.
c     nrow    -- row dimension of the matrix
c     a,ja,ia -- input matrix in CSR format
c
c     On return:
c     a,ja,ia -- cleaned matrix
c     indu    -- pointers to the beginning of the upper triangular
c                portion if job > 1
c
c     Work space:
c     iwk     -- integer work space of size nrow+1
c
c     .. Local Scalars ..
      integer i,j,k,ko,ipos,kfirst,klast
      real*4  tmp
c     ..
c
      if (job.le.0) return
c
c     .. eliminate duplicate entries --
c     array INDU is used as marker for existing indices, it is also the
c     location of the entry.
c     IWK is used to stored the old IA array.
c     matrix is copied to squeeze out the space taken by the duplicated
c     entries.
c
      do 90 i = 1, nrow
         indu(i) = 0
         iwk(i) = ia(i)
 90   continue
      iwk(nrow+1) = ia(nrow+1)
      k = 1
      do 120 i = 1, nrow
         ia(i) = k
         ipos = iwk(i)
         klast = iwk(i+1)
 100     if (ipos.lt.klast) then
            j = ja(ipos)
            if (indu(j).eq.0) then
c     .. new entry ..
               if (value2.ne.0) then
                  if (a(ipos) .ne. 0.0D0) then
                     indu(j) = k
                     ja(k) = ja(ipos)
                     a(k) = a(ipos)
                     k = k + 1
                  endif
               else
                  indu(j) = k
                  ja(k) = ja(ipos)
                  k = k + 1
               endif
            else if (value2.ne.0) then
c     .. duplicate entry ..
               a(indu(j)) = a(indu(j)) + a(ipos)
            endif
            ipos = ipos + 1
            go to 100
         endif
c     .. remove marks before working on the next row ..
         do 110 ipos = ia(i), k - 1
            indu(ja(ipos)) = 0
 110     continue
 120  continue
      ia(nrow+1) = k
      if (job.le.1) return
c
c     .. partial ordering ..
c     split the matrix into strict upper/lower triangular
c     parts, INDU points to the the beginning of the upper part.
c
      do 140 i = 1, nrow
         klast = ia(i+1) - 1
         kfirst = ia(i)
 130     if (klast.gt.kfirst) then
            if (ja(klast).lt.i .and. ja(kfirst).ge.i) then
c     .. swap klast with kfirst ..
               j = ja(klast)
               ja(klast) = ja(kfirst)
               ja(kfirst) = j
               if (value2.ne.0) then
                  tmp = a(klast)
                  a(klast) = a(kfirst)
                  a(kfirst) = tmp
               endif
            endif
            if (ja(klast).ge.i)
     &         klast = klast - 1
            if (ja(kfirst).lt.i)
     &         kfirst = kfirst + 1
            go to 130
         endif
c
         if (ja(klast).lt.i) then
            indu(i) = klast + 1
         else
            indu(i) = klast
         endif
 140  continue
      if (job.le.2) return
c
c     .. order the entries according to column indices
c     burble-sort is used
c
      do 190 i = 1, nrow
         do 160 ipos = ia(i), indu(i)-1
            do 150 j = indu(i)-1, ipos+1, -1
               k = j - 1
               if (ja(k).gt.ja(j)) then
                  ko = ja(k)
                  ja(k) = ja(j)
                  ja(j) = ko
                  if (value2.ne.0) then
                     tmp = a(k)
                     a(k) = a(j)
                     a(j) = tmp
                  endif
               endif
 150        continue
 160     continue
         do 180 ipos = indu(i), ia(i+1)-1
            do 170 j = ia(i+1)-1, ipos+1, -1
               k = j - 1
               if (ja(k).gt.ja(j)) then
                  ko = ja(k)
                  ja(k) = ja(j)
                  ja(j) = ko
                  if (value2.ne.0) then
                     tmp = a(k)
                     a(k) = a(j)
                     a(j) = tmp
                  endif
               endif
 170        continue
 180     continue
 190  continue
      return
c---- end of clncsr ----------------------------------------------------
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine scopmat (nrow,a,ja,ia,ao,jao,iao,ipos,job) 
      real*4 a(*),ao(*) 
      integer nrow, ia(*),ja(*),jao(*),iao(*), ipos, job 
c----------------------------------------------------------------------
c copies the matrix a, ja, ia, into the matrix ao, jao, iao. 
c----------------------------------------------------------------------
c on entry:
c---------
c nrow	= row dimension of the matrix 
c a,
c ja,
c ia    = input matrix in compressed sparse row format. 
c ipos  = integer. indicates the position in the array ao, jao
c         where the first element should be copied. Thus 
c         iao(1) = ipos on return. 
c job   = job indicator. if (job .ne. 1) the values are not copies 
c         (i.e., pattern only is copied in the form of arrays ja, ia).
c
c on return:
c----------
c ao,
c jao,
c iao   = output matrix containing the same data as a, ja, ia.
c-----------------------------------------------------------------------
c           Y. Saad, March 1990. 
c-----------------------------------------------------------------------
c local variables
      integer kst, i, k 
c
      kst    = ipos -ia(1) 
      do 100 i = 1, nrow+1
         iao(i) = ia(i) + kst
 100  continue
c     
      do 200 k=ia(1), ia(nrow+1)-1
         jao(kst+k)= ja(k)
 200  continue
c
      if (job .ne. 1) return
      do 201 k=ia(1), ia(nrow+1)-1
         ao(kst+k) = a(k)
 201  continue
c
      return
c--------end-of-copmat -------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine smsrcop (nrow,a,ja,ao,jao,job) 
      real*4 a(*),ao(*) 
      integer nrow, ja(*),jao(*), job 
c----------------------------------------------------------------------
c copies the MSR matrix a, ja, into the MSR matrix ao, jao 
c----------------------------------------------------------------------
c on entry:
c---------
c nrow	= row dimension of the matrix 
c a,ja  = input matrix in Modified compressed sparse row format. 
c job   = job indicator. Values are not copied if job .ne. 1 
c       
c on return:
c----------
c ao, jao   = output matrix containing the same data as a, ja.
c-----------------------------------------------------------------------
c           Y. Saad, 
c-----------------------------------------------------------------------
c local variables
      integer i, k 
c
      do 100 i = 1, nrow+1
         jao(i) = ja(i) 
 100  continue
c     
      do 200 k=ja(1), ja(nrow+1)-1
         jao(k)= ja(k)
 200  continue
c
      if (job .ne. 1) return
      do 201 k=ja(1), ja(nrow+1)-1
         ao(k) = a(k)
 201  continue
      do 202 k=1,nrow
         ao(k) = a(k)
 202  continue
c
      return
c--------end-of-msrcop -------------------------------------------------
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      real*4 function sgetelm (i,j,a,ja,ia,iadd,sorted) 
c-----------------------------------------------------------------------
c     purpose:
c     -------- 
c     this function returns the element a(i,j) of a matrix a, 
c     for any pair (i,j).  the matrix is assumed to be stored 
c     in compressed sparse row (csr) format. getelm performs a
c     binary search in the case where it is known that the elements 
c     are sorted so that the column indices are in increasing order. 
c     also returns (in iadd) the address of the element a(i,j) in 
c     arrays a and ja when the search is successsful (zero if not).
c----- 
c     first contributed by noel nachtigal (mit). 
c     recoded jan. 20, 1991, by y. saad [in particular
c     added handling of the non-sorted case + the iadd output] 
c-----------------------------------------------------------------------
c     parameters:
c     ----------- 
c on entry: 
c---------- 
c     i      = the row index of the element sought (input).
c     j      = the column index of the element sought (input).
c     a      = the matrix a in compressed sparse row format (input).
c     ja     = the array of column indices (input).
c     ia     = the array of pointers to the rows' data (input).
c     sorted = logical indicating whether the matrix is knonw to 
c              have its column indices sorted in increasing order 
c              (sorted=.true.) or not (sorted=.false.).
c              (input). 
c on return:
c----------- 
c     getelm = value of a(i,j). 
c     iadd   = address of element a(i,j) in arrays a, ja if found,
c              zero if not found. (output) 
c
c     note: the inputs i and j are not checked for validity. 
c-----------------------------------------------------------------------
c     noel m. nachtigal october 28, 1990 -- youcef saad jan 20, 1991.
c----------------------------------------------------------------------- 
      integer i, ia(*), iadd, j, ja(*)
      real*4 a(*)
      logical sorted 
c
c     local variables.
c
      integer ibeg, iend, imid, k
c
c     initialization 
c
      iadd = 0 
      getelm = 0.0
      ibeg = ia(i)
      iend = ia(i+1)-1
c
c     case where matrix is not necessarily sorted
c     
      if (.not. sorted) then 
c
c scan the row - exit as soon as a(i,j) is found
c
         do 5  k=ibeg, iend
            if (ja(k) .eq.  j) then
               iadd = k 
               goto 20 
            endif
 5       continue
c     
c     end unsorted case. begin sorted case
c     
      else
c     
c     begin binary search.   compute the middle index.
c     
 10      imid = ( ibeg + iend ) / 2
c     
c     test if  found
c     
         if (ja(imid).eq.j) then
            iadd = imid 
            goto 20
         endif
         if (ibeg .ge. iend) goto 20
c     
c     else     update the interval bounds. 
c     
         if (ja(imid).gt.j) then
            iend = imid -1
         else 
            ibeg = imid +1
         endif
         goto 10  
c     
c     end both cases
c     
      endif
c     
 20   if (iadd .ne. 0) getelm = a(iadd) 
c
      return
c--------end-of-getelm--------------------------------------------------
c-----------------------------------------------------------------------
      end 
c-----------------------------------------------------------------------
      subroutine sgetdia (nrow,ncol,job,a,ja,ia,len,diag,idiag,ioff)
      real*4 diag(*),a(*)
      integer nrow, ncol, job, len, ioff, ia(*), ja(*), idiag(*)
c-----------------------------------------------------------------------
c this subroutine extracts a given diagonal from a matrix stored in csr 
c format. the output matrix may be transformed with the diagonal removed
c from it if desired (as indicated by job.) 
c----------------------------------------------------------------------- 
c our definition of a diagonal of matrix is a vector of length nrow
c (always) which contains the elements in rows 1 to nrow of
c the matrix that are contained in the diagonal offset by ioff
c with respect to the main diagonal. if the diagonal element
c falls outside the matrix then it is defined as a zero entry.
c thus the proper definition of diag(*) with offset ioff is 
c
c     diag(i) = a(i,ioff+i) i=1,2,...,nrow
c     with elements falling outside the matrix being defined as zero.
c 
c----------------------------------------------------------------------- 
c 
c on entry:
c---------- 
c
c nrow	= integer. the row dimension of the matrix a.
c ncol	= integer. the column dimension of the matrix a.
c job   = integer. job indicator.  if job = 0 then
c         the matrix a, ja, ia, is not altered on return.
c         if job.ne.0  then getdia will remove the entries
c         collected in diag from the original matrix.
c         this is done in place.
c
c a,ja,
c    ia = matrix stored in compressed sparse row a,ja,ia,format
c ioff  = integer,containing the offset of the wanted diagonal
c	  the diagonal extracted is the one corresponding to the
c	  entries a(i,j) with j-i = ioff.
c	  thus ioff = 0 means the main diagonal
c
c on return:
c----------- 
c len   = number of nonzero elements found in diag.
c         (len .le. min(nrow,ncol-ioff)-max(1,1-ioff) + 1 )
c
c diag  = real*4 array of length nrow containing the wanted diagonal.
c	  diag contains the diagonal (a(i,j),j-i = ioff ) as defined 
c         above. 
c
c idiag = integer array of  length len, containing the poisitions 
c         in the original arrays a and ja of the diagonal elements
c         collected in diag. a zero entry in idiag(i) means that 
c         there was no entry found in row i belonging to the diagonal.
c         
c a, ja,
c    ia = if job .ne. 0 the matrix is unchanged. otherwise the nonzero
c         diagonal entries collected in diag are removed from the 
c         matrix and therefore the arrays a, ja, ia will change.
c	  (the matrix a, ja, ia will contain len fewer elements) 
c 
c----------------------------------------------------------------------c
c     Y. Saad, sep. 21 1989 - modified and retested Feb 17, 1996.      c 
c----------------------------------------------------------------------c
c     local variables
      integer istart, max, iend, i, kold, k, kdiag, ko
c     
      istart = max(0,-ioff)
      iend = min(nrow,ncol-ioff)
      len = 0
      do 1 i=1,nrow
         idiag(i) = 0
	 diag(i) = 0.0d0 
 1    continue
c     
c     extract  diagonal elements
c     
      do 6 i=istart+1, iend
         do 51 k= ia(i),ia(i+1) -1
            if (ja(k)-i .eq. ioff) then
               diag(i)= a(k)
               idiag(i) = k
               len = len+1
               goto 6
            endif
 51      continue
 6    continue
      if (job .eq. 0 .or. len .eq.0) return
c
c     remove diagonal elements and rewind structure
c 
      ko = 0
      do  7 i=1, nrow 
         kold = ko
         kdiag = idiag(i) 
         do 71 k= ia(i), ia(i+1)-1 
            if (k .ne. kdiag) then
               ko = ko+1
               a(ko) = a(k)
               ja(ko) = ja(k)
            endif
 71      continue
         ia(i) = kold+1
 7    continue
c
c     redefine ia(nrow+1)
c
      ia(nrow+1) = ko+1
      return
c------------end-of-getdia----------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine stransp (nrow,ncol,a,ja,ia,iwk,ierr)
      integer nrow, ncol, ia(*), ja(*), iwk(*), ierr
      real*4 a(*) 
c------------------------------------------------------------------------
c In-place transposition routine.
c------------------------------------------------------------------------
c this subroutine transposes a matrix stored in compressed sparse row 
c format. the transposition is done in place in that the arrays a,ja,ia
c of the transpose are overwritten onto the original arrays.
c------------------------------------------------------------------------
c on entry:
c--------- 
c nrow	= integer. The row dimension of A.
c ncol	= integer. The column dimension of A.
c a	= real array of size nnz (number of nonzero elements in A).
c         containing the nonzero elements 
c ja	= integer array of length nnz containing the column positions
c 	  of the corresponding elements in a.
c ia	= integer of size n+1, where n = max(nrow,ncol). On entry
c         ia(k) contains the position in a,ja of  the beginning of 
c         the k-th row.
c
c iwk	= integer work array of same length as ja.
c
c on return:
c----------
c
c ncol	= actual row dimension of the transpose of the input matrix.
c         Note that this may be .le. the input value for ncol, in
c         case some of the last columns of the input matrix are zero
c         columns. In the case where the actual number of rows found
c         in transp(A) exceeds the input value of ncol, transp will
c         return without completing the transposition. see ierr.
c a,
c ja,
c ia	= contains the transposed matrix in compressed sparse
c         row format. The row dimension of a, ja, ia is now ncol.
c
c ierr	= integer. error message. If the number of rows for the
c         transposed matrix exceeds the input value of ncol,
c         then ierr is  set to that number and transp quits.
c         Otherwise ierr is set to 0 (normal return).
c
c Note: 
c----- 1) If you do not need the transposition to be done in place
c         it is preferrable to use the conversion routine csrcsc 
c         (see conversion routines in formats).
c      2) the entries of the output matrix are not sorted (the column
c         indices in each are not in increasing order) use csrcsc
c         if you want them sorted.
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c  modified Oct. 11, 1989.                                             c
c----------------------------------------------------------------------c
c local variables
      real*4 t, t1
      ierr = 0
      nnz = ia(nrow+1)-1
c
c     determine column dimension
c
      jcol = 0
      do 1 k=1, nnz
         jcol = max(jcol,ja(k))
 1    continue
      if (jcol .gt. ncol) then
         ierr = jcol
         return
      endif
c     
c     convert to coordinate format. use iwk for row indices.
c     
      ncol = jcol
c     
      do 3 i=1,nrow
         do 2 k=ia(i),ia(i+1)-1
            iwk(k) = i
 2       continue 
 3    continue
c     find pointer array for transpose. 
      do 35 i=1,ncol+1
         ia(i) = 0
 35   continue
      do 4 k=1,nnz
         i = ja(k)
         ia(i+1) = ia(i+1)+1
 4    continue 
      ia(1) = 1 
c------------------------------------------------------------------------
      do 44 i=1,ncol
         ia(i+1) = ia(i) + ia(i+1)
 44   continue 
c     
c     loop for a cycle in chasing process. 
c     
      init = 1
      k = 0
 5    t = a(init)
      i = ja(init)
      j = iwk(init)
      iwk(init) = -1
c------------------------------------------------------------------------
 6    k = k+1 		
c     current row number is i.  determine  where to go. 
      l = ia(i)
c     save the chased element. 
      t1 = a(l)
      inext = ja(l)
c     then occupy its location.
      a(l)  = t
      ja(l) = j
c     update pointer information for next element to be put in row i. 
      ia(i) = l+1
c     determine  next element to be chased
      if (iwk(l) .lt. 0) goto 65
      t = t1
      i = inext
      j = iwk(l)
      iwk(l) = -1
      if (k .lt. nnz) goto 6
      goto 70
 65   init = init+1
      if (init .gt. nnz) goto 70
      if (iwk(init) .lt. 0) goto 65
c     restart chasing --	
      goto 5
 70   continue
      do 80 i=ncol,1,-1 
         ia(i+1) = ia(i)
 80   continue
      ia(1) = 1
c
      return
c------------------end-of-transp ----------------------------------------
c------------------------------------------------------------------------
      end 
c------------------------------------------------------------------------ 
      subroutine sgetl (n,a,ja,ia,ao,jao,iao)
      integer n, ia(*), ja(*), iao(*), jao(*)
      real*4 a(*), ao(*)
c------------------------------------------------------------------------
c this subroutine extracts the lower triangular part of a matrix 
c and writes the result ao, jao, iao. The routine is in place in
c that ao, jao, iao can be the same as a, ja, ia if desired.
c-----------
c on input:
c
c n     = dimension of the matrix a.
c a, ja, 
c    ia = matrix stored in compressed sparse row format.
c On return:
c ao, jao, 
c    iao = lower triangular matrix (lower part of a) 
c	stored in a, ja, ia, format
c note: the diagonal element is the last element in each row.
c i.e. in  a(ia(i+1)-1 ) 
c ao, jao, iao may be the same as a, ja, ia on entry -- in which case
c getl will overwrite the result on a, ja, ia.
c
c------------------------------------------------------------------------
c local variables
      real*4 t
      integer ko, kold, kdiag, k, i
c
c inititialize ko (pointer for output matrix)
c
      ko = 0
      do  7 i=1, n
         kold = ko
         kdiag = 0
         do 71 k = ia(i), ia(i+1) -1
            if (ja(k)  .gt. i) goto 71
            ko = ko+1
            ao(ko) = a(k)
            jao(ko) = ja(k)
            if (ja(k)  .eq. i) kdiag = ko
 71      continue
         if (kdiag .eq. 0 .or. kdiag .eq. ko) goto 72
c
c     exchange
c
         t = ao(kdiag)
         ao(kdiag) = ao(ko)
         ao(ko) = t
c     
         k = jao(kdiag)
         jao(kdiag) = jao(ko)
         jao(ko) = k
 72      iao(i) = kold+1
 7    continue
c     redefine iao(n+1)
      iao(n+1) = ko+1
      return
c----------end-of-getl ------------------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine sgetu (n,a,ja,ia,ao,jao,iao)
      integer n, ia(*), ja(*), iao(*), jao(*)
      real*4 a(*), ao(*)
c------------------------------------------------------------------------
c this subroutine extracts the upper triangular part of a matrix 
c and writes the result ao, jao, iao. The routine is in place in
c that ao, jao, iao can be the same as a, ja, ia if desired.
c-----------
c on input:
c
c n     = dimension of the matrix a.
c a, ja, 
c    ia = matrix stored in a, ja, ia, format
c On return:
c ao, jao, 
c    iao = upper triangular matrix (upper part of a) 
c	stored in compressed sparse row format
c note: the diagonal element is the last element in each row.
c i.e. in  a(ia(i+1)-1 ) 
c ao, jao, iao may be the same as a, ja, ia on entry -- in which case
c getu will overwrite the result on a, ja, ia.
c
c------------------------------------------------------------------------
c local variables
      real*4 t 
      integer ko, k, i, kdiag, kfirst 
      ko = 0
      do  7 i=1, n
         kfirst = ko+1
         kdiag = 0
         do 71 k = ia(i), ia(i+1) -1
            if (ja(k)  .lt. i) goto 71
            ko = ko+1
            ao(ko) = a(k)
            jao(ko) = ja(k)
            if (ja(k)  .eq. i) kdiag = ko
 71      continue
         if (kdiag .eq. 0 .or. kdiag .eq. kfirst) goto 72
c     exchange
         t = ao(kdiag)
         ao(kdiag) = ao(kfirst)
         ao(kfirst) = t
c     
         k = jao(kdiag)
         jao(kdiag) = jao(kfirst)
         jao(kfirst) = k
 72      iao(i) = kfirst
 7    continue
c     redefine iao(n+1)
      iao(n+1) = ko+1
      return
c----------end-of-getu ------------------------------------------------- 
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine samask (nrow,ncol,a,ja,ia,jmask,imask,
     *                  c,jc,ic,iw,nzmax,ierr)
c---------------------------------------------------------------------
      real*4 a(*),c(*) 
      integer ia(nrow+1),ja(*),jc(*),ic(nrow+1),jmask(*),imask(nrow+1) 
      logical iw(ncol)
c-----------------------------------------------------------------------
c This subroutine builds a sparse matrix from an input matrix by 
c extracting only elements in positions defined by the mask jmask, imask
c-----------------------------------------------------------------------
c On entry:
c---------
c nrow  = integer. row dimension of input matrix 
c ncol	= integer. Column dimension of input matrix.
c
c a,
c ja,
c ia	= matrix in Compressed Sparse Row format
c
c jmask,
c imask = matrix defining mask (pattern only) stored in compressed
c         sparse row format.
c
c nzmax = length of arrays c and jc. see ierr.
c 
c On return:
c-----------
c
c a, ja, ia and jmask, imask are unchanged.
c
c c
c jc, 
c ic	= the output matrix in Compressed Sparse Row format.
c 
c ierr  = integer. serving as error message.c
c         ierr = 1  means normal return
c         ierr .gt. 1 means that amask stopped when processing
c         row number ierr, because there was not enough space in
c         c, jc according to the value of nzmax.
c
c work arrays:
c------------- 
c iw	= logical work array of length ncol.
c
c note: 
c------ the  algorithm is in place: c, jc, ic can be the same as 
c a, ja, ia in which cas the code will overwrite the matrix c
c on a, ja, ia
c
c-----------------------------------------------------------------------
      ierr = 0
      len = 0
      do 1 j=1, ncol
         iw(j) = .false.
 1    continue
c     unpack the mask for row ii in iw
      do 100 ii=1, nrow
c     save pointer in order to be able to do things in place
         do 2 k=imask(ii), imask(ii+1)-1
            iw(jmask(k)) = .true.
 2       continue
c     add umasked elemnts of row ii
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         ic(ii) = len+1
         do 200 k=k1,k2 
            j = ja(k)
            if (iw(j)) then
               len = len+1
               if (len .gt. nzmax) then
                  ierr = ii
                  return
               endif
               jc(len) = j
               c(len) = a(k)
            endif
 200     continue	      
c     
         do 3 k=imask(ii), imask(ii+1)-1
            iw(jmask(k)) = .false.
 3       continue
 100  continue	  
      ic(nrow+1)=len+1
c
      return
c-----end-of-amask -----------------------------------------------------
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine srperm (nrow,a,ja,ia,ao,jao,iao,perm,job)
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow),job
      real*4 a(*),ao(*) 
c-----------------------------------------------------------------------
c this subroutine permutes the rows of a matrix in CSR format. 
c rperm  computes B = P A  where P is a permutation matrix.  
c the permutation P is defined through the array perm: for each j, 
c perm(j) represents the destination row number of row number j. 
c Youcef Saad -- recoded Jan 28, 1991.
c-----------------------------------------------------------------------
c on entry:
c----------
c n 	= dimension of the matrix
c a, ja, ia = input matrix in csr format
c perm 	= integer array of length nrow containing the permutation arrays
c	  for the rows: perm(i) is the destination of row i in the
c         permuted matrix. 
c         ---> a(i,j) in the original matrix becomes a(perm(i),j) 
c         in the output  matrix.
c
c job	= integer indicating the work to be done:
c 		job = 1	permute a, ja, ia into ao, jao, iao 
c                       (including the copying of real values ao and
c                       the array iao).
c 		job .ne. 1 :  ignore real values.
c                     (in which case arrays a and ao are not needed nor
c                      used).
c
c------------
c on return: 
c------------ 
c ao, jao, iao = input matrix in a, ja, ia format
c note : 
c        if (job.ne.1)  then the arrays a and ao are not used.
c----------------------------------------------------------------------c
c           Y. Saad, May  2, 1990                                      c
c----------------------------------------------------------------------c
      logical values
      values = (job .eq. 1) 
c     
c     determine pointers for output matix. 
c     
      do 50 j=1,nrow
         i = perm(j)
         iao(i+1) = ia(j+1) - ia(j)
 50   continue
c
c get pointers from lengths
c
      iao(1) = 1
      do 51 j=1,nrow
         iao(j+1)=iao(j+1)+iao(j)
 51   continue
c
c copying 
c
      do 100 ii=1,nrow
c
c old row = ii  -- new row = iperm(ii) -- ko = new pointer
c        
         ko = iao(perm(ii)) 
         do 60 k=ia(ii), ia(ii+1)-1 
            jao(ko) = ja(k) 
            if (values) ao(ko) = a(k)
            ko = ko+1
 60      continue
 100  continue
c
      return
c---------end-of-rperm ------------------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine scperm (nrow,a,ja,ia,ao,jao,iao,perm,job) 
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(*), job
      real*4 a(*), ao(*) 
c-----------------------------------------------------------------------
c this subroutine permutes the columns of a matrix a, ja, ia.
c the result is written in the output matrix  ao, jao, iao.
c cperm computes B = A P, where  P is a permutation matrix
c that maps column j into column perm(j), i.e., on return 
c      a(i,j) becomes a(i,perm(j)) in new matrix 
c Y. Saad, May 2, 1990 / modified Jan. 28, 1991. 
c-----------------------------------------------------------------------
c on entry:
c----------
c nrow 	= row dimension of the matrix
c
c a, ja, ia = input matrix in csr format. 
c
c perm	= integer array of length ncol (number of columns of A
c         containing the permutation array  the columns: 
c         a(i,j) in the original matrix becomes a(i,perm(j))
c         in the output matrix.
c
c job	= integer indicating the work to be done:
c 		job = 1	permute a, ja, ia into ao, jao, iao 
c                       (including the copying of real values ao and
c                       the array iao).
c 		job .ne. 1 :  ignore real values ao and ignore iao.
c
c------------
c on return: 
c------------ 
c ao, jao, iao = input matrix in a, ja, ia format (array ao not needed)
c
c Notes:
c------- 
c 1. if job=1 then ao, iao are not used.
c 2. This routine is in place: ja, jao can be the same. 
c 3. If the matrix is initially sorted (by increasing column number) 
c    then ao,jao,iao  may not be on return. 
c 
c----------------------------------------------------------------------c
c local parameters:
      integer k, i, nnz
c
      nnz = ia(nrow+1)-1
      do 100 k=1,nnz
         jao(k) = perm(ja(k)) 
 100  continue
c
c     done with ja array. return if no need to touch values.
c
      if (job .ne. 1) return
c
c else get new pointers -- and copy values too.
c 
      do 1 i=1, nrow+1
         iao(i) = ia(i)
 1    continue
c
      do 2 k=1, nnz
         ao(k) = a(k)
 2    continue
c
      return
c---------end-of-cperm-------------------------------------------------- 
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine sdperm (nrow,a,ja,ia,ao,jao,iao,perm,qperm,job)
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow),
     +        qperm(*),job
      real*4 a(*),ao(*) 
c-----------------------------------------------------------------------
c This routine permutes the rows and columns of a matrix stored in CSR
c format. i.e., it computes P A Q, where P, Q are permutation matrices. 
c P maps row i into row perm(i) and Q maps column j into column qperm(j): 
c      a(i,j)    becomes   a(perm(i),qperm(j)) in new matrix
c In the particular case where Q is the transpose of P (symmetric 
c permutation of A) then qperm is not needed. 
c note that qperm should be of length ncol (number of columns) but this
c is not checked. 
c-----------------------------------------------------------------------
c Y. Saad, Sep. 21 1989 / recoded Jan. 28 1991. 
c-----------------------------------------------------------------------
c on entry: 
c---------- 
c n 	= dimension of the matrix
c a, ja, 
c    ia = input matrix in a, ja, ia format
c perm 	= integer array of length n containing the permutation arrays
c	  for the rows: perm(i) is the destination of row i in the
c         permuted matrix -- also the destination of column i in case
c         permutation is symmetric (job .le. 2) 
c
c qperm	= same thing for the columns. This should be provided only
c         if job=3 or job=4, i.e., only in the case of a nonsymmetric
c	  permutation of rows and columns. Otherwise qperm is a dummy
c
c job	= integer indicating the work to be done:
c * job = 1,2 permutation is symmetric  Ao :== P * A * transp(P)
c 		job = 1	permute a, ja, ia into ao, jao, iao 
c 		job = 2 permute matrix ignoring real values.
c * job = 3,4 permutation is non-symmetric  Ao :== P * A * Q 
c 		job = 3	permute a, ja, ia into ao, jao, iao 
c 		job = 4 permute matrix ignoring real values.
c		
c on return: 
c-----------
c ao, jao, iao = input matrix in a, ja, ia format
c
c in case job .eq. 2 or job .eq. 4, a and ao are never referred to 
c and can be dummy arguments. 
c Notes:
c------- 
c  1) algorithm is in place 
c  2) column indices may not be sorted on return even  though they may be 
c     on entry.
c----------------------------------------------------------------------c
c local variables 
      integer locjob, mod
c
c     locjob indicates whether or not real values must be copied. 
c     
      locjob = mod(job,2) 
c
c permute rows first 
c 
      call srperm (nrow,a,ja,ia,ao,jao,iao,perm,locjob)
c
c then permute columns
c
      locjob = 0
c
      if (job .le. 2) then
         call scperm (nrow,ao,jao,iao,ao,jao,iao,perm,locjob) 
      else 
         call scperm (nrow,ao,jao,iao,ao,jao,iao,qperm,locjob) 
      endif 
c     
      return
c-------end-of-dperm----------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine sdperm1 (i1,i2,a,ja,ia,b,jb,ib,perm,ipos,job)
      integer i1,i2,job,ja(*),ia(*),jb(*),ib(*),perm(*)
      real*4 a(*),b(*)
c----------------------------------------------------------------------- 
c     general submatrix extraction routine.
c----------------------------------------------------------------------- 
c     extracts rows perm(i1), perm(i1+1), ..., perm(i2) (in this order) 
c     from a matrix (doing nothing in the column indices.) The resulting 
c     submatrix is constructed in b, jb, ib. A pointer ipos to the
c     beginning of arrays b,jb,is also allowed (i.e., nonzero elements
c     are accumulated starting in position ipos of b, jb). 
c-----------------------------------------------------------------------
c Y. Saad,Sep. 21 1989 / recoded Jan. 28 1991 / modified for PSPARSLIB 
c Sept. 1997.. 
c-----------------------------------------------------------------------
c on entry: 
c---------- 
c n 	= dimension of the matrix
c a,ja,
c   ia  = input matrix in CSR format
c perm 	= integer array of length n containing the indices of the rows
c         to be extracted. 
c
c job   = job indicator. if (job .ne.1) values are not copied (i.e.,
c         only pattern is copied).
c
c on return: 
c-----------
c b,ja,
c ib   = matrix in csr format. b(ipos:ipos+nnz-1),jb(ipos:ipos+nnz-1) 
c     contain the value and column indices respectively of the nnz
c     nonzero elements of the permuted matrix. thus ib(1)=ipos.
c
c Notes:
c------- 
c  algorithm is NOT in place 
c----------------------------------------------------------------------- 
c local variables
c
      integer ko,irow,k 
      logical values
c-----------------------------------------------------------------------
      values = (job .eq. 1) 
      ko = ipos 
      ib(1) = ko
      do 900 i=i1,i2
         irow = perm(i) 
         do 800 k=ia(irow),ia(irow+1)-1
            if (values) b(ko) = a(k)
            jb(ko) = ja(k)
            ko=ko+1
 800     continue
         ib(i-i1+2) = ko
 900  continue
      return
c--------end-of-dperm1--------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine sdperm2 (i1,i2,a,ja,ia,b,jb,ib,cperm,rperm,istart,
     *        ipos,job)
      integer i1,i2,job,istart,ja(*),ia(*),jb(*),ib(*),cperm(*),rperm(*) 
      real*4 a(*),b(*)
c----------------------------------------------------------------------- 
c     general submatrix permutation/ extraction routine.
c----------------------------------------------------------------------- 
c     extracts rows rperm(i1), rperm(i1+1), ..., rperm(i2) and does an
c     associated column permutation (using array cperm). The resulting
c     submatrix is constructed in b, jb, ib. For added flexibility, the
c     extracted elements are put in sequence starting from row 'istart' 
c     of B. In addition a pointer ipos to the beginning of arrays b,jb,
c     is also allowed (i.e., nonzero elements are accumulated starting in
c     position ipos of b, jb). In most applications istart and ipos are 
c     equal to one. However, the generality adds substantial flexiblity.
c     EXPLE: (1) to permute msr to msr (excluding diagonals) 
c     call dperm2 (1,n,a,ja,ja,b,jb,jb,rperm,rperm,1,n+2) 
c            (2) To extract rows 1 to 10: define rperm and cperm to be
c     identity permutations (rperm(i)=i, i=1,n) and then
c            call dperm2 (1,10,a,ja,ia,b,jb,ib,rperm,rperm,1,1) 
c            (3) to achieve a symmetric permutation as defined by perm: 
c            call dperm2 (1,10,a,ja,ia,b,jb,ib,perm,perm,1,1) 
c            (4) to get a symmetric permutation of A and append the
c            resulting data structure to A's data structure (useful!) 
c            call dperm2 (1,10,a,ja,ia,a,ja,ia(n+1),perm,perm,1,ia(n+1))
c-----------------------------------------------------------------------
c Y. Saad,Sep. 21 1989 / recoded Jan. 28 1991. 
c-----------------------------------------------------------------------
c on entry: 
c---------- 
c n 	= dimension of the matrix
c i1,i2 = extract rows rperm(i1) to rperm(i2) of A, with i1<i2.
c
c a,ja,
c   ia  = input matrix in CSR format
c cperm = integer array of length n containing the permutation arrays
c	  for the columns: cperm(i) is the destination of column j, 
c         i.e., any column index ja(k) is transformed into cperm(ja(k)) 
c
c rperm	=  permutation array for the rows. rperm(i) = origin (in A) of
c          row i in B. This is the reverse permutation relative to the
c          ones used in routines cperm, dperm,.... 
c          rows rperm(i1), rperm(i1)+1, ... rperm(i2) are 
c          extracted from A and stacked into B, starting in row istart
c          of B. 
c istart= starting row for B where extracted matrix is to be added.
c         this is also only a pointer of the be beginning address for
c         ib , on return. 
c ipos  = beginning position in arrays b and jb where to start copying 
c         elements. Thus, ib(istart) = ipos. 
c
c job   = job indicator. if (job .ne.1) values are not copied (i.e.,
c         only pattern is copied).
c
c on return: 
c-----------
c b,ja,
c ib   = matrix in csr format. positions 1,2,...,istart-1 of ib 
c     are not touched. b(ipos:ipos+nnz-1),jb(ipos:ipos+nnz-1) 
c     contain the value and column indices respectively of the nnz
c     nonzero elements of the permuted matrix. thus ib(istart)=ipos.
c
c Notes:
c------- 
c  1) algorithm is NOT in place 
c  2) column indices may not be sorted on return even  though they 
c     may be on entry.
c----------------------------------------------------------------------- 
c local variables
c
      integer ko,irow,k 
      logical values
c-----------------------------------------------------------------------
      values = (job .eq. 1) 
      ko = ipos 
      ib(istart) = ko
      do 900 i=i1,i2
         irow = rperm(i) 
         do 800 k=ia(irow),ia(irow+1)-1
            if (values) b(ko) = a(k)
            jb(ko) = cperm(ja(k))
            ko=ko+1
 800     continue
         ib(istart+i-i1+1) = ko
 900  continue
      return
c--------end-of-dperm2--------------------------------------------------
c-----------------------------------------------------------------------
      end 
c-----------------------------------------------------------------------
      subroutine sdmperm (nrow,a,ja,ao,jao,perm,job)
      integer nrow,ja(*),jao(*),perm(nrow),job
      real*4 a(*),ao(*) 
c-----------------------------------------------------------------------
c This routine performs a symmetric permutation of the rows and 
c columns of a matrix stored in MSR format. i.e., it computes 
c B = P A transp(P), where P, is  a permutation matrix. 
c P maps row i into row perm(i) and column j into column perm(j): 
c      a(i,j)    becomes   a(perm(i),perm(j)) in new matrix
c (i.e.  ao(perm(i),perm(j)) = a(i,j) ) 
c calls dperm. 
c-----------------------------------------------------------------------
c Y. Saad, Nov 15, 1991. 
c-----------------------------------------------------------------------
c on entry: 
c---------- 
c n 	= dimension of the matrix
c a, ja = input matrix in MSR format. 
c perm 	= integer array of length n containing the permutation arrays
c	  for the rows: perm(i) is the destination of row i in the
c         permuted matrix -- also the destination of column i in case
c         permutation is symmetric (job .le. 2) 
c
c job	= integer indicating the work to be done:
c 		job = 1	permute a, ja, ia into ao, jao, iao 
c 		job = 2 permute matrix ignoring real values.
c		
c on return: 
c-----------
c ao, jao = output matrix in MSR. 
c
c in case job .eq. 2 a and ao are never referred to and can be dummy 
c arguments. 
c
c Notes:
c------- 
c  1) algorithm is NOT in place 
c  2) column indices may not be sorted on return even  though they may be 
c     on entry.
c----------------------------------------------------------------------c
c     local variables
c     
      integer n1, n2
      n1 = nrow+1
      n2 = n1+1
c      
      call sdperm (nrow,a,ja,ja,ao(n2),jao(n2),jao,perm,perm,job) 
c     
      jao(1) = n2
      do 101 j=1, nrow 
         ao(perm(j)) = a(j) 
         jao(j+1) = jao(j+1)+n1
 101  continue
c
c done
c     
      return
c-----------------------------------------------------------------------
c--------end-of-dmperm--------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine sdvperm (n, x, perm) 
      integer n, perm(n) 
      real*4 x(n)
c-----------------------------------------------------------------------
c this subroutine performs an in-place permutation of a real vector x 
c according to the permutation array perm(*), i.e., on return, 
c the vector x satisfies,
c
c	x(perm(j)) :== x(j), j=1,2,.., n
c
c-----------------------------------------------------------------------
c on entry:
c---------
c n 	= length of vector x.
c perm 	= integer array of length n containing the permutation  array.
c x	= input vector
c
c on return:
c---------- 
c x	= vector x permuted according to x(perm(*)) :=  x(*)
c
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c
c local variables 
      real*4 tmp, tmp1
c
      init      = 1
      tmp	= x(init)	
      ii        = perm(init)
      perm(init)= -perm(init)
      k         = 0
c     
c loop
c 
 6    k = k+1
c
c save the chased element --
c 
      tmp1	  = x(ii) 
      x(ii)     = tmp
      next	  = perm(ii) 
      if (next .lt. 0 ) goto 65
c     
c test for end 
c
      if (k .gt. n) goto 101
      tmp       = tmp1
      perm(ii)  = - perm(ii)
      ii        = next 
c
c end loop 
c
      goto 6
c
c reinitilaize cycle --
c
 65   init      = init+1
      if (init .gt. n) goto 101
      if (perm(init) .lt. 0) goto 65
      tmp	= x(init)
      ii	= perm(init)
      perm(init)=-perm(init)
      goto 6
c     
 101  continue
      do 200 j=1, n
         perm(j) = -perm(j)
 200  continue 
c     
      return
c-------------------end-of-dvperm--------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine sretmx (n,a,ja,ia,dd)
      real*4 a(*),dd(*)
      integer n,ia(*),ja(*)
c-----------------------------------------------------------------------
c returns in dd(*) the max absolute value of elements in row *.
c used for scaling purposes. superseded by rnrms  .
c
c on entry:
c n	= dimension of A
c a,ja,ia
c	= matrix stored in compressed sparse row format
c dd	= real*4 array of length n. On output,entry dd(i) contains
c	  the element of row i that has the largest absolute value.
c	  Moreover the sign of dd is modified such that it is the
c	  same as that of the diagonal element in row i.
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c
c local variables
      integer k2, i, k1, k
      real*4 t, t1, t2
c
c initialize 
c
      k2 = 1
      do 11 i=1,n
         k1 = k2
         k2 = ia(i+1) - 1
         t = 0.0d0
         do 101  k=k1,k2
            t1 = abs(a(k))
            if (t1 .gt. t) t = t1
            if (ja(k) .eq. i) then 
               if (a(k) .ge. 0.0) then 
                  t2 = a(k) 
               else 
                  t2 = - a(k)
               endif
            endif
 101     continue		
         dd(i) =  t2*t
c     we do not invert diag
 11   continue
      return
c---------end of retmx -------------------------------------------------
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine sdscaldg (n,a,ja,ia,diag,job)
      real*4 a(*), diag(*),t
      integer ia(*),ja(*)
c----------------------------------------------------------------------- 
c scales rows by diag where diag is either given (job=0)
c or to be computed:
c  job = 1 ,scale row i by  by  +/- max |a(i,j) | and put inverse of 
c       scaling factor in diag(i),where +/- is the sign of a(i,i).
c  job = 2 scale by 2-norm of each row..
c if diag(i) = 0,then diag(i) is replaced by one
c (no scaling)..
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c
      goto (12,11,10) job+1
 10   do 110 j=1,n
         k1= ia(j)
         k2 = ia(j+1)-1
         t = 0.0d0
         do 111 k = k1,k2
 111        t = t+a(k)*a(k)
 110        diag(j) = sqrt(t)
            goto 12
 11   continue
      call sretmx (n,a,ja,ia,diag)
c------
 12   do 1 j=1,n
         if (diag(j) .ne. 0.0d0) then 
            diag(j) = 1.0d0/diag(j)
         else 
            diag(j) = 1.0d0
         endif
 1    continue
      do 2 i=1,n
         t = diag(i)
         do 21 k=ia(i),ia(i+1) -1
            a(k) = a(k)*t
 21      continue
 2    continue
      return 
c--------end of dscaldg -----------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine sextbdg (n,a,ja,ia,bdiag,nblk,ao,jao,iao)
      implicit real*4 (a-h,o-z)
      real*4 bdiag(*),a(*),ao(*)
      integer ia(*),ja(*),jao(*),iao(*) 
c-----------------------------------------------------------------------
c this subroutine extracts the main diagonal blocks of a 
c matrix stored in compressed sparse row format and puts the result
c into the array bdiag and the remainder in ao,jao,iao.
c-----------------------------------------------------------------------
c on entry:
c----------
c n	= integer. The row dimension of the matrix a.
c a,
c ja,
c ia    = matrix stored in csr format
c nblk  = dimension of each diagonal block. The diagonal blocks are
c         stored in compressed format rowwise,i.e.,we store in 
c	  succession the i nonzeros of the i-th row after those of
c	  row number i-1..
c
c on return:
c----------
c bdiag = real*4 array of size (n x nblk) containing the diagonal
c	  blocks of A on return
c ao,
c jao,
C iao   = remainder of the matrix stored in csr format.
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c
      m = 1 + (n-1)/nblk
c this version is sequential -- there is a more parallel version
c that goes through the structure twice ....
      ltr =  ((nblk-1)*nblk)/2 
      l = m * ltr
      do 1 i=1,l
         bdiag(i) = 0.0d0
 1    continue
      ko = 0
      kb = 1
      iao(1) = 1
c-------------------------
      do 11 jj = 1,m
         j1 = (jj-1)*nblk+1
         j2 =  min0 (n,j1+nblk-1)
         do 12 j=j1,j2
            do 13 i=ia(j),ia(j+1) -1
               k = ja(i)
               if (k .lt. j1) then
                  ko = ko+1
                  ao(ko) = a(i)
                  jao(ko) = k
               else if (k .lt. j) then
c     kb = (jj-1)*ltr+((j-j1)*(j-j1-1))/2+k-j1+1
c     bdiag(kb) = a(i)
                  bdiag(kb+k-j1) = a(i)
               endif
 13         continue 	
            kb = kb + j-j1
            iao(j+1) = ko+1
 12      continue
 11   continue
      return
c---------end-of-extbdg------------------------------------------------- 
c----------------------------------------------------------------------- 
      end
c----------------------------------------------------------------------- 
      subroutine sgetbwd (n,a,ja,ia,ml,mu)
c-----------------------------------------------------------------------
c gets the bandwidth of lower part and upper part of A.
c does not assume that A is sorted.
c-----------------------------------------------------------------------
c on entry:
c----------
c n	= integer = the row dimension of the matrix
c a, ja,
c    ia = matrix in compressed sparse row format.
c 
c on return:
c----------- 
c ml	= integer. The bandwidth of the strict lower part of A
c mu	= integer. The bandwidth of the strict upper part of A 
c
c Notes:
c ===== ml and mu are allowed to be negative or return. This may be 
c       useful since it will tell us whether a band is confined 
c       in the strict  upper/lower triangular part. 
c       indeed the definitions of ml and mu are
c
c       ml = max ( (i-j)  s.t. a(i,j) .ne. 0  )
c       mu = max ( (j-i)  s.t. a(i,j) .ne. 0  )
c----------------------------------------------------------------------c
c Y. Saad, Sep. 21 1989                                                c
c----------------------------------------------------------------------c
      real*4 a(*) 
      integer ja(*),ia(n+1),ml,mu,ldist,i,k 
      ml = - n
      mu = - n
      do 3 i=1,n
         do 31 k=ia(i),ia(i+1)-1 
            ldist = i-ja(k)
            ml = max(ml,ldist)
            mu = max(mu,-ldist)
 31      continue
 3    continue
      return
c---------------end-of-getbwd ------------------------------------------ 
c----------------------------------------------------------------------- 
      end
c----------------------------------------------------------------------- 
      subroutine srnrms (nrow, nrm, a, ja, ia, diag) 
      real*4 a(*), diag(nrow), scal 
      integer ja(*), ia(nrow+1) 
c-----------------------------------------------------------------------
c gets the norms of each row of A. (choice of three norms)
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c on return:
c----------
c
c diag = real vector of length nrow containing the norms
c
c-----------------------------------------------------------------
      do 1 ii=1,nrow
c
c     compute the norm if each element.
c     
         scal = 0.0d0
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         if (nrm .eq. 0) then
            do 2 k=k1, k2
               scal = max(scal,abs(a(k) ) ) 
 2          continue
         elseif (nrm .eq. 1) then
            do 3 k=k1, k2
               scal = scal + abs(a(k) ) 
 3          continue
         else
            do 4 k=k1, k2
               scal = scal+a(k)**2
 4          continue
         endif 
         if (nrm .eq. 2) scal = sqrt(scal) 
         diag(ii) = scal
 1    continue
      return
c-----------------------------------------------------------------------
c-------------end-of-rnrms----------------------------------------------
      end 
c----------------------------------------------------------------------- 
      subroutine scnrms (nrow, nrm, a, ja, ia, diag) 
      real*4 a(*), diag(nrow) 
      integer ja(*), ia(nrow+1) 
c-----------------------------------------------------------------------
c gets the norms of each column of A. (choice of three norms)
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c on return:
c----------
c
c diag = real vector of length nrow containing the norms
c
c-----------------------------------------------------------------
      do 10 k=1, nrow 
         diag(k) = 0.0d0
 10   continue
      do 1 ii=1,nrow
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         do 2 k=k1, k2
            j = ja(k) 
c     update the norm of each column
            if (nrm .eq. 0) then
               diag(j) = max(diag(j),abs(a(k) ) ) 
            elseif (nrm .eq. 1) then
               diag(j) = diag(j) + abs(a(k) ) 
            else
               diag(j) = diag(j)+a(k)**2
            endif 
 2       continue
 1    continue
      if (nrm .ne. 2) return
      do 3 k=1, nrow
         diag(k) = sqrt(diag(k))
 3    continue
      return
c-----------------------------------------------------------------------
c------------end-of-cnrms-----------------------------------------------
      end 
c----------------------------------------------------------------------- 
      subroutine sroscal (nrow,job,nrm,a,ja,ia,diag,b,jb,ib,ierr) 
      real*4 a(*), b(*), diag(nrow) 
      integer nrow,job,nrm,ja(*),jb(*),ia(nrow+1),ib(nrow+1),ierr 
c-----------------------------------------------------------------------
c scales the rows of A such that their norms are one on return
c 3 choices of norms: 1-norm, 2-norm, max-norm.
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c on return:
c----------
c
c diag = diagonal matrix stored as a vector containing the matrix
c        by which the rows have been scaled, i.e., on return 
c        we have B = Diag*A.
c
c b, 
c jb, 
c ib	= resulting matrix B in compressed sparse row sparse format.
c	    
c ierr  = error message. ierr=0     : Normal return 
c                        ierr=i > 0 : Row number i is a zero row.
c Notes:
c-------
c 1)        The column dimension of A is not needed. 
c 2)        algorithm in place (B can take the place of A).
c-----------------------------------------------------------------
      call srnrms (nrow,nrm,a,ja,ia,diag)
      ierr = 0
      do 1 j=1, nrow
         if (diag(j) .eq. 0.0d0) then
            ierr = j 
            return
         else
            diag(j) = 1.0d0/diag(j)
         endif
 1    continue
      call sdiamua (nrow,job,a,ja,ia,diag,b,jb,ib)
      return
c-------end-of-roscal---------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine scoscal (nrow,job,nrm,a,ja,ia,diag,b,jb,ib,ierr) 
c----------------------------------------------------------------------- 
      real*4 a(*),b(*),diag(nrow) 
      integer nrow,job,ja(*),jb(*),ia(nrow+1),ib(nrow+1),ierr 
c-----------------------------------------------------------------------
c scales the columns of A such that their norms are one on return
c result matrix written on b, or overwritten on A.
c 3 choices of norms: 1-norm, 2-norm, max-norm. in place.
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c on return:
c----------
c
c diag = diagonal matrix stored as a vector containing the matrix
c        by which the columns have been scaled, i.e., on return 
c        we have B = A * Diag
c
c b, 
c jb, 
c ib	= resulting matrix B in compressed sparse row sparse format.
c
c ierr  = error message. ierr=0     : Normal return 
c                        ierr=i > 0 : Column number i is a zero row.
c Notes:
c-------
c 1)     The column dimension of A is not needed. 
c 2)     algorithm in place (B can take the place of A).
c-----------------------------------------------------------------
      call scnrms (nrow,nrm,a,ja,ia,diag)
      ierr = 0
      do 1 j=1, nrow
         if (diag(j) .eq. 0.0) then
            ierr = j 
            return
         else
            diag(j) = 1.0d0/diag(j)
         endif
 1    continue
      call samudia (nrow,job,a,ja,ia,diag,b,jb,ib)
      return
c--------end-of-coscal-------------------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine saddblk (nrowa, ncola, a, ja, ia, ipos, jpos, job,
     & nrowb, ncolb, b, jb, ib, nrowc, ncolc, c, jc, ic, nzmx, ierr)
c      implicit none
      integer nrowa, nrowb, nrowc, ncola, ncolb, ncolc, ipos, jpos
      integer nzmx, ierr, job
      integer ja(1:*), ia(1:*), jb(1:*), ib(1:*), jc(1:*), ic(1:*)
      real*4 a(1:*), b(1:*), c(1:*)
c-----------------------------------------------------------------------
c     This subroutine adds a matrix B into a submatrix of A whose 
c     (1,1) element is located in the starting position (ipos, jpos). 
c     The resulting matrix is allowed to be larger than A (and B), 
c     and the resulting dimensions nrowc, ncolc will be redefined 
c     accordingly upon return.  
c     The input matrices are assumed to be sorted, i.e. in each row
c     the column indices appear in ascending order in the CSR format.
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrowa    = number of rows in A.
c bcola    = number of columns in A.
c a,ja,ia  = Matrix A in compressed sparse row format with entries sorted
c nrowb    = number of rows in B.
c ncolb    = number of columns in B.
c b,jb,ib  = Matrix B in compressed sparse row format with entries sorted
c
c nzmax	   = integer. The  length of the arrays c and jc. addblk will 
c            stop if the number of nonzero elements in the matrix C
c            exceeds nzmax. See ierr.
c 
c on return:
c----------
c nrowc    = number of rows in C.
c ncolc    = number of columns in C.
c c,jc,ic  = resulting matrix C in compressed sparse row sparse format
c            with entries sorted ascendly in each row. 
c	    
c ierr	   = integer. serving as error message. 
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that addblk stopped while computing the
c         i-th row  of C with i=ierr, because the number 
c         of elements in C exceeds nzmax.
c
c Notes: 
c-------
c     this will not work if any of the two input matrices is not sorted
c-----------------------------------------------------------------------
      logical values
      integer i,j1,j2,ka,kb,kc,kamax,kbmax
      values = (job .ne. 0) 
      ierr = 0
      nrowc = max(nrowa, nrowb+ipos-1)
      ncolc = max(ncola, ncolb+jpos-1)
      kc = 1
      kbmax = 0
      ic(1) = kc
c
      do 10 i=1, nrowc
         if (i.le.nrowa) then
            ka = ia(i)
            kamax = ia(i+1)-1
         else
            ka = ia(nrowa+1)
         end if
         if ((i.ge.ipos).and.((i-ipos).le.nrowb)) then
            kb = ib(i-ipos+1)
            kbmax = ib(i-ipos+2)-1 
         else
            kb = ib(nrowb+1)
         end if
c
c     a do-while type loop -- goes through all the elements in a row.
c
 20      continue 
         if (ka .le. kamax) then
            j1 = ja(ka)
         else
            j1 = ncolc+1
         endif
         if (kb .le. kbmax) then 
            j2 = jb(kb) + jpos - 1
         else 
            j2 = ncolc+1
         endif
c
c     if there are more elements to be added.
c
         if ((ka .le. kamax .or. kb .le. kbmax) .and.
     &        (j1 .le. ncolc .or. j2 .le. ncolc)) then
c
c     three cases
c
            if (j1 .eq. j2) then 
               if (values) c(kc) = a(ka)+b(kb)
               jc(kc) = j1
               ka = ka+1
               kb = kb+1
               kc = kc+1
            else if (j1 .lt. j2) then
               jc(kc) = j1
               if (values) c(kc) = a(ka)
               ka = ka+1
               kc = kc+1
            else if (j1 .gt. j2) then
               jc(kc) = j2
               if (values) c(kc) = b(kb)
               kb = kb+1
               kc = kc+1
            endif
            if (kc .gt. nzmx) goto 999
            goto 20
         end if
         ic(i+1) = kc
 10   continue
      return
 999  ierr = i 
      return
c---------end-of-addblk------------------------------------------------- 
      end
c----------------------------------------------------------------------- 
      subroutine sxtrows (i1,i2,a,ja,ia,ao,jao,iao,iperm,job)
      integer i1,i2,ja(*),ia(*),jao(*),iao(*),iperm(*),job
      real*4 a(*),ao(*) 
c-----------------------------------------------------------------------
c this subroutine extracts given rows from a matrix in CSR format. 
c Specifically, rows number iperm(i1), iperm(i1+1), ...., iperm(i2)
c are extracted and put in the output matrix ao, jao, iao, in CSR
c format.  NOT in place. 
c Youcef Saad -- coded Feb 15, 1992. 
c-----------------------------------------------------------------------
c on entry:
c----------
c i1,i2   = two integers indicating the rows to be extracted.
c           xtrows will extract rows iperm(i1), iperm(i1+1),..,iperm(i2),
c           from original matrix and stack them in output matrix 
c           ao, jao, iao in csr format
c
c a, ja, ia = input matrix in csr format
c
c iperm	= integer array of length nrow containing the reverse permutation 
c         array for the rows. row number iperm(j) in permuted matrix PA
c         used to be row number j in unpermuted matrix.
c         ---> a(i,j) in the permuted matrix was a(iperm(i),j) 
c         in the inout matrix.
c
c job	= integer indicating the work to be done:
c 		job .ne. 1 : get structure only of output matrix,,
c               i.e., ignore real values. (in which case arrays a 
c               and ao are not used nor accessed).
c 		job = 1	get complete data structure of output matrix. 
c               (i.e., including arrays ao and iao).
c------------
c on return: 
c------------ 
c ao, jao, iao = input matrix in a, ja, ia format
c note : 
c        if (job.ne.1)  then the arrays a and ao are not used.
c----------------------------------------------------------------------c
c           Y. Saad, revised May  2, 1990                              c
c----------------------------------------------------------------------c
      logical values
      values = (job .eq. 1) 
c
c copying 
c
      ko = 1
      iao(1) = ko
      do 100 j=i1,i2 
c
c ii=iperm(j) is the index of old row to be copied.
c        
         ii = iperm(j) 
         do 60 k=ia(ii), ia(ii+1)-1 
            jao(ko) = ja(k) 
            if (values) ao(ko) = a(k)
            ko = ko+1
 60      continue
         iao(j-i1+2) = ko
 100  continue
c
      return
c---------end-of-xtrows------------------------------------------------- 
c-----------------------------------------------------------------------
      end

