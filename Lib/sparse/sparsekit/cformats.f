c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c                    FORMAT CONVERSION MODULE                          c
c----------------------------------------------------------------------c
c contents:                                                            c
c----------                                                            c
c csrdns  : converts a row-stored sparse matrix into the dense format. c
c dnscsr  : converts a dense matrix to a sparse storage format.        c
c coocsr  : converts coordinate to  to csr format                      c
c coicsr  : in-place conversion of coordinate to csr format            c
c csrcoo  : converts compressed sparse row to coordinate.              c
c csrssr  : converts compressed sparse row to symmetric sparse row     c
c ssrcsr  : converts symmetric sparse row to compressed sparse row     c
c csrell  : converts compressed sparse row to ellpack format           c
c ellcsr  : converts ellpack format to compressed sparse row format    c
c csrmsr  : converts compressed sparse row format to modified sparse   c
c           row format                                                 c
c msrcsr  : converts modified sparse row format to compressed sparse   c
c           row format.                                                c
c csrcsc  : converts compressed sparse row format to compressed sparse c
c           column format (transposition)                              c
c csrcsc2 : rectangular version of csrcsc                              c
c csrlnk  : converts compressed sparse row to linked list format       c
c lnkcsr  : converts linked list format to compressed sparse row fmt   c
c csrdia  : converts a compressed sparse row format into a diagonal    c
c           format.                                                    c
c diacsr  : converts a diagonal format into a compressed sparse row    c
c           format.                                                    c
c bsrcsr  : converts a block-row sparse format into a compressed       c
c           sparse row format.                                         c
c csrbsr  : converts a compressed sparse row format into a block-row   c
c           sparse format.                                             c
c csrbnd  : converts a compressed sparse row format into a banded      c
c           format (linpack style).                                    c
c bndcsr  : converts a banded format (linpack style) into a compressed c
c           sparse row storage.                                        c
c csrssk  : converts the compressed sparse row format to the symmetric c
c           skyline format                                             c
c sskssr  : converts symmetric skyline format to symmetric  sparse row c
c           format.                                                    c
c csrjad  : converts the csr format into the jagged diagonal format    c
c jadcsr  : converts the jagged-diagonal format into the csr format    c
c csruss  : Compressed Sparse Row to Unsymmetric Sparse Skyline        c
c           format                                                     c
c usscsr  : Unsymmetric Sparse Skyline format to Compressed Sparse Row c
c csrsss  : Compressed Sparse Row to Symmetric Sparse Skyline format   c
c ssscsr  : Symmetric Sparse Skyline format to Compressed Sparse Row   c
c csrvbr  : Converts compressed sparse row to var block row format     c
c vbrcsr  : Converts var block row to compressed sparse row format     c
c csorted : Checks if matrix in CSR format is sorted by columns        c
c--------- miscalleneous additions not involving the csr format--------c
c cooell  : converts coordinate to Ellpack/Itpack format               c
c dcsort  : sorting routine used by crsjad                             c
c----------------------------------------------------------------------c
      subroutine ccsrdns (nrow,ncol,a,ja,ia,dns,ndns,ierr) 
      complex*8 dns(ndns,*),a(*)
      integer ja(*),ia(*)
c-----------------------------------------------------------------------
c Compressed Sparse Row    to    Dense 
c-----------------------------------------------------------------------
c
c converts a row-stored sparse matrix into a densely stored one
c
c On entry:
c---------- 
c
c nrow	= row-dimension of a
c ncol	= column dimension of a
c a, 
c ja, 
c ia    = input matrix in compressed sparse row format. 
c         (a=value array, ja=column array, ia=pointer array)
c dns   = array where to store dense matrix
c ndns	= first dimension of array dns 
c
c on return: 
c----------- 
c dns   = the sparse matrix a, ja, ia has been stored in dns(ndns,*)
c 
c ierr  = integer error indicator. 
c         ierr .eq. 0  means normal return
c         ierr .eq. i  means that the code has stopped when processing
c         row number i, because it found a column number .gt. ncol.
c 
c----------------------------------------------------------------------- 
      ierr = 0
      do 1 i=1, nrow
         do 2 j=1,ncol
	    dns(i,j) = 0.0d0
 2       continue
 1    continue
c     
      do 4 i=1,nrow
         do 3 k=ia(i),ia(i+1)-1
            j = ja(k) 
	    if (j .gt. ncol) then
               ierr = i
               return
	    endif
	    dns(i,j) = a(k)
 3       continue	   
 4    continue
      return
c---- end of csrdns ----------------------------------------------------
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine cdnscsr (nrow,ncol,nzmax,dns,ndns,a,ja,ia,ierr)
      complex*8 dns(ndns,*),a(*)
      integer ia(*),ja(*)
c-----------------------------------------------------------------------
c Dense		to    Compressed Row Sparse 
c----------------------------------------------------------------------- 
c
c converts a densely stored matrix into a row orientied
c compactly sparse matrix. ( reverse of csrdns )
c Note: this routine does not check whether an element 
c is small. It considers that a(i,j) is zero if it is exactly
c equal to zero: see test below.
c-----------------------------------------------------------------------
c on entry:
c---------
c
c nrow	= row-dimension of a
c ncol	= column dimension of a
c nzmax = maximum number of nonzero elements allowed. This
c         should be set to be the lengths of the arrays a and ja.
c dns   = input nrow x ncol (dense) matrix.
c ndns	= first dimension of dns. 
c
c on return:
c---------- 
c 
c a, ja, ia = value, column, pointer  arrays for output matrix 
c
c ierr	= integer error indicator: 
c         ierr .eq. 0 means normal retur
c         ierr .eq. i means that the the code stopped while
c         processing row number i, because there was no space left in
c         a, and ja (as defined by parameter nzmax).
c----------------------------------------------------------------------- 
      ierr = 0
      next = 1
      ia(1) = 1
      do 4 i=1,nrow
         do 3 j=1, ncol 
            if (dns(i,j) .eq. 0.0d0) goto 3
            if (next .gt. nzmax) then
               ierr = i
               return
            endif
            ja(next) = j
            a(next) = dns(i,j)
            next = next+1
 3       continue	   
         ia(i+1) = next
 4    continue
      return
c---- end of dnscsr ---------------------------------------------------- 
c----------------------------------------------------------------------- 
      end
c----------------------------------------------------------------------- 
      subroutine ccoocsr (nrow,nnz,a,ir,jc,ao,jao,iao)
c----------------------------------------------------------------------- 
      complex*8 a(*),ao(*),x
      integer ir(*),jc(*),jao(*),iao(*)
c-----------------------------------------------------------------------
c  Coordinate     to   Compressed Sparse Row 
c----------------------------------------------------------------------- 
c converts a matrix that is stored in coordinate format
c  a, ir, jc into a row general sparse ao, jao, iao format.
c
c on entry:
c--------- 
c nrow	= dimension of the matrix 
c nnz	= number of nonzero elements in matrix
c a,
c ir, 
c jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
c         nonzero elements of the matrix with a(k) = actual real value of
c 	  the elements, ir(k) = its row number and jc(k) = its column 
c	  number. The order of the elements is arbitrary. 
c
c on return:
c----------- 
c ir 	is destroyed
c
c ao, jao, iao = matrix in general sparse matrix format with ao 
c 	continung the real values, jao containing the column indices, 
c	and iao being the pointer to the beginning of the row, 
c	in arrays ao, jao.
c
c Notes:
c------ This routine is NOT in place.  See coicsr
c
c------------------------------------------------------------------------
      do 1 k=1,nrow+1
         iao(k) = 0
 1    continue
c determine row-lengths.
      do 2 k=1, nnz
         iao(ir(k)) = iao(ir(k))+1
 2    continue
c starting position of each row..
      k = 1
      do 3 j=1,nrow+1
         k0 = iao(j)
         iao(j) = k
         k = k+k0
 3    continue
c go through the structure  once more. Fill in output matrix.
      do 4 k=1, nnz
         i = ir(k)
         j = jc(k)
         x = a(k)
         iad = iao(i)
         ao(iad) =  x
         jao(iad) = j
         iao(i) = iad+1
 4    continue
c shift back iao
      do 5 j=nrow,1,-1
         iao(j+1) = iao(j)
 5    continue
      iao(1) = 1
      return
c------------- end of coocsr ------------------------------------------- 
c----------------------------------------------------------------------- 
      end
c----------------------------------------------------------------------- 
      subroutine ccoicsr (n,nnz,job,a,ja,ia,iwk)
      integer ia(nnz),ja(nnz),iwk(n+1) 
      complex*8 a(*)
c------------------------------------------------------------------------
c IN-PLACE coo-csr conversion routine.
c------------------------------------------------------------------------
c this subroutine converts a matrix stored in coordinate format into 
c the csr format. The conversion is done in place in that the arrays 
c a,ja,ia of the result are overwritten onto the original arrays.
c------------------------------------------------------------------------
c on entry:
c--------- 
c n	= integer. row dimension of A.
c nnz	= integer. number of nonzero elements in A.
c job   = integer. Job indicator. when job=1, the real values in a are
c         filled. Otherwise a is not touched and the structure of the
c         array only (i.e. ja, ia)  is obtained.
c a	= real array of size nnz (number of nonzero elements in A)
c         containing the nonzero elements 
c ja	= integer array of length nnz containing the column positions
c 	  of the corresponding elements in a.
c ia	= integer array of length nnz containing the row positions
c 	  of the corresponding elements in a.
c iwk	= integer work array of length n+1 
c on return:
c----------
c a
c ja 
c ia	= contains the compressed sparse row data structure for the 
c         resulting matrix.
c Note: 
c-------
c         the entries of the output matrix are not sorted (the column
c         indices in each are not in increasing order) use coocsr
c         if you want them sorted.
c----------------------------------------------------------------------c
c  Coded by Y. Saad, Sep. 26 1989                                      c
c----------------------------------------------------------------------c
      complex*8 t,tnext
      logical values
c----------------------------------------------------------------------- 
      values = (job .eq. 1) 
c find pointer array for resulting matrix. 
      do 35 i=1,n+1
         iwk(i) = 0
 35   continue
      do 4 k=1,nnz
         i = ia(k)
         iwk(i+1) = iwk(i+1)+1
 4    continue 
c------------------------------------------------------------------------
      iwk(1) = 1 
      do 44 i=2,n
         iwk(i) = iwk(i-1) + iwk(i)
 44   continue 
c
c     loop for a cycle in chasing process. 
c
      init = 1
      k = 0
 5    if (values) t = a(init)
      i = ia(init)
      j = ja(init)
      ia(init) = -1
c------------------------------------------------------------------------
 6    k = k+1 		
c     current row number is i.  determine  where to go. 
      ipos = iwk(i)
c     save the chased element. 
      if (values) tnext = a(ipos)
      inext = ia(ipos)
      jnext = ja(ipos)
c     then occupy its location.
      if (values) a(ipos)  = t
      ja(ipos) = j
c     update pointer information for next element to come in row i. 
      iwk(i) = ipos+1
c     determine  next element to be chased,
      if (ia(ipos) .lt. 0) goto 65
      t = tnext
      i = inext
      j = jnext 
      ia(ipos) = -1
      if (k .lt. nnz) goto 6
      goto 70
 65   init = init+1
      if (init .gt. nnz) goto 70
      if (ia(init) .lt. 0) goto 65
c     restart chasing --	
      goto 5
 70   do 80 i=1,n 
         ia(i+1) = iwk(i)
 80   continue
      ia(1) = 1
      return
c----------------- end of coicsr ----------------------------------------
c------------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine ccsrcoo (nrow,job,nzmax,a,ja,ia,nnz,ao,ir,jc,ierr)
c-----------------------------------------------------------------------
      complex*8 a(*),ao(*) 
      integer ir(*),jc(*),ja(*),ia(nrow+1) 
c----------------------------------------------------------------------- 
c  Compressed Sparse Row      to      Coordinate 
c----------------------------------------------------------------------- 
c converts a matrix that is stored in coordinate format
c  a, ir, jc into a row general sparse ao, jao, iao format.
c
c on entry: 
c---------
c nrow	= dimension of the matrix.
c job   = integer serving as a job indicator. 
c         if job = 1 fill in only the array ir, ignore jc, and ao.
c         if job = 2 fill in ir, and jc but not ao 
c         if job = 3 fill in everything.
c         The reason why these options are provided is that on return 
c         ao and jc are the same as a, ja. So when job = 3, a and ja are
c         simply copied into ao, jc.  When job=2, only jc and ir are
c         returned. With job=1 only the array ir is returned. Moreover,
c         the algorithm is in place:
c	     call csrcoo (nrow,1,nzmax,a,ja,ia,nnz,a,ia,ja,ierr) 
c         will write the output matrix in coordinate format on a, ja,ia.
c
c a,
c ja,
c ia    = matrix in compressed sparse row format.
c nzmax = length of space available in ao, ir, jc.
c         the code will stop immediatly if the number of
c         nonzero elements found in input matrix exceeds nzmax.
c 
c on return:
c----------- 
c ao, ir, jc = matrix in coordinate format.
c
c nnz        = number of nonzero elements in matrix.
c ierr       = integer error indicator.
c         ierr .eq. 0 means normal retur
c         ierr .eq. 1 means that the the code stopped 
c         because there was no space in ao, ir, jc 
c         (according to the value of  nzmax).
c 
c NOTES: 1)This routine is PARTIALLY in place: csrcoo can be called with 
c         ao being the same array as as a, and jc the same array as ja. 
c         but ir CANNOT be the same as ia. 
c         2) note the order in the output arrays, 
c------------------------------------------------------------------------
      ierr = 0
      nnz = ia(nrow+1)-1
      if (nnz .gt. nzmax) then
         ierr = 1
         return
      endif
c------------------------------------------------------------------------
      goto (3,2,1) job
 1    do 10 k=1,nnz
         ao(k) = a(k)
 10   continue
 2    do 11 k=1,nnz
         jc(k) = ja(k)
 11   continue
c
c     copy backward to allow for in-place processing. 
c
 3    do 13 i=nrow,1,-1
         k1 = ia(i+1)-1
         k2 = ia(i)
         do 12 k=k1,k2,-1
            ir(k) = i
 12      continue
 13   continue
      return
c------------- end-of-csrcoo ------------------------------------------- 
c----------------------------------------------------------------------- 
      end
c----------------------------------------------------------------------- 
      subroutine ccsrssr (nrow,a,ja,ia,nzmax,ao,jao,iao,ierr)
      complex*8 a(*), ao(*), t
      integer ia(*), ja(*), iao(*), jao(*)
c-----------------------------------------------------------------------
c Compressed Sparse Row     to     Symmetric Sparse Row
c----------------------------------------------------------------------- 
c this subroutine extracts the lower triangular part of a matrix.
c It can used as a means for converting a symmetric matrix for 
c which all the entries are stored in sparse format into one
c in which only the lower part is stored. The routine is in place in 
c that the output matrix ao, jao, iao can be overwritten on 
c the  input matrix  a, ja, ia if desired. Csrssr has been coded to
c put the diagonal elements of the matrix in the last position in
c each row (i.e. in position  ao(ia(i+1)-1   of ao and jao) 
c----------------------------------------------------------------------- 
c On entry
c-----------
c nrow  = dimension of the matrix a.
c a, ja, 
c    ia = matrix stored in compressed row sparse format
c
c nzmax = length of arrays ao,  and jao. 
c
c On return:
c----------- 
c ao, jao, 
c     iao = lower part of input matrix (a,ja,ia) stored in compressed sparse 
c          row format format.
c  
c ierr   = integer error indicator. 
c          ierr .eq. 0  means normal return
c          ierr .eq. i  means that the code has stopped when processing
c          row number i, because there is not enough space in ao, jao
c          (according to the value of nzmax) 
c
c----------------------------------------------------------------------- 
      ierr = 0
      ko = 0
c-----------------------------------------------------------------------
      do  7 i=1, nrow
         kold = ko
         kdiag = 0
         do 71 k = ia(i), ia(i+1) -1
            if (ja(k)  .gt. i) goto 71
            ko = ko+1
            if (ko .gt. nzmax) then
               ierr = i
               return
            endif
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
      iao(nrow+1) = ko+1
      return
c--------- end of csrssr ----------------------------------------------- 
c----------------------------------------------------------------------- 
      end
c----------------------------------------------------------------------- 
      subroutine cssrcsr (job, value2, nrow, a, ja, ia, nzmax,
     &                  ao, jao, iao, indu, iwk, ierr)
c     .. Scalar Arguments ..
      integer            ierr, job, nrow, nzmax, value2
c     ..
c     .. Array Arguments ..
      integer            ia(nrow+1), iao(nrow+1), indu(nrow),
     &                   iwk(nrow+1), ja(*), jao(nzmax)
      complex*8             a(*), ao(nzmax)
c     ..
c-----------------------------------------------------------------------
c     Symmetric Sparse Row to Compressed Sparse Row format
c-----------------------------------------------------------------------
c     This subroutine converts a given matrix in SSR format to regular
c     CSR format by computing Ao = A + A' - diag(A), where A' is A
c     transpose.
c
c     Typically this routine is used to expand the SSR matrix of
c     Harwell Boeing matrices, or to obtain a symmetrized graph of
c     unsymmetric matrices.
c
c     This routine is inplace, i.e., (Ao,jao,iao) may be same as
c     (a,ja,ia).
c
c     It is possible to input an arbitrary CSR matrix to this routine,
c     since there is no syntactical difference between CSR and SSR
c     format. It also removes duplicate entries and perform a partial
c     ordering. The output matrix has an order of lower half, main
c     diagonal and upper half after the partial ordering.
c-----------------------------------------------------------------------
c on entry:
c---------
c
c job   = options
c         0 -- duplicate entries are not removed. If the input matrix is
c              SSR (not an arbitary CSR) matrix, no duplicate entry should
c              arise from this routine.
c         1 -- eliminate duplicate entries, zero entries.
c         2 -- eliminate duplicate entries and perform partial ordering.
c         3 -- eliminate duplicate entries, sort the entries in the
c              increasing order of clumn indices.
c              
c value2= will the values of A be copied?
c         0 -- only expand the graph (a, ao are not touched)
c         1 -- expand the matrix with the values.
c
c nrow  = column dimension of inout matrix
c a,
c ia,
c ja    = matrix in compressed sparse row format.
c
c nzmax = size of arrays ao and jao. SSRCSR will abort if the storage
c          provided in ao, jao is not sufficient to store A. See ierr.
c
c on return:
c----------
c ao, jao, iao
c       = output matrix in compressed sparse row format. The resulting
c         matrix is symmetric and is equal to A+A'-D. ao, jao, iao,
c         can be the same as a, ja, ia in the calling sequence.
c
c indu  = integer array of length nrow. INDU will contain pointers
c         to the beginning of upper traigular part if job > 1.
c         Otherwise it is also used as a work array (size nrow).
c
c iwk   = integer work space (size nrow+1).
c
c ierr  = integer. Serving as error message. If the length of the arrays
c         ao, jao exceeds nzmax, ierr returns the minimum value
c         needed for nzmax. otherwise ierr=0 (normal return).
c
c-----------------------------------------------------------------------
c     .. Local Scalars ..
      integer            i, ipos, j, k, kfirst, klast, ko, kosav, nnz
      complex*8             tmp
c     ..
c     .. Executable Statements ..
      ierr = 0
      do 10 i = 1, nrow
         indu(i) = 0
         iwk(i) = 0
 10   continue
      iwk(nrow+1) = 0
c
c     .. compute number of elements in each row of (A'-D)
c     put result in iwk(i+1)  for row i.
c
      do 30 i = 1, nrow
         do 20 k = ia(i), ia(i+1) - 1
            j = ja(k)
            if (j.ne.i)
     &         iwk(j+1) = iwk(j+1) + 1
 20      continue
 30   continue
c
c     .. find addresses of first elements of ouput matrix. result in iwk
c
      iwk(1) = 1
      do 40 i = 1, nrow
         indu(i) = iwk(i) + ia(i+1) - ia(i)
         iwk(i+1) = iwk(i+1) + indu(i)
         indu(i) = indu(i) - 1
 40   continue
c.....Have we been given enough storage in ao, jao ?
      nnz = iwk(nrow+1) - 1
      if (nnz.gt.nzmax) then
         ierr = nnz
         return
      endif
c
c     .. copy the existing matrix (backwards).
c
      kosav = iwk(nrow+1)
      do 60 i = nrow, 1, -1
         klast = ia(i+1) - 1
         kfirst = ia(i)
         iao(i+1) = kosav
         kosav = iwk(i)
         ko = iwk(i) - kfirst
         iwk(i) = ko + klast + 1
         do 50 k = klast, kfirst, -1
            if (value2.ne.0)
     &         ao(k+ko) = a(k)
            jao(k+ko) = ja(k)
 50      continue
 60   continue
      iao(1) = 1
c
c     now copy (A'-D). Go through the structure of ao, jao, iao
c     that has already been copied. iwk(i) is the address
c     of the next free location in row i for ao, jao.
c
      do 80 i = 1, nrow
         do 70 k = iao(i), indu(i)
            j = jao(k)
            if (j.ne.i) then
               ipos = iwk(j)
               if (value2.ne.0)
     &            ao(ipos) = ao(k)
               jao(ipos) = i
               iwk(j) = ipos + 1
            endif
 70      continue
 80   continue
      if (job.le.0) return
c
c     .. eliminate duplicate entries --
c     array INDU is used as marker for existing indices, it is also the
c     location of the entry.
c     IWK is used to stored the old IAO array.
c     matrix is copied to squeeze out the space taken by the duplicated
c     entries.
c
      do 90 i = 1, nrow
         indu(i) = 0
         iwk(i) = iao(i)
 90   continue
      iwk(nrow+1) = iao(nrow+1)
      k = 1
      do 120 i = 1, nrow
         iao(i) = k
         ipos = iwk(i)
         klast = iwk(i+1)
 100     if (ipos.lt.klast) then
            j = jao(ipos)
            if (indu(j).eq.0) then
c     .. new entry ..
               if (value2.ne.0) then
                  if (ao(ipos) .ne. 0.0D0) then
                     indu(j) = k
                     jao(k) = jao(ipos)
                     ao(k) = ao(ipos)
                     k = k + 1
                  endif
               else
                  indu(j) = k
                  jao(k) = jao(ipos)
                  k = k + 1
               endif
            else if (value2.ne.0) then
c     .. duplicate entry ..
               ao(indu(j)) = ao(indu(j)) + ao(ipos)
            endif
            ipos = ipos + 1
            go to 100
         endif
c     .. remove marks before working on the next row ..
         do 110 ipos = iao(i), k - 1
            indu(jao(ipos)) = 0
 110     continue
 120  continue
      iao(nrow+1) = k
      if (job.le.1) return
c
c     .. partial ordering ..
c     split the matrix into strict upper/lower triangular
c     parts, INDU points to the the beginning of the strict upper part.
c
      do 140 i = 1, nrow
         klast = iao(i+1) - 1
         kfirst = iao(i)
 130     if (klast.gt.kfirst) then
            if (jao(klast).lt.i .and. jao(kfirst).ge.i) then
c     .. swap klast with kfirst ..
               j = jao(klast)
               jao(klast) = jao(kfirst)
               jao(kfirst) = j
               if (value2.ne.0) then
                  tmp = ao(klast)
                  ao(klast) = ao(kfirst)
                  ao(kfirst) = tmp
               endif
            endif
            if (jao(klast).ge.i)
     &         klast = klast - 1
            if (jao(kfirst).lt.i)
     &         kfirst = kfirst + 1
            go to 130
         endif
c
         if (jao(klast).lt.i) then
            indu(i) = klast + 1
         else
            indu(i) = klast
         endif
 140  continue
      if (job.le.2) return
c
c     .. order the entries according to column indices
c     bubble-sort is used
c
      do 190 i = 1, nrow
         do 160 ipos = iao(i), indu(i)-1
            do 150 j = indu(i)-1, ipos+1, -1
               k = j - 1
               if (jao(k).gt.jao(j)) then
                  ko = jao(k)
                  jao(k) = jao(j)
                  jao(j) = ko
                  if (value2.ne.0) then
                     tmp = ao(k)
                     ao(k) = ao(j)
                     ao(j) = tmp
                  endif
               endif
 150        continue
 160     continue
         do 180 ipos = indu(i), iao(i+1)-1
            do 170 j = iao(i+1)-1, ipos+1, -1
               k = j - 1
               if (jao(k).gt.jao(j)) then
                  ko = jao(k)
                  jao(k) = jao(j)
                  jao(j) = ko
                  if (value2.ne.0) then
                     tmp = ao(k)
                     ao(k) = ao(j)
                     ao(j) = tmp
                  endif
               endif
 170        continue
 180     continue
 190  continue
c
      return
c---- end of ssrcsr ----------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine cxssrcsr (nrow,a,ja,ia,nzmax,ao,jao,iao,indu,ierr)
      integer ia(nrow+1),iao(nrow+1),ja(*),jao(nzmax),indu(nrow+1)
      complex*8 a(*),ao(nzmax)
c-----------------------------------------------------------------------
c Symmetric Sparse Row   to    (regular) Compressed Sparse Row
c----------------------------------------------------------------------- 
c this subroutine converts  a symmetric  matrix in which only the lower 
c part is  stored in compressed sparse row format, i.e.,
c a matrix stored in symmetric sparse format, into a fully stored matrix
c i.e., a matrix where both the lower and upper parts are stored in 
c compressed sparse row format. the algorithm is in place (i.e. result 
c may be overwritten onto the input matrix a, ja, ia ----- ). 
c the output matrix delivered by ssrcsr is such that each row starts with
c the elements of the lower part followed by those of the upper part.
c----------------------------------------------------------------------- 
c on entry:
c--------- 
c	
c nrow  = row dimension of inout matrix
c a, 
c ia, 
c ja    = matrix in compressed sparse row format. This is assumed to be
c         a lower triangular matrix. 
c
c nzmax	= size of arrays ao and jao. ssrcsr will abort if the storage 
c	   provided in a, ja is not sufficient to store A. See ierr. 
c	
c on return:
c----------
c ao, iao, 
c   jao = output matrix in compressed sparse row format. The resulting 
c         matrix is symmetric and is equal to A+A**T - D, if
c         A is the original lower triangular matrix. ao, jao, iao,
c         can be the same as a, ja, ia in the calling sequence.
c      
c indu  = integer array of length nrow+1. If the input matrix is such 
c         that the last element in each row is its diagonal element then
c         on return, indu will contain the pointers to the diagonal 
c         element in each row of the output matrix. Otherwise used as
c         work array.
c ierr  = integer. Serving as error message. If the length of the arrays
c         ao, jao exceeds nzmax, ierr returns the minimum value
c         needed for nzmax. otherwise ierr=0 (normal return).
c 
c----------------------------------------------------------------------- 
      ierr = 0
      do 1 i=1,nrow+1
         indu(i) = 0     
 1    continue
c     
c     compute  number of elements in each row of strict upper part. 
c     put result in indu(i+1)  for row i. 
c     
      do 3 i=1, nrow
         do 2 k=ia(i),ia(i+1)-1 
            j = ja(k)
            if (j .lt. i) indu(j+1) = indu(j+1)+1
 2       continue 
 3    continue
c-----------
c     find addresses of first elements of ouput matrix. result in indu
c-----------
      indu(1) = 1 
      do 4 i=1,nrow
         lenrow = ia(i+1)-ia(i)
         indu(i+1) = indu(i) + indu(i+1) + lenrow
 4    continue
c--------------------- enough storage in a, ja ? --------
      nnz = indu(nrow+1)-1 
      if (nnz .gt. nzmax) then
         ierr = nnz
         return
      endif
c
c     now copy lower part (backwards).
c     
      kosav = indu(nrow+1)
      do 6 i=nrow,1,-1
         klast = ia(i+1)-1
         kfirst = ia(i)
         iao(i+1) = kosav
         ko = indu(i) 
         kosav = ko
         do 5 k = kfirst, klast
            ao(ko) = a(k)
            jao(ko) = ja(k)
	    ko = ko+1
 5       continue
         indu(i) = ko 
 6    continue
      iao(1) = 1
c
c     now copy upper part. Go through the structure of ao, jao, iao
c     that has already been copied (lower part). indu(i) is the address
c     of the next free location in row i for ao, jao.
c     
      do 8 i=1,nrow
c     i-th row is now in ao, jao, iao structure -- lower half part
         do 9 k=iao(i), iao(i+1)-1 
            j = jao(k)
            if (j .ge. i)  goto 8
            ipos = indu(j)
            ao(ipos) = ao(k)
            jao(ipos) = i
            indu(j) = indu(j) + 1 
 9       continue
 8    continue
      return
c----- end of xssrcsr -------------------------------------------------- 
c----------------------------------------------------------------------- 
      end
c-----------------------------------------------------------------------
      subroutine ccsrell (nrow,a,ja,ia,maxcol,coef,jcoef,ncoef,
     *                   ndiag,ierr)
      integer ia(nrow+1), ja(*), jcoef(ncoef,1)  
      complex*8 a(*), coef(ncoef,1)
c----------------------------------------------------------------------- 
c Compressed Sparse Row	    to    Ellpack - Itpack format 
c----------------------------------------------------------------------- 
c this subroutine converts  matrix stored in the general a, ja, ia 
c format into the coef, jcoef itpack format.
c
c----------------------------------------------------------------------- 
c on entry:
c---------- 
c nrow 	  = row dimension of the matrix A.
c
c a, 
c ia, 
c ja      = input matrix in compressed sparse row format. 
c
c ncoef  = first dimension of arrays coef, and jcoef.
c 
c maxcol = integer equal to the number of columns available in coef.
c
c on return: 
c----------
c coef	= real array containing the values of the matrix A in 
c         itpack-ellpack format.
c jcoef = integer array containing the column indices of coef(i,j) 
c         in A.
c ndiag = number of active 'diagonals' found. 
c
c ierr 	= error message. 0 = correct return. If ierr .ne. 0 on
c	  return this means that the number of diagonals found
c         (ndiag) exceeds maxcol.
c
c----------------------------------------------------------------------- 
c first determine the length of each row of lower-part-of(A)
      ierr = 0
      ndiag = 0
      do 3 i=1, nrow
         k = ia(i+1)-ia(i)
         ndiag = max0(ndiag,k) 
 3    continue
c----- check whether sufficient columns are available. ----------------- 
      if (ndiag .gt. maxcol) then
         ierr = 1 
         return
      endif
c
c fill coef with zero elements and jcoef with row numbers.------------ 
c
      do 4 j=1,ndiag 
         do 41 i=1,nrow
            coef(i,j) = 0.0d0
            jcoef(i,j) = i
 41      continue
 4    continue
c     
c------- copy elements row by row.-------------------------------------- 
c     
      do 6 i=1, nrow
         k1 = ia(i)
         k2 = ia(i+1)-1
         do 5 k=k1,k2
            coef(i,k-k1+1) = a(k)
            jcoef(i,k-k1+1) = ja(k)
 5       continue
 6    continue
      return
c--- end of csrell------------------------------------------------------ 
c----------------------------------------------------------------------- 
      end
c-----------------------------------------------------------------------
      subroutine cellcsr (nrow,coef,jcoef,ncoef,ndiag,a,ja,ia,nzmax,
     *                    ierr)
      integer ia(nrow+1), ja(*), jcoef(ncoef,1) 
      complex*8 a(*), coef(ncoef,1)
c----------------------------------------------------------------------- 
c  Ellpack - Itpack format  to  Compressed Sparse Row
c----------------------------------------------------------------------- 
c this subroutine converts a matrix stored in ellpack-itpack format 
c coef-jcoef into the compressed sparse row format. It actually checks
c whether an entry in the input matrix is a nonzero element before
c putting it in the output matrix. The test does not account for small
c values but only for exact zeros. 
c----------------------------------------------------------------------- 
c on entry:
c---------- 
c
c nrow 	= row dimension of the matrix A.
c coef	= array containing the values of the matrix A in ellpack format.
c jcoef = integer arraycontains the column indices of coef(i,j) in A.
c ncoef = first dimension of arrays coef, and jcoef.
c ndiag = number of active columns in coef, jcoef.
c 
c ndiag = on entry the number of columns made available in coef.
c
c on return: 
c----------
c a, ia, 
c    ja = matrix in a, ia, ja format where. 
c 
c nzmax	= size of arrays a and ja. ellcsr will abort if the storage 
c	   provided in a, ja is not sufficient to store A. See ierr. 
c
c ierr 	= integer. serves are output error message. 
c         ierr = 0 means normal return. 
c         ierr = 1 means that there is not enough space in
c         a and ja to store output matrix.
c----------------------------------------------------------------------- 
c first determine the length of each row of lower-part-of(A)
      ierr = 0
c-----check whether sufficient columns are available. ----------------- 
c
c------- copy elements row by row.-------------------------------------- 
      kpos = 1
      ia(1) = kpos
      do 6 i=1, nrow
         do 5 k=1,ndiag
            if (coef(i,k) .ne. 0.0d0) then
               if (kpos .gt. nzmax) then
                  ierr = kpos
                  return
               endif
               a(kpos) = coef(i,k)
               ja(kpos) = jcoef(i,k)
               kpos = kpos+1
	    endif
 5       continue
         ia(i+1) = kpos
 6    continue	
      return
c--- end of ellcsr ----------------------------------------------------- 
c----------------------------------------------------------------------- 
      end
c-----------------------------------------------------------------------
      subroutine ccsrmsr (n,a,ja,ia,ao,jao,wk,iwk)
      complex*8 a(*),ao(*),wk(n)
      integer ia(n+1),ja(*),jao(*),iwk(n+1)
c----------------------------------------------------------------------- 
c Compressed Sparse Row   to      Modified - Sparse Row 
c                                 Sparse row with separate main diagonal
c----------------------------------------------------------------------- 
c converts a general sparse matrix a, ja, ia into 
c a compressed matrix using a separated diagonal (referred to as
c the bell-labs format as it is used by bell labs semi conductor
c group. We refer to it here as the modified sparse row format.
c Note: this has been coded in such a way that one can overwrite
c the output matrix onto the input matrix if desired by a call of
c the form 
c
c     call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
c
c In case ao, jao, are different from a, ja, then one can
c use ao, jao as the work arrays in the calling sequence:
c
c     call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
c
c----------------------------------------------------------------------- 
c
c on entry :
c---------
c a, ja, ia = matrix in csr format. note that the 
c	     algorithm is in place: ao, jao can be the same
c            as a, ja, in which case it will be overwritten on it
c            upon return.
c	 
c on return :
c-----------
c
c ao, jao  = sparse matrix in modified sparse row storage format:
c	   +  ao(1:n) contains the diagonal of the matrix. 
c	   +  ao(n+2:nnz) contains the nondiagonal elements of the
c             matrix, stored rowwise.
c	   +  jao(n+2:nnz) : their column indices
c	   +  jao(1:n+1) contains the pointer array for the nondiagonal
c             elements in ao(n+1:nnz) and jao(n+2:nnz).
c             i.e., for i .le. n+1 jao(i) points to beginning of row i 
c	      in arrays ao, jao.
c	       here nnz = number of nonzero elements+1 
c work arrays:
c------------
c wk	= real work array of length n
c iwk   = integer work array of length n+1
c
c notes: 
c------- 
c        Algorithm is in place.  i.e. both:
c
c          call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
c          (in which  ao, jao, are different from a, ja)
c           and
c          call csrmsr (n, a, ja, ia, a, ja, wk,iwk) 
c          (in which  wk, jwk, are different from a, ja)
c        are OK.
c--------
c coded by Y. Saad Sep. 1989. Rechecked Feb 27, 1990.
c-----------------------------------------------------------------------
      icount = 0
c
c store away diagonal elements and count nonzero diagonal elements.
c
      do 1 i=1,n
         wk(i) = 0.0d0
         iwk(i+1) = ia(i+1)-ia(i)
         do 2 k=ia(i),ia(i+1)-1
            if (ja(k) .eq. i) then
               wk(i) = a(k)
               icount = icount + 1 
               iwk(i+1) = iwk(i+1)-1
            endif
 2       continue
 1    continue
c     
c compute total length
c     
      iptr = n + ia(n+1) - icount
c     
c     copy backwards (to avoid collisions)
c     
      do 500 ii=n,1,-1
         do 100 k=ia(ii+1)-1,ia(ii),-1
            j = ja(k)
            if (j .ne. ii) then
               ao(iptr) = a(k)
               jao(iptr) = j 
               iptr = iptr-1
            endif
 100     continue
 500  continue
c
c compute pointer values and copy wk(*)
c
      jao(1) = n+2
      do 600 i=1,n
         ao(i) = wk(i) 
         jao(i+1) = jao(i)+iwk(i+1)
 600  continue
      return	
c------------ end of subroutine csrmsr ---------------------------------
c----------------------------------------------------------------------- 
      end
c-----------------------------------------------------------------------
      subroutine cmsrcsr (n,a,ja,ao,jao,iao,wk,iwk)
      complex*8 a(*),ao(*),wk(n)
      integer ja(*),jao(*),iao(n+1),iwk(n+1)
c----------------------------------------------------------------------- 
c       Modified - Sparse Row  to   Compressed Sparse Row   
c
c----------------------------------------------------------------------- 
c converts a compressed matrix using a separated diagonal 
c (modified sparse row format) in the Compressed Sparse Row   
c format.
c does not check for zero elements in the diagonal.
c
c
c on entry :
c---------
c n          = row dimension of matrix
c a, ja      = sparse matrix in msr sparse storage format
c              see routine csrmsr for details on data structure 
c        
c on return :
c-----------
c
c ao,jao,iao = output matrix in csr format.  
c
c work arrays:
c------------
c wk       = real work array of length n
c iwk      = integer work array of length n+1
c
c notes:
c   The original version of this was NOT in place, but has
c   been modified by adding the vector iwk to be in place.
c   The original version had ja instead of iwk everywhere in
c   loop 500.  Modified  Sun 29 May 1994 by R. Bramley (Indiana).
c   
c----------------------------------------------------------------------- 
      logical added
      do 1 i=1,n
         wk(i) = a(i)
         iwk(i) = ja(i)
 1    continue
      iwk(n+1) = ja(n+1)
      iao(1) = 1
      iptr = 1
c---------
      do 500 ii=1,n 
         added = .false.
         idiag = iptr + (iwk(ii+1)-iwk(ii)) 
         do 100 k=iwk(ii),iwk(ii+1)-1
            j = ja(k)
            if (j .lt. ii) then
               ao(iptr) = a(k)
               jao(iptr) = j 
               iptr = iptr+1
            elseif (added) then
               ao(iptr) = a(k)
               jao(iptr) = j 
               iptr = iptr+1
            else 
c add diag element - only reserve a position for it. 
               idiag = iptr
               iptr = iptr+1
               added = .true.
c     then other element
               ao(iptr) = a(k)
               jao(iptr) = j 
               iptr = iptr+1
            endif
 100     continue
         ao(idiag) = wk(ii)
         jao(idiag) = ii
         if (.not. added) iptr = iptr+1
         iao(ii+1) = iptr 
 500  continue
      return    
c------------ end of subroutine msrcsr --------------------------------- 
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine ccsrcsc (n,job,ipos,a,ja,ia,ao,jao,iao)
      integer ia(n+1),iao(n+1),ja(*),jao(*)
      complex*8  a(*),ao(*)
c-----------------------------------------------------------------------
c Compressed Sparse Row     to      Compressed Sparse Column
c
c (transposition operation)   Not in place. 
c----------------------------------------------------------------------- 
c -- not in place --
c this subroutine transposes a matrix stored in a, ja, ia format.
c ---------------
c on entry:
c----------
c n	= dimension of A.
c job	= integer to indicate whether to fill the values (job.eq.1) of the
c         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
c
c ipos  = starting position in ao, jao of the transposed matrix.
c         the iao array takes this into account (thus iao(1) is set to ipos.)
c         Note: this may be useful if one needs to append the data structure
c         of the transpose to that of A. In this case use for example
c                call csrcsc (n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2)) 
c	  for any other normal usage, enter ipos=1.
c a	= real array of length nnz (nnz=number of nonzero elements in input 
c         matrix) containing the nonzero elements.
c ja	= integer array of length nnz containing the column positions
c 	  of the corresponding elements in a.
c ia	= integer of size n+1. ia(k) contains the position in a, ja of
c	  the beginning of the k-th row.
c
c on return:
c ---------- 
c output arguments:
c ao	= real array of size nzz containing the "a" part of the transpose
c jao	= integer array of size nnz containing the column indices.
c iao	= integer array of size n+1 containing the "ia" index array of
c	  the transpose. 
c
c----------------------------------------------------------------------- 
      call ccsrcsc2 (n,n,job,ipos,a,ja,ia,ao,jao,iao)
      end
c-----------------------------------------------------------------------
      subroutine ccsrcsc2 (n,n2,job,ipos,a,ja,ia,ao,jao,iao)
      integer ia(n+1),iao(n2+1),ja(*),jao(*)
      complex*8  a(*),ao(*)
c-----------------------------------------------------------------------
c Compressed Sparse Row     to      Compressed Sparse Column
c
c (transposition operation)   Not in place. 
c----------------------------------------------------------------------- 
c Rectangular version.  n is number of rows of CSR matrix,
c                       n2 (input) is number of columns of CSC matrix.
c----------------------------------------------------------------------- 
c -- not in place --
c this subroutine transposes a matrix stored in a, ja, ia format.
c ---------------
c on entry:
c----------
c n	= number of rows of CSR matrix.
c n2    = number of columns of CSC matrix.
c job	= integer to indicate whether to fill the values (job.eq.1) of the
c         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
c
c ipos  = starting position in ao, jao of the transposed matrix.
c         the iao array takes this into account (thus iao(1) is set to ipos.)
c         Note: this may be useful if one needs to append the data structure
c         of the transpose to that of A. In this case use for example
c                call csrcsc2 (n,n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2)) 
c	  for any other normal usage, enter ipos=1.
c a	= real array of length nnz (nnz=number of nonzero elements in input 
c         matrix) containing the nonzero elements.
c ja	= integer array of length nnz containing the column positions
c 	  of the corresponding elements in a.
c ia	= integer of size n+1. ia(k) contains the position in a, ja of
c	  the beginning of the k-th row.
c
c on return:
c ---------- 
c output arguments:
c ao	= real array of size nzz containing the "a" part of the transpose
c jao	= integer array of size nnz containing the column indices.
c iao	= integer array of size n+1 containing the "ia" index array of
c	  the transpose. 
c
c----------------------------------------------------------------------- 
c----------------- compute lengths of rows of transp(A) ----------------
      do 1 i=1,n2+1
         iao(i) = 0
 1    continue
      do 3 i=1, n
         do 2 k=ia(i), ia(i+1)-1 
            j = ja(k)+1
            iao(j) = iao(j)+1
 2       continue 
 3    continue
c---------- compute pointers from lengths ------------------------------
      iao(1) = ipos 
      do 4 i=1,n2
         iao(i+1) = iao(i) + iao(i+1)
 4    continue
c--------------- now do the actual copying ----------------------------- 
      do 6 i=1,n
         do 62 k=ia(i),ia(i+1)-1 
            j = ja(k) 
            next = iao(j)
            if (job .eq. 1)  ao(next) = a(k)
            jao(next) = i
            iao(j) = next+1
 62      continue
 6    continue
c-------------------------- reshift iao and leave ---------------------- 
      do 7 i=n2,1,-1
         iao(i+1) = iao(i)
 7    continue
      iao(1) = ipos
c--------------- end of csrcsc2 ---------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine ccsrlnk (n,a,ja,ia,link) 
      complex*8 a(*) 
      integer n, ja(*), ia(n+1), link(*)
c----------------------------------------------------------------------- 
c      Compressed Sparse Row         to    Linked storage format. 
c----------------------------------------------------------------------- 
c this subroutine translates a matrix stored in compressed sparse
c row into one with a linked list storage format. Only the link
c array needs to be obtained since the arrays a, ja, and ia may
c be unchanged and  carry the same meaning for the output matrix.
c in  other words a, ja, ia, link   is the output linked list data
c structure with a, ja, unchanged from input, and ia possibly 
c altered (in case therea re null rows in matrix). Details on
c the output array link are given below.
c----------------------------------------------------------------------- 
c Coded by Y. Saad, Feb 21, 1991.
c----------------------------------------------------------------------- 
c
c on entry:
c----------
c n	= integer equal to the dimension of A.	
c         
c a	= real array of size nna containing the nonzero elements
c ja	= integer array of size	nnz containing the column positions
c 	  of the corresponding elements in a.
c ia	= integer of size n+1 containing the pointers to the beginning 
c         of each row. ia(k) contains the position in a, ja of the 
c         beginning of the k-th row.
c
c on return:
c---------- 
c a, ja, are not changed.
c ia    may be changed if there are null rows.
c 
c a     = nonzero elements.
c ja    = column positions. 
c ia    = ia(i) points to the first element of row i in linked structure.
c link	= integer array of size containing the linked list information.
c         link(k) points to the next element of the row after element 
c         a(k), ja(k). if link(k) = 0, then there is no next element,
c         i.e., a(k), jcol(k) is the last element of the current row.
c
c  Thus row number i can be accessed as follows:
c     next = ia(i) 
c     while(next .ne. 0) do 
c          value = a(next)      ! value a(i,j) 
c          jcol  = ja(next)     ! column index j
c          next  = link(next)   ! address of next element in row
c     endwhile
c notes:
c ------ ia may be altered on return.
c----------------------------------------------------------------------- 
c local variables
      integer i, k
c
c loop through all rows
c
      do 100 i =1, n
         istart = ia(i) 
         iend = ia(i+1)-1
         if (iend .gt. istart) then
            do 99  k=istart, iend-1 
               link(k) = k+1
 99         continue
            link(iend) = 0
         else
            ia(i) = 0
         endif
 100  continue
c     
      return
c-------------end-of-csrlnk --------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine clnkcsr (n, a, jcol, istart, link, ao, jao, iao) 
      complex*8 a(*), ao(*) 
      integer n, jcol(*), istart(n), link(*), jao(*), iao(*) 
c----------------------------------------------------------------------- 
c     Linked list storage format   to      Compressed Sparse Row  format
c----------------------------------------------------------------------- 
c this subroutine translates a matrix stored in linked list storage 
c format into the compressed sparse row format. 
c----------------------------------------------------------------------- 
c Coded by Y. Saad, Feb 21, 1991.
c----------------------------------------------------------------------- 
c
c on entry:
c----------
c n	= integer equal to the dimension of A.	
c         
c a	= real array of size nna containing the nonzero elements
c jcol	= integer array of size	nnz containing the column positions
c 	  of the corresponding elements in a.
c istart= integer array of size n poiting to the beginning of the rows.
c         istart(i) contains the position of the first element of 
c         row i in data structure. (a, jcol, link).
c         if a row is empty istart(i) must be zero.
c link	= integer array of size nnz containing the links in the linked 
c         list data structure. link(k) points to the next element 
c         of the row after element ao(k), jcol(k). if link(k) = 0, 
c         then there is no next element, i.e., ao(k), jcol(k) is 
c         the last element of the current row.
c
c on return:
c-----------
c ao, jao, iao = matrix stored in csr format:
c
c ao    = real array containing the values of the nonzero elements of 
c         the matrix stored row-wise. 
c jao	= integer array of size nnz containing the column indices.
c iao	= integer array of size n+1 containing the pointers array to the 
c         beginning of each row. iao(i) is the address in ao,jao of
c         first element of row i.
c
c----------------------------------------------------------------------- 
c first determine individial bandwidths and pointers.
c----------------------------------------------------------------------- 
c local variables
      integer irow, ipos, next
c-----------------------------------------------------------------------
      ipos = 1
      iao(1) = ipos
c     
c     loop through all rows
c     
      do 100 irow =1, n
c     
c     unroll i-th row.
c     
         next = istart(irow)
 10      if (next .eq. 0) goto 99
         jao(ipos) = jcol(next)
         ao(ipos)  = a(next)
         ipos = ipos+1
         next = link(next) 
         goto 10
 99      iao(irow+1) = ipos 
 100  continue
c     
      return
c-------------end-of-lnkcsr ------------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine ccsrdia (n,idiag,job,a,ja,ia,ndiag,
     *                   diag,ioff,ao,jao,iao,ind)
      complex*8 diag(ndiag,idiag), a(*), ao(*)
      integer ia(*), ind(*), ja(*), jao(*), iao(*), ioff(*)
c----------------------------------------------------------------------- 
c Compressed sparse row     to    diagonal format
c----------------------------------------------------------------------- 
c this subroutine extracts  idiag diagonals  from the  input matrix a,
c a, ia, and puts the rest of  the matrix  in the  output matrix ao,
c jao, iao.  The diagonals to be extracted depend  on the  value of job
c (see below for details.)  In  the first  case, the  diagonals to be
c extracted are simply identified by  their offsets  provided in ioff
c by the caller.  In the second case, the  code internally determines
c the idiag most significant diagonals, i.e., those  diagonals of the
c matrix which  have  the  largest  number  of  nonzero elements, and
c extracts them.
c----------------------------------------------------------------------- 
c on entry:
c---------- 
c n	= dimension of the matrix a.
c idiag = integer equal to the number of diagonals to be extracted. 
c         Note: on return idiag may be modified.
c a, ja, 			
c    ia = matrix stored in a, ja, ia, format
c job	= integer. serves as a job indicator.  Job is better thought 
c         of as a two-digit number job=xy. If the first (x) digit
c         is one on entry then the diagonals to be extracted are 
c         internally determined. In this case csrdia exctracts the
c         idiag most important diagonals, i.e. those having the largest
c         number on nonzero elements. If the first digit is zero
c         then csrdia assumes that ioff(*) contains the offsets 
c         of the diagonals to be extracted. there is no verification 
c         that ioff(*) contains valid entries.
c         The second (y) digit of job determines whether or not
c         the remainder of the matrix is to be written on ao,jao,iao.
c         If it is zero  then ao, jao, iao is not filled, i.e., 
c         the diagonals are found  and put in array diag and the rest is
c         is discarded. if it is one, ao, jao, iao contains matrix
c         of the remaining elements.
c         Thus:
c         job= 0 means do not select diagonals internally (pick those
c                defined by ioff) and do not fill ao,jao,iao
c         job= 1 means do not select diagonals internally 
c                      and fill ao,jao,iao
c         job=10 means  select diagonals internally 
c                      and do not fill ao,jao,iao
c         job=11 means select diagonals internally 
c                      and fill ao,jao,iao
c 
c ndiag = integer equal to the first dimension of array diag.
c
c on return:
c----------- 
c
c idiag = number of diagonals found. This may be smaller than its value 
c         on entry. 
c diag  = real array of size (ndiag x idiag) containing the diagonals
c         of A on return
c          
c ioff  = integer array of length idiag, containing the offsets of the
c   	  diagonals to be extracted.
c ao, jao
c  iao  = remainder of the matrix in a, ja, ia format.
c work arrays:
c------------ 
c ind   = integer array of length 2*n-1 used as integer work space.
c         needed only when job.ge.10 i.e., in case the diagonals are to
c         be selected internally.
c
c Notes:
c-------
c    1) The algorithm is in place: ao, jao, iao can be overwritten on 
c       a, ja, ia if desired 
c    2) When the code is required to select the diagonals (job .ge. 10) 
c       the selection of the diagonals is done from left to right 
c       as a result if several diagonals have the same weight (number 
c       of nonzero elemnts) the leftmost one is selected first.
c-----------------------------------------------------------------------
      job1 = job/10
      job2 = job-job1*10
      if (job1 .eq. 0) goto 50
      n2 = n+n-1
      call infdia (n,ja,ia,ind,idum)
c----------- determine diagonals to  accept.---------------------------- 
c----------------------------------------------------------------------- 
      ii = 0
 4    ii=ii+1
      jmax = 0
      do 41 k=1, n2
         j = ind(k)
         if (j .le. jmax) goto 41
         i = k
         jmax = j
 41   continue
      if (jmax .le. 0) then
         ii = ii-1
         goto 42
      endif
      ioff(ii) = i-n
      ind(i) = - jmax
      if (ii .lt.  idiag) goto 4
 42   idiag = ii
c---------------- initialize diago to zero ----------------------------- 
 50   continue
      do 55 j=1,idiag
         do 54 i=1,n
            diag(i,j) = 0.0d0
 54      continue
 55   continue
c----------------------------------------------------------------------- 
      ko = 1
c----------------------------------------------------------------------- 
c extract diagonals and accumulate remaining matrix.
c----------------------------------------------------------------------- 
      do 6 i=1, n
         do 51 k=ia(i),ia(i+1)-1 
            j = ja(k)
            do 52 l=1,idiag
               if (j-i .ne. ioff(l)) goto 52
               diag(i,l) = a(k)
               goto 51
 52         continue
c--------------- append element not in any diagonal to ao,jao,iao ----- 
            if (job2 .eq. 0) goto 51
            ao(ko) = a(k)
            jao(ko) = j
            ko = ko+1
 51      continue
         if (job2 .ne. 0 ) ind(i+1) = ko
 6    continue
      if (job2 .eq. 0) return
c     finish with iao
      iao(1) = 1
      do 7 i=2,n+1
         iao(i) = ind(i)
 7    continue
      return
c----------- end of csrdia ---------------------------------------------
c-----------------------------------------------------------------------
      end

c----------------------------------------------------------------------- 
c  The following routine has been modified by Travis Oliphant (1999)
c      to allow for description of diagonal along non-square arrays
      subroutine cdiacsr (m,n,job,idiag,diag,ndiag,ioff,a,ja,ia)
      complex*8 diag(ndiag,idiag), a(*), t
      integer ia(*), ja(*), ioff(*)
c----------------------------------------------------------------------- 
c    diagonal format     to     compressed sparse row     
c----------------------------------------------------------------------- 
c this subroutine extract the idiag most important diagonals from the 
c input matrix a, ja, ia, i.e, those diagonals of the matrix which have
c the largest number of nonzero elements. If requested (see job),
c the rest of the matrix is put in a the output matrix ao, jao, iao
c----------------------------------------------------------------------- 
c on entry:
c---------- 
c m     = integer. number of rows in A
c n	= integer. number of columns in A.  Sparse array is (m x n)
c job	= integer. job indicator with the following meaning.
c         if (job .eq. 0) then check for each entry in diag
c         whether this entry is zero. If it is then do not include
c         in the output matrix. Note that the test is a test for
c         an exact arithmetic zero. Be sure that the zeros are
c         actual zeros in double precision otherwise this would not
c         work.
c         
c idiag = integer equal to the number of diagonals to be extracted. 
c         Note: on return idiag may be modified.
c
c diag  = real array of size (ndiag x idiag) containing the diagonals
c         of A on return. 
c 
c ndiag = integer equal to the first dimension of array diag.
c
c ioff  = integer array of length idiag, containing the offsets of the
c   	  diagonals to be extracted.
c
c on return:
c----------- 
c a, 
c ja, 			
c ia    = matrix stored in a, ja, ia, format
c
c Note:
c ----- the arrays a and ja should be of length n*idiag.
c
c----------------------------------------------------------------------- 
      ia(1) = 1
      ko = 1
      do 80 i=1, m
         do 70 jj = 1, idiag
            j = i+ioff(jj) 
            if (j .lt. 1 .or. j .gt. n) goto 70
            t = diag(i,jj) 
            if (job .eq. 0 .and. t .eq. 0.0d0) goto 70
            a(ko) = t
            ja(ko) = j
            ko = ko+1
 70      continue
         ia(i+1) = ko
 80   continue
      return
c----------- end of diacsr ---------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine cbsrcsr (job, n, m, na, a, ja, ia, ao, jao, iao)
      implicit none
      integer job, n, m, na, ia(*), ja(*), jao(*), iao(n+1)
      complex*8 a(na,*), ao(*)
c-----------------------------------------------------------------------
c             Block Sparse Row  to Compressed Sparse Row.
c----------------------------------------------------------------------- 
c NOTE: ** meanings of parameters may have changed wrt earlier versions
c FORMAT DEFINITION HAS CHANGED WRT TO EARLIER VERSIONS... 
c-----------------------------------------------------------------------
c
c converts a  matrix stored in block-reduced   a, ja, ia  format to the
c general  sparse row a,  ja, ia  format.  A matrix   that has  a block
c structure is a matrix whose entries are blocks  of the same size m
c (e.g.  3 x 3).   Then it is often preferred  to work with the reduced
c graph of the matrix. Instead of storing one element at a time one can
c store a whole block at a time.  In this storage scheme  an entry is a
c square array holding the m**2 elements of a block.
c 
c-----------------------------------------------------------------------
c on entry:
c----------
c job   = if job.eq.0 on entry, values are not copied (pattern only)
c
c n	= the block row dimension of the matrix.
c
c m     = the dimension of each block. Thus, the actual row dimension 
c         of A is n x m.  
c
c na	= first dimension of array a as declared in calling program.
c         This should be .ge. m**2.
c
c a	= real array containing the real entries of the matrix. Recall
c         that each entry is in fact an m x m block. These entries 
c         are stored column-wise in locations a(1:m*m,k) for each k-th
c         entry. See details below.
c 
c ja	= integer array of length n. ja(k) contains the column index 
c         of the leading element, i.e., the element (1,1) of the block
c         that is held in the column a(*,k) of the value array. 
c
c ia    = integer array of length n+1. ia(i) points to the beginning
c         of block row number i in the arrays a and ja. 
c 
c on return:
c-----------
c ao, jao, 
c iao   = matrix stored in compressed sparse row format. The number of
c         rows in the new matrix is n x m. 
c 
c Notes: THIS CODE IS NOT IN PLACE.
c 
c-----------------------------------------------------------------------
c BSR FORMAT.
c---------- 
c Each row of A contains the m x m block matrix unpacked column-
c wise (this allows the user to declare the array a as a(m,m,*) on entry 
c if desired). The block rows are stored in sequence just as for the
c compressed sparse row format. 
c
c-----------------------------------------------------------------------
c     example  with m = 2:
c                                                       1  2 3   
c    +-------|--------|--------+                       +-------+
c    | 1   2 |  0   0 |  3   4 |     Block             | x 0 x | 1
c    | 5   6 |  0   0 |  7   8 |     Representation:   | 0 x x | 2 
c    +-------+--------+--------+                       | x 0 0 | 3 
c    | 0   0 |  9  10 | 11  12 |                       +-------+ 
c    | 0   0 | 13  14 | 15  16 |  
c    +-------+--------+--------+   
c    | 17 18 |  0   0 |  0   0 |
c    | 22 23 |  0   0 |  0   0 |
c    +-------+--------+--------+
c
c    For this matrix:     n    = 3
c                         m    = 2
c                         nnz  = 5 
c-----------------------------------------------------------------------
c Data structure in Block Sparse Row format:      
c-------------------------------------------
c Array A:
c------------------------- 
c     1   3   9   11   17   <<--each m x m block is stored column-wise 
c     5   7   13  15   22       in a  column of the array A.
c     2   4   10  12   18      
c     6   8   14  16   23
c------------------------- 
c JA  1   3   2    3    1   <<-- column indices for each block. Note that
c-------------------------       these indices are wrt block matrix.
c IA  1   3   5    6        <<-- pointers to beginning of each block row 
c-------------------------       in arrays A and JA. 
c-----------------------------------------------------------------------
c locals 
c 
      integer i, i1, i2, ij, ii, irow, j, jstart, k, krow, no
      logical val
c     
      val = (job.ne.0)
      no = n * m 
      irow = 1	
      krow = 1
      iao(irow) = 1      
c-----------------------------------------------------------------------
      do 2 ii=1, n
c
c     recall: n is the block-row dimension
c
         i1 = ia(ii)
         i2 = ia(ii+1)-1
c
c     create m rows for each block row -- i.e., each k. 
c     
         do 23 i=1,m
            do 21 k=i1, i2
               jstart = m*(ja(k)-1)
               do 22  j=1,m
                  ij = (j-1)*m + i
                  if (val) ao(krow) = a(ij,k) 
                  jao(krow) = jstart+j
                  krow = krow+1
 22            continue	    
 21         continue
            irow = irow+1 
            iao(irow) = krow
 23      continue
 2    continue
      return
c-------------end-of-bsrcsr -------------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine ccsrbsr (job,nrow,m,na,a,ja,ia,ao,jao,iao,iw,ierr)
      implicit none
      integer job,ierr,nrow,m,na,ia(nrow+1),ja(*),jao(na),iao(*),iw(*)
      complex*8 a(*),ao(na,*)
c-----------------------------------------------------------------------
c     Compressed Sparse Row  to    Block Sparse Row
c-----------------------------------------------------------------------
c
c This  subroutine converts a matrix stored  in a general compressed a,
c ja, ia format into a a block  sparse row format a(m,m,*),ja(*),ia(*).
c See routine  bsrcsr  for more  details on  data   structure for block
c matrices. 
c
c NOTES: 1) the initial matrix does not have to have a block structure. 
c zero padding is done for general sparse matrices. 
c        2) For most practical purposes, na should be the same as m*m.
c 
c-----------------------------------------------------------------------
c 
c In what follows nr=1+(nrow-1)/m = block-row dimension of output matrix 
c 
c on entry:
c----------
c 
c job   =  job indicator.
c          job =  0 -> only the pattern of output matrix is generated 
c          job >  0 -> both pattern and values are generated. 
c          job = -1 -> iao(1) will return the number of nonzero blocks,
c            in the output matrix. In this case jao(1:nr) is used as 
c            workspace, ao is untouched, iao is untouched except iao(1)
c
c nrow	= integer, the actual row dimension of the matrix.
c 
c m     = integer equal to the dimension of each block. m should be > 0. 
c 
c na	= first dimension of array ao as declared in calling program.
c         na should be .ge. m*m. 
c
c a, ja, 
c    ia = input matrix stored in compressed sparse row format.
c
c on return:
c-----------
c 
c ao    = real  array containing the  values of the matrix. For details
c         on the format  see below. Each  row of  a contains the  m x m
c         block matrix  unpacked column-wise (this  allows the  user to
c         declare the  array a as ao(m,m,*) on  entry if desired).  The
c         block rows are stored in sequence  just as for the compressed
c         sparse row format. The block  dimension  of the output matrix
c         is  nr = 1 + (nrow-1) / m.
c         
c jao   = integer array. containing the block-column indices of the 
c         block-matrix. Each jao(k) is an integer between 1 and nr
c         containing the block column index of the block ao(*,k).  
c
c iao   = integer array of length nr+1. iao(i) points to the beginning
c         of block row number i in the arrays ao and jao. When job=-1
c         iao(1) contains the number of nonzero blocks of the output
c         matrix and the rest of iao is unused. This is useful for
c         determining the lengths of ao and jao. 
c
c ierr  = integer, error code. 
c              0 -- normal termination
c              1 -- m is equal to zero 
c              2 -- NA too small to hold the blocks (should be .ge. m**2)
c
c Work arrays:
c------------- 
c iw    = integer work array of dimension  nr = 1 + (nrow-1) / m
c
c NOTES: 
c-------
c     1) this code is not in place.
c     2) see routine bsrcsr for details on data sctructure for block 
c        sparse row format.
c     
c-----------------------------------------------------------------------
c     nr is the block-dimension of the output matrix.
c     
      integer nr, m2, io, ko, ii, len, k, jpos, j, i, ij, jr, irow   
      logical vals  
c----- 
      ierr = 0 
      if (m*m .gt. na) ierr = 2 
      if (m .eq. 0) ierr = 1 
      if (ierr .ne. 0) return
c----------------------------------------------------------------------- 
      vals = (job .gt. 0) 
      nr = 1 + (nrow-1) / m
      m2 = m*m 
      ko = 1 
      io = 1 
      iao(io) = 1 
      len = 0 
c     
c     iw determines structure of block-row (nonzero indicator) 
c 
         do j=1, nr
            iw(j) = 0
         enddo
c     
c     big loop -- leap by m rows each time.
c     
      do ii=1, nrow, m
         irow = 0
c
c     go through next m rows -- make sure not to go beyond nrow. 
c
         do while (ii+irow .le. nrow .and. irow .le. m-1) 
            do k=ia(ii+irow),ia(ii+irow+1)-1
c     
c     block column index = (scalar column index -1) / m + 1 
c
               j = ja(k)-1 
               jr = j/m + 1                
               j = j - (jr-1)*m 
               jpos = iw(jr) 
               if (jpos .eq. 0) then
c
c     create a new block
c     
                  iw(jr) = ko 
                  jao(ko) = jr 
                  if (vals) then
c     
c     initialize new block to zero -- then copy nonzero element
c     
                     do i=1, m2
                        ao(i,ko) = 0.0d0
                     enddo
                     ij = j*m + irow + 1 
                     ao(ij,ko) = a(k) 
                  endif
                  ko = ko+1
               else
c
c     copy column index and nonzero element 
c     
                  jao(jpos) = jr 
                  ij = j*m + irow + 1 
                  if (vals) ao(ij,jpos) = a(k) 
               endif 
            enddo  
            irow = irow+1
         enddo
c     
c     refresh iw
c                      
         do j = iao(io),ko-1 
            iw(jao(j)) = 0
         enddo
         if (job .eq. -1) then
            len = len + ko-1
            ko = 1
         else
            io = io+1 
            iao(io) = ko
         endif
      enddo
      if (job .eq. -1) iao(1) = len
c
      return
c--------------end-of-csrbsr-------------------------------------------- 
c----------------------------------------------------------------------- 
      end
c-----------------------------------------------------------------------
      subroutine ccsrbnd (n,a,ja,ia,job,abd,nabd,lowd,ml,mu,ierr)
      complex*8 a(*),abd(nabd,n)
      integer ia(n+1),ja(*)
c----------------------------------------------------------------------- 
c   Compressed Sparse Row  to  Banded (Linpack ) format.
c----------------------------------------------------------------------- 
c this subroutine converts a general sparse matrix stored in
c compressed sparse row format into the banded format. for the
c banded format,the Linpack conventions are assumed (see below).
c----------------------------------------------------------------------- 
c on entry:
c----------
c n	= integer,the actual row dimension of the matrix.
c
c a,
c ja,
c ia    = input matrix stored in compressed sparse row format.
c
c job	= integer. if job=1 then the values of the lower bandwith ml 
c         and the upper bandwidth mu are determined internally. 
c         otherwise it is assumed that the values of ml and mu 
c         are the correct bandwidths on input. See ml and mu below.
c
c nabd  = integer. first dimension of array abd.
c
c lowd  = integer. this should be set to the row number in abd where
c         the lowest diagonal (leftmost) of A is located. 
c         lowd should be  ( 1  .le.  lowd  .le. nabd).
c         if it is not known in advance what lowd should be
c         enter lowd = 0 and the default value lowd = ml+mu+1
c         will be chosen. Alternative: call routine getbwd from unary
c         first to detrermione ml and mu then define lowd accordingly.
c         (Note: the banded solvers in linpack use lowd=2*ml+mu+1. )
c
c ml	= integer. equal to the bandwidth of the strict lower part of A
c mu	= integer. equal to the bandwidth of the strict upper part of A
c         thus the total bandwidth of A is ml+mu+1.
c         if ml+mu+1 is found to be larger than lowd then an error 
c         flag is raised (unless lowd = 0). see ierr.
c
c note:   ml and mu are assumed to have	 the correct bandwidth values
c         as defined above if job is set to zero on entry.
c
c on return:
c-----------
c
c abd   = real array of dimension abd(nabd,n).
c         on return contains the values of the matrix stored in
c         banded form. The j-th column of abd contains the elements
c         of the j-th column of  the original matrix comprised in the
c         band ( i in (j-ml,j+mu) ) with the lowest diagonal at
c         the bottom row (row lowd). See details below for this format.
c
c ml	= integer. equal to the bandwidth of the strict lower part of A
c mu	= integer. equal to the bandwidth of the strict upper part of A
c         if job=1 on entry then these two values are internally computed.
c
c lowd  = integer. row number in abd where the lowest diagonal 
c         (leftmost) of A is located on return. In case lowd = 0
c         on return, then it is defined to ml+mu+1 on return and the
c         lowd will contain this value on return. `
c
c ierr  = integer. used for error messages. On return:
c         ierr .eq. 0  :means normal return
c         ierr .eq. -1 : means invalid value for lowd. (either .lt. 0
c         or larger than nabd).
c         ierr .eq. -2 : means that lowd is not large enough and as 
c         result the matrix cannot be stored in array abd. 
c         lowd should be at least ml+mu+1, where ml and mu are as
c         provided on output.
c
c----------------------------------------------------------------------* 
c Additional details on banded format.  (this closely follows the      *
c format used in linpack. may be useful for converting a matrix into   *
c this storage format in order to use the linpack  banded solvers).    * 
c----------------------------------------------------------------------*
c             ---  band storage format  for matrix abd ---             * 
c uses ml+mu+1 rows of abd(nabd,*) to store the diagonals of           *
c a in rows of abd starting from the lowest (sub)-diagonal  which  is  *
c stored in row number lowd of abd. the minimum number of rows needed  *
c in abd is ml+mu+1, i.e., the minimum value for lowd is ml+mu+1. the  *
c j-th  column  of  abd contains the elements of the j-th column of a, *
c from bottom to top: the element a(j+ml,j) is stored in  position     *
c abd(lowd,j), then a(j+ml-1,j) in position abd(lowd-1,j) and so on.   *
c Generally, the element a(j+k,j) of original matrix a is stored in    *
c position abd(lowd+k-ml,j), for k=ml,ml-1,..,0,-1, -mu.               *
c The first dimension nabd of abd must be .ge. lowd                    *
c                                                                      *
c     example [from linpack ]:   if the original matrix is             *
c                                                                      *
c              11 12 13  0  0  0                                       *
c              21 22 23 24  0  0                                       *
c               0 32 33 34 35  0     original banded matrix            *
c               0  0 43 44 45 46                                       *
c               0  0  0 54 55 56                                       *
c               0  0  0  0 65 66                                       *
c                                                                      *
c then  n = 6, ml = 1, mu = 2. lowd should be .ge. 4 (=ml+mu+1)  and   *
c if lowd = 5 for example, abd  should be:                             *
c                                                                      *
c untouched --> x  x  x  x  x  x                                       *
c               *  * 13 24 35 46                                       *
c               * 12 23 34 45 56    resulting abd matrix in banded     *
c              11 22 33 44 55 66    format                             *
c  row lowd--> 21 32 43 54 65  *                                       *
c                                                                      *
c * = not used                                                         *
c                                                                      
*
c----------------------------------------------------------------------*
c first determine ml and mu.
c----------------------------------------------------------------------- 
      ierr = 0
c-----------
      if (job .eq. 1) call cgetbwd(n,a,ja,ia,ml,mu)
      m = ml+mu+1
      if (lowd .eq. 0) lowd = m
      if (m .gt. lowd)  ierr = -2
      if (lowd .gt. nabd .or. lowd .lt. 0) ierr = -1
      if (ierr .lt. 0) return
c------------
      do 15  i=1,m
         ii = lowd -i+1
         do 10 j=1,n
	    abd(ii,j) = 0.0d0
 10      continue
 15   continue
c---------------------------------------------------------------------	   
      mdiag = lowd-ml
      do 30 i=1,n
         do 20 k=ia(i),ia(i+1)-1
            j = ja(k)
            abd(i-j+mdiag,j) = a(k) 
 20      continue
 30   continue
      return
c------------- end of csrbnd ------------------------------------------- 
c----------------------------------------------------------------------- 
      end
c-----------------------------------------------------------------------
      subroutine cbndcsr (n,abd,nabd,lowd,ml,mu,a,ja,ia,len,ierr)
      complex*8 a(*),abd(nabd,*), t
      integer ia(n+1),ja(*)
c----------------------------------------------------------------------- 
c Banded (Linpack ) format   to    Compressed Sparse Row  format.
c----------------------------------------------------------------------- 
c on entry:
c----------
c n	= integer,the actual row dimension of the matrix.
c
c nabd  = first dimension of array abd.
c
c abd   = real array containing the values of the matrix stored in
c         banded form. The j-th column of abd contains the elements
c         of the j-th column of  the original matrix,comprised in the
c         band ( i in (j-ml,j+mu) ) with the lowest diagonal located
c         in row lowd (see below). 
c    
c lowd  = integer. this should be set to the row number in abd where
c         the lowest diagonal (leftmost) of A is located. 
c         lowd should be s.t.  ( 1  .le.  lowd  .le. nabd).
c         The subroutines dgbco, ... of linpack use lowd=2*ml+mu+1.
c
c ml	= integer. equal to the bandwidth of the strict lower part of A
c mu	= integer. equal to the bandwidth of the strict upper part of A
c         thus the total bandwidth of A is ml+mu+1.
c         if ml+mu+1 is found to be larger than nabd then an error 
c         message is set. see ierr.
c
c len   = integer. length of arrays a and ja. bndcsr will stop if the
c         length of the arrays a and ja is insufficient to store the 
c         matrix. see ierr.
c
c on return:
c-----------
c a,
c ja,
c ia    = input matrix stored in compressed sparse row format.
c
c lowd  = if on entry lowd was zero then lowd is reset to the default
c         value ml+mu+l. 
c
c ierr  = integer. used for error message output. 
c         ierr .eq. 0 :means normal return
c         ierr .eq. -1 : means invalid value for lowd. 
c	  ierr .gt. 0 : means that there was not enough storage in a and ja
c         for storing the ourput matrix. The process ran out of space 
c         (as indicated by len) while trying to fill row number ierr. 
c         This should give an idea of much more storage might be required. 
c         Moreover, the first irow-1 rows are correctly filled. 
c
c notes:  the values in abd found to be equal to zero
c -----   (actual test: if (abd(...) .eq. 0.0d0) are removed.
c         The resulting may not be identical to a csr matrix
c         originally transformed to a bnd format.
c          
c----------------------------------------------------------------------- 
      ierr = 0
c-----------
      if (lowd .gt. nabd .or. lowd .le. 0) then 
         ierr = -1
         return
      endif
c-----------
      ko = 1
      ia(1) = 1
      do 30 irow=1,n
c-----------------------------------------------------------------------
         i = lowd 
          do  20 j=irow-ml,irow+mu
             if (j .le. 0 ) goto 19
             if (j .gt. n) goto 21
             t = abd(i,j) 
             if (t .eq. 0.0d0) goto 19
             if (ko .gt. len) then 
               ierr = irow 
               return
            endif
            a(ko) = t
            ja(ko) = j
            ko = ko+1
 19         i = i-1
 20      continue
c     end for row irow
 21      ia(irow+1) = ko
 30   continue
      return
c------------- end of bndcsr ------------------------------------------- 
c----------------------------------------------------------------------- 
      end
c----------------------------------------------------------------------- 
      subroutine ccsrssk (n,imod,a,ja,ia,asky,isky,nzmax,ierr)
      complex*8 a(*),asky(nzmax) 
      integer n, imod, nzmax, ierr, ia(n+1), isky(n+1), ja(*)
c----------------------------------------------------------------------- 
c      Compressed Sparse Row         to     Symmetric Skyline Format 
c  or  Symmetric Sparse Row        
c----------------------------------------------------------------------- 
c this subroutine translates a compressed sparse row or a symmetric
c sparse row format into a symmetric skyline format.
c the input matrix can be in either compressed sparse row or the 
c symmetric sparse row format. The output matrix is in a symmetric
c skyline format: a real array containing the (active portions) of the
c rows in  sequence and a pointer to the beginning of each row.
c
c This module is NOT  in place.
c----------------------------------------------------------------------- 
c Coded by Y. Saad, Oct 5, 1989. Revised Feb. 18, 1991.
c----------------------------------------------------------------------- 
c
c on entry:
c----------
c n	= integer equal to the dimension of A.	
c imod  = integer indicating the variant of skyline format wanted:
c         imod = 0 means the pointer isky points to the `zeroth' 
c         element of the row, i.e., to the position of the diagonal
c         element of previous row (for i=1, isky(1)= 0)
c         imod = 1 means that itpr points to the beginning of the row. 
c         imod = 2 means that isky points to the end of the row (diagonal
c                  element) 
c         
c a	= real array of size nna containing the nonzero elements
c ja	= integer array of size	nnz containing the column positions
c 	  of the corresponding elements in a.
c ia	= integer of size n+1. ia(k) contains the position in a, ja of
c	  the beginning of the k-th row.
c nzmax = integer. must be set to the number of available locations
c         in the output array asky. 
c
c on return:
c---------- 
c
c asky    = real array containing the values of the matrix stored in skyline
c         format. asky contains the sequence of active rows from 
c         i=1, to n, an active row being the row of elemnts of 
c         the matrix contained between the leftmost nonzero element 
c         and the diagonal element. 
c isky	= integer array of size n+1 containing the pointer array to 
c         each row. The meaning of isky depends on the input value of
c         imod (see above). 
c ierr  =  integer.  Error message. If the length of the 
c         output array asky exceeds nzmax. ierr returns the minimum value
c         needed for nzmax. otherwise ierr=0 (normal return).
c 
c Notes:
c         1) This module is NOT  in place.
c         2) even when imod = 2, length of  isky is  n+1, not n.
c
c----------------------------------------------------------------------- 
c first determine individial bandwidths and pointers.
c----------------------------------------------------------------------- 
      ierr = 0
      isky(1) = 0
      do 3 i=1,n
         ml = 0
         do 31 k=ia(i),ia(i+1)-1 
            ml = max(ml,i-ja(k)+1) 
 31      continue
         isky(i+1) = isky(i)+ml
 3    continue
c
c     test if there is enough space  asky to do the copying.  
c
      nnz = isky(n+1) 
      if (nnz .gt. nzmax) then
         ierr = nnz
         return
      endif
c    
c   fill asky with zeros.
c     
      do 1 k=1, nnz 
         asky(k) = 0.0d0
 1    continue
c     
c     copy nonzero elements.
c     
      do 4 i=1,n
         kend = isky(i+1) 
         do 41 k=ia(i),ia(i+1)-1 
            j = ja(k)
            if (j .le. i) asky(kend+j-i) = a(k)
 41      continue
 4    continue
c 
c modify pointer according to imod if necessary.
c
      if (imod .eq. 0) return
      if (imod .eq. 1) then 
         do 50 k=1, n+1
            isky(k) = isky(k)+1
 50      continue
      endif
      if (imod .eq. 2) then
         do 60 k=1, n
            isky(k) = isky(k+1) 
 60      continue
      endif
c
      return
c------------- end of csrssk ------------------------------------------- 
c----------------------------------------------------------------------- 
      end
c-----------------------------------------------------------------------
      subroutine csskssr (n,imod,asky,isky,ao,jao,iao,nzmax,ierr)
      complex*8 asky(*),ao(nzmax) 
      integer n, imod,nzmax,ierr, isky(n+1),iao(n+1),jao(nzmax) 
c----------------------------------------------------------------------- 
c     Symmetric Skyline Format  to  Symmetric Sparse Row format.
c----------------------------------------------------------------------- 
c  tests for exact zeros in skyline matrix (and ignores them in
c  output matrix).  In place routine (a, isky :: ao, iao)
c----------------------------------------------------------------------- 
c this subroutine translates a  symmetric skyline format into a 
c symmetric sparse row format. Each element is tested to see if it is
c a zero element. Only the actual nonzero elements are retained. Note 
c that the test used is simple and does take into account the smallness 
c of a value. the subroutine filter (see unary module) can be used
c for this purpose. 
c----------------------------------------------------------------------- 
c Coded by Y. Saad, Oct 5, 1989. Revised Feb 18, 1991./
c----------------------------------------------------------------------- 
c
c on entry:
c----------
c n	= integer equal to the dimension of A.	
c imod  = integer indicating the variant of skyline format used:
c         imod = 0 means the pointer iao points to the `zeroth' 
c         element of the row, i.e., to the position of the diagonal
c         element of previous row (for i=1, iao(1)= 0)
c         imod = 1 means that itpr points to the beginning of the row. 
c         imod = 2 means that iao points to the end of the row 
c                  (diagonal element) 
c asky  = real array containing the values of the matrix. asky contains 
c         the sequence of active rows from i=1, to n, an active row 
c         being the row of elemnts of the matrix contained between the 
c         leftmost nonzero element and the diagonal element. 
c isky 	= integer array of size n+1 containing the pointer array to 
c         each row. isky (k) contains the address of the beginning of the
c         k-th active row in the array asky. 
c nzmax = integer. equal to the number of available locations in the 
c         output array ao.  
c
c on return:
c ---------- 
c ao	= real array of size nna containing the nonzero elements
c jao	= integer array of size	nnz containing the column positions
c 	  of the corresponding elements in a.
c iao	= integer of size n+1. iao(k) contains the position in a, ja of
c	  the beginning of the k-th row.
c ierr  = integer. Serving as error message. If the length of the 
c         output arrays ao, jao exceeds nzmax then ierr returns 
c         the row number where the algorithm stopped: rows
c         i, to ierr-1 have been processed succesfully.
c         ierr = 0 means normal return.
c         ierr = -1  : illegal value for imod
c Notes:
c------- 
c This module is in place: ao and iao can be the same as asky, and isky.
c-----------------------------------------------------------------------
c local variables
      integer next, kend, kstart, i, j 
      ierr = 0
c
c check for validity of imod
c 
      if (imod.ne.0 .and. imod.ne.1 .and. imod .ne. 2) then
         ierr =-1
         return
      endif 
c
c next  = pointer to next available position in output matrix
c kend  = pointer to end of current row in skyline matrix. 
c
      next = 1
c
c set kend = start position -1 in  skyline matrix.
c 
      kend = 0 
      if (imod .eq. 1) kend = isky(1)-1
      if (imod .eq. 0) kend = isky(1) 
c
c loop through all rows
c     
      do 50 i=1,n
c
c save value of pointer to ith row in output matrix
c
         iao(i) = next
c
c get beginnning and end of skyline  row 
c
         kstart = kend+1
         if (imod .eq. 0) kend = isky(i+1)
         if (imod .eq. 1) kend = isky(i+1)-1
         if (imod .eq. 2) kend = isky(i) 
c 
c copy element into output matrix unless it is a zero element.
c 
         do 40 k=kstart,kend
            if (asky(k) .eq. 0.0d0) goto 40
            j = i-(kend-k) 
            jao(next) = j
            ao(next)  = asky(k)
            next=next+1
            if (next .gt. nzmax+1) then
               ierr = i
               return
            endif 
 40      continue
 50    continue
      iao(n+1) = next
      return
c-------------end-of-sskssr -------------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine ccsrjad (nrow, a, ja, ia, idiag, iperm, ao, jao, iao) 
      integer ja(*), jao(*), ia(nrow+1), iperm(nrow), iao(nrow) 
      complex*8 a(*), ao(*)
c-----------------------------------------------------------------------
c    Compressed Sparse Row  to   JAgged Diagonal storage. 
c----------------------------------------------------------------------- 
c this subroutine converts  matrix stored in the compressed sparse
c row format to the jagged diagonal format. The data structure
c for the JAD (Jagged Diagonal storage) is as follows. The rows of 
c the matrix are (implicitly) permuted so that their lengths are in
c decreasing order. The real entries ao(*) and their column indices 
c jao(*) are stored in succession. The number of such diagonals is idiag.
c the lengths of each of these diagonals is stored in iao(*).
c For more details see [E. Anderson and Y. Saad,
c ``Solving sparse triangular systems on parallel computers'' in
c Inter. J. of High Speed Computing, Vol 1, pp. 73-96 (1989).]
c or  [Y. Saad, ``Krylov Subspace Methods on Supercomputers''
c SIAM J. on  Stat. Scient. Comput., volume 10, pp. 1200-1232 (1989).]
c----------------------------------------------------------------------- 
c on entry:
c---------- 
c nrow 	  = row dimension of the matrix A.
c
c a, 
c ia, 
c ja      = input matrix in compressed sparse row format. 
c
c on return: 
c----------
c 
c idiag = integer. The number of jagged diagonals in the matrix.
c
c iperm = integer array of length nrow containing the permutation
c         of the rows that leads to a decreasing order of the
c         number of nonzero elements.
c
c ao    = real array containing the values of the matrix A in 
c         jagged diagonal storage. The j-diagonals are stored
c         in ao in sequence. 
c
c jao   = integer array containing the column indices of the 
c         entries in ao.
c
c iao   = integer array containing pointers to the beginning 
c         of each j-diagonal in ao, jao. iao is also used as 
c         a work array and it should be of length n at least.
c
c----------------------------------------------------------------------- 
c     ---- define initial iperm and get lengths of each row
c     ---- jao is used a work vector to store tehse lengths
c     
      idiag = 0
      ilo = nrow 
      do 10 j=1, nrow
         iperm(j) = j 
         len = ia(j+1) - ia(j)
         ilo = min(ilo,len) 
         idiag = max(idiag,len) 
         jao(j) = len
 10   continue 
c     
c     call sorter to get permutation. use iao as work array.
c    
      call dcsort (jao, nrow, iao, iperm, ilo, idiag) 
c     
c     define output data structure. first lengths of j-diagonals
c     
      do 20 j=1, nrow
         iao(j) = 0
 20   continue
      do 40 k=1, nrow
         len = jao(iperm(k)) 
         do 30 i=1,len
            iao(i) = iao(i)+1
 30      continue
 40   continue
c     
c     get the output matrix itself
c     
      k1 = 1
      k0 = k1
      do 60 jj=1, idiag
         len = iao(jj)
         do 50 k=1,len
            i = ia(iperm(k))+jj-1
            ao(k1) = a(i)
            jao(k1) = ja(i) 
            k1 = k1+1
 50      continue
         iao(jj) = k0
         k0 = k1
 60   continue
      iao(idiag+1) = k1
      return
c----------end-of-csrjad------------------------------------------------
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine cjadcsr (nrow, idiag, a, ja, ia, iperm, ao, jao, iao) 
      integer ja(*), jao(*), ia(idiag+1), iperm(nrow), iao(nrow+1) 
      complex*8 a(*), ao(*)
c-----------------------------------------------------------------------
c     Jagged Diagonal Storage   to     Compressed Sparse Row  
c----------------------------------------------------------------------- 
c this subroutine converts a matrix stored in the jagged diagonal format
c to the compressed sparse row format.
c----------------------------------------------------------------------- 
c on entry:
c---------- 
c nrow 	  = integer. the row dimension of the matrix A.
c 
c idiag   = integer. The  number of jagged diagonals in the data
c           structure a, ja, ia.
c 
c a, 
c ja,
c ia      = input matrix in jagged diagonal format. 
c
c iperm   = permutation of the rows used to obtain the JAD ordering. 
c        
c on return: 
c----------
c 
c ao, jao,
c iao     = matrix in CSR format.
c-----------------------------------------------------------------------
c determine first the pointers for output matrix. Go through the
c structure once:
c
      do 137 j=1,nrow
         jao(j) = 0
 137  continue
c     
c     compute the lengths of each row of output matrix - 
c     
      do 140 i=1, idiag
         len = ia(i+1)-ia(i) 
         do 138 k=1,len
            jao(iperm(k)) = jao(iperm(k))+1
 138     continue
 140  continue
c     
c     remember to permute
c     
      kpos = 1
      iao(1) = 1
      do 141 i=1, nrow 
         kpos = kpos+jao(i) 
         iao(i+1) = kpos
 141  continue
c     
c     copy elemnts one at a time.
c     
      do 200 jj = 1, idiag
         k1 = ia(jj)-1
         len = ia(jj+1)-k1-1 
         do 160 k=1,len
            kpos = iao(iperm(k))
            ao(kpos) = a(k1+k) 
            jao(kpos) = ja(k1+k) 
            iao(iperm(k)) = kpos+1
 160     continue
 200  continue
c     
c     rewind pointers
c     
      do 5 j=nrow,1,-1
         iao(j+1) = iao(j)
 5    continue
      iao(1) = 1
      return
c----------end-of-jadcsr------------------------------------------------
c-----------------------------------------------------------------------
      end
      subroutine ccooell (job,n,nnz,a,ja,ia,ao,jao,lda,ncmax,nc,ierr)
      implicit none
      integer job,n,nnz,lda,ncmax,nc,ierr
      integer ja(nnz),ia(nnz),jao(lda,ncmax)
      complex*8  a(nnz),ao(lda,ncmax)
c-----------------------------------------------------------------------
c     COOrdinate format to ELLpack format
c-----------------------------------------------------------------------
c     On entry:
c     job     -- 0 if only pattern is to be processed(AO is not touched)
c     n       -- number of rows in the matrix
c     a,ja,ia -- input matix in COO format
c     lda     -- leading dimension of array AO and JAO
c     ncmax   -- size of the second dimension of array AO and JAO
c
c     On exit:
c     ao,jao  -- the matrix in ELL format
c     nc      -- maximum number of nonzeros per row
c     ierr    -- 0 if convertion succeeded
c                -1 if LDA < N
c                nc if NC > ncmax
c
c     NOTE: the last column of JAO is used as work space!!
c-----------------------------------------------------------------------
      integer i,j,k,ip
      complex*8  zero
      logical copyval
      parameter (zero=0.0D0)
c     .. first executable statement ..
      copyval = (job.ne.0)
      if (lda .lt. n) then
         ierr = -1
         return
      endif
c     .. use the last column of JAO as workspace
c     .. initialize the work space
      do i = 1, n
         jao(i,ncmax) = 0
      enddo
      nc = 0
c     .. go through ia and ja to find out number nonzero per row
      do k = 1, nnz
         i = ia(k)
         jao(i,ncmax) = jao(i,ncmax) + 1
      enddo
c     .. maximum number of nonzero per row
      nc = 0
      do i = 1, n
         if (nc.lt.jao(i,ncmax)) nc = jao(i,ncmax)
         jao(i,ncmax) = 0
      enddo
c     .. if nc > ncmax retrun now
      if (nc.gt.ncmax) then
         ierr = nc
         return
      endif
c     .. go through ia and ja to copy the matrix to AO and JAO
      do k = 1, nnz
         i = ia(k)
         j = ja(k)
         jao(i,ncmax) = jao(i,ncmax) + 1
         ip = jao(i,ncmax)
         if (ip.gt.nc) nc = ip
         if (copyval) ao(i,ip) = a(k)
         jao(i,ip) = j
      enddo
c     .. fill the unspecified elements of AO and JAO with zero diagonals
      do i = 1, n
         do j = ia(i+1)-ia(i)+1, nc
            jao(i,j)=i
            if(copyval) ao(i,j) = zero
         enddo
      enddo
      ierr = 0
c
      return
      end
c-----end-of-cooell-----------------------------------------------------
c-----------------------------------------------------------------------
      subroutine cxcooell (n,nnz,a,ja,ia,ac,jac,nac,ner,ncmax,ierr)
C-----------------------------------------------------------------------
C   coordinate format to ellpack format.
C-----------------------------------------------------------------------
C
C   DATE WRITTEN: June 4, 1989. 
C
C   PURPOSE
C   -------
C  This subroutine takes a sparse matrix in coordinate format and
C  converts it into the Ellpack-Itpack storage.
C
C  Example:
C  -------
C       (   11   0   13    0     0     0  )
C       |   21  22    0   24     0     0  |
C       |    0  32   33    0    35     0  |
C   A = |    0   0   43   44     0    46  |
C       |   51   0    0   54    55     0  |
C       (   61  62    0    0    65    66  )
C
C   Coordinate storage scheme:
C
C    A  = (11,22,33,44,55,66,13,21,24,32,35,43,46,51,54,61,62,65)
C    IA = (1, 2, 3, 4, 5, 6, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 6 )
C    JA = ( 1, 2, 3, 4, 5, 6, 3, 1, 4, 2, 5, 3, 6, 1, 4, 1, 2, 5)
C
C   Ellpack-Itpack storage scheme:
C
C       (   11  13    0    0   )          (   1   3   *    *  )
C       |   22  21   24    0   |          |   2   1   4    *  |
C  AC = |   33  32   35    0   |    JAC = |   3   2   5    *  |
C       |   44  43   46    0   |          |   4   3   6    *  |
C       |   55  51   54    0   |          |   5   1   4    *  |
C       (   66  61   62   65   )          (   6   1   2    5  )
C
C   Note: * means that you can store values from 1 to 6 (1 to n, where
C         n is the order of the matrix) in that position in the array.
C
C   Contributed by:
C   --------------- 
C   Ernest E. Rothman
C   Cornell Thoery Center/Cornell National Supercomputer Facility
C   e-mail address: BITNET:   EER@CORNELLF.BITNET
C                   INTERNET: eer@cornellf.tn.cornell.edu
C   
C   checked and modified  04/13/90 Y.Saad.
C
C   REFERENCES
C   ----------
C   Kincaid, D. R.; Oppe, T. C.; Respess, J. R.; Young, D. M. 1984.
C   ITPACKV 2C User's Guide, CNA-191. Center for Numerical Analysis,
C   University of Texas at Austin.
C
C   "Engineering and Scientific Subroutine Library; Guide and
C   Reference; Release 3 (SC23-0184-3). Pp. 79-86.
C
C-----------------------------------------------------------------------
C
C   INPUT PARAMETERS
C   ----------------
C  N       - Integer. The size of the square matrix.
C
C  NNZ     - Integer. Must be greater than or equal to the number of
C            nonzero elements in the sparse matrix. Dimension of A, IA 
C            and JA.
C
C  NCA     - Integer. First dimension of output arrays ca and jac.
C
C  A(NNZ)  - Real array. (Double precision)
C            Stored entries of the sparse matrix A.
C            NNZ is the number of nonzeros.
C
C  IA(NNZ) - Integer array.
C            Pointers to specify rows for the stored nonzero entries
C            in A.
C
C  JA(NNZ) - Integer array.
C            Pointers to specify columns for the stored nonzero
C            entries in A.
C
C  NER     - Integer. Must be set greater than or equal to the maximum
C            number of nonzeros in any row of the sparse matrix.
C
C  OUTPUT PARAMETERS
C  -----------------
C  AC(NAC,*)  - Real array. (Double precision)
C               Stored entries of the sparse matrix A in compressed
C               storage mode.
C
C  JAC(NAC,*) - Integer array.
C               Contains the column numbers of the sparse matrix
C               elements stored in the corresponding positions in
C               array AC.
C
C  NCMAX   -  Integer. Equals the maximum number of nonzeros in any
C             row of the sparse matrix.
C
C  IERR    - Error parameter is returned as zero on successful
C             execution of the subroutin<e.
C             Error diagnostics are given by means of positive values
C             of this parameter as follows:
C
C             IERR = -1   -  NER is too small and should be set equal
C                            to NCMAX. The array AC may not be large
C                            enough to accomodate all the non-zeros of
C                            of the sparse matrix.
C             IERR =  1   -  The array AC has a zero column. (Warning) 
C             IERR =  2   -  The array AC has a zero row.    (Warning)
C
C---------------------------------------------------------------------
      complex*8 a(nnz), ac(nac,ner)
      integer ja(nnz), ia(nnz), jac(nac,ner), ierr, ncmax, icount
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   Initial error parameter to zero:
c
      ierr = 0
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   Initial output arrays to zero:
c
      do 4 in = 1,ner
         do 4 innz =1,n
            jac(innz,in) = n
            ac(innz,in) = 0.0d0
 4    continue
c     
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c   Assign nonzero elements of the sparse matrix (stored in the one
c   dimensional array A to the two dimensional array AC.
c   Also, assign the correct values with information about their
c   column indices to the two dimensional array KA. And at the same
c   time count the number of nonzeros in each row so that the
c   parameter NCMAX equals the maximum number of nonzeros in any row
c   of the sparse matrix.
c
      ncmax = 1
      do 10 is = 1,n
         k = 0
         do 30 ii = 1,nnz
            if(ia(ii).eq.is)then
               k = k + 1
               if (k .le. ner) then
                  ac(is,k) = a(ii)
                  jac(is,k) = ja(ii)
               endif 
            endif
 30      continue
         if (k.ge.ncmax) ncmax = k
 10   continue
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     
c     Perform some simple error checks:
c     
check maximum number of nonzeros in each row:
      if (ncmax.eq.ner) ierr = 0
      if (ncmax.gt.ner) then
         ierr = -1
         return
      endif
c     
check if there are any zero columns in AC:
c     
      do 45 in = 1,ncmax
         icount = 0
         do 44 inn =1,n
            if (ac(inn,in).ne.0.0d0) icount = 1
 44      continue
         if (icount.eq.0) then
            ierr = 1
            return
         endif
 45   continue
c     
check if there are any zero rows in AC:
c     
      do 55 inn = 1,n
         icount = 0
         do 54 in =1,ncmax
            if (ac(inn,in).ne.0.0d0) icount = 1
 54      continue
         if (icount.eq.0) then
            ierr = 2
            return
         endif
 55   continue
      return
c------------- end of xcooell ------------------------------------------- 
      end
c----------------------------------------------------------------------- 
      subroutine ccsruss (nrow,a,ja,ia,diag,al,jal,ial,au,jau,iau) 
      complex*8 a(*),al(*),diag(*),au(*) 
      integer nrow,ja(*),ia(nrow+1),jal(*),ial(nrow+1),jau(*),
     *     iau(nrow+1)
c-----------------------------------------------------------------------
c Compressed Sparse Row     to     Unsymmetric Sparse Skyline format
c----------------------------------------------------------------------- 
c this subroutine converts a matrix stored in csr format into a nonsym. 
c sparse skyline format. This latter format does not assume
c that the matrix has a symmetric pattern and consists of the following 
c * the diagonal of A stored separately in diag(*);
c * The strict lower part of A is stored  in CSR format in al,jal,ial 
c * The strict upper part is stored in CSC format in au,jau,iau.
c----------------------------------------------------------------------- 
c On entry
c---------
c nrow  = dimension of the matrix a.
c a     = real array containing the nonzero values of the matrix 
c         stored rowwise.
c ja    = column indices of the values in array a
c ia    = integer array of length n+1 containing the pointers to
c         beginning of each row in arrays a, ja.
c 
c On return
c----------
c diag  = array containing the diagonal entries of A
c al,jal,ial = matrix in CSR format storing the strict lower 
c              trangular part of A.
c au,jau,iau = matrix in CSC format storing the strict upper
c              triangular part of A. 
c----------------------------------------------------------------------- 
      integer i, j, k, kl, ku 
c
c determine U's data structure first
c 
      do 1 i=1,nrow+1
         iau(i) = 0
 1    continue
      do 3 i=1, nrow
         do 2 k=ia(i), ia(i+1)-1 
            j = ja(k)
            if (j .gt. i) iau(j+1) = iau(j+1)+1
 2       continue 
 3    continue
c
c     compute pointers from lengths
c
      iau(1) = 1
      do 4 i=1,nrow
         iau(i+1) = iau(i)+iau(i+1)
         ial(i+1) = ial(i)+ial(i+1)
 4    continue
c
c     now do the extractions. scan all rows.
c
      kl = 1
      ial(1) = kl
      do  7 i=1, nrow
c
c     scan all elements in a row
c 
         do 71 k = ia(i), ia(i+1)-1
            j = ja(k) 
c
c     if in upper part, store in row j (of transp(U) )
c     
            if (j  .gt. i) then
               ku = iau(j) 
               au(ku) = a(k)
               jau(ku) = i
               iau(j) = ku+1
            elseif (j  .eq. i) then
               diag(i) = a(k) 
            elseif (j .lt. i) then
               al(kl) = a(k)
               jal(kl) = j
               kl = kl+1
            endif
 71      continue
         ial(i+1) = kl 
 7    continue
c
c readjust iau
c
      do 8 i=nrow,1,-1
         iau(i+1) = iau(i)
 8    continue
      iau(1) = 1
c--------------- end-of-csruss ----------------------------------------- 
c-----------------------------------------------------------------------
      end 
c-----------------------------------------------------------------------
      subroutine cusscsr (nrow,a,ja,ia,diag,al,jal,ial,au,jau,iau) 
      complex*8 a(*),al(*),diag(*),au(*) 
      integer ja(*),ia(nrow+1),jal(*),ial(nrow+1),jau(*),iau(nrow+1)
c-----------------------------------------------------------------------
c Unsymmetric Sparse Skyline   format   to Compressed Sparse Row 
c----------------------------------------------------------------------- 
c this subroutine converts a matrix stored in nonsymmetric sparse
c skyline format into csr format. The sparse skyline format is 
c described in routine csruss. 
c----------------------------------------------------------------------- 
c----------------------------------------------------------------------- 
c On entry
c-----------------------------------------------------------------------
c nrow  = dimension of the matrix a.
c diag  = array containing the diagonal entries of A
c al,jal,ial = matrix in CSR format storing the strict lower 
c              trangular part of A.
c au,jau,iau = matrix in CSC format storing the strict upper
c              trangular part of A.
c On return
c --------- 
c a     = real array containing the nonzero values of the matrix 
c         stored rowwise.
c ja    = column indices of the values in array a
c ia    = integer array of length n+1 containing the pointers to
c         beginning of each row in arrays a, ja.
c 
c-----------------------------------------------------------------------
c
c count elements in lower part + diagonal 
c 
      do 1 i=1, nrow
         ia(i+1) = ial(i+1)-ial(i)+1
 1    continue
c
c count elements in upper part
c 
      do 3 i=1, nrow
         do 2 k=iau(i), iau(i+1)-1 
            j = jau(k)
            ia(j+1) = ia(j+1)+1
 2       continue 
 3    continue
c---------- compute pointers from lengths ------------------------------
      ia(1) = 1
      do 4 i=1,nrow
         ia(i+1) = ia(i)+ia(i+1)
 4    continue
c
c copy lower part + diagonal 
c 
      do 6 i=1, nrow
         ka = ia(i) 
         do 5 k=ial(i), ial(i+1)-1
            a(ka) = al(k) 
            ja(ka) = jal(k) 
            ka = ka+1
 5       continue
         a(ka) = diag(i) 
         ja(ka) = i
         ia(i) = ka+1
 6    continue
c     
c     copy upper part
c     
      do 8 i=1, nrow
         do 7 k=iau(i), iau(i+1)-1
c
c row number
c
            jak = jau(k) 
c
c where element goes
c
            ka = ia(jak) 
            a(ka) = au(k) 
            ja(ka) = i
            ia(jak) = ka+1
 7       continue
 8    continue
c
c readjust ia
c
      do 9 i=nrow,1,-1
         ia(i+1) = ia(i)
 9    continue
      ia(1) = 1
c----------end-of-usscsr------------------------------------------------
      end 
c----------------------------------------------------------------------- 
      subroutine ccsrsss (nrow,a,ja,ia,sorted,diag,al,jal,ial,au)
      complex*8 a(*),al(*),diag(*),au(*) 
      integer ja(*),ia(nrow+1),jal(*),ial(nrow+1)
      logical sorted 
c-----------------------------------------------------------------------
c Compressed Sparse Row     to     Symmetric Sparse Skyline   format 
c----------------------------------------------------------------------- 
c this subroutine converts a matrix stored in csr format into the 
c Symmetric sparse skyline   format. This latter format assumes that 
c that the matrix has a symmetric pattern. It consists of the following 
c * the diagonal of A stored separately in diag(*);
c * The strict lower part of A is stored  in csr format in al,jal,ial 
c * The values only of strict upper part as stored in csc format in au. 
c----------------------------------------------------------------------- 
c On entry
c-----------
c nrow  = dimension of the matrix a.
c a     = real array containing the nonzero values of the matrix 
c         stored rowwise.
c ja    = column indices of the values in array a
c ia    = integer array of length n+1 containing the pointers to
c         beginning of each row in arrays a, ja.
c sorted= a logical indicating whether or not the elements in a,ja,ia
c         are sorted. 
c 
c On return
c --------- 
c diag  = array containing the diagonal entries of A
c al,jal,ial = matrix in csr format storing the strict lower 
c              trangular part of A.
c au    = values of the strict upper trangular part of A, column wise.
c----------------------------------------------------------------------- 
c 
c     extract lower part and diagonal.
c
      kl = 1
      ial(1) = kl
      do  7 i=1, nrow
c
c scan all elements in a row
c 
         do 71 k = ia(i), ia(i+1)-1
            jak = ja(k) 
            if (jak  .eq. i) then
               diag(i) = a(k) 
            elseif (jak .lt. i) then
               al(kl) = a(k)
               jal(kl) = jak
               kl = kl+1
            endif
 71      continue
         ial(i+1) = kl 
 7    continue
c
c sort if not sorted
c 
      if (.not. sorted) then
         call ccsort (nrow, al, jal, ial, au, .true.) 
      endif
c
c copy u
c 
      do  8 i=1, nrow
c
c scan all elements in a row
c 
         do 81 k = ia(i), ia(i+1)-1
            jak = ja(k) 
            if (jak  .gt. i) then
               ku = ial(jak) 
               au(ku) = a(k)
               ial(jak) = ku+1
            endif
 81      continue
 8    continue
c   
c readjust ial
c
      do 9 i=nrow,1,-1
         ial(i+1) = ial(i)
 9    continue
      ial(1) = 1
c--------------- end-of-csrsss ----------------------------------------- 
c-----------------------------------------------------------------------
      end 
c
      subroutine cssscsr (nrow,a,ja,ia,diag,al,jal,ial,au) 
      complex*8 a(*),al(*),diag(*),au(*) 
      integer ja(*),ia(nrow+1),jal(*),ial(nrow+1) 
c-----------------------------------------------------------------------
c Unsymmetric Sparse Skyline   format   to Compressed Sparse Row 
c----------------------------------------------------------------------- 
c this subroutine converts a matrix stored in nonsymmetric sparse 
c skyline format into csr format. The sparse skyline format is 
c described in routine csruss. 
c----------------------------------------------------------------------- 
c On entry
c--------- 
c diag  = array containing the diagonal entries of A
c al,jal,ial = matrix in csr format storing the strict lower 
c              trangular part of A.
c au    = values of strict upper part. 
c
c On return
c --------- 
c nrow  = dimension of the matrix a.
c a     = real array containing the nonzero values of the matrix 
c         stored rowwise.
c ja    = column indices of the values in array a
c ia    = integer array of length n+1 containing the pointers to
c         beginning of each row in arrays a, ja.
c 
c-----------------------------------------------------------------------
c
c count elements in lower part + diagonal 
c 
      do 1 i=1, nrow
         ia(i+1) = ial(i+1)-ial(i)+1
 1    continue
c
c count elements in upper part
c 
      do 3 i=1, nrow
         do 2 k=ial(i), ial(i+1)-1 
            j = jal(k)
            ia(j+1) = ia(j+1)+1
 2       continue 
 3    continue
c---------- compute pointers from lengths ------------------------------
      ia(1) = 1
      do 4 i=1,nrow
         ia(i+1) = ia(i)+ia(i+1)
 4    continue
c
c copy lower part + diagonal 
c 
      do 6 i=1, nrow
         ka = ia(i) 
         do 5 k=ial(i), ial(i+1)-1
            a(ka) = al(k) 
            ja(ka) = jal(k) 
            ka = ka+1
 5       continue
         a(ka) = diag(i) 
         ia(i) = ka+1
 6    continue
c     
c     copy upper part
c     
      do 8 i=1, nrow
         do 7 k=ial(i), ial(i+1)-1
c
c row number
c
            jak = jal(k) 
c
c where element goes
c
            ka = ia(jak) 
            a(ka) = au(k) 
            ja(ka) = i
            ia(jak) = ka+1
 7       continue
 8    continue
c
c readjust ia
c
      do 9 i=nrow,1,-1
         ia(i+1) = ia(i)
 9    continue
      ia(1) = 1
c----------end-of-ssscsr------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine ccsrvbr (n,ia,ja,a,nr,nc,kvstr,kvstc,ib,jb,kb,
     &     b, job, iwk, nkmax, nzmax, ierr )
c-----------------------------------------------------------------------
      integer n, ia(n+1), ja(*), nr, nc, ib(*), jb(nkmax-1), kb(nkmax)
      integer kvstr(*), kvstc(*), job, iwk(*), nkmax, nzmax, ierr
      complex*8  a(*), b(nzmax)
c-----------------------------------------------------------------------
c     Converts compressed sparse row to variable block row format.
c-----------------------------------------------------------------------
c     On entry:
c--------------
c     n       = number of matrix rows
c     ia,ja,a = input matrix in CSR format
c
c     job     = job indicator.
c               If job=0, kvstr and kvstc are used as supplied.
c               If job=1, kvstr and kvstc are determined by the code.
c               If job=2, a conformal row/col partitioning is found and
c               returned in both kvstr and kvstc.  In the latter two cases,
c               an optimized algorithm can be used to perform the
c               conversion because all blocks are full.
c
c     nkmax   = size of supplied jb and kb arrays
c     nzmax   = size of supplied b array
c
c     If job=0 then the following are input:
c     nr,nc   = matrix block row and block column dimension
c     kvstr   = first row number for each block row
c     kvstc   = first column number for each block column.
c               (kvstr and kvstc may be the same array)
c
c     On return:
c---------------
c
c     ib,jb,kb,b = output matrix in VBR format
c
c     ierr    = error message
c               ierr = 0 means normal return
c               ierr = 1 out of space in jb and/or kb arrays
c               ierr = 2 out of space in b array
c               ierr = 3 nonsquare matrix used with job=2
c
c     If job=1,2 then the following are output:
c     nr,nc   = matrix block row and block column dimension
c     kvstr   = first row number for each block row
c     kvstc   = first column number for each block column
c               If job=2, then kvstr and kvstc contain the same info.
c
c     Work space:
c----------------
c     iwk(1:ncol) = inverse kvstc array.  If job=1,2 then we also need:
c     iwk(ncol+1:ncol+nr) = used to help determine sparsity of each block row.
c     The workspace is not assumed to be initialized to zero, nor is it
c     left that way.
c
c     Algorithms:
c----------------
c     There are two conversion codes in this routine.  The first assumes
c     that all blocks are full (there is a nonzero in the CSR data
c     structure for each entry in the block), and is used if the routine
c     determines the block partitioning itself.  The second code makes
c     no assumptions about the block partitioning, and is used if the
c     caller provides the partitioning.  The second code is much less
c     efficient than the first code.
c
c     In the first code, the CSR data structure is traversed sequentially
c     and entries are placed into the VBR data structure with stride
c     equal to the row dimension of the block row.  The columns of the
c     CSR data structure are sorted first if necessary.
c
c     In the second code, the block sparsity pattern is first determined.
c     This is done by traversing the CSR data structure and using an
c     implied linked list to determine which blocks are nonzero.  Then
c     the VBR data structure is filled by mapping each individual entry
c     in the CSR data structure into the VBR data structure.  The columns
c     of the CSR data structure are sorted first if necessary.
c
c-----------------------------------------------------------------------
c     Local variables:
c---------------------
      integer ncol, nb, neqr, numc, a0, b0, b1, k0, i, ii, j, jj, jnew
      logical sorted
c
c     ncol = number of scalar columns in matrix
c     nb = number of blocks in conformal row/col partitioning
c     neqr = number of rows in block row
c     numc = number of nonzero columns in row
c     a0 = index for entries in CSR a array
c     b0 = index for entries in VBR b array
c     b1 = temp
c     k0 = index for entries in VBR kb array
c     i  = loop index for block rows
c     ii = loop index for scalar rows in block row
c     j  = loop index for block columns
c     jj = loop index for scalar columns in block column
c     jnew = block column number
c     sorted = used to indicate if matrix already sorted by columns
c
c-----------------------------------------------------------------------
      ierr = 0
c-----sort matrix by column indices
      call csorted (n, ia, ja, sorted)
      if (.not. sorted) then
         call ccsort (n, a, ja, ia, b, .true.)
      endif
      if (job .eq. 1 .or. job .eq. 2) then
c--------need to zero workspace; first find ncol
         ncol = 0
         do i = 2, n
            ncol = max0(ncol, ja(ia(i)-1))
         enddo
         do i = 1, ncol
            iwk(i) = 0
         enddo
         call csrkvstr (n, ia, ja, nr, kvstr)
         call csrkvstc (n, ia, ja, nc, kvstc, iwk)
      endif
c-----check if want conformal partitioning
      if (job .eq. 2) then
         if (kvstr(nr+1) .ne. kvstc(nc+1)) then
            ierr = 3
            return
         endif
c        use iwk temporarily
         call kvstmerge (nr, kvstr, nc, kvstc, nb, iwk)
         nr = nb
         nc = nb
         do i = 1, nb+1
            kvstr(i) = iwk(i)
            kvstc(i) = iwk(i)
         enddo
      endif
c-----------------------------------------------------------------------
c     inverse kvst (scalar col number) = block col number
c     stored in iwk(1:n)
c-----------------------------------------------------------------------
      do i = 1, nc
         do j = kvstc(i), kvstc(i+1)-1
            iwk(j) = i
         enddo
      enddo
      ncol = kvstc(nc+1)-1
c-----jump to conversion routine
      if (job .eq. 0) goto 400
c-----------------------------------------------------------------------
c     Fast conversion for computed block partitioning
c-----------------------------------------------------------------------
      a0 = 1
      b0 = 1
      k0 = 1
      kb(1) = 1
c-----loop on block rows
      do i = 1, nr
         neqr = kvstr(i+1) - kvstr(i)
         numc = ia(kvstr(i)+1) - ia(kvstr(i))
         ib(i) = k0
c--------loop on first row in block row to determine block sparsity
         j = 0
         do jj = ia(kvstr(i)), ia(kvstr(i)+1)-1
            jnew = iwk(ja(jj))
            if (jnew .ne. j) then
c--------------check there is enough space in kb and jb arrays
               if (k0+1 .gt. nkmax) then
                  ierr = 1
                  write (*,*) 'csrvbr: no space in kb for block row ', i
                  return
               endif
c--------------set entries for this block
               j = jnew
               b0 = b0 + neqr * (kvstc(j+1) - kvstc(j))
               kb(k0+1) = b0
               jb(k0) = j
               k0 = k0 + 1
            endif
         enddo
c--------loop on scalar rows in block row
         do ii = 0, neqr-1
            b1 = kb(ib(i))+ii
c-----------loop on elements in a scalar row
            do jj = 1, numc
c--------------check there is enough space in b array
               if (b1 .gt. nzmax) then
                  ierr = 2
                  write (*,*) 'csrvbr: no space in b for block row ', i
                  return
               endif
               b(b1) = a(a0)
               b1 = b1 + neqr
               a0 = a0 + 1
            enddo
         enddo
      enddo
      ib(nr+1) = k0
      return
c-----------------------------------------------------------------------
c     Conversion for user supplied block partitioning
c-----------------------------------------------------------------------
 400  continue
c-----initialize workspace for sparsity indicator
      do i = ncol+1, ncol+nc
         iwk(i) = 0
      enddo
      k0 = 1
      kb(1) = 1
c-----find sparsity of block rows
      do i = 1, nr
         neqr = kvstr(i+1) - kvstr(i)
         numc = ia(kvstr(i)+1) - ia(kvstr(i))
         ib(i) = k0
c--------loop on all the elements in the block row to determine block sparsity
         do jj = ia(kvstr(i)), ia(kvstr(i+1))-1
            iwk(iwk(ja(jj))+ncol) = 1
         enddo
c--------use sparsity to set jb and kb arrays
         do j = 1, nc
            if (iwk(j+ncol) .ne. 0) then
c--------------check there is enough space in kb and jb arrays
               if (k0+1 .gt. nkmax) then
                  ierr = 1
                  write (*,*) 'csrvbr: no space in kb for block row ', i
                  return
               endif
               kb(k0+1) = kb(k0) + neqr * (kvstc(j+1) - kvstc(j))
               jb(k0) = j
               k0 = k0 + 1
               iwk(j+ncol) = 0
            endif
         enddo
      enddo
      ib(nr+1) = k0
c-----Fill b with entries from a by traversing VBR data structure.
      a0 = 1
c-----loop on block rows
      do i = 1, nr
         neqr = kvstr(i+1) - kvstr(i)
c--------loop on scalar rows in block row
         do ii = 0, neqr-1
            b0 = kb(ib(i)) + ii
c-----------loop on block columns
            do j = ib(i), ib(i+1)-1
c--------------loop on scalar columns within block column
               do jj = kvstc(jb(j)), kvstc(jb(j)+1)-1
c-----------------check there is enough space in b array
                  if (b0 .gt. nzmax) then
                     ierr = 2
                     write (*,*)'csrvbr: no space in b for blk row',i
                     return
                  endif
                  if (a0 .ge. ia(kvstr(i)+ii+1)) then
                     b(b0) = 0.d0
                  else
                     if (jj .eq. ja(a0)) then
                        b(b0) = a(a0)
                        a0 = a0 + 1
                     else
                        b(b0) = 0.d0
                     endif
                  endif
                  b0 = b0 + neqr
c--------------endloop on scalar columns
               enddo
c-----------endloop on block columns
            enddo
 2020       continue
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
c----------------------------end-of-csrvbr------------------------------
c----------------------------------------------------------------------c
      subroutine cvbrcsr (ia, ja, a, nr, kvstr, kvstc, ib, jb, kb,
     &   b, nzmax, ierr)
c-----------------------------------------------------------------------
      integer ia(*), ja(*), nr, ib(nr+1), jb(*), kb(*)
      integer kvstr(nr+1), kvstc(*), nzmax, ierr
      complex*8  a(*), b(nzmax)
c-----------------------------------------------------------------------
c     Converts variable block row to compressed sparse row format.
c-----------------------------------------------------------------------
c     On entry:
c--------------
c     nr      = number of block rows
c     kvstr   = first row number for each block row
c     kvstc   = first column number for each block column
c     ib,jb,kb,b = input matrix in VBR format
c     nzmax   = size of supplied ja and a arrays
c
c     On return:
c---------------
c     ia,ja,a = output matrix in CSR format
c
c     ierr    = error message
c               ierr = 0 means normal return
c               ierr = negative row number when out of space in
c                      ja and a arrays
c
c     Work space:
c----------------
c     None
c
c     Algorithm:
c---------------
c     The VBR data structure is traversed in the order that is required
c     to fill the CSR data structure.  In a given block row, consecutive
c     entries in the CSR data structure are entries in the VBR data
c     structure with stride equal to the row dimension of the block.
c     The VBR data structure is assumed to be sorted by block columns.
c
c-----------------------------------------------------------------------
c     Local variables:
c---------------------
      integer neqr, numc, a0, b0, i, ii, j, jj
c
c     neqr = number of rows in block row
c     numc = number of nonzero columns in row
c     a0 = index for entries in CSR a array
c     b0 = index for entries in VBR b array
c     i  = loop index for block rows
c     ii = loop index for scalar rows in block row
c     j  = loop index for block columns
c     jj = loop index for scalar columns in block column
c
c-----------------------------------------------------------------------
      ierr = 0
      a0 = 1
      b0 = 1
c-----loop on block rows
      do i = 1, nr
c--------set num of rows in block row, and num of nonzero cols in row
         neqr = kvstr(i+1) - kvstr(i)
         numc = ( kb(ib(i+1)) - kb(ib(i)) ) / neqr
c--------construct ja for a scalar row
         do j = ib(i), ib(i+1)-1
            do jj = kvstc(jb(j)), kvstc(jb(j)+1)-1
               ja(a0) = jj
               a0 = a0 + 1
            enddo
         enddo
c--------construct neqr-1 additional copies of ja for the block row
         do ii = 1, neqr-1
            do j = 1, numc
               ja(a0) = ja(a0-numc)
               a0 = a0 + 1
            enddo
         enddo
c--------reset a0 back to beginning of block row
         a0 = kb(ib(i))
c--------loop on scalar rows in block row
         do ii = 0, neqr-1
            ia(kvstr(i)+ii) = a0
            b0 = kb(ib(i)) + ii
c-----------loop on elements in a scalar row
            do jj = 1, numc
c--------------check there is enough space in a array
               if (a0 .gt. nzmax) then
                  ierr = -(kvstr(i)+ii)
                  write (*,*) 'vbrcsr: no space for row ', -ierr
                  return
               endif
               a(a0) = b(b0)
               a0 = a0 + 1
               b0 = b0 + neqr
            enddo
         enddo
c-----endloop on block rows
      enddo
      ia(kvstr(nr+1)) = a0
      return
      end
c-----------------------------------------------------------------------
c---------------------------end-of-vbrcsr-------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
