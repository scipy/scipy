c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c        BASIC LINEAR ALGEBRA FOR SPARSE MATRICES. BLASSM MODULE       c
c----------------------------------------------------------------------c
c amub   :   computes     C = A*B                                      c
c aplb   :   computes     C = A+B                                      c
c aplb1  :   computes     C = A+B  [Sorted version: A, B, C sorted]    c
c aplsb  :   computes     C = A + s B                                  c
c aplsb1 :   computes     C = A+sB  [Sorted version: A, B, C sorted]   c
c apmbt  :   Computes     C = A +/- transp(B)                          c
c aplsbt :   Computes     C = A + s * transp(B)                        c
c diamua :   Computes     C = Diag * A                                 c
c amudia :   Computes     C = A* Diag                                  c
c aplsca :   Computes     A:= A + s I    (s = scalar)                  c 
c apldia :   Computes     C = A + Diag.                                c
c----------------------------------------------------------------------c 
c Note: this module still incomplete.                                  c
c----------------------------------------------------------------------c
       subroutine damub (nrow,ncol,job,a,ja,ia,b,jb,ib,
     *                  c,jc,ic,nzmax,iw,ierr) 
      real*8 a(*), b(*), c(*) 
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(*),ic(*),iw(ncol)
c-----------------------------------------------------------------------
c performs the matrix by matrix product C = A B 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow  = integer. The row dimension of A = row dimension of C
c ncol  = integer. The column dimension of B = column dimension of C
c job   = integer. Job indicator. When job = 0, only the structure
c                  (i.e. the arrays jc, ic) is computed and the
c                  real values are ignored.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c b, 
c jb, 
c ib    =  Matrix B in compressed sparse row format.
c
c nzmax = integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number 
c         of elements that exceeds exceeds nzmax. See ierr.
c 
c on return:
c----------
c c, 
c jc, 
c ic    = resulting matrix C in compressed sparse row sparse format.
c           
c ierr  = integer. serving as error message. 
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number 
c         of elements in C exceeds nzmax.
c
c work arrays:
c------------
c iw    = integer work array of length equal to the number of
c         columns in A.
c Note: 
c-------
c   The row dimension of B is not needed. However there is no checking 
c   on the condition that ncol(A) = nrow(B). 
c
c----------------------------------------------------------------------- 
      real*8 scal 
      logical values
      values = (job .ne. 0) 
      len = 0
      ic(1) = 1 
      ierr = 0
c     initialize array iw.
      do 1 j=1, ncol
         iw(j) = 0
 1    continue
c
      do 500 ii=1, nrow 
c     row i 
         do 200 ka=ia(ii), ia(ii+1)-1 
	    if (values) scal = a(ka)
	    jj   = ja(ka)
	    do 100 kb=ib(jj),ib(jj+1)-1
               jcol = jb(kb)
               jpos = iw(jcol)
               if (jpos .eq. 0) then
                  len = len+1
                  if (len .gt. nzmax) then
                     ierr = ii
                     return
                  endif
                  jc(len) = jcol
                  iw(jcol)= len
                  if (values) c(len)  = scal*b(kb)
               else
                  if (values) c(jpos) = c(jpos) + scal*b(kb)
               endif
 100	    continue
 200     continue
         do 201 k=ic(ii), len
	    iw(jc(k)) = 0
 201     continue
         ic(ii+1) = len+1
 500  continue
      return
c-------------end-of-amub-----------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine daplb (nrow,ncol,job,a,ja,ia,b,jb,ib,
     *     c,jc,ic,nzmax,iw,ierr)
      real*8 a(*), b(*), c(*) 
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1),
     *     iw(ncol)
c-----------------------------------------------------------------------
c performs the matrix sum  C = A+B. 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A and B
c ncol  = integer. The column dimension of A and B.
c job   = integer. Job indicator. When job = 0, only the structure
c                  (i.e. the arrays jc, ic) is computed and the
c                  real values are ignored.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c b, 
c jb, 
c ib	=  Matrix B in compressed sparse row format.
c
c nzmax	= integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number 
c         of elements that exceeds exceeds nzmax. See ierr.
c 
c on return:
c----------
c c, 
c jc, 
c ic	= resulting matrix C in compressed sparse row sparse format.
c	    
c ierr	= integer. serving as error message. 
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number 
c         of elements in C exceeds nzmax.
c
c work arrays:
c------------
c iw	= integer work array of length equal to the number of
c         columns in A.
c
c-----------------------------------------------------------------------
      logical values
      values = (job .ne. 0) 
      ierr = 0
      len = 0
      ic(1) = 1 
      do 1 j=1, ncol
         iw(j) = 0
 1    continue
c     
      do 500 ii=1, nrow
c     row i 
         do 200 ka=ia(ii), ia(ii+1)-1 
            len = len+1
            jcol    = ja(ka)
            if (len .gt. nzmax) goto 999
            jc(len) = jcol 
            if (values) c(len)  = a(ka) 
            iw(jcol)= len
 200     continue
c     
         do 300 kb=ib(ii),ib(ii+1)-1
            jcol = jb(kb)
            jpos = iw(jcol)
            if (jpos .eq. 0) then
               len = len+1
               if (len .gt. nzmax) goto 999
               jc(len) = jcol
               if (values) c(len)  = b(kb)
               iw(jcol)= len
            else
               if (values) c(jpos) = c(jpos) + b(kb)
            endif
 300     continue
         do 301 k=ic(ii), len
	    iw(jc(k)) = 0
 301     continue
         ic(ii+1) = len+1
 500  continue
      return
 999  ierr = ii
      return
c------------end of aplb ----------------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine daplb1 (nrow,ncol,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,
     *                   ierr)
      real*8 a(*), b(*), c(*) 
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1)
c-----------------------------------------------------------------------
c performs the matrix sum  C = A+B for matrices in sorted CSR format.
c the difference with aplb  is that the resulting matrix is such that
c the elements of each row are sorted with increasing column indices in
c each row, provided the original matrices are sorted in the same way. 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A and B
c ncol  = integer. The column dimension of A and B.
c job   = integer. Job indicator. When job = 0, only the structure
c                  (i.e. the arrays jc, ic) is computed and the
c                  real values are ignored.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format with entries sorted
c 
c b, 
c jb, 
c ib	=  Matrix B in compressed sparse row format with entries sorted
c        ascendly in each row   
c
c nzmax	= integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number 
c         of elements that exceeds exceeds nzmax. See ierr.
c 
c on return:
c----------
c c, 
c jc, 
c ic	= resulting matrix C in compressed sparse row sparse format
c         with entries sorted ascendly in each row. 
c	    
c ierr	= integer. serving as error message. 
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number 
c         of elements in C exceeds nzmax.
c
c Notes: 
c-------
c     this will not work if any of the two input matrices is not sorted
c-----------------------------------------------------------------------
      logical values
      values = (job .ne. 0) 
      ierr = 0
      kc = 1
      ic(1) = kc 
c
      do 6 i=1, nrow
         ka = ia(i)
         kb = ib(i)
         kamax = ia(i+1)-1
         kbmax = ib(i+1)-1 
 5       continue 
         if (ka .le. kamax) then
            j1 = ja(ka)
         else
            j1 = ncol+1
         endif
         if (kb .le. kbmax) then 
            j2 = jb(kb)         
         else 
            j2 = ncol+1
         endif
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
         if (kc .gt. nzmax) goto 999
         if (ka .le. kamax .or. kb .le. kbmax) goto 5
         ic(i+1) = kc
 6    continue
      return
 999  ierr = i 
      return
c------------end-of-aplb1----------------------------------------------- 
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine daplsb (nrow,ncol,a,ja,ia,s,b,jb,ib,c,jc,ic,
     *     nzmax,ierr)
      real*8 a(*), b(*), c(*), s
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1)
c-----------------------------------------------------------------------
c performs the operation C = A+s B for matrices in sorted CSR format.
c the difference with aplsb is that the resulting matrix is such that
c the elements of each row are sorted with increasing column indices in
c each row, provided the original matrices are sorted in the same way. 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A and B
c ncol  = integer. The column dimension of A and B.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format with entries sorted
c
c s	= real. scalar factor for B.
c 
c b, 
c jb, 
c ib	=  Matrix B in compressed sparse row format with entries sorted
c        ascendly in each row   
c
c nzmax	= integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number 
c         of elements that exceeds exceeds nzmax. See ierr.
c 
c on return:
c----------
c c, 
c jc, 
c ic	= resulting matrix C in compressed sparse row sparse format
c         with entries sorted ascendly in each row. 
c	    
c ierr	= integer. serving as error message. 
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number 
c         of elements in C exceeds nzmax.
c
c Notes: 
c-------
c     this will not work if any of the two input matrices is not sorted
c-----------------------------------------------------------------------
      ierr = 0
      kc = 1
      ic(1) = kc 
c
c     the following loop does a merge of two sparse rows + adds  them.
c 
      do 6 i=1, nrow
         ka = ia(i)
         kb = ib(i)
         kamax = ia(i+1)-1
         kbmax = ib(i+1)-1 
 5       continue 
c
c     this is a while  -- do loop -- 
c 
         if (ka .le. kamax .or. kb .le. kbmax) then 
c     
            if (ka .le. kamax) then
               j1 = ja(ka)
            else
c     take j1 large enough  that always j2 .lt. j1
               j1 = ncol+1
            endif
            if (kb .le. kbmax) then 
               j2 = jb(kb)         
            else 
c     similarly take j2 large enough  that always j1 .lt. j2 
               j2 = ncol+1
            endif
c     
c     three cases
c     
            if (j1 .eq. j2) then 
               c(kc) = a(ka)+s*b(kb)
               jc(kc) = j1
               ka = ka+1
               kb = kb+1
               kc = kc+1
            else if (j1 .lt. j2) then
               jc(kc) = j1
               c(kc) = a(ka)
               ka = ka+1
               kc = kc+1
            else if (j1 .gt. j2) then
               jc(kc) = j2
               c(kc) = s*b(kb)
               kb = kb+1
               kc = kc+1
            endif
            if (kc .gt. nzmax) goto 999
            goto 5
c
c     end while loop
c
         endif
         ic(i+1) = kc
 6    continue
      return
 999  ierr = i 
      return
c------------end-of-aplsb --------------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine daplsb1 (nrow,ncol,a,ja,ia,s,b,jb,ib,c,jc,ic,
     *     nzmax,ierr)
      real*8 a(*), b(*), c(*), s
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1)
c-----------------------------------------------------------------------
c performs the operation C = A+s B for matrices in sorted CSR format.
c the difference with aplsb is that the resulting matrix is such that
c the elements of each row are sorted with increasing column indices in
c each row, provided the original matrices are sorted in the same way. 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A and B
c ncol  = integer. The column dimension of A and B.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format with entries sorted
c
c s	= real. scalar factor for B.
c 
c b, 
c jb, 
c ib	=  Matrix B in compressed sparse row format with entries sorted
c        ascendly in each row   
c
c nzmax	= integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number 
c         of elements that exceeds exceeds nzmax. See ierr.
c 
c on return:
c----------
c c, 
c jc, 
c ic	= resulting matrix C in compressed sparse row sparse format
c         with entries sorted ascendly in each row. 
c	    
c ierr	= integer. serving as error message. 
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number 
c         of elements in C exceeds nzmax.
c
c Notes: 
c-------
c     this will not work if any of the two input matrices is not sorted
c-----------------------------------------------------------------------
      ierr = 0
      kc = 1
      ic(1) = kc 
c
c     the following loop does a merge of two sparse rows + adds  them.
c 
      do 6 i=1, nrow
         ka = ia(i)
         kb = ib(i)
         kamax = ia(i+1)-1
         kbmax = ib(i+1)-1 
 5       continue 
c
c     this is a while  -- do loop -- 
c 
         if (ka .le. kamax .or. kb .le. kbmax) then 
c     
            if (ka .le. kamax) then
               j1 = ja(ka)
            else
c     take j1 large enough  that always j2 .lt. j1
               j1 = ncol+1
            endif
            if (kb .le. kbmax) then 
               j2 = jb(kb)         
            else 
c     similarly take j2 large enough  that always j1 .lt. j2 
               j2 = ncol+1
            endif
c     
c     three cases
c     
            if (j1 .eq. j2) then 
               c(kc) = a(ka)+s*b(kb)
               jc(kc) = j1
               ka = ka+1
               kb = kb+1
               kc = kc+1
            else if (j1 .lt. j2) then
               jc(kc) = j1
               c(kc) = a(ka)
               ka = ka+1
               kc = kc+1
            else if (j1 .gt. j2) then
               jc(kc) = j2
               c(kc) = s*b(kb)
               kb = kb+1
               kc = kc+1
            endif
            if (kc .gt. nzmax) goto 999
            goto 5
c
c     end while loop
c
         endif
         ic(i+1) = kc
 6    continue
      return
 999  ierr = i 
      return
c------------end-of-aplsb1 --------------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine dapmbt (nrow,ncol,job,a,ja,ia,b,jb,ib,
     *     c,jc,ic,nzmax,iw,ierr)
      real*8 a(*), b(*), c(*) 
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(ncol+1),ic(*),iw(*) 
c-----------------------------------------------------------------------
c performs the matrix sum  C = A + transp(B) or C = A - transp(B) 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A and transp(B)
c ncol  = integer. The column dimension of A. Also the row 
c                  dimension of B. 
c
c job	= integer. if job = -1, apmbt will compute C= A - transp(B)
c         (structure + values) 
c         if (job .eq. 1)  it will compute C=A+transp(A) 
c         (structure+ values) 
c         if (job .eq. 0) it will compute the structure of
c         C= A+/-transp(B) only (ignoring all real values).
c         any other value of job will be treated as  job=1
c a,
c ja,
c ia    = Matrix A in compressed sparse row format.
c
c b, 
c jb, 
c ib	=  Matrix B in compressed sparse row format.
c
c nzmax	= integer. The  length of the arrays c, jc, and ic.
c         amub will stop if the result matrix C  has a number 
c         of elements that exceeds exceeds nzmax. See ierr.
c 
c on return:
c----------
c c, 
c jc, 
c ic	= resulting matrix C in compressed sparse row format.
c	    
c ierr	= integer. serving as error message. 
c         ierr = 0 means normal return.
c         ierr = -1 means that nzmax was .lt. either the number of
c         nonzero elements of A or the number of nonzero elements in B.
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number 
c         of elements in C exceeds nzmax.
c
c work arrays:
c------------
c iw	= integer work array of length at least max(ncol,nrow) 
c
c Notes:
c------- It is important to note that here all of three arrays c, ic, 
c        and jc are assumed to be of length nnz(c). This is because 
c        the matrix is internally converted in coordinate format.
c        
c-----------------------------------------------------------------------
      logical values
      values = (job .ne. 0) 
c
      ierr = 0
      do 1 j=1, ncol
         iw(j) = 0
 1    continue
c     
      nnza = ia(nrow+1)-1
      nnzb = ib(ncol+1)-1
      len = nnzb
      if (nzmax .lt. nnzb .or. nzmax .lt. nnza) then
         ierr = -1
         return
      endif
c     
c trasnpose matrix b into c
c
      ljob = 0
      if (values) ljob = 1
      ipos = 1
      call scsrcsc (ncol,ljob,ipos,b,jb,ib,c,jc,ic) 
c----------------------------------------------------------------------- 
      if (job .eq. -1) then
         do 2 k=1,len
	    c(k) = -c(k)
 2       continue
      endif
c
c--------------- main loop --------------------------------------------
c
      do 500 ii=1, nrow
         do 200 k = ic(ii),ic(ii+1)-1
            iw(jc(k)) = k
 200     continue
c-----------------------------------------------------------------------     
         do 300 ka = ia(ii), ia(ii+1)-1 
            jcol = ja(ka)
            jpos = iw(jcol)
            if (jpos .eq. 0) then
c
c     if fill-in append in coordinate format to matrix.
c 
               len = len+1
               if (len .gt. nzmax) goto 999
               jc(len) = jcol
               
               ic(len) = ii
               if (values) c(len)  = a(ka)
            else
c     else do addition.
               if (values) c(jpos) = c(jpos) + a(ka)
            endif
 300     continue
         do 301 k=ic(ii), ic(ii+1)-1
	    iw(jc(k)) = 0
 301     continue
 500  continue
c     
c     convert first part of matrix (without fill-ins) into coo format
c     
      ljob = 2
      if (values) ljob = 3
      do 501 i=1, nrow+1
         iw(i) = ic(i) 
 501  continue
      call scsrcoo (nrow,ljob,nnzb,c,jc,iw,nnzb,c,ic,jc,ierr)
c
c     convert the whole thing back to csr format. 
c 
      ljob = 0
      if (values) ljob = 1
      call scoicsr (nrow,len,ljob,c,jc,ic,iw)
      return
 999  ierr = ii
      return
c--------end-of-apmbt---------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine daplsbt(nrow,ncol,a,ja,ia,s,b,jb,ib,
     *     c,jc,ic,nzmax,iw,ierr)
      real*8 a(*), b(*), c(*), s
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(ncol+1),ic(*),iw(*)
c-----------------------------------------------------------------------
c performs the matrix sum  C = A + transp(B).
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A and transp(B)
c ncol  = integer. The column dimension of A. Also the row 
c                  dimension of B. 
c
c a,
c ja,
c ia    = Matrix A in compressed sparse row format.
c
c s	= real. scalar factor for B.
c
c 
c b, 
c jb, 
c ib	=  Matrix B in compressed sparse row format.
c
c nzmax	= integer. The  length of the arrays c, jc, and ic.
c         amub will stop if the result matrix C  has a number 
c         of elements that exceeds exceeds nzmax. See ierr.
c 
c on return:
c----------
c c, 
c jc, 
c ic	= resulting matrix C in compressed sparse row format.
c	    
c ierr	= integer. serving as error message. 
c         ierr = 0 means normal return.
c         ierr = -1 means that nzmax was .lt. either the number of
c         nonzero elements of A or the number of nonzero elements in B.
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number 
c         of elements in C exceeds nzmax.
c
c work arrays:
c------------
c iw	= integer work array of length at least max(nrow,ncol) 
c
c Notes:
c------- It is important to note that here all of three arrays c, ic, 
c        and jc are assumed to be of length nnz(c). This is because 
c        the matrix is internally converted in coordinate format.
c        
c-----------------------------------------------------------------------
      ierr = 0
      do 1 j=1, ncol
         iw(j) = 0
 1    continue
c     
      nnza = ia(nrow+1)-1
      nnzb = ib(ncol+1)-1
      len = nnzb
      if (nzmax .lt. nnzb .or. nzmax .lt. nnza) then
         ierr = -1
         return
      endif
c     
c     transpose matrix b into c
c
      ljob = 1
      ipos = 1
      call scsrcsc (ncol,ljob,ipos,b,jb,ib,c,jc,ic) 
      do 2 k=1,len
 2       c(k) = c(k)*s
c     
c     main loop. add rows from ii = 1 to nrow.
c     
         do 500 ii=1, nrow
c     iw is used as a system to recognize whether there
c     was a nonzero element in c. 
            do 200 k = ic(ii),ic(ii+1)-1
               iw(jc(k)) = k
 200        continue
c     
            do 300 ka = ia(ii), ia(ii+1)-1 
               jcol = ja(ka)
               jpos = iw(jcol)
           if (jpos .eq. 0) then
c     
c     if fill-in append in coordinate format to matrix.
c     
              len = len+1
              if (len .gt. nzmax) goto 999
              jc(len) = jcol              
              ic(len) = ii
              c(len)  = a(ka)
           else
c     else do addition.
              c(jpos) = c(jpos) + a(ka)
           endif
 300    continue
        do 301 k=ic(ii), ic(ii+1)-1
           iw(jc(k)) = 0
 301    continue
 500  continue
c     
c     convert first part of matrix (without fill-ins) into coo format
c     
      ljob = 3
      do 501 i=1, nrow+1
         iw(i) = ic(i) 
 501  continue
      call scsrcoo (nrow,ljob,nnzb,c,jc,iw,nnzb,c,ic,jc,ierr)
c
c     convert the whole thing back to csr format. 
c 
      ljob = 1
      call scoicsr (nrow,len,ljob,c,jc,ic,iw)
      return
 999  ierr = ii
      return
c--------end-of-aplsbt--------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine ddiamua (nrow,job, a, ja, ia, diag, b, jb, ib)
      real*8 a(*), b(*), diag(nrow), scal
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1) 
c-----------------------------------------------------------------------
c performs the matrix by matrix product B = Diag * A  (in place) 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c diag = diagonal matrix stored as a vector dig(1:n)
c
c on return:
c----------
c
c b, 
c jb, 
c ib	= resulting matrix B in compressed sparse row sparse format.
c	    
c Notes:
c-------
c 1)        The column dimension of A is not needed. 
c 2)        algorithm in place (B can take the place of A).
c           in this case use job=0.
c-----------------------------------------------------------------
      do 1 ii=1,nrow
c     
c     normalize each row 
c     
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         scal = diag(ii) 
         do 2 k=k1, k2
            b(k) = a(k)*scal
 2       continue
 1    continue
c     
      if (job .eq. 0) return
c     
      do 3 ii=1, nrow+1
         ib(ii) = ia(ii)
 3    continue
      do 31 k=ia(1), ia(nrow+1) -1 
         jb(k) = ja(k)
 31   continue
      return
c----------end-of-diamua------------------------------------------------
c-----------------------------------------------------------------------
      end 
c----------------------------------------------------------------------- 
      subroutine damudia (nrow,job, a, ja, ia, diag, b, jb, ib)
      real*8 a(*), b(*), diag(nrow) 
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1) 
c-----------------------------------------------------------------------
c performs the matrix by matrix product B = A * Diag  (in place) 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c diag = diagonal matrix stored as a vector dig(1:n)
c
c on return:
c----------
c
c b, 
c jb, 
c ib	= resulting matrix B in compressed sparse row sparse format.
c	    
c Notes:
c-------
c 1)        The column dimension of A is not needed. 
c 2)        algorithm in place (B can take the place of A).
c-----------------------------------------------------------------
      do 1 ii=1,nrow
c     
c     scale each element 
c     
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         do 2 k=k1, k2
            b(k) = a(k)*diag(ja(k)) 
 2       continue
 1    continue
c     
      if (job .eq. 0) return
c     
      do 3 ii=1, nrow+1
         ib(ii) = ia(ii)
 3    continue
      do 31 k=ia(1), ia(nrow+1) -1 
         jb(k) = ja(k)
 31   continue
      return
c-----------------------------------------------------------------------
c-----------end-of-amudiag----------------------------------------------
      end 
c----------------------------------------------------------------------- 
      subroutine daplsca (nrow, a, ja, ia, scal,iw) 
      real*8 a(*), scal
      integer ja(*), ia(nrow+1),iw(*)
c-----------------------------------------------------------------------
c Adds a scalar to the diagonal entries of a sparse matrix A :=A + s I 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c a,
c ja,
c ia    = Matrix A in compressed sparse row format.
c 
c scal  = real. scalar to add to the diagonal entries. 
c
c on return:
c----------
c
c a, 
c ja, 
c ia	= matrix A with diagonal elements shifted (or created).
c	    
c iw    = integer work array of length n. On return iw will
c         contain  the positions of the diagonal entries in the 
c         output matrix. (i.e., a(iw(k)), ja(iw(k)), k=1,...n,
c         are the values/column indices of the diagonal elements 
c         of the output matrix. ). 
c
c Notes:
c-------
c     The column dimension of A is not needed. 
c     important: the matrix a may be expanded slightly to allow for
c     additions of nonzero elements to previously nonexisting diagonals.
c     The is no checking as to whether there is enough space appended
c     to the arrays a and ja. if not sure allow for n additional 
c     elemnts. 
c coded by Y. Saad. Latest version July, 19, 1990
c-----------------------------------------------------------------------
      logical test
c
      call diapos (nrow,ja,ia,iw)
      icount = 0
      do 1 j=1, nrow
         if (iw(j) .eq. 0) then
            icount = icount+1
         else
            a(iw(j)) = a(iw(j)) + scal 
         endif
 1    continue
c
c     if no diagonal elements to insert in data structure return.
c
      if (icount .eq. 0) return
c
c shift the nonzero elements if needed, to allow for created 
c diagonal elements. 
c
      ko = ia(nrow+1)+icount
c
c     copy rows backward
c
      do 5 ii=nrow, 1, -1 
c     
c     go through  row ii
c     
         k1 = ia(ii)
         k2 = ia(ii+1)-1 
         ia(ii+1) = ko
         test = (iw(ii) .eq. 0) 
         do 4 k = k2,k1,-1 
            j = ja(k)
            if (test .and. (j .lt. ii)) then 
               test = .false. 
               ko = ko - 1
               a(ko) = scal 
               ja(ko) = ii
               iw(ii) = ko
            endif
            ko = ko-1
            a(ko) = a(k) 
            ja(ko) = j
 4       continue
c     diagonal element has not been added yet.
         if (test) then
            ko = ko-1
            a(ko) = scal 
            ja(ko) = ii
            iw(ii) = ko
         endif
 5    continue
      ia(1) = ko 
      return
c-----------------------------------------------------------------------
c----------end-of-aplsca------------------------------------------------ 
      end
c-----------------------------------------------------------------------
      subroutine dapldia (nrow, job, a, ja, ia, diag, b, jb, ib, iw) 
      real*8 a(*), b(*), diag(nrow) 
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1), iw(*)
c-----------------------------------------------------------------------
c Adds a diagonal matrix to a general sparse matrix:  B = A + Diag 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         (i.e. assume that a has already been copied into array b,
c         or that algorithm is used in place. ) For all practical
c         puposes enter job=0 for an in-place call and job=1 otherwise.
c 
c         Note: in case there are missing diagonal elements in A, 
c         then the option job =0 will be ignored, since the algorithm 
c         must modify the data structure (i.e. jb, ib) in this 
c         situation.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c     
c diag = diagonal matrix stored as a vector dig(1:n)
c
c on return:
c----------
c
c b, 
c jb, 
c ib	= resulting matrix B in compressed sparse row sparse format.
c
c
c iw    = integer work array of length n. On return iw will
c         contain  the positions of the diagonal entries in the 
c         output matrix. (i.e., a(iw(k)), ja(iw(k)), k=1,...n,
c         are the values/column indices of the diagonal elements 
c         of the output matrix. ). 
c
c Notes:
c-------
c 1)        The column dimension of A is not needed. 
c 2)        algorithm in place (b, jb, ib, can be the same as
c           a, ja, ia, on entry). See comments for parameter job.
c
c coded by Y. Saad. Latest version July, 19, 1990
c-----------------------------------------------------------------
      logical test
c
c     copy integer arrays into b's data structure if required
c
      if (job .ne. 0) then 
         nnz = ia(nrow+1)-1
         do 2  k=1, nnz
            jb(k) = ja(k)
 2       continue
         do 3 k=1, nrow+1
            ib(k) = ia(k)
 3       continue
      endif 
c
c     get positions of diagonal elements in data structure.
c     
      call diapos (nrow,ja,ia,iw)
c     
c     count number of holes in diagonal and add diag(*) elements to
c     valid diagonal entries.
c     
      icount = 0
      do 1 j=1, nrow
         if (iw(j) .eq. 0) then
            icount = icount+1
         else
            b(iw(j)) = a(iw(j)) + diag(j) 
         endif
 1    continue
c     
c     if no diagonal elements to insert return
c     
      if (icount .eq. 0) return
c     
c     shift the nonzero elements if needed, to allow for created 
c     diagonal elements. 
c     
      ko = ib(nrow+1)+icount
c     
c     copy rows backward
c     
      do 5 ii=nrow, 1, -1 
c     
c     go through  row ii
c     
         k1 = ib(ii)
         k2 = ib(ii+1)-1 
         ib(ii+1) = ko
         test = (iw(ii) .eq. 0) 
         do 4 k = k2,k1,-1 
            j = jb(k)
            if (test .and. (j .lt. ii)) then 
               test = .false. 
               ko = ko - 1
               b(ko) = diag(ii) 
               jb(ko) = ii
               iw(ii) = ko
            endif
            ko = ko-1
            b(ko) = a(k) 
            jb(ko) = j
 4       continue
c     diagonal element has not been added yet.
         if (test) then
            ko = ko-1
            b(ko) =  diag(ii) 
            jb(ko) = ii
            iw(ii) = ko
         endif
 5    continue
      ib(1) = ko 
      return
c-----------------------------------------------------------------------
c------------end-of-apldiag---------------------------------------------
      end 
 
