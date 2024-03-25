c       this file contains the following user-callable routines:
c
c
c       routine iddp_id computes the ID of a matrix,
c       to a specified precision.
c
c       routine iddr_id computes the ID of a matrix,
c       to a specified rank.
c
c       routine idd_reconid reconstructs a matrix from its ID.
c
c       routine idd_copycols collects together selected columns
c       of a matrix.
c
c       routine idd_getcols collects together selected columns
c       of a matrix specified by a routine for applying the matrix
c       to arbitrary vectors.
c
c       routine idd_reconint constructs p in the ID a = b p,
c       where the columns of b are a subset of the columns of a,
c       and p is the projection coefficient matrix,
c       given list, krank, and proj output by routines iddr_id
c       or iddp_id.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine iddp_id(eps,m,n,a,krank,list,rnorms)
c
c       computes the ID of a, i.e., lists in list the indices
c       of krank columns of a such that
c
c       a(j,list(k))  =  a(j,list(k))
c
c       for all j = 1, ..., m; k = 1, ..., krank, and
c
c                        krank
c       a(j,list(k))  =  Sigma  a(j,list(l)) * proj(l,k-krank)       (*)
c                         l=1
c
c                     +  epsilon(j,k-krank)
c
c       for all j = 1, ..., m; k = krank+1, ..., n,
c
c       for some matrix epsilon dimensioned epsilon(m,n-krank)
c       such that the greatest singular value of epsilon
c       <= the greatest singular value of a * eps.
c       The present routine stores the krank x (n-krank) matrix proj
c       in the memory initially occupied by a.
c
c       input:
c       eps -- relative precision of the resulting ID
c       m -- first dimension of a
c       n -- second dimension of a, as well as the dimension required
c            of list
c       a -- matrix to be ID'd
c
c       output:
c       a -- the first krank*(n-krank) elements of a constitute
c            the krank x (n-krank) interpolation matrix proj
c       krank -- numerical rank
c       list -- list of the indices of the krank columns of a
c               through which the other columns of a are expressed;
c               also, list describes the permutation of proj
c               required to reconstruct a as indicated in (*) above
c       rnorms -- absolute values of the entries on the diagonal
c                 of the triangular matrix used to compute the ID
c                 (these may be used to check the stability of the ID)
c
c       _N.B._: This routine changes a.
c
c       reference:
c       Cheng, Gimbutas, Martinsson, Rokhlin, "On the compression of
c            low-rank matrices," SIAM Journal on Scientific Computing,
c            26 (4): 1389-1404, 2005.
c
        implicit none
        integer m,n,krank,k,list(n),iswap
        real*8 a(m,n),eps,rnorms(n)
c
c
c       QR decompose a.
c
        call iddp_qrpiv(eps,m,n,a,krank,list,rnorms)
c
c
c       Build the list of columns chosen in a
c       by multiplying together the permutations in list,
c       with the permutation swapping 1 and list(1) taken rightmost
c       in the product, that swapping 2 and list(2) taken next
c       rightmost, ..., that swapping krank and list(krank) taken
c       leftmost.
c
        do k = 1,n
          rnorms(k) = k
        enddo ! k
c
        if(krank .gt. 0) then
          do k = 1,krank
c
c           Swap rnorms(k) and rnorms(list(k)).
c
            iswap = rnorms(k)
            rnorms(k) = rnorms(list(k))
            rnorms(list(k)) = iswap
c
          enddo ! k
        endif
c
        do k = 1,n
          list(k) = rnorms(k)
        enddo ! k
c
c
c       Fill rnorms for the output.
c
        if(krank .gt. 0) then
c
          do k = 1,krank
            rnorms(k) = a(k,k)
          enddo ! k
c
        endif
c
c
c       Backsolve for proj, storing it at the beginning of a.
c
        if(krank .gt. 0) then
          call idd_lssolve(m,n,a,krank)
        endif
c
c
        return
        end
c
c
c
c
        subroutine iddr_id(m,n,a,krank,list,rnorms)
c
c       computes the ID of a, i.e., lists in list the indices
c       of krank columns of a such that
c
c       a(j,list(k))  =  a(j,list(k))
c
c       for all j = 1, ..., m; k = 1, ..., krank, and
c
c                        krank
c       a(j,list(k))  =  Sigma  a(j,list(l)) * proj(l,k-krank)       (*)
c                         l=1
c
c                     +  epsilon(j,k-krank)
c
c       for all j = 1, ..., m; k = krank+1, ..., n,
c
c       for some matrix epsilon, dimensioned epsilon(m,n-krank),
c       whose norm is (hopefully) minimized by the pivoting procedure.
c       The present routine stores the krank x (n-krank) matrix proj
c       in the memory initially occupied by a.
c
c       input:
c       m -- first dimension of a
c       n -- second dimension of a, as well as the dimension required
c            of list
c       a -- matrix to be ID'd
c       krank -- desired rank of the output matrix
c                (please note that if krank > m or krank > n,
c                then the rank of the output matrix will be
c                less than krank)
c
c       output:
c       a -- the first krank*(n-krank) elements of a constitute
c            the krank x (n-krank) interpolation matrix proj
c       list -- list of the indices of the krank columns of a
c               through which the other columns of a are expressed;
c               also, list describes the permutation of proj
c               required to reconstruct a as indicated in (*) above
c       rnorms -- absolute values of the entries on the diagonal
c                 of the triangular matrix used to compute the ID
c                 (these may be used to check the stability of the ID)
c
c       _N.B._: This routine changes a.
c
c       reference:
c       Cheng, Gimbutas, Martinsson, Rokhlin, "On the compression of
c            low-rank matrices," SIAM Journal on Scientific Computing,
c            26 (4): 1389-1404, 2005.
c
        implicit none
        integer m,n,krank,j,k,list(n),iswap
        real*8 a(m,n),rnorms(n),ss
c
c
c       QR decompose a.
c
        call iddr_qrpiv(m,n,a,krank,list,rnorms)
c
c
c       Build the list of columns chosen in a
c       by multiplying together the permutations in list,
c       with the permutation swapping 1 and list(1) taken rightmost
c       in the product, that swapping 2 and list(2) taken next
c       rightmost, ..., that swapping krank and list(krank) taken
c       leftmost.
c
        do k = 1,n
          rnorms(k) = k
        enddo ! k
c
        if(krank .gt. 0) then
          do k = 1,krank
c
c           Swap rnorms(k) and rnorms(list(k)).
c
            iswap = rnorms(k)
            rnorms(k) = rnorms(list(k))
            rnorms(list(k)) = iswap
c
          enddo ! k
        endif
c
        do k = 1,n
          list(k) = rnorms(k)
        enddo ! k
c
c
c       Fill rnorms for the output.
c
        ss = 0
c
        do k = 1,krank
          rnorms(k) = a(k,k)
          ss = ss+rnorms(k)**2
        enddo ! k
c
c
c       Backsolve for proj, storing it at the beginning of a.
c
        if(krank .gt. 0 .and. ss .gt. 0) then
          call idd_lssolve(m,n,a,krank)
        endif
c
        if(ss .eq. 0) then
c
          do k = 1,n
            do j = 1,m
c
              a(j,k) = 0
c
            enddo ! j
          enddo ! k
c
        endif
c
c
        return
        end
c
c
c
c
        subroutine idd_reconid(m,krank,col,n,list,proj,approx)
c
c       reconstructs the matrix that the routine iddp_id
c       or iddr_id has decomposed, using the columns col
c       of the reconstructed matrix whose indices are listed in list,
c       in addition to the interpolation matrix proj.
c
c       input:
c       m -- first dimension of cols and approx
c       krank -- first dimension of cols and proj; also,
c                n-krank is the second dimension of proj
c       col -- columns of the matrix to be reconstructed
c       n -- second dimension of approx; also,
c            n-krank is the second dimension of proj
c       list(k) -- index of col(1:m,k) in the reconstructed matrix
c                  when k <= krank; in general, list describes
c                  the permutation required for reconstruction
c                  via cols and proj
c       proj -- interpolation matrix
c
c       output:
c       approx -- reconstructed matrix
c
        implicit none
        integer m,n,krank,j,k,l,list(n)
        real*8 col(m,krank),proj(krank,n-krank),approx(m,n)
c
c
        do j = 1,m
          do k = 1,n
c
            approx(j,list(k)) = 0
c
c           Add in the contributions due to the identity matrix.
c
            if(k .le. krank) then
              approx(j,list(k)) = approx(j,list(k)) + col(j,k)
            endif
c
c           Add in the contributions due to proj.
c
            if(k .gt. krank) then
              if(krank .gt. 0) then
c
                do l = 1,krank
                  approx(j,list(k)) = approx(j,list(k))
     1                              + col(j,l)*proj(l,k-krank)
                enddo ! l
c
              endif
            endif
c
          enddo ! k
        enddo ! j
c
c
        return
        end
c
c
c
c
        subroutine idd_lssolve(m,n,a,krank)
c
c       backsolves for proj satisfying R_11 proj ~ R_12,
c       where R_11 = a(1:krank,1:krank)
c       and R_12 = a(1:krank,krank+1:n).
c       This routine overwrites the beginning of a with proj.
c
c       input:
c       m -- first dimension of a
c       n -- second dimension of a; also,
c            n-krank is the second dimension of proj
c       a -- trapezoidal input matrix
c       krank -- first dimension of proj; also,
c                n-krank is the second dimension of proj
c
c       output:
c       a -- the first krank*(n-krank) elements of a constitute
c            the krank x (n-krank) matrix proj
c
        implicit none
        integer m,n,krank,j,k,l
        real*8 a(m,n),sum
c
c
c       Overwrite a(1:krank,krank+1:n) with proj.
c
        do k = 1,n-krank
          do j = krank,1,-1
c
            sum = 0
c
            do l = j+1,krank
              sum = sum+a(j,l)*a(l,krank+k)
            enddo ! l
c
            a(j,krank+k) = a(j,krank+k)-sum
c
c           Make sure that the entry in proj won't be too big;
c           set the entry to 0 when roundoff would make it too big
c           (in which case a(j,j) is so small that the contribution
c           from this entry in proj to the overall matrix approximation
c           is supposed to be negligible).
c
            if(abs(a(j,krank+k)) .lt. 2**20*abs(a(j,j))) then
              a(j,krank+k) = a(j,krank+k)/a(j,j)
            else
              a(j,krank+k) = 0
            endif
c
          enddo ! j
        enddo ! k
c
c
c       Move proj from a(1:krank,krank+1:n) to the beginning of a.
c
        call idd_moverup(m,n,krank,a)
c
c
        return
        end
c
c
c
c
        subroutine idd_moverup(m,n,krank,a)
c
c       moves the krank x (n-krank) matrix in a(1:krank,krank+1:n),
c       where a is initially dimensioned m x n, to the beginning of a.
c       (This is not the most natural way to code the move,
c       but one of my usually well-behaved compilers chokes
c       on more natural ways.)
c
c       input:
c       m -- initial first dimension of a
c       n -- initial second dimension of a
c       krank -- number of rows to move
c       a -- m x n matrix whose krank x (n-krank) block
c            a(1:krank,krank+1:n) is to be moved
c
c       output:
c       a -- array starting with the moved krank x (n-krank) block
c
        implicit none
        integer m,n,krank,j,k
        real*8 a(m*n)
c
c
        do k = 1,n-krank
          do j = 1,krank
            a(j+krank*(k-1)) = a(j+m*(krank+k-1))
          enddo ! j
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine idd_getcols(m,n,matvec,p1,p2,p3,p4,krank,list,
     1                         col,x)
c
c       collects together the columns of the matrix a indexed by list
c       into the matrix col, where routine matvec applies a
c       to an arbitrary vector.
c
c       input:
c       m -- first dimension of a
c       n -- second dimension of a
c       matvec -- routine which applies a to an arbitrary vector;
c                 this routine must have a calling sequence of the form
c
c                 matvec(m,x,n,y,p1,p2,p3,p4)
c
c                 where m is the length of x,
c                 x is the vector to which the matrix is to be applied,
c                 n is the length of y,
c                 y is the product of the matrix and x,
c                 and p1, p2, p3, and p4 are user-specified parameters
c       p1 -- parameter to be passed to routine matvec
c       p2 -- parameter to be passed to routine matvec
c       p3 -- parameter to be passed to routine matvec
c       p4 -- parameter to be passed to routine matvec
c       krank -- number of columns to be extracted
c       list -- indices of the columns to be extracted
c
c       output:
c       col -- columns of a indexed by list
c
c       work:
c       x -- must be at least n real*8 elements long
c
        implicit none
        integer m,n,krank,list(krank),j,k
        real*8 col(m,krank),x(n),p1,p2,p3,p4
        external matvec
c
c
        do j = 1,krank
c
          do k = 1,n
            x(k) = 0
          enddo ! k
c
          x(list(j)) = 1
c
          call matvec(n,x,m,col(1,j),p1,p2,p3,p4)
c
        enddo ! j
c
c
        return
        end
c
c
c
c
        subroutine idd_reconint(n,list,krank,proj,p)
c
c       constructs p in the ID a = b p,
c       where the columns of b are a subset of the columns of a,
c       and p is the projection coefficient matrix,
c       given list, krank, and proj output
c       by routines iddp_id or iddr_id.
c
c       input:
c       n -- part of the second dimension of proj and p
c       list -- list of columns retained from the original matrix
c               in the ID
c       krank -- rank of the ID
c       proj -- matrix of projection coefficients in the ID
c
c       output:
c       p -- projection matrix in the ID
c
        implicit none
        integer n,krank,list(n),j,k
        real*8 proj(krank,n-krank),p(krank,n)
c
c
        do k = 1,krank
          do j = 1,n
c
            if(j .le. krank) then
              if(j .eq. k) p(k,list(j)) = 1
              if(j .ne. k) p(k,list(j)) = 0
            endif
c
            if(j .gt. krank) then
              p(k,list(j)) = proj(k,j-krank)
            endif
c
          enddo ! j
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine idd_copycols(m,n,a,krank,list,col)
c
c       collects together the columns of the matrix a indexed by list
c       into the matrix col.
c
c       input:
c       m -- first dimension of a
c       n -- second dimension of a
c       a -- matrix whose columns are to be extracted
c       krank -- number of columns to be extracted
c       list -- indices of the columns to be extracted
c
c       output:
c       col -- columns of a indexed by list
c
        implicit none
        integer m,n,krank,list(krank),j,k
        real*8 a(m,n),col(m,krank)
c
c
        do k = 1,krank
          do j = 1,m
c
            col(j,k) = a(j,list(k))
c
          enddo ! j
        enddo ! k
c
c
        return
        end
