c       this file contains the following user-callable routines:
c
c
c       routine idzr_aid computes the ID, to a specified rank,
c       of an arbitrary matrix. This routine is randomized.
c
c       routine idzr_aidi initializes routine idzr_aid.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine idzr_aid(m,n,a,krank,w,list,proj)
c
c       computes the ID of the matrix a, i.e., lists in list
c       the indices of krank columns of a such that
c
c       a(j,list(k))  =  a(j,list(k))
c
c       for all j = 1, ..., m; k = 1, ..., krank, and
c
c                       min(m,n,krank)
c       a(j,list(k))  =     Sigma      a(j,list(l)) * proj(l,k-krank)(*)
c                            l=1
c
c                     +  epsilon(j,k-krank)
c
c       for all j = 1, ..., m; k = krank+1, ..., n,
c
c       for some matrix epsilon, dimensioned epsilon(m,n-krank),
c       whose norm is (hopefully) minimized by the pivoting procedure.
c
c       input:
c       m -- number of rows in a
c       n -- number of columns in a
c       a -- matrix to be ID'd; the present routine does not alter a
c       krank -- rank of the ID to be constructed
c       w -- initialization array that routine idzr_aidi
c            has constructed
c
c       output:
c       list -- indices of the columns in the ID
c       proj -- matrix of coefficients needed to interpolate
c               from the selected columns to the other columns
c               in the original matrix being ID'd
c
c       _N.B._: The algorithm used by this routine is randomized.
c
c       reference:
c       Halko, Martinsson, Tropp, "Finding structure with randomness:
c            probabilistic algorithms for constructing approximate
c            matrix decompositions," SIAM Review, 53 (2): 217-288,
c            2011.
c
        implicit none
        integer m,n,krank,list(n),lw,ir,lr,lw2,iw
        complex*16 a(m,n),proj(krank*(n-krank)),
     1             w((2*krank+17)*n+21*m+80)
c
c
c       Allocate memory in w.
c
        lw = 0
c
        iw = lw+1
        lw2 = 21*m+80+n
        lw = lw+lw2
c
        ir = lw+1
        lr = (krank+8)*2*n
        lw = lw+lr
c
c
        call idzr_aid0(m,n,a,krank,w(iw),list,proj,w(ir))
c
c
        return
        end
c
c
c
c
        subroutine idzr_aid0(m,n,a,krank,w,list,proj,r)
c
c       routine idzr_aid serves as a memory wrapper
c       for the present routine
c       (see idzr_aid for further documentation).
c
        implicit none
        integer k,l,m,n2,n,krank,list(n),mn,lproj
        complex*16 a(m,n),r(krank+8,2*n),proj(krank,n-krank),
     1             w(21*m+80+n)
c
c       Please note that the second dimension of r is 2*n
c       (instead of n) so that if krank+8 >= m/2, then
c       we can copy the whole of a into r.
c
c
c       Retrieve the number of random test vectors
c       and the greatest integer less than m that is
c       a positive integer power of two.
c
        l = w(1)
        n2 = w(2)
c
c
        if(l .lt. n2 .and. l .le. m) then
c
c         Apply the random matrix.
c
          do k = 1,n
            call idz_sfrm(l,m,n2,w(11),a(1,k),r(1,k))
          enddo ! k
c
c         ID r.
c
          call idzr_id(l,n,r,krank,list,w(20*m+81))
c
c         Retrieve proj from r.
c
          lproj = krank*(n-krank)
          call idzr_copyzarr(lproj,r,proj)
c
        endif
c
c
        if(l .ge. n2 .or. l .gt. m) then
c
c         ID a directly.
c
          mn = m*n
          call idzr_copyzarr(mn,a,r)
          call idzr_id(m,n,r,krank,list,w(20*m+81))
c
c         Retrieve proj from r.
c
          lproj = krank*(n-krank)
          call idzr_copyzarr(lproj,r,proj)
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
        subroutine idzr_copyzarr(n,a,b)
c
c       copies a into b.
c
c       input:
c       n -- length of a and b
c       a -- array to copy into b
c
c       output:
c       b -- copy of a
c
        implicit none
        integer n,k
        complex*16 a(n),b(n)
c
c
        do k = 1,n
          b(k) = a(k)
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine idzr_aidi(m,n,krank,w)
c
c       initializes the array w for using routine idzr_aid.
c
c       input:
c       m -- number of rows in the matrix to be ID'd
c       n -- number of columns in the matrix to be ID'd
c       krank -- rank of the ID to be constructed
c
c       output:
c       w -- initialization array for using routine idzr_aid
c
        implicit none
        integer m,n,krank,l,n2
        complex*16 w((2*krank+17)*n+21*m+80)
c
c
c       Set the number of random test vectors to 8 more than the rank.
c
        l = krank+8
        w(1) = l
c
c
c       Initialize the rest of the array w.
c
        n2 = 0
        if(l .le. m) call idz_sfrmi(l,m,n2,w(11))
        w(2) = n2
c
c
        return
        end
