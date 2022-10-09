c       this file contains the following user-callable routines:
c
c
c       routine idzr_rid computes the ID, to a specified rank,
c       of a matrix specified by a routine for applying its adjoint
c       to arbitrary vectors. This routine is randomized.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine idzr_rid(m,n,matveca,p1,p2,p3,p4,krank,list,proj)
c
c       computes the ID of a matrix "a" specified by
c       the routine matveca -- matveca must apply the adjoint
c       of the matrix being ID'd to an arbitrary vector --
c       i.e., the present routine lists in list the indices
c       of krank columns of a such that
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
c       m -- number of rows in the matrix to be ID'd
c       n -- number of columns in the matrix to be ID'd
c       matveca -- routine which applies the adjoint
c                  of the matrix to be ID'd to an arbitrary vector;
c                  this routine must have a calling sequence
c                  of the form
c
c                  matveca(m,x,n,y,p1,p2,p3,p4),
c
c                  where m is the length of x,
c                  x is the vector to which the adjoint
c                  of the matrix is to be applied,
c                  n is the length of y,
c                  y is the product of the adjoint of the matrix and x,
c                  and p1, p2, p3, and p4 are user-specified parameters
c       p1 -- parameter to be passed to routine matveca
c       p2 -- parameter to be passed to routine matveca
c       p3 -- parameter to be passed to routine matveca
c       p4 -- parameter to be passed to routine matveca
c       krank -- rank of the ID to be constructed
c
c       output:
c       list -- indices of the columns in the ID
c       proj -- matrix of coefficients needed to interpolate
c               from the selected columns to the other columns
c               in the original matrix being ID'd;
c               proj doubles as a work array in the present routine, so
c               proj must be at least m+(krank+3)*n complex*16 elements
c               long
c
c       _N.B._: The algorithm used by this routine is randomized.
c               proj must be at least m+(krank+3)*n complex*16 elements
c               long.
c
c       reference:
c       Halko, Martinsson, Tropp, "Finding structure with randomness:
c            probabilistic algorithms for constructing approximate
c            matrix decompositions," SIAM Review, 53 (2): 217-288,
c            2011.
c
        implicit none
        integer m,n,krank,list(n),lw,ix,lx,iy,ly,ir,lr
        complex*16 p1,p2,p3,p4,proj(m+(krank+3)*n)
        external matveca
c
c
c       Allocate memory in w.
c
        lw = 0
c
        ir = lw+1
        lr = (krank+2)*n
        lw = lw+lr
c
        ix = lw+1
        lx = m
        lw = lw+lx
c
        iy = lw+1
        ly = n
        lw = lw+ly
c
c
        call idzr_ridall0(m,n,matveca,p1,p2,p3,p4,krank,
     1                    list,proj(ir),proj(ix),proj(iy))
c
c
        return
        end
c
c
c
c
        subroutine idzr_ridall0(m,n,matveca,p1,p2,p3,p4,krank,
     1                          list,r,x,y)
c
c       routine idzr_ridall serves as a memory wrapper
c       for the present routine
c       (see idzr_ridall for further documentation).
c
        implicit none
        integer j,k,l,m,n,krank,list(n),m2
        complex*16 x(m),y(n),p1,p2,p3,p4,r(krank+2,n)
        external matveca
c
c
c       Set the number of random test vectors to 2 more than the rank.
c
        l = krank+2
c
c       Apply the adjoint of the original matrix to l random vectors.
c
        do j = 1,l
c
c         Generate a random vector.
c
          m2 = m*2
          call id_srand(m2,x)
c
c         Apply the adjoint of the matrix to x, obtaining y.
c
          call matveca(m,x,n,y,p1,p2,p3,p4)
c
c         Copy the conjugate of y into row j of r.
c
          do k = 1,n
            r(j,k) = conjg(y(k))
          enddo ! k
c
        enddo ! j
c
c
c       ID r.
c
        call idzr_id(l,n,r,krank,list,y)
c
c
        return
        end
