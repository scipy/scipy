c       this file contains the following user-callable routines:
c
c
c       routine iddr_rsvd computes the SVD, to a specified rank,
c       of a matrix specified by routines for applying the matrix
c       and its transpose to arbitrary vectors.
c       This routine is randomized.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine iddr_rsvd(m,n,matvect,p1t,p2t,p3t,p4t,
     1                       matvec,p1,p2,p3,p4,krank,u,v,s,ier,w)
c
c       constructs a rank-krank SVD  u diag(s) v^T  approximating a,
c       where matvect is a routine which applies a^T
c       to an arbitrary vector, and matvec is a routine
c       which applies a to an arbitrary vector;
c       u is an m x krank matrix whose columns are orthonormal,
c       v is an n x krank matrix whose columns are orthonormal,
c       and diag(s) is a diagonal krank x krank matrix whose entries
c       are all nonnegative. This routine uses a randomized algorithm.
c
c       input:
c       m -- number of rows in a
c       n -- number of columns in a
c       matvect -- routine which applies the transpose
c                  of the matrix to be SVD'd
c                  to an arbitrary vector; this routine must have
c                  a calling sequence of the form
c
c                  matvect(m,x,n,y,p1t,p2t,p3t,p4t),
c
c                  where m is the length of x,
c                  x is the vector to which the transpose
c                  of the matrix is to be applied,
c                  n is the length of y,
c                  y is the product of the transposed matrix and x,
c                  and p1t, p2t, p3t, and p4t are user-specified
c                  parameters
c       p1t -- parameter to be passed to routine matvect
c       p2t -- parameter to be passed to routine matvect
c       p3t -- parameter to be passed to routine matvect
c       p4t -- parameter to be passed to routine matvect
c       matvec -- routine which applies the matrix to be SVD'd
c                 to an arbitrary vector; this routine must have
c                 a calling sequence of the form
c
c                 matvec(n,x,m,y,p1,p2,p3,p4),
c
c                 where n is the length of x,
c                 x is the vector to which the matrix is to be applied,
c                 m is the length of y,
c                 y is the product of the matrix and x,
c                 and p1, p2, p3, and p4 are user-specified parameters
c       p1 -- parameter to be passed to routine matvec
c       p2 -- parameter to be passed to routine matvec
c       p3 -- parameter to be passed to routine matvec
c       p4 -- parameter to be passed to routine matvec
c       krank -- rank of the SVD being constructed
c
c       output:
c       u -- matrix of orthonormal left singular vectors of a
c       v -- matrix of orthonormal right singular vectors of a
c       s -- array of singular values of a
c       ier -- 0 when the routine terminates successfully;
c              nonzero otherwise
c
c       work:
c       w -- must be at least (krank+1)*(2*m+4*n)+25*krank**2
c            real*8 elements long
c
c       _N.B._: The algorithm used by this routine is randomized.
c
        implicit none
        integer m,n,krank,lw,ilist,llist,iproj,lproj,icol,lcol,
     1          iwork,lwork,ier
        real*8 p1t,p2t,p3t,p4t,p1,p2,p3,p4,u(m,krank),v(n,krank),
     1         s(krank),w((krank+1)*(2*m+4*n)+25*krank**2)
        external matvect,matvec
c
c
c       Allocate memory in w.
c
        lw = 0
c
        ilist = lw+1
        llist = n
        lw = lw+llist
c
        iproj = lw+1
        lproj = krank*(n-krank)
        lw = lw+lproj
c
        icol = lw+1
        lcol = m*krank
        lw = lw+lcol
c
        iwork = lw+1
        lwork = (krank+1)*(m+3*n)+26*krank**2
        lw = lw+lwork
c
c
        call iddr_rsvd0(m,n,matvect,p1t,p2t,p3t,p4t,
     1                  matvec,p1,p2,p3,p4,krank,u,v,s,ier,
     2                  w(ilist),w(iproj),w(icol),w(iwork))
c
c
        return
        end
c
c
c
c
        subroutine iddr_rsvd0(m,n,matvect,p1t,p2t,p3t,p4t,
     1                        matvec,p1,p2,p3,p4,krank,u,v,s,ier,
     2                        list,proj,col,work)
c
c       routine iddr_rsvd serves as a memory wrapper
c       for the present routine (please see routine iddr_rsvd
c       for further documentation).
c
        implicit none
        integer m,n,krank,list(n),ier,k
        real*8 p1t,p2t,p3t,p4t,p1,p2,p3,p4,u(m,krank),v(n,krank),
     1         s(krank),proj(krank*(n-krank)),col(m*krank),
     2         work((krank+1)*(m+3*n)+26*krank**2)
        external matvect,matvec
c
c
c       ID a.
c
        call iddr_rid(m,n,matvect,p1t,p2t,p3t,p4t,krank,list,work)
c
c
c       Retrieve proj from work.
c
        do k = 1,krank*(n-krank)
          proj(k) = work(k)
        enddo ! k
c
c
c       Collect together the columns of a indexed by list into col.
c
        call idd_getcols(m,n,matvec,p1,p2,p3,p4,krank,list,col,work)
c
c
c       Convert the ID to an SVD.
c
        call idd_id2svd(m,krank,col,n,list,proj,u,v,s,ier,work)
c
c
        return
        end
