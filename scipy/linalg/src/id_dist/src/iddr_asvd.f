c       this file contains the following user-callable routines:
c
c
c       routine iddr_aid computes the SVD, to a specified rank,
c       of an arbitrary matrix. This routine is randomized.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine iddr_asvd(m,n,a,krank,w,u,v,s,ier)
c
c       constructs a rank-krank SVD  u diag(s) v^T  approximating a,
c       where u is an m x krank matrix whose columns are orthonormal,
c       v is an n x krank matrix whose columns are orthonormal,
c       and diag(s) is a diagonal krank x krank matrix whose entries
c       are all nonnegative. This routine uses a randomized algorithm.
c
c       input:
c       m -- number of rows in a
c       n -- number of columns in a 
c       a -- matrix to be decomposed; the present routine does not
c            alter a
c       krank -- rank of the SVD being constructed
c       w -- initialization array that routine iddr_aidi
c            has constructed (for use in the present routine, w must
c            be at least (2*krank+28)*m+(6*krank+21)*n+25*krank**2+100
c            real*8 elements long)
c
c       output:
c       u -- matrix of orthonormal left singular vectors of a
c       v -- matrix of orthonormal right singular vectors of a
c       s -- array of singular values of a
c       ier -- 0 when the routine terminates successfully;
c              nonzero otherwise
c
c       _N.B._: The algorithm used by this routine is randomized.
c
        implicit none
        integer m,n,krank,lw,ilist,llist,iproj,lproj,icol,lcol,
     1          iwork,lwork,iwinit,lwinit,ier
        real*8 a(m,n),u(m,krank),v(n,krank),s(krank),
     1         w((2*krank+28)*m+(6*krank+21)*n+25*krank**2+100)
c
c
c       Allocate memory in w.
c
        lw = 0
c
        iwinit = lw+1
        lwinit = (2*krank+17)*n+27*m+100
        lw = lw+lwinit
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
        call iddr_asvd0(m,n,a,krank,w(iwinit),u,v,s,ier,
     1                  w(ilist),w(iproj),w(icol),w(iwork))
c
c
        return
        end
c
c
c
c
        subroutine iddr_asvd0(m,n,a,krank,winit,u,v,s,ier,
     1                        list,proj,col,work)
c
c       routine iddr_asvd serves as a memory wrapper
c       for the present routine (please see routine iddr_asvd
c       for further documentation).
c
        implicit none
        integer m,n,krank,list(n),ier
        real*8 a(m,n),u(m,krank),v(n,krank),s(krank),
     1         proj(krank,n-krank),col(m*krank),
     2         winit((2*krank+17)*n+27*m+100),
     3         work((krank+1)*(m+3*n)+26*krank**2)
c
c
c       ID a.
c
        call iddr_aid(m,n,a,krank,winit,list,proj)
c
c
c       Collect together the columns of a indexed by list into col.
c
        call idd_copycols(m,n,a,krank,list,col)
c
c
c       Convert the ID to an SVD.
c
        call idd_id2svd(m,krank,col,n,list,proj,u,v,s,ier,work)
c
c
        return
        end
