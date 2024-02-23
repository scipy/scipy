c       this file contains the following user-callable routines:
c
c
c       routine idzp_rsvd computes the SVD, to a specified precision,
c       of a matrix specified by routines for applying the matrix
c       and its adjoint to arbitrary vectors.
c       This routine is randomized.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine idzp_rsvd(lw,eps,m,n,matveca,p1t,p2t,p3t,p4t,
     1                       matvec,p1,p2,p3,p4,krank,iu,iv,is,w,ier)
c
c       constructs a rank-krank SVD  U Sigma V^*  approximating a
c       to precision eps, where matveca is a routine which applies a^*
c       to an arbitrary vector, and matvec is a routine
c       which applies a to an arbitrary vector; U is an m x krank
c       matrix whose columns are orthonormal, V is an n x krank
c       matrix whose columns are orthonormal, and Sigma is a diagonal
c       krank x krank matrix whose entries are all nonnegative.
c       The entries of U are stored in w, starting at w(iu);
c       the entries of V are stored in w, starting at w(iv).
c       The diagonal entries of Sigma are stored in w,
c       starting at w(is). This routine uses a randomized algorithm.
c
c       input:
c       lw -- maximum usable length (in complex*16 elements)
c             of the array w
c       eps -- precision of the desired approximation
c       m -- number of rows in a
c       n -- number of columns in a
c       matveca -- routine which applies the adjoint
c                  of the matrix to be SVD'd
c                  to an arbitrary vector; this routine must have
c                  a calling sequence of the form
c
c                  matveca(m,x,n,y,p1t,p2t,p3t,p4t),
c
c                  where m is the length of x,
c                  x is the vector to which the adjoint
c                  of the matrix is to be applied,
c                  n is the length of y,
c                  y is the product of the adjoint of the matrix and x,
c                  and p1t, p2t, p3t, and p4t are user-specified
c                  parameters
c       p1t -- parameter to be passed to routine matveca
c       p2t -- parameter to be passed to routine matveca
c       p3t -- parameter to be passed to routine matveca
c       p4t -- parameter to be passed to routine matveca
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
c
c       output:
c       krank -- rank of the SVD constructed
c       iu -- index in w of the first entry of the matrix
c             of orthonormal left singular vectors of a
c       iv -- index in w of the first entry of the matrix
c             of orthonormal right singular vectors of a
c       is -- index in w of the first entry of the array
c             of singular values of a; the singular values are stored
c             as complex*16 numbers whose imaginary parts are zeros
c       w -- array containing the singular values and singular vectors
c            of a; w doubles as a work array, and so must be at least
c            (krank+1)*(3*m+5*n+11)+8*krank^2 complex*16 elements long,
c            where krank is the rank returned by the present routine
c       ier -- 0 when the routine terminates successfully;
c              -1000 when lw is too small;
c              other nonzero values when idz_id2svd bombs
c
c       _N.B._: w must be at least (krank+1)*(3*m+5*n+11)+8*krank**2
c               complex*16 elements long, where krank is the rank
c               returned by the present routine. Also, the algorithm
c               used by the present routine is randomized.
c
        implicit none
        integer m,n,krank,lw,lw2,ilist,llist,iproj,icol,lcol,lp,
     1          iwork,lwork,ier,lproj,iu,iv,is,lu,lv,ls,iui,ivi,isi,k
        real*8 eps
        complex*16 p1t,p2t,p3t,p4t,p1,p2,p3,p4,w(*)
        external matveca,matvec
c
c
c       Allocate some memory.
c
        lw2 = 0
c
        ilist = lw2+1
        llist = n
        lw2 = lw2+llist
c
        iproj = lw2+1
c
c
c       ID a.
c
        lp = lw-lw2
        call idzp_rid(lp,eps,m,n,matveca,p1t,p2t,p3t,p4t,krank,
     1                w(ilist),w(iproj),ier)
        if(ier .ne. 0) return
c
c
        if(krank .gt. 0) then
c
c
c         Allocate more memory.
c
          lproj = krank*(n-krank)
          lw2 = lw2+lproj
c
          icol = lw2+1
          lcol = m*krank
          lw2 = lw2+lcol
c
          iui = lw2+1
          lu = m*krank
          lw2 = lw2+lu
c
          ivi = lw2+1
          lv = n*krank
          lw2 = lw2+lv
c
          isi = lw2+1
          ls = krank
          lw2 = lw2+ls
c
          iwork = lw2+1
          lwork = (krank+1)*(m+3*n+10)+9*krank**2
          lw2 = lw2+lwork
c
c
          if(lw .lt. lw2) then
            ier = -1000
            return
          endif
c
c
          call idzp_rsvd0(m,n,matveca,p1t,p2t,p3t,p4t,
     1                    matvec,p1,p2,p3,p4,krank,w(iui),w(ivi),
     2                    w(isi),ier,w(ilist),w(iproj),w(icol),
     3                    w(iwork))
          if(ier .ne. 0) return
c
c
          iu = 1
          iv = iu+lu
          is = iv+lv
c
c
c         Copy the singular values and singular vectors
c         into their proper locations.
c
          do k = 1,lu
            w(iu+k-1) = w(iui+k-1)
          enddo ! k
c
          do k = 1,lv
            w(iv+k-1) = w(ivi+k-1)
          enddo ! k
c
          call idz_reco(ls,w(isi),w(is))
c
c
        endif ! krank .gt. 0
c
c
        return
        end
c
c
c
c
        subroutine idzp_rsvd0(m,n,matveca,p1t,p2t,p3t,p4t,
     1                        matvec,p1,p2,p3,p4,krank,u,v,s,ier,
     2                        list,proj,col,work)
c
c       routine idzp_rsvd serves as a memory wrapper
c       for the present routine (please see routine idzp_rsvd
c       for further documentation).
c
        implicit none
        integer m,n,krank,list(n),ier
        real*8 s(krank)
        complex*16 p1t,p2t,p3t,p4t,p1,p2,p3,p4,u(m,krank),v(n,krank),
     1             proj(krank,n-krank),col(m*krank),
     2             work((krank+1)*(m+3*n+10)+9*krank**2)
        external matveca,matvec
c
c
c       Collect together the columns of a indexed by list into col.
c
        call idz_getcols(m,n,matvec,p1,p2,p3,p4,krank,list,col,work)
c
c
c       Convert the ID to an SVD.
c
        call idz_id2svd(m,krank,col,n,list,proj,u,v,s,ier,work)
c
c
        return
        end
c
c
c
c
        subroutine idz_reco(n,a,b)
c
c       copies the real*8 array a into the complex*16 array b.
c
c       input:
c       n -- length of a and b
c       a -- real*8 array to be copied into b
c
c       output:
c       b -- complex*16 copy of a
c
        integer n,k
        real*8 a(n)
        complex*16 b(n)
c
c
        do k = 1,n
          b(k) = a(k)
        enddo ! k
c
c
        return
        end
