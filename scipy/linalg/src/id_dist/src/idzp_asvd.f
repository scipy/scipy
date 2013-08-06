c       this file contains the following user-callable routines:
c
c
c       routine idzp_asvd computes the SVD, to a specified precision,
c       of an arbitrary matrix. This routine is randomized.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine idzp_asvd(lw,eps,m,n,a,winit,krank,iu,iv,is,w,ier)
c
c       constructs a rank-krank SVD  U Sigma V^*  approximating a
c       to precision eps, where U is an m x krank matrix whose
c       columns are orthonormal, V is an n x krank matrix whose
c       columns are orthonormal, and Sigma is a diagonal krank x krank
c       matrix whose entries are all nonnegative.
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
c       a -- matrix to be approximated; the present routine does not
c            alter a
c       winit -- initialization array that has been constructed
c                by routine idz_frmi
c
c       output:
c       krank -- rank of the SVD constructed
c       iu -- index in w of the first entry of the matrix
c             of orthonormal left singular vectors of a
c       iv -- index in w of the first entry of the matrix
c             of orthonormal right singular vectors of a
c       is -- index in w of the first entry of the array
c             of singular values of a
c       w -- array containing the singular values and singular vectors
c            of a; w doubles as a work array, and so must be at least
c            max( (krank+1)*(3*m+5*n+11)+8*krank**2, (2*n+1)*(n2+1) )
c            complex*16 elements long, where n2 is the greatest integer
c            less than or equal to m, such that n2 is
c            a positive integer power of two; krank is the rank output
c            by this routine
c       ier -- 0 when the routine terminates successfully;
c              -1000 when lw is too small;
c              other nonzero values when idz_id2svd bombs
c
c       _N.B._: w must be at least
c               max( (krank+1)*(3*m+5*n+11)+8*krank^2, (2*n+1)*(n2+1) )
c               complex*16 elements long, where n2 is
c               the greatest integer less than or equal to m,
c               such that n2 is a positive integer power of two;
c               krank is the rank output by this routine.
c               Also, the algorithm used by this routine is randomized.
c
        implicit none
        integer m,n,krank,lw,ilist,llist,iproj,lproj,icol,lcol,
     1          iwork,lwork,k,ier,lw2,iu,iv,is,iui,ivi,isi,lu,lv,ls
        real*8 eps
        complex*16 a(m,n),winit(17*m+70),w(*)
c
c
c       Allocate memory in w.
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
        call idzp_aid(eps,m,n,a,winit,krank,w(ilist),w(iproj))
c
c
        if(krank .gt. 0) then
c
c
c         Allocate more memory in w.
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
          call idzp_asvd0(m,n,a,krank,w(ilist),w(iproj),
     1                    w(iui),w(ivi),w(isi),ier,w(icol),w(iwork))
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
          call idz_realcomplex(ls,w(isi),w(is))
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
        subroutine idzp_asvd0(m,n,a,krank,list,proj,u,v,s,ier,
     1                        col,work)
c
c       routine idzp_asvd serves as a memory wrapper
c       for the present routine (please see routine idzp_asvd
c       for further documentation).
c
        implicit none
        integer m,n,krank,list(n),ier
        real*8 s(krank)
        complex*16 a(m,n),u(m,krank),v(n,krank),
     1             proj(krank,n-krank),col(m,krank),
     2             work((krank+1)*(m+3*n+10)+9*krank**2)
c
c
c       Collect together the columns of a indexed by list into col.
c
        call idz_copycols(m,n,a,krank,list,col)
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
        subroutine idz_realcomplex(n,a,b)
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
