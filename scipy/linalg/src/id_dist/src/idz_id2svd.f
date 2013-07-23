c       this file contains the following user-callable routines:
c
c
c       routine idz_id2svd converts an approximation to a matrix
c       in the form of an ID to an approximation in the form of an SVD.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine idz_id2svd(m,krank,b,n,list,proj,u,v,s,ier,w)
c
c       converts an approximation to a matrix in the form of an ID
c       to an approximation in the form of an SVD.
c
c       input:
c       m -- first dimension of b
c       krank -- rank of the ID
c       b -- columns of the original matrix in the ID
c       list -- list of columns chosen from the original matrix
c               in the ID
c       n -- length of list and part of the second dimension of proj
c       proj -- projection coefficients in the ID
c
c       output:
c       u -- left singular vectors
c       v -- right singular vectors
c       s -- singular values
c       ier -- 0 when the routine terminates successfully;
c              nonzero otherwise
c
c       work:
c       w -- must be at least (krank+1)*(m+3*n+10)+9*krank**2
c            complex*16 elements long
c
c       _N.B._: This routine destroys b.
c
        implicit none
        integer m,krank,n,list(n),iwork,lwork,ip,lp,it,lt,ir,lr,
     1          ir2,lr2,ir3,lr3,iind,lind,iindt,lindt,lw,ier
        real*8 s(krank)
        complex*16 b(m,krank),proj(krank,n-krank),u(m,krank),
     1             v(n,krank),w((krank+1)*(m+3*n+10)+9*krank**2)
c
c
c       Allocate memory for idz_id2svd0.
c
        lw = 0
c
        iwork = lw+1
        lwork = 8*krank**2+10*krank
        lw = lw+lwork
c
        ip = lw+1
        lp = krank*n
        lw = lw+lp
c
        it = lw+1
        lt = n*krank
        lw = lw+lt
c
        ir = lw+1
        lr = krank*n
        lw = lw+lr
c
        ir2 = lw+1
        lr2 = krank*m
        lw = lw+lr2
c
        ir3 = lw+1
        lr3 = krank*krank
        lw = lw+lr3
c
        iind = lw+1
        lind = n/4+1
        lw = lw+1
c
        iindt = lw+1
        lindt = m/4+1
        lw = lw+1
c
c
        call idz_id2svd0(m,krank,b,n,list,proj,u,v,s,ier,
     1                   w(iwork),w(ip),w(it),w(ir),w(ir2),w(ir3),
     2                   w(iind),w(iindt))
c
c
        return
        end
c
c
c
c
        subroutine idz_id2svd0(m,krank,b,n,list,proj,u,v,s,ier,
     1                         work,p,t,r,r2,r3,ind,indt)
c
c       routine idz_id2svd serves as a memory wrapper
c       for the present routine (please see routine idz_id2svd
c       for further documentation).
c
        implicit none
c
        character*1 jobz
        integer m,n,krank,list(n),ind(n),indt(m),ifadjoint,
     1          lwork,ldu,ldvt,ldr,info,j,k,ier
        real*8 s(krank)
        complex*16 b(m,krank),proj(krank,n-krank),p(krank,n),
     1             r(krank,n),r2(krank,m),t(n,krank),r3(krank,krank),
     2             u(m,krank),v(n,krank),work(8*krank**2+10*krank)
c
c
c
        ier = 0
c
c
c
c       Construct the projection matrix p from the ID.
c
        call idz_reconint(n,list,krank,proj,p)
c
c
c
c       Compute a pivoted QR decomposition of b.
c
        call idzr_qrpiv(m,krank,b,krank,ind,r)
c
c
c       Extract r from the QR decomposition.
c
        call idz_rinqr(m,krank,b,krank,r)
c
c
c       Rearrange r according to ind.
c
        call idz_rearr(krank,ind,krank,krank,r)
c
c
c
c       Take the adjoint of p to obtain t.
c
        call idz_matadj(krank,n,p,t)
c
c
c       Compute a pivoted QR decomposition of t.
c
        call idzr_qrpiv(n,krank,t,krank,indt,r2)
c
c
c       Extract r2 from the QR decomposition.
c
        call idz_rinqr(n,krank,t,krank,r2)
c
c
c       Rearrange r2 according to indt.
c
        call idz_rearr(krank,indt,krank,krank,r2)
c
c
c
c       Multiply r and r2^* to obtain r3.
c
        call idz_matmulta(krank,krank,r,krank,r2,r3)
c
c
c
c       Use LAPACK to SVD r3.
c
        jobz = 'S'
        ldr = krank
        lwork = 8*krank**2+10*krank
     1        - (krank**2+2*krank+3*krank**2+4*krank)
        ldu = krank
        ldvt = krank
c
        call zgesdd(jobz,krank,krank,r3,ldr,s,work,ldu,r,ldvt,
     1              work(krank**2+2*krank+3*krank**2+4*krank+1),lwork,
     2              work(krank**2+2*krank+1),work(krank**2+1),info)
c
        if(info .ne. 0) then
          ier = info
          return
        endif
c
c
c
c       Multiply the u from r3 from the left by the q from b
c       to obtain the u for a.
c
        do k = 1,krank
c
          do j = 1,krank
            u(j,k) = work(j+krank*(k-1))
          enddo ! j
c
          do j = krank+1,m
            u(j,k) = 0
          enddo ! j
c
        enddo ! k
c
        ifadjoint = 0
        call idz_qmatmat(ifadjoint,m,krank,b,krank,krank,u,r2)
c
c
c
c       Take the adjoint of r to obtain r2.
c
        call idz_matadj(krank,krank,r,r2)
c
c
c       Multiply the v from r3 from the left by the q from p^*
c       to obtain the v for a.
c
        do k = 1,krank
c
          do j = 1,krank
            v(j,k) = r2(j,k)
          enddo ! j
c
          do j = krank+1,n
            v(j,k) = 0
          enddo ! j
c
        enddo ! k
c
        ifadjoint = 0
        call idz_qmatmat(ifadjoint,n,krank,t,krank,krank,v,r2)
c
c
        return
        end
c
c
c
c
        subroutine idz_matadj(m,n,a,aa)
c
c       Takes the adjoint of a to obtain aa.
c
c       input:
c       m -- first dimension of a, and second dimension of aa
c       n -- second dimension of a, and first dimension of aa
c       a -- matrix whose adjoint is to be taken
c
c       output:
c       aa -- adjoint of a
c
        implicit none
        integer m,n,j,k
        complex*16 a(m,n),aa(n,m)
c
c
        do k = 1,n
          do j = 1,m
            aa(k,j) = conjg(a(j,k))
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
        subroutine idz_matmulta(l,m,a,n,b,c)
c
c       multiplies a and b^* to obtain c.
c
c       input:
c       l -- first dimension of a and c
c       m -- second dimension of a and b
c       a -- leftmost matrix in the product c = a b^*
c       n -- first dimension of b and second dimension of c
c       b -- rightmost matrix in the product c = a b^*
c
c       output:
c       c -- product of a and b^*
c
        implicit none
        integer l,m,n,i,j,k
        complex*16 a(l,m),b(n,m),c(l,n),sum
c
c
        do i = 1,l
          do k = 1,n
c
            sum = 0
c
            do j = 1,m
              sum = sum+a(i,j)*conjg(b(k,j))
            enddo ! j
c
            c(i,k) = sum
c
          enddo ! k
        enddo ! i
c
c
        return
        end
c
c
c
c
        subroutine idz_rearr(krank,ind,m,n,a)
c
c       rearranges a according to ind obtained
c       from routines idzr_qrpiv or idzp_qrpiv,
c       assuming that a = q r, where q and r are from idzr_qrpiv
c       or idzp_qrpiv.
c
c       input:
c       krank -- rank obtained from routine idzp_qrpiv,
c                or provided to routine idzr_qrpiv
c       ind -- indexing array obtained from routine idzr_qrpiv
c              or idzp_qrpiv
c       m -- first dimension of a
c       n -- second dimension of a
c       a -- matrix to be rearranged
c
c       output:
c       a -- rearranged matrix
c
        implicit none
        integer k,krank,m,n,j,ind(krank)
        complex*16 cswap,a(m,n)
c
c
        do k = krank,1,-1
          do j = 1,m
c
            cswap = a(j,k)
            a(j,k) = a(j,ind(k))
            a(j,ind(k)) = cswap
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
        subroutine idz_rinqr(m,n,a,krank,r)
c
c       extracts R in the QR decomposition specified by the output a
c       of the routine idzr_qrpiv or idzp_qrpiv.
c
c       input:
c       m -- first dimension of a
c       n -- second dimension of a and r
c       a -- output of routine idzr_qrpiv or idzp_qrpiv
c       krank -- rank output by routine idzp_qrpiv (or specified
c                to routine idzr_qrpiv)
c
c       output:
c       r -- triangular factor in the QR decomposition specified
c            by the output a of the routine idzr_qrpiv or idzp_qrpiv
c
        implicit none
        integer m,n,j,k,krank
        complex*16 a(m,n),r(krank,n)
c
c
c       Copy a into r and zero out the appropriate
c       Householder vectors that are stored in one triangle of a.
c
        do k = 1,n
          do j = 1,krank
            r(j,k) = a(j,k)
          enddo ! j
        enddo ! k
c
        do k = 1,n
          if(k .lt. krank) then
            do j = k+1,krank
              r(j,k) = 0
            enddo ! j
          endif
        enddo ! k
c
c
        return
        end
