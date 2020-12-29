c
c     Rasmus Munk Larsen, Stanford University, 1999, 2004.
c

      subroutine dreorth(n,k,V,ldv,vnew,normvnew,index,alpha,work,
     c     iflag)
c
c     Orthogonalize the N-vector VNEW against a subset of the columns of
c     the N-by-K matrix V(1:N,1:K) using iterated classical or modified
c     Gram-Schmidt. LDV is the leading dimension of the array containing
c     V.
c     
c     Which columns to orthogonalize against is decided by the integer
c     array INDEX = [s_1,e_1, s_2,e_2,..., s_k,e_l, s_{l+1}], which
c     selects the columns V(:,[s_1:e_1 s_2:e_2 ... s_l:e_l]). s_{l+1}
c     must be larger than k and marks the end of INDEX.
c
c     The reorthogonalization is repeated until
c
c       ||VNEW'|| > ALPHA * ||VNEW|| , 
c
c     where VNEW' is the vector obtained by orthogonalizing VNEW.  If
c     VNEW' fails to satisfy this after 4 tries, VNEW is deemed to lie
c     numerically in the span of V(:,[s_1:e_1 s_2:e_2 ... s_l:e_l]), and
c     is set to the zero vector.
c
c     On return NORMVNEW contains ||VNEW||.
c
c     WORK is a workspace array of length at least 
c
c       max e_i-s_i+1, i=1...l.
c     
c     WORK is only used if IFLAG==1.
c
c     If IFLAG==0 then iterated modified Gram-Schmidt is used.
c     If IFLAG==1 then iterated classical Gram-Schmidt is used.
c

c     References: 
c       Aake Bjorck, "Numerical Methods for Least Squares Problems",
c       SIAM, Philadelphia, 1996, pp. 68-69.
c     
c       J.~W. Daniel, W.~B. Gragg, L. Kaufman and G.~W. Stewart, 
c       ``Reorthogonalization and Stable Algorithms Updating the
c       Gram-Schmidt QR Factorization'', Math. Comp.,  30 (1976), no.
c       136, pp. 772-795.
c
c       B. N. Parlett, ``The Symmetric Eigenvalue Problem'', 
c       Prentice-Hall, Englewood Cliffs, NJ, 1980. pp. 105-109

c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      integer n,k,ldv,iflag,index(*)
      double precision V(ldv,*),vnew(*),work(*),normvnew

c     %------------%
c     | Parameters |
c     %------------%
      integer NTRY
      double precision one, zero
      parameter(one = 1.0, zero = 0.0, NTRY=5)
      
c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer itry
      double precision alpha,normvnew_0
      real t2,t3
      
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
      external dmgs,dcgs
      
c     %--------------------%
c     | External Functions |
c     %--------------------%
      double precision pdnrm2
      external pdnrm2,pdzero

      if (k.le.0 .or. n.le.0) return

      call second(t2)

      do itry=1,NTRY
         normvnew_0 = normvnew         
         if ( iflag.eq.1 ) then
            call dcgs(n,k,V,ldv,vnew,index,work)
         else
            call dmgs(n,k,V,ldv,vnew,index)
         endif
         ndot = ndot + k
         normvnew = pdnrm2(n,vnew,1)
         if (normvnew.gt.alpha*normvnew_0) goto 9999
      enddo
      normvnew = zero
c     vnew is numerically in span(V) => return vnew = (0,0,...,0)^T
      call pdzero(n,vnew,1)
 9999 call second(t3)
      treorth = treorth + (t3-t2)
      nreorth = nreorth + 1
      return
      end
c
c****************************************************************************
c

      subroutine dcgs(n,k,V,ldv,vnew,index,work)

c     Block  Gram-Schmidt orthogonalization:
c     FOR i= 1:l
c         vnew = vnew - V(:,[s_i:e_i])*(V(:,[s_i:e_i])'*vnew)
c      
c     If l=1 and s_1=1 and e_1=k then this becomes classical Gram-Schmidt.

 
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      integer n,k,ldv,index(*)
      double precision V(ldv,*),vnew(*),work(*)
c     %------------%
c     | Parameters |
c     %------------%
      double precision one, zero      
      parameter(one = 1.0, zero = 0.0)
c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,j,p,q,l, ld
      integer tid, nt, cnk,st
      double precision ylocal(n)
c     %--------------------%
c     | External Functions |
c     %--------------------%
      external dgemv
#ifdef _OPENMP
      integer omp_get_thread_num, omp_get_num_threads
      external omp_get_thread_num, omp_get_num_threads
#endif

c The local variable ld was introduced to circumvent an apparent
c bug in Intel's OpenMP implementation.
      ld = ldv

c$OMP PARALLEL private(i,p,q,l,tid,nt,cnk,st,j,ylocal) 
c$OMP& firstprivate(ld) shared(ndot)
#ifdef _OPENMP
      tid = omp_get_thread_num()
      nt = omp_get_num_threads()
#else
      tid = 0
      nt = 1
#endif
      cnk = n/nt
      st = tid*cnk+1

      i=1
      do while(index(i).le.k .and. index(i).gt.0)
c
c     Select the next block of columns from V
c         
         p = index(i)
         q = index(i+1)
         l = q-p+1
         if (tid.eq.0) then
            ndot = ndot + l
         endif
c     Classical Gram-Schmidt: vnew = vnew - V(:,p:q)*(V(:,p:q)'*vnew)
         if (l.gt.0) then
            if (tid.eq.nt-1) then
               cnk = n-st+1
            endif
            call dgemv('T',cnk,l,one,V(st,p),ld,vnew(st),1,zero,
     c           ylocal,1)

            if (tid.eq.0) then
               do j=1,l
                  work(j) = ylocal(j)
               enddo
            endif
c$OMP BARRIER
            if (.not.tid.eq.0) then
c$OMP CRITICAL
               do j=1,l
                  work(j) = work(j) + ylocal(j)
               enddo
c$OMP END CRITICAL
            endif
c$OMP BARRIER
            call dgemv('N',cnk,l,-one,V(st,p),ld,work,1,zero,
     c           ylocal,1)
            do j=1,cnk
               vnew(st+j-1) = vnew(st+j-1) + ylocal(j)
            enddo
         endif
         i = i+2
      enddo
c$OMP END PARALLEL
      end      

