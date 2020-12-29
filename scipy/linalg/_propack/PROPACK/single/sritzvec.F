      subroutine sritzvec(which,jobu,jobv,m,n,k,dim,D,E,S,U,ldu,
     c     V,ldv,work,in_lwrk,iwork)

c
c     DRITZVEC: Compute Ritz-vectors corresponding to the K largest 
c               or smallest (depending on the value of WHICH) 
c               Ritz-values for A from the Lanczos bidiagonalization 
c               A*V_{dim} = U_{dim+1}*B_{dim}.
c
c     Parameters:
c
c     WHICH: CHARACTER*1. Decides which singular triplets to compute. 
c            If WHICH.EQ.'L' then compute triplets corresponding to the K
c            largest singular values. 
c            If WHICH.EQ.'S' then compute triplets corresponding to the K
c            smallest singular values. 
c     JOBU: CHARACTER*1. If JOBU.EQ.'Y' then compute the left singular vectors.
c           Otherwise the array U is not touched.
c     JOBV: CHARACTER*1. If JOBV.EQ.'Y' then compute the right singular 
c           vectors. Otherwise the array V is not touched.
c     M:    INTEGER. Number of rows of A.
c     N:    INTEGER. Number of columns of A.
c     K:    INTEGER. Number of desired singular triplets. K <= MIN(DIM,M,N)
c     DIM:  INTEGER. Dimension of the Krylov subspace.
c     D(DIM): REAL array. Contains the diagonal of B.
c     E(DIM): REAL array. Contains the first sub-diagonal of B.
c     S(K): REAL array. On return S contains approximation
c               to the K largest or smallest (depending on the 
c               value of WHICH) singular values of A.
c     U(LDU,DIM+1): REAL array. On return the first K columns of U
c               will contain approximations to the left singular vectors 
c               corresponding to the K largest or smallest (depending on the 
c               value of WHICH)  singular values of A.
c               On entry the first column of U contains the starting vector
c               for the Lanczos bidiagonalization. A random starting vector
c               is used if U is zero.
c     LDU: INTEGER. Leading dimension of the array U. LDV >= M.
c     V(LDV,DIM): REAL array. On return the first K columns of V
c               will contain approximations to the right singular vectors 
c               corresponding to the K largest or smallest (depending on the 
c               value of WHICH) singular values of A.
c     LDV: INTEGER. Leading dimension of the array V. LDV >= N.
c     WORK(LWORK): REAL array. Workspace of dimension LWORK.
c     IN_LWORK: INTEGER. Dimension of WORK. 
c            LWORK should be at least 3*DIM**2 + 
c            MAX(3*DIM**2+4*DIM+4, NB*MAX(M,N)), where NB>0 is a block 
c            size, which determines how large a fraction of the work in
c            setting up the singular vectors is done using fast BLAS-3 
c            operations. NB should probably be at least 32 to achieve good
c            performance.
c     IWORK(8*DIM): INTEGER array. Integer workspace of dimension >= 8*DIM. 
c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c

      
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      character*1 which, jobu,jobv
      integer m,n,k,dim,ldu,ldv,in_lwrk,iwork(*)
      real U(ldu,*),V(ldv,*),D(*),E(*),S(*),work(*)

c     %------------%
c     | Parameters |
c     %------------%
      real one, zero
      parameter(one = 1.0, zero = 0.0)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer lwrk, mstart
      integer ip,iqt,imt,iwrk,id(1),info
      real c1,c2,dd(1)
      real t0,t1
      integer st,cnk,wst,wcnk, tid, nt

c     %--------------------%
c     | External Functions |
c     %--------------------%
      logical lsame
      external lsame, sbdqr, sbdsdc, sgemm_ovwr_left
#ifdef _OPENMP
      integer omp_get_thread_num, omp_get_num_threads
      external omp_get_thread_num, omp_get_num_threads
#endif
    
c-------------------- Here begins executable code ---------------------

c     %----------------------------------------------------------------------%
c     | The bidiagonal SVD is computed in a two-stage procedure:             
c     |                                                   
c     | 1. Compute a QR-factorization M^T*B = [R; 0] of the (k+1)-by-k lower 
c     |    bidiagonal matrix B.
c     | 2. Compute the SVD of the k-by-k upper bidiagonal matrix 
c     |    R = P*S*Q^T. The SVD of B is then (M*P)*S*Q^T.
c     %----------------------------------------------------------------------%


      call second(t0)

c     %-----------------------------------------%
c     | Set pointers into workspace array
c     %-----------------------------------------%
      iwrk = 1
      lwrk = in_lwrk
      imt = 1
      iqt = imt + (dim+1)**2 
      ip = iqt + dim**2
      iwrk = ip + dim**2
      lwrk = lwrk - iwrk + 1  
c     %-----------------------------------------%
c     | Compute QR-factorization
c     |   B = M * [R; 0]
c     %-----------------------------------------%
      call sbdqr((dim.eq.min(m,n)),jobu,dim,D,E,c1,c2,work(imt),dim+1)

c     %-----------------------------------------%
c     | Compute SVD of R:
c     |      R = P * S * Q^T, 
c     | using the Divide-and-conquer SVD
c     %-----------------------------------------%
      call sbdsdc('u','I',dim,D,E,work(ip),dim,work(iqt),dim,dd,id,
     c     work(iwrk),iwork,info)

c     %-----------------------------------------%
c     | Compute left singular vectors for B
c     |    X = P^T * M^T 
c     %-----------------------------------------%
      call sgemm_ovwr('t',dim,dim+1,dim,1.0,work(ip),dim,0.0,
     c     work(imt),dim+1,work(iwrk),lwrk)
      
      
      if (lsame(jobu,'y')) then
c     %-----------------------------------------%
c     | Form left Ritz-vectors                  |
c     |   U = U * X^T                           |
c     %-----------------------------------------%
         if (lsame(which,'s')) then
            mstart = dim-k+1
         else
            mstart = 1
         endif
c$OMP PARALLEL private(tid,nt,cnk,st,wcnk,wst) 
#ifdef _OPENMP
         tid = omp_get_thread_num()
         nt = omp_get_num_threads()
#else
         tid = 0
         nt = 1
#endif
         wcnk = lwrk/nt
         wst = tid*wcnk+1
         cnk = m/nt
         st = tid*cnk+1
         if (tid.eq.nt-1) then
            wcnk = lwrk-wst+1
            cnk = m-st+1
         endif
         call sgemm_ovwr_left('t',cnk,k,dim+1,1.0,U(st,1),
     c        ldu,0.0,work(imt+mstart-1),dim+1,work(iwrk+wst-1),wcnk)
c$OMP END PARALLEL 
      endif

      if (lsame(jobv,'y')) then
         if (lsame(which,'s')) then
            mstart = dim-k+1
         else
            mstart = 1
         endif
c     %-----------------------------------------%
c     | Form right Ritz-vectors
c     |   V = V * Q
c     %-----------------------------------------%
c$OMP PARALLEL private(tid,nt,cnk,st,wcnk,wst) 
#ifdef _OPENMP
         tid = omp_get_thread_num()
         nt = omp_get_num_threads()
#else
         tid = 0
         nt = 1
#endif
         wcnk = lwrk/nt
         wst = tid*wcnk+1
         cnk = n/nt
         st = tid*cnk+1
         if (tid.eq.nt-1) then
            wcnk = lwrk-wst+1
            cnk = n-st+1
         endif
         call sgemm_ovwr_left('t',cnk,k,dim,1.0,V(st,1),ldv,0.0,
     c        work(iqt+mstart-1),dim,work(iwrk+wst-1),wcnk)
c$OMP END PARALLEL 
      endif    
      
      call  second(t1)           
      tritzvec = t1-t0
      end
      
