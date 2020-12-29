      subroutine slanbpro( m, n, k0, k, APROD, U, ldu, V, ldv, B, ldb,
     c     rnorm, doption, ioption, work, iwork, dparm, iparm, ierr)

c     DLANBPRO: Computes K steps of the Lanczos bidiagonalization (LBD) 
c     algorithm with partial reorthogonalization (BPRO) with M-by-1 starting 
c     vector U(:,k0+1), producing a lower bidiagonal K+1-by-K matrix B_k, an 
c     N-by-K matrix V_k, an M-by-K+1 matrix U_{k+1} such that
c           A*V_k = U_{k+1}*B_k
c     Partial reorthogonalization is used to keep the columns of V_K and U_k
c     semiorthogonal to a level prescribed in DOPTION(1), i.e.
c           MAX(DIAG((EYE(K) - V_K'*V_K))) <= DOPTION(1)
c     and
c           MAX(DIAG((EYE(K) - U_K'*U_K))) <= DOPTION(1).
c
c     If K0>0 and K>K0 an existing K0-step LBD of stored in U, V and B is 
c     extended to a K-step LBD.
c
c     Parameters:
c
c     M: INTEGER. Number of rows of A.
c     N: INTEGER. Number of columns of A.
c     K0: INTEGER. The dimension of the previously computed Lanczos
c                  bidiagonalization stored in U, V, and B.
c     K: INTEGER. On entry: The desired dimension of the Lanczos 
c           bidiagonalization. On exit: the actual size of the LBD computed.
c           This can be smaller than the input value if an invariant subspace
c           is computed.
c     APROD: Subroutine defining the linear operator A. 
c            APROD should be of the form:
c
c           SUBROUTINE DAPROD(TRANSA,M,N,X,Y,DPARM,IPARM)
c           CHARACTER*1 TRANSA
c           INTEGER M,N,IPARM(*)
c           REAL X(*),Y(*),DPARM(*)
c
c           If TRANSA.EQ.'N' then the function should compute the matrix-vector
c           product Y = A * X.
c           If TRANSA.EQ.'T' then the function should compute the matrix-vector
c           product Y = A^T * X.
c           The arrays IPARM and DPARM are a means to pass user supplied
c           data to APROD without the use of common blocks.
c     U(LDU,K+1): REAL array. On return the first K+1 columns of
c               U will contain the left Lanczos vectors.
c               On entry: 
c                  If K0==0 the first column of U contains the starting
c                  vector for the Lanczos bidiagonalization. A random 
c                  starting vector is used if the first column is U is zero.
c                  If K0>0 the first K0+1 columns of U are assumed to 
c                  contain the first K0+1 left Lanczos vectors of an
c                  existing LBD.
c
c     LDU: INTEGER. Leading dimension of the array U. LDU >= M.
c     V(LDV,K): REAL array. On return the first K columns of
c               V will contain the right Lanczos vectors.
c               On entry: 
c                  If K0>0 the first K0 columns of V are assumed to 
c                  contain the first K0 right Lanczos vectors of an
c                  existing LBD.
c     LDV: INTEGER. Leading dimension of the array V. LDV >= N.
c     B(K,2): REAL array. On return the first columns of
c               B will contain the K diagonal elements of B_k, and
c               the second column of B will contain the K elements
c               of the first sub-diagonal of B_k. 
c     LDB: INTEGER. Leading dimension of the array B. LDB >= K.
c     RNORM: REAL. On entry RNORM must contain the norm of
c               the K0+1st column of U.
c               On exit RNORM contains the value of the (K+1,K) element
c               of B_k.
c     DOPTION: REAL array. 
c        doption(1) = delta. Level of orthogonality to maintain among
c          Lanczos vectors.
c        doption(2) = eta. During reorthogonalization, all vectors with
c          with components larger than eta along the latest Lanczos vector
c          will be purged.
c        doption(3) = anorm. Estimate of || A ||.
c     IOPTION: INTEGER array. 
c        ioption(1) = CGS.  If CGS.EQ.1 then reorthogonalization is done
c          using iterated classical GRAM-SCHMIDT. IF CGS.EQ.0 then 
c          reorthogonalization is done using iterated modified Gram-Schmidt.
c        ioption(2) = ELR. If ELR.EQ.1 then extended local orthogonality is
c          enforced among u_{k}, u_{k+1} and v_{k} and v_{k+1} respectively.
c     WORK(LWORK): REAL array of dimension LWORK. Workspace.
c     LWORK: INTEGER. LWORK >= 2*(m+n+k+1)
c     IWORK(2*K+1): INTEGER ARRAY. Integer workspace.
c     DPARM: REAL array. Array used for passing data to the APROD
c         function.   
c     IPARM: INTEGER array. Array used for passing data to the APROD
c         function.   
c     IERR: INTEGER. Error status code.
c         IERR < 0  : An invariant subspace of dimension -J was found.
c         IERR == 0 : The computation succeeded.
c         IERR > 0  : The computation succeeded, but the algorithm 
c                     came close to computing an invariant subspace after 
c                     IERR steps. In a previous version this would have caused 
c                     the algorithm to switch to full reorthogonalization 
c                     after IERR steps, but that is no longer the case. 
c                     It is probably safe to ignore.
c
c     (C) Rasmus Munk Larsen, Stanford University, 1999, 2004
c


c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      integer m, n, k0, k, ldb, ldu, ldv, ierr
      integer ioption(*), iwork(*), iparm(*)
      real rnorm,B(ldb,*), doption(*), work(*), dparm(*)
      real U(ldu,*),V(ldv,*)
      external APROD

c     %------------%
c     | Parameters |
c     %------------%
      real one, zero, FUDGE, kappa
      parameter(one = 1.0, zero = 0.0, FUDGE = 1.01)
      parameter (kappa = 0.717)
      
c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,j,inu,imu,is,iidx,j0
      real eps,eps34,epsn2,epsn,sfmin,delta,eta,anorm,s
      real mumax,numax,alpha,beta,a1,b1,amax,anormest
      logical force_reorth,full_reorth
      real t1,t2,t3

c     %----------------------%
c     | External Subroutines |
c     %----------------------%
      external sgetu0,sreorth,ssafescal,szero,izero,pscopy,psaxpy
      external scompute_int,supdate_nu,supdate_mu
      
c     %--------------------%
c     | External Functions |
c     %--------------------%
      real slamch,psnrm2,psdot,slapy2
      external psnrm2,psdot
      external slamch,slapy2

c-------------------- Here begins executable code ---------------------
      call second(t1)

c     %---------------------------------%
c     | Set machine dependent constants |
c     %---------------------------------%
      eps = slamch('e')
      eps34 = eps**(3.0/4.0)
      epsn = real(max(m,n))*eps
      epsn2 = sqrt(real(max(m,n)))*eps
      sfmin = slamch('s')
      
c     %------------------------%
c     | Set default parameters |
c     %------------------------%
      if (doption(1).lt.zero) then
         delta = sqrt(eps/k)
      else
         delta = doption(1)
      endif
      if (doption(2).lt.zero) then
         eta = eps34/sqrt(real(k))
      else
         eta = doption(2)
      endif
      if (delta.le.eta .or. delta.eq.zero) then
         full_reorth = .true.
      else
         full_reorth = .false.
      endif
      if (doption(3).gt.zero) then
         anorm = doption(3)
      else if (k0.gt.0) then
         anorm = slapy2(B(1,1),B(1,2))
         if (anorm.le.zero) then
            ierr = -1
            goto 9999         
         endif
      else
         anorm = zero
      endif
      ierr = 0

c     %---------------------%
c     | Get starting vector |
c     %---------------------%
      if (rnorm .eq. zero) then
         call sgetu0('n',m, n, k0, 3, U(1,k0+1), rnorm, U,
     c        ldu, aprod, dparm, iparm, ierr, ioption(1), anormest,
     c        work)
         anorm = max(anorm,anormest)
      endif

c     %------------------------------%
c     | Set pointers into work array |
c     %------------------------------%
      imu = 1
      inu = imu+k+1
      is = inu+k+1
      iidx = 1
      call pszero(max(m,n)+2*k+2,work,1)
      call pizero(2*k+1,iwork,1)

c     %---------------------------%
c     | Prepare Lanczos iteration |
c     %---------------------------%
      if (k0.eq.0) then
         amax = zero
         alpha = zero
         beta = rnorm
         force_reorth = .false.

c     %---------------------------------------------------%
c     | Compute ||A x|| / ||x|| for a random vector x     |
c     | to make it less likely that ||A|| is grossly      |
c     | underestimated at the beginning of the iteration. |
c     %---------------------------------------------------%
         if (n.gt.m) then            
            call sgetu0('n',m,n,0,1,work(is),s,U,ldu,aprod,dparm,iparm,
     c           ierr,ioption(1),anormest,work(is+m))
         else            
            call sgetu0('t',m,n,0,1,work(is),s,V,ldv,aprod,dparm,iparm,
     c           ierr,ioption(1),anormest,work(is+n))
         endif
         anorm = max(anorm,FUDGE*anormest)
         j0 = 1
         if (beta.ne.zero) then
            call ssafescal(m,beta,U(1,1))
         endif
         work(imu) = one
         work(inu) = one
      else
         force_reorth = .true.
         alpha = B(k0,1)
         beta = rnorm
         if (k0.lt.k .and. beta*delta.lt.anorm*eps) then
            full_reorth = .true.
            ierr = k0
         endif

         iwork(iidx) = 1
         iwork(iidx+1) = k0
         iwork(iidx+2) = k0+1

         call psscal(m,rnorm,U(1,k0+1),1)
         call second(t2)         
         call sreorth(m,k0,U,ldu,U(1,k0+1),rnorm,iwork(iidx),kappa,
     c        work(is),ioption(1))
         call second(t3)
         treorthu = treorthu+(t3-t2)         
         call ssafescal(m,rnorm,U(1,k0+1))
         call sset_mu(k0,work(imu),iwork(iidx),epsn2)         
         call sset_mu(k0,work(inu),iwork(iidx),epsn2)         
         beta = rnorm

c     %--------------------------------------%
c     | Estimate ||B||_2^2 as ||B^T * B||_1  |   
c     %--------------------------------------%
         B(k0,2) = beta
         amax = zero
         do j=1,k0
            amax = max(amax,B(j,1),B(j,2))
            if (j.eq.1) then
               anorm = max(anorm,FUDGE*alpha)
            else if (j.eq.2) then
               a1 = B(1,2)/amax
               a1 = FUDGE*amax*sqrt((B(1,1)/amax)**2 + a1**2 +
     c              B(2,1)/amax*a1)
               anorm = max(anorm,a1)
            else
               a1 = B(j-1,1)/amax
               b1 = B(j-1,2)/amax 
               a1 = FUDGE*amax*sqrt( a1**2 + b1**2 +
     c              a1*B(j-2,2)/amax + B(j,1)/amax*b1)
               anorm = max(anorm,a1)
            endif
         enddo
         j0 = k0+1
      endif
      numax = zero
      mumax = zero
      
c     %-------------------------------------------%
c     | Start Lanczos bidiagonalization iteration |
c     %-------------------------------------------%            

      do j=j0,k

c     %---------------------------------------------%
c     | alpha_{j} v_{j} = A'*u_{j} - beta_{j} v_{j} |
c     %---------------------------------------------%
         call second(t2)
         call aprod('t',m,n,U(1,j),V(1,j),dparm,iparm)
         call second(t3)
         tmvopx = tmvopx + (t3-t2)
         nopx = nopx+1

         if (j.eq.1) then
            alpha = psnrm2(n,V(1,j),1)
            anorm = max(anorm,FUDGE*alpha)
         else            
            call psaxpy(n,-beta,V(1,j-1),1,V(1,j),1)
            alpha = psnrm2(n,V(1,j),1)

c     %------------------------------------%
c     | Extended local reorthogonalization |
c     %------------------------------------%
            call second(t2)
            if (j.gt.1 .and. ioption(2).gt.0 .and.
     c           alpha.lt.kappa*beta) then
               do i=1,ioption(2)
                  s = psdot(n,V(1,j-1),1,V(1,j),1)
                  call psaxpy(n,-s,V(1,j-1),1,V(1,j),1)
                  if (beta .ne. zero) then
                     beta = beta + s
                     B(j-1,2) = beta
                  endif
                  s = psnrm2(n,V(1,j),1)
                  if (s .ge. kappa*alpha) goto 10                  
                  alpha = s
               enddo
 10            work(inu+j-2) = eps
               alpha = s
            endif
            call second(t3)
            telrv = telrv + (t3-t2)

            B(j,1) = alpha
            amax = max(amax,alpha)

c     %----------------------------%
c     | Update estimate of ||A||_2 |         
c     %----------------------------%
            if (j.eq.2) then
               a1 = B(1,2)/amax
               a1 = FUDGE*amax*sqrt((B(1,1)/amax)**2 + a1**2 +
     c              B(2,1)/amax*a1)
            else
               a1 = B(j-1,1)/amax
               b1 = B(j-1,2)/amax 
               a1 = FUDGE*amax*sqrt( a1**2 + b1**2 +
     c              a1*B(j-2,2)/amax + B(j,1)/amax*b1)
            endif
            anorm = max(anorm,a1)
         endif

c     %--------------------------%
c     | Update the nu recurrence |
c     %--------------------------%
         if (.not.full_reorth .and. alpha.ne.zero) then
            call supdate_nu(numax,work(imu),work(inu),j,
     c           B(1,1),B(1,2),anorm,epsn2)
         endif

c     %------------------------------%
c     | Reorthogonalize if necessary |
c     %------------------------------%
         if ( (full_reorth .or. numax.gt.delta .or. force_reorth) 
     c        .and. alpha.ne.zero) then
            if (full_reorth .or. eta.eq.zero) then
               iwork(iidx) = 1
               iwork(iidx+1) = j-1
               iwork(iidx+2) = j
            else if (.not. force_reorth) then
               call scompute_int(work(inu),j-1,delta,eta,iwork(iidx))
            endif
            call second(t2)
            call sreorth(n,j-1,V,ldv,V(1,j),alpha,iwork(iidx),
     c           kappa,work(is),ioption(1))
            call second(t3)
            treorthv = treorthv+(t3-t2)

            call sset_mu(j-1,work(inu),iwork(iidx),eps)
            numax = eta
            if (force_reorth) then
               force_reorth = .false.
            else
               force_reorth = .true.
            endif
         endif

c     %-----------------------------------------------%
c     | Check whether an invariant subspace was found |
c     %-----------------------------------------------% 
         if (alpha .lt. anorm*epsn .and. j.lt.k) then
            rnorm = alpha
            alpha = zero
c     %------------------------------------------------%
c     | Try to build an orthogonal subspace, starting  |
c     | with a random vector.                          |
c     %------------------------------------------------%
            call sgetu0('t', m, n, j-1, 3, V(1,j), alpha, V, ldv,
     c           aprod, dparm, iparm, ierr,ioption(1),anormest,
     c           work(is))
            if (alpha .eq. zero) then
c     %------------------------------------------------%
c     | We failed to generate a new random vector      |
c     | in span(A^T) orthogonal to span(V(:,1:j-1)).   |
c     | Most likely span(V(:,1:j-1)) is an invariant   |
c     | subspace.                                      |
c     %------------------------------------------------%
               k = j-1
               ierr = -j
               goto 9999
            else
c     %-------------------------------------------------%
c     | We have managed to generate a random vector     |
c     | in span(A^T) orthogonal to V(:,1:j-1), so we    |
c     | can continue the LBD and "deflate" the subspace |
c     | by setting alpha_{j} = 0.                       |
c     %-------------------------------------------------%
              call ssafescal(n,alpha,V(1,j))
               alpha = zero
               force_reorth = .true.
               if (delta.gt.zero) then
                  full_reorth = .false.
               endif
            endif
         else if (j.gt.1 .and. .not. full_reorth .and. j.lt.k .and.
     c           (delta*alpha .lt. anorm*eps)) then
            ierr = j
         endif            
         B(j,1) = alpha

         if (alpha.ne.zero) then
            call ssafescal(n,alpha,V(1,j))
         endif


c     %------------------------------------------------%
c     | beta_{j+1} u_{j+1} = A*v_{j} - alpha_{j} u_{j} |
c     %------------------------------------------------%
         call second(t2)
         call aprod('n',m,n,V(1,j),U(1,j+1),dparm,iparm)
         call second(t3)
         tmvopx = tmvopx + (t3-t2)
         nopx = nopx+1

         call psaxpy(m,-alpha,U(1,j),1,U(1,j+1),1)
         beta = psnrm2(m,U(1,j+1),1)
         
c     %------------------------------------%
c     | Extended local reorthogonalization |
c     %------------------------------------%
         call second(t2)
         if (ioption(2).gt.0 .and. beta.lt.kappa*alpha) then
            do i=1,ioption(2)
               s = psdot(m,U(1,j),1,U(1,j+1),1)
               call psaxpy(m,-s,U(1,j),1,U(1,j+1),1)
               if (alpha .ne. zero) then
                  alpha = alpha + s
                  B(j,1) = alpha
               endif
               s = psnrm2(m,U(1,j+1),1)       
               if (s .ge. kappa*beta) goto 20
               beta = s
            enddo
 20         work(imu+j-1) = eps
            beta = s
         endif
         call second(t3)
         telru = telru + (t3-t2)

         B(j,2) = beta
         amax = max(amax,beta)
      
c     %----------------------------%
c     | Update estimate of ||A||_2 |         
c     %----------------------------%
         if (j.le.1) then
            a1 = slapy2(B(1,1), B(1,2))
         else
            a1 = B(j,1)/amax
            a1 = amax*sqrt(a1**2 + (B(j,2)/amax)**2 +
     c           a1*B(j-1,2)/amax)
         endif      
         anorm = max(anorm,a1)
         
c     %--------------------------%
c     | Update the mu recurrence |
c     %--------------------------%
         if (.not.full_reorth .and. beta.ne.zero) then
            call supdate_mu(mumax,work(imu),work(inu),j,B(1,1),
     c           B(1,2),anorm,epsn2)
         endif

c     %--------------------------------------%
c     | Reorthogonalize u_{j+1} if necessary |
c     %--------------------------------------%
         if ( (full_reorth .or. mumax.gt.delta .or. force_reorth) 
     c        .and. beta.ne.zero) then
            if (full_reorth .or. eta.eq.zero) then
               iwork(iidx) = 1
               iwork(iidx+1) = j
               iwork(iidx+2) = j+1
            else if (.not. force_reorth) then
               call scompute_int(work(imu),j,delta,eta,iwork(iidx))
            else
               do i=1,2*j+1
                  if (iwork(iidx+i-1).eq.j) then
                     iwork(iidx+i-1) = j+1
                     goto 25
                  endif
               enddo
            endif
 
 25         call second(t2)        
            call sreorth(m,j,U,ldu,U(1,j+1),beta,iwork(iidx),
     c              kappa, work(is),ioption(1))
            call second(t3)
            treorthu = treorthu+(t3-t2)

            call sset_mu(j,work(imu),iwork(iidx),eps)
            mumax = eta
            if (force_reorth) then
               force_reorth = .false.
            else
               force_reorth = .true.
            endif
         endif

c     %-----------------------------------------------%
c     | Check whether an invariant subspace was found |
c     %-----------------------------------------------%
         if (beta .lt. anorm*epsn .and. j.lt.k) then
            rnorm = beta
            beta = zero
c     %-----------------------------------------------%
c     | Try to build an orthogonal subspace, starting |
c     | with a random vector.                         |
c     %-----------------------------------------------%
            call sgetu0('n', m, n, j, 3, U(1,j+1), beta, U, ldu, aprod,
     c           dparm, iparm, ierr,ioption(1),anormest,work(is))
            if (beta .eq. zero) then
c     %-----------------------------------------------%
c     | We failed to generate a new random vector     |
c     | in span(A) orthogonal to span(U(:,1:j)).      |
c     | Most likely span(U(:,1:j)) is an invariant    |
c     | subspace.                                     |
c     %-----------------------------------------------%
               k = j
               ierr = -j
               goto 9999
            else
c     %------------------------------------------------%
c     | We have managed to generate a random vector    |
c     | in span(A) orthogonal to U(:,1:j), so we can   |
c     | continue the LBD and "deflate" the subspace by |
c     | setting beta_{j+1} = 0.                        |
c     %------------------------------------------------%
               call ssafescal(n,beta,U(1,j+1))
               beta = zero
               force_reorth = .true.
               if (delta .gt. zero) then
                  full_reorth = .false.
               endif
            endif
         else if (.not.full_reorth .and. j.lt.k .and. 
     c           (delta*beta .lt. anorm*eps)) then
            ierr = j
         endif            

         B(j,2) = beta
         if (beta.ne.zero .and. beta.ne.one) then
            call ssafescal(m,beta,U(1,j+1))
         endif
         rnorm = beta
         call second(t2)
      enddo
 9999 doption(3) = anorm      
      call second(t2)
      tlanbpro = tlanbpro + (t2-t1)
      return
      end

c
c**********************************************************************
c

      subroutine sset_mu(k,mu,index,val)
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      integer k,index(*)
      real mu(*), val

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,j,p,q

      i=1
      do while(index(i).le.k .and. index(i).gt.0)
         p = index(i)
         q = index(i+1)
         do j=p,q
            mu(j) = val
         enddo
         i = i+2
      enddo      
      end
c
c**********************************************************************
c
      subroutine scompute_int(mu,j,delta,eta,index)
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      integer j,index(*)
      real mu(*)
      real delta,eta

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,k,s,ip
      real t1,t2

      call second(t1)
      if (delta.lt.eta) then
         write (*,*) 'Warning delta<eta in scompute_int'
         return
      endif

      ip = 0
      index(1) = 0
      i=0
      do while(i.lt.j)
c     find the next mu(k), k>i where abs(mu(k)) > delta
         do k=i+1,j
            if (abs(mu(k)).gt.delta) goto 10
         enddo         
         goto 40
c     find smallest i<k such that for all j=i,..,k, m(j) >= eta
 10      do s=k,max(i,1),-1
             if (abs(mu(s)).lt.eta) goto 20
         enddo
 20      ip= ip+1
         index(ip) = s+1
         do i=s+1,j            
            if (abs(mu(i)).lt.eta) goto 30
         enddo         
 30      ip= ip+1
         index(ip) = i-1
      enddo
 40   ip = ip+1
      index(ip) = j+1
      call second(t2)
      tintv = tintv + (t2-t1)
      end
c
c**********************************************************************
c
      subroutine supdate_mu(mumax,mu,nu,j,alpha,beta,anorm,eps1)
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      integer j
      real mumax,eps1,anorm
      real mu(*),nu(*),alpha(*),beta(*)

c     %------------%
c     | Parameters |
c     %------------%
      real one, zero, FUDGE
      parameter(one = 1.0, zero = 0.0, FUDGE = 1.01)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      real d
      integer k
      real t1,t2
      
c     %--------------------%
c     | External Functions |
c     %--------------------%
      real slapy2
      external slapy2

      call second(t1)
      if (j.eq.1) then
         d = eps1*(slapy2(alpha(j), beta(j)) + alpha(1)) + eps1*anorm
         mu(1) = eps1/beta(1)
         mumax = abs(mu(1))
      else
         mu(1) = alpha(1)*nu(1)-alpha(j)*mu(1)
         d = eps1*(slapy2(alpha(j), beta(j)) + alpha(1)) + eps1*anorm
         mu(1) = (mu(1) + sign(d,mu(1))) / beta(j)
         mumax = abs(mu(1))
         do k=2,j-1
            mu(k) = alpha(k)*nu(k) +beta(k-1)*nu(k-1)-alpha(j)*mu(k)
            d = eps1*(slapy2(alpha(j), beta(j)) + 
     c           slapy2(alpha(k), beta(k-1))) + eps1*anorm
            mu(k) = (mu(k) + sign(d,mu(k))) / beta(j)
            mumax = max(mumax,abs(mu(k)))
         enddo
         mu(j) = beta(j-1)*nu(j-1)
         d = eps1*(slapy2(alpha(j), beta(j)) + 
     c        slapy2(alpha(j), beta(j-1))) + eps1*anorm
         mu(j) = (mu(j) + sign(d,mu(j))) / beta(j)
         mumax = max(mumax,abs(mu(j)))
      endif
      mu(j+1) = one
      call second(t2)
      tupdmu = tupdmu + (t2-t1)
      end
c
c**********************************************************************
c

      subroutine supdate_nu(numax,mu,nu,j,alpha,beta,anorm,eps1)
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      integer j
      real numax,eps1,anorm
      real mu(*),nu(*),alpha(*),beta(*)

c     %------------%
c     | Parameters |
c     %------------%
      real one, zero, FUDGE
      parameter(one = 1.0, zero = 0.0, FUDGE = 1.01)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      real d
      integer k
      real t1,t2

c     %--------------------%
c     | External Functions |
c     %--------------------%
      real slapy2
      external slapy2

      call second(t1)
      if (j.gt.1) then
         numax = zero
         do k=1,j-1
            nu(k) = beta(k)*mu(k+1) + alpha(k)*mu(k) -beta(j-1)*nu(k)
            d = eps1*(slapy2(alpha(k),beta(k)) +
     c           slapy2(alpha(j),beta(j-1))) + eps1*anorm
            nu(k) = (nu(k) + sign(d,nu(k))) / alpha(j)
            numax = max(numax,abs(nu(k)))
         enddo
         nu(j) = one
      endif
      call second(t2)
      tupdnu = tupdnu + (t2-t1)
      end
