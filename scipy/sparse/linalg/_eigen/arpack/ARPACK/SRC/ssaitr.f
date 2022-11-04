c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: ssaitr
c
c\Description: 
c  Reverse communication interface for applying NP additional steps to 
c  a K step symmetric Arnoldi factorization.
c
c  Input:  OP*V_{k}  -  V_{k}*H = r_{k}*e_{k}^T
c
c          with (V_{k}^T)*B*V_{k} = I, (V_{k}^T)*B*r_{k} = 0.
c
c  Output: OP*V_{k+p}  -  V_{k+p}*H = r_{k+p}*e_{k+p}^T
c
c          with (V_{k+p}^T)*B*V_{k+p} = I, (V_{k+p}^T)*B*r_{k+p} = 0.
c
c  where OP and B are as in ssaupd.  The B-norm of r_{k+p} is also
c  computed and returned.
c
c\Usage:
c  call ssaitr
c     ( IDO, BMAT, N, K, NP, MODE, RESID, RNORM, V, LDV, H, LDH, 
c       IPNTR, WORKD, INFO )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y.
c                    This is for the restart phase to force the new
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y,
c                    IPNTR(3) is the pointer into WORK for B * X.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y.
c          IDO = 99: done
c          -------------------------------------------------------------
c          When the routine is used in the "shift-and-invert" mode, the
c          vector B * Q is already available and does not need to be
c          recomputed in forming OP * Q.
c
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of matrix B that defines the
c          semi-inner product for the operator OP.  See ssaupd.
c          B = 'I' -> standard eigenvalue problem A*x = lambda*x
c          B = 'G' -> generalized eigenvalue problem A*x = lambda*M*x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  K       Integer.  (INPUT)
c          Current order of H and the number of columns of V.
c
c  NP      Integer.  (INPUT)
c          Number of additional Arnoldi steps to take.
c
c  MODE    Integer.  (INPUT)
c          Signifies which form for "OP". If MODE=2 then
c          a reduction in the number of B matrix vector multiplies
c          is possible since the B-norm of OP*x is equivalent to
c          the inv(B)-norm of A*x.
c
c  RESID   Real array of length N.  (INPUT/OUTPUT)
c          On INPUT:  RESID contains the residual vector r_{k}.
c          On OUTPUT: RESID contains the residual vector r_{k+p}.
c
c  RNORM   Real scalar.  (INPUT/OUTPUT)
c          On INPUT the B-norm of r_{k}.
c          On OUTPUT the B-norm of the updated residual r_{k+p}.
c
c  V       Real N by K+NP array.  (INPUT/OUTPUT)
c          On INPUT:  V contains the Arnoldi vectors in the first K 
c          columns.
c          On OUTPUT: V contains the new NP Arnoldi vectors in the next
c          NP columns.  The first K columns are unchanged.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling 
c          program.
c
c  H       Real (K+NP) by 2 array.  (INPUT/OUTPUT)
c          H is used to store the generated symmetric tridiagonal matrix
c          with the subdiagonal in the first column starting at H(2,1)
c          and the main diagonal in the second column.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling 
c          program.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c          Pointer to mark the starting locations in the WORK for 
c          vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X.
c          IPNTR(2): pointer to the current result vector Y.
c          IPNTR(3): pointer to the vector B * X when used in the 
c                    shift-and-invert mode.  X is the current operand.
c          -------------------------------------------------------------
c          
c  WORKD   Real work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The calling program should not 
c          use WORKD as temporary workspace during the iteration !!!!!!
c          On INPUT, WORKD(1:N) = B*RESID where RESID is associated
c          with the K step Arnoldi factorization. Used to save some 
c          computation at the first step. 
c          On OUTPUT, WORKD(1:N) = B*RESID where RESID is associated
c          with the K+NP step Arnoldi factorization.
c
c  INFO    Integer.  (OUTPUT)
c          = 0: Normal exit.
c          > 0: Size of an invariant subspace of OP is found that is
c               less than K + NP.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     sgetv0  ARPACK routine to generate the initial vector.
c     ivout   ARPACK utility routine that prints integers.
c     smout   ARPACK utility routine that prints matrices.
c     svout   ARPACK utility routine that prints vectors.
c     slamch  LAPACK routine that determines machine constants.
c     slascl  LAPACK routine for careful scaling of a matrix.
c     sgemv   Level 2 BLAS routine for matrix vector multiplication.
c     saxpy   Level 1 BLAS that computes a vector triad.
c     sscal   Level 1 BLAS that scales a vector.
c     scopy   Level 1 BLAS that copies one vector to another .
c     sdot    Level 1 BLAS that computes the scalar product of two vectors. 
c     snrm2   Level 1 BLAS that computes the norm of a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c 
c\Revision history:
c     xx/xx/93: Version ' 2.4'
c
c\SCCS Information: @(#) 
c FILE: saitr.F   SID: 2.6   DATE OF SID: 8/28/96   RELEASE: 2
c
c\Remarks
c  The algorithm implemented is:
c  
c  restart = .false.
c  Given V_{k} = [v_{1}, ..., v_{k}], r_{k}; 
c  r_{k} contains the initial residual vector even for k = 0;
c  Also assume that rnorm = || B*r_{k} || and B*r_{k} are already 
c  computed by the calling program.
c
c  betaj = rnorm ; p_{k+1} = B*r_{k} ;
c  For  j = k+1, ..., k+np  Do
c     1) if ( betaj < tol ) stop or restart depending on j.
c        if ( restart ) generate a new starting vector.
c     2) v_{j} = r(j-1)/betaj;  V_{j} = [V_{j-1}, v_{j}];  
c        p_{j} = p_{j}/betaj
c     3) r_{j} = OP*v_{j} where OP is defined as in ssaupd
c        For shift-invert mode p_{j} = B*v_{j} is already available.
c        wnorm = || OP*v_{j} ||
c     4) Compute the j-th step residual vector.
c        w_{j} =  V_{j}^T * B * OP * v_{j}
c        r_{j} =  OP*v_{j} - V_{j} * w_{j}
c        alphaj <- j-th component of w_{j}
c        rnorm = || r_{j} ||
c        betaj+1 = rnorm
c        If (rnorm > 0.717*wnorm) accept step and go back to 1)
c     5) Re-orthogonalization step:
c        s = V_{j}'*B*r_{j}
c        r_{j} = r_{j} - V_{j}*s;  rnorm1 = || r_{j} ||
c        alphaj = alphaj + s_{j};   
c     6) Iterative refinement step:
c        If (rnorm1 > 0.717*rnorm) then
c           rnorm = rnorm1
c           accept step and go back to 1)
c        Else
c           rnorm = rnorm1
c           If this is the first time in step 6), go to 5)
c           Else r_{j} lies in the span of V_{j} numerically.
c              Set r_{j} = 0 and rnorm = 0; go to 1)
c        EndIf 
c  End Do
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine ssaitr
     &   (ido, bmat, n, k, np, mode, resid, rnorm, v, ldv, h, ldh, 
     &    ipntr, workd, info)
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1
      integer    ido, info, k, ldh, ldv, n, mode, np
      Real
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(3)
      Real
     &           h(ldh,2), resid(n), v(ldv,k+np), workd(3*n)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Real
     &           one, zero
      parameter (one = 1.0E+0, zero = 0.0E+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      logical    first, orth1, orth2, rstart, step3, step4
      integer    i, ierr, ipj, irj, ivj, iter, itry, j, msglvl, 
     &           infol, jj
      Real
     &           rnorm1, wnorm, safmin, temp1
      save       orth1, orth2, rstart, step3, step4,
     &           ierr, ipj, irj, ivj, iter, itry, j, msglvl,
     &           rnorm1, safmin, wnorm
c
c     %-----------------------%
c     | Local Array Arguments | 
c     %-----------------------%
c
      Real
     &           xtemp(2)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   saxpy, scopy, sscal, sgemv, sgetv0, svout, smout,
     &           slascl, ivout, arscnd
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Real
     &           sdot, snrm2, slamch
      external   sdot, snrm2, slamch
c
c     %-----------------%
c     | Data statements |
c     %-----------------%
c
      data      first / .true. /
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (first) then
         first = .false.
c
c        %--------------------------------%
c        | safmin = safe minimum is such  |
c        | that 1/sfmin does not overflow |
c        %--------------------------------%
c
         safmin = slamch('safmin')
      end if
c
      if (ido .eq. 0) then
c 
c        %-------------------------------%
c        | Initialize timing statistics  |
c        | & message level for debugging |
c        %-------------------------------%
c
         call arscnd (t0)
         msglvl = msaitr
c 
c        %------------------------------%
c        | Initial call to this routine |
c        %------------------------------%
c
         info   = 0
         step3  = .false.
         step4  = .false.
         rstart = .false.
         orth1  = .false.
         orth2  = .false.
c 
c        %--------------------------------%
c        | Pointer to the current step of |
c        | the factorization to build     |
c        %--------------------------------%
c
         j      = k + 1
c 
c        %------------------------------------------%
c        | Pointers used for reverse communication  |
c        | when using WORKD.                        |
c        %------------------------------------------%
c
         ipj    = 1
         irj    = ipj   + n
         ivj    = irj   + n
      end if
c 
c     %-------------------------------------------------%
c     | When in reverse communication mode one of:      |
c     | STEP3, STEP4, ORTH1, ORTH2, RSTART              |
c     | will be .true.                                  |
c     | STEP3: return from computing OP*v_{j}.          |
c     | STEP4: return from computing B-norm of OP*v_{j} |
c     | ORTH1: return from computing B-norm of r_{j+1}  |
c     | ORTH2: return from computing B-norm of          |
c     |        correction to the residual vector.       |
c     | RSTART: return from OP computations needed by   |
c     |         sgetv0.                                 |
c     %-------------------------------------------------%
c
      if (step3)  go to 50
      if (step4)  go to 60
      if (orth1)  go to 70
      if (orth2)  go to 90
      if (rstart) go to 30
c
c     %------------------------------%
c     | Else this is the first step. |
c     %------------------------------%
c 
c     %--------------------------------------------------------------%
c     |                                                              |
c     |        A R N O L D I     I T E R A T I O N     L O O P       |
c     |                                                              |
c     | Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) |
c     %--------------------------------------------------------------%
c
 1000 continue
c
         if (msglvl .gt. 2) then
            call ivout (logfil, 1, j, ndigit, 
     &                  '_saitr: generating Arnoldi vector no.')
            call svout (logfil, 1, rnorm, ndigit, 
     &                  '_saitr: B-norm of the current residual =')
         end if
c 
c        %---------------------------------------------------------%
c        | Check for exact zero. Equivalent to determing whether a |
c        | j-step Arnoldi factorization is present.                |
c        %---------------------------------------------------------%
c
         if (rnorm .gt. zero) go to 40
c
c           %---------------------------------------------------%
c           | Invariant subspace found, generate a new starting |
c           | vector which is orthogonal to the current Arnoldi |
c           | basis and continue the iteration.                 |
c           %---------------------------------------------------%
c
            if (msglvl .gt. 0) then
               call ivout (logfil, 1, j, ndigit,
     &                     '_saitr: ****** restart at step ******')
            end if
c 
c           %---------------------------------------------%
c           | ITRY is the loop variable that controls the |
c           | maximum amount of times that a restart is   |
c           | attempted. NRSTRT is used by stat.h         |
c           %---------------------------------------------%
c
            nrstrt = nrstrt + 1
            itry   = 1
   20       continue
            rstart = .true.
            ido    = 0
   30       continue
c
c           %--------------------------------------%
c           | If in reverse communication mode and |
c           | RSTART = .true. flow returns here.   |
c           %--------------------------------------%
c
            call sgetv0 (ido, bmat, itry, .false., n, j, v, ldv, 
     &                   resid, rnorm, ipntr, workd, ierr)
            if (ido .ne. 99) go to 9000
            if (ierr .lt. 0) then
               itry = itry + 1
               if (itry .le. 3) go to 20
c
c              %------------------------------------------------%
c              | Give up after several restart attempts.        |
c              | Set INFO to the size of the invariant subspace |
c              | which spans OP and exit.                       |
c              %------------------------------------------------%
c
               info = j - 1
               call arscnd (t1)
               tsaitr = tsaitr + (t1 - t0)
               ido = 99
               go to 9000
            end if
c 
   40    continue
c
c        %---------------------------------------------------------%
c        | STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  |
c        | Note that p_{j} = B*r_{j-1}. In order to avoid overflow |
c        | when reciprocating a small RNORM, test against lower    |
c        | machine bound.                                          |
c        %---------------------------------------------------------%
c
         call scopy (n, resid, 1, v(1,j), 1)
         if (rnorm .ge. safmin) then
             temp1 = one / rnorm
             call sscal (n, temp1, v(1,j), 1)
             call sscal (n, temp1, workd(ipj), 1)
         else
c
c            %-----------------------------------------%
c            | To scale both v_{j} and p_{j} carefully |
c            | use LAPACK routine SLASCL               |
c            %-----------------------------------------%
c
             call slascl ('General', i, i, rnorm, one, n, 1, 
     &                    v(1,j), n, infol)
             call slascl ('General', i, i, rnorm, one, n, 1, 
     &                    workd(ipj), n, infol)
         end if
c 
c        %------------------------------------------------------%
c        | STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} |
c        | Note that this is not quite yet r_{j}. See STEP 4    |
c        %------------------------------------------------------%
c
         step3 = .true.
         nopx  = nopx + 1
         call arscnd (t2)
         call scopy (n, v(1,j), 1, workd(ivj), 1)
         ipntr(1) = ivj
         ipntr(2) = irj
         ipntr(3) = ipj
         ido = 1
c 
c        %-----------------------------------%
c        | Exit in order to compute OP*v_{j} |
c        %-----------------------------------%
c 
         go to 9000
   50    continue
c 
c        %-----------------------------------%
c        | Back from reverse communication;  |
c        | WORKD(IRJ:IRJ+N-1) := OP*v_{j}.   |
c        %-----------------------------------%
c
         call arscnd (t3)
         tmvopx = tmvopx + (t3 - t2)
c 
         step3 = .false.
c
c        %------------------------------------------%
c        | Put another copy of OP*v_{j} into RESID. |
c        %------------------------------------------%
c
         call scopy (n, workd(irj), 1, resid, 1)
c 
c        %-------------------------------------------%
c        | STEP 4:  Finish extending the symmetric   |
c        |          Arnoldi to length j. If MODE = 2 |
c        |          then B*OP = B*inv(B)*A = A and   |
c        |          we don't need to compute B*OP.   |
c        | NOTE: If MODE = 2 WORKD(IVJ:IVJ+N-1) is   |
c        | assumed to have A*v_{j}.                  |
c        %-------------------------------------------%
c
         if (mode .eq. 2) go to 65
         call arscnd (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            step4 = .true.
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c 
c           %-------------------------------------%
c           | Exit in order to compute B*OP*v_{j} |
c           %-------------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
              call scopy(n, resid, 1 , workd(ipj), 1)
         end if
   60    continue
c 
c        %-----------------------------------%
c        | Back from reverse communication;  |
c        | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j}. |
c        %-----------------------------------%
c
         if (bmat .eq. 'G') then
            call arscnd (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if 
c
         step4 = .false.
c
c        %-------------------------------------%
c        | The following is needed for STEP 5. |
c        | Compute the B-norm of OP*v_{j}.     |
c        %-------------------------------------%
c
   65    continue
         if (mode .eq. 2) then
c
c           %----------------------------------%
c           | Note that the B-norm of OP*v_{j} |
c           | is the inv(B)-norm of A*v_{j}.   |
c           %----------------------------------%
c
            wnorm = sdot (n, resid, 1, workd(ivj), 1)
            wnorm = sqrt(abs(wnorm))
         else if (bmat .eq. 'G') then         
            wnorm = sdot (n, resid, 1, workd(ipj), 1)
            wnorm = sqrt(abs(wnorm))
         else if (bmat .eq. 'I') then
            wnorm = snrm2(n, resid, 1)
         end if
c
c        %-----------------------------------------%
c        | Compute the j-th residual corresponding |
c        | to the j step factorization.            |
c        | Use Classical Gram Schmidt and compute: |
c        | w_{j} <-  V_{j}^T * B * OP * v_{j}      |
c        | r_{j} <-  OP*v_{j} - V_{j} * w_{j}      |
c        %-----------------------------------------%
c
c
c        %------------------------------------------%
c        | Compute the j Fourier coefficients w_{j} |
c        | WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  |
c        %------------------------------------------%
c
         if (mode .ne. 2 ) then
            call sgemv('T', n, j, one, v, ldv, workd(ipj), 1, zero, 
     &                  workd(irj), 1)
         else if (mode .eq. 2) then
            call sgemv('T', n, j, one, v, ldv, workd(ivj), 1, zero, 
     &                  workd(irj), 1)
         end if
c
c        %--------------------------------------%
c        | Orthgonalize r_{j} against V_{j}.    |
c        | RESID contains OP*v_{j}. See STEP 3. | 
c        %--------------------------------------%
c
         call sgemv('N', n, j, -one, v, ldv, workd(irj), 1, one, 
     &               resid, 1)
c
c        %--------------------------------------%
c        | Extend H to have j rows and columns. |
c        %--------------------------------------%
c
         h(j,2) = workd(irj + j - 1)
         if (j .eq. 1  .or.  rstart) then
            h(j,1) = zero
         else
            h(j,1) = rnorm
         end if
         call arscnd (t4)
c 
         orth1 = .true.
         iter  = 0
c 
         call arscnd (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call scopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c 
c           %----------------------------------%
c           | Exit in order to compute B*r_{j} |
c           %----------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call scopy (n, resid, 1, workd(ipj), 1)
         end if
   70    continue
c 
c        %---------------------------------------------------%
c        | Back from reverse communication if ORTH1 = .true. |
c        | WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    |
c        %---------------------------------------------------%
c
         if (bmat .eq. 'G') then
            call arscnd (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
c 
         orth1 = .false.
c
c        %------------------------------%
c        | Compute the B-norm of r_{j}. |
c        %------------------------------%
c
         if (bmat .eq. 'G') then         
            rnorm = sdot (n, resid, 1, workd(ipj), 1)
            rnorm = sqrt(abs(rnorm))
         else if (bmat .eq. 'I') then
            rnorm = snrm2(n, resid, 1)
         end if
c
c        %-----------------------------------------------------------%
c        | STEP 5: Re-orthogonalization / Iterative refinement phase |
c        | Maximum NITER_ITREF tries.                                |
c        |                                                           |
c        |          s      = V_{j}^T * B * r_{j}                     |
c        |          r_{j}  = r_{j} - V_{j}*s                         |
c        |          alphaj = alphaj + s_{j}                          |
c        |                                                           |
c        | The stopping criteria used for iterative refinement is    |
c        | discussed in Parlett's book SEP, page 107 and in Gragg &  |
c        | Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         |
c        | Determine if we need to correct the residual. The goal is |
c        | to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||  |
c        %-----------------------------------------------------------%
c
         if (rnorm .gt. 0.717*wnorm) go to 100
         nrorth = nrorth + 1
c 
c        %---------------------------------------------------%
c        | Enter the Iterative refinement phase. If further  |
c        | refinement is necessary, loop back here. The loop |
c        | variable is ITER. Perform a step of Classical     |
c        | Gram-Schmidt using all the Arnoldi vectors V_{j}  |
c        %---------------------------------------------------%
c
   80    continue
c
         if (msglvl .gt. 2) then
            xtemp(1) = wnorm
            xtemp(2) = rnorm
            call svout (logfil, 2, xtemp, ndigit, 
     &           '_saitr: re-orthonalization ; wnorm and rnorm are')
         end if
c
c        %----------------------------------------------------%
c        | Compute V_{j}^T * B * r_{j}.                       |
c        | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). |
c        %----------------------------------------------------%
c
         call sgemv ('T', n, j, one, v, ldv, workd(ipj), 1, 
     &               zero, workd(irj), 1)
c
c        %----------------------------------------------%
c        | Compute the correction to the residual:      |
c        | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1).  |
c        | The correction to H is v(:,1:J)*H(1:J,1:J) + |
c        | v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j, but only   |
c        | H(j,j) is updated.                           |
c        %----------------------------------------------%
c
         call sgemv ('N', n, j, -one, v, ldv, workd(irj), 1, 
     &               one, resid, 1)
c
         if (j .eq. 1  .or.  rstart) h(j,1) = zero
         h(j,2) = h(j,2) + workd(irj + j - 1)
c 
         orth2 = .true.
         call arscnd (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call scopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c 
c           %-----------------------------------%
c           | Exit in order to compute B*r_{j}. |
c           | r_{j} is the corrected residual.  |
c           %-----------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call scopy (n, resid, 1, workd(ipj), 1)
         end if
   90    continue
c
c        %---------------------------------------------------%
c        | Back from reverse communication if ORTH2 = .true. |
c        %---------------------------------------------------%
c
         if (bmat .eq. 'G') then
            call arscnd (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
c
c        %-----------------------------------------------------%
c        | Compute the B-norm of the corrected residual r_{j}. |
c        %-----------------------------------------------------%
c 
         if (bmat .eq. 'G') then         
             rnorm1 = sdot (n, resid, 1, workd(ipj), 1)
             rnorm1 = sqrt(abs(rnorm1))
         else if (bmat .eq. 'I') then
             rnorm1 = snrm2(n, resid, 1)
         end if
c
         if (msglvl .gt. 0 .and. iter .gt. 0) then
            call ivout (logfil, 1, j, ndigit,
     &           '_saitr: Iterative refinement for Arnoldi residual')
            if (msglvl .gt. 2) then
                xtemp(1) = rnorm
                xtemp(2) = rnorm1
                call svout (logfil, 2, xtemp, ndigit,
     &           '_saitr: iterative refinement ; rnorm and rnorm1 are')
            end if
         end if
c 
c        %-----------------------------------------%
c        | Determine if we need to perform another |
c        | step of re-orthogonalization.           |
c        %-----------------------------------------%
c
         if (rnorm1 .gt. 0.717*rnorm) then
c
c           %--------------------------------%
c           | No need for further refinement |
c           %--------------------------------%
c
            rnorm = rnorm1
c 
         else
c
c           %-------------------------------------------%
c           | Another step of iterative refinement step |
c           | is required. NITREF is used by stat.h     |
c           %-------------------------------------------%
c
            nitref = nitref + 1
            rnorm  = rnorm1
            iter   = iter + 1
            if (iter .le. 1) go to 80
c
c           %-------------------------------------------------%
c           | Otherwise RESID is numerically in the span of V |
c           %-------------------------------------------------%
c
            do 95 jj = 1, n
               resid(jj) = zero
  95        continue
            rnorm = zero
         end if
c 
c        %----------------------------------------------%
c        | Branch here directly if iterative refinement |
c        | wasn't necessary or after at most NITER_REF  |
c        | steps of iterative refinement.               |
c        %----------------------------------------------%
c
  100    continue
c 
         rstart = .false.
         orth2  = .false.
c 
         call arscnd (t5)
         titref = titref + (t5 - t4)
c 
c        %----------------------------------------------------------%
c        | Make sure the last off-diagonal element is non negative  |
c        | If not perform a similarity transformation on H(1:j,1:j) |
c        | and scale v(:,j) by -1.                                  |
c        %----------------------------------------------------------%
c
         if (h(j,1) .lt. zero) then
            h(j,1) = -h(j,1)
            if ( j .lt. k+np) then 
               call sscal(n, -one, v(1,j+1), 1)
            else
               call sscal(n, -one, resid, 1)
            end if
         end if
c 
c        %------------------------------------%
c        | STEP 6: Update  j = j+1;  Continue |
c        %------------------------------------%
c
         j = j + 1
         if (j .gt. k+np) then
            call arscnd (t1)
            tsaitr = tsaitr + (t1 - t0)
            ido = 99
c
            if (msglvl .gt. 1) then
               call svout (logfil, k+np, h(1,2), ndigit, 
     &         '_saitr: main diagonal of matrix H of step K+NP.')
               if (k+np .gt. 1) then
               call svout (logfil, k+np-1, h(2,1), ndigit, 
     &         '_saitr: sub diagonal of matrix H of step K+NP.')
               end if
            end if
c
            go to 9000
         end if
c
c        %--------------------------------------------------------%
c        | Loop back to extend the factorization by another step. |
c        %--------------------------------------------------------%
c
      go to 1000
c 
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
 9000 continue
      return
c
c     %---------------%
c     | End of ssaitr |
c     %---------------%
c
      end
