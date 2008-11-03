c\BeginDoc
c
c\Name: cnaup2
c
c\Description: 
c  Intermediate level interface called by cnaupd.
c
c\Usage:
c  call cnaup2
c     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
c       ISHIFT, MXITER, V, LDV, H, LDH, RITZ, BOUNDS, 
c       Q, LDQ, WORKL, IPNTR, WORKD, RWORK, INFO )
c
c\Arguments
c
c  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in cnaupd.
c  MODE, ISHIFT, MXITER: see the definition of IPARAM in cnaupd.
c
c  NP      Integer.  (INPUT/OUTPUT)
c          Contains the number of implicit shifts to apply during
c          each Arnoldi iteration.
c          If ISHIFT=1, NP is adjusted dynamically at each iteration
c          to accelerate convergence and prevent stagnation.
c          This is also roughly equal to the number of matrix-vector
c          products (involving the operator OP) per Arnoldi iteration.
c          The logic for adjusting is contained within the current
c          subroutine.
c          If ISHIFT=0, NP is the number of shifts the user needs
c          to provide via reverse comunication. 0 < NP < NCV-NEV.
c          NP may be less than NCV-NEV since a leading block of the current
c          upper Hessenberg matrix has split off and contains "unwanted"
c          Ritz values.
c          Upon termination of the IRA iteration, NP contains the number
c          of "converged" wanted Ritz values.
c
c  IUPD    Integer.  (INPUT)
c          IUPD .EQ. 0: use explicit restart instead implicit update.
c          IUPD .NE. 0: use implicit update.
c
c  V       Complex  N by (NEV+NP) array.  (INPUT/OUTPUT)
c          The Arnoldi basis vectors are returned in the first NEV 
c          columns of V.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling 
c          program.
c
c  H       Complex  (NEV+NP) by (NEV+NP) array.  (OUTPUT)
c          H is used to store the generated upper Hessenberg matrix
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling 
c          program.
c
c  RITZ    Complex  array of length NEV+NP.  (OUTPUT)
c          RITZ(1:NEV)  contains the computed Ritz values of OP.
c
c  BOUNDS  Complex  array of length NEV+NP.  (OUTPUT)
c          BOUNDS(1:NEV) contain the error bounds corresponding to 
c          the computed Ritz values.
c          
c  Q       Complex  (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
c          Private (replicated) work array used to accumulate the
c          rotation in the shift application step.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Complex  work array of length at least 
c          (NEV+NP)**2 + 3*(NEV+NP).  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  It is used in shifts calculation, shifts
c          application and convergence checking.
c
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD for 
c          vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X.
c          IPNTR(2): pointer to the current result vector Y.
c          IPNTR(3): pointer to the vector B * X when used in the 
c                    shift-and-invert mode.  X is the current operand.
c          -------------------------------------------------------------
c          
c  WORKD   Complex  work array of length 3*N.  (WORKSPACE)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD
c          as temporary workspace during the iteration !!!!!!!!!!
c          See Data Distribution Note in CNAUPD.
c
c  RWORK   Real    work array of length  NEV+NP ( WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =     0: Normal return.
c          =     1: Maximum number of iterations taken.
c                   All possible eigenvalues of OP has been found.  
c                   NP returns the number of converged Ritz values.
c          =     2: No shifts could be applied.
c          =    -8: Error return from LAPACK eigenvalue calculation;
c                   This should never happen.
c          =    -9: Starting vector is zero.
c          = -9999: Could not build an Arnoldi factorization.
c                   Size that was built in returned in NP.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  Complex 
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c
c\Routines called:
c     cgetv0  ARPACK initial vector generation routine. 
c     cnaitr  ARPACK Arnoldi factorization routine.
c     cnapps  ARPACK application of implicit shifts routine.
c     cneigh  ARPACK compute Ritz values and error bounds routine. 
c     cngets  ARPACK reorder Ritz values and error bounds routine.
c     csortc  ARPACK sorting routine.
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     cmout   ARPACK utility routine that prints matrices
c     cvout   ARPACK utility routine that prints vectors.
c     svout   ARPACK utility routine that prints vectors.
c     slamch  LAPACK routine that determines machine constants.
c     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     ccopy   Level 1 BLAS that copies one vector to another .
c     wcdotc   Level 1 BLAS that computes the scalar product of two vectors. 
c     cswap   Level 1 BLAS that swaps two vectors.
c     scnrm2  Level 1 BLAS that computes the norm of a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice Universitya
c     Chao Yang                    Houston, Texas
c     Dept. of Computational &
c     Applied Mathematics 
c     Rice University           
c     Houston, Texas 
c 
c\SCCS Information: @(#)
c FILE: naup2.F   SID: 2.6   DATE OF SID: 06/01/00   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine cnaup2
     &   ( ido, bmat, n, which, nev, np, tol, resid, mode, iupd, 
     &     ishift, mxiter, v, ldv, h, ldh, ritz, bounds, 
     &     q, ldq, workl, ipntr, workd, rwork, info )
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
      character  bmat*1, which*2
      integer    ido, info, ishift, iupd, mode, ldh, ldq, ldv, mxiter,
     &           n, nev, np
      Real   
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(13)
      Complex 
     &           bounds(nev+np), h(ldh,nev+np), q(ldq,nev+np), 
     &           resid(n), ritz(nev+np),  v(ldv,nev+np), 
     &           workd(3*n), workl( (nev+np)*(nev+np+3) )
       Real   
     &           rwork(nev+np)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex 
     &           one, zero
      Real 
     &           rzero
      parameter (one = (1.0E+0, 0.0E+0) , zero = (0.0E+0, 0.0E+0) ,
     &           rzero = 0.0E+0 )
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      logical    cnorm , getv0, initv , update, ushift
      integer    ierr  , iter , kplusp, msglvl, nconv, 
     &           nevbef, nev0 , np0   , nptemp, i    ,
     &           j    
      Complex 
     &           cmpnorm
      Real 
     &           rnorm , eps23, rtemp
      character  wprime*2
c
      save       cnorm,  getv0, initv , update, ushift, 
     &           rnorm,  iter , kplusp, msglvl, nconv ,
     &           nevbef, nev0 , np0   , eps23
c
c
c     %-----------------------%
c     | Local array arguments |
c     %-----------------------%
c
      integer    kp(3)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   ccopy, cgetv0, cnaitr, cneigh, cngets, cnapps,
     &           csortc, cswap, cmout, cvout, ivout, second
c
c     %--------------------%
c     | External functions |
c     %--------------------%
c
      Complex 
     &           wcdotc
      Real   
     &           scnrm2, slamch, slapy2
      external   wcdotc, scnrm2, slamch, slapy2
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic  aimag, real , min, max
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (ido .eq. 0) then
c 
         call second (t0)
c 
         msglvl = mcaup2
c 
         nev0   = nev
         np0    = np
c
c        %-------------------------------------%
c        | kplusp is the bound on the largest  |
c        |        Lanczos factorization built. |
c        | nconv is the current number of      |
c        |        "converged" eigenvalues.     |
c        | iter is the counter on the current  |
c        |      iteration step.                |
c        %-------------------------------------%
c
         kplusp = nev + np
         nconv  = 0
         iter   = 0
c 
c        %---------------------------------%
c        | Get machine dependent constant. |
c        %---------------------------------%
c
         eps23 = slamch('Epsilon-Machine')
         eps23 = eps23**(2.0E+0  / 3.0E+0 )
c
c        %---------------------------------------%
c        | Set flags for computing the first NEV |
c        | steps of the Arnoldi factorization.   |
c        %---------------------------------------%
c
         getv0    = .true.
         update   = .false.
         ushift   = .false.
         cnorm    = .false.
c
         if (info .ne. 0) then
c
c           %--------------------------------------------%
c           | User provides the initial residual vector. |
c           %--------------------------------------------%
c
            initv = .true.
            info  = 0
         else
            initv = .false.
         end if
      end if
c 
c     %---------------------------------------------%
c     | Get a possibly random starting vector and   |
c     | force it into the range of the operator OP. |
c     %---------------------------------------------%
c
   10 continue
c
      if (getv0) then
         call cgetv0 (ido, bmat, 1, initv, n, 1, v, ldv, resid, rnorm,
     &                ipntr, workd, info)
c
         if (ido .ne. 99) go to 9000
c
         if (rnorm .eq. rzero) then
c
c           %-----------------------------------------%
c           | The initial vector is zero. Error exit. | 
c           %-----------------------------------------%
c
            info = -9
            go to 1100
         end if
         getv0 = .false.
         ido  = 0
      end if
c 
c     %-----------------------------------%
c     | Back from reverse communication : |
c     | continue with update step         |
c     %-----------------------------------%
c
      if (update) go to 20
c
c     %-------------------------------------------%
c     | Back from computing user specified shifts |
c     %-------------------------------------------%
c
      if (ushift) go to 50
c
c     %-------------------------------------%
c     | Back from computing residual norm   |
c     | at the end of the current iteration |
c     %-------------------------------------%
c
      if (cnorm)  go to 100
c 
c     %----------------------------------------------------------%
c     | Compute the first NEV steps of the Arnoldi factorization |
c     %----------------------------------------------------------%
c
      call cnaitr (ido, bmat, n, 0, nev, mode, resid, rnorm, v, ldv, 
     &             h, ldh, ipntr, workd, info)
c
      if (ido .ne. 99) go to 9000
c
      if (info .gt. 0) then
         np   = info
         mxiter = iter
         info = -9999
         go to 1200
      end if
c 
c     %--------------------------------------------------------------%
c     |                                                              |
c     |           M A I N  ARNOLDI  I T E R A T I O N  L O O P       |
c     |           Each iteration implicitly restarts the Arnoldi     |
c     |           factorization in place.                            |
c     |                                                              |
c     %--------------------------------------------------------------%
c 
 1000 continue
c
         iter = iter + 1
c
         if (msglvl .gt. 0) then
            call ivout (logfil, 1, iter, ndigit, 
     &           '_naup2: **** Start of major iteration number ****')
         end if
c 
c        %-----------------------------------------------------------%
c        | Compute NP additional steps of the Arnoldi factorization. |
c        | Adjust NP since NEV might have been updated by last call  |
c        | to the shift application routine cnapps.                  |
c        %-----------------------------------------------------------%
c
         np  = kplusp - nev
c
         if (msglvl .gt. 1) then
            call ivout (logfil, 1, nev, ndigit, 
     &     '_naup2: The length of the current Arnoldi factorization')
            call ivout (logfil, 1, np, ndigit, 
     &           '_naup2: Extend the Arnoldi factorization by')
         end if
c
c        %-----------------------------------------------------------%
c        | Compute NP additional steps of the Arnoldi factorization. |
c        %-----------------------------------------------------------%
c
         ido = 0
   20    continue
         update = .true.
c
         call cnaitr(ido, bmat, n, nev, np,    mode,  resid, rnorm,
     &               v  , ldv , h, ldh, ipntr, workd, info)
c
         if (ido .ne. 99) go to 9000
c
         if (info .gt. 0) then
            np = info
            mxiter = iter
            info = -9999
            go to 1200
         end if
         update = .false.
c
         if (msglvl .gt. 1) then
            call svout (logfil, 1, rnorm, ndigit, 
     &           '_naup2: Corresponding B-norm of the residual')
         end if
c 
c        %--------------------------------------------------------%
c        | Compute the eigenvalues and corresponding error bounds |
c        | of the current upper Hessenberg matrix.                |
c        %--------------------------------------------------------%
c
         call cneigh (rnorm, kplusp, h, ldh, ritz, bounds,
     &                q, ldq, workl, rwork,  ierr)
c
         if (ierr .ne. 0) then
            info = -8
            go to 1200
         end if
c
c        %---------------------------------------------------%
c        | Select the wanted Ritz values and their bounds    |
c        | to be used in the convergence test.               |
c        | The wanted part of the spectrum and corresponding |
c        | error bounds are in the last NEV loc. of RITZ,    |
c        | and BOUNDS respectively.                          | 
c        %---------------------------------------------------%
c
         nev = nev0
         np = np0
c
c        %--------------------------------------------------%
c        | Make a copy of Ritz values and the corresponding |
c        | Ritz estimates obtained from cneigh.             |
c        %--------------------------------------------------%
c
         call ccopy(kplusp,ritz,1,workl(kplusp**2+1),1)
         call ccopy(kplusp,bounds,1,workl(kplusp**2+kplusp+1),1)
c
c        %---------------------------------------------------%
c        | Select the wanted Ritz values and their bounds    |
c        | to be used in the convergence test.               |
c        | The wanted part of the spectrum and corresponding |
c        | bounds are in the last NEV loc. of RITZ           |
c        | BOUNDS respectively.                              |
c        %---------------------------------------------------%
c
         call cngets (ishift, which, nev, np, ritz, bounds)
c 
c        %------------------------------------------------------------%
c        | Convergence test: currently we use the following criteria. |
c        | The relative accuracy of a Ritz value is considered        |
c        | acceptable if:                                             |
c        |                                                            |
c        | error_bounds(i) .le. tol*max(eps23, magnitude_of_ritz(i)). |
c        |                                                            |
c        %------------------------------------------------------------%
c
         nconv  = 0
c
         do 25 i = 1, nev
            rtemp = max( eps23, slapy2( real (ritz(np+i)),
     &                                  aimag(ritz(np+i)) ) ) 
            if ( slapy2(real (bounds(np+i)),aimag(bounds(np+i))) 
     &                 .le. tol*rtemp ) then
               nconv = nconv + 1
            end if
   25    continue
c 
         if (msglvl .gt. 2) then
            kp(1) = nev
            kp(2) = np
            kp(3) = nconv
            call ivout (logfil, 3, kp, ndigit, 
     &                  '_naup2: NEV, NP, NCONV are')
            call cvout (logfil, kplusp, ritz, ndigit,
     &           '_naup2: The eigenvalues of H')
            call cvout (logfil, kplusp, bounds, ndigit, 
     &          '_naup2: Ritz estimates of the current NCV Ritz values')
         end if
c
c        %---------------------------------------------------------%
c        | Count the number of unwanted Ritz values that have zero |
c        | Ritz estimates. If any Ritz estimates are equal to zero |
c        | then a leading block of H of order equal to at least    |
c        | the number of Ritz values with zero Ritz estimates has  |
c        | split off. None of these Ritz values may be removed by  |
c        | shifting. Decrease NP the number of shifts to apply. If |
c        | no shifts may be applied, then prepare to exit          |
c        %---------------------------------------------------------%
c
         nptemp = np
         do 30 j=1, nptemp
            if (bounds(j) .eq. zero) then
               np = np - 1
               nev = nev + 1
            end if
 30      continue
c     
         if ( (nconv .ge. nev0) .or. 
     &        (iter .gt. mxiter) .or.
     &        (np .eq. 0) ) then
c
            if (msglvl .gt. 4) then
               call cvout(logfil, kplusp, workl(kplusp**2+1), ndigit,
     &             '_naup2: Eigenvalues computed by _neigh:')
               call cvout(logfil, kplusp, workl(kplusp**2+kplusp+1),
     &                     ndigit,
     &             '_naup2: Ritz estimates computed by _neigh:')
            end if
c     
c           %------------------------------------------------%
c           | Prepare to exit. Put the converged Ritz values |
c           | and corresponding bounds in RITZ(1:NCONV) and  |
c           | BOUNDS(1:NCONV) respectively. Then sort. Be    |
c           | careful when NCONV > NP                        |
c           %------------------------------------------------%
c
c           %------------------------------------------%
c           |  Use h( 3,1 ) as storage to communicate  |
c           |  rnorm to cneupd if needed               |
c           %------------------------------------------%

            h(3,1) = cmplx(rnorm,rzero)
c
c           %----------------------------------------------%
c           | Sort Ritz values so that converged Ritz      |
c           | values appear within the first NEV locations |
c           | of ritz and bounds, and the most desired one |
c           | appears at the front.                        |
c           %----------------------------------------------%
c
            if (which .eq. 'LM') wprime = 'SM'
            if (which .eq. 'SM') wprime = 'LM'
            if (which .eq. 'LR') wprime = 'SR'
            if (which .eq. 'SR') wprime = 'LR'
            if (which .eq. 'LI') wprime = 'SI'
            if (which .eq. 'SI') wprime = 'LI'
c
            call csortc(wprime, .true., kplusp, ritz, bounds)
c
c           %--------------------------------------------------%
c           | Scale the Ritz estimate of each Ritz value       |
c           | by 1 / max(eps23, magnitude of the Ritz value).  |
c           %--------------------------------------------------%
c
            do 35 j = 1, nev0 
                rtemp = max( eps23, slapy2( real (ritz(j)),
     &                                       aimag(ritz(j)) ) )
                bounds(j) = bounds(j)/rtemp
 35         continue
c
c           %---------------------------------------------------%
c           | Sort the Ritz values according to the scaled Ritz |
c           | estimates.  This will push all the converged ones |
c           | towards the front of ritz, bounds (in the case    |
c           | when NCONV < NEV.)                                |
c           %---------------------------------------------------%
c
            wprime = 'LM'
            call csortc(wprime, .true., nev0, bounds, ritz)
c
c           %----------------------------------------------%
c           | Scale the Ritz estimate back to its original |
c           | value.                                       |
c           %----------------------------------------------%
c
            do 40 j = 1, nev0
                rtemp = max( eps23, slapy2( real (ritz(j)),
     &                                       aimag(ritz(j)) ) )
                bounds(j) = bounds(j)*rtemp
 40         continue
c
c           %-----------------------------------------------%
c           | Sort the converged Ritz values again so that  |
c           | the "threshold" value appears at the front of |
c           | ritz and bound.                               |
c           %-----------------------------------------------%
c
            call csortc(which, .true., nconv, ritz, bounds)
c
            if (msglvl .gt. 1) then
               call cvout (logfil, kplusp, ritz, ndigit,
     &            '_naup2: Sorted eigenvalues')
               call cvout (logfil, kplusp, bounds, ndigit,
     &            '_naup2: Sorted ritz estimates.')
            end if
c
c           %------------------------------------%
c           | Max iterations have been exceeded. | 
c           %------------------------------------%
c
            if (iter .gt. mxiter .and. nconv .lt. nev0) info = 1
c
c           %---------------------%
c           | No shifts to apply. | 
c           %---------------------%
c
            if (np .eq. 0 .and. nconv .lt. nev0)  info = 2
c
            np = nconv
            go to 1100
c
         else if ( (nconv .lt. nev0) .and. (ishift .eq. 1) ) then
c     
c           %-------------------------------------------------%
c           | Do not have all the requested eigenvalues yet.  |
c           | To prevent possible stagnation, adjust the size |
c           | of NEV.                                         |
c           %-------------------------------------------------%
c
            nevbef = nev
            nev = nev + min(nconv, np/2)
            if (nev .eq. 1 .and. kplusp .ge. 6) then
               nev = kplusp / 2
            else if (nev .eq. 1 .and. kplusp .gt. 3) then
               nev = 2
            end if
            np = kplusp - nev
c     
c           %---------------------------------------%
c           | If the size of NEV was just increased |
c           | resort the eigenvalues.               |
c           %---------------------------------------%
c     
            if (nevbef .lt. nev) 
     &         call cngets (ishift, which, nev, np, ritz, bounds)
c
         end if              
c     
         if (msglvl .gt. 0) then
            call ivout (logfil, 1, nconv, ndigit, 
     &           '_naup2: no. of "converged" Ritz values at this iter.')
            if (msglvl .gt. 1) then
               kp(1) = nev
               kp(2) = np
               call ivout (logfil, 2, kp, ndigit, 
     &              '_naup2: NEV and NP are')
               call cvout (logfil, nev, ritz(np+1), ndigit,
     &              '_naup2: "wanted" Ritz values ')
               call cvout (logfil, nev, bounds(np+1), ndigit,
     &              '_naup2: Ritz estimates of the "wanted" values ')
            end if
         end if
c
         if (ishift .eq. 0) then
c
c           %-------------------------------------------------------%
c           | User specified shifts: pop back out to get the shifts |
c           | and return them in the first 2*NP locations of WORKL. |
c           %-------------------------------------------------------%
c
            ushift = .true.
            ido = 3
            go to 9000
         end if
   50    continue
         ushift = .false.
c
         if ( ishift .ne. 1 ) then
c 
c            %----------------------------------%
c            | Move the NP shifts from WORKL to |
c            | RITZ, to free up WORKL           |
c            | for non-exact shift case.        |
c            %----------------------------------%
c
             call ccopy (np, workl, 1, ritz, 1)
         end if
c
         if (msglvl .gt. 2) then 
            call ivout (logfil, 1, np, ndigit, 
     &                  '_naup2: The number of shifts to apply ')
            call cvout (logfil, np, ritz, ndigit,
     &                  '_naup2: values of the shifts')
            if ( ishift .eq. 1 ) 
     &          call cvout (logfil, np, bounds, ndigit,
     &                  '_naup2: Ritz estimates of the shifts')
         end if
c
c        %---------------------------------------------------------%
c        | Apply the NP implicit shifts by QR bulge chasing.       |
c        | Each shift is applied to the whole upper Hessenberg     |
c        | matrix H.                                               |
c        | The first 2*N locations of WORKD are used as workspace. |
c        %---------------------------------------------------------%
c
         call cnapps (n, nev, np, ritz, v, ldv, 
     &                h, ldh, resid, q, ldq, workl, workd)
c
c        %---------------------------------------------%
c        | Compute the B-norm of the updated residual. |
c        | Keep B*RESID in WORKD(1:N) to be used in    |
c        | the first step of the next call to cnaitr.  |
c        %---------------------------------------------%
c
         cnorm = .true.
         call second (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call ccopy (n, resid, 1, workd(n+1), 1)
            ipntr(1) = n + 1
            ipntr(2) = 1
            ido = 2
c 
c           %----------------------------------%
c           | Exit in order to compute B*RESID |
c           %----------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call ccopy (n, resid, 1, workd, 1)
         end if
c 
  100    continue
c 
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(1:N) := B*RESID            |
c        %----------------------------------%
c
         if (bmat .eq. 'G') then
            call second (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
c 
         if (bmat .eq. 'G') then         
            cmpnorm = wcdotc (n, resid, 1, workd, 1)
            rnorm = sqrt(slapy2(real (cmpnorm),aimag(cmpnorm)))
         else if (bmat .eq. 'I') then
            rnorm = scnrm2(n, resid, 1)
         end if
         cnorm = .false.
c
         if (msglvl .gt. 2) then
            call svout (logfil, 1, rnorm, ndigit, 
     &      '_naup2: B-norm of residual for compressed factorization')
            call cmout (logfil, nev, nev, h, ldh, ndigit,
     &        '_naup2: Compressed upper Hessenberg matrix H')
         end if
c 
      go to 1000
c
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
 1100 continue
c
      mxiter = iter
      nev = nconv
c     
 1200 continue
      ido = 99
c
c     %------------%
c     | Error Exit |
c     %------------%
c
      call second (t1)
      tcaup2 = t1 - t0
c     
 9000 continue
c
c     %---------------%
c     | End of cnaup2 |
c     %---------------%
c
      return
      end
