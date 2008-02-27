c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsaup2
c
c\Description: 
c  Intermediate level interface called by dsaupd.
c
c\Usage:
c  call dsaup2 
c     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
c       ISHIFT, MXITER, V, LDV, H, LDH, RITZ, BOUNDS, Q, LDQ, WORKL, 
c       IPNTR, WORKD, INFO )
c
c\Arguments
c
c  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in dsaupd.
c  MODE, ISHIFT, MXITER: see the definition of IPARAM in dsaupd.
c  
c  NP      Integer.  (INPUT/OUTPUT)
c          Contains the number of implicit shifts to apply during 
c          each Arnoldi/Lanczos iteration.  
c          If ISHIFT=1, NP is adjusted dynamically at each iteration 
c          to accelerate convergence and prevent stagnation.
c          This is also roughly equal to the number of matrix-vector 
c          products (involving the operator OP) per Arnoldi iteration.
c          The logic for adjusting is contained within the current
c          subroutine.
c          If ISHIFT=0, NP is the number of shifts the user needs
c          to provide via reverse comunication. 0 < NP < NCV-NEV.
c          NP may be less than NCV-NEV since a leading block of the current
c          upper Tridiagonal matrix has split off and contains "unwanted"
c          Ritz values.
c          Upon termination of the IRA iteration, NP contains the number 
c          of "converged" wanted Ritz values.
c
c  IUPD    Integer.  (INPUT)
c          IUPD .EQ. 0: use explicit restart instead implicit update.
c          IUPD .NE. 0: use implicit update.
c
c  V       Double precision N by (NEV+NP) array.  (INPUT/OUTPUT)
c          The Lanczos basis vectors.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling 
c          program.
c
c  H       Double precision (NEV+NP) by 2 array.  (OUTPUT)
c          H is used to store the generated symmetric tridiagonal matrix
c          The subdiagonal is stored in the first column of H starting 
c          at H(2,1).  The main diagonal is stored in the second column
c          of H starting at H(1,2). If dsaup2 converges store the 
c          B-norm of the final residual vector in H(1,1).
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling 
c          program.
c
c  RITZ    Double precision array of length NEV+NP.  (OUTPUT)
c          RITZ(1:NEV) contains the computed Ritz values of OP.
c
c  BOUNDS  Double precision array of length NEV+NP.  (OUTPUT)
c          BOUNDS(1:NEV) contain the error bounds corresponding to RITZ.
c
c  Q       Double precision (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
c          Private (replicated) work array used to accumulate the 
c          rotation in the shift application step.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c          
c  WORKL   Double precision array of length at least 3*(NEV+NP).  (INPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  It is used in the computation of the 
c          tridiagonal eigenvalue problem, the calculation and
c          application of the shifts and convergence checking.
c          If ISHIFT .EQ. O and IDO .EQ. 3, the first NP locations
c          of WORKL are used in reverse communication to hold the user 
c          supplied shifts.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD for 
c          vectors used by the Lanczos iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X.
c          IPNTR(2): pointer to the current result vector Y.
c          IPNTR(3): pointer to the vector B * X when used in one of  
c                    the spectral transformation modes.  X is the current
c                    operand.
c          -------------------------------------------------------------
c          
c  WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Lanczos iteration
c          for reverse communication.  The user should not use WORKD
c          as temporary workspace during the iteration !!!!!!!!!!
c          See Data Distribution Note in dsaupd.
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =     0: Normal return.
c          =     1: All possible eigenvalues of OP has been found.  
c                   NP returns the size of the invariant subspace
c                   spanning the operator OP. 
c          =     2: No shifts could be applied.
c          =    -8: Error return from trid. eigenvalue calculation;
c                   This should never happen.
c          =    -9: Starting vector is zero.
c          = -9999: Could not build an Lanczos factorization.
c                   Size that was built in returned in NP.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B.N. Parlett, "The Symmetric Eigenvalue Problem". Prentice-Hall,
c     1980.
c  4. B.N. Parlett, B. Nour-Omid, "Towards a Black Box Lanczos Program",
c     Computer Physics Communications, 53 (1989), pp 169-179.
c  5. B. Nour-Omid, B.N. Parlett, T. Ericson, P.S. Jensen, "How to
c     Implement the Spectral Transformation", Math. Comp., 48 (1987),
c     pp 663-673.
c  6. R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos 
c     Algorithm for Solving Sparse Symmetric Generalized Eigenproblems", 
c     SIAM J. Matr. Anal. Apps.,  January (1993).
c  7. L. Reichel, W.B. Gragg, "Algorithm 686: FORTRAN Subroutines
c     for Updating the QR decomposition", ACM TOMS, December 1990,
c     Volume 16 Number 4, pp 369-377.
c
c\Routines called:
c     dgetv0  ARPACK initial vector generation routine. 
c     dsaitr  ARPACK Lanczos factorization routine.
c     dsapps  ARPACK application of implicit shifts routine.
c     dsconv  ARPACK convergence of Ritz values routine.
c     dseigt  ARPACK compute Ritz values and error bounds routine.
c     dsgets  ARPACK reorder Ritz values and error bounds routine.
c     dsortr  ARPACK sorting routine.
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     dvout   ARPACK utility routine that prints vectors.
c     dlamch  LAPACK routine that determines machine constants.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     ddot    Level 1 BLAS that computes the scalar product of two vectors. 
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     dscal   Level 1 BLAS that scales a vector.
c     dswap   Level 1 BLAS that swaps two vectors.
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
c     12/15/93: Version ' 2.4'
c     xx/xx/95: Version ' 2.4'.  (R.B. Lehoucq)
c
c\SCCS Information: @(#) 
c FILE: saup2.F   SID: 2.7   DATE OF SID: 5/19/98   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsaup2
     &   ( ido, bmat, n, which, nev, np, tol, resid, mode, iupd, 
     &     ishift, mxiter, v, ldv, h, ldh, ritz, bounds, 
     &     q, ldq, workl, ipntr, workd, info )
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
      integer    ido, info, ishift, iupd, ldh, ldq, ldv, mxiter,
     &           n, mode, nev, np
      Double precision
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(3)
      Double precision
     &           bounds(nev+np), h(ldh,2), q(ldq,nev+np), resid(n), 
     &           ritz(nev+np), v(ldv,nev+np), workd(3*n), 
     &           workl(3*(nev+np))
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character  wprime*2
      logical    cnorm, getv0, initv, update, ushift
      integer    ierr, iter, j, kplusp, msglvl, nconv, nevbef, nev0, 
     &           np0, nptemp, nevd2, nevm2, kp(3) 
      Double precision
     &           rnorm, temp, eps23
      save       cnorm, getv0, initv, update, ushift,
     &           iter, kplusp, msglvl, nconv, nev0, np0,
     &           rnorm, eps23
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, dgetv0, dsaitr, dscal, dsconv, dseigt, dsgets, 
     &           dsapps, dsortr, dvout, ivout, second, dswap
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           ddot, dnrm2, dlamch
      external   ddot, dnrm2, dlamch
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    min
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (ido .eq. 0) then
c 
c        %-------------------------------%
c        | Initialize timing statistics  |
c        | & message level for debugging |
c        %-------------------------------%
c
         call second (t0)
         msglvl = msaup2
c
c        %---------------------------------%
c        | Set machine dependent constant. |
c        %---------------------------------%
c
         eps23 = dlamch('Epsilon-Machine')
         eps23 = eps23**(2.0D+0/3.0D+0)
c
c        %-------------------------------------%
c        | nev0 and np0 are integer variables  |
c        | hold the initial values of NEV & NP |
c        %-------------------------------------%
c
         nev0   = nev
         np0    = np
c
c        %-------------------------------------%
c        | kplusp is the bound on the largest  |
c        |        Lanczos factorization built. |
c        | nconv is the current number of      |
c        |        "converged" eigenvlues.      |
c        | iter is the counter on the current  |
c        |      iteration step.                |
c        %-------------------------------------%
c
         kplusp = nev0 + np0
         nconv  = 0
         iter   = 0
c 
c        %--------------------------------------------%
c        | Set flags for computing the first NEV steps |
c        | of the Lanczos factorization.              |
c        %--------------------------------------------%
c
         getv0    = .true.
         update   = .false.
         ushift   = .false.
         cnorm    = .false.
c
         if (info .ne. 0) then
c
c        %--------------------------------------------%
c        | User provides the initial residual vector. |
c        %--------------------------------------------%
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
         call dgetv0 (ido, bmat, 1, initv, n, 1, v, ldv, resid, rnorm,
     &                ipntr, workd, info)
c
         if (ido .ne. 99) go to 9000
c
         if (rnorm .eq. zero) then
c
c           %-----------------------------------------%
c           | The initial vector is zero. Error exit. | 
c           %-----------------------------------------%
c
            info = -9
            go to 1200
         end if
         getv0 = .false.
         ido  = 0
      end if
c 
c     %------------------------------------------------------------%
c     | Back from reverse communication: continue with update step |
c     %------------------------------------------------------------%
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
c     | Compute the first NEV steps of the Lanczos factorization |
c     %----------------------------------------------------------%
c
      call dsaitr (ido, bmat, n, 0, nev0, mode, resid, rnorm, v, ldv, 
     &             h, ldh, ipntr, workd, info)
c 
c     %---------------------------------------------------%
c     | ido .ne. 99 implies use of reverse communication  |
c     | to compute operations involving OP and possibly B |
c     %---------------------------------------------------%
c
      if (ido .ne. 99) go to 9000
c
      if (info .gt. 0) then
c
c        %-----------------------------------------------------%
c        | dsaitr was unable to build an Lanczos factorization |
c        | of length NEV0. INFO is returned with the size of   |
c        | the factorization built. Exit main loop.            |
c        %-----------------------------------------------------%
c
         np   = info
         mxiter = iter
         info = -9999
         go to 1200
      end if
c 
c     %--------------------------------------------------------------%
c     |                                                              |
c     |           M A I N  LANCZOS  I T E R A T I O N  L O O P       |
c     |           Each iteration implicitly restarts the Lanczos     |
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
     &           '_saup2: **** Start of major iteration number ****')
         end if
         if (msglvl .gt. 1) then
            call ivout (logfil, 1, nev, ndigit, 
     &     '_saup2: The length of the current Lanczos factorization')
            call ivout (logfil, 1, np, ndigit, 
     &           '_saup2: Extend the Lanczos factorization by')
         end if
c 
c        %------------------------------------------------------------%
c        | Compute NP additional steps of the Lanczos factorization. |
c        %------------------------------------------------------------%
c
         ido = 0
   20    continue
         update = .true.
c
         call dsaitr (ido, bmat, n, nev, np, mode, resid, rnorm, v, 
     &                ldv, h, ldh, ipntr, workd, info)
c 
c        %---------------------------------------------------%
c        | ido .ne. 99 implies use of reverse communication  |
c        | to compute operations involving OP and possibly B |
c        %---------------------------------------------------%
c
         if (ido .ne. 99) go to 9000
c
         if (info .gt. 0) then
c
c           %-----------------------------------------------------%
c           | dsaitr was unable to build an Lanczos factorization |
c           | of length NEV0+NP0. INFO is returned with the size  |  
c           | of the factorization built. Exit main loop.         |
c           %-----------------------------------------------------%
c
            np = info
            mxiter = iter
            info = -9999
            go to 1200
         end if
         update = .false.
c
         if (msglvl .gt. 1) then
            call dvout (logfil, 1, rnorm, ndigit, 
     &           '_saup2: Current B-norm of residual for factorization')
         end if
c 
c        %--------------------------------------------------------%
c        | Compute the eigenvalues and corresponding error bounds |
c        | of the current symmetric tridiagonal matrix.           |
c        %--------------------------------------------------------%
c
         call dseigt (rnorm, kplusp, h, ldh, ritz, bounds, workl, ierr)
c
         if (ierr .ne. 0) then
            info = -8
            go to 1200
         end if
c
c        %----------------------------------------------------%
c        | Make a copy of eigenvalues and corresponding error |
c        | bounds obtained from _seigt.                       |
c        %----------------------------------------------------%
c
         call dcopy(kplusp, ritz, 1, workl(kplusp+1), 1)
         call dcopy(kplusp, bounds, 1, workl(2*kplusp+1), 1)
c
c        %---------------------------------------------------%
c        | Select the wanted Ritz values and their bounds    |
c        | to be used in the convergence test.               |
c        | The selection is based on the requested number of |
c        | eigenvalues instead of the current NEV and NP to  |
c        | prevent possible misconvergence.                  |
c        | * Wanted Ritz values := RITZ(NP+1:NEV+NP)         |
c        | * Shifts := RITZ(1:NP) := WORKL(1:NP)             |
c        %---------------------------------------------------%
c
         nev = nev0
         np = np0
         call dsgets (ishift, which, nev, np, ritz, bounds, workl)
c 
c        %-------------------%
c        | Convergence test. |
c        %-------------------%
c
         call dcopy (nev, bounds(np+1), 1, workl(np+1), 1)
         call dsconv (nev, ritz(np+1), workl(np+1), tol, nconv)
c
         if (msglvl .gt. 2) then
            kp(1) = nev
            kp(2) = np
            kp(3) = nconv
            call ivout (logfil, 3, kp, ndigit,
     &                  '_saup2: NEV, NP, NCONV are')
            call dvout (logfil, kplusp, ritz, ndigit,
     &           '_saup2: The eigenvalues of H')
            call dvout (logfil, kplusp, bounds, ndigit,
     &          '_saup2: Ritz estimates of the current NCV Ritz values')
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
c           %------------------------------------------------%
c           | Prepare to exit. Put the converged Ritz values |
c           | and corresponding bounds in RITZ(1:NCONV) and  |
c           | BOUNDS(1:NCONV) respectively. Then sort. Be    |
c           | careful when NCONV > NP since we don't want to |
c           | swap overlapping locations.                    |
c           %------------------------------------------------%
c
            if (which .eq. 'BE') then
c
c              %-----------------------------------------------------%
c              | Both ends of the spectrum are requested.            |
c              | Sort the eigenvalues into algebraically decreasing  |
c              | order first then swap low end of the spectrum next  |
c              | to high end in appropriate locations.               |
c              | NOTE: when np < floor(nev/2) be careful not to swap |
c              | overlapping locations.                              |
c              %-----------------------------------------------------%
c
               wprime = 'SA'
               call dsortr (wprime, .true., kplusp, ritz, bounds)
               nevd2 = nev0 / 2
               nevm2 = nev0 - nevd2 
               if ( nev .gt. 1 ) then
                  call dswap ( min(nevd2,np), ritz(nevm2+1), 1,
     &                 ritz( max(kplusp-nevd2+1,kplusp-np+1) ), 1)
                  call dswap ( min(nevd2,np), bounds(nevm2+1), 1,
     &                 bounds( max(kplusp-nevd2+1,kplusp-np+1)), 1)
               end if
c
            else
c
c              %--------------------------------------------------%
c              | LM, SM, LA, SA case.                             |
c              | Sort the eigenvalues of H into the an order that |
c              | is opposite to WHICH, and apply the resulting    |
c              | order to BOUNDS.  The eigenvalues are sorted so  |
c              | that the wanted part are always within the first |
c              | NEV locations.                                   |
c              %--------------------------------------------------%
c
               if (which .eq. 'LM') wprime = 'SM'
               if (which .eq. 'SM') wprime = 'LM'
               if (which .eq. 'LA') wprime = 'SA'
               if (which .eq. 'SA') wprime = 'LA'
c
               call dsortr (wprime, .true., kplusp, ritz, bounds)
c
            end if
c
c           %--------------------------------------------------%
c           | Scale the Ritz estimate of each Ritz value       |
c           | by 1 / max(eps23,magnitude of the Ritz value).   |
c           %--------------------------------------------------%
c
            do 35 j = 1, nev0
               temp = max( eps23, abs(ritz(j)) )
               bounds(j) = bounds(j)/temp
 35         continue
c
c           %----------------------------------------------------%
c           | Sort the Ritz values according to the scaled Ritz  |
c           | esitmates.  This will push all the converged ones  |
c           | towards the front of ritzr, ritzi, bounds          |
c           | (in the case when NCONV < NEV.)                    |
c           %----------------------------------------------------%
c
            wprime = 'LA'
            call dsortr(wprime, .true., nev0, bounds, ritz)
c
c           %----------------------------------------------%
c           | Scale the Ritz estimate back to its original |
c           | value.                                       |
c           %----------------------------------------------%
c
            do 40 j = 1, nev0
                temp = max( eps23, abs(ritz(j)) )
                bounds(j) = bounds(j)*temp
 40         continue
c
c           %--------------------------------------------------%
c           | Sort the "converged" Ritz values again so that   |
c           | the "threshold" values and their associated Ritz |
c           | estimates appear at the appropriate position in  |
c           | ritz and bound.                                  |
c           %--------------------------------------------------%
c
            if (which .eq. 'BE') then
c
c              %------------------------------------------------%
c              | Sort the "converged" Ritz values in increasing |
c              | order.  The "threshold" values are in the      |
c              | middle.                                        |
c              %------------------------------------------------%
c
               wprime = 'LA'
               call dsortr(wprime, .true., nconv, ritz, bounds)
c
            else
c
c              %----------------------------------------------%
c              | In LM, SM, LA, SA case, sort the "converged" |
c              | Ritz values according to WHICH so that the   |
c              | "threshold" value appears at the front of    |
c              | ritz.                                        |
c              %----------------------------------------------%

               call dsortr(which, .true., nconv, ritz, bounds)
c
            end if
c
c           %------------------------------------------%
c           |  Use h( 1,1 ) as storage to communicate  |
c           |  rnorm to _seupd if needed               |
c           %------------------------------------------%
c
            h(1,1) = rnorm
c
            if (msglvl .gt. 1) then
               call dvout (logfil, kplusp, ritz, ndigit,
     &            '_saup2: Sorted Ritz values.')
               call dvout (logfil, kplusp, bounds, ndigit,
     &            '_saup2: Sorted ritz estimates.')
            end if
c
c           %------------------------------------%
c           | Max iterations have been exceeded. | 
c           %------------------------------------%
c
            if (iter .gt. mxiter .and. nconv .lt. nev) info = 1
c
c           %---------------------%
c           | No shifts to apply. | 
c           %---------------------%
c
            if (np .eq. 0 .and. nconv .lt. nev0) info = 2
c
            np = nconv
            go to 1100
c
         else if (nconv .lt. nev .and. ishift .eq. 1) then
c
c           %---------------------------------------------------%
c           | Do not have all the requested eigenvalues yet.    |
c           | To prevent possible stagnation, adjust the number |
c           | of Ritz values and the shifts.                    |
c           %---------------------------------------------------%
c
            nevbef = nev
            nev = nev + min (nconv, np/2)
            if (nev .eq. 1 .and. kplusp .ge. 6) then
               nev = kplusp / 2
            else if (nev .eq. 1 .and. kplusp .gt. 2) then
               nev = 2
            end if
            np  = kplusp - nev
c     
c           %---------------------------------------%
c           | If the size of NEV was just increased |
c           | resort the eigenvalues.               |
c           %---------------------------------------%
c     
            if (nevbef .lt. nev) 
     &         call dsgets (ishift, which, nev, np, ritz, bounds,
     &              workl)
c
         end if
c
         if (msglvl .gt. 0) then
            call ivout (logfil, 1, nconv, ndigit,
     &           '_saup2: no. of "converged" Ritz values at this iter.')
            if (msglvl .gt. 1) then
               kp(1) = nev
               kp(2) = np
               call ivout (logfil, 2, kp, ndigit,
     &              '_saup2: NEV and NP are')
               call dvout (logfil, nev, ritz(np+1), ndigit,
     &              '_saup2: "wanted" Ritz values.')
               call dvout (logfil, nev, bounds(np+1), ndigit,
     &              '_saup2: Ritz estimates of the "wanted" values ')
            end if
         end if

c 
         if (ishift .eq. 0) then
c
c           %-----------------------------------------------------%
c           | User specified shifts: reverse communication to     |
c           | compute the shifts. They are returned in the first  |
c           | NP locations of WORKL.                              |
c           %-----------------------------------------------------%
c
            ushift = .true.
            ido = 3
            go to 9000
         end if
c
   50    continue
c
c        %------------------------------------%
c        | Back from reverse communication;   |
c        | User specified shifts are returned |
c        | in WORKL(1:*NP)                   |
c        %------------------------------------%
c
         ushift = .false.
c 
c 
c        %---------------------------------------------------------%
c        | Move the NP shifts to the first NP locations of RITZ to |
c        | free up WORKL.  This is for the non-exact shift case;   |
c        | in the exact shift case, dsgets already handles this.   |
c        %---------------------------------------------------------%
c
         if (ishift .eq. 0) call dcopy (np, workl, 1, ritz, 1)
c
         if (msglvl .gt. 2) then
            call ivout (logfil, 1, np, ndigit,
     &                  '_saup2: The number of shifts to apply ')
            call dvout (logfil, np, workl, ndigit,
     &                  '_saup2: shifts selected')
            if (ishift .eq. 1) then
               call dvout (logfil, np, bounds, ndigit,
     &                  '_saup2: corresponding Ritz estimates')
             end if
         end if
c 
c        %---------------------------------------------------------%
c        | Apply the NP0 implicit shifts by QR bulge chasing.      |
c        | Each shift is applied to the entire tridiagonal matrix. |
c        | The first 2*N locations of WORKD are used as workspace. |
c        | After dsapps is done, we have a Lanczos                 |
c        | factorization of length NEV.                            |
c        %---------------------------------------------------------%
c
         call dsapps (n, nev, np, ritz, v, ldv, h, ldh, resid, q, ldq,
     &        workd)
c
c        %---------------------------------------------%
c        | Compute the B-norm of the updated residual. |
c        | Keep B*RESID in WORKD(1:N) to be used in    |
c        | the first step of the next call to dsaitr.  |
c        %---------------------------------------------%
c
         cnorm = .true.
         call second (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy (n, resid, 1, workd(n+1), 1)
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
            call dcopy (n, resid, 1, workd, 1)
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
            rnorm = ddot (n, resid, 1, workd, 1)
            rnorm = sqrt(abs(rnorm))
         else if (bmat .eq. 'I') then
            rnorm = dnrm2(n, resid, 1)
         end if
         cnorm = .false.
  130    continue
c
         if (msglvl .gt. 2) then
            call dvout (logfil, 1, rnorm, ndigit, 
     &      '_saup2: B-norm of residual for NEV factorization')
            call dvout (logfil, nev, h(1,2), ndigit,
     &           '_saup2: main diagonal of compressed H matrix')
            call dvout (logfil, nev-1, h(2,1), ndigit,
     &           '_saup2: subdiagonal of compressed H matrix')
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
c     | Error exit |
c     %------------%
c
      call second (t1)
      tsaup2 = t1 - t0
c 
 9000 continue
      return
c
c     %---------------%
c     | End of dsaup2 |
c     %---------------%
c
      end
