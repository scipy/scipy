c\BeginDoc
c 
c\Name: cneupd 
c 
c\Description: 
c  This subroutine returns the converged approximations to eigenvalues 
c  of A*z = lambda*B*z and (optionally): 
c 
c      (1) The corresponding approximate eigenvectors; 
c 
c      (2) An orthonormal basis for the associated approximate 
c          invariant subspace; 
c 
c      (3) Both.  
c
c  There is negligible additional cost to obtain eigenvectors.  An orthonormal 
c  basis is always computed.  There is an additional storage cost of n*nev
c  if both are requested (in this case a separate array Z must be supplied). 
c
c  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
c  are derived from approximate eigenvalues and eigenvectors of
c  of the linear operator OP prescribed by the MODE selection in the
c  call to CNAUPD.  CNAUPD must be called before this routine is called.
c  These approximate eigenvalues and vectors are commonly called Ritz
c  values and Ritz vectors respectively.  They are referred to as such 
c  in the comments that follow.   The computed orthonormal basis for the 
c  invariant subspace corresponding to these Ritz values is referred to as a 
c  Schur basis. 
c 
c  The definition of OP as well as other terms and the relation of computed
c  Ritz values and vectors of OP with respect to the given problem
c  A*z = lambda*B*z may be found in the header of CNAUPD.  For a brief 
c  description, see definitions of IPARAM(7), MODE and WHICH in the
c  documentation of CNAUPD.
c
c\Usage:
c  call cneupd 
c     ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, WORKEV, BMAT, 
c       N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, 
c       WORKL, LWORKL, RWORK, INFO )
c
c\Arguments:
c  RVEC    LOGICAL  (INPUT)
c          Specifies whether a basis for the invariant subspace corresponding
c          to the converged Ritz value approximations for the eigenproblem 
c          A*z = lambda*B*z is computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute Ritz vectors or Schur vectors.
c                                See Remarks below.
c
c  HOWMNY  Character*1  (INPUT)
c          Specifies the form of the basis for the invariant subspace 
c          corresponding to the converged Ritz values that is to be computed.
c
c          = 'A': Compute NEV Ritz vectors;
c          = 'P': Compute NEV Schur vectors;
c          = 'S': compute some of the Ritz vectors, specified
c                 by the logical array SELECT.
c
c  SELECT  Logical array of dimension NCV.  (INPUT)
c          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c          computed. To select the  Ritz vector corresponding to a
c          Ritz value D(j), SELECT(j) must be set to .TRUE.. 
c          If HOWMNY = 'A' or 'P', SELECT need not be initialized 
c          but it is used as internal workspace.
c
c  D       Complex array of dimension NEV+1.  (OUTPUT)
c          On exit, D contains the  Ritz  approximations 
c          to the eigenvalues lambda for A*z = lambda*B*z.
c
c  Z       Complex N by NEV array    (OUTPUT)
c          On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns of 
c          Z represents approximate eigenvectors (Ritz vectors) corresponding 
c          to the NCONV=IPARAM(5) Ritz values for eigensystem
c          A*z = lambda*B*z.
c
c          If RVEC = .FALSE. or HOWMNY = 'P', then Z is NOT REFERENCED.
c
c          NOTE: If if RVEC = .TRUE. and a Schur basis is not required, 
c          the array Z may be set equal to first NEV+1 columns of the Arnoldi 
c          basis array V computed by CNAUPD.  In this case the Arnoldi basis 
c          will be destroyed and overwritten with the eigenvector basis.
c
c  LDZ     Integer.  (INPUT)
c          The leading dimension of the array Z.  If Ritz vectors are
c          desired, then  LDZ .ge.  max( 1, N ) is required.  
c          In any case,  LDZ .ge. 1 is required.
c
c  SIGMA   Complex  (INPUT)
c          If IPARAM(7) = 3 then SIGMA represents the shift. 
c          Not referenced if IPARAM(7) = 1 or 2.
c
c  WORKEV  Complex work array of dimension 2*NCV.  (WORKSPACE)
c
c  **** The remaining arguments MUST be the same as for the   ****
c  **** call to CNAUPD that was just completed.               ****
c
c  NOTE: The remaining arguments 
c
c           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, 
c           WORKD, WORKL, LWORKL, RWORK, INFO 
c
c         must be passed directly to CNEUPD following the last call 
c         to CNAUPD.  These arguments MUST NOT BE MODIFIED between
c         the the last call to CNAUPD and the call to CNEUPD.
c
c  Three of these parameters (V, WORKL and INFO) are also output parameters:
c
c  V       Complex N by NCV array.  (INPUT/OUTPUT)
c
c          Upon INPUT: the NCV columns of V contain the Arnoldi basis
c                      vectors for OP as constructed by CNAUPD .
c
c          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
c                       contain approximate Schur vectors that span the
c                       desired invariant subspace.
c
c          NOTE: If the array Z has been set equal to first NEV+1 columns
c          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then the
c          Arnoldi basis held by V has been overwritten by the desired
c          Ritz vectors.  If a separate array Z has been passed then
c          the first NCONV=IPARAM(5) columns of V will contain approximate
c          Schur vectors that span the desired invariant subspace.
c
c  WORKL   Real work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          WORKL(1:ncv*ncv+2*ncv) contains information obtained in
c          cnaupd.  They are not changed by cneupd.
c          WORKL(ncv*ncv+2*ncv+1:3*ncv*ncv+4*ncv) holds the
c          untransformed Ritz values, the untransformed error estimates of 
c          the Ritz values, the upper triangular matrix for H, and the
c          associated matrix representation of the invariant subspace for H.
c
c          Note: IPNTR(9:13) contains the pointer into WORKL for addresses
c          of the above information computed by cneupd.
c          -------------------------------------------------------------
c          IPNTR(9):  pointer to the NCV RITZ values of the
c                     original system.
c          IPNTR(10): Not used
c          IPNTR(11): pointer to the NCV corresponding error estimates.
c          IPNTR(12): pointer to the NCV by NCV upper triangular
c                     Schur matrix for H.
c          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
c                     of the upper Hessenberg matrix H. Only referenced by
c                     cneupd if RVEC = .TRUE. See Remark 2 below.
c          -------------------------------------------------------------
c
c  INFO    Integer.  (OUTPUT)
c          Error flag on output.
c          =  0: Normal exit.
c
c          =  1: The Schur form computed by LAPACK routine csheqr
c                could not be reordered by LAPACK routine ctrsen.
c                Re-enter subroutine cneupd with IPARAM(5)=NCV and
c                increase the size of the array D to have
c                dimension at least dimension NCV and allocate at least NCV
c                columns for Z. NOTE: Not necessary if Z and V share
c                the same space. Please notify the authors if this error
c                occurs.
c
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from LAPACK eigenvalue calculation.
c                This should never happened.
c          = -9: Error return from calculation of eigenvectors.
c                Informational error from LAPACK routine ctrevc.
c          = -10: IPARAM(7) must be 1,2,3
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: HOWMNY = 'S' not yet implemented
c          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.
c          = -14: CNAUPD did not find any eigenvalues to sufficient
c                 accuracy.
c          = -15: CNEUPD got a different count of the number of converged
c                 Ritz values than CNAUPD got.  This indicates the user
c                 probably made an error in passing data from CNAUPD to
c                 CNEUPD or that the data was modified before entering
c                 CNEUPD
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
c  3. B. Nour-Omid, B. N. Parlett, T. Ericsson and P. S. Jensen,
c     "How to Implement the Spectral Transformation", Math Comp.,
c     Vol. 48, No. 178, April, 1987 pp. 664-673. 
c
c\Routines called:
c     ivout   ARPACK utility routine that prints integers.
c     cmout   ARPACK utility routine that prints matrices
c     cvout   ARPACK utility routine that prints vectors.
c     cgeqr2  LAPACK routine that computes the QR factorization of 
c             a matrix.
c     clacpy  LAPACK matrix copy routine.
c     clahqr  LAPACK routine that computes the Schur form of a
c             upper Hessenberg matrix.
c     claset  LAPACK matrix initialization routine.
c     ctrevc  LAPACK routine to compute the eigenvectors of a matrix
c             in upper triangular form.
c     ctrsen  LAPACK routine that re-orders the Schur form.
c     cunm2r  LAPACK routine that applies an orthogonal matrix in 
c             factored form.
c     slamch  LAPACK routine that determines machine constants.
c     ctrmm   Level 3 BLAS matrix times an upper triangular matrix.
c     cgeru   Level 2 BLAS rank one update to a matrix.
c     ccopy   Level 1 BLAS that copies one vector to another .
c     cscal   Level 1 BLAS that scales a vector.
c     csscal  Level 1 BLAS that scales a complex vector by a real number.
c     scnrm2  Level 1 BLAS that computes the norm of a complex vector.
c
c\Remarks
c
c  1. Currently only HOWMNY = 'A' and 'P' are implemented. 
c
c  2. Schur vectors are an orthogonal representation for the basis of
c     Ritz vectors. Thus, their numerical properties are often superior.
c     If RVEC = .true. then the relationship
c             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
c       transpose( V(:,1:IPARAM(5)) ) * V(:,1:IPARAM(5)) = I
c     are approximately satisfied.
c     Here T is the leading submatrix of order IPARAM(5) of the 
c     upper triangular matrix stored workl(ipntr(12)). 
c
c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Chao Yang                    Houston, Texas 
c     Dept. of Computational & 
c     Applied Mathematics 
c     Rice University 
c     Houston, Texas
c
c\SCCS Information: @(#)
c FILE: neupd.F   SID: 2.7   DATE OF SID: 09/20/00   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
      subroutine cneupd(rvec , howmny, select, d     ,
     &                   z    , ldz   , sigma , workev,
     &                   bmat , n     , which , nev   ,
     &                   tol  , resid , ncv   , v     ,
     &                   ldv  , iparam, ipntr , workd ,
     &                   workl, lworkl, rwork , info  )
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
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Complex     
     &           sigma
      Real 
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(11), ipntr(14)
      logical    select(ncv)
      Real
     &           rwork(ncv)
      Complex
     &           d(nev)     , resid(n)     , v(ldv,ncv),
     &           z(ldz, nev), 
     &           workd(3*n) , workl(lworkl), workev(2*ncv)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex
     &           one, zero
      parameter  (one = (1.0E+0, 0.0E+0), zero = (0.0E+0, 0.0E+0))
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character  type*6
      integer    bounds, ierr  , ih    , ihbds, iheig , nconv ,
     &           invsub, iuptri, iwev  , j    , ldh   , ldq   ,
     &           mode  , msglvl, ritz  , wr   , k     , irz   ,
     &           ibd   , outncv, iq    , np   , numcnv, jj    ,
     &           ishift, nconv2
      Complex 
     &           rnorm, temp, vl(1)
      Real
     &           conds, sep, rtemp, eps23
      logical    reord
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   ccopy , cgeru, cgeqr2, clacpy, cmout,
     &           cunm2r, ctrmm, cvout, ivout,
     &           clahqr
c  
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Real
     &           scnrm2, slamch, slapy2
      external   scnrm2, slamch, slapy2
c
      Complex
     &           wcdotc
      external   wcdotc
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
c     %------------------------%
c     | Set default parameters |
c     %------------------------%
c
      msglvl = mceupd
      mode = iparam(7)
      nconv = iparam(5)
      info = 0
c
c
c     %---------------------------------%
c     | Get machine dependent constant. |
c     %---------------------------------%
c
      eps23 = slamch('Epsilon-Machine')
      eps23 = eps23**(2.0E+0 / 3.0E+0)
c
c     %-------------------------------%
c     | Quick return                  |
c     | Check for incompatible input  |
c     %-------------------------------%
c
      ierr = 0
c
      if (nconv .le. 0) then
         ierr = -14
      else if (n .le. 0) then
         ierr = -1
      else if (nev .le. 0) then
         ierr = -2
      else if (ncv .le. nev+1 .or.  ncv .gt. n) then
         ierr = -3
      else if (which .ne. 'LM' .and.
     &        which .ne. 'SM' .and.
     &        which .ne. 'LR' .and.
     &        which .ne. 'SR' .and.
     &        which .ne. 'LI' .and.
     &        which .ne. 'SI') then
         ierr = -5
      else if (bmat .ne. 'I' .and. bmat .ne. 'G') then
         ierr = -6
      else if (lworkl .lt. 3*ncv**2 + 4*ncv) then
         ierr = -7
      else if ( (howmny .ne. 'A' .and.
     &           howmny .ne. 'P' .and.
     &           howmny .ne. 'S') .and. rvec ) then
         ierr = -13
      else if (howmny .eq. 'S' ) then
         ierr = -12
      end if
c     
      if (mode .eq. 1 .or. mode .eq. 2) then
         type = 'REGULR'
      else if (mode .eq. 3 ) then
         type = 'SHIFTI'
      else 
                                              ierr = -10
      end if
      if (mode .eq. 1 .and. bmat .eq. 'G')    ierr = -11
c
c     %------------%
c     | Error Exit |
c     %------------%
c
      if (ierr .ne. 0) then
         info = ierr
         go to 9000
      end if
c 
c     %--------------------------------------------------------%
c     | Pointer into WORKL for address of H, RITZ, WORKEV, Q   |
c     | etc... and the remaining workspace.                    |
c     | Also update pointer to be used on output.              |
c     | Memory is laid out as follows:                         |
c     | workl(1:ncv*ncv) := generated Hessenberg matrix        |
c     | workl(ncv*ncv+1:ncv*ncv+ncv) := ritz values            |
c     | workl(ncv*ncv+ncv+1:ncv*ncv+2*ncv) := error bounds     |
c     %--------------------------------------------------------%
c
c     %-----------------------------------------------------------%
c     | The following is used and set by CNEUPD.                 |
c     | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := The untransformed |
c     |                                      Ritz values.         |
c     | workl(ncv*ncv+3*ncv+1:ncv*ncv+4*ncv) := The untransformed |
c     |                                      error bounds of      |
c     |                                      the Ritz values      |
c     | workl(ncv*ncv+4*ncv+1:2*ncv*ncv+4*ncv) := Holds the upper |
c     |                                      triangular matrix    |
c     |                                      for H.               |
c     | workl(2*ncv*ncv+4*ncv+1: 3*ncv*ncv+4*ncv) := Holds the    |
c     |                                      associated matrix    |
c     |                                      representation of    |
c     |                                      the invariant        |
c     |                                      subspace for H.      |
c     | GRAND total of NCV * ( 3 * NCV + 4 ) locations.           |
c     %-----------------------------------------------------------%
c     
      ih     = ipntr(5)
      ritz   = ipntr(6)
      iq     = ipntr(7)
      bounds = ipntr(8)
      ldh    = ncv
      ldq    = ncv
      iheig  = bounds + ldh
      ihbds  = iheig  + ldh
      iuptri = ihbds  + ldh
      invsub = iuptri + ldh*ncv
      ipntr(9)  = iheig
      ipntr(11) = ihbds
      ipntr(12) = iuptri
      ipntr(13) = invsub
      wr = 1
      iwev = wr + ncv
c
c     %-----------------------------------------%
c     | irz points to the Ritz values computed  |
c     |     by _neigh before exiting _naup2.    |
c     | ibd points to the Ritz estimates        |
c     |     computed by _neigh before exiting   |
c     |     _naup2.                             |
c     %-----------------------------------------%
c
      irz = ipntr(14) + ncv*ncv
      ibd = irz + ncv
c
c     %------------------------------------%
c     | RNORM is B-norm of the RESID(1:N). |
c     %------------------------------------%
c
      rnorm = workl(ih+2)
      workl(ih+2) = zero
c
      if (msglvl .gt. 2) then
         call cvout(logfil, ncv, workl(irz), ndigit,
     &   '_neupd: Ritz values passed in from _NAUPD.')
         call cvout(logfil, ncv, workl(ibd), ndigit,
     &   '_neupd: Ritz estimates passed in from _NAUPD.')
      end if
c
      if (rvec) then
c
         reord = .false.
c
c        %---------------------------------------------------%
c        | Use the temporary bounds array to store indices   |
c        | These will be used to mark the select array later |
c        %---------------------------------------------------%
c
         do 10 j = 1,ncv
            workl(bounds+j-1) = j
            select(j) = .false.
   10    continue
c
c        %-------------------------------------%
c        | Select the wanted Ritz values.      |
c        | Sort the Ritz values so that the    |
c        | wanted ones appear at the tailing   |
c        | NEV positions of workl(irr) and     |
c        | workl(iri).  Move the corresponding |
c        | error estimates in workl(ibd)       |
c        | accordingly.                        |
c        %-------------------------------------%
c
         np     = ncv - nev
         ishift = 0
         call cngets(ishift, which     , nev          ,
     &                np    , workl(irz), workl(bounds))
c
         if (msglvl .gt. 2) then
            call cvout (logfil, ncv, workl(irz), ndigit,
     &      '_neupd: Ritz values after calling _NGETS.')
            call cvout (logfil, ncv, workl(bounds), ndigit,
     &      '_neupd: Ritz value indices after calling _NGETS.')
         end if
c
c        %-----------------------------------------------------%
c        | Record indices of the converged wanted Ritz values  |
c        | Mark the select array for possible reordering       |
c        %-----------------------------------------------------%
c
         numcnv = 0
         do 11 j = 1,ncv
            rtemp = max(eps23,
     &                 slapy2 ( real(workl(irz+ncv-j)),
     &                          aimag(workl(irz+ncv-j)) ))
            jj = workl(bounds + ncv - j)
            if (numcnv .lt. nconv .and.
     &          slapy2( real(workl(ibd+jj-1)),
     &          aimag(workl(ibd+jj-1)) )
     &          .le. tol*rtemp) then
               select(jj) = .true.
               numcnv = numcnv + 1
               if (jj .gt. nconv) reord = .true.
            endif
   11    continue
c
c        %-----------------------------------------------------------%
c        | Check the count (numcnv) of converged Ritz values with    |
c        | the number (nconv) reported by dnaupd.  If these two      |
c        | are different then there has probably been an error       |
c        | caused by incorrect passing of the dnaupd data.           |
c        %-----------------------------------------------------------%
c
         if (msglvl .gt. 2) then
             call ivout(logfil, 1, numcnv, ndigit,
     &            '_neupd: Number of specified eigenvalues')
             call ivout(logfil, 1, nconv, ndigit,
     &            '_neupd: Number of "converged" eigenvalues')
         end if
c
         if (numcnv .ne. nconv) then
            info = -15
            go to 9000
         end if
c
c        %-------------------------------------------------------%
c        | Call LAPACK routine clahqr to compute the Schur form |
c        | of the upper Hessenberg matrix returned by CNAUPD.   |
c        | Make a copy of the upper Hessenberg matrix.           |
c        | Initialize the Schur vector matrix Q to the identity. |
c        %-------------------------------------------------------%
c
         call ccopy(ldh*ncv, workl(ih), 1, workl(iuptri), 1)
         call claset('All', ncv, ncv          , 
     &                zero , one, workl(invsub),
     &                ldq)
         call clahqr(.true., .true.       , ncv          , 
     &                1     , ncv          , workl(iuptri),
     &                ldh   , workl(iheig) , 1            ,
     &                ncv   , workl(invsub), ldq          ,
     &                ierr)
         call ccopy(ncv         , workl(invsub+ncv-1), ldq,
     &               workl(ihbds), 1)
c
         if (ierr .ne. 0) then
            info = -8
            go to 9000
         end if
c
         if (msglvl .gt. 1) then
            call cvout (logfil, ncv, workl(iheig), ndigit,
     &           '_neupd: Eigenvalues of H')
            call cvout (logfil, ncv, workl(ihbds), ndigit,
     &           '_neupd: Last row of the Schur vector matrix')
            if (msglvl .gt. 3) then
               call cmout (logfil       , ncv, ncv   , 
     &                     workl(iuptri), ldh, ndigit,
     &              '_neupd: The upper triangular matrix ')
            end if
         end if
c
         if (reord) then
c
c           %-----------------------------------------------%
c           | Reorder the computed upper triangular matrix. |
c           %-----------------------------------------------%
c
            call ctrsen('None'       , 'V'          , select      ,
     &                   ncv          , workl(iuptri), ldh         ,
     &                   workl(invsub), ldq          , workl(iheig),
     &                   nconv2       , conds        , sep         , 
     &                   workev       , ncv          , ierr)
c
            if (nconv2 .lt. nconv) then
               nconv = nconv2
            end if

            if (ierr .eq. 1) then
               info = 1
               go to 9000
            end if
c
            if (msglvl .gt. 2) then
                call cvout (logfil, ncv, workl(iheig), ndigit,
     &           '_neupd: Eigenvalues of H--reordered')
                if (msglvl .gt. 3) then
                   call cmout(logfil       , ncv, ncv   ,
     &                         workl(iuptri), ldq, ndigit,
     &              '_neupd: Triangular matrix after re-ordering')
                end if
            end if
c
         end if
c
c        %---------------------------------------------%
c        | Copy the last row of the Schur basis matrix |
c        | to workl(ihbds).  This vector will be used  |
c        | to compute the Ritz estimates of converged  |
c        | Ritz values.                                |
c        %---------------------------------------------%
c
         call ccopy(ncv         , workl(invsub+ncv-1), ldq,
     &               workl(ihbds), 1)
c 
c        %--------------------------------------------%
c        | Place the computed eigenvalues of H into D |
c        | if a spectral transformation was not used. |
c        %--------------------------------------------%
c
         if (type .eq. 'REGULR') then
            call ccopy(nconv, workl(iheig), 1, d, 1)
         end if
c
c        %----------------------------------------------------------%
c        | Compute the QR factorization of the matrix representing  |
c        | the wanted invariant subspace located in the first NCONV |
c        | columns of workl(invsub,ldq).                            |
c        %----------------------------------------------------------%
c
         call cgeqr2(ncv , nconv , workl(invsub),
     &                ldq , workev, workev(ncv+1),
     &                ierr)
c
c        %--------------------------------------------------------%
c        | * Postmultiply V by Q using cunm2r.                    |
c        | * Copy the first NCONV columns of VQ into Z.           |
c        | * Postmultiply Z by R.                                 |
c        | The N by NCONV matrix Z is now a matrix representation |
c        | of the approximate invariant subspace associated with  |
c        | the Ritz values in workl(iheig). The first NCONV       | 
c        | columns of V are now approximate Schur vectors         |
c        | associated with the upper triangular matrix of order   |
c        | NCONV in workl(iuptri).                                |
c        %--------------------------------------------------------%
c
         call cunm2r('Right', 'Notranspose', n            ,
     &                ncv    , nconv        , workl(invsub),
     &                ldq    , workev       , v            ,
     &                ldv    , workd(n+1)   , ierr)
         call clacpy('All', n, nconv, v, ldv, z, ldz)
c
         do 20 j=1, nconv
c
c           %---------------------------------------------------%
c           | Perform both a column and row scaling if the      |
c           | diagonal element of workl(invsub,ldq) is negative |
c           | I'm lazy and don't take advantage of the upper    |
c           | triangular form of workl(iuptri,ldq).             |
c           | Note that since Q is orthogonal, R is a diagonal  |
c           | matrix consisting of plus or minus ones.          |
c           %---------------------------------------------------%
c
            if ( real( workl(invsub+(j-1)*ldq+j-1) ) .lt. 
     &                  real(zero) ) then
               call cscal(nconv, -one, workl(iuptri+j-1), ldq)
               call cscal(nconv, -one, workl(iuptri+(j-1)*ldq), 1)
            end if
c
 20      continue
c
         if (howmny .eq. 'A') then
c
c           %--------------------------------------------%
c           | Compute the NCONV wanted eigenvectors of T |
c           | located in workl(iuptri,ldq).              |
c           %--------------------------------------------%
c
            do 30 j=1, ncv
               if (j .le. nconv) then
                  select(j) = .true.
               else
                  select(j) = .false.
               end if
 30         continue
c
            call ctrevc('Right', 'Select'     , select       ,
     &                   ncv    , workl(iuptri), ldq          ,
     &                   vl     , 1            , workl(invsub),
     &                   ldq    , ncv          , outncv       ,
     &                   workev , rwork        , ierr)
c
            if (ierr .ne. 0) then
                info = -9
                go to 9000
            end if
c
c           %------------------------------------------------%
c           | Scale the returning eigenvectors so that their |
c           | Euclidean norms are all one. LAPACK subroutine |
c           | ctrevc returns each eigenvector normalized so  |
c           | that the element of largest magnitude has      |
c           | magnitude 1.                                   |
c           %------------------------------------------------%
c
            do 40 j=1, nconv
                  rtemp = scnrm2(ncv, workl(invsub+(j-1)*ldq), 1)
                  rtemp = real(one) / rtemp
                  call csscal ( ncv, rtemp,
     &                 workl(invsub+(j-1)*ldq), 1 )
c
c                 %------------------------------------------%
c                 | Ritz estimates can be obtained by taking |
c                 | the inner product of the last row of the |
c                 | Schur basis of H with eigenvectors of T. |
c                 | Note that the eigenvector matrix of T is |
c                 | upper triangular, thus the length of the |
c                 | inner product can be set to j.           |
c                 %------------------------------------------%
c 
                  workev(j) = wcdotc(j, workl(ihbds), 1,
     &                        workl(invsub+(j-1)*ldq), 1)
 40         continue
c
            if (msglvl .gt. 2) then
               call ccopy(nconv, workl(invsub+ncv-1), ldq,
     &                    workl(ihbds), 1)
               call cvout (logfil, nconv, workl(ihbds), ndigit,
     &            '_neupd: Last row of the eigenvector matrix for T')
               if (msglvl .gt. 3) then
                  call cmout(logfil       , ncv, ncv   ,
     &                        workl(invsub), ldq, ndigit,
     &               '_neupd: The eigenvector matrix for T')
               end if
            end if
c
c           %---------------------------------------%
c           | Copy Ritz estimates into workl(ihbds) |
c           %---------------------------------------%
c 
            call ccopy(nconv, workev, 1, workl(ihbds), 1)
c
c           %----------------------------------------------%
c           | The eigenvector matrix Q of T is triangular. |
c           | Form Z*Q.                                    |
c           %----------------------------------------------%
c
            call ctrmm('Right'   , 'Upper'      , 'No transpose',
     &                  'Non-unit', n            , nconv         ,
     &                  one       , workl(invsub), ldq           ,
     &                  z         , ldz)
         end if 
c
      else
c
c        %--------------------------------------------------%
c        | An approximate invariant subspace is not needed. |
c        | Place the Ritz values computed CNAUPD into D.    |
c        %--------------------------------------------------%
c
         call ccopy(nconv, workl(ritz), 1, d, 1)
         call ccopy(nconv, workl(ritz), 1, workl(iheig), 1)
         call ccopy(nconv, workl(bounds), 1, workl(ihbds), 1)
c
      end if
c
c     %------------------------------------------------%
c     | Transform the Ritz values and possibly vectors |
c     | and corresponding error bounds of OP to those  |
c     | of A*x = lambda*B*x.                           |
c     %------------------------------------------------%
c
      if (type .eq. 'REGULR') then
c
         if (rvec) 
     &      call cscal(ncv, rnorm, workl(ihbds), 1)
c      
      else
c     
c        %---------------------------------------%
c        |   A spectral transformation was used. |
c        | * Determine the Ritz estimates of the |
c        |   Ritz values in the original system. |
c        %---------------------------------------%
c
         if (rvec) 
     &      call cscal(ncv, rnorm, workl(ihbds), 1)
c    
         do 50 k=1, ncv
            temp = workl(iheig+k-1)
            workl(ihbds+k-1) = workl(ihbds+k-1) / temp / temp
  50     continue
c  
      end if
c
c     %-----------------------------------------------------------%
c     | *  Transform the Ritz values back to the original system. |
c     |    For TYPE = 'SHIFTI' the transformation is              |
c     |             lambda = 1/theta + sigma                      |
c     | NOTES:                                                    |
c     | *The Ritz vectors are not affected by the transformation. |
c     %-----------------------------------------------------------%
c    
      if (type .eq. 'SHIFTI') then
         do 60 k=1, nconv
            d(k) = one / workl(iheig+k-1) + sigma
  60     continue
      end if
c
      if (type .ne. 'REGULR' .and. msglvl .gt. 1) then
         call cvout (logfil, nconv, d, ndigit,
     &     '_neupd: Untransformed Ritz values.')
         call cvout (logfil, nconv, workl(ihbds), ndigit,
     &     '_neupd: Ritz estimates of the untransformed Ritz values.')
      else if ( msglvl .gt. 1) then
         call cvout (logfil, nconv, d, ndigit,
     &     '_neupd: Converged Ritz values.')
         call cvout (logfil, nconv, workl(ihbds), ndigit,
     &     '_neupd: Associated Ritz estimates.')
      end if
c
c     %-------------------------------------------------%
c     | Eigenvector Purification step. Formally perform |
c     | one of inverse subspace iteration. Only used    |
c     | for MODE = 3. See reference 3.                  |
c     %-------------------------------------------------%
c
      if (rvec .and. howmny .eq. 'A' .and. type .eq. 'SHIFTI') then
c
c        %------------------------------------------------%
c        | Purify the computed Ritz vectors by adding a   |
c        | little bit of the residual vector:             |
c        |                      T                         |
c        |          resid(:)*( e    s ) / theta           |
c        |                      NCV                       |
c        | where H s = s theta.                           |
c        %------------------------------------------------%
c
         do 100 j=1, nconv
            if (workl(iheig+j-1) .ne. zero) then
               workev(j) =  workl(invsub+(j-1)*ldq+ncv-1) /
     &                      workl(iheig+j-1)
            endif
 100     continue

c        %---------------------------------------%
c        | Perform a rank one update to Z and    |
c        | purify all the Ritz vectors together. |
c        %---------------------------------------%
c
         call cgeru (n, nconv, one, resid, 1, workev, 1, z, ldz)
c
      end if
c
 9000 continue
c
      return
c     
c     %---------------%
c     | End of cneupd|
c     %---------------%
c
      end
