c\BeginDoc
c
c\Name: dneupd
c
c\Description:
c
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
c  call to DNAUPD .  DNAUPD  must be called before this routine is called.
c  These approximate eigenvalues and vectors are commonly called Ritz
c  values and Ritz vectors respectively.  They are referred to as such
c  in the comments that follow.  The computed orthonormal basis for the
c  invariant subspace corresponding to these Ritz values is referred to as a
c  Schur basis.
c
c  See documentation in the header of the subroutine DNAUPD  for
c  definition of OP as well as other terms and the relation of computed
c  Ritz values and Ritz vectors of OP with respect to the given problem
c  A*z = lambda*B*z.  For a brief description, see definitions of
c  IPARAM(7), MODE and WHICH in the documentation of DNAUPD .
c
c\Usage:
c  call dneupd
c     ( RVEC, HOWMNY, SELECT, DR, DI, Z, LDZ, SIGMAR, SIGMAI, WORKEV, BMAT,
c       N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL,
c       LWORKL, INFO )
c
c\Arguments:
c  RVEC    LOGICAL  (INPUT)
c          Specifies whether a basis for the invariant subspace corresponding
c          to the converged Ritz value approximations for the eigenproblem
c          A*z = lambda*B*z is computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute the Ritz vectors or Schur vectors.
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
c          computed. To select the Ritz vector corresponding to a
c          Ritz value (DR(j), DI(j)), SELECT(j) must be set to .TRUE..
c          If HOWMNY = 'A' or 'P', SELECT is used as internal workspace.
c
c  DR      Double precision  array of dimension NEV+1.  (OUTPUT)
c          If IPARAM(7) = 1,2 or 3 and SIGMAI=0.0  then on exit: DR contains
c          the real part of the Ritz  approximations to the eigenvalues of
c          A*z = lambda*B*z.
c          If IPARAM(7) = 3, 4 and SIGMAI is not equal to zero, then on exit:
c          DR contains the real part of the Ritz values of OP computed by
c          DNAUPD . A further computation must be performed by the user
c          to transform the Ritz values computed for OP by DNAUPD  to those
c          of the original system A*z = lambda*B*z. See remark 3 below.
c
c  DI      Double precision  array of dimension NEV+1.  (OUTPUT)
c          On exit, DI contains the imaginary part of the Ritz value
c          approximations to the eigenvalues of A*z = lambda*B*z associated
c          with DR.
c
c          NOTE: When Ritz values are complex, they will come in complex
c                conjugate pairs.  If eigenvectors are requested, the
c                corresponding Ritz vectors will also come in conjugate
c                pairs and the real and imaginary parts of these are
c                represented in two consecutive columns of the array Z
c                (see below).
c
c  Z       Double precision  N by NEV+1 array if RVEC = .TRUE. and HOWMNY = 'A'. (OUTPUT)
c          On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns of
c          Z represent approximate eigenvectors (Ritz vectors) corresponding
c          to the NCONV=IPARAM(5) Ritz values for eigensystem
c          A*z = lambda*B*z.
c
c          The complex Ritz vector associated with the Ritz value
c          with positive imaginary part is stored in two consecutive
c          columns.  The first column holds the real part of the Ritz
c          vector and the second column holds the imaginary part.  The
c          Ritz vector associated with the Ritz value with negative
c          imaginary part is simply the complex conjugate of the Ritz vector
c          associated with the positive imaginary part.
c
c          If  RVEC = .FALSE. or HOWMNY = 'P', then Z is not referenced.
c
c          NOTE: If if RVEC = .TRUE. and a Schur basis is not required,
c          the array Z may be set equal to first NEV+1 columns of the Arnoldi
c          basis array V computed by DNAUPD .  In this case the Arnoldi basis
c          will be destroyed and overwritten with the eigenvector basis.
c
c  LDZ     Integer.  (INPUT)
c          The leading dimension of the array Z.  If Ritz vectors are
c          desired, then  LDZ >= max( 1, N ).  In any case,  LDZ >= 1.
c
c  SIGMAR  Double precision   (INPUT)
c          If IPARAM(7) = 3 or 4, represents the real part of the shift.
c          Not referenced if IPARAM(7) = 1 or 2.
c
c  SIGMAI  Double precision   (INPUT)
c          If IPARAM(7) = 3 or 4, represents the imaginary part of the shift.
c          Not referenced if IPARAM(7) = 1 or 2. See remark 3 below.
c
c  WORKEV  Double precision  work array of dimension 3*NCV.  (WORKSPACE)
c
c  **** The remaining arguments MUST be the same as for the   ****
c  **** call to DNAUPD  that was just completed.               ****
c
c  NOTE: The remaining arguments
c
c           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
c           WORKD, WORKL, LWORKL, INFO
c
c         must be passed directly to DNEUPD  following the last call
c         to DNAUPD .  These arguments MUST NOT BE MODIFIED between
c         the the last call to DNAUPD  and the call to DNEUPD .
c
c  Three of these parameters (V, WORKL, INFO) are also output parameters:
c
c  V       Double precision  N by NCV array.  (INPUT/OUTPUT)
c
c          Upon INPUT: the NCV columns of V contain the Arnoldi basis
c                      vectors for OP as constructed by DNAUPD  .
c
c          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
c                       contain approximate Schur vectors that span the
c                       desired invariant subspace.  See Remark 2 below.
c
c          NOTE: If the array Z has been set equal to first NEV+1 columns
c          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then the
c          Arnoldi basis held by V has been overwritten by the desired
c          Ritz vectors.  If a separate array Z has been passed then
c          the first NCONV=IPARAM(5) columns of V will contain approximate
c          Schur vectors that span the desired invariant subspace.
c
c  WORKL   Double precision  work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          WORKL(1:ncv*ncv+3*ncv) contains information obtained in
c          dnaupd .  They are not changed by dneupd .
c          WORKL(ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) holds the
c          real and imaginary part of the untransformed Ritz values,
c          the upper quasi-triangular matrix for H, and the
c          associated matrix representation of the invariant subspace for H.
c
c          Note: IPNTR(9:13) contains the pointer into WORKL for addresses
c          of the above information computed by dneupd .
c          -------------------------------------------------------------
c          IPNTR(9):  pointer to the real part of the NCV RITZ values of the
c                     original system.
c          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of
c                     the original system.
c          IPNTR(11): pointer to the NCV corresponding error bounds.
c          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
c                     Schur matrix for H.
c          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
c                     of the upper Hessenberg matrix H. Only referenced by
c                     dneupd  if RVEC = .TRUE. See Remark 2 below.
c          -------------------------------------------------------------
c
c  INFO    Integer.  (OUTPUT)
c          Error flag on output.
c
c          =  0: Normal exit.
c
c          =  1: The Schur form computed by LAPACK routine dlahqr
c                could not be reordered by LAPACK routine dtrsen .
c                Re-enter subroutine dneupd  with IPARAM(5)=NCV and
c                increase the size of the arrays DR and DI to have
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
c          = -8: Error return from calculation of a real Schur form.
c                Informational error from LAPACK routine dlahqr .
c          = -9: Error return from calculation of eigenvectors.
c                Informational error from LAPACK routine dtrevc .
c          = -10: IPARAM(7) must be 1,2,3,4.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: HOWMNY = 'S' not yet implemented
c          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.
c          = -14: DNAUPD  did not find any eigenvalues to sufficient
c                 accuracy.
c          = -15: DNEUPD got a different count of the number of converged
c                 Ritz values than DNAUPD got.  This indicates the user
c                 probably made an error in passing data from DNAUPD to
c                 DNEUPD or that the data was modified before entering
c                 DNEUPD
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
c  3. B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for
c     Real Matrices", Linear Algebra and its Applications, vol 88/89,
c     pp 575-595, (1987).
c
c\Routines called:
c     ivout   ARPACK utility routine that prints integers.
c     dmout    ARPACK utility routine that prints matrices
c     dvout    ARPACK utility routine that prints vectors.
c     dgeqr2   LAPACK routine that computes the QR factorization of
c             a matrix.
c     dlacpy   LAPACK matrix copy routine.
c     dlahqr   LAPACK routine to compute the real Schur form of an
c             upper Hessenberg matrix.
c     dlamch   LAPACK routine that determines machine constants.
c     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dlaset   LAPACK matrix initialization routine.
c     dorm2r   LAPACK routine that applies an orthogonal matrix in
c             factored form.
c     dtrevc   LAPACK routine to compute the eigenvectors of a matrix
c             in upper quasi-triangular form.
c     dtrsen   LAPACK routine that re-orders the Schur form.
c     dtrmm    Level 3 BLAS matrix times an upper triangular matrix.
c     dger     Level 2 BLAS rank one update to a matrix.
c     dcopy    Level 1 BLAS that copies one vector to another .
c     ddot     Level 1 BLAS that computes the scalar product of two vectors.
c     dnrm2    Level 1 BLAS that computes the norm of a vector.
c     dscal    Level 1 BLAS that scales a vector.
c
c\Remarks
c
c  1. Currently only HOWMNY = 'A' and 'P' are implemented.
c
c     Let trans(X) denote the transpose of X.
c
c  2. Schur vectors are an orthogonal representation for the basis of
c     Ritz vectors. Thus, their numerical properties are often superior.
c     If RVEC = .TRUE. then the relationship
c             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
c     trans(V(:,1:IPARAM(5))) * V(:,1:IPARAM(5)) = I are approximately
c     satisfied. Here T is the leading submatrix of order IPARAM(5) of the
c     real upper quasi-triangular matrix stored workl(ipntr(12)). That is,
c     T is block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
c     each 2-by-2 diagonal block has its diagonal elements equal and its
c     off-diagonal elements of opposite sign.  Corresponding to each 2-by-2
c     diagonal block is a complex conjugate pair of Ritz values. The real
c     Ritz values are stored on the diagonal of T.
c
c  3. If IPARAM(7) = 3 or 4 and SIGMAI is not equal zero, then the user must
c     form the IPARAM(5) Rayleigh quotients in order to transform the Ritz
c     values computed by DNAUPD  for OP to those of A*z = lambda*B*z.
c     Set RVEC = .true. and HOWMNY = 'A', and
c     compute
c           trans(Z(:,I)) * A * Z(:,I) if DI(I) = 0.
c     If DI(I) is not equal to zero and DI(I+1) = - D(I),
c     then the desired real and imaginary parts of the Ritz value are
c           trans(Z(:,I)) * A * Z(:,I) +  trans(Z(:,I+1)) * A * Z(:,I+1),
c           trans(Z(:,I)) * A * Z(:,I+1) -  trans(Z(:,I+1)) * A * Z(:,I),
c     respectively.
c     Another possibility is to set RVEC = .true. and HOWMNY = 'P' and
c     compute trans(V(:,1:IPARAM(5))) * A * V(:,1:IPARAM(5)) and then an upper
c     quasi-triangular matrix of order IPARAM(5) is computed. See remark
c     2 above.
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
      subroutine dneupd (rvec , howmny, select, dr    , di,
     &                   z    , ldz   , sigmar, sigmai, workev,
     &                   bmat , n     , which , nev   , tol,
     &                   resid, ncv   , v     , ldv   , iparam,
     &                   ipntr, workd , workl , lworkl, info)
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
      Double precision
     &           sigmar, sigmai, tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(11), ipntr(14)
      logical    select(ncv)
      Double precision
     &           dr(nev+1)    , di(nev+1), resid(n)  ,
     &           v(ldv,ncv)   , z(ldz,*) , workd(3*n),
     &           workl(lworkl), workev(3*ncv)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0 , zero = 0.0D+0 )
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character  type*6
      integer    bounds, ierr  , ih    , ihbds   ,
     &           iheigr, iheigi, iconj , nconv   ,
     &           invsub, iuptri, iwev  , iwork(1),
     &           j     , k     , ldh   , ldq     ,
     &           mode  , msglvl, outncv, ritzr   ,
     &           ritzi , wri   , wrr   , irr     ,
     &           iri   , ibd   , ishift, numcnv  ,
     &           np    , jj    , nconv2
      logical    reord
      Double precision
     &           conds  , rnorm, sep  , temp,
     &           vl(1,1), temp1, eps23
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy  , dger   , dgeqr2 , dlacpy ,
     &           dlahqr , dlaset , dmout  , dorm2r ,
     &           dtrevc , dtrmm  , dtrsen , dscal  ,
     &           dvout  , ivout
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dlapy2 , dnrm2 , dlamch , ddot
      external   dlapy2 , dnrm2 , dlamch , ddot
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    abs, min, sqrt
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %------------------------%
c     | Set default parameters |
c     %------------------------%
c
      msglvl = mneupd
      mode = iparam(7)
      nconv = iparam(5)
      info = 0
c
c     %---------------------------------%
c     | Get machine dependent constant. |
c     %---------------------------------%
c
      eps23 = dlamch ('Epsilon-Machine')
      eps23 = eps23**(2.0D+0  / 3.0D+0 )
c
c     %--------------%
c     | Quick return |
c     %--------------%
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
      else if (lworkl .lt. 3*ncv**2 + 6*ncv) then
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
      else if (mode .eq. 3 .and. sigmai .eq. zero) then
         type = 'SHIFTI'
      else if (mode .eq. 3 ) then
         type = 'REALPT'
      else if (mode .eq. 4 ) then
         type = 'IMAGPT'
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
c     | Pointer into WORKL for address of H, RITZ, BOUNDS, Q   |
c     | etc... and the remaining workspace.                    |
c     | Also update pointer to be used on output.              |
c     | Memory is laid out as follows:                         |
c     | workl(1:ncv*ncv) := generated Hessenberg matrix        |
c     | workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary   |
c     |                                   parts of ritz values |
c     | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds   |
c     %--------------------------------------------------------%
c
c     %-----------------------------------------------------------%
c     | The following is used and set by DNEUPD .                  |
c     | workl(ncv*ncv+3*ncv+1:ncv*ncv+4*ncv) := The untransformed |
c     |                             real part of the Ritz values. |
c     | workl(ncv*ncv+4*ncv+1:ncv*ncv+5*ncv) := The untransformed |
c     |                        imaginary part of the Ritz values. |
c     | workl(ncv*ncv+5*ncv+1:ncv*ncv+6*ncv) := The untransformed |
c     |                           error bounds of the Ritz values |
c     | workl(ncv*ncv+6*ncv+1:2*ncv*ncv+6*ncv) := Holds the upper |
c     |                             quasi-triangular matrix for H |
c     | workl(2*ncv*ncv+6*ncv+1: 3*ncv*ncv+6*ncv) := Holds the    |
c     |       associated matrix representation of the invariant   |
c     |       subspace for H.                                     |
c     | GRAND total of NCV * ( 3 * NCV + 6 ) locations.           |
c     %-----------------------------------------------------------%
c
      ih     = ipntr(5)
      ritzr  = ipntr(6)
      ritzi  = ipntr(7)
      bounds = ipntr(8)
      ldh    = ncv
      ldq    = ncv
      iheigr = bounds + ldh
      iheigi = iheigr + ldh
      ihbds  = iheigi + ldh
      iuptri = ihbds  + ldh
      invsub = iuptri + ldh*ncv
      ipntr(9)  = iheigr
      ipntr(10) = iheigi
      ipntr(11) = ihbds
      ipntr(12) = iuptri
      ipntr(13) = invsub
      wrr = 1
      wri = ncv + 1
      iwev = wri + ncv
c
c     %-----------------------------------------%
c     | irr points to the REAL part of the Ritz |
c     |     values computed by _neigh before    |
c     |     exiting _naup2.                     |
c     | iri points to the IMAGINARY part of the |
c     |     Ritz values computed by _neigh      |
c     |     before exiting _naup2.              |
c     | ibd points to the Ritz estimates        |
c     |     computed by _neigh before exiting   |
c     |     _naup2.                             |
c     %-----------------------------------------%
c
      irr = ipntr(14)+ncv*ncv
      iri = irr+ncv
      ibd = iri+ncv
c
c     %------------------------------------%
c     | RNORM is B-norm of the RESID(1:N). |
c     %------------------------------------%
c
      rnorm = workl(ih+2)
      workl(ih+2) = zero
c
      if (msglvl .gt. 2) then
         call dvout (logfil, ncv, workl(irr), ndigit,
     &   '_neupd: Real part of Ritz values passed in from _NAUPD.')
         call dvout (logfil, ncv, workl(iri), ndigit,
     &   '_neupd: Imag part of Ritz values passed in from _NAUPD.')
         call dvout (logfil, ncv, workl(ibd), ndigit,
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
c        | error estimates in workl(bound)     |
c        | accordingly.                        |
c        %-------------------------------------%
c
         np     = ncv - nev
         ishift = 0
         call dngets (ishift       , which     , nev       ,
     &                np           , workl(irr), workl(iri),
     &                workl(bounds), workl     , workl(np+1))
c
         if (msglvl .gt. 2) then
            call dvout (logfil, ncv, workl(irr), ndigit,
     &      '_neupd: Real part of Ritz values after calling _NGETS.')
            call dvout (logfil, ncv, workl(iri), ndigit,
     &      '_neupd: Imag part of Ritz values after calling _NGETS.')
            call dvout (logfil, ncv, workl(bounds), ndigit,
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
            temp1 = max(eps23,
     &                 dlapy2 ( workl(irr+ncv-j), workl(iri+ncv-j) ))
            jj = workl(bounds + ncv - j)
            if (numcnv .lt. nconv .and.
     &          workl(ibd+jj-1) .le. tol*temp1) then
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
             call ivout(logfil, 1, [numcnv], ndigit,
     &            '_neupd: Number of specified eigenvalues')
             call ivout(logfil, 1, [nconv], ndigit,
     &            '_neupd: Number of "converged" eigenvalues')
         end if
c
         if (numcnv .ne. nconv) then
            info = -15
            go to 9000
         end if
c
c        %-----------------------------------------------------------%
c        | Call LAPACK routine dlahqr  to compute the real Schur form |
c        | of the upper Hessenberg matrix returned by DNAUPD .        |
c        | Make a copy of the upper Hessenberg matrix.               |
c        | Initialize the Schur vector matrix Q to the identity.     |
c        %-----------------------------------------------------------%
c
         call dcopy (ldh*ncv, workl(ih), 1, workl(iuptri), 1)
         call dlaset ('All', ncv, ncv,
     &                zero , one, workl(invsub),
     &                ldq)
         call dlahqr (.true., .true.       , ncv,
     &                1     , ncv          , workl(iuptri),
     &                ldh   , workl(iheigr), workl(iheigi),
     &                1     , ncv          , workl(invsub),
     &                ldq   , ierr)
         call dcopy (ncv         , workl(invsub+ncv-1), ldq,
     &               workl(ihbds), 1)
c
         if (ierr .ne. 0) then
            info = -8
            go to 9000
         end if
c
         if (msglvl .gt. 1) then
            call dvout (logfil, ncv, workl(iheigr), ndigit,
     &           '_neupd: Real part of the eigenvalues of H')
            call dvout (logfil, ncv, workl(iheigi), ndigit,
     &           '_neupd: Imaginary part of the Eigenvalues of H')
            call dvout (logfil, ncv, workl(ihbds), ndigit,
     &           '_neupd: Last row of the Schur vector matrix')
            if (msglvl .gt. 3) then
               call dmout (logfil       , ncv, ncv   ,
     &                     workl(iuptri), ldh, ndigit,
     &              '_neupd: The upper quasi-triangular matrix ')
            end if
         end if
c
         if (reord) then
c
c           %-----------------------------------------------------%
c           | Reorder the computed upper quasi-triangular matrix. |
c           %-----------------------------------------------------%
c
            call dtrsen ('None'       , 'V'          ,
     &                   select       , ncv          ,
     &                   workl(iuptri), ldh          ,
     &                   workl(invsub), ldq          ,
     &                   workl(iheigr), workl(iheigi),
     &                   nconv2       , conds        ,
     &                   sep          , workl(ihbds) ,
     &                   ncv          , iwork        ,
     &                   1            , ierr)
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
                call dvout (logfil, ncv, workl(iheigr), ndigit,
     &           '_neupd: Real part of the eigenvalues of H--reordered')
                call dvout (logfil, ncv, workl(iheigi), ndigit,
     &           '_neupd: Imag part of the eigenvalues of H--reordered')
                if (msglvl .gt. 3) then
                   call dmout (logfil       , ncv, ncv   ,
     &                         workl(iuptri), ldq, ndigit,
     &             '_neupd: Quasi-triangular matrix after re-ordering')
                end if
            end if
c
         end if
c
c        %---------------------------------------%
c        | Copy the last row of the Schur vector |
c        | into workl(ihbds).  This will be used |
c        | to compute the Ritz estimates of      |
c        | converged Ritz values.                |
c        %---------------------------------------%
c
         call dcopy (ncv, workl(invsub+ncv-1), ldq, workl(ihbds), 1)
c
c        %----------------------------------------------------%
c        | Place the computed eigenvalues of H into DR and DI |
c        | if a spectral transformation was not used.         |
c        %----------------------------------------------------%
c
         if (type .eq. 'REGULR') then
            call dcopy (nconv, workl(iheigr), 1, dr, 1)
            call dcopy (nconv, workl(iheigi), 1, di, 1)
         end if
c
c        %----------------------------------------------------------%
c        | Compute the QR factorization of the matrix representing  |
c        | the wanted invariant subspace located in the first NCONV |
c        | columns of workl(invsub,ldq).                            |
c        %----------------------------------------------------------%
c
         call dgeqr2 (ncv, nconv , workl(invsub),
     &               ldq, workev, workev(ncv+1),
     &               ierr)
c
c        %---------------------------------------------------------%
c        | * Postmultiply V by Q using dorm2r .                     |
c        | * Copy the first NCONV columns of VQ into Z.            |
c        | * Postmultiply Z by R.                                  |
c        | The N by NCONV matrix Z is now a matrix representation  |
c        | of the approximate invariant subspace associated with   |
c        | the Ritz values in workl(iheigr) and workl(iheigi)      |
c        | The first NCONV columns of V are now approximate Schur  |
c        | vectors associated with the real upper quasi-triangular |
c        | matrix of order NCONV in workl(iuptri)                  |
c        %---------------------------------------------------------%
c
         call dorm2r ('Right', 'Notranspose', n            ,
     &                ncv   , nconv        , workl(invsub),
     &                ldq   , workev       , v            ,
     &                ldv   , workd(n+1)   , ierr)
         call dlacpy ('All', n, nconv, v, ldv, z, ldz)
c
         do 20 j=1, nconv
c
c           %---------------------------------------------------%
c           | Perform both a column and row scaling if the      |
c           | diagonal element of workl(invsub,ldq) is negative |
c           | I'm lazy and don't take advantage of the upper    |
c           | quasi-triangular form of workl(iuptri,ldq)        |
c           | Note that since Q is orthogonal, R is a diagonal  |
c           | matrix consisting of plus or minus ones           |
c           %---------------------------------------------------%
c
            if (workl(invsub+(j-1)*ldq+j-1) .lt. zero) then
               call dscal (nconv, -one, workl(iuptri+j-1), ldq)
               call dscal (nconv, -one, workl(iuptri+(j-1)*ldq), 1)
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
            call dtrevc ('Right', 'Select'     , select       ,
     &                   ncv    , workl(iuptri), ldq          ,
     &                   vl     , 1            , workl(invsub),
     &                   ldq    , ncv          , outncv       ,
     &                   workev , ierr)
c
            if (ierr .ne. 0) then
                info = -9
                go to 9000
            end if
c
c           %------------------------------------------------%
c           | Scale the returning eigenvectors so that their |
c           | Euclidean norms are all one. LAPACK subroutine |
c           | dtrevc  returns each eigenvector normalized so  |
c           | that the element of largest magnitude has      |
c           | magnitude 1;                                   |
c           %------------------------------------------------%
c
            iconj = 0
            do 40 j=1, nconv
c
               if ( workl(iheigi+j-1) .eq. zero ) then
c
c                 %----------------------%
c                 | real eigenvalue case |
c                 %----------------------%
c
                  temp = dnrm2 ( ncv, workl(invsub+(j-1)*ldq), 1 )
                  call dscal ( ncv, one / temp,
     &                 workl(invsub+(j-1)*ldq), 1 )
c
               else
c
c                 %-------------------------------------------%
c                 | Complex conjugate pair case. Note that    |
c                 | since the real and imaginary part of      |
c                 | the eigenvector are stored in consecutive |
c                 | columns, we further normalize by the      |
c                 | square root of two.                       |
c                 %-------------------------------------------%
c
                  if (iconj .eq. 0) then
                     temp = dlapy2 (dnrm2 (ncv,
     &                                   workl(invsub+(j-1)*ldq),
     &                                   1),
     &                             dnrm2 (ncv,
     &                                   workl(invsub+j*ldq),
     &                                   1))
                     call dscal (ncv, one/temp,
     &                           workl(invsub+(j-1)*ldq), 1 )
                     call dscal (ncv, one/temp,
     &                           workl(invsub+j*ldq), 1 )
                     iconj = 1
                  else
                     iconj = 0
                  end if
c
               end if
c
 40         continue
c
            call dgemv ('T', ncv, nconv, one, workl(invsub),
     &                 ldq, workl(ihbds), 1, zero,  workev, 1)
c
            iconj = 0
            do 45 j=1, nconv
               if (workl(iheigi+j-1) .ne. zero) then
c
c                 %-------------------------------------------%
c                 | Complex conjugate pair case. Note that    |
c                 | since the real and imaginary part of      |
c                 | the eigenvector are stored in consecutive |
c                 %-------------------------------------------%
c
                  if (iconj .eq. 0) then
                     workev(j) = dlapy2 (workev(j), workev(j+1))
                     workev(j+1) = workev(j)
                     iconj = 1
                  else
                     iconj = 0
                  end if
               end if
 45         continue
c
            if (msglvl .gt. 2) then
               call dcopy (ncv, workl(invsub+ncv-1), ldq,
     &                    workl(ihbds), 1)
               call dvout (logfil, ncv, workl(ihbds), ndigit,
     &              '_neupd: Last row of the eigenvector matrix for T')
               if (msglvl .gt. 3) then
                  call dmout (logfil, ncv, ncv, workl(invsub), ldq,
     &                 ndigit, '_neupd: The eigenvector matrix for T')
               end if
            end if
c
c           %---------------------------------------%
c           | Copy Ritz estimates into workl(ihbds) |
c           %---------------------------------------%
c
            call dcopy (nconv, workev, 1, workl(ihbds), 1)
c
c           %---------------------------------------------------------%
c           | Compute the QR factorization of the eigenvector matrix  |
c           | associated with leading portion of T in the first NCONV |
c           | columns of workl(invsub,ldq).                           |
c           %---------------------------------------------------------%
c
            call dgeqr2 (ncv, nconv , workl(invsub),
     &                   ldq, workev, workev(ncv+1),
     &                   ierr)
c
c           %----------------------------------------------%
c           | * Postmultiply Z by Q.                       |
c           | * Postmultiply Z by R.                       |
c           | The N by NCONV matrix Z is now contains the  |
c           | Ritz vectors associated with the Ritz values |
c           | in workl(iheigr) and workl(iheigi).          |
c           %----------------------------------------------%
c
            call dorm2r ('Right', 'Notranspose', n            ,
     &                   ncv  , nconv        , workl(invsub),
     &                   ldq  , workev       , z            ,
     &                   ldz  , workd(n+1)   , ierr)
c
            call dtrmm ('Right'   , 'Upper'       , 'No transpose',
     &                  'Non-unit', n            , nconv         ,
     &                  one       , workl(invsub), ldq           ,
     &                  z         , ldz)
c
         end if
c
      else
c
c        %------------------------------------------------------%
c        | An approximate invariant subspace is not needed.     |
c        | Place the Ritz values computed DNAUPD  into DR and DI |
c        %------------------------------------------------------%
c
         call dcopy (nconv, workl(ritzr), 1, dr, 1)
         call dcopy (nconv, workl(ritzi), 1, di, 1)
         call dcopy (nconv, workl(ritzr), 1, workl(iheigr), 1)
         call dcopy (nconv, workl(ritzi), 1, workl(iheigi), 1)
         call dcopy (nconv, workl(bounds), 1, workl(ihbds), 1)
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
     &      call dscal (ncv, rnorm, workl(ihbds), 1)
c
      else
c
c        %---------------------------------------%
c        |   A spectral transformation was used. |
c        | * Determine the Ritz estimates of the |
c        |   Ritz values in the original system. |
c        %---------------------------------------%
c
         if (type .eq. 'SHIFTI') then
c
            if (rvec)
     &         call dscal (ncv, rnorm, workl(ihbds), 1)
c
            do 50 k=1, ncv
               temp = dlapy2 ( workl(iheigr+k-1),
     &                        workl(iheigi+k-1) )
               workl(ihbds+k-1) = abs( workl(ihbds+k-1) )
     &                          / temp / temp
 50         continue
c
         else if (type .eq. 'REALPT') then
c
            do 60 k=1, ncv
 60         continue
c
         else if (type .eq. 'IMAGPT') then
c
            do 70 k=1, ncv
 70         continue
c
         end if
c
c        %-----------------------------------------------------------%
c        | *  Transform the Ritz values back to the original system. |
c        |    For TYPE = 'SHIFTI' the transformation is              |
c        |             lambda = 1/theta + sigma                      |
c        |    For TYPE = 'REALPT' or 'IMAGPT' the user must from     |
c        |    Rayleigh quotients or a projection. See remark 3 above.|
c        | NOTES:                                                    |
c        | *The Ritz vectors are not affected by the transformation. |
c        %-----------------------------------------------------------%
c
         if (type .eq. 'SHIFTI') then
c
            do 80 k=1, ncv
               temp = dlapy2 ( workl(iheigr+k-1),
     &                        workl(iheigi+k-1) )
               workl(iheigr+k-1) = workl(iheigr+k-1)/temp/temp
     &                           + sigmar
               workl(iheigi+k-1) = -workl(iheigi+k-1)/temp/temp
     &                           + sigmai
 80         continue
c
            call dcopy (nconv, workl(iheigr), 1, dr, 1)
            call dcopy (nconv, workl(iheigi), 1, di, 1)
c
         else if (type .eq. 'REALPT' .or. type .eq. 'IMAGPT') then
c
            call dcopy (nconv, workl(iheigr), 1, dr, 1)
            call dcopy (nconv, workl(iheigi), 1, di, 1)
c
         end if
c
      end if
c
      if (type .eq. 'SHIFTI' .and. msglvl .gt. 1) then
         call dvout (logfil, nconv, dr, ndigit,
     &   '_neupd: Untransformed real part of the Ritz values.')
         call dvout  (logfil, nconv, di, ndigit,
     &   '_neupd: Untransformed imag part of the Ritz values.')
         call dvout (logfil, nconv, workl(ihbds), ndigit,
     &   '_neupd: Ritz estimates of untransformed Ritz values.')
      else if (type .eq. 'REGULR' .and. msglvl .gt. 1) then
         call dvout (logfil, nconv, dr, ndigit,
     &   '_neupd: Real parts of converged Ritz values.')
         call dvout  (logfil, nconv, di, ndigit,
     &   '_neupd: Imag parts of converged Ritz values.')
         call dvout (logfil, nconv, workl(ihbds), ndigit,
     &   '_neupd: Associated Ritz estimates.')
      end if
c
c     %-------------------------------------------------%
c     | Eigenvector Purification step. Formally perform |
c     | one of inverse subspace iteration. Only used    |
c     | for MODE = 2.                                   |
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
c        | where H s = s theta. Remember that when theta  |
c        | has nonzero imaginary part, the corresponding  |
c        | Ritz vector is stored across two columns of Z. |
c        %------------------------------------------------%
c
         iconj = 0
         do 110 j=1, nconv
            if ((workl(iheigi+j-1) .eq. zero) .and.
     &           (workl(iheigr+j-1) .ne. zero)) then
               workev(j) =  workl(invsub+(j-1)*ldq+ncv-1) /
     &                      workl(iheigr+j-1)
            else if (iconj .eq. 0) then
               temp = dlapy2 ( workl(iheigr+j-1), workl(iheigi+j-1) )
               if (temp .ne. zero) then
                  workev(j) = ( workl(invsub+(j-1)*ldq+ncv-1) *
     &                          workl(iheigr+j-1) +
     &                          workl(invsub+j*ldq+ncv-1) *
     &                          workl(iheigi+j-1) ) / temp / temp
                  workev(j+1) = ( workl(invsub+j*ldq+ncv-1) *
     &                            workl(iheigr+j-1) -
     &                            workl(invsub+(j-1)*ldq+ncv-1) *
     &                            workl(iheigi+j-1) ) / temp / temp
               end if
               iconj = 1
            else
               iconj = 0
            end if
 110     continue
c
c        %---------------------------------------%
c        | Perform a rank one update to Z and    |
c        | purify all the Ritz vectors together. |
c        %---------------------------------------%
c
         call dger (n, nconv, one, resid, 1, workev, 1, z, ldz)
c
      end if
c
 9000 continue
c
      return
c
c     %---------------%
c     | End of DNEUPD  |
c     %---------------%
c
      end
