c\BeginDoc
c
c\Name: dseupd 
c
c\Description: 
c
c  This subroutine returns the converged approximations to eigenvalues
c  of A*z = lambda*B*z and (optionally):
c
c      (1) the corresponding approximate eigenvectors,
c
c      (2) an orthonormal (Lanczos) basis for the associated approximate
c          invariant subspace,
c
c      (3) Both.
c
c  There is negligible additional cost to obtain eigenvectors.  An orthonormal
c  (Lanczos) basis is always computed.  There is an additional storage cost 
c  of n*nev if both are requested (in this case a separate array Z must be 
c  supplied).
c
c  These quantities are obtained from the Lanczos factorization computed
c  by DSAUPD  for the linear operator OP prescribed by the MODE selection
c  (see IPARAM(7) in DSAUPD  documentation.)  DSAUPD  must be called before
c  this routine is called. These approximate eigenvalues and vectors are 
c  commonly called Ritz values and Ritz vectors respectively.  They are 
c  referred to as such in the comments that follow.   The computed orthonormal 
c  basis for the invariant subspace corresponding to these Ritz values is 
c  referred to as a Lanczos basis.
c
c  See documentation in the header of the subroutine DSAUPD  for a definition 
c  of OP as well as other terms and the relation of computed Ritz values 
c  and vectors of OP with respect to the given problem  A*z = lambda*B*z.  
c
c  The approximate eigenvalues of the original problem are returned in
c  ascending algebraic order.  The user may elect to call this routine
c  once for each desired Ritz vector and store it peripherally if desired.
c  There is also the option of computing a selected set of these vectors
c  with a single call.
c
c\Usage:
c  call dseupd  
c     ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, BMAT, N, WHICH, NEV, TOL,
c       RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, INFO )
c
c  RVEC    LOGICAL  (INPUT) 
c          Specifies whether Ritz vectors corresponding to the Ritz value 
c          approximations to the eigenproblem A*z = lambda*B*z are computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute Ritz vectors.
c
c  HOWMNY  Character*1  (INPUT) 
c          Specifies how many Ritz vectors are wanted and the form of Z
c          the matrix of Ritz vectors. See remark 1 below.
c          = 'A': compute NEV Ritz vectors;
c          = 'S': compute some of the Ritz vectors, specified
c                 by the logical array SELECT.
c
c  SELECT  Logical array of dimension NCV.  (INPUT/WORKSPACE)
c          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c          computed. To select the Ritz vector corresponding to a
c          Ritz value D(j), SELECT(j) must be set to .TRUE.. 
c          If HOWMNY = 'A' , SELECT is used as a workspace for
c          reordering the Ritz values.
c
c  D       Double precision  array of dimension NEV.  (OUTPUT)
c          On exit, D contains the Ritz value approximations to the
c          eigenvalues of A*z = lambda*B*z. The values are returned
c          in ascending order. If IPARAM(7) = 3,4,5 then D represents
c          the Ritz values of OP computed by dsaupd  transformed to
c          those of the original eigensystem A*z = lambda*B*z. If 
c          IPARAM(7) = 1,2 then the Ritz values of OP are the same 
c          as the those of A*z = lambda*B*z.
c
c  Z       Double precision  N by NEV array if HOWMNY = 'A'.  (OUTPUT)
c          On exit, Z contains the B-orthonormal Ritz vectors of the
c          eigensystem A*z = lambda*B*z corresponding to the Ritz
c          value approximations.
c          If  RVEC = .FALSE. then Z is not referenced.
c          NOTE: The array Z may be set equal to first NEV columns of the 
c          Arnoldi/Lanczos basis array V computed by DSAUPD .
c
c  LDZ     Integer.  (INPUT)
c          The leading dimension of the array Z.  If Ritz vectors are
c          desired, then  LDZ .ge.  max( 1, N ).  In any case,  LDZ .ge. 1.
c
c  SIGMA   Double precision   (INPUT)
c          If IPARAM(7) = 3,4,5 represents the shift. Not referenced if
c          IPARAM(7) = 1 or 2.
c
c
c  **** The remaining arguments MUST be the same as for the   ****
c  **** call to DSAUPD  that was just completed.               ****
c
c  NOTE: The remaining arguments
c
c           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
c           WORKD, WORKL, LWORKL, INFO
c
c         must be passed directly to DSEUPD  following the last call
c         to DSAUPD .  These arguments MUST NOT BE MODIFIED between
c         the the last call to DSAUPD  and the call to DSEUPD .
c
c  Two of these parameters (WORKL, INFO) are also output parameters:
c
c  WORKL   Double precision  work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          WORKL(1:4*ncv) contains information obtained in
c          dsaupd .  They are not changed by dseupd .
c          WORKL(4*ncv+1:ncv*ncv+8*ncv) holds the
c          untransformed Ritz values, the computed error estimates,
c          and the associated eigenvector matrix of H.
c
c          Note: IPNTR(8:10) contains the pointer into WORKL for addresses
c          of the above information computed by dseupd .
c          -------------------------------------------------------------
c          IPNTR(8): pointer to the NCV RITZ values of the original system.
c          IPNTR(9): pointer to the NCV corresponding error bounds.
c          IPNTR(10): pointer to the NCV by NCV matrix of eigenvectors
c                     of the tridiagonal matrix T. Only referenced by
c                     dseupd  if RVEC = .TRUE. See Remarks.
c          -------------------------------------------------------------
c
c  INFO    Integer.  (OUTPUT)
c          Error flag on output.
c          =  0: Normal exit.
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV must be greater than NEV and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from trid. eigenvalue calculation;
c                Information error from LAPACK routine dsteqr .
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4,5.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: NEV and WHICH = 'BE' are incompatible.
c          = -14: DSAUPD  did not find any eigenvalues to sufficient
c                 accuracy.
c          = -15: HOWMNY must be one of 'A' or 'S' if RVEC = .true.
c          = -16: HOWMNY = 'S' not yet implemented
c          = -17: DSEUPD  got a different count of the number of converged
c                 Ritz values than DSAUPD  got.  This indicates the user
c                 probably made an error in passing data from DSAUPD  to
c                 DSEUPD  or that the data was modified before entering 
c                 DSEUPD .
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
c\Remarks
c  1. The converged Ritz values are always returned in increasing 
c     (algebraic) order.
c
c  2. Currently only HOWMNY = 'A' is implemented. It is included at this
c     stage for the user who wants to incorporate it. 
c
c\Routines called:
c     dsesrt   ARPACK routine that sorts an array X, and applies the
c             corresponding permutation to a matrix A.
c     dsortr   dsortr   ARPACK sorting routine.
c     ivout   ARPACK utility routine that prints integers.
c     dvout    ARPACK utility routine that prints vectors.
c     dgeqr2   LAPACK routine that computes the QR factorization of
c             a matrix.
c     dlacpy   LAPACK matrix copy routine.
c     dlamch   LAPACK routine that determines machine constants.
c     dorm2r   LAPACK routine that applies an orthogonal matrix in
c             factored form.
c     dsteqr   LAPACK routine that computes eigenvalues and eigenvectors
c             of a tridiagonal matrix.
c     dger     Level 2 BLAS rank one update to a matrix.
c     dcopy    Level 1 BLAS that copies one vector to another .
c     dnrm2    Level 1 BLAS that computes the norm of a vector.
c     dscal    Level 1 BLAS that scales a vector.
c     dswap    Level 1 BLAS that swaps the contents of two vectors.

c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Chao Yang                    Houston, Texas
c     Dept. of Computational & 
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c 
c\Revision history:
c     12/15/93: Version ' 2.1'
c
c\SCCS Information: @(#) 
c FILE: seupd.F   SID: 2.11   DATE OF SID: 04/10/01   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
      subroutine dseupd (rvec  , howmny, select, d    ,
     &                   z     , ldz   , sigma , bmat ,
     &                   n     , which , nev   , tol  ,
     &                   resid , ncv   , v     , ldv  ,
     &                   iparam, ipntr , workd , workl,
     &                   lworkl, info )
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
     &           sigma, tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(7), ipntr(11)
      logical    select(ncv)
      Double precision 
     &           d(nev)     , resid(n)  , v(ldv,ncv),
     &           z(ldz, nev), workd(2*n), workl(lworkl)
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
      integer    bounds , ierr   , ih    , ihb   , ihd   ,
     &           iq     , iw     , j     , k     , ldh   ,
     &           ldq    , mode   , msglvl, nconv , next  ,
     &           ritz   , irz    , ibd   , np    , ishift,
     &           leftptr, rghtptr, numcnv, jj
      Double precision 
     &           bnorm2 , rnorm, temp, temp1, eps23
      logical    reord
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy  , dger   , dgeqr2 , dlacpy , dorm2r , dscal , 
     &           dsesrt , dsteqr , dswap  , dvout  , ivout , dsortr 
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision 
     &           dnrm2 , dlamch 
      external   dnrm2 , dlamch 
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
c     %------------------------%
c     | Set default parameters |
c     %------------------------%
c
      msglvl = mseupd
      mode = iparam(7)
      nconv = iparam(5)
      info = 0
c
c     %--------------%
c     | Quick return |
c     %--------------%
c
      if (nconv .eq. 0) go to 9000
      ierr = 0
c
      if (nconv .le. 0)                        ierr = -14 
      if (n .le. 0)                            ierr = -1
      if (nev .le. 0)                          ierr = -2
      if (ncv .le. nev .or.  ncv .gt. n)       ierr = -3
      if (which .ne. 'LM' .and.
     &    which .ne. 'SM' .and.
     &    which .ne. 'LA' .and.
     &    which .ne. 'SA' .and.
     &    which .ne. 'BE')                     ierr = -5
      if (bmat .ne. 'I' .and. bmat .ne. 'G')   ierr = -6
      if ( (howmny .ne. 'A' .and.
     &           howmny .ne. 'P' .and.
     &           howmny .ne. 'S') .and. rvec ) 
     &                                         ierr = -15
      if (rvec .and. howmny .eq. 'S')           ierr = -16
c
      if (rvec .and. lworkl .lt. ncv**2+8*ncv) ierr = -7
c     
      if (mode .eq. 1 .or. mode .eq. 2) then
         type = 'REGULR'
      else if (mode .eq. 3 ) then
         type = 'SHIFTI'
      else if (mode .eq. 4 ) then
         type = 'BUCKLE'
      else if (mode .eq. 5 ) then
         type = 'CAYLEY'
      else 
                                               ierr = -10
      end if
      if (mode .eq. 1 .and. bmat .eq. 'G')     ierr = -11
      if (nev .eq. 1 .and. which .eq. 'BE')    ierr = -12
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
c     %-------------------------------------------------------%
c     | Pointer into WORKL for address of H, RITZ, BOUNDS, Q  |
c     | etc... and the remaining workspace.                   |
c     | Also update pointer to be used on output.             |
c     | Memory is laid out as follows:                        |
c     | workl(1:2*ncv) := generated tridiagonal matrix H      |
c     |       The subdiagonal is stored in workl(2:ncv).      |
c     |       The dead spot is workl(1) but upon exiting      |
c     |       dsaupd  stores the B-norm of the last residual   |
c     |       vector in workl(1). We use this !!!             |
c     | workl(2*ncv+1:2*ncv+ncv) := ritz values               |
c     |       The wanted values are in the first NCONV spots. |
c     | workl(3*ncv+1:3*ncv+ncv) := computed Ritz estimates   |
c     |       The wanted values are in the first NCONV spots. |
c     | NOTE: workl(1:4*ncv) is set by dsaupd  and is not      |
c     |       modified by dseupd .                             |
c     %-------------------------------------------------------%
c
c     %-------------------------------------------------------%
c     | The following is used and set by dseupd .              |
c     | workl(4*ncv+1:4*ncv+ncv) := used as workspace during  |
c     |       computation of the eigenvectors of H. Stores    |
c     |       the diagonal of H. Upon EXIT contains the NCV   |
c     |       Ritz values of the original system. The first   |
c     |       NCONV spots have the wanted values. If MODE =   |
c     |       1 or 2 then will equal workl(2*ncv+1:3*ncv).    |
c     | workl(5*ncv+1:5*ncv+ncv) := used as workspace during  |
c     |       computation of the eigenvectors of H. Stores    |
c     |       the subdiagonal of H. Upon EXIT contains the    |
c     |       NCV corresponding Ritz estimates of the         |
c     |       original system. The first NCONV spots have the |
c     |       wanted values. If MODE = 1,2 then will equal    |
c     |       workl(3*ncv+1:4*ncv).                           |
c     | workl(6*ncv+1:6*ncv+ncv*ncv) := orthogonal Q that is  |
c     |       the eigenvector matrix for H as returned by     |
c     |       dsteqr . Not referenced if RVEC = .False.        |
c     |       Ordering follows that of workl(4*ncv+1:5*ncv)   |
c     | workl(6*ncv+ncv*ncv+1:6*ncv+ncv*ncv+2*ncv) :=         |
c     |       Workspace. Needed by dsteqr  and by dseupd .      |
c     | GRAND total of NCV*(NCV+8) locations.                 |
c     %-------------------------------------------------------%
c
c
      ih     = ipntr(5)
      ritz   = ipntr(6)
      bounds = ipntr(7)
      ldh    = ncv
      ldq    = ncv
      ihd    = bounds + ldh
      ihb    = ihd    + ldh
      iq     = ihb    + ldh
      iw     = iq     + ldh*ncv
      next   = iw     + 2*ncv
      ipntr(4)  = next
      ipntr(8)  = ihd
      ipntr(9)  = ihb
      ipntr(10) = iq
c
c     %----------------------------------------%
c     | irz points to the Ritz values computed |
c     |     by _seigt before exiting _saup2.   |
c     | ibd points to the Ritz estimates       |
c     |     computed by _seigt before exiting  |
c     |     _saup2.                            |
c     %----------------------------------------%
c
      irz = ipntr(11)+ncv
      ibd = irz+ncv
c
c
c     %---------------------------------%
c     | Set machine dependent constant. |
c     %---------------------------------%
c
      eps23 = dlamch ('Epsilon-Machine') 
      eps23 = eps23**(2.0D+0  / 3.0D+0 )
c
c     %---------------------------------------%
c     | RNORM is B-norm of the RESID(1:N).    |
c     | BNORM2 is the 2 norm of B*RESID(1:N). |
c     | Upon exit of dsaupd  WORKD(1:N) has    |
c     | B*RESID(1:N).                         |
c     %---------------------------------------%
c
      rnorm = workl(ih)
      if (bmat .eq. 'I') then
         bnorm2 = rnorm
      else if (bmat .eq. 'G') then
         bnorm2 = dnrm2 (n, workd, 1)
      end if
c
      if (msglvl .gt. 2) then
         call dvout (logfil, ncv, workl(irz), ndigit,
     &   '_seupd: Ritz values passed in from _SAUPD.')
         call dvout (logfil, ncv, workl(ibd), ndigit,
     &   '_seupd: Ritz estimates passed in from _SAUPD.')
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
         call dsgets (ishift, which       , nev          ,
     &                np    , workl(irz)  , workl(bounds),
     &                workl)
c
         if (msglvl .gt. 2) then
            call dvout (logfil, ncv, workl(irz), ndigit,
     &      '_seupd: Ritz values after calling _SGETS.')
            call dvout (logfil, ncv, workl(bounds), ndigit,
     &      '_seupd: Ritz value indices after calling _SGETS.')
         end if
c
c        %-----------------------------------------------------%
c        | Record indices of the converged wanted Ritz values  |
c        | Mark the select array for possible reordering       |
c        %-----------------------------------------------------%
c
         numcnv = 0
         do 11 j = 1,ncv
            temp1 = max(eps23, abs(workl(irz+ncv-j)) )
            jj = workl(bounds + ncv - j)
            if (numcnv .lt. nconv .and.
     &          workl(ibd+jj-1) .le. tol*temp1) then
               select(jj) = .true.
               numcnv = numcnv + 1
               if (jj .gt. nev) reord = .true.
            endif
   11    continue
c
c        %-----------------------------------------------------------%
c        | Check the count (numcnv) of converged Ritz values with    |
c        | the number (nconv) reported by _saupd.  If these two      |
c        | are different then there has probably been an error       |
c        | caused by incorrect passing of the _saupd data.           |
c        %-----------------------------------------------------------%
c
         if (msglvl .gt. 2) then
             call ivout(logfil, 1, numcnv, ndigit,
     &            '_seupd: Number of specified eigenvalues')
             call ivout(logfil, 1, nconv, ndigit,
     &            '_seupd: Number of "converged" eigenvalues')
         end if
c
         if (numcnv .ne. nconv) then
            info = -17
            go to 9000
         end if
c
c        %-----------------------------------------------------------%
c        | Call LAPACK routine _steqr to compute the eigenvalues and |
c        | eigenvectors of the final symmetric tridiagonal matrix H. |
c        | Initialize the eigenvector matrix Q to the identity.      |
c        %-----------------------------------------------------------%
c
         call dcopy (ncv-1, workl(ih+1), 1, workl(ihb), 1)
         call dcopy (ncv, workl(ih+ldh), 1, workl(ihd), 1)
c
         call dsteqr ('Identity', ncv, workl(ihd), workl(ihb),
     &                workl(iq) , ldq, workl(iw), ierr)
c
         if (ierr .ne. 0) then
            info = -8
            go to 9000
         end if
c
         if (msglvl .gt. 1) then
            call dcopy (ncv, workl(iq+ncv-1), ldq, workl(iw), 1)
            call dvout (logfil, ncv, workl(ihd), ndigit,
     &          '_seupd: NCV Ritz values of the final H matrix')
            call dvout (logfil, ncv, workl(iw), ndigit,
     &           '_seupd: last row of the eigenvector matrix for H')
         end if
c
         if (reord) then
c
c           %---------------------------------------------%
c           | Reordered the eigenvalues and eigenvectors  |
c           | computed by _steqr so that the "converged"  |
c           | eigenvalues appear in the first NCONV       |
c           | positions of workl(ihd), and the associated |
c           | eigenvectors appear in the first NCONV      |
c           | columns.                                    |
c           %---------------------------------------------%
c
            leftptr = 1
            rghtptr = ncv
c
            if (ncv .eq. 1) go to 30
c
 20         if (select(leftptr)) then
c
c              %-------------------------------------------%
c              | Search, from the left, for the first Ritz |
c              | value that has not converged.             |
c              %-------------------------------------------%
c
               leftptr = leftptr + 1
c
            else if ( .not. select(rghtptr)) then
c
c              %----------------------------------------------%
c              | Search, from the right, the first Ritz value |
c              | that has converged.                          |
c              %----------------------------------------------%
c
               rghtptr = rghtptr - 1
c
            else
c
c              %----------------------------------------------%
c              | Swap the Ritz value on the left that has not |
c              | converged with the Ritz value on the right   |
c              | that has converged.  Swap the associated     |
c              | eigenvector of the tridiagonal matrix H as   |
c              | well.                                        |
c              %----------------------------------------------%
c
               temp = workl(ihd+leftptr-1)
               workl(ihd+leftptr-1) = workl(ihd+rghtptr-1)
               workl(ihd+rghtptr-1) = temp
               call dcopy (ncv, workl(iq+ncv*(leftptr-1)), 1,
     &                    workl(iw), 1)
               call dcopy (ncv, workl(iq+ncv*(rghtptr-1)), 1,
     &                    workl(iq+ncv*(leftptr-1)), 1)
               call dcopy (ncv, workl(iw), 1,
     &                    workl(iq+ncv*(rghtptr-1)), 1)
               leftptr = leftptr + 1
               rghtptr = rghtptr - 1
c
            end if
c
            if (leftptr .lt. rghtptr) go to 20
c
 30      end if
c
         if (msglvl .gt. 2) then
             call dvout  (logfil, ncv, workl(ihd), ndigit,
     &       '_seupd: The eigenvalues of H--reordered')
         end if
c
c        %----------------------------------------%
c        | Load the converged Ritz values into D. |
c        %----------------------------------------%
c
         call dcopy (nconv, workl(ihd), 1, d, 1)
c
      else
c
c        %-----------------------------------------------------%
c        | Ritz vectors not required. Load Ritz values into D. |
c        %-----------------------------------------------------%
c
         call dcopy (nconv, workl(ritz), 1, d, 1)
         call dcopy (ncv, workl(ritz), 1, workl(ihd), 1)
c
      end if
c
c     %------------------------------------------------------------------%
c     | Transform the Ritz values and possibly vectors and corresponding |
c     | Ritz estimates of OP to those of A*x=lambda*B*x. The Ritz values |
c     | (and corresponding data) are returned in ascending order.        |
c     %------------------------------------------------------------------%
c
      if (type .eq. 'REGULR') then
c
c        %---------------------------------------------------------%
c        | Ascending sort of wanted Ritz values, vectors and error |
c        | bounds. Not necessary if only Ritz values are desired.  |
c        %---------------------------------------------------------%
c
         if (rvec) then
            call dsesrt ('LA', rvec , nconv, d, ncv, workl(iq), ldq)
         else
            call dcopy (ncv, workl(bounds), 1, workl(ihb), 1)
         end if
c
      else 
c 
c        %-------------------------------------------------------------%
c        | *  Make a copy of all the Ritz values.                      |
c        | *  Transform the Ritz values back to the original system.   |
c        |    For TYPE = 'SHIFTI' the transformation is                |
c        |             lambda = 1/theta + sigma                        |
c        |    For TYPE = 'BUCKLE' the transformation is                |
c        |             lambda = sigma * theta / ( theta - 1 )          |
c        |    For TYPE = 'CAYLEY' the transformation is                |
c        |             lambda = sigma * (theta + 1) / (theta - 1 )     |
c        |    where the theta are the Ritz values returned by dsaupd .  |
c        | NOTES:                                                      |
c        | *The Ritz vectors are not affected by the transformation.   |
c        |  They are only reordered.                                   |
c        %-------------------------------------------------------------%
c
         call dcopy  (ncv, workl(ihd), 1, workl(iw), 1)
         if (type .eq. 'SHIFTI') then 
            do 40 k=1, ncv
               workl(ihd+k-1) = one / workl(ihd+k-1) + sigma
  40        continue
         else if (type .eq. 'BUCKLE') then
            do 50 k=1, ncv
               workl(ihd+k-1) = sigma * workl(ihd+k-1) / 
     &                          (workl(ihd+k-1) - one)
  50        continue
         else if (type .eq. 'CAYLEY') then
            do 60 k=1, ncv
               workl(ihd+k-1) = sigma * (workl(ihd+k-1) + one) /
     &                          (workl(ihd+k-1) - one)
  60        continue
         end if
c 
c        %-------------------------------------------------------------%
c        | *  Store the wanted NCONV lambda values into D.             |
c        | *  Sort the NCONV wanted lambda in WORKL(IHD:IHD+NCONV-1)   |
c        |    into ascending order and apply sort to the NCONV theta   |
c        |    values in the transformed system. We will need this to   |
c        |    compute Ritz estimates in the original system.           |
c        | *  Finally sort the lambda`s into ascending order and apply |
c        |    to Ritz vectors if wanted. Else just sort lambda`s into  |
c        |    ascending order.                                         |
c        | NOTES:                                                      |
c        | *workl(iw:iw+ncv-1) contain the theta ordered so that they  |
c        |  match the ordering of the lambda. We`ll use them again for |
c        |  Ritz vector purification.                                  |
c        %-------------------------------------------------------------%
c
         call dcopy (nconv, workl(ihd), 1, d, 1)
         call dsortr ('LA', .true., nconv, workl(ihd), workl(iw))
         if (rvec) then
            call dsesrt ('LA', rvec , nconv, d, ncv, workl(iq), ldq)
         else
            call dcopy (ncv, workl(bounds), 1, workl(ihb), 1)
            call dscal (ncv, bnorm2/rnorm, workl(ihb), 1)
            call dsortr ('LA', .true., nconv, d, workl(ihb))
         end if
c
      end if 
c 
c     %------------------------------------------------%
c     | Compute the Ritz vectors. Transform the wanted |
c     | eigenvectors of the symmetric tridiagonal H by |
c     | the Lanczos basis matrix V.                    |
c     %------------------------------------------------%
c
      if (rvec .and. howmny .eq. 'A') then
c    
c        %----------------------------------------------------------%
c        | Compute the QR factorization of the matrix representing  |
c        | the wanted invariant subspace located in the first NCONV |
c        | columns of workl(iq,ldq).                                |
c        %----------------------------------------------------------%
c     
         call dgeqr2 (ncv, nconv        , workl(iq) ,
     &                ldq, workl(iw+ncv), workl(ihb),
     &                ierr)
c
c        %--------------------------------------------------------%
c        | * Postmultiply V by Q.                                 |   
c        | * Copy the first NCONV columns of VQ into Z.           |
c        | The N by NCONV matrix Z is now a matrix representation |
c        | of the approximate invariant subspace associated with  |
c        | the Ritz values in workl(ihd).                         |
c        %--------------------------------------------------------%
c     
         call dorm2r ('Right', 'Notranspose', n        ,
     &                ncv    , nconv        , workl(iq),
     &                ldq    , workl(iw+ncv), v        ,
     &                ldv    , workd(n+1)   , ierr)
         call dlacpy ('All', n, nconv, v, ldv, z, ldz)
c
c        %-----------------------------------------------------%
c        | In order to compute the Ritz estimates for the Ritz |
c        | values in both systems, need the last row of the    |
c        | eigenvector matrix. Remember, it`s in factored form |
c        %-----------------------------------------------------%
c
         do 65 j = 1, ncv-1
            workl(ihb+j-1) = zero 
  65     continue
         workl(ihb+ncv-1) = one
         call dorm2r ('Left', 'Transpose'  , ncv       ,
     &                1     , nconv        , workl(iq) ,
     &                ldq   , workl(iw+ncv), workl(ihb),
     &                ncv   , temp         , ierr)
c
      else if (rvec .and. howmny .eq. 'S') then
c
c     Not yet implemented. See remark 2 above.
c
      end if
c
      if (type .eq. 'REGULR' .and. rvec) then
c
            do 70 j=1, ncv
               workl(ihb+j-1) = rnorm * abs( workl(ihb+j-1) )
 70         continue
c
      else if (type .ne. 'REGULR' .and. rvec) then
c
c        %-------------------------------------------------%
c        | *  Determine Ritz estimates of the theta.       |
c        |    If RVEC = .true. then compute Ritz estimates |
c        |               of the theta.                     |
c        |    If RVEC = .false. then copy Ritz estimates   |
c        |              as computed by dsaupd .             |
c        | *  Determine Ritz estimates of the lambda.      |
c        %-------------------------------------------------%
c
         call dscal  (ncv, bnorm2, workl(ihb), 1)
         if (type .eq. 'SHIFTI') then 
c
            do 80 k=1, ncv
               workl(ihb+k-1) = abs( workl(ihb+k-1) ) 
     &                        / workl(iw+k-1)**2
 80         continue
c
         else if (type .eq. 'BUCKLE') then
c
            do 90 k=1, ncv
               workl(ihb+k-1) = sigma * abs( workl(ihb+k-1) )
     &                        / (workl(iw+k-1)-one )**2
 90         continue
c
         else if (type .eq. 'CAYLEY') then
c
            do 100 k=1, ncv
               workl(ihb+k-1) = abs( workl(ihb+k-1)
     &                        / workl(iw+k-1)*(workl(iw+k-1)-one) )
 100        continue
c
         end if
c
      end if
c
      if (type .ne. 'REGULR' .and. msglvl .gt. 1) then
         call dvout (logfil, nconv, d, ndigit,
     &          '_seupd: Untransformed converged Ritz values')
         call dvout (logfil, nconv, workl(ihb), ndigit, 
     &     '_seupd: Ritz estimates of the untransformed Ritz values')
      else if (msglvl .gt. 1) then
         call dvout (logfil, nconv, d, ndigit,
     &          '_seupd: Converged Ritz values')
         call dvout (logfil, nconv, workl(ihb), ndigit, 
     &     '_seupd: Associated Ritz estimates')
      end if
c 
c     %-------------------------------------------------%
c     | Ritz vector purification step. Formally perform |
c     | one of inverse subspace iteration. Only used    |
c     | for MODE = 3,4,5. See reference 7               |
c     %-------------------------------------------------%
c
      if (rvec .and. (type .eq. 'SHIFTI' .or. type .eq. 'CAYLEY')) then
c
         do 110 k=0, nconv-1
            workl(iw+k) = workl(iq+k*ldq+ncv-1)
     &                  / workl(iw+k)
 110     continue
c
      else if (rvec .and. type .eq. 'BUCKLE') then
c
         do 120 k=0, nconv-1
            workl(iw+k) = workl(iq+k*ldq+ncv-1)
     &                  / (workl(iw+k)-one)
 120     continue
c
      end if 
c
      if (type .ne. 'REGULR')
     &   call dger  (n, nconv, one, resid, 1, workl(iw), 1, z, ldz)
c
 9000 continue
c
      return
c
c     %---------------%
c     | End of dseupd |
c     %---------------%
c
      end
