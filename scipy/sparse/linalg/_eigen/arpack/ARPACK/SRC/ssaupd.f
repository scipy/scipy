c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: ssaupd
c
c\Description:
c
c  Reverse communication interface for the Implicitly Restarted Arnoldi
c  Iteration.  For symmetric problems this reduces to a variant of the Lanczos
c  method.  This method has been designed to compute approximations to a
c  few eigenpairs of a linear operator OP that is real and symmetric
c  with respect to a real positive semi-definite symmetric matrix B,
c  i.e.
c
c       B*OP = (OP`)*B.
c
c  Another way to express this condition is
c
c       < x,OPy > = < OPx,y >  where < z,w > = z`Bw  .
c
c  In the standard eigenproblem B is the identity matrix.
c  ( A` denotes transpose of A)
c
c  The computed approximate eigenvalues are called Ritz values and
c  the corresponding approximate eigenvectors are called Ritz vectors.
c
c  ssaupd is usually called iteratively to solve one of the
c  following problems:
c
c  Mode 1:  A*x = lambda*x, A symmetric
c           ===> OP = A  and  B = I.
c
c  Mode 2:  A*x = lambda*M*x, A symmetric, M symmetric positive definite
c           ===> OP = inv[M]*A  and  B = M.
c           ===> (If M can be factored see remark 3 below)
c
c  Mode 3:  K*x = lambda*M*x, K symmetric, M symmetric positive semi-definite
c           ===> OP = (inv[K - sigma*M])*M  and  B = M.
c           ===> Shift-and-Invert mode
c
c  Mode 4:  K*x = lambda*KG*x, K symmetric positive semi-definite,
c           KG symmetric indefinite
c           ===> OP = (inv[K - sigma*KG])*K  and  B = K.
c           ===> Buckling mode
c
c  Mode 5:  A*x = lambda*M*x, A symmetric, M symmetric positive semi-definite
c           ===> OP = inv[A - sigma*M]*[A + sigma*M]  and  B = M.
c           ===> Cayley transformed mode
c
c  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
c        should be accomplished either by a direct method
c        using a sparse matrix factorization and solving
c
c           [A - sigma*M]*w = v  or M*w = v,
c
c        or through an iterative method for solving these
c        systems.  If an iterative method is used, the
c        convergence test must be more stringent than
c        the accuracy requirements for the eigenvalue
c        approximations.
c
c\Usage:
c  call ssaupd
c     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
c       IPNTR, WORKD, WORKL, LWORKL, INFO )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.  IDO must be zero on the first
c          call to ssaupd.  IDO will be set internally to
c          indicate the type of operation to be performed.  Control is
c          then given back to the calling routine which has the
c          responsibility to carry out the requested operation and call
c          ssaupd with the result.  The operand is given in
c          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
c          (If Mode = 2 see remark 5 below)
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    This is for the initialization phase to force the
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    In mode 3,4 and 5, the vector B * X is already
c                    available in WORKD(ipntr(3)).  It does not
c                    need to be recomputed in forming OP * X.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c          IDO =  3: compute the IPARAM(8) shifts where
c                    IPNTR(11) is the pointer into WORKL for
c                    placing the shifts. See remark 6 below.
c          IDO = 99: done
c          -------------------------------------------------------------
c
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.
c          B = 'I' -> standard eigenvalue problem A*x = lambda*x
c          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  WHICH   Character*2.  (INPUT)
c          Specify which of the Ritz values of OP to compute.
c
c          'LA' - compute the NEV largest (algebraic) eigenvalues.
c          'SA' - compute the NEV smallest (algebraic) eigenvalues.
c          'LM' - compute the NEV largest (in magnitude) eigenvalues.
c          'SM' - compute the NEV smallest (in magnitude) eigenvalues.
c          'BE' - compute NEV eigenvalues, half from each end of the
c                 spectrum.  When NEV is odd, compute one more from the
c                 high end than from the low end.
c           (see remark 1 below)
c
c  NEV     Integer.  (INPUT)
c          Number of eigenvalues of OP to be computed. 0 < NEV < N.
c
c  TOL     Real  scalar.  (INPUT)
c          Stopping criterion: the relative accuracy of the Ritz value
c          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
c          If TOL .LE. 0. is passed a default is set:
c          DEFAULT = SLAMCH('EPS')  (machine precision as computed
c                    by the LAPACK auxiliary subroutine SLAMCH).
c
c  RESID   Real  array of length N.  (INPUT/OUTPUT)
c          On INPUT:
c          If INFO .EQ. 0, a random initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          On OUTPUT:
c          RESID contains the final residual vector.
c
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V (less than or equal to N).
c          This will indicate how many Lanczos vectors are generated
c          at each iteration.  After the startup phase in which NEV
c          Lanczos vectors are generated, the algorithm generates
c          NCV-NEV Lanczos vectors at each subsequent update iteration.
c          Most of the cost in generating each Lanczos vector is in the
c          matrix-vector product OP*x. (See remark 4 below).
c
c  V       Real  N by NCV array.  (OUTPUT)
c          The NCV columns of V contain the Lanczos basis vectors.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
c          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
c          The shifts selected at each iteration are used to restart
c          the Arnoldi iteration in an implicit fashion.
c          -------------------------------------------------------------
c          ISHIFT = 0: the shifts are provided by the user via
c                      reverse communication.  The NCV eigenvalues of
c                      the current tridiagonal matrix T are returned in
c                      the part of WORKL array corresponding to RITZ.
c                      See remark 6 below.
c          ISHIFT = 1: exact shifts with respect to the reduced
c                      tridiagonal matrix T.  This is equivalent to
c                      restarting the iteration with a starting vector
c                      that is a linear combination of Ritz vectors
c                      associated with the "wanted" Ritz values.
c          -------------------------------------------------------------
c
c          IPARAM(2) = LEVEC
c          No longer referenced. See remark 2 below.
c
c          IPARAM(3) = MXITER
c          On INPUT:  maximum number of Arnoldi update iterations allowed.
c          On OUTPUT: actual number of Arnoldi update iterations taken.
c
c          IPARAM(4) = NB: blocksize to be used in the recurrence.
c          The code currently works only for NB = 1.
c
c          IPARAM(5) = NCONV: number of "converged" Ritz values.
c          This represents the number of Ritz values that satisfy
c          the convergence criterion.
c
c          IPARAM(6) = IUPD
c          No longer referenced. Implicit restarting is ALWAYS used.
c
c          IPARAM(7) = MODE
c          On INPUT determines what type of eigenproblem is being solved.
c          Must be 1,2,3,4,5; See under \Description of ssaupd for the
c          five modes available.
c
c          IPARAM(8) = NP
c          When ido = 3 and the user provides shifts through reverse
c          communication (IPARAM(1)=0), ssaupd returns NP, the number
c          of shifts the user is to provide. 0 < NP <=NCV-NEV. See Remark
c          6 below.
c
c          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
c          OUTPUT: NUMOP  = total number of OP*x operations,
c                  NUMOPB = total number of B*x operations if BMAT='G',
c                  NUMREO = total number of steps of re-orthogonalization.
c
c  IPNTR   Integer array of length 11.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD and WORKL
c          arrays for matrices/vectors used by the Lanczos iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X in WORKD.
c          IPNTR(2): pointer to the current result vector Y in WORKD.
c          IPNTR(3): pointer to the vector B * X in WORKD when used in
c                    the shift-and-invert mode.
c          IPNTR(4): pointer to the next available location in WORKL
c                    that is untouched by the program.
c          IPNTR(5): pointer to the NCV by 2 tridiagonal matrix T in WORKL.
c          IPNTR(6): pointer to the NCV RITZ values array in WORKL.
c          IPNTR(7): pointer to the Ritz estimates in array WORKL associated
c                    with the Ritz values located in RITZ in WORKL.
c          IPNTR(11): pointer to the NP shifts in WORKL. See Remark 6 below.
c
c          Note: IPNTR(8:10) is only referenced by sseupd. See Remark 2.
c          IPNTR(8): pointer to the NCV RITZ values of the original system.
c          IPNTR(9): pointer to the NCV corresponding error bounds.
c          IPNTR(10): pointer to the NCV by NCV matrix of eigenvectors
c                     of the tridiagonal matrix T. Only referenced by
c                     sseupd if RVEC = .TRUE. See Remarks.
c          -------------------------------------------------------------
c
c  WORKD   Real  work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD
c          as temporary workspace during the iteration. Upon termination
c          WORKD(1:N) contains B*RESID(1:N). If the Ritz vectors are desired
c          subroutine sseupd uses this output.
c          See Data Distribution Note below.
c
c  WORKL   Real  work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  See Data Distribution Note below.
c
c  LWORKL  Integer.  (INPUT)
c          LWORKL must be at least NCV**2 + 8*NCV .
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: Maximum number of iterations taken.
c                All possible eigenvalues of OP has been found. IPARAM(5)
c                returns the number of wanted converged Ritz values.
c          =  2: No longer an informational error. Deprecated starting
c                with release 2 of ARPACK.
c          =  3: No shifts could be applied during a cycle of the
c                Implicitly restarted Arnoldi iteration. One possibility
c                is to increase the size of NCV relative to NEV.
c                See remark 4 below.
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV must be greater than NEV and less than or equal to N.
c          = -4: The maximum number of Arnoldi update iterations allowed
c                must be greater than zero.
c          = -5: WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work array WORKL is not sufficient.
c          = -8: Error return from trid. eigenvalue calculation;
c                Informatinal error from LAPACK routine ssteqr.
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4,5.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: IPARAM(1) must be equal to 0 or 1.
c          = -13: NEV and WHICH = 'BE' are incompatible.
c          = -9999: Could not build an Arnoldi factorization.
c                   IPARAM(5) returns the size of the current Arnoldi
c                   factorization. The user is advised to check that
c                   enough workspace and array storage has been allocated.
c
c
c\Remarks
c  1. The converged Ritz values are always returned in ascending
c     algebraic order.  The computed Ritz values are approximate
c     eigenvalues of OP.  The selection of WHICH should be made
c     with this in mind when Mode = 3,4,5.  After convergence,
c     approximate eigenvalues of the original problem may be obtained
c     with the ARPACK subroutine sseupd.
c
c  2. If the Ritz vectors corresponding to the converged Ritz values
c     are needed, the user must call sseupd immediately following completion
c     of ssaupd. This is new starting with version 2.1 of ARPACK.
c
c  3. If M can be factored into a Cholesky factorization M = LL`
c     then Mode = 2 should not be selected.  Instead one should use
c     Mode = 1 with  OP = inv(L)*A*inv(L`).  Appropriate triangular
c     linear systems should be solved with L and L` rather
c     than computing inverses.  After convergence, an approximate
c     eigenvector z of the original problem is recovered by solving
c     L`z = x  where x is a Ritz vector of OP.
c
c  4. At present there is no a-priori analysis to guide the selection
c     of NCV relative to NEV.  The only formal requirement is that NCV > NEV.
c     However, it is recommended that NCV .ge. 2*NEV.  If many problems of
c     the same type are to be solved, one should experiment with increasing
c     NCV while keeping NEV fixed for a given test problem.  This will
c     usually decrease the required number of OP*x operations but it
c     also increases the work and storage required to maintain the orthogonal
c     basis vectors.   The optimal "cross-over" with respect to CPU time
c     is problem dependent and must be determined empirically.
c
c  5. If IPARAM(7) = 2 then in the Reverse communication interface the user
c     must do the following. When IDO = 1, Y = OP * X is to be computed.
c     When IPARAM(7) = 2 OP = inv(B)*A. After computing A*X the user
c     must overwrite X with A*X. Y is then the solution to the linear set
c     of equations B*Y = A*X.
c
c  6. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the
c     NP = IPARAM(8) shifts in locations:
c     1   WORKL(IPNTR(11))
c     2   WORKL(IPNTR(11)+1)
c                        .
c                        .
c                        .
c     NP  WORKL(IPNTR(11)+NP-1).
c
c     The eigenvalues of the current tridiagonal matrix are located in
c     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1). They are in the
c     order defined by WHICH. The associated Ritz estimates are located in
c     WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , WORKL(IPNTR(8)+NCV-1).
c
c-----------------------------------------------------------------------
c
c\Data Distribution Note:
c
c  Fortran-D syntax:
c  ================
c  REAL       RESID(N), V(LDV,NCV), WORKD(3*N), WORKL(LWORKL)
c  DECOMPOSE  D1(N), D2(N,NCV)
c  ALIGN      RESID(I) with D1(I)
c  ALIGN      V(I,J)   with D2(I,J)
c  ALIGN      WORKD(I) with D1(I)     range (1:N)
c  ALIGN      WORKD(I) with D1(I-N)   range (N+1:2*N)
c  ALIGN      WORKD(I) with D1(I-2*N) range (2*N+1:3*N)
c  DISTRIBUTE D1(BLOCK), D2(BLOCK,:)
c  REPLICATED WORKL(LWORKL)
c
c  Cray MPP syntax:
c  ===============
c  REAL       RESID(N), V(LDV,NCV), WORKD(N,3), WORKL(LWORKL)
c  SHARED     RESID(BLOCK), V(BLOCK,:), WORKD(BLOCK,:)
c  REPLICATED WORKL(LWORKL)
c
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
c  8. R.B. Lehoucq, D.C. Sorensen, "Implementation of Some Spectral
c     Transformations in a k-Step Arnoldi Method". In Preparation.
c
c\Routines called:
c     ssaup2  ARPACK routine that implements the Implicitly Restarted
c             Arnoldi Iteration.
c     sstats  ARPACK routine that initialize timing and other statistics
c             variables.
c     ivout   ARPACK utility routine that prints integers.
c     arscnd  ARPACK utility routine for timing.
c     svout   ARPACK utility routine that prints vectors.
c     slamch  LAPACK routine that determines machine constants.
c
c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Revision history:
c     12/15/93: Version ' 2.4'
c
c\SCCS Information: @(#)
c FILE: saupd.F   SID: 2.8   DATE OF SID: 04/10/01   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine ssaupd
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam,
     &     ipntr, workd, workl, lworkl, info )
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
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Real
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(11), ipntr(11)
      Real
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Real
     &           one, zero
      parameter (one = 1.0E+0 , zero = 0.0E+0 )
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    bounds, ierr, ih, iq, ishift, iupd, iw,
     &           ldh, ldq, msglvl, mxiter, mode, nb,
     &           nev0, next, np, ritz, j
      save       bounds, ierr, ih, iq, ishift, iupd, iw,
     &           ldh, ldq, msglvl, mxiter, mode, nb,
     &           nev0, next, np, ritz
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   ssaup2,  svout, ivout, arscnd, sstats
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Real
     &           slamch
      external   slamch
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
         call sstats
         call arscnd (t0)
         msglvl = msaupd
c
         ierr   = 0
         ishift = iparam(1)
         mxiter = iparam(3)
c         nb     = iparam(4)
         nb     = 1
c
c        %--------------------------------------------%
c        | Revision 2 performs only implicit restart. |
c        %--------------------------------------------%
c
         iupd   = 1
         mode   = iparam(7)
c
c        %----------------%
c        | Error checking |
c        %----------------%
c
         if (n .le. 0) then
            ierr = -1
         else if (nev .le. 0) then
            ierr = -2
         else if (ncv .le. nev .or.  ncv .gt. n) then
            ierr = -3
         end if
c
c        %----------------------------------------------%
c        | NP is the number of additional steps to      |
c        | extend the length NEV Lanczos factorization. |
c        %----------------------------------------------%
c
         np     = ncv - nev
c
         if (mxiter .le. 0)                     ierr = -4
         if (which .ne. 'LM' .and.
     &       which .ne. 'SM' .and.
     &       which .ne. 'LA' .and.
     &       which .ne. 'SA' .and.
     &       which .ne. 'BE')                   ierr = -5
         if (bmat .ne. 'I' .and. bmat .ne. 'G') ierr = -6
c
         if (lworkl .lt. ncv**2 + 8*ncv)        ierr = -7
         if (mode .lt. 1 .or. mode .gt. 5) then
                                                ierr = -10
         else if (mode .eq. 1 .and. bmat .eq. 'G') then
                                                ierr = -11
         else if (ishift .lt. 0 .or. ishift .gt. 1) then
                                                ierr = -12
         else if (nev .eq. 1 .and. which .eq. 'BE') then
                                                ierr = -13
         end if
c
c        %------------%
c        | Error Exit |
c        %------------%
c
         if (ierr .ne. 0) then
            info = ierr
            ido  = 99
            go to 9000
         end if
c
c        %------------------------%
c        | Set default parameters |
c        %------------------------%
c
         if (nb .le. 0)                         nb = 1
         if (tol .le. zero)                     tol = slamch('EpsMach')
c
c        %----------------------------------------------%
c        | NP is the number of additional steps to      |
c        | extend the length NEV Lanczos factorization. |
c        | NEV0 is the local variable designating the   |
c        | size of the invariant subspace desired.      |
c        %----------------------------------------------%
c
         np     = ncv - nev
         nev0   = nev
c
c        %-----------------------------%
c        | Zero out internal workspace |
c        %-----------------------------%
c
         do 10 j = 1, ncv**2 + 8*ncv
            workl(j) = zero
 10      continue
c
c        %-------------------------------------------------------%
c        | Pointer into WORKL for address of H, RITZ, BOUNDS, Q  |
c        | etc... and the remaining workspace.                   |
c        | Also update pointer to be used on output.             |
c        | Memory is laid out as follows:                        |
c        | workl(1:2*ncv) := generated tridiagonal matrix        |
c        | workl(2*ncv+1:2*ncv+ncv) := ritz values               |
c        | workl(3*ncv+1:3*ncv+ncv) := computed error bounds     |
c        | workl(4*ncv+1:4*ncv+ncv*ncv) := rotation matrix Q     |
c        | workl(4*ncv+ncv*ncv+1:7*ncv+ncv*ncv) := workspace     |
c        %-------------------------------------------------------%
c
         ldh    = ncv
         ldq    = ncv
         ih     = 1
         ritz   = ih     + 2*ldh
         bounds = ritz   + ncv
         iq     = bounds + ncv
         iw     = iq     + ncv**2
         next   = iw     + 3*ncv
c
         ipntr(4) = next
         ipntr(5) = ih
         ipntr(6) = ritz
         ipntr(7) = bounds
         ipntr(11) = iw
      end if
c
c     %-------------------------------------------------------%
c     | Carry out the Implicitly restarted Lanczos Iteration. |
c     %-------------------------------------------------------%
c
      call ssaup2
     &   ( ido, bmat, n, which, nev0, np, tol, resid, mode, iupd,
     &     ishift, mxiter, v, ldv, workl(ih), ldh, workl(ritz),
     &     workl(bounds), workl(iq), ldq, workl(iw), ipntr, workd,
     &     info )
c
c     %--------------------------------------------------%
c     | ido .ne. 99 implies use of reverse communication |
c     | to compute operations involving OP or shifts.    |
c     %--------------------------------------------------%
c
      if (ido .eq. 3) iparam(8) = np
      if (ido .ne. 99) go to 9000
c
      iparam(3) = mxiter
      iparam(5) = np
      iparam(9) = nopx
      iparam(10) = nbx
      iparam(11) = nrorth
c
c     %------------------------------------%
c     | Exit if there was an informational |
c     | error within ssaup2.               |
c     %------------------------------------%
c
      if (info .lt. 0) go to 9000
      if (info .eq. 2) info = 3
c
      if (msglvl .gt. 0) then
         call ivout (logfil, 1, [mxiter], ndigit,
     &               '_saupd: number of update iterations taken')
         call ivout (logfil, 1, [np], ndigit,
     &               '_saupd: number of "converged" Ritz values')
         call svout (logfil, np, workl(Ritz), ndigit,
     &               '_saupd: final Ritz values')
         call svout (logfil, np, workl(Bounds), ndigit,
     &               '_saupd: corresponding error bounds')
      end if
c
      call arscnd (t1)
      tsaupd = t1 - t0
c
      if (msglvl .gt. 0) then
c
c        %--------------------------------------------------------%
c        | Version Number & Version Date are defined in version.h |
c        %--------------------------------------------------------%
c
         write (6,1000)
         write (6,1100) mxiter, nopx, nbx, nrorth, nitref, nrstrt,
     &                  tmvopx, tmvbx, tsaupd, tsaup2, tsaitr, titref,
     &                  tgetv0, tseigt, tsgets, tsapps, tsconv
 1000    format (//,
     &      5x, '==========================================',/
     &      5x, '= Symmetric implicit Arnoldi update code =',/
     &      5x, '= Version Number:', ' 2.4' , 19x, ' =',/
     &      5x, '= Version Date:  ', ' 07/31/96' , 14x, ' =',/
     &      5x, '==========================================',/
     &      5x, '= Summary of timing statistics           =',/
     &      5x, '==========================================',//)
 1100    format (
     &      5x, 'Total number update iterations             = ', i5,/
     &      5x, 'Total number of OP*x operations            = ', i5,/
     &      5x, 'Total number of B*x operations             = ', i5,/
     &      5x, 'Total number of reorthogonalization steps  = ', i5,/
     &      5x, 'Total number of iterative refinement steps = ', i5,/
     &      5x, 'Total number of restart steps              = ', i5,/
     &      5x, 'Total time in user OP*x operation          = ', f12.6,/
     &      5x, 'Total time in user B*x operation           = ', f12.6,/
     &      5x, 'Total time in Arnoldi update routine       = ', f12.6,/
     &      5x, 'Total time in saup2 routine                = ', f12.6,/
     &      5x, 'Total time in basic Arnoldi iteration loop = ', f12.6,/
     &      5x, 'Total time in reorthogonalization phase    = ', f12.6,/
     &      5x, 'Total time in (re)start vector generation  = ', f12.6,/
     &      5x, 'Total time in trid eigenvalue subproblem   = ', f12.6,/
     &      5x, 'Total time in getting the shifts           = ', f12.6,/
     &      5x, 'Total time in applying the shifts          = ', f12.6,/
     &      5x, 'Total time in convergence testing          = ', f12.6)
      end if
c
 9000 continue
c
      return
c
c     %---------------%
c     | End of ssaupd |
c     %---------------%
c
      end
