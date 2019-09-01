c\BeginDoc
c
c\Name: snaupd
c
c\Description: 
c  Reverse communication interface for the Implicitly Restarted Arnoldi
c  iteration. This subroutine computes approximations to a few eigenpairs 
c  of a linear operator "OP" with respect to a semi-inner product defined by 
c  a symmetric positive semi-definite real matrix B. B may be the identity 
c  matrix. NOTE: If the linear operator "OP" is real and symmetric 
c  with respect to the real positive semi-definite symmetric matrix B, 
c  i.e. B*OP = (OP`)*B, then subroutine ssaupd should be used instead.
c
c  The computed approximate eigenvalues are called Ritz values and
c  the corresponding approximate eigenvectors are called Ritz vectors.
c
c  snaupd is usually called iteratively to solve one of the 
c  following problems:
c
c  Mode 1:  A*x = lambda*x.
c           ===> OP = A  and  B = I.
c
c  Mode 2:  A*x = lambda*M*x, M symmetric positive definite
c           ===> OP = inv[M]*A  and  B = M.
c           ===> (If M can be factored see remark 3 below)
c
c  Mode 3:  A*x = lambda*M*x, M symmetric semi-definite
c           ===> OP = Real_Part{ inv[A - sigma*M]*M }  and  B = M. 
c           ===> shift-and-invert mode (in real arithmetic)
c           If OP*x = amu*x, then 
c           amu = 1/2 * [ 1/(lambda-sigma) + 1/(lambda-conjg(sigma)) ].
c           Note: If sigma is real, i.e. imaginary part of sigma is zero;
c                 Real_Part{ inv[A - sigma*M]*M } == inv[A - sigma*M]*M 
c                 amu == 1/(lambda-sigma). 
c  
c  Mode 4:  A*x = lambda*M*x, M symmetric semi-definite
c           ===> OP = Imaginary_Part{ inv[A - sigma*M]*M }  and  B = M. 
c           ===> shift-and-invert mode (in real arithmetic)
c           If OP*x = amu*x, then 
c           amu = 1/2i * [ 1/(lambda-sigma) - 1/(lambda-conjg(sigma)) ].
c
c  Both mode 3 and 4 give the same enhancement to eigenvalues close to
c  the (complex) shift sigma.  However, as lambda goes to infinity,
c  the operator OP in mode 4 dampens the eigenvalues more strongly than
c  does OP defined in mode 3.
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
c  call snaupd
c     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
c       IPNTR, WORKD, WORKL, LWORKL, INFO )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.  IDO must be zero on the first 
c          call to snaupd.  IDO will be set internally to
c          indicate the type of operation to be performed.  Control is
c          then given back to the calling routine which has the
c          responsibility to carry out the requested operation and call
c          snaupd with the result.  The operand is given in
c          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    This is for the initialization phase to force the
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    In mode 3 and 4, the vector B * X is already
c                    available in WORKD(ipntr(3)).  It does not
c                    need to be recomputed in forming OP * X.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c          IDO =  3: compute the IPARAM(8) real and imaginary parts 
c                    of the shifts where INPTR(14) is the pointer
c                    into WORKL for placing the shifts. See Remark
c                    5 below.
c          IDO = 99: done
c          -------------------------------------------------------------
c             
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.
c          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
c          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  WHICH   Character*2.  (INPUT)
c          'LM' -> want the NEV eigenvalues of largest magnitude.
c          'SM' -> want the NEV eigenvalues of smallest magnitude.
c          'LR' -> want the NEV eigenvalues of largest real part.
c          'SR' -> want the NEV eigenvalues of smallest real part.
c          'LI' -> want the NEV eigenvalues of largest imaginary part.
c          'SI' -> want the NEV eigenvalues of smallest imaginary part.
c
c  NEV     Integer.  (INPUT)
c          Number of eigenvalues of OP to be computed. 0 < NEV < N-1.
c
c  TOL     Real  scalar.  (INPUT)
c          Stopping criterion: the relative accuracy of the Ritz value 
c          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I))
c          where ABS(RITZ(I)) is the magnitude when RITZ(I) is complex.
c          DEFAULT = slamch('EPS')  (machine precision as computed
c                    by the LAPACK auxiliary subroutine slamch).
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
c          Number of columns of the matrix V. NCV must satisfy the two
c          inequalities 2 <= NCV-NEV and NCV <= N.
c          This will indicate how many Arnoldi vectors are generated 
c          at each iteration.  After the startup phase in which NEV 
c          Arnoldi vectors are generated, the algorithm generates 
c          approximately NCV-NEV Arnoldi vectors at each subsequent update 
c          iteration. Most of the cost in generating each Arnoldi vector is 
c          in the matrix-vector operation OP*x. 
c          NOTE: 2 <= NCV-NEV in order that complex conjugate pairs of Ritz 
c          values are kept together. (See remark 4 below)
c
c  V       Real  array N by NCV.  (OUTPUT)
c          Contains the final set of Arnoldi basis vectors. 
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling program.
c
c  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
c          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
c          The shifts selected at each iteration are used to restart
c          the Arnoldi iteration in an implicit fashion.
c          -------------------------------------------------------------
c          ISHIFT = 0: the shifts are provided by the user via
c                      reverse communication.  The real and imaginary
c                      parts of the NCV eigenvalues of the Hessenberg
c                      matrix H are returned in the part of the WORKL 
c                      array corresponding to RITZR and RITZI. See remark 
c                      5 below.
c          ISHIFT = 1: exact shifts with respect to the current
c                      Hessenberg matrix H.  This is equivalent to 
c                      restarting the iteration with a starting vector
c                      that is a linear combination of approximate Schur
c                      vectors associated with the "wanted" Ritz values.
c          -------------------------------------------------------------
c
c          IPARAM(2) = No longer referenced.
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
c          Must be 1,2,3,4; See under \Description of snaupd for the 
c          four modes available.
c
c          IPARAM(8) = NP
c          When ido = 3 and the user provides shifts through reverse
c          communication (IPARAM(1)=0), snaupd returns NP, the number
c          of shifts the user is to provide. 0 < NP <=NCV-NEV. See Remark
c          5 below.
c
c          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
c          OUTPUT: NUMOP  = total number of OP*x operations,
c                  NUMOPB = total number of B*x operations if BMAT='G',
c                  NUMREO = total number of steps of re-orthogonalization.        
c
c  IPNTR   Integer array of length 14.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD and WORKL
c          arrays for matrices/vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X in WORKD.
c          IPNTR(2): pointer to the current result vector Y in WORKD.
c          IPNTR(3): pointer to the vector B * X in WORKD when used in 
c                    the shift-and-invert mode.
c          IPNTR(4): pointer to the next available location in WORKL
c                    that is untouched by the program.
c          IPNTR(5): pointer to the NCV by NCV upper Hessenberg matrix
c                    H in WORKL.
c          IPNTR(6): pointer to the real part of the ritz value array 
c                    RITZR in WORKL.
c          IPNTR(7): pointer to the imaginary part of the ritz value array
c                    RITZI in WORKL.
c          IPNTR(8): pointer to the Ritz estimates in array WORKL associated
c                    with the Ritz values located in RITZR and RITZI in WORKL.
c
c          IPNTR(14): pointer to the NP shifts in WORKL. See Remark 5 below.
c
c          Note: IPNTR(9:13) is only referenced by sneupd. See Remark 2 below.
c
c          IPNTR(9):  pointer to the real part of the NCV RITZ values of the 
c                     original system.
c          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of 
c                     the original system.
c          IPNTR(11): pointer to the NCV corresponding error bounds.
c          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
c                     Schur matrix for H.
c          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
c                     of the upper Hessenberg matrix H. Only referenced by
c                     sneupd if RVEC = .TRUE. See Remark 2 below.
c          -------------------------------------------------------------
c          
c  WORKD   Real  work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD 
c          as temporary workspace during the iteration. Upon termination
c          WORKD(1:N) contains B*RESID(1:N). If an invariant subspace
c          associated with the converged Ritz values is desired, see remark
c          2 below, subroutine sneupd uses this output.
c          See Data Distribution Note below.  
c
c  WORKL   Real  work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  See Data Distribution Note below.
c
c  LWORKL  Integer.  (INPUT)
c          LWORKL must be at least 3*NCV**2 + 6*NCV.
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
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -4: The maximum number of Arnoldi update iteration 
c                must be greater than zero.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work array is not sufficient.
c          = -8: Error return from LAPACK eigenvalue calculation;
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
c          = -12: IPARAM(1) must be equal to 0 or 1.
c          = -9999: Could not build an Arnoldi factorization.
c                   IPARAM(5) returns the size of the current Arnoldi
c                   factorization.
c
c\Remarks
c  1. The computed Ritz values are approximate eigenvalues of OP. The
c     selection of WHICH should be made with this in mind when
c     Mode = 3 and 4.  After convergence, approximate eigenvalues of the
c     original problem may be obtained with the ARPACK subroutine sneupd.
c
c  2. If a basis for the invariant subspace corresponding to the converged Ritz 
c     values is needed, the user must call sneupd immediately following 
c     completion of snaupd. This is new starting with release 2 of ARPACK.
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
c     of NCV relative to NEV.  The only formal requrement is that NCV > NEV + 2.
c     However, it is recommended that NCV .ge. 2*NEV+1.  If many problems of
c     the same type are to be solved, one should experiment with increasing
c     NCV while keeping NEV fixed for a given test problem.  This will 
c     usually decrease the required number of OP*x operations but it
c     also increases the work and storage required to maintain the orthogonal
c     basis vectors.  The optimal "cross-over" with respect to CPU time
c     is problem dependent and must be determined empirically. 
c     See Chapter 8 of Reference 2 for further information.
c
c  5. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the 
c     NP = IPARAM(8) real and imaginary parts of the shifts in locations 
c         real part                  imaginary part
c         -----------------------    --------------
c     1   WORKL(IPNTR(14))           WORKL(IPNTR(14)+NP)
c     2   WORKL(IPNTR(14)+1)         WORKL(IPNTR(14)+NP+1)
c                        .                          .
c                        .                          .
c                        .                          .
c     NP  WORKL(IPNTR(14)+NP-1)      WORKL(IPNTR(14)+2*NP-1).
c
c     Only complex conjugate pairs of shifts may be applied and the pairs 
c     must be placed in consecutive locations. The real part of the 
c     eigenvalues of the current upper Hessenberg matrix are located in 
c     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1) and the imaginary part 
c     in WORKL(IPNTR(7)) through WORKL(IPNTR(7)+NCV-1). They are ordered
c     according to the order defined by WHICH. The complex conjugate
c     pairs are kept together and the associated Ritz estimates are located in
c     WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , WORKL(IPNTR(8)+NCV-1).
c
c-----------------------------------------------------------------------
c
c\Data Distribution Note: 
c
c  Fortran-D syntax:
c  ================
c  Real  resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
c  decompose  d1(n), d2(n,ncv)
c  align      resid(i) with d1(i)
c  align      v(i,j)   with d2(i,j)
c  align      workd(i) with d1(i)     range (1:n)
c  align      workd(i) with d1(i-n)   range (n+1:2*n)
c  align      workd(i) with d1(i-2*n) range (2*n+1:3*n)
c  distribute d1(block), d2(block,:)
c  replicated workl(lworkl)
c
c  Cray MPP syntax:
c  ===============
c  Real   resid(n), v(ldv,ncv), workd(n,3), workl(lworkl)
c  shared     resid(block), v(block,:), workd(block,:)
c  replicated workl(lworkl)
c  
c  CM2/CM5 syntax:
c  ==============
c  
c-----------------------------------------------------------------------
c
c     include   'ex-nonsym.doc'
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
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
c     snaup2  ARPACK routine that implements the Implicitly Restarted
c             Arnoldi Iteration.
c     ivout   ARPACK utility routine that prints integers.
c     arscnd  ARPACK utility routine for timing.
c     svout   ARPACK utility routine that prints vectors.
c     slamch  LAPACK routine that determines machine constants.
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
c     12/16/93: Version '1.1'
c
c\SCCS Information: @(#) 
c FILE: naupd.F   SID: 2.8   DATE OF SID: 04/10/01   RELEASE: 2
c
c\Remarks
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine snaupd
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
      integer    iparam(11), ipntr(14)
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
     &           ldh, ldq, levec, mode, msglvl, mxiter, nb,
     &           nev0, next, np, ritzi, ritzr, j
      save       bounds, ih, iq, ishift, iupd, iw, ldh, ldq,
     &           levec, mode, msglvl, mxiter, nb, nev0, next,
     &           np, ritzi, ritzr
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   snaup2, svout, ivout, arscnd, sstatn
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
         call sstatn
         call arscnd (t0)
         msglvl = mnaupd
c
c        %----------------%
c        | Error checking |
c        %----------------%
c
         ierr   = 0
         ishift = iparam(1)
c         levec  = iparam(2)
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
         if (n .le. 0) then
             ierr = -1
         else if (nev .le. 0) then
             ierr = -2
         else if (ncv .le. nev+1 .or.  ncv .gt. n) then
             ierr = -3
         else if (mxiter .le. 0) then
             ierr = -4
         else if (which .ne. 'LM' .and.
     &       which .ne. 'SM' .and.
     &       which .ne. 'LR' .and.
     &       which .ne. 'SR' .and.
     &       which .ne. 'LI' .and.
     &       which .ne. 'SI') then
            ierr = -5
         else if (bmat .ne. 'I' .and. bmat .ne. 'G') then
            ierr = -6
         else if (lworkl .lt. 3*ncv**2 + 6*ncv) then
            ierr = -7
         else if (mode .lt. 1 .or. mode .gt. 4) then
                                                ierr = -10
         else if (mode .eq. 1 .and. bmat .eq. 'G') then
                                                ierr = -11
         else if (ishift .lt. 0 .or. ishift .gt. 1) then
                                                ierr = -12
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
         if (nb .le. 0)				nb = 1
         if (tol .le. zero)			tol = slamch('EpsMach')
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
         do 10 j = 1, 3*ncv**2 + 6*ncv
            workl(j) = zero
  10     continue
c 
c        %-------------------------------------------------------------%
c        | Pointer into WORKL for address of H, RITZ, BOUNDS, Q        |
c        | etc... and the remaining workspace.                         |
c        | Also update pointer to be used on output.                   |
c        | Memory is laid out as follows:                              |
c        | workl(1:ncv*ncv) := generated Hessenberg matrix             |
c        | workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary        |
c        |                                   parts of ritz values      |
c        | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds        |
c        | workl(ncv*ncv+3*ncv+1:2*ncv*ncv+3*ncv) := rotation matrix Q |
c        | workl(2*ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) := workspace       |
c        | The final workspace is needed by subroutine sneigh called   |
c        | by snaup2. Subroutine sneigh calls LAPACK routines for      |
c        | calculating eigenvalues and the last row of the eigenvector |
c        | matrix.                                                     |
c        %-------------------------------------------------------------%
c
         ldh    = ncv
         ldq    = ncv
         ih     = 1
         ritzr  = ih     + ldh*ncv
         ritzi  = ritzr  + ncv
         bounds = ritzi  + ncv
         iq     = bounds + ncv
         iw     = iq     + ldq*ncv
         next   = iw     + ncv**2 + 3*ncv
c
         ipntr(4) = next
         ipntr(5) = ih
         ipntr(6) = ritzr
         ipntr(7) = ritzi
         ipntr(8) = bounds
         ipntr(14) = iw 
c
      end if
c
c     %-------------------------------------------------------%
c     | Carry out the Implicitly restarted Arnoldi Iteration. |
c     %-------------------------------------------------------%
c
      call snaup2 
     &   ( ido, bmat, n, which, nev0, np, tol, resid, mode, iupd,
     &     ishift, mxiter, v, ldv, workl(ih), ldh, workl(ritzr), 
     &     workl(ritzi), workl(bounds), workl(iq), ldq, workl(iw), 
     &     ipntr, workd, info )
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
c     | error within snaup2.               |
c     %------------------------------------%
c
      if (info .lt. 0) go to 9000
      if (info .eq. 2) info = 3
c
      if (msglvl .gt. 0) then
         call ivout (logfil, 1, mxiter, ndigit,
     &               '_naupd: Number of update iterations taken')
         call ivout (logfil, 1, np, ndigit,
     &               '_naupd: Number of wanted "converged" Ritz values')
         call svout (logfil, np, workl(ritzr), ndigit, 
     &               '_naupd: Real part of the final Ritz values')
         call svout (logfil, np, workl(ritzi), ndigit, 
     &               '_naupd: Imaginary part of the final Ritz values')
         call svout (logfil, np, workl(bounds), ndigit, 
     &               '_naupd: Associated Ritz estimates')
      end if
c
      call arscnd (t1)
      tnaupd = t1 - t0
c
      if (msglvl .gt. 0) then
c
c        %--------------------------------------------------------%
c        | Version Number & Version Date are defined in version.h |
c        %--------------------------------------------------------%
c
         write (6,1000)
         write (6,1100) mxiter, nopx, nbx, nrorth, nitref, nrstrt,
     &                  tmvopx, tmvbx, tnaupd, tnaup2, tnaitr, titref,
     &                  tgetv0, tneigh, tngets, tnapps, tnconv, trvec
 1000    format (//,
     &      5x, '=============================================',/
     &      5x, '= Nonsymmetric implicit Arnoldi update code =',/
     &      5x, '= Version Number: ', ' 2.4' , 21x, ' =',/
     &      5x, '= Version Date:   ', ' 07/31/96' , 16x,   ' =',/
     &      5x, '=============================================',/
     &      5x, '= Summary of timing statistics              =',/
     &      5x, '=============================================',//)
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
     &      5x, 'Total time in naup2 routine                = ', f12.6,/
     &      5x, 'Total time in basic Arnoldi iteration loop = ', f12.6,/
     &      5x, 'Total time in reorthogonalization phase    = ', f12.6,/
     &      5x, 'Total time in (re)start vector generation  = ', f12.6,/
     &      5x, 'Total time in Hessenberg eig. subproblem   = ', f12.6,/
     &      5x, 'Total time in getting the shifts           = ', f12.6,/
     &      5x, 'Total time in applying the shifts          = ', f12.6,/
     &      5x, 'Total time in convergence testing          = ', f12.6,/
     &      5x, 'Total time in computing final Ritz vectors = ', f12.6/)
      end if
c
 9000 continue
c
      return
c
c     %---------------%
c     | End of snaupd |
c     %---------------%
c
      end
