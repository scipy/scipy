c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: sseigt
c
c\Description: 
c  Compute the eigenvalues of the current symmetric tridiagonal matrix
c  and the corresponding error bounds given the current residual norm.
c
c\Usage:
c  call sseigt
c     ( RNORM, N, H, LDH, EIG, BOUNDS, WORKL, IERR )
c
c\Arguments
c  RNORM   Real scalar.  (INPUT)
c          RNORM contains the residual norm corresponding to the current
c          symmetric tridiagonal matrix H.
c
c  N       Integer.  (INPUT)
c          Size of the symmetric tridiagonal matrix H.
c
c  H       Real N by 2 array.  (INPUT)
c          H contains the symmetric tridiagonal matrix with the 
c          subdiagonal in the first column starting at H(2,1) and the 
c          main diagonal in second column.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling 
c          program.
c
c  EIG     Real array of length N.  (OUTPUT)
c          On output, EIG contains the N eigenvalues of H possibly 
c          unsorted.  The BOUNDS arrays are returned in the
c          same sorted order as EIG.
c
c  BOUNDS  Real array of length N.  (OUTPUT)
c          On output, BOUNDS contains the error estimates corresponding
c          to the eigenvalues EIG.  This is equal to RNORM times the
c          last components of the eigenvectors corresponding to the
c          eigenvalues in EIG.
c
c  WORKL   Real work array of length 3*N.  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.
c
c  IERR    Integer.  (OUTPUT)
c          Error exit flag from sstqrb.
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
c     sstqrb  ARPACK routine that computes the eigenvalues and the
c             last components of the eigenvectors of a symmetric
c             and tridiagonal matrix.
c     second  ARPACK utility routine for timing.
c     svout   ARPACK utility routine that prints vectors.
c     scopy   Level 1 BLAS that copies one vector to another.
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
c     xx/xx/92: Version ' 2.4'
c
c\SCCS Information: @(#) 
c FILE: seigt.F   SID: 2.4   DATE OF SID: 8/27/96   RELEASE: 2
c
c\Remarks
c     None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine sseigt 
     &   ( rnorm, n, h, ldh, eig, bounds, workl, ierr )
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
      integer    ierr, ldh, n
      Real
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Real
     &           eig(n), bounds(n), h(ldh,2), workl(3*n)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Real
     &           zero
      parameter (zero = 0.0E+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i, k, msglvl
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   scopy, sstqrb, svout, second
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------% 
c
      call second (t0)
      msglvl = mseigt
c
      if (msglvl .gt. 0) then
         call svout (logfil, n, h(1,2), ndigit,
     &              '_seigt: main diagonal of matrix H')
         if (n .gt. 1) then
         call svout (logfil, n-1, h(2,1), ndigit,
     &              '_seigt: sub diagonal of matrix H')
         end if
      end if
c
      call scopy  (n, h(1,2), 1, eig, 1)
      call scopy  (n-1, h(2,1), 1, workl, 1)
      call sstqrb (n, eig, workl, bounds, workl(n+1), ierr)
      if (ierr .ne. 0) go to 9000
      if (msglvl .gt. 1) then
         call svout (logfil, n, bounds, ndigit,
     &              '_seigt: last row of the eigenvector matrix for H')
      end if
c
c     %-----------------------------------------------%
c     | Finally determine the error bounds associated |
c     | with the n Ritz values of H.                  |
c     %-----------------------------------------------%
c
      do 30 k = 1, n
         bounds(k) = rnorm*abs(bounds(k))
   30 continue
c 
      call second (t1)
      tseigt = tseigt + (t1 - t0)
c
 9000 continue
      return
c
c     %---------------%
c     | End of sseigt |
c     %---------------%
c
      end
