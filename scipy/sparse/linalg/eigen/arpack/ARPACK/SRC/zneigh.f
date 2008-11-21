c\BeginDoc
c
c\Name: zneigh
c
c\Description:
c  Compute the eigenvalues of the current upper Hessenberg matrix
c  and the corresponding Ritz estimates given the current residual norm.
c
c\Usage:
c  call zneigh
c     ( RNORM, N, H, LDH, RITZ, BOUNDS, Q, LDQ, WORKL, RWORK, IERR )
c
c\Arguments
c  RNORM   Double precision scalar.  (INPUT)
c          Residual norm corresponding to the current upper Hessenberg 
c          matrix H.
c
c  N       Integer.  (INPUT)
c          Size of the matrix H.
c
c  H       Complex*16 N by N array.  (INPUT)
c          H contains the current upper Hessenberg matrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RITZ    Complex*16 array of length N.  (OUTPUT)
c          On output, RITZ(1:N) contains the eigenvalues of H.
c
c  BOUNDS  Complex*16 array of length N.  (OUTPUT)
c          On output, BOUNDS contains the Ritz estimates associated with
c          the eigenvalues held in RITZ.  This is equal to RNORM 
c          times the last components of the eigenvectors corresponding 
c          to the eigenvalues in RITZ.
c
c  Q       Complex*16 N by N array.  (WORKSPACE)
c          Workspace needed to store the eigenvectors of H.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Complex*16 work array of length N**2 + 3*N.  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  This is needed to keep the full Schur form
c          of H and also in the calculation of the eigenvectors of H.
c
c  RWORK   Double precision  work array of length N (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end. 
c
c  IERR    Integer.  (OUTPUT)
c          Error exit flag from zlahqr or ztrevc.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  Complex*16
c
c\Routines called:
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     zmout   ARPACK utility routine that prints matrices
c     zvout   ARPACK utility routine that prints vectors.
c     dvout   ARPACK utility routine that prints vectors.
c     zlacpy  LAPACK matrix copy routine.
c     zlahqr  LAPACK routine to compute the Schur form of an
c             upper Hessenberg matrix.
c     zlaset  LAPACK matrix initialization routine.
c     ztrevc  LAPACK routine to compute the eigenvectors of a matrix
c             in upper triangular form
c     zcopy   Level 1 BLAS that copies one vector to another. 
c     zdscal  Level 1 BLAS that scales a complex vector by a real number.
c     dznrm2  Level 1 BLAS that computes the norm of a vector.
c     
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics 
c     Rice University           
c     Houston, Texas 
c
c\SCCS Information: @(#)
c FILE: neigh.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2
c
c\Remarks
c     None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine zneigh (rnorm, n, h, ldh, ritz, bounds, 
     &                   q, ldq, workl, rwork, ierr)
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
      integer    ierr, n, ldh, ldq
      Double precision     
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Complex*16     
     &           bounds(n), h(ldh,n), q(ldq,n), ritz(n),
     &           workl(n*(n+3)) 
      Double precision 
     &           rwork(n)
c 
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex*16     
     &           one, zero
      Double precision
     &           rone
      parameter  (one = (1.0D+0, 0.0D+0), zero = (0.0D+0, 0.0D+0),
     &           rone = 1.0D+0)
c 
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
      logical    select(1)
      integer    j,  msglvl
      Complex*16     
     &           vl(1)
      Double precision
     &           temp
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   zlacpy, zlahqr, ztrevc, zcopy, 
     &           zdscal, zmout, zvout, second
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision 
     &           dznrm2
      external   dznrm2
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
      msglvl = mceigh
c 
      if (msglvl .gt. 2) then
          call zmout (logfil, n, n, h, ldh, ndigit, 
     &         '_neigh: Entering upper Hessenberg matrix H ')
      end if
c 
c     %----------------------------------------------------------%
c     | 1. Compute the eigenvalues, the last components of the   |
c     |    corresponding Schur vectors and the full Schur form T |
c     |    of the current upper Hessenberg matrix H.             |
c     |    zlahqr returns the full Schur form of H               | 
c     |    in WORKL(1:N**2), and the Schur vectors in q.         |
c     %----------------------------------------------------------%
c
      call zlacpy ('All', n, n, h, ldh, workl, n)
      call zlaset ('All', n, n, zero, one, q, ldq)
      call zlahqr (.true., .true., n, 1, n, workl, ldh, ritz,
     &             1, n, q, ldq, ierr)
      if (ierr .ne. 0) go to 9000
c
      call zcopy (n, q(n-1,1), ldq, bounds, 1)
      if (msglvl .gt. 1) then
         call zvout (logfil, n, bounds, ndigit,
     &              '_neigh: last row of the Schur matrix for H')
      end if
c
c     %----------------------------------------------------------%
c     | 2. Compute the eigenvectors of the full Schur form T and |
c     |    apply the Schur vectors to get the corresponding      |
c     |    eigenvectors.                                         |
c     %----------------------------------------------------------%
c
      call ztrevc ('Right', 'Back', select, n, workl, n, vl, n, q, 
     &             ldq, n, n, workl(n*n+1), rwork, ierr)
c
      if (ierr .ne. 0) go to 9000
c
c     %------------------------------------------------%
c     | Scale the returning eigenvectors so that their |
c     | Euclidean norms are all one. LAPACK subroutine |
c     | ztrevc returns each eigenvector normalized so  |
c     | that the element of largest magnitude has      |
c     | magnitude 1; here the magnitude of a complex   |
c     | number (x,y) is taken to be |x| + |y|.         |
c     %------------------------------------------------%
c
      do 10 j=1, n
            temp = dznrm2( n, q(1,j), 1 )
            call zdscal ( n, rone / temp, q(1,j), 1 )
   10 continue
c
      if (msglvl .gt. 1) then
         call zcopy(n, q(n,1), ldq, workl, 1)
         call zvout (logfil, n, workl, ndigit,
     &              '_neigh: Last row of the eigenvector matrix for H')
      end if
c
c     %----------------------------%
c     | Compute the Ritz estimates |
c     %----------------------------%
c
      call zcopy(n, q(n,1), n, bounds, 1)
      call zdscal(n, rnorm, bounds, 1)
c
      if (msglvl .gt. 2) then
         call zvout (logfil, n, ritz, ndigit,
     &              '_neigh: The eigenvalues of H')
         call zvout (logfil, n, bounds, ndigit,
     &              '_neigh: Ritz estimates for the eigenvalues of H')
      end if
c
      call second(t1)
      tceigh = tceigh + (t1 - t0)
c
 9000 continue
      return
c
c     %---------------%
c     | End of zneigh |
c     %---------------%
c
      end
