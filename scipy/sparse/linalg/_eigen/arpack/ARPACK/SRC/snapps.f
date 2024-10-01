c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: snapps
c
c\Description:
c  Given the Arnoldi factorization
c
c     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
c
c  apply NP implicit shifts resulting in
c
c     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
c
c  where Q is an orthogonal matrix which is the product of rotations
c  and reflections resulting from the NP bulge change sweeps.
c  The updated Arnoldi factorization becomes:
c
c     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
c
c\Usage:
c  call snapps
c     ( N, KEV, NP, SHIFTR, SHIFTI, V, LDV, H, LDH, RESID, Q, LDQ,
c       WORKL, WORKD )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Problem size, i.e. size of matrix A.
c
c  KEV     Integer.  (INPUT/OUTPUT)
c          KEV+NP is the size of the input matrix H.
c          KEV is the size of the updated matrix HNEW.  KEV is only
c          updated on output when fewer than NP shifts are applied in
c          order to keep the conjugate pair together.
c
c  NP      Integer.  (INPUT)
c          Number of implicit shifts to be applied.
c
c  SHIFTR, Real array of length NP.  (INPUT)
c  SHIFTI  Real and imaginary part of the shifts to be applied.
c          Upon, entry to snapps, the shifts must be sorted so that the
c          conjugate pairs are in consecutive locations.
c
c  V       Real N by (KEV+NP) array.  (INPUT/OUTPUT)
c          On INPUT, V contains the current KEV+NP Arnoldi vectors.
c          On OUTPUT, V contains the updated KEV Arnoldi vectors
c          in the first KEV columns of V.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  H       Real (KEV+NP) by (KEV+NP) array.  (INPUT/OUTPUT)
c          On INPUT, H contains the current KEV+NP by KEV+NP upper
c          Hessenber matrix of the Arnoldi factorization.
c          On OUTPUT, H contains the updated KEV by KEV upper Hessenberg
c          matrix in the KEV leading submatrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RESID   Real array of length N.  (INPUT/OUTPUT)
c          On INPUT, RESID contains the the residual vector r_{k+p}.
c          On OUTPUT, RESID is the update residual vector rnew_{k}
c          in the first KEV locations.
c
c  Q       Real KEV+NP by KEV+NP work array.  (WORKSPACE)
c          Work array used to accumulate the rotations and reflections
c          during the bulge chase sweep.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Real work array of length (KEV+NP).  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.
c
c  WORKD   Real work array of length 2*N.  (WORKSPACE)
c          Distributed array used in the application of the accumulated
c          orthogonal matrix Q.
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
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c
c\Routines called:
c     ivout   ARPACK utility routine that prints integers.
c     arscnd  ARPACK utility routine for timing.
c     smout   ARPACK utility routine that prints matrices.
c     svout   ARPACK utility routine that prints vectors.
c     slabad  LAPACK routine that computes machine constants.
c     slacpy  LAPACK matrix copy routine.
c     slamch  LAPACK routine that determines machine constants.
c     slanhs  LAPACK routine that computes various norms of a matrix.
c     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     slarf   LAPACK routine that applies Householder reflection to
c             a matrix.
c     slarfg  LAPACK Householder reflection construction routine.
c     slartg  LAPACK Givens rotation construction routine.
c     slaset  LAPACK matrix initialization routine.
c     sgemv   Level 2 BLAS routine for matrix vector multiplication.
c     saxpy   Level 1 BLAS that computes a vector triad.
c     scopy   Level 1 BLAS that copies one vector to another .
c     sscal   Level 1 BLAS that scales a vector.
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
c FILE: napps.F   SID: 2.4   DATE OF SID: 3/28/97   RELEASE: 2
c
c\Remarks
c  1. In this version, each shift is applied to all the sublocks of
c     the Hessenberg matrix H and not just to the submatrix that it
c     comes from. Deflation as in LAPACK routine slahqr (QR algorithm
c     for upper Hessenberg matrices ) is used.
c     The subdiagonals of H are enforced to be non-negative.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine snapps
     &   ( n, kev, np, shiftr, shifti, v, ldv, h, ldh, resid, q, ldq,
     &     workl, workd )
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
      integer    kev, ldh, ldq, ldv, n, np
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Real
     &           h(ldh,kev+np), resid(n), shifti(np), shiftr(np),
     &           v(ldv,kev+np), q(ldq,kev+np), workd(2*n), workl(kev+np)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Real
     &           one, zero
      parameter (one = 1.0E+0, zero = 0.0E+0)
c
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
      integer    i, iend, ir, istart, j, jj, kplusp, msglvl, nr
      logical    cconj, first
      Real
     &           c, f, g, h11, h12, h21, h22, h32, ovfl, r, s, sigmai,
     &           sigmar, smlnum, ulp, unfl, u(3), t, tau, tst1
      save       first, ovfl, smlnum, ulp, unfl
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   saxpy, scopy, sscal, slacpy, slarfg, slarf,
     &           slaset, slabad, arscnd, slartg
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Real
     &           slamch, slanhs, slapy2
      external   slamch, slanhs, slapy2
c
c     %----------------------%
c     | Intrinsics Functions |
c     %----------------------%
c
      intrinsic  abs, max, min
c
c     %----------------%
c     | Data statements |
c     %----------------%
c
      data       first / .true. /
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (first) then
c
c        %-----------------------------------------------%
c        | Set machine-dependent constants for the       |
c        | stopping criterion. If norm(H) <= sqrt(OVFL), |
c        | overflow should not occur.                    |
c        | REFERENCE: LAPACK subroutine slahqr           |
c        %-----------------------------------------------%
c
         unfl = slamch( 'safe minimum' )
         ovfl = one / unfl
         call slabad( unfl, ovfl )
         ulp = slamch( 'precision' )
         smlnum = unfl*( n / ulp )
         first = .false.
      end if
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c
      call arscnd (t0)
      msglvl = mnapps
      kplusp = kev + np
c
c     %--------------------------------------------%
c     | Initialize Q to the identity to accumulate |
c     | the rotations and reflections              |
c     %--------------------------------------------%
c
      call slaset ('All', kplusp, kplusp, zero, one, q, ldq)
c
c     %----------------------------------------------%
c     | Quick return if there are no shifts to apply |
c     %----------------------------------------------%
c
      if (np .eq. 0) go to 9000
c
c     %----------------------------------------------%
c     | Chase the bulge with the application of each |
c     | implicit shift. Each shift is applied to the |
c     | whole matrix including each block.           |
c     %----------------------------------------------%
c
      cconj = .false.
      do 110 jj = 1, np
         sigmar = shiftr(jj)
         sigmai = shifti(jj)
c
         if (msglvl .gt. 2 ) then
            call ivout (logfil, 1, [jj], ndigit,
     &               '_napps: shift number.')
            call svout (logfil, 1, [sigmar], ndigit,
     &               '_napps: The real part of the shift ')
            call svout (logfil, 1, [sigmai], ndigit,
     &               '_napps: The imaginary part of the shift ')
         end if
c
c        %-------------------------------------------------%
c        | The following set of conditionals is necessary  |
c        | in order that complex conjugate pairs of shifts |
c        | are applied together or not at all.             |
c        %-------------------------------------------------%
c
         if ( cconj ) then
c
c           %-----------------------------------------%
c           | cconj = .true. means the previous shift |
c           | had non-zero imaginary part.            |
c           %-----------------------------------------%
c
            cconj = .false.
            go to 110
         else if ( jj .lt. np .and. abs( sigmai ) .gt. zero ) then
c
c           %------------------------------------%
c           | Start of a complex conjugate pair. |
c           %------------------------------------%
c
            cconj = .true.
         else if ( jj .eq. np .and. abs( sigmai ) .gt. zero ) then
c
c           %----------------------------------------------%
c           | The last shift has a nonzero imaginary part. |
c           | Don't apply it; thus the order of the        |
c           | compressed H is order KEV+1 since only np-1  |
c           | were applied.                                |
c           %----------------------------------------------%
c
            kev = kev + 1
            go to 110
         end if
         istart = 1
   20    continue
c
c        %--------------------------------------------------%
c        | if sigmai = 0 then                               |
c        |    Apply the jj-th shift ...                     |
c        | else                                             |
c        |    Apply the jj-th and (jj+1)-th together ...    |
c        |    (Note that jj < np at this point in the code) |
c        | end                                              |
c        | to the current block of H. The next do loop      |
c        | determines the current block ;                   |
c        %--------------------------------------------------%
c
         do 30 i = istart, kplusp-1
c
c           %----------------------------------------%
c           | Check for splitting and deflation. Use |
c           | a standard test as in the QR algorithm |
c           | REFERENCE: LAPACK subroutine slahqr    |
c           %----------------------------------------%
c
            tst1 = abs( h( i, i ) ) + abs( h( i+1, i+1 ) )
            if( tst1.eq.zero )
     &         tst1 = slanhs( '1', kplusp-jj+1, h, ldh, workl )
            if( abs( h( i+1,i ) ).le.max( ulp*tst1, smlnum ) ) then
               if (msglvl .gt. 0) then
                  call ivout (logfil, 1, [i], ndigit,
     &                 '_napps: matrix splitting at row/column no.')
                  call ivout (logfil, 1, [jj], ndigit,
     &                 '_napps: matrix splitting with shift number.')
                  call svout (logfil, 1, h(i+1,i), ndigit,
     &                 '_napps: off diagonal element.')
               end if
               iend = i
               h(i+1,i) = zero
               go to 40
            end if
   30    continue
         iend = kplusp
   40    continue
c
         if (msglvl .gt. 2) then
             call ivout (logfil, 1, [istart], ndigit,
     &                   '_napps: Start of current block ')
             call ivout (logfil, 1, [iend], ndigit,
     &                   '_napps: End of current block ')
         end if
c
c        %------------------------------------------------%
c        | No reason to apply a shift to block of order 1 |
c        %------------------------------------------------%
c
         if ( istart .eq. iend ) go to 100
c
c        %------------------------------------------------------%
c        | If istart + 1 = iend then no reason to apply a       |
c        | complex conjugate pair of shifts on a 2 by 2 matrix. |
c        %------------------------------------------------------%
c
         if ( istart + 1 .eq. iend .and. abs( sigmai ) .gt. zero )
     &      go to 100
c
         h11 = h(istart,istart)
         h21 = h(istart+1,istart)
         if ( abs( sigmai ) .le. zero ) then
c
c           %---------------------------------------------%
c           | Real-valued shift ==> apply single shift QR |
c           %---------------------------------------------%
c
            f = h11 - sigmar
            g = h21
c
            do 80 i = istart, iend-1
c
c              %-----------------------------------------------------%
c              | Construct the plane rotation G to zero out the bulge |
c              %-----------------------------------------------------%
c
               call slartg (f, g, c, s, r)
               if (i .gt. istart) then
c
c                 %-------------------------------------------%
c                 | The following ensures that h(1:iend-1,1), |
c                 | the first iend-2 off diagonal of elements |
c                 | H, remain non negative.                   |
c                 %-------------------------------------------%
c
                  if (r .lt. zero) then
                     r = -r
                     c = -c
                     s = -s
                  end if
                  h(i,i-1) = r
                  h(i+1,i-1) = zero
               end if
c
c              %---------------------------------------------%
c              | Apply rotation to the left of H;  H <- G'*H |
c              %---------------------------------------------%
c
               do 50 j = i, kplusp
                  t        =  c*h(i,j) + s*h(i+1,j)
                  h(i+1,j) = -s*h(i,j) + c*h(i+1,j)
                  h(i,j)   = t
   50          continue
c
c              %---------------------------------------------%
c              | Apply rotation to the right of H;  H <- H*G |
c              %---------------------------------------------%
c
               do 60 j = 1, min(i+2,iend)
                  t        =  c*h(j,i) + s*h(j,i+1)
                  h(j,i+1) = -s*h(j,i) + c*h(j,i+1)
                  h(j,i)   = t
   60          continue
c
c              %----------------------------------------------------%
c              | Accumulate the rotation in the matrix Q;  Q <- Q*G |
c              %----------------------------------------------------%
c
               do 70 j = 1, min( i+jj, kplusp )
                  t        =   c*q(j,i) + s*q(j,i+1)
                  q(j,i+1) = - s*q(j,i) + c*q(j,i+1)
                  q(j,i)   = t
   70          continue
c
c              %---------------------------%
c              | Prepare for next rotation |
c              %---------------------------%
c
               if (i .lt. iend-1) then
                  f = h(i+1,i)
                  g = h(i+2,i)
               end if
   80       continue
c
c           %-----------------------------------%
c           | Finished applying the real shift. |
c           %-----------------------------------%
c
         else
c
c           %----------------------------------------------------%
c           | Complex conjugate shifts ==> apply double shift QR |
c           %----------------------------------------------------%
c
            h12 = h(istart,istart+1)
            h22 = h(istart+1,istart+1)
            h32 = h(istart+2,istart+1)
c
c           %---------------------------------------------------------%
c           | Compute 1st column of (H - shift*I)*(H - conj(shift)*I) |
c           %---------------------------------------------------------%
c
            s    = 2.0*sigmar
            t = slapy2 ( sigmar, sigmai )
            u(1) = ( h11 * (h11 - s) + t * t ) / h21 + h12
            u(2) = h11 + h22 - s
            u(3) = h32
c
            do 90 i = istart, iend-1
c
               nr = min ( 3, iend-i+1 )
c
c              %-----------------------------------------------------%
c              | Construct Householder reflector G to zero out u(1). |
c              | G is of the form I - tau*( 1 u )' * ( 1 u' ).       |
c              %-----------------------------------------------------%
c
               call slarfg ( nr, u(1), u(2), 1, tau )
c
               if (i .gt. istart) then
                  h(i,i-1)   = u(1)
                  h(i+1,i-1) = zero
                  if (i .lt. iend-1) h(i+2,i-1) = zero
               end if
               u(1) = one
c
c              %--------------------------------------%
c              | Apply the reflector to the left of H |
c              %--------------------------------------%
c
               call slarf ('Left', nr, kplusp-i+1, u, 1, tau,
     &                     h(i,i), ldh, workl)
c
c              %---------------------------------------%
c              | Apply the reflector to the right of H |
c              %---------------------------------------%
c
               ir = min ( i+3, iend )
               call slarf ('Right', ir, nr, u, 1, tau,
     &                     h(1,i), ldh, workl)
c
c              %-----------------------------------------------------%
c              | Accumulate the reflector in the matrix Q;  Q <- Q*G |
c              %-----------------------------------------------------%
c
               call slarf ('Right', kplusp, nr, u, 1, tau,
     &                     q(1,i), ldq, workl)
c
c              %----------------------------%
c              | Prepare for next reflector |
c              %----------------------------%
c
               if (i .lt. iend-1) then
                  u(1) = h(i+1,i)
                  u(2) = h(i+2,i)
                  if (i .lt. iend-2) u(3) = h(i+3,i)
               end if
c
   90       continue
c
c           %--------------------------------------------%
c           | Finished applying a complex pair of shifts |
c           | to the current block                       |
c           %--------------------------------------------%
c
         end if
c
  100    continue
c
c        %---------------------------------------------------------%
c        | Apply the same shift to the next block if there is any. |
c        %---------------------------------------------------------%
c
         istart = iend + 1
         if (iend .lt. kplusp) go to 20
c
c        %---------------------------------------------%
c        | Loop back to the top to get the next shift. |
c        %---------------------------------------------%
c
  110 continue
c
c     %--------------------------------------------------%
c     | Perform a similarity transformation that makes   |
c     | sure that H will have non negative sub diagonals |
c     %--------------------------------------------------%
c
      do 120 j=1,kev
         if ( h(j+1,j) .lt. zero ) then
              call sscal( kplusp-j+1, -one, h(j+1,j), ldh )
              call sscal( min(j+2, kplusp), -one, h(1,j+1), 1 )
              call sscal( min(j+np+1,kplusp), -one, q(1,j+1), 1 )
         end if
 120  continue
c
      do 130 i = 1, kev
c
c        %--------------------------------------------%
c        | Final check for splitting and deflation.   |
c        | Use a standard test as in the QR algorithm |
c        | REFERENCE: LAPACK subroutine slahqr        |
c        %--------------------------------------------%
c
         tst1 = abs( h( i, i ) ) + abs( h( i+1, i+1 ) )
         if( tst1.eq.zero )
     &       tst1 = slanhs( '1', kev, h, ldh, workl )
         if( h( i+1,i ) .le. max( ulp*tst1, smlnum ) )
     &       h(i+1,i) = zero
 130  continue
c
c     %-------------------------------------------------%
c     | Compute the (kev+1)-st column of (V*Q) and      |
c     | temporarily store the result in WORKD(N+1:2*N). |
c     | This is needed in the residual update since we  |
c     | cannot GUARANTEE that the corresponding entry   |
c     | of H would be zero as in exact arithmetic.      |
c     %-------------------------------------------------%
c
      if (h(kev+1,kev) .gt. zero)
     &    call sgemv ('N', n, kplusp, one, v, ldv, q(1,kev+1), 1, zero,
     &                workd(n+1), 1)
c
c     %----------------------------------------------------------%
c     | Compute column 1 to kev of (V*Q) in backward order       |
c     | taking advantage of the upper Hessenberg structure of Q. |
c     %----------------------------------------------------------%
c
      do 140 i = 1, kev
         call sgemv ('N', n, kplusp-i+1, one, v, ldv,
     &               q(1,kev-i+1), 1, zero, workd, 1)
         call scopy (n, workd, 1, v(1,kplusp-i+1), 1)
  140 continue
c
c     %-------------------------------------------------%
c     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
c     %-------------------------------------------------%
c
      call slacpy ('A', n, kev, v(1,kplusp-kev+1), ldv, v, ldv)
c
c     %--------------------------------------------------------------%
c     | Copy the (kev+1)-st column of (V*Q) in the appropriate place |
c     %--------------------------------------------------------------%
c
      if (h(kev+1,kev) .gt. zero)
     &   call scopy (n, workd(n+1), 1, v(1,kev+1), 1)
c
c     %-------------------------------------%
c     | Update the residual vector:         |
c     |    r <- sigmak*r + betak*v(:,kev+1) |
c     | where                               |
c     |    sigmak = (e_{kplusp}'*Q)*e_{kev} |
c     |    betak = e_{kev+1}'*H*e_{kev}     |
c     %-------------------------------------%
c
      call sscal (n, q(kplusp,kev), resid, 1)
      if (h(kev+1,kev) .gt. zero)
     &   call saxpy (n, h(kev+1,kev), v(1,kev+1), 1, resid, 1)
c
      if (msglvl .gt. 1) then
         call svout (logfil, 1, q(kplusp,kev), ndigit,
     &        '_napps: sigmak = (e_{kev+p}^T*Q)*e_{kev}')
         call svout (logfil, 1, h(kev+1,kev), ndigit,
     &        '_napps: betak = e_{kev+1}^T*H*e_{kev}')
         call ivout (logfil, 1, [kev], ndigit,
     &               '_napps: Order of the final Hessenberg matrix ')
         if (msglvl .gt. 2) then
            call smout (logfil, kev, kev, h, ldh, ndigit,
     &      '_napps: updated Hessenberg matrix H for next iteration')
         end if
c
      end if
c
 9000 continue
      call arscnd (t1)
      tnapps = tnapps + (t1 - t0)
c
      return
c
c     %---------------%
c     | End of snapps |
c     %---------------%
c
      end
