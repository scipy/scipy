c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: ssapps
c
c\Description:
c  Given the Arnoldi factorization
c
c     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
c
c  apply NP shifts implicitly resulting in
c
c     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
c
c  where Q is an orthogonal matrix of order KEV+NP. Q is the product of 
c  rotations resulting from the NP bulge chasing sweeps.  The updated Arnoldi 
c  factorization becomes:
c
c     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
c
c\Usage:
c  call ssapps
c     ( N, KEV, NP, SHIFT, V, LDV, H, LDH, RESID, Q, LDQ, WORKD )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Problem size, i.e. dimension of matrix A.
c
c  KEV     Integer.  (INPUT)
c          INPUT: KEV+NP is the size of the input matrix H.
c          OUTPUT: KEV is the size of the updated matrix HNEW.
c
c  NP      Integer.  (INPUT)
c          Number of implicit shifts to be applied.
c
c  SHIFT   Real array of length NP.  (INPUT)
c          The shifts to be applied.
c
c  V       Real N by (KEV+NP) array.  (INPUT/OUTPUT)
c          INPUT: V contains the current KEV+NP Arnoldi vectors.
c          OUTPUT: VNEW = V(1:n,1:KEV); the updated Arnoldi vectors
c          are in the first KEV columns of V.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  H       Real (KEV+NP) by 2 array.  (INPUT/OUTPUT)
c          INPUT: H contains the symmetric tridiagonal matrix of the
c          Arnoldi factorization with the subdiagonal in the 1st column
c          starting at H(2,1) and the main diagonal in the 2nd column.
c          OUTPUT: H contains the updated tridiagonal matrix in the 
c          KEV leading submatrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RESID   Real array of length (N).  (INPUT/OUTPUT)
c          INPUT: RESID contains the the residual vector r_{k+p}.
c          OUTPUT: RESID is the updated residual vector rnew_{k}.
c
c  Q       Real KEV+NP by KEV+NP work array.  (WORKSPACE)
c          Work array used to accumulate the rotations during the bulge
c          chase sweep.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
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
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c
c\Routines called:
c     ivout   ARPACK utility routine that prints integers. 
c     second  ARPACK utility routine for timing.
c     svout   ARPACK utility routine that prints vectors.
c     slamch  LAPACK routine that determines machine constants.
c     slartg  LAPACK Givens rotation construction routine.
c     slacpy  LAPACK matrix copy routine.
c     slaset  LAPACK matrix initialization routine.
c     sgemv   Level 2 BLAS routine for matrix vector multiplication.
c     saxpy   Level 1 BLAS that computes a vector triad.
c     scopy   Level 1 BLAS that copies one vector to another.
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
c     12/16/93: Version ' 2.4'
c
c\SCCS Information: @(#) 
c FILE: sapps.F   SID: 2.6   DATE OF SID: 3/28/97   RELEASE: 2
c
c\Remarks
c  1. In this version, each shift is applied to all the subblocks of
c     the tridiagonal matrix H and not just to the submatrix that it 
c     comes from. This routine assumes that the subdiagonal elements 
c     of H that are stored in h(1:kev+np,1) are nonegative upon input
c     and enforce this condition upon output. This version incorporates
c     deflation. See code for documentation.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine ssapps
     &   ( n, kev, np, shift, v, ldv, h, ldh, resid, q, ldq, workd )
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
     &           h(ldh,2), q(ldq,kev+np), resid(n), shift(np), 
     &           v(ldv,kev+np), workd(2*n)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Real
     &           one, zero
      parameter (one = 1.0E+0, zero = 0.0E+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i, iend, istart, itop, j, jj, kplusp, msglvl
      logical    first
      Real
     &           a1, a2, a3, a4, big, c, epsmch, f, g, r, s
      save       epsmch, first
c
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   saxpy, scopy, sscal, slacpy, slartg, slaset, svout, 
     &           ivout, second, sgemv
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Real
     &           slamch
      external   slamch
c
c     %----------------------%
c     | Intrinsics Functions |
c     %----------------------%
c
      intrinsic  abs
c
c     %----------------%
c     | Data statments |
c     %----------------%
c
      data       first / .true. /
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (first) then
         epsmch = slamch('Epsilon-Machine')
         first = .false.
      end if
      itop = 1
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c
      call second (t0)
      msglvl = msapps
c 
      kplusp = kev + np 
c 
c     %----------------------------------------------%
c     | Initialize Q to the identity matrix of order |
c     | kplusp used to accumulate the rotations.     |
c     %----------------------------------------------%
c
      call slaset ('All', kplusp, kplusp, zero, one, q, ldq)
c
c     %----------------------------------------------%
c     | Quick return if there are no shifts to apply |
c     %----------------------------------------------%
c
      if (np .eq. 0) go to 9000
c 
c     %----------------------------------------------------------%
c     | Apply the np shifts implicitly. Apply each shift to the  |
c     | whole matrix and not just to the submatrix from which it |
c     | comes.                                                   |
c     %----------------------------------------------------------%
c
      do 90 jj = 1, np
c 
         istart = itop
c
c        %----------------------------------------------------------%
c        | Check for splitting and deflation. Currently we consider |
c        | an off-diagonal element h(i+1,1) negligible if           |
c        |         h(i+1,1) .le. epsmch*( |h(i,2)| + |h(i+1,2)| )   |
c        | for i=1:KEV+NP-1.                                        |
c        | If above condition tests true then we set h(i+1,1) = 0.  |
c        | Note that h(1:KEV+NP,1) are assumed to be non negative.  |
c        %----------------------------------------------------------%
c
   20    continue
c
c        %------------------------------------------------%
c        | The following loop exits early if we encounter |
c        | a negligible off diagonal element.             |
c        %------------------------------------------------%
c
         do 30 i = istart, kplusp-1
            big   = abs(h(i,2)) + abs(h(i+1,2))
            if (h(i+1,1) .le. epsmch*big) then
               if (msglvl .gt. 0) then
                  call ivout (logfil, 1, i, ndigit, 
     &                 '_sapps: deflation at row/column no.')
                  call ivout (logfil, 1, jj, ndigit, 
     &                 '_sapps: occured before shift number.')
                  call svout (logfil, 1, h(i+1,1), ndigit, 
     &                 '_sapps: the corresponding off diagonal element')
               end if
               h(i+1,1) = zero
               iend = i
               go to 40
            end if
   30    continue
         iend = kplusp
   40    continue
c
         if (istart .lt. iend) then
c 
c           %--------------------------------------------------------%
c           | Construct the plane rotation G'(istart,istart+1,theta) |
c           | that attempts to drive h(istart+1,1) to zero.          |
c           %--------------------------------------------------------%
c
             f = h(istart,2) - shift(jj)
             g = h(istart+1,1)
             call slartg (f, g, c, s, r)
c 
c            %-------------------------------------------------------%
c            | Apply rotation to the left and right of H;            |
c            | H <- G' * H * G,  where G = G(istart,istart+1,theta). |
c            | This will create a "bulge".                           |
c            %-------------------------------------------------------%
c
             a1 = c*h(istart,2)   + s*h(istart+1,1)
             a2 = c*h(istart+1,1) + s*h(istart+1,2)
             a4 = c*h(istart+1,2) - s*h(istart+1,1)
             a3 = c*h(istart+1,1) - s*h(istart,2) 
             h(istart,2)   = c*a1 + s*a2
             h(istart+1,2) = c*a4 - s*a3
             h(istart+1,1) = c*a3 + s*a4
c 
c            %----------------------------------------------------%
c            | Accumulate the rotation in the matrix Q;  Q <- Q*G |
c            %----------------------------------------------------%
c
             do 60 j = 1, min(istart+jj,kplusp)
                a1            =   c*q(j,istart) + s*q(j,istart+1)
                q(j,istart+1) = - s*q(j,istart) + c*q(j,istart+1)
                q(j,istart)   = a1
   60        continue
c
c
c            %----------------------------------------------%
c            | The following loop chases the bulge created. |
c            | Note that the previous rotation may also be  |
c            | done within the following loop. But it is    |
c            | kept separate to make the distinction among  |
c            | the bulge chasing sweeps and the first plane |
c            | rotation designed to drive h(istart+1,1) to  |
c            | zero.                                        |
c            %----------------------------------------------%
c
             do 70 i = istart+1, iend-1
c 
c               %----------------------------------------------%
c               | Construct the plane rotation G'(i,i+1,theta) |
c               | that zeros the i-th bulge that was created   |
c               | by G(i-1,i,theta). g represents the bulge.   |
c               %----------------------------------------------%
c
                f = h(i,1)
                g = s*h(i+1,1)
c
c               %----------------------------------%
c               | Final update with G(i-1,i,theta) |
c               %----------------------------------%
c
                h(i+1,1) = c*h(i+1,1)
                call slartg (f, g, c, s, r)
c
c               %-------------------------------------------%
c               | The following ensures that h(1:iend-1,1), |
c               | the first iend-2 off diagonal of elements |
c               | H, remain non negative.                   |
c               %-------------------------------------------%
c
                if (r .lt. zero) then
                   r = -r
                   c = -c
                   s = -s
                end if
c 
c               %--------------------------------------------%
c               | Apply rotation to the left and right of H; |
c               | H <- G * H * G',  where G = G(i,i+1,theta) |
c               %--------------------------------------------%
c
                h(i,1) = r
c 
                a1 = c*h(i,2)   + s*h(i+1,1)
                a2 = c*h(i+1,1) + s*h(i+1,2)
                a3 = c*h(i+1,1) - s*h(i,2)
                a4 = c*h(i+1,2) - s*h(i+1,1)
c 
                h(i,2)   = c*a1 + s*a2
                h(i+1,2) = c*a4 - s*a3
                h(i+1,1) = c*a3 + s*a4
c 
c               %----------------------------------------------------%
c               | Accumulate the rotation in the matrix Q;  Q <- Q*G |
c               %----------------------------------------------------%
c
                do 50 j = 1, min( i+jj, kplusp )
                   a1       =   c*q(j,i) + s*q(j,i+1)
                   q(j,i+1) = - s*q(j,i) + c*q(j,i+1)
                   q(j,i)   = a1
   50           continue
c
   70        continue
c
         end if
c
c        %--------------------------%
c        | Update the block pointer |
c        %--------------------------%
c
         istart = iend + 1
c
c        %------------------------------------------%
c        | Make sure that h(iend,1) is non-negative |
c        | If not then set h(iend,1) <-- -h(iend,1) |
c        | and negate the last column of Q.         |
c        | We have effectively carried out a        |
c        | similarity on transformation H           |
c        %------------------------------------------%
c
         if (h(iend,1) .lt. zero) then
             h(iend,1) = -h(iend,1)
             call sscal(kplusp, -one, q(1,iend), 1)
         end if
c
c        %--------------------------------------------------------%
c        | Apply the same shift to the next block if there is any |
c        %--------------------------------------------------------%
c
         if (iend .lt. kplusp) go to 20
c
c        %-----------------------------------------------------%
c        | Check if we can increase the the start of the block |
c        %-----------------------------------------------------%
c
         do 80 i = itop, kplusp-1
            if (h(i+1,1) .gt. zero) go to 90
            itop  = itop + 1
   80    continue
c
c        %-----------------------------------%
c        | Finished applying the jj-th shift |
c        %-----------------------------------%
c
   90 continue
c
c     %------------------------------------------%
c     | All shifts have been applied. Check for  |
c     | more possible deflation that might occur |
c     | after the last shift is applied.         |                               
c     %------------------------------------------%
c
      do 100 i = itop, kplusp-1
         big   = abs(h(i,2)) + abs(h(i+1,2))
         if (h(i+1,1) .le. epsmch*big) then
            if (msglvl .gt. 0) then
               call ivout (logfil, 1, i, ndigit, 
     &              '_sapps: deflation at row/column no.')
               call svout (logfil, 1, h(i+1,1), ndigit, 
     &              '_sapps: the corresponding off diagonal element')
            end if
            h(i+1,1) = zero
         end if
 100  continue
c
c     %-------------------------------------------------%
c     | Compute the (kev+1)-st column of (V*Q) and      |
c     | temporarily store the result in WORKD(N+1:2*N). |
c     | This is not necessary if h(kev+1,1) = 0.         |
c     %-------------------------------------------------%
c
      if ( h(kev+1,1) .gt. zero ) 
     &   call sgemv ('N', n, kplusp, one, v, ldv,
     &                q(1,kev+1), 1, zero, workd(n+1), 1)
c 
c     %-------------------------------------------------------%
c     | Compute column 1 to kev of (V*Q) in backward order    |
c     | taking advantage that Q is an upper triangular matrix |    
c     | with lower bandwidth np.                              |
c     | Place results in v(:,kplusp-kev:kplusp) temporarily.  |
c     %-------------------------------------------------------%
c
      do 130 i = 1, kev
         call sgemv ('N', n, kplusp-i+1, one, v, ldv,
     &               q(1,kev-i+1), 1, zero, workd, 1)
         call scopy (n, workd, 1, v(1,kplusp-i+1), 1)
  130 continue
c
c     %-------------------------------------------------%
c     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
c     %-------------------------------------------------%
c
      call slacpy ('All', n, kev, v(1,np+1), ldv, v, ldv)
c 
c     %--------------------------------------------%
c     | Copy the (kev+1)-st column of (V*Q) in the |
c     | appropriate place if h(kev+1,1) .ne. zero. |
c     %--------------------------------------------%
c
      if ( h(kev+1,1) .gt. zero ) 
     &     call scopy (n, workd(n+1), 1, v(1,kev+1), 1)
c 
c     %-------------------------------------%
c     | Update the residual vector:         |
c     |    r <- sigmak*r + betak*v(:,kev+1) |
c     | where                               |
c     |    sigmak = (e_{kev+p}'*Q)*e_{kev}  |
c     |    betak = e_{kev+1}'*H*e_{kev}     |
c     %-------------------------------------%
c
      call sscal (n, q(kplusp,kev), resid, 1)
      if (h(kev+1,1) .gt. zero) 
     &   call saxpy (n, h(kev+1,1), v(1,kev+1), 1, resid, 1)
c
      if (msglvl .gt. 1) then
         call svout (logfil, 1, q(kplusp,kev), ndigit, 
     &      '_sapps: sigmak of the updated residual vector')
         call svout (logfil, 1, h(kev+1,1), ndigit, 
     &      '_sapps: betak of the updated residual vector')
         call svout (logfil, kev, h(1,2), ndigit, 
     &      '_sapps: updated main diagonal of H for next iteration')
         if (kev .gt. 1) then
         call svout (logfil, kev-1, h(2,1), ndigit, 
     &      '_sapps: updated sub diagonal of H for next iteration')
         end if
      end if
c
      call second (t1)
      tsapps = tsapps + (t1 - t0)
c 
 9000 continue 
      return
c
c     %---------------%
c     | End of ssapps |
c     %---------------%
c
      end
