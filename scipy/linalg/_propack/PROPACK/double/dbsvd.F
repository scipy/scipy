c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c

      subroutine dbsvdstep(jobu,jobv,m,n,k,sigma,D,E,U,ldu,V,ldv)
c
c     Perform one implicit LQ SVD sweep with shift SIGMA.
c

c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      character*1 jobu,jobv
      integer m,n,k,ldu,ldv
      double precision D(*),E(*),U(ldu,*),V(ldv,*),sigma

c     %------------%
c     | Parameters |
c     %------------%
      double precision one, zero
      parameter(one = 1.0, zero = 0.0)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i
      double precision c,s,x,y,r
      logical dou,dov

      logical lsame
      external lsame,dlartg,drot
      
c-------------------- Here begins executable code ---------------------

      if (k.le.1) return

      dou = lsame(jobu,'y')
      dov = lsame(jobv,'y')

c     Compute the initial rotation based on B*B^T-sigma^2
      x = D(1)*D(1) - sigma*sigma
      y = E(1)*D(1)

c     Chase the "bulge" down the lower bidiagonal with Givens rotations.
c     Below 'y' is the "bulge" and 'x' is the element used to eliminate it.
      do i=1,k-1
        if (i.gt.1) then
           call dlartg(x,y,c,s,E(i-1))
        else
           call dlartg(x,y,c,s,r)
        endif
        x = c*D(i) + s*E(i)
        E(i) = -s*D(i) + c*E(i)
        D(i) = x
        y = s*D(i+1)
        D(i+1) = c*D(i+1)

        if (dou .and. m.gt.0) then
           call drot(m,U(1,i),1,U(1,i+1),1,c,s)
        endif

        call dlartg(x,y,c,s,D(i))
        x = c*E(i) + s*D(i+1)
        D(i+1) = -s*E(i) + c*D(i+1)
        E(i) = x
        y = s*E(i+1)
        E(i+1) = c*E(i+1)

        if (dov .and. n.gt.0) then
           call drot(n,V(1,i),1,V(1,i+1),1,c,s)
        endif
      enddo
      call dlartg(x,y,c,s,E(k-1))
      x = c*D(k) + s*E(k)
      E(k) = -s*D(k) + c*E(k)
      D(k) = x
      if (dou .and. m.gt.0) then
         call drot(m,U(1,k),1,U(1,k+1),1,c,s)
      endif
      return
      end




      subroutine dbdqr(ignorelast, jobq, n, D, E, c1, c2, Qt, ldq)
      implicit none

c Compute QR factorization B = Q*R of (n+1) x n lower bidiagonal matrix 
c with diagonal elements d(1)...d(n) and first subdiagonal elements
c e(1)...e(n). On return [0 ... 0 c1 c2]' = Q'*[0 ... 0 1]'.
c If ignorelast.eq..true. then e(n) is assumed to be zero.
c
c If jobq=='Y' then on return Qt contains Q^T.

c     %------------%
c     | Parameters |
c     %------------%
      character*1 jobq
      logical ignorelast
      integer n,ldq
      double precision D(*),E(*),c1,c2,Qt(ldq,*)
      
c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,j
      double precision cs,sn,r

c     %------------------------------------%
c     | External Functions and Subroutines |
c     %------------------------------------%
      logical lsame
      external lsame

c-------------------- Here begins executable code ---------------------

      if (n.lt.1) return
      if (lsame(jobq,'Y')) then
         do j=1,n+1
            do i=1,n+1
               Qt(i,j) = 0.0
            enddo
            Qt(j,j) = 1.0
         enddo
      endif
      do i=1,n-1
         call dlartg(d(i),e(i),cs,sn,r)
         d(i) = r
         e(i) = sn*d(i+1)
         d(i+1) = cs*d(i+1)
         if (lsame(jobq,'Y')) then
            do j=1,i
               Qt(i+1,j) = -sn*Qt(i,j)
               Qt(i,j) = cs*Qt(i,j)
            enddo
            Qt(i,i+1) = sn
            Qt(i+1,i+1) = cs
         endif
      enddo
      if (.not.ignorelast) then
         call dlartg(d(n),e(n),cs,sn,r)
         d(n) = r
         e(n) = 0.0
         c1 = sn
         c2 = cs
         if (lsame(jobq,'Y')) then
            do j=1,i
               Qt(i+1,j) = -sn*Qt(i,j)
               Qt(i,j) = cs*Qt(i,j)
            enddo
            Qt(i,i+1) = sn
            Qt(i+1,i+1) = cs
         endif
      endif
      end




      subroutine drefinebounds(n,k,theta,bound,tol,eps34)
c
c     Refine Lanczos error bounds using the gap theorem.
c     
c     Input arguments: 
c              n:     smallest dimension of original matrix
c              k:     number of Ritz values to refine
c              theta: array of Ritz values
c              bound: array of unrefined error bounds
c              tol:   clustering tolerance
c              eps34: machine epsilon to the power 3/4.

c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      integer n,k
      double precision theta(*), bound(*), tol, eps34

c     %-----------------%
c     | Local variables |
c     %-----------------%
      double precision gap
      integer i,l

c     %------------------------------------%
c     | External Functions and Subroutines |
c     %------------------------------------%
      double precision dlapy2
      external dlapy2

c-------------------- Here begins executable code ---------------------
      if (k.le.1) return
      do i=1,k
         do l=-1,1,2
            if ((l.eq.1.and.i.lt.k) .or. (l.eq.-1.and.i.gt.1)) then
               if (abs(theta(i)-theta(i+l)) .lt. eps34*(theta(i))) then
                  if (bound(i).gt.tol .and. bound(i+l).gt.tol) then
                     bound(i+l) = dlapy2(bound(i),bound(i+l))
                     bound(i) = 0.0
                  endif
               endif
            endif
         enddo         
      enddo
      do i=1,k
         if (i.lt.k .or. k.eq.n) then
c
c     We cannot compute a reliable value for the gap of the last
c     Ritz value unless we know it is an approximation to the 
c     smallest singular value (k.eq.n). In this case we can take the 
c     distance to the next bigger one as the gap, which can really 
c     save us from getting stuck on matrices with a single isolated tiny 
c     singular value.
c
            if (i.eq.1) then
               gap = abs(theta(i)-theta(i+1))-max(bound(i),bound(i+1))
            else if (i.eq.n) then
               gap = abs(theta(i-1)-theta(i))-max(bound(i-1),bound(i))
            else
               gap = abs(theta(i)-theta(i+1))-max(bound(i),bound(i+1)) 
               gap = min(gap,abs(theta(i-1) - theta(i)) - 
     c              max(bound(i-1),bound(i)))
            endif
            if (gap.gt.bound(i)) then
               bound(i) = bound(i) * (bound(i)/gap)
            endif
         endif
      enddo
      end

