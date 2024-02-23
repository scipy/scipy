      subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,stpmin,
     +                  stpmax)
      logical brackt
      double precision stx, fx, dx, sty, fy, dy, stp, fp, dp, stpmin,
     +                 stpmax
c     **********
c
c     Subroutine dcstep
c
c     This subroutine computes a safeguarded step for a search
c     procedure and updates an interval that contains a step that
c     satisfies a sufficient decrease and a curvature condition.
c
c     The parameter stx contains the step with the least function
c     value. If brackt is set to .true. then a minimizer has
c     been bracketed in an interval with endpoints stx and sty.
c     The parameter stp contains the current step.
c     The subroutine assumes that if brackt is set to .true. then
c
c           min(stx,sty) < stp < max(stx,sty),
c
c     and that the derivative at stx is negative in the direction
c     of the step.
c
c     The subroutine statement is
c
c       subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
c                         stpmin,stpmax)
c
c     where
c
c       stx is a double precision variable.
c         On entry stx is the best step obtained so far and is an
c            endpoint of the interval that contains the minimizer.
c         On exit stx is the updated best step.
c
c       fx is a double precision variable.
c         On entry fx is the function at stx.
c         On exit fx is the function at stx.
c
c       dx is a double precision variable.
c         On entry dx is the derivative of the function at
c            stx. The derivative must be negative in the direction of
c            the step, that is, dx and stp - stx must have opposite
c            signs.
c         On exit dx is the derivative of the function at stx.
c
c       sty is a double precision variable.
c         On entry sty is the second endpoint of the interval that
c            contains the minimizer.
c         On exit sty is the updated endpoint of the interval that
c            contains the minimizer.
c
c       fy is a double precision variable.
c         On entry fy is the function at sty.
c         On exit fy is the function at sty.
c
c       dy is a double precision variable.
c         On entry dy is the derivative of the function at sty.
c         On exit dy is the derivative of the function at the exit sty.
c
c       stp is a double precision variable.
c         On entry stp is the current step. If brackt is set to .true.
c            then on input stp must be between stx and sty.
c         On exit stp is a new trial step.
c
c       fp is a double precision variable.
c         On entry fp is the function at stp
c         On exit fp is unchanged.
c
c       dp is a double precision variable.
c         On entry dp is the the derivative of the function at stp.
c         On exit dp is unchanged.
c
c       brackt is an logical variable.
c         On entry brackt specifies if a minimizer has been bracketed.
c            Initially brackt must be set to .false.
c         On exit brackt specifies if a minimizer has been bracketed.
c            When a minimizer is bracketed brackt is set to .true.
c
c       stpmin is a double precision variable.
c         On entry stpmin is a lower bound for the step.
c         On exit stpmin is unchanged.
c
c       stpmax is a double precision variable.
c         On entry stpmax is an upper bound for the step.
c         On exit stpmax is unchanged.
c
c     MINPACK-1 Project. June 1983
c     Argonne National Laboratory.
c     Jorge J. More' and David J. Thuente.
c
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick and Jorge J. More'.
c
c     **********
      double precision zero, p66, two, three
      parameter (zero=0.0d0,p66=0.66d0,two=2.0d0,three=3.0d0)

      double precision gamma, p, q, r, s, sgnd, stpc, stpf, stpq, theta

      sgnd = dp*(dx/abs(dx))

c     First case: A higher function value. The minimum is bracketed.
c     If the cubic step is closer to stx than the quadratic step, the
c     cubic step is taken, otherwise the average of the cubic and
c     quadratic steps is taken.

      if (fp .gt. fx) then
         theta = three*(fx-fp)/(stp-stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))
         gamma = s*sqrt((theta/s)**2-(dx/s)*(dp/s))
         if (stp .lt. stx) gamma = -gamma
         p = (gamma-dx) + theta
         q = ((gamma-dx)+gamma) + dp
         r = p/q
         stpc = stx + r*(stp-stx)
         stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))/two)*(stp-stx)
         if (abs(stpc-stx) .lt. abs(stpq-stx)) then
            stpf = stpc
         else
            stpf = stpc + (stpq-stpc)/two
         end if
         brackt = .true.

c     Second case: A lower function value and derivatives of opposite
c     sign. The minimum is bracketed. If the cubic step is farther from
c     stp than the secant step, the cubic step is taken, otherwise the
c     secant step is taken.

      else if (sgnd .lt. zero) then
         theta = three*(fx-fp)/(stp-stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))
         gamma = s*sqrt((theta/s)**2-(dx/s)*(dp/s))
         if (stp .gt. stx) gamma = -gamma
         p = (gamma-dp) + theta
         q = ((gamma-dp)+gamma) + dx
         r = p/q
         stpc = stp + r*(stx-stp)
         stpq = stp + (dp/(dp-dx))*(stx-stp)
         if (abs(stpc-stp) .gt. abs(stpq-stp)) then
            stpf = stpc
         else
            stpf = stpq
         end if
         brackt = .true.

c     Third case: A lower function value, derivatives of the same sign,
c     and the magnitude of the derivative decreases.

      else if (abs(dp) .lt. abs(dx)) then

c        The cubic step is computed only if the cubic tends to infinity
c        in the direction of the step or if the minimum of the cubic
c        is beyond stp. Otherwise the cubic step is defined to be the
c        secant step.

         theta = three*(fx-fp)/(stp-stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))

c        The case gamma = 0 only arises if the cubic does not tend
c        to infinity in the direction of the step.

         gamma = s*sqrt(max(zero,(theta/s)**2-(dx/s)*(dp/s)))
         if (stp .gt. stx) gamma = -gamma
         p = (gamma-dp) + theta
         q = (gamma+(dx-dp)) + gamma
         r = p/q
         if (r .lt. zero .and. gamma .ne. zero) then
            stpc = stp + r*(stx-stp)
         else if (stp .gt. stx) then
            stpc = stpmax
         else
            stpc = stpmin
         end if
         stpq = stp + (dp/(dp-dx))*(stx-stp)

         if (brackt) then

c           A minimizer has been bracketed. If the cubic step is
c           closer to stp than the secant step, the cubic step is
c           taken, otherwise the secant step is taken.

            if (abs(stpc-stp) .lt. abs(stpq-stp)) then
               stpf = stpc
            else
               stpf = stpq
            end if
            if (stp .gt. stx) then
               stpf = min(stp+p66*(sty-stp),stpf)
            else
               stpf = max(stp+p66*(sty-stp),stpf)
            end if
         else

c           A minimizer has not been bracketed. If the cubic step is
c           farther from stp than the secant step, the cubic step is
c           taken, otherwise the secant step is taken.

            if (abs(stpc-stp) .gt. abs(stpq-stp)) then
               stpf = stpc
            else
               stpf = stpq
            end if
            stpf = min(stpmax,stpf)
            stpf = max(stpmin,stpf)
         end if

c     Fourth case: A lower function value, derivatives of the same sign,
c     and the magnitude of the derivative does not decrease. If the
c     minimum is not bracketed, the step is either stpmin or stpmax,
c     otherwise the cubic step is taken.

      else
         if (brackt) then
            theta = three*(fp-fy)/(sty-stp) + dy + dp
            s = max(abs(theta),abs(dy),abs(dp))
            gamma = s*sqrt((theta/s)**2-(dy/s)*(dp/s))
            if (stp .gt. sty) gamma = -gamma
            p = (gamma-dp) + theta
            q = ((gamma-dp)+gamma) + dy
            r = p/q
            stpc = stp + r*(sty-stp)
            stpf = stpc
         else if (stp .gt. stx) then
            stpf = stpmax
         else
            stpf = stpmin
         end if
      end if

c     Update the interval which contains a minimizer.

      if (fp .gt. fx) then
         sty = stp
         fy = fp
         dy = dp
      else
         if (sgnd .lt. zero) then
            sty = stx
            fy = fx
            dy = dx
         end if
         stx = stp
         fx = fp
         dx = dp
      end if

c     Compute the new step.

      stp = stpf

      end
