! Fortran 90 wrapper around lsoda
! Author: Juan Luis Cano, Aug 2013

module pyodepack

    implicit none

    contains

    subroutine odeint(func, neq, y0, t, rtol, atol, h0, hmax, hmin, &
                      dfunc, jt, &
                      ixpr, mxstep, mxhnil, mxordn, mxords, &
                      yout, iostate)

        interface odeint_user_interface

            subroutine func(n, t, y, ydot)

                !f2py integer, intent(hide), check(len(y) >= n), depend(y) :: n = len(y)
                integer, intent(in) :: n
                double precision, intent(in) :: t
                double precision, intent(in), dimension(n) :: y
                double precision, intent(out), dimension(n) :: ydot

            end subroutine func

            subroutine dfunc(n, t, y, ml, mu, jac, nrowjac)

                !f2py integer, intent(hide), check(len(y) >= n), depend(y) :: n = len(y)
                integer, intent(in) :: n
                !f2py integer, intent(hide) :: nrowjac
                integer, intent(in) :: nrowjac
                double precision, intent(in) :: t
                double precision, intent(in), dimension(n) :: y
                !f2py intent(hide) :: ml
                integer, intent(in) :: ml
                !f2py intent(hide) :: mu
                integer, intent(in) :: mu
                double precision, intent(out), dimension(nrowjac, n) :: jac

            end subroutine dfunc

        end interface odeint_user_interface

!             interface lsoda_interface
!             end interface lsoda_interface

        !f2py intent(hide) :: neq = len(y0)
        integer, intent(in) :: neq
        double precision, intent(in), dimension(neq) :: y0
        double precision, intent(in), dimension(:) :: t
        integer, intent(in) :: jt

        double precision, intent(in), dimension(:) :: rtol
        double precision, intent(in), dimension(:) :: atol

        double precision, intent(in) :: h0
        double precision, intent(in) :: hmax
        ! Optionals are handled almost correctly by f2py. In _pyodepackmodule.c:
        ! if (hmin_capi == Py_None) hmin = 0.0d0; else
        ! but 0.0d0 is not a valid literal in C, so this fails in compile time.
!         double precision, intent(in), optional :: hmin = 0.0d0
        double precision, intent(in) :: hmin

        !f2py check((ixpr == 0) || (ixpr == 1)) :: ixpr
        integer, intent(in) :: ixpr
        !f2py check(mxstep >= 0) :: mxstep
        integer, intent(in) :: mxstep
        !f2py check(mxhnil >= 0) :: mxhnil
        integer, intent(in) :: mxhnil
        !f2py check(mxordn >= 0) :: mxordn
        integer, intent(in) :: mxordn
        !f2py check(mxords >= 0) :: mxords
        integer, intent(in) :: mxords

        ! Output
        ! The same variable is reused for lsoda istate and exit state of the
        ! subroutine. On output, one code is added:
        ! iostate = 0 -> memory error
        ! The rest of the codes mean the same as in lsoda. On output:
        ! iostate = 1 -> nothing was done
        ! iostate = 2 -> integration was successful
        ! iostate = -3 -> illegal input detected
        ! and so forth.
        double precision, intent(out), dimension(size(t), neq) :: yout
        integer, intent(out) :: iostate

        ! Parameters used for only input in lsoda
        integer :: itol
        integer :: itask

        ! Parameters used for both input and output in lsoda
        double precision, dimension(neq) :: y
        double precision :: t_
        double precision, dimension(:), allocatable :: rwork
        integer, dimension(:), allocatable :: iwork

        ! Internal variables
        integer :: nrtol, natol
        integer :: ierr
        integer :: lrn, lrs
        integer :: lrw, liw
        integer :: ii

        ! Verify input
        iostate = 1
        nrtol = size(rtol)
        natol = size(atol)
        itol = 0
        if (natol /= 1) then
            if (natol == neq) then
                itol = 1
            else
                iostate = -3
                return
            end if
        end if
        if (nrtol /= 1) then
            if (nrtol == neq) then
                itol = itol + 2
            else
                iostate = -3
                return
            end if
        end if
        itol = itol + 1

        y = y0
        yout(1, :) = y0
        t_ = t(1)

        itask = 1

        ! Compute size of the work arrays and allocate them
        lrn = 20 + 16 * neq
        lrs = 22 + 9 * neq + neq * neq  ! jt = 2
        lrw = max(lrn, lrs)
        liw = 20 + neq

        allocate(rwork(lrw), stat=ierr)
        if (ierr > 0) then
            iostate = 0
            return
        end if
        allocate(iwork(liw), stat=ierr)
        if (ierr > 0) then
            iostate = 0
            return
        end if

        ! Optional inputs
        rwork(5) = h0
        rwork(6) = hmax
        rwork(7) = hmin
        iwork(5) = ixpr
        iwork(6) = mxstep
        iwork(7) = mxhnil
        iwork(8) = mxordn
        iwork(9) = mxords

        do ii = 2, size(t)
            call lsoda(func, neq, y, t_, t(ii), itol, rtol, atol, itask, &
                        iostate, 1, rwork, lrw, iwork, liw, &
                        dfunc, jt)

            if (iostate < 0) exit
            ! TODO: iostate 3
            ! changes are allowed in
            ! neq, itol, rtol, atol, iopt, lrw, liw, jt, ml, mu,
            ! and any optional inputs except h0, mxordn, and mxords.
            yout(ii, :) = y
        end do

        deallocate(rwork)
        deallocate(iwork)

    end subroutine odeint

end module pyodepack
