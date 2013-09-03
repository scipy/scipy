! Fortran 90 wrapper around lsoda

module pyodepack

    implicit none

    contains

    subroutine odeint(func, neq, y0, t, rtol, atol, tcrit, &
                      h0, hmax, hmin, &
                      dfunc, jt, ml, mu, &
                      ixpr, mxstep, mxhnil, mxordn, mxords, &
                      yout, iostate, rout, iout)

        interface odeint_user_interface

            subroutine func(n, t, y, ydot)

                integer, intent(in) :: n
                double precision, intent(in) :: t
                double precision, intent(in), dimension(n) :: y
                double precision, intent(out), dimension(n) :: ydot

            end subroutine func

            subroutine dfunc(n, t, y, ml, mu, jac, nrowjac)

                integer, intent(in) :: n
                integer, intent(in) :: nrowjac
                double precision, intent(in) :: t
                double precision, intent(in), dimension(n) :: y
                integer, intent(in) :: ml
                integer, intent(in) :: mu
                double precision, intent(out), dimension(nrowjac, n) :: jac

            end subroutine dfunc

        end interface odeint_user_interface

!             interface lsoda_interface
!             end interface lsoda_interface

        integer, intent(in) :: neq
        double precision, intent(in), dimension(neq) :: y0
        double precision, intent(in), dimension(:) :: t
        integer, intent(in) :: jt
        integer, intent(in) :: ml
        integer, intent(in) :: mu

        double precision, intent(in), dimension(:) :: rtol
        double precision, intent(in), dimension(:) :: atol

        double precision, intent(in), dimension(:) :: tcrit

        double precision, intent(in) :: h0
        double precision, intent(in) :: hmax
        ! Optionals are handled almost correctly by f2py. In _pyodepackmodule.c:
        ! if (hmin_capi == Py_None) hmin = 0.0d0; else
        ! but 0.0d0 is not a valid literal in C, so this fails in compile time.
!         double precision, intent(in), optional :: hmin = 0.0d0
        double precision, intent(in) :: hmin

        integer, intent(in) :: ixpr
        integer, intent(in) :: mxstep
        integer, intent(in) :: mxhnil
        integer, intent(in) :: mxordn
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
        double precision, intent(out), dimension(size(t), 5) :: rout
        integer, intent(out), dimension(size(t), 10) :: iout

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
        integer :: ncrit, icrit
        integer :: ierr
        integer :: lrn, lrs, lmat
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

        ! Compute size of the work arrays and allocate them
        if ((jt == 1) .or. (jt == 2)) then
            lmat = neq ** 2 + 2
        else if ((jt == 4) .or. (jt == 5)) then
            lmat = (2 * ml + mu + 1) * neq + 2
        else
            iostate = -3
            return
        end if
        lrn = 20 + (mxordn + 4) * neq
        lrs = 20 + (mxords + 1) * neq + 3 * neq + lmat
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
        rwork(5:7) = (/h0, hmax, hmin/)
        iwork(1:2) = (/ml, mu/)
        iwork(5:9) = (/ixpr, mxstep, mxhnil, mxordn, mxords/)

        itask = 1
        ! This block will run if tcrit is not NaN (because NaN /= NaN)
        if (all(tcrit == tcrit)) then
            itask = 4
            ncrit = size(tcrit)
            icrit = 1
            rwork(1) = tcrit(icrit)
        end if

        do ii = 2, size(t)
            if ((itask == 4) .and. (t(ii) > tcrit(icrit))) then
                if (icrit < ncrit) then
                    icrit = icrit + 1
                    rwork(1) = tcrit(icrit)
                else
                    itask = 1
                end if
            end if
            call lsoda(func, neq, y, t_, t(ii), itol, rtol, atol, itask, &
                       iostate, 1, rwork, lrw, iwork, liw, &
                       dfunc, jt)

            yout(ii, :) = y
            rout(ii, :) = rwork(11:15)
            iout(ii, :) = iwork(11:20)

            if (iostate < 0) exit
            ! TODO: iostate 3
            ! changes are allowed in
            ! neq, itol, rtol, atol, iopt, lrw, liw, jt, ml, mu,
            ! and any optional inputs except h0, mxordn, and mxords.
        end do

        deallocate(rwork)
        deallocate(iwork)

    end subroutine odeint

end module pyodepack
