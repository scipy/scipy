      subroutine cdotcwrapper(ret, n, cx, incx, cy, incy)
        external wcdotc
        complex wcdotc
        complex ret(1)
        integer n
        complex cx(n)
        integer incx
        complex cy(n)
        integer incy
        ret(1) = wcdotc(n, cx, incx, cy, incy)
      end

      subroutine cdotuwrapper(ret, n, cx, incx, cy, incy)
        external wcdotu
        complex wcdotu
        complex ret(1)
        integer n
        complex cx(n)
        integer incx
        complex cy(n)
        integer incy
        ret(1) = wcdotu(n, cx, incx, cy, incy)
      end

      subroutine dasumwrapper(ret, n, dx, incx)
        external dasum
        double precision dasum
        double precision ret(1)
        integer n
        double precision dx(n)
        integer incx
        ret(1) = dasum(n, dx, incx)
      end

      subroutine ddotwrapper(ret, n, dx, incx, dy, incy)
        external ddot
        double precision ddot
        double precision ret(1)
        integer n
        double precision dx(n)
        integer incx
        double precision dy(n)
        integer incy
        ret(1) = ddot(n, dx, incx, dy, incy)
      end

      subroutine dnrm2wrapper(ret, n, x, incx)
        external dnrm2
        double precision dnrm2
        double precision ret(1)
        integer n
        double precision x(n)
        integer incx
        ret(1) = dnrm2(n, x, incx)
      end

      subroutine dzasumwrapper(ret, n, zx, incx)
        external dzasum
        double precision dzasum
        double precision ret(1)
        integer n
        complex*16 zx(n)
        integer incx
        ret(1) = dzasum(n, zx, incx)
      end

      subroutine dznrm2wrapper(ret, n, x, incx)
        external dznrm2
        double precision dznrm2
        double precision ret(1)
        integer n
        complex*16 x(n)
        integer incx
        ret(1) = dznrm2(n, x, incx)
      end

      subroutine icamaxwrapper(ret, n, cx, incx)
        external icamax
        integer icamax
        integer ret(1)
        integer n
        complex cx(n)
        integer incx
        ret(1) = icamax(n, cx, incx)
      end

      subroutine idamaxwrapper(ret, n, dx, incx)
        external idamax
        integer idamax
        integer ret(1)
        integer n
        double precision dx(n)
        integer incx
        ret(1) = idamax(n, dx, incx)
      end

      subroutine isamaxwrapper(ret, n, sx, incx)
        external isamax
        integer isamax
        integer ret(1)
        integer n
        real sx(n)
        integer incx
        ret(1) = isamax(n, sx, incx)
      end

      subroutine izamaxwrapper(ret, n, zx, incx)
        external izamax
        integer izamax
        integer ret(1)
        integer n
        complex*16 zx(n)
        integer incx
        ret(1) = izamax(n, zx, incx)
      end

      subroutine sasumwrapper(ret, n, sx, incx)
        external wsasum
        real wsasum
        real ret(1)
        integer n
        real sx(n)
        integer incx
        ret(1) = wsasum(n, sx, incx)
      end

      subroutine scasumwrapper(ret, n, cx, incx)
        external wscasum
        real wscasum
        real ret(1)
        integer n
        complex cx(n)
        integer incx
        ret(1) = wscasum(n, cx, incx)
      end

      subroutine scnrm2wrapper(ret, n, x, incx)
        external wscnrm2
        real wscnrm2
        real ret(1)
        integer n
        complex x(n)
        integer incx
        ret(1) = wscnrm2(n, x, incx)
      end

      subroutine sdotwrapper(ret, n, sx, incx, sy, incy)
        external wsdot
        real wsdot
        real ret(1)
        integer n
        real sx(n)
        integer incx
        real sy(n)
        integer incy
        ret(1) = wsdot(n, sx, incx, sy, incy)
      end

      subroutine snrm2wrapper(ret, n, x, incx)
        external wsnrm2
        real wsnrm2
        real ret(1)
        integer n
        real x(n)
        integer incx
        ret(1) = wsnrm2(n, x, incx)
      end

      subroutine zdotcwrapper(ret, n, zx, incx, zy, incy)
        external wzdotc
        complex*16 wzdotc
        complex*16 ret(1)
        integer n
        complex*16 zx(n)
        integer incx
        complex*16 zy(n)
        integer incy
        ret(1) = wzdotc(n, zx, incx, zy, incy)
      end

      subroutine zdotuwrapper(ret, n, zx, incx, zy, incy)
        external wzdotu
        complex*16 wzdotu
        complex*16 ret(1)
        integer n
        complex*16 zx(n)
        integer incx
        complex*16 zy(n)
        integer incy
        ret(1) = wzdotu(n, zx, incx, zy, incy)
      end
