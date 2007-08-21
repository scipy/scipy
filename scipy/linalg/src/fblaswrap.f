      subroutine wcdotu (r, n, cx, incx, cy, incy)
      external cdotu
      complex cdotu, r
      integer n
      complex cx (*)
      integer incx
      complex cy (*)
      integer incy
      r = cdotu (n, cx, incx, cy, incy)
      end

      subroutine wzdotu (r, n, zx, incx, zy, incy)
      external zdotu
      double complex zdotu, r
      integer n
      double complex zx (*)
      integer incx
      double complex zy (*)
      integer incy
      r = zdotu (n, zx, incx, zy, incy)
      end

      subroutine wcdotc (r, n, cx, incx, cy, incy)
      external cdotc
      complex cdotc, r
      integer n
      complex cx (*)
      integer incx
      complex cy (*)
      integer incy
      r = cdotc (n, cx, incx, cy, incy)
      end

      subroutine wzdotc (r, n, zx, incx, zy, incy)
      external zdotc
      double complex zdotc, r
      integer n
      double complex zx (*)
      integer incx
      double complex zy (*)
      integer incy
      r = zdotc (n, zx, incx, zy, incy)
      end
