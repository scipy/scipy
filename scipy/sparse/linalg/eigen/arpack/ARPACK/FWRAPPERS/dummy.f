      double complex function wzdotc(n, zx, incx, zy, incy)
      double complex zx(*), zy(*), z
      double complex zdotc
      integer n, incx, incy
      
      z = zdotc(n, zx, incx, zy, incy)
      wzdotc = z
      
      end
      
      double complex function wzdotu(n, zx, incx, zy, incy)
      double complex zx(*), zy(*), z, zdotu
      integer n, incx, incy
      
      z = zdotu(n, zx, incx, zy, incy)
      wzdotu = z
      
      return
      end
      
      complex function wcdotc(n, cx, incx, cy, incy)
      complex cx(*), cy(*), c, cdotc
      integer n, incx, incy
      
      c = cdotc(n, cx, incx, cy, incy)
      wcdotc = c
      
      return
      end
      
      complex function wcdotu(n, cx, incx, cy, incy)
      complex cx(*), cy(*), c, cdotu
      integer n, incx, incy
      
      c = cdotu(n, cx, incx, cy, incy)
      wcdotu = c
      
      return
      end

      real function wsdot(n, x, incx, y, incy)
      real x(*), y(*), s, sdot
      integer n, incx, incy

      s = sdot(n, x, incx, y, incy)
      wsdot = s

      return
      end

      complex function wcladiv(x, y)
      complex x, y, z
      complex cladiv
      
      z = cladiv(x, y)
      wcladiv = z
      return
      end

      double complex function wzladiv(x, y)
      double complex x, y, z
      double complex zladiv
      
      z = zladiv(x, y)
      wzladiv = z
      return
      end
