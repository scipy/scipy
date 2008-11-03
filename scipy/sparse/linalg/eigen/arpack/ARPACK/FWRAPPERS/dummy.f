      double complex function wzdotc(n, zx, incx, zy, incy)
      double complex zx(*), zy(*), z, zdotc
      integer n, incx, incy
      
      zdotc(n, zx, incx, zy, incy, z)
      wzdotc = zdotc
      
      end
      
      double complex function wzdotu(n, zx, incx, zy, incy)
      double complex zx(*), zy(*), z, zdotu
      integer n, incx, incy
      
      zdotu(n, zx, incx, zy, incy, z)
      wzdotu = zdotu
      
      return
      end
      
      complex function wcdotc(n, cx, incx, cy, incy)
      complex cx(*), cy(*), c, cdotc
      integer n, incx, incy
      
      cdotc(n, cx, incx, cy, incy, c)
      wzdotc = zdotc
      
      return
      end
      
      complex function wcdotu(n, cx, incx, cy, incy)
      complex cx(*), cy(*), c, cdotu
      integer n, incx, incy
      
      cdotu(n, cx, incx, cy, incy, c)
      wzdotu = zdotu
      
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
