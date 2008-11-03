      double complex function wzdotc(n, zx, incx, zy, incy)
      double complex zx(*), zy(*), z
      integer n, incx, incy
      
      call veclib_zdotc(n, zx, incx, zy, incy, z)
      
      wzdotc = z
      return
      end
      
      double complex function wzdotu(n, zx, incx, zy, incy)
      double complex zx(*), zy(*), z
      integer n, incx, incy
      
      call veclib_zdotu(n, zx, incx, zy, incy, z)
      
      wzdotu = z
      return
      end
      
      complex function wcdotc(n, cx, incx, cy, incy)
      complex cx(*), cy(*), c
      integer n, incx, incy
      
      call veclib_cdotc(n, cx, incx, cy, incy, c)
      
      wcdotc = c
      return
      end
      
      complex function wcdotu(n, cx, incx, cy, incy)
      complex cx(*), cy(*), c
      integer n, incx, incy
      
      call veclib_cdotu(n, cx, incx, cy, incy, c)
      
      wcdotu = c
      return
      end

      complex function wcladiv(x, y)
      complex x, y, z
      
      call cladiv(z, x, y)
      wcladiv = z
      return
      end

      double complex function wzladiv(x, y)
      double complex x, y, z
      
      call zladiv(z, x, y)
      wzladiv = z
      return
      end
